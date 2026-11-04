#!/usr/bin/env python3
"""S12 — exon/intron models for every gene drawn in the toolbox (data/regions/{ref}.exons.tsv).

GRCh38: GENCODE v49 basic, the transcript tagged Ensembl_canonical (fallback: longest).
CHM13:  RefSeq Liftoff v5.2, the transcript whose structure best matches the GRCh38 canonical one (same exon count,
        closest exon lengths; NM_ preferred) — "longest NM_" had picked an SMN1 isoform skipping exon 7.
        Genes missing from the Liftoff file (CBS, U2AF1, CRYAA) get their GRCh38 exons lifted with the chain (CrossMap).
GRCh38 false-duplication copies (CBSL, U2AF1L5, CRYAA2; found by sequence in 05) get the parent's exons projected
        by offset (gene lengths agree within 16 bp), flagged `projected`.
SMN1/SMN2 exons are labelled with the clinical names 1, 2a, 2b, 3–8 (canonical transcript has 9 exons).

Columns: gene, transcript, source, chrom, start (0-based), end, strand, exon_label, cds_start, cds_end ('' if non-coding)
"""
import gzip
import os
import re
import subprocess
import tempfile
from collections import defaultdict

from cfg import CFG, path
from common import log, status

SMN_LABELS = ["1", "2a", "2b", "3", "4", "5", "6", "7", "8"]


def genes_bed(ref):
    return {l.split("\t")[3]: l.rstrip("\n").split("\t") for l in open(path("regions", f"{ref}.genes.bed"))}


def gencode(names):
    tx = defaultdict(lambda: dict(exons=[], cds=[], canon=False, gene=None, chrom=None, strand=None))
    with gzip.open(path("refs", "gencode.basic.annotation.gtf.gz"), "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("transcript", "exon", "CDS"):
                continue
            g = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
            if g not in names:
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            d = tx[t]; d.update(gene=g, chrom=f[0], strand=f[6])
            if f[2] == "transcript":
                d["canon"] = "Ensembl_canonical" in f[8]
            elif f[2] == "exon":
                d["exons"].append((int(f[3]) - 1, int(f[4])))
            else:
                d["cds"].append((int(f[3]) - 1, int(f[4])))
    return pick(tx, prefer=lambda t, d: (d["canon"], sum(e - s for s, e in d["exons"])))


def liftoff_gff(names, template=None):
    tx = defaultdict(lambda: dict(exons=[], cds=[], gene=None, chrom=None, strand=None))
    with gzip.open(path("refs", "chm13v2.0.genes.gff3.gz"), "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            m = re.search(r"(?:^|;)gene=([^;]+)", f[8]); p = re.search(r"(?:^|;)Parent=([^;]+)", f[8])
            if not m or not p or m.group(1) not in names:
                continue
            d = tx[p.group(1)]; d.update(gene=m.group(1), chrom=f[0], strand=f[6])
            d["exons" if f[2] == "exon" else "cds"].append((int(f[3]) - 1, int(f[4])))
    # prefer the transcript whose structure matches the GRCh38 canonical one (same exon count, closest exon lengths):
    # "longest NM_" picked e.g. an SMN1 isoform that skips exon 7
    def score(t, d):
        tm = (template or {}).get(d["gene"])
        L = sorted(e - s for s, e in set(d["exons"]))
        if tm:
            T = sorted(e - s for s, e in tm["exons"])
            same = len(L) == len(T)
            diff = sum(abs(a - b) for a, b in zip(L, T)) if same else 10**9
            return (same, -diff, t.startswith("NM_"), sum(L))
        return (True, 0, t.startswith("NM_"), sum(L))
    return pick(tx, prefer=score)


def pick(tx, prefer):
    best = {}
    for t, d in tx.items():
        if d["exons"] and (d["gene"] not in best or prefer(t, d) > prefer(*best[d["gene"]])):
            best[d["gene"]] = (t, d)
    return {g: dict(tx=t, chrom=d["chrom"], strand=d["strand"], exons=sorted(set(d["exons"])), cds=sorted(d["cds"]))
            for g, (t, d) in best.items()}


def rows_for(g, m, source):
    ex = m["exons"]
    order = ex if m["strand"] == "+" else ex[::-1]
    labels = SMN_LABELS if g in ("SMN1", "SMN2") and len(ex) == 9 else [str(i + 1) for i in range(len(ex))]
    lab = {e: labels[i] for i, e in enumerate(order)}
    cs = min((s for s, _ in m["cds"]), default=""); ce = max((e for _, e in m["cds"]), default="")
    return [(g, m["tx"], source, m["chrom"], s, e, m["strand"], lab[(s, e)], cs, ce) for s, e in ex]


def crossmap_lift(model):
    """Lift one GRCh38 model to CHM13 exon by exon; None unless every exon lifts 1:1 to one chromosome/strand."""
    chain = path("refs", "grch38-chm13v2.chain")
    with tempfile.TemporaryDirectory() as td:
        bi, bo = os.path.join(td, "in.bed"), os.path.join(td, "out.bed")
        with open(bi, "w") as fh:
            for i, (s, e) in enumerate(model["exons"]):
                fh.write(f"{model['chrom']}\t{s}\t{e}\tx{i}\t0\t{model['strand']}\n")
            if model["cds"]:
                fh.write(f"{model['chrom']}\t{min(s for s, _ in model['cds'])}\t{min(s for s, _ in model['cds']) + 1}\tcs\t0\t{model['strand']}\n")
                fh.write(f"{model['chrom']}\t{max(e for _, e in model['cds']) - 1}\t{max(e for _, e in model['cds'])}\tce\t0\t{model['strand']}\n")
        subprocess.run(["CrossMap", "bed", chain, bi, bo], capture_output=True, text=True)
        if not os.path.exists(bo):
            return None
        got = {}
        for l in open(bo):
            f = l.rstrip("\n").split("\t")
            got[f[3]] = (f[0], int(f[1]), int(f[2]), f[5])
    n = len(model["exons"])
    if any(f"x{i}" not in got for i in range(n)) or len({v[0] for v in got.values()}) != 1:
        return None
    ex = sorted((got[f"x{i}"][1], got[f"x{i}"][2]) for i in range(n))
    cds = []
    if "cs" in got and "ce" in got:
        a, b = sorted([got["cs"][1], got["ce"][2]]) if got["cs"][1] < got["ce"][1] else sorted([got["ce"][1], got["cs"][2]])
        cds = [(a, b)]
    return dict(tx=model["tx"], chrom=got["x0"][0], strand=got["x0"][3], exons=ex, cds=cds)


def main():
    paralog = CFG["loci"]["CHR21_FALSEDUP"].get("paralog_of", {})
    for ref in ["grch38", "chm13"]:
        gb = genes_bed(ref)
        names = set(gb)
        models = gencode(names | set(paralog.values())) if ref == "grch38" else liftoff_gff(names, template=gencode(names))
        out = []
        for g in sorted(names):
            if g in models:
                out += rows_for(g, models[g], "GENCODE v49 canonical" if ref == "grch38" else "RefSeq Liftoff v5.2")
        if ref == "grch38":
            for copy, parent in paralog.items():
                if copy in gb and parent in models and copy not in models:
                    pm, c = models[parent], gb[copy]
                    pg = gb[parent]
                    off = int(c[1]) - int(pg[1])
                    m = dict(tx=pm["tx"] + "(projected)", chrom=c[0], strand=c[5],
                             exons=[(s + off, e + off) for s, e in pm["exons"]], cds=[(s + off, e + off) for s, e in pm["cds"]])
                    out += rows_for(copy, m, f"projected from {parent} by offset")
        else:
            g38 = gencode(names)
            for g in sorted(names - set(models)):
                if g in g38:
                    m = crossmap_lift(g38[g])
                    if m:
                        out += rows_for(g, m, "GENCODE v49 canonical lifted with chain")
                    else:
                        log(f"exons: {g} could not be lifted to CHM13")
        with open(path("regions", f"{ref}.exons.tsv"), "w") as fh:
            fh.write("gene\ttranscript\tsource\tchrom\tstart\tend\tstrand\texon_label\tcds_start\tcds_end\n")
            for r in out:
                fh.write("\t".join(map(str, r)) + "\n")
        got = sorted({r[0] for r in out})
        log(f"exons {ref}: {len(got)}/{len(names)} genes; missing {sorted(names - set(got))}")
    status("S12_exons", "DONE", "data/regions/{grch38,chm13}.exons.tsv")


if __name__ == "__main__":
    main()
