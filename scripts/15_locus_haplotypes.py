#!/usr/bin/env python3
"""S7 — per-locus haplotype FASTAs for a local pangenome graph (pggb) of SMN and HBA.

For a diverse subset of complete (usable) HPRC R2 haplotypes, pull the locus from the inner unique anchor on one
side to the inner anchor on the other side (SMN: MARVELD2..BDP1, HBA: NPRL3..LUC7L) with remote samtools faidx
(same code as 08: local .fai/.gzi, only the locus bytes are fetched), orient it like GRCh38, and add the GRCh38 and
CHM13 reference segments. Also writes the gene layout of every path (path-relative coordinates) so the graph
figures can be annotated.

Outputs: data/pangraph/{SMN,HBA}.fa(.fai), data/pangraph/{SMN,HBA}.genes.tsv
"""
import csv
import importlib.util
import os
import random
import sys
from collections import defaultdict

from cfg import CFG, path
from common import log, read_fasta_region, revcomp, s3_to_https, status

spec = importlib.util.spec_from_file_location("t8", os.path.join(os.path.dirname(__file__), "08_truth_from_assemblies.py"))
t8 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t8)

OUT = path("pangraph")
os.makedirs(OUT, exist_ok=True)
N_HAPS = dict(CFG.get("pangraph", {}).get("n_haps", {"SMN": 16, "HBA": 48}))
INNER = {fam: (CFG["loci"][fam]["anchor_sides"][0][-1], CFG["loci"][fam]["anchor_sides"][1][0]) for fam in N_HAPS}


def ref_segment(ref, fam):
    """(chrom, start, end, strand_of_left_anchor, genes) on a reference, anchors included."""
    genes = [l.rstrip("\n").split("\t") for l in open(path("regions", f"{ref}.genes.bed"))]
    want = set(CFG["loci"][fam]["genes"]) | set(INNER[fam])
    gs = [g for g in genes if g[3] in want]
    chrom = max({g[0] for g in gs}, key=lambda c: sum(g[0] == c for g in gs))
    gs = [g for g in gs if g[0] == chrom]
    lo, hi = min(int(g[1]) for g in gs), max(int(g[2]) for g in gs)
    strand = next(g[5] for g in gs if g[3] == INNER[fam][0])
    return chrom, lo, hi, strand, [(g[3], int(g[1]) - lo, int(g[2]) - lo, g[5]) for g in gs]


def pick(fam, haps):
    """Stratified pick: every copy-number class first, then random fill (seeded)."""
    cols = {"SMN": ["SMN1", "SMN2", "SMN_del7-8"], "HBA": ["HBA_total"]}[fam]
    ok = [h for h in haps if not h[f"{fam}_flags"] and h[cols[0]] not in ("", None)]
    by = defaultdict(list)
    for h in ok:
        by[tuple(h[c] for c in cols)].append(h)
    rng = random.Random(7)
    chosen = []
    for k in sorted(by, key=lambda k: len(by[k])):
        chosen.append(rng.choice(by[k]))
    rest = [h for h in ok if h not in chosen]
    rng.shuffle(rest)
    chosen += rest[: max(0, N_HAPS[fam] - len(chosen))]
    return chosen[: max(N_HAPS[fam], len(by))]


def main():
    haps = list(csv.DictReader(open(path("truth", "haplotypes.tsv")), delimiter="\t"))
    asm = {r["assembly_name"]: r for r in csv.DictReader(open(path("meta", "hprc_r2_assemblies.csv")))}
    cat = {r["assembly_name"]: r["location"] for r in csv.DictReader(open(path("meta", "hprc_r2_cat_annotation.csv")))}
    for fam in N_HAPS:
        if os.path.exists(os.path.join(OUT, f"{fam}.fa.ok")):
            continue
        recs, layout = [], []
        for ref in ["grch38", "chm13"]:
            chrom, lo, hi, strand, gl = ref_segment(ref, fam)
            name = f"{ref.upper()}#0#{chrom}:{lo + 1}-{hi}"
            recs.append((name, read_fasta_region(path(CFG["ref_fasta"][ref]), f"{chrom}:{lo + 1}-{hi}")[1]))
            layout += [(name, *g) for g in gl]
            if ref == "grch38":
                ref_strand = strand
        for h in pick(fam, haps):
            a = asm[h["assembly"]]
            genes = t8.cat_genes(h["assembly"], s3_to_https(cat[h["assembly"]]))
            fai = t8.local_idx(s3_to_https(a["assembly_fai"]), h["assembly"], ".fai")
            gzi = t8.local_idx(s3_to_https(a["assembly_gzi"]), h["assembly"], ".gzi")
            lens = t8.fai_lengths(fai)
            alias = {n.split("#")[-1]: n for n in lens}
            core = [g for g in genes if g[5] == fam and g[4] in CFG["loci"][fam]["genes"]]
            contig = max({g[0] for g in core}, key=lambda c: sum(g[0] == c for g in core))
            gs = [g for g in genes if g[0] == contig and (g[4] in CFG["loci"][fam]["genes"] or g[4] in INNER[fam])]
            if not all(any(g[4] == x for g in gs) for x in INNER[fam]):
                log(f"pangraph {fam}: {h['assembly']} lacks an inner anchor on {contig}; skipped")
                continue
            lo, hi = min(int(g[1]) for g in gs), max(int(g[2]) for g in gs)
            seq = t8.extract(s3_to_https(a["assembly"]), fai, gzi, f"{alias.get(contig, contig)}:{lo + 1}-{hi}")
            flip = next(g[3] for g in gs if g[4] == INNER[fam][0]) != ref_strand
            L = hi - lo
            if flip:
                seq = revcomp(seq)
            tag = {"SMN": f"SMN1={h['SMN1']},SMN2={h['SMN2']},d78={h['SMN_del7-8']}", "HBA": f"HBA={h['HBA_total']}"}[fam]
            name = f"{h['sample']}#{h['hap']}#{contig}:{lo + 1}-{hi}{'(-)' if flip else ''}"
            recs.append((name, seq))
            for g in gs:
                s, e = int(g[1]) - lo, int(g[2]) - lo
                if flip:
                    s, e = L - e, L - s
                layout.append((name, g[4], s, e, ("-" if g[3] == "+" else "+") if flip else g[3]))
            layout.append((name, "#truth", 0, 0, tag))
            log(f"pangraph {fam}: {name} {len(seq) / 1e3:.0f} kb {tag}")
        with open(os.path.join(OUT, f"{fam}.fa"), "w") as fh:
            for n, s in recs:
                fh.write(f">{n}\n")
                for i in range(0, len(s), 80):
                    fh.write(s[i:i + 80] + "\n")
        with open(os.path.join(OUT, f"{fam}.genes.tsv"), "w") as fh:
            fh.write("path\tgene\tstart\tend\tstrand\n")
            for r in layout:
                fh.write("\t".join(map(str, r)) + "\n")
        open(os.path.join(OUT, f"{fam}.fa.ok"), "w").close()
        status(f"S7_pangraph_fa_{fam}", "DONE", f"{len(recs)} paths ({sum(len(s) for _, s in recs) / 1e6:.1f} Mb)")


if __name__ == "__main__":
    main()
