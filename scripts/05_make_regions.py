#!/usr/bin/env python3
"""Phase 5 — regions for both references (written to data/regions/).

  {ref}.genes.bed      gene copies (chrom, start0, end, name, locus)        name gets #k for extra copies
  {ref}.loci.bed       locus windows (span of genes ± pad)
  {ref}.bins.bed       fixed-width bins tiling each locus window (for depth profiles)
  {ref}.controls.bed   2-kb sub-windows of contiguous unique control blocks (same IDs on both refs)
  {ref}.slice.bed      merged union of loci + control blocks: what 06_slice_crams.sh extracts
  {ref}.psv.tsv        paralog-specific sites for SMN1/SMN2 and HBA1/HBA2, per gene copy
  projection_check.tsv GRCh38 SMN sites re-projected onto GRCh38 (self-test of the projection code)
"""
import gzip
import os
import re
import subprocess
import tempfile
from collections import defaultdict

from cfg import CFG, path
from common import log, project_positions, read_fasta_region, revcomp, run, status

R, OUT = path("refs"), path("regions")
os.makedirs(OUT, exist_ok=True)
REF_FA = {k: path(v) for k, v in CFG["ref_fasta"].items()}


def opener(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)


def parse_genes(fn, wanted, fmt):
    """Return {symbol: [(chrom, start0, end, strand)]} for gene features named in `wanted`."""
    hits = defaultdict(list)
    pat_gtf = re.compile(r'gene_name "([^"]+)"')
    pat_gff = re.compile(r"(?:^|;)(?:Name|gene|gene_name)=([^;]+)")
    with opener(fn) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            names = pat_gtf.findall(f[8]) if fmt == "gtf" else pat_gff.findall(f[8])
            for n in set(names):
                if n in wanted:
                    hits[n].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
    for n in hits:  # de-duplicate identical records (Name= and gene= on the same line)
        hits[n] = sorted(set(hits[n]), key=lambda x: (x[0], x[1]))
    return hits


def rescue_missing_genes(genes, loci):
    """Fill gaps in the annotations by sequence/chain so every reference gets the same gene set.

    - grch38_only_genes (GRCh38 false-duplication copies) that GENCODE no longer names: align the parent
      gene (config `paralog_of`) to its chromosome with minimap2 and take the best non-self hit.
    - genes missing from the CHM13 GFF3 (e.g. CBS is absent from the RefSeq Liftoff v5.2 file):
      lift the GRCh38 record with the GRCh38->CHM13 chain (CrossMap), keep a lift covering >= 90%.
    """
    G38 = genes["grch38"]
    for lname, L in loci.items():
        for child, parent in (L.get("paralog_of") or {}).items():
            if G38.get(child) or not G38.get(parent):
                continue
            c, s, e, strand = G38[parent][0]
            chrom_fa = os.path.join(OUT, f"_{c}.grch38.fa")
            if not os.path.exists(chrom_fa):
                with open(chrom_fa, "w") as fh:
                    fh.write(run(["samtools", "faidx", REF_FA["grch38"], c]))
            q = read_fasta_region(REF_FA["grch38"], f"{c}:{s + 1}-{e}")[1]
            paf = run(["minimap2", "-x", "asm20", "-N", "10", "-p", "0.3", chrom_fa, "/dev/stdin"], input=f">q\n{q}\n")
            best = None
            for line in paf.splitlines():
                f = line.split("\t")
                ts, te, nm, al = int(f[7]), int(f[8]), int(f[9]), int(f[10])
                if ts < e and s < te:  # the gene itself
                    continue
                if al >= 0.8 * len(q) and (best is None or nm > best[2]):
                    best = (ts, te, nm, f[4])
            if best:
                G38[child] = [(c, best[0], best[1], strand if best[3] == "+" else ("-" if strand == "+" else "+"))]
                log(f"grch38: {child} not in GENCODE; located by aligning {parent}: {c}:{best[0] + 1}-{best[1]}")
            else:
                log(f"grch38: {child} not found by sequence either")
    if "chm13" not in genes:
        return
    G13 = genes["chm13"]
    for lname, L in loci.items():
        for sym in L["genes"] + L.get("anchors", []):
            if any(r[0] == L["chrom"] for r in G13.get(sym, [])) or not G38.get(sym):
                continue
            c, s, e, strand = next((r for r in G38[sym] if r[0] == L["chrom"]), G38[sym][0])
            tmp = tempfile.mkdtemp(dir=OUT)
            src, dst = os.path.join(tmp, "g.bed"), os.path.join(tmp, "g.chm13.bed")
            write_bed(src, [(c, s, e, sym, 0, strand)])
            subprocess.run(["CrossMap", "bed", os.path.join(R, "grch38-chm13v2.chain"), src, dst],
                           check=True, capture_output=True)
            parts = [l.split("\t") for l in open(dst)] if os.path.exists(dst) else []
            parts = [p for p in parts if p[0] == L["chrom"]]
            if parts:
                lo, hi = min(int(p[1]) for p in parts), max(int(p[2]) for p in parts)
                if sum(int(p[2]) - int(p[1]) for p in parts) >= 0.9 * (e - s) and hi - lo <= 1.5 * (e - s):
                    G13[sym] = [(L["chrom"], lo, hi, parts[0][5].strip())]
                    log(f"chm13: {sym} not in GFF3; lifted from GRCh38 by chain: {L['chrom']}:{lo + 1}-{hi}")
                    continue
            log(f"chm13: {sym} not in GFF3 and chain lift failed")


def write_bed(fn, rows):
    with open(fn, "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


def merge(rows):
    rows = sorted(rows, key=lambda r: (r[0], r[1]))
    out = []
    for c, s, e, *_ in rows:
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out


def n_frac(fa, chrom, s, e):
    seq = read_fasta_region(fa, f"{chrom}:{s + 1}-{e}")[1]
    return seq.count("N") / max(1, len(seq))


def load_superdups():
    fn = os.path.join(R, "hg38.genomicSuperDups.txt.gz")
    sd = defaultdict(list)
    if os.path.exists(fn):
        with gzip.open(fn, "rt") as fh:
            for line in fh:
                f = line.split("\t")
                sd[f[1]].append((int(f[2]), int(f[3])))
    return sd


def overlaps(iv, s, e):
    return any(a < e and s < b for a, b in iv)


def main():
    loci = CFG["loci"]
    gene_sets = {"grch38": set(), "chm13": set()}
    for name, L in loci.items():
        gene_sets["grch38"] |= set(L["genes"]) | set(L.get("grch38_only_genes", [])) | set(L.get("anchors", []))
        gene_sets["chm13"] |= set(L["genes"]) | set(L.get("anchors", []))

    ann = {"grch38": (os.path.join(R, "gencode.basic.annotation.gtf.gz"), "gtf")}
    for cand in ["chm13v2.0.genes.gff3.gz", "chm13v2.0.genes.gff3"]:
        if os.path.exists(os.path.join(R, cand)):
            ann["chm13"] = (os.path.join(R, cand), "gff")
    genes = {ref: parse_genes(fn, gene_sets[ref], fmt) for ref, (fn, fmt) in ann.items()}
    rescue_missing_genes(genes, loci)

    # ---- genes / loci / bins
    for ref, G in genes.items():
        gene_rows, loci_rows, bin_rows = [], [], []
        for lname, L in loci.items():
            syms = L["genes"] + (L.get("grch38_only_genes", []) if ref == "grch38" else [])
            recs = [(sym, r) for sym in syms for r in G.get(sym, []) if r[0] == L["chrom"]]
            missing = [s for s in syms if not any(x[0] == s for x in recs)]
            if missing:
                log(f"{ref} {lname}: genes not found on {L['chrom']}: {missing}")
            if not recs:
                if not L.get("optional"):
                    log(f"ERROR {ref} {lname}: no genes found")
                continue
            for sym in syms:
                copies = [r for s, r in recs if s == sym]
                for k, (c, s, e, strand) in enumerate(copies):
                    gene_rows.append((c, s, e, sym if k == 0 else f"{sym}#{k}", lname, strand))
            c = L["chrom"]
            # one window per cluster of genes (genes > 2 Mb apart, e.g. the chr21 false duplication, get
            # separate windows instead of one huge span)
            spans = merge([(c, max(0, r[1] - L["pad_bp"]), r[2] + L["pad_bp"]) for _, r in recs])
            spans = merge([(c, s - 1_000_000, e + 1_000_000) for c, s, e in spans])
            spans = [(c, s + 1_000_000, e - 1_000_000) for c, s, e in spans]
            i = 0
            for _, lo, hi in spans:
                lo = max(0, lo)
                loci_rows.append((c, lo, hi, lname))
                for s in range(lo, hi, L["bin_bp"]):
                    bin_rows.append((c, s, min(hi, s + L["bin_bp"]), f"{lname}:bin{i:05d}", lname))
                    i += 1
            # anchors kept separately for assembly-completeness checks downstream
            for a in L.get("anchors", []):
                for k, (ac, s, e, strand) in enumerate(G.get(a, [])):
                    if ac == c:
                        gene_rows.append((ac, s, e, a if k == 0 else f"{a}#{k}", f"{lname}_anchor", strand))
        write_bed(os.path.join(OUT, f"{ref}.genes.bed"), sorted(gene_rows, key=lambda r: (r[0], r[1])))
        write_bed(os.path.join(OUT, f"{ref}.loci.bed"), loci_rows)
        write_bed(os.path.join(OUT, f"{ref}.bins.bed"), bin_rows)
        log(f"{ref}: {len(gene_rows)} gene records, {len(loci_rows)} loci, {len(bin_rows)} bins")

    # ---- control blocks (GRCh38), chosen around SMNCopyNumberCaller's normalisation windows
    C = CFG["controls"]
    sd = load_superdups()
    norm = []
    for line in open(path("meta", "SMN_region_38.bed")):
        f = line.split("\t")
        if len(f) > 4 and f[4].strip() == "norm" and f[0] not in C["exclude_chroms"]:
            mid = (int(f[1]) + int(f[2])) // 2
            norm.append((f[0], mid))
    norm.sort(key=lambda x: (int(re.sub(r"\D", "", x[0]) or 99), x[1]))
    step = max(1, len(norm) // (C["n_blocks"] * 3))
    blocks = []
    for c, mid in norm[::step]:
        s, e = mid - C["block_bp"] // 2, mid + C["block_bp"] // 2
        if s < 0 or overlaps(sd.get(c, []), s, e) or any(b[0] == c and b[1] < e + 1_000_000 and s < b[2] + 1_000_000 for b in blocks):
            continue
        if n_frac(REF_FA["grch38"], c, s, e) > 0.005:
            continue
        blocks.append((c, s, e))
    # keep n_blocks spread over the genome
    if len(blocks) > C["n_blocks"]:
        idx = [round(i * (len(blocks) - 1) / (C["n_blocks"] - 1)) for i in range(C["n_blocks"])]
        blocks = [blocks[i] for i in sorted(set(idx))]
    sub = []
    for bi, (c, s, e) in enumerate(blocks):
        for wi, ws in enumerate(range(s, e, C["subwindow_bp"])):
            sub.append((c, ws, min(e, ws + C["subwindow_bp"]), f"ctrl{bi:02d}_{wi:03d}"))
    log(f"control blocks: {len(blocks)} ({len(sub)} sub-windows)")

    # lift sub-windows to CHM13 with CrossMap; keep 1:1 lifts of near-equal length
    tmp = tempfile.mkdtemp(dir=OUT)
    src, dst = os.path.join(tmp, "ctrl.grch38.bed"), os.path.join(tmp, "ctrl.chm13.bed")
    write_bed(src, sub)
    subprocess.run(["CrossMap", "bed", os.path.join(R, "grch38-chm13v2.chain"), src, dst],
                   check=True, capture_output=True)
    lifted = defaultdict(list)
    for line in open(dst):
        f = line.rstrip("\n").split("\t")
        lifted[f[3].split("(")[0]].append((f[0], int(f[1]), int(f[2])))
    keep_g, keep_c = [], []
    for c, s, e, wid in sub:
        L = lifted.get(wid, [])
        if len(L) == 1 and abs((L[0][2] - L[0][1]) - (e - s)) <= 0.05 * (e - s):
            keep_g.append((c, s, e, wid))
            keep_c.append((*L[0], wid))
    write_bed(os.path.join(OUT, "grch38.controls.bed"), keep_g)
    write_bed(os.path.join(OUT, "chm13.controls.bed"), sorted(keep_c, key=lambda r: (r[0], r[1])))
    log(f"control sub-windows lifted 1:1 to CHM13: {len(keep_c)}/{len(sub)}")

    # ---- slice BEDs (merged loci + controls)
    for ref in ["grch38", "chm13"]:
        rows = [tuple(l.split("\t")[:3]) for l in open(os.path.join(OUT, f"{ref}.loci.bed"))]
        rows = [(c, int(s), int(e)) for c, s, e in rows]
        rows += [(r[0], int(r[1]), int(r[2])) for r in (l.split("\t") for l in open(os.path.join(OUT, f"{ref}.controls.bed")))]
        m = merge(rows)
        write_bed(os.path.join(OUT, f"{ref}.slice.bed"), m)
        log(f"{ref}.slice.bed: {len(m)} intervals, {sum(e - s for _, s, e in m) / 1e6:.2f} Mb")

    build_psv_tables()
    status("P5_regions", "DONE", "genes/loci/bins/controls/slice/psv for grch38 + chm13")


# ------------------------------------------------------------------ PSV tables
def smn_sites():
    rows = []
    for line in open(path("meta", "SMN_SNP_38.txt")):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        rows.append(dict(chrom=f[0], pos1=int(f[1]), b1=f[2], pos2=int(f[3]), b2=f[4], ann=f[5], sel=f[6]))
    return rows


def genes_bed(ref):
    out = defaultdict(list)
    for line in open(os.path.join(OUT, f"{ref}.genes.bed")):
        c, s, e, n, locus, strand = line.rstrip("\n").split("\t")
        out[n.split("#")[0]].append((c, int(s), int(e), n, strand))
    return out


def project_smn(ref, sites, q_start1, q_seq, fh):
    """Project GRCh38 SMN1 sites onto every SMN-like copy in `ref`'s SMN locus window."""
    locus = next(l.split("\t") for l in open(os.path.join(OUT, f"{ref}.loci.bed")) if l.split("\t")[3].strip() == "SMN")
    c, lo, hi = locus[0], int(locus[1]), int(locus[2])
    tmpfa = os.path.join(OUT, f"_{ref}_smn_locus.fa")
    with open(tmpfa, "w") as t:
        t.write(f">{c}_{lo}\n{read_fasta_region(REF_FA[ref], f'{c}:{lo + 1}-{hi}')[1]}\n")
    offs = {i: s["pos1"] - q_start1 for i, s in enumerate(sites)}
    hits = project_positions(q_seq, tmpfa, list(offs.values()))
    hits = [h for h in hits if h["qcov"] >= 0.5 and h["identity"] >= 0.9]
    hits.sort(key=lambda h: -(h["identity"] * (h["tend"] - h["tstart"])))
    kept = []
    for h in hits:  # drop hits overlapping a better hit by >50%
        if all(min(h["tend"], k["tend"]) - max(h["tstart"], k["tstart"]) < 0.5 * (h["tend"] - h["tstart"]) for k in kept):
            kept.append(h)
    i840 = next(i for i, s in enumerate(sites) if "splice" in s["ann"])
    out = []
    for k, h in enumerate(sorted(kept, key=lambda h: h["tstart"])):
        def base(i):
            p = h["sites"].get(offs[i])
            if p is None:
                return None, "-"
            b = p[1] if h["strand"] == "+" else revcomp(p[1])
            return lo + p[0] + 1, b
        _, b840 = base(i840)
        label = "SMN1-like" if b840 == sites[i840]["b1"] else "SMN2-like" if b840 == sites[i840]["b2"] else "SMN-other"
        copy = f"SMNcopy{k + 1}_{label}@{c}:{lo + h['tstart'] + 1}"
        for i, s in enumerate(sites):
            pos, b = base(i)
            if pos is None:
                continue
            site_id = f"SMN_{i + 1:02d}_{s['ann']}"
            fh.write(f"SMN\t{site_id}\t{copy}\t{c}\t{pos}\t{b}\t{s['b1']}\t{s['b2']}\t{s['sel']}\t{h['strand']}\n")
            out.append((copy, site_id, pos, b))
    return out


def hba_psvs(ref, fh):
    """HBA1 vs HBA2 single-base differences in `ref`, written for both copies."""
    G = genes_bed(ref)
    if not G.get("HBA1") or not G.get("HBA2"):
        log(f"{ref}: HBA1/HBA2 not both annotated; skipping HBA PSVs")
        return 0
    pad = 300
    c1, s1, e1, _, st1 = G["HBA1"][0]
    c2, s2, e2, _, st2 = G["HBA2"][0]
    tfa = os.path.join(OUT, f"_{ref}_hba1.fa")
    t_seq = read_fasta_region(REF_FA[ref], f"{c1}:{s1 - pad + 1}-{e1 + pad}")[1]
    with open(tfa, "w") as t:
        t.write(f">HBA1\n{t_seq}\n")
    q_seq = read_fasta_region(REF_FA[ref], f"{c2}:{s2 - pad + 1}-{e2 + pad}")[1]
    paf = run(["minimap2", "-x", "asm20", "-c", "--cs=long", tfa, "/dev/stdin"], input=f">HBA2\n{q_seq}\n")
    best = max((l.split("\t") for l in paf.splitlines()), key=lambda f: int(f[9]), default=None)
    if not best:
        log(f"{ref}: HBA2->HBA1 alignment failed")
        return 0
    cs = next(x[5:] for x in best[12:] if x.startswith("cs:Z:"))
    q, t, n = int(best[2]), int(best[7]), 0
    strand = best[4]
    if strand != "+":
        log(f"{ref}: HBA2 aligned to HBA1 on '-' strand (unexpected); skipping HBA PSVs — check annotation")
        return 0
    for op in re.findall(r"(=[A-Z]+|\*[a-z]{2}|\+[a-z]+|-[a-z]+)", cs):
        if op[0] == "=":
            q += len(op) - 1; t += len(op) - 1
        elif op[0] == "*":
            b_hba1, b_hba2q = op[1].upper(), op[2].upper()
            pos_hba1 = s1 - pad + t + 1
            pos_hba2 = (s2 - pad + q + 1) if strand == "+" else (s2 - pad + (len(q_seq) - 1 - q) + 1)
            n += 1
            sid = f"HBA_{n:02d}"
            fh.write(f"HBA\t{sid}\tHBA1\t{c1}\t{pos_hba1}\t{b_hba1}\t{b_hba1}\t{b_hba2q}\tyes\t+\n")
            fh.write(f"HBA\t{sid}\tHBA2\t{c2}\t{pos_hba2}\t{b_hba2q}\t{b_hba1}\t{b_hba2q}\tyes\t+\n")
            q += 1; t += 1
        elif op[0] == "+":
            q += len(op) - 1
        elif op[0] == "-":
            t += len(op) - 1
    log(f"{ref}: {n} HBA1/HBA2 single-base PSVs")
    return n


def build_psv_tables():
    sites = smn_sites()
    q_start1 = 70924941  # SMN1 exon1-8 span in GRCh38 (SMN_region_38.bed: 70924940-70954173, 0-based)
    q_seq = read_fasta_region(REF_FA["grch38"], f"chr5:{q_start1}-70954173")[1]
    header = "family\tsite_id\tcopy\tchrom\tpos\tbase\tbase_SMN1_or_HBA1\tbase_SMN2_or_HBA2\tselected\tstrand\n"

    # self-test: projecting onto GRCh38 must recover SMN1 and SMN2 coordinates from SMN_SNP_38.txt
    chk = os.path.join(OUT, "projection_check.tsv")
    with open(chk, "w") as fh:
        fh.write(header)
        got = project_smn("grch38", sites, q_start1, q_seq, fh)
    expected = {s["pos1"] for s in sites} | {s["pos2"] for s in sites}
    found = {p for _, _, p, _ in got}
    ok = len(expected & found)
    log(f"projection self-test on GRCh38: recovered {ok}/{len(expected)} SMN1/SMN2 site coordinates")
    if ok < 0.9 * len(expected):
        log("WARNING projection self-test below 90% — inspect data/regions/projection_check.tsv before trusting chm13.psv.tsv")

    for ref in ["grch38", "chm13"]:
        with open(os.path.join(OUT, f"{ref}.psv.tsv"), "w") as fh:
            fh.write(header)
            if ref == "grch38":
                for i, s in enumerate(sites):
                    sid = f"SMN_{i + 1:02d}_{s['ann']}"
                    fh.write(f"SMN\t{sid}\tSMN1\t{s['chrom']}\t{s['pos1']}\t{s['b1']}\t{s['b1']}\t{s['b2']}\t{s['sel']}\t+\n")
                    fh.write(f"SMN\t{sid}\tSMN2\t{s['chrom']}\t{s['pos2']}\t{s['b2']}\t{s['b1']}\t{s['b2']}\t{s['sel']}\t+\n")
            else:
                got = project_smn("chm13", sites, q_start1, q_seq, fh)
                log(f"chm13: {len({c for c, *_ in got})} SMN-like copies, {len(got)} projected sites")
            hba_psvs(ref, fh)


if __name__ == "__main__":
    main()
