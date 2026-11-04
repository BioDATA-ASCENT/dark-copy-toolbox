#!/usr/bin/env python3
"""Stretch 2 — Locityper tidy tables (data/tables/locityper_*.tsv).

  locityper_our_loci.tsv   Locityper db loci whose GRCh38 ref.bed overlaps any of our locus windows
  locityper_loo.tsv        leave-one-out benchmark rows (illumina) for those loci and our truth samples
  locityper_1kgp.tsv       1kGP genotype (haplotype pair) for every selected sample at those loci;
                           `circular` = yes when the sample's own haplotypes are in the Locityper panel
                           (non-LOO calls for panel samples are not an independent test)
  locityper_cn.tsv         Locityper-implied gene copy numbers: the two predicted haplotypes are typed with
                           type_copies() from 08 for every family whose query gene lies inside the locus.

NOTE (verified 2026-11-04): the Locityper v0.17 db (326 loci) has NO locus covering SMN1/SMN2 or HBA1/HBA2.
HBM (chr16:152,245-166,764) stops 6 kb before HBA2; CBS covers the chr21 false-duplication demo.
"""
import csv
import glob
import gzip
import importlib.util
import os

from cfg import path
from common import log, status

L = path("precomputed", "locityper")
DB = next((d for d in [os.path.join(L, "db", "db", "loci"), os.path.join(L, "db", "loci")] if os.path.isdir(d)), None)
TB = path("tables")


def overlaps(a, b):
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


def read_csv_gz(fn):
    with gzip.open(fn, "rt") as fh:
        rows = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def main():
    windows = [(l.split("\t")[0], int(l.split("\t")[1]), int(l.split("\t")[2]), l.split("\t")[3].strip())
               for l in open(path("regions", "grch38.loci.bed"))]
    genes = [(l.split("\t")[0], int(l.split("\t")[1]), int(l.split("\t")[2]), l.split("\t")[3])
             for l in open(path("regions", "grch38.genes.bed"))]
    ours = []
    for d in sorted(glob.glob(os.path.join(DB, "*"))):
        f = open(os.path.join(d, "ref.bed")).readline().split("\t")
        iv = (f[0], int(f[1]), int(f[2]))
        hit = [w[3] for w in windows if overlaps(iv, w)]
        if hit:
            g = [x[3] for x in genes if overlaps(iv, x)]
            ours.append(dict(locus=os.path.basename(d), chrom=iv[0], start=iv[1], end=iv[2],
                             our_window=",".join(sorted(set(hit))), genes_inside=",".join(g)))
    with open(os.path.join(TB, "locityper_our_loci.tsv"), "w") as fh:
        cols = ["locus", "chrom", "start", "end", "our_window", "genes_inside"]
        fh.write("\t".join(cols) + "\n")
        for r in ours:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    loci = {r["locus"] for r in ours}
    log(f"Locityper loci overlapping our windows: {sorted(loci)}")

    sel = list(csv.DictReader(open(path("meta", "samples_selected.tsv")), delimiter="\t"))
    truth = {r["sample"] for r in sel if r["role"] == "truth"}
    selset = {r["sample"] for r in sel}

    bench = glob.glob(os.path.join(L, "locityper-benchmarking", "**", "QV", "loo_illumina.csv.gz"), recursive=True)
    loo = [r for r in read_csv_gz(bench[0]) if r["locus"] in loci] if bench else []
    with open(os.path.join(TB, "locityper_loo.tsv"), "w") as fh:
        if loo:
            cols = ["in_our_truth"] + list(loo[0].keys())
            fh.write("\t".join(cols) + "\n")
            for r in loo:
                r["in_our_truth"] = "yes" if r["sample"] in truth else "no"
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    log(f"LOO rows: {len(loo)} ({sum(r['sample'] in truth for r in loo)} for our truth samples)")

    panel = {}
    for lc in loci:
        with gzip.open(os.path.join(DB, lc, "haplotypes.fa.gz"), "rt") as fh:
            panel[lc] = {l[1:].split()[0] for l in fh if l.startswith(">")}
    gts = [r for r in read_csv_gz(os.path.join(L, "1kgp_genotypes.csv.gz")) if r["locus"] in loci and r["sample"] in selset]
    with open(os.path.join(TB, "locityper_1kgp.tsv"), "w") as fh:
        cols = ["sample", "locus", "genotype", "quality", "total_reads", "unexpl_reads", "weight_dist", "warnings", "circular"]
        fh.write("\t".join(cols) + "\n")
        for r in gts:
            r["circular"] = "yes" if any(h.split(".")[0] == r["sample"] for h in panel[r["locus"]]) else "no"
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    log(f"1kGP genotypes for selected samples: {len(gts)}")

    # Locityper-implied CN: only families whose GRCh38 query gene (SMN1 / HBA2) lies inside a Locityper locus
    spec = importlib.util.spec_from_file_location("truth08", os.path.join(os.path.dirname(__file__), "08_truth_from_assemblies.py"))
    t8 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t8)
    fam_gene = {"SMN": "SMN1", "HBA": "HBA2"}
    typable = {lc: [fam for fam, g in fam_gene.items() if g in r["genes_inside"].split(",")] for lc, r in
               ((r["locus"], r) for r in ours)}
    n_typed = 0
    with open(os.path.join(TB, "locityper_cn.tsv"), "w") as fh:
        fh.write("sample\tlocus\tfamily\thap1\thap2\tcalls_hap1\tcalls_hap2\tcircular\tnote\n")
        if not any(typable.values()):
            fh.write("#\tNA\tSMN,HBA\t\t\t\t\t\tno Locityper locus contains SMN1 or HBA2: Locityper-implied CN not available\n")
        else:
            queries = t8.load_queries()
            seqs = {}
            for lc, fams in typable.items():
                if not fams:
                    continue
                with gzip.open(os.path.join(DB, lc, "haplotypes.fa.gz"), "rt") as h:
                    name = None
                    for line in h:
                        if line.startswith(">"):
                            name = line[1:].split()[0]
                            seqs[(lc, name)] = []
                        else:
                            seqs[(lc, name)].append(line.strip())
            for r in gts:
                for fam in typable.get(r["locus"], []):
                    haps = r["genotype"].split(",")
                    calls = []
                    for h in haps:
                        s = "".join(seqs.get((r["locus"], h), []))
                        calls.append("+".join(c["call"] for c in t8.type_copies(fam, s.upper(), queries[fam])) if s else "NA")
                    fh.write(f"{r['sample']}\t{r['locus']}\t{fam}\t{haps[0]}\t{haps[-1]}\t{calls[0]}\t{calls[-1]}\t{r['circular']}\t\n")
                    n_typed += 1
    status("S2_locityper", "DONE", f"loci {sorted(loci)}; LOO rows {len(loo)}; 1kGP gts {len(gts)}; CN typed {n_typed} (no SMN/HBA locus in db)")


if __name__ == "__main__":
    main()
