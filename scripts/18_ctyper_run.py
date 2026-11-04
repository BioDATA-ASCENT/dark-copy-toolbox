#!/usr/bin/env python3
"""S8b — run ctyper (k-mer genotyping against pangenome alleles; no read mapping) on every region slice, then turn
the called pangenome alleles into gene copy numbers.

Per sample × reference:
  slice CRAM -> temporary BAM (ctyper reads every record in the file; it does not need the reference coordinates,
  so GRCh38 and CHM13 slices are two independent read sets for the same k-mer matrix)
  ctyper -i tmp.bam -m data/ctyper/subset_matrix.txt -d <haploid 31-mer depth>
     haploid 31-mer depth = median control-window read depth (MAPQ>=0, from 09) × (L-30)/L ÷ 2
     (calibrated on 5 truth samples: with the diploid value every copy number came out halved)
Allele -> genes: transcript column of PangenomeAlleles_annotationfix (one HBA allele = the HBA2+HBA1 block, so
2 genes normally, 1 for -alpha3.7/-alpha4.2 hybrids, 3 for triplications; one SMN allele = SMN1 or SMN2).

Caveat: slices hold only reads aligned inside our windows (+ their unmapped mates); reads aligned to GRCh38 alt
contigs (~5% at SMN, checked on HG00097) or elsewhere are missing, so this is "ctyper on a slice", not a full run.

Outputs: data/ctyper/out/{sample}.{ref}.out(.ok), data/tables/ctyper_calls.tsv
"""
import csv
import gzip
import os
import re
import statistics
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from cfg import CFG, path
from common import log, status

C = CFG["ctyper"]
OUTD = path("ctyper", "out")
TMP = path("ctyper", "tmp")
os.makedirs(OUTD, exist_ok=True)
os.makedirs(TMP, exist_ok=True)
MATRIX = path("ctyper", "subset_matrix.txt")
L = C["read_length"]


def control_depth(ref):
    d = defaultdict(list)
    with gzip.open(path("tables", f"depth.{ref}.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] == "control" and r["min_mapq"] == "0":
                d[r["sample"]].append(float(r["mean_depth"]))
    return {s: statistics.median(v) for s, v in d.items() if len(v) > 10}


def run_one(sample, ref, depth):
    out = os.path.join(OUTD, f"{sample}.{ref}.out")
    if os.path.exists(out + ".ok"):
        return out, "cached"
    cram = path("slices", ref, f"{sample}.{ref}.cram")
    fa = path(CFG["ref_fasta"][ref])
    kdep = round(depth * (L - 30) / L / 2, 3)
    with tempfile.TemporaryDirectory(dir=TMP) as td:
        bam = os.path.join(td, "s.bam")
        subprocess.run(["samtools", "view", "-b", "-T", fa, "-o", bam, cram], check=True)
        subprocess.run(["samtools", "index", bam], check=True)
        r = subprocess.run(["conda", "run", "-n", C["conda_env"], "ctyper", "-i", bam, "-m", MATRIX, "-o", out,
                            "-d", str(kdep), "-n", str(C["threads_per_job"])], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out):
        return out, f"FAIL {r.stderr[-300:]}"
    with open(out + ".ok", "w") as fh:
        fh.write(f"depth={depth} kmer_depth_haploid={kdep}\n")
    return out, "ok"


def allele_genes():
    """allele name -> list of genes it carries (from its transcript column), e.g. ['HBA1','HBA2']."""
    g = {}
    for line in open(path("ctyper", "annotation_subset.tsv")):
        f = line.rstrip("\n").split("\t")
        g[f[0]] = [m.group(2) for m in re.finditer(r"ENST\d+(\([^)]*\))?:([A-Za-z0-9-]+):[\d.]+", f[4])]
        if any("(" in (m.group(1) or "") for m in re.finditer(r"ENST\d+(\([^)]*\))?:", f[4])):
            g[f[0]] = [x + ("_partial" if x.startswith("SMN") else "") for x in g[f[0]]]
    return g


def parse(out, genes_of):
    res, cur = {}, None
    for line in open(out):
        if line.startswith(">"):
            cur = line[2:].split()[0]
        elif line.startswith("result:") and cur:
            res[cur] = [a for a in line.split(":", 1)[1].strip().strip(",").split(",") if a]
    cn = Counter()
    for grp, alleles in res.items():
        fam = grp.split("_")[0]
        cn[f"{fam}_alleles"] += len(alleles)
        for a in alleles:
            for gname in genes_of.get(a, []):
                cn[gname] += 1
    return cn, res


def main():
    if not os.path.exists(MATRIX + ".ok"):
        raise SystemExit("run scripts/17_ctyper_subset.py first")
    sel = [r["sample"] for r in csv.DictReader(open(path("meta", "samples_selected.tsv")), delimiter="\t")]
    jobs = []
    for ref in ["grch38", "chm13"]:
        dep = control_depth(ref)
        for s in sel:
            if s in dep and os.path.exists(path("slices", ref, f"{s}.{ref}.cram.ok")):
                jobs.append((s, ref, dep[s]))
    log(f"ctyper: {len(jobs)} sample×ref jobs")
    n_fail = 0
    with ThreadPoolExecutor(max_workers=C["jobs"]) as ex:
        futs = {ex.submit(run_one, *j): j for j in jobs}
        for i, f in enumerate(as_completed(futs), 1):
            out, st = f.result()
            if st.startswith("FAIL"):
                n_fail += 1
                log(f"ctyper {futs[f][:2]}: {st}")
            if i % 100 == 0:
                log(f"ctyper: {i}/{len(jobs)}")
    genes_of = allele_genes()
    cols = ["sample", "ref", "SMN1", "SMN2", "SMN_partial", "SMN_alleles", "HBA1", "HBA2", "HBA_total", "HBA_alleles",
            "CBS", "U2AF1", "alleles"]
    with open(path("tables", "ctyper_calls.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s, ref, _ in jobs:
            out = os.path.join(OUTD, f"{s}.{ref}.out")
            if not os.path.exists(out + ".ok"):
                continue
            cn, res = parse(out, genes_of)
            row = dict(sample=s, ref=ref, SMN1=cn["SMN1"], SMN2=cn["SMN2"],
                       SMN_partial=cn["SMN1_partial"] + cn["SMN2_partial"], SMN_alleles=cn["SMN_alleles"],
                       HBA1=cn["HBA1"], HBA2=cn["HBA2"], HBA_total=cn["HBA1"] + cn["HBA2"], HBA_alleles=cn["HBA_alleles"],
                       CBS=cn["CBS"], U2AF1=cn["U2AF1"],
                       alleles=";".join(f"{g}={','.join(a)}" for g, a in sorted(res.items())))
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    status("S8_ctyper", "DONE", f"{len(jobs) - n_fail}/{len(jobs)} sample×ref runs; data/tables/ctyper_calls.tsv")


if __name__ == "__main__":
    main()
