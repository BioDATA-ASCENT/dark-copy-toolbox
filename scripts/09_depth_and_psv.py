#!/usr/bin/env python3
"""Phase 9 — per-slice depth and PSV allele counts (inputs for tomorrow's callers).

  data/tables/depth.{ref}.tsv.gz       sample, ref, id, kind(bin|gene|control), chrom, start, end, min_mapq,
                                       mean_depth, norm_cn (= 2 * mean / median control mean at that MAPQ)
  data/tables/psv_counts.{ref}.tsv.gz  sample, ref, family, site_id, copy, chrom, pos, fwd_base_A, fwd_base_B,
                                       nA, nC, nG, nT, n_baseA, n_baseB, depth   (A = SMN1/HBA1, B = SMN2/HBA2)
  data/tables/control_depth.tsv        per sample/ref/min_mapq control median (sanity: ~30x)

MAPQ thresholds: 0 (everything), 1, 20 — the 0 vs 20 gap is the "multi-mapping loss" signal.
"""
import csv
import gzip
import os
import statistics
import subprocess
import sys
from multiprocessing import Pool

import pysam

from cfg import CFG, path
from common import log, revcomp, status

REF_FA = {k: path(v) for k, v in CFG["ref_fasta"].items()}
G, S, TB = path("regions"), path("slices"), path("tables")
os.makedirs(TB, exist_ok=True)
MAPQS = [0, 1, 20]
EXT = CFG["slice_format"]


def region_bed(ref):
    fn = os.path.join(G, f"{ref}.depth_regions.bed")
    rows = []
    for kind, f in [("bin", "bins"), ("gene", "genes"), ("control", "controls")]:
        for line in open(os.path.join(G, f"{ref}.{f}.bed")):
            x = line.rstrip("\n").split("\t")
            rows.append((x[0], int(x[1]), int(x[2]), x[3], kind))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(fn, "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return fn


def depth_one(args):
    """Depth at 3 MAPQ thresholds + PSV counts for one slice.

    The CRAM slice is first decoded once into a temporary uncompressed BAM: samtools bedcov / pysam fetches
    on CRAM re-decode a whole container per region (~3 min per slice for ~4k regions), on BAM it is ~1 s.
    """
    ref, sample, slice_fn, bed, sites = args
    import tempfile
    tmpd = tempfile.mkdtemp(dir=TB, prefix=".tmp_")
    bam = os.path.join(tmpd, "s.bam")
    try:
        p = subprocess.run(["samtools", "view", "-u", "-T", REF_FA[ref], "-o", bam, slice_fn], capture_output=True, text=True)
        if p.returncode == 0:
            p = subprocess.run(["samtools", "index", bam], capture_output=True, text=True)
        if p.returncode:
            return sample, None, None, p.stderr[-500:]
        out = []
        for q in MAPQS:
            p = subprocess.run(["samtools", "bedcov", "-Q", str(q), bed, bam], capture_output=True, text=True)
            if p.returncode:
                return sample, None, None, p.stderr[-500:]
            rows = [l.split("\t") for l in p.stdout.splitlines()]
            means = [(r[0], int(r[1]), int(r[2]), r[3], r[4], int(r[5]) / max(1, int(r[2]) - int(r[1]))) for r in rows]
            ctrl = [m[5] for m in means if m[4] == "control"]
            med = statistics.median(ctrl) if ctrl else float("nan")
            for c, s, e, rid, kind, mean in means:
                out.append((sample, ref, rid, kind, c, s, e, q, round(mean, 3), round(2 * mean / med, 4) if med else "nan"))
            out.append(("__ctrl__", sample, ref, q, round(med, 3)))
        psv = psv_one((ref, sample, bam, sites))
        return sample, out, psv, ""
    finally:
        for f in os.listdir(tmpd):
            os.unlink(os.path.join(tmpd, f))
        os.rmdir(tmpd)


def psv_one(args):
    ref, sample, slice_fn, sites = args
    out = []
    with pysam.AlignmentFile(slice_fn, reference_filename=REF_FA[ref]) as af:
        for r in sites:
            pos0 = int(r["pos"]) - 1
            cov = af.count_coverage(r["chrom"], pos0, pos0 + 1, quality_threshold=13, read_callback="all")
            nA, nC, nG, nT = (int(c[0]) for c in cov)
            a, b = r["base_SMN1_or_HBA1"], r["base_SMN2_or_HBA2"]
            if r["strand"] == "-":
                a, b = revcomp(a), revcomp(b)
            idx = {"A": nA, "C": nC, "G": nG, "T": nT}
            out.append((sample, ref, r["family"], r["site_id"], r["copy"], r["chrom"], r["pos"], a, b,
                        nA, nC, nG, nT, idx.get(a, 0), idx.get(b, 0), nA + nC + nG + nT))
    return out


def main():
    refs = [r for r in ["grch38", "chm13"] if os.path.isdir(os.path.join(S, r))]
    ctrl_rows = []
    for ref in refs:
        bed = region_bed(ref)
        slices = sorted(f for f in os.listdir(os.path.join(S, ref)) if f.endswith(f".{EXT}") and os.path.exists(os.path.join(S, ref, f + ".ok")))
        sites = list(csv.DictReader(open(os.path.join(G, f"{ref}.psv.tsv")), delimiter="\t"))
        jobs = [(ref, f.split(".")[0], os.path.join(S, ref, f), bed, sites) for f in slices]
        log(f"P9 {ref}: {len(jobs)} slices")
        with Pool(CFG["threads"]) as pool, gzip.open(os.path.join(TB, f"depth.{ref}.tsv.gz"), "wt") as fh, \
                gzip.open(os.path.join(TB, f"psv_counts.{ref}.tsv.gz"), "wt") as fp:
            fh.write("sample\tref\tid\tkind\tchrom\tstart\tend\tmin_mapq\tmean_depth\tnorm_cn\n")
            fp.write("sample\tref\tfamily\tsite_id\tcopy\tchrom\tpos\tfwd_base_A\tfwd_base_B\tnA\tnC\tnG\tnT\tn_baseA\tn_baseB\tdepth\n")
            for i, (sample, rows, psv, err) in enumerate(pool.imap_unordered(depth_one, jobs), 1):
                if rows is None:
                    log(f"depth failed {sample} {ref}: {err}")
                    continue
                for r in rows:
                    if r[0] == "__ctrl__":
                        ctrl_rows.append(r[1:])
                    else:
                        fh.write("\t".join(map(str, r)) + "\n")
                for r in psv:
                    fp.write("\t".join(map(str, r)) + "\n")
                if i % 50 == 0:
                    log(f"P9 {ref}: {i}/{len(jobs)}")
    with open(os.path.join(TB, "control_depth.tsv"), "w") as fh:
        fh.write("sample\tref\tmin_mapq\tcontrol_median_depth\n")
        for r in sorted(ctrl_rows):
            fh.write("\t".join(map(str, r)) + "\n")
    meds = [r[3] for r in ctrl_rows if r[2] == 0]
    status("P9_tables", "DONE", f"{len(meds)} sample×ref depth tables; median control depth {statistics.median(meds) if meds else 'NA'}x")


if __name__ == "__main__":
    sys.exit(main())
