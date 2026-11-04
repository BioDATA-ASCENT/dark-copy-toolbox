#!/usr/bin/env python3
"""S11 — depth at a ladder of MAPQ thresholds (config mapq_ladder: 0, 5, …, 60) for the toolbox's MAPQ slider.

Same measurement as 09 (samtools bedcov -Q q over bins, genes and control windows of a slice decoded once to a
temporary BAM), at every threshold; normalised CN = 2 × mean depth / median control depth *at the same threshold*.
The 0/1/20 tables from 09 are not touched.

Output: data/tables/mapq_ladder.{ref}.npz  (samples, qs, region ids/kinds, float32 CN array [sample, q, region])
"""
import os
import statistics
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from cfg import CFG, path
from common import log, status

QS = list(CFG.get("mapq_ladder", list(range(0, 61, 5))))
REF_FA = {k: path(v) for k, v in CFG["ref_fasta"].items()}
TB = path("tables")


def one(args):
    ref, sample, slice_fn, bed = args
    tmpd = tempfile.mkdtemp(dir=TB, prefix=".tmpq_")
    bam = os.path.join(tmpd, "s.bam")
    try:
        p = subprocess.run(["samtools", "view", "-u", "-T", REF_FA[ref], "-o", bam, slice_fn], capture_output=True, text=True)
        if p.returncode == 0:
            p = subprocess.run(["samtools", "index", bam], capture_output=True, text=True)
        if p.returncode:
            return sample, None, p.stderr[-300:]
        mat = []
        for q in QS:
            p = subprocess.run(["samtools", "bedcov", "-Q", str(q), bed, bam], capture_output=True, text=True)
            if p.returncode:
                return sample, None, p.stderr[-300:]
            rows = [l.split("\t") for l in p.stdout.splitlines()]
            mean = np.array([int(r[5]) / max(1, int(r[2]) - int(r[1])) for r in rows], dtype=np.float64)
            kinds = [r[4] for r in rows]
            med = statistics.median([m for m, k in zip(mean, kinds) if k == "control"])
            mat.append((2 * mean / med) if med else np.full(len(mean), np.nan))
        return sample, np.array(mat, dtype=np.float32), ""
    finally:
        for f in os.listdir(tmpd):
            os.unlink(os.path.join(tmpd, f))
        os.rmdir(tmpd)


def main():
    for ref in ["grch38", "chm13"]:
        out = path("tables", f"mapq_ladder.{ref}.npz")
        if os.path.exists(out):
            continue
        bed = path("regions", f"{ref}.depth_regions.bed")          # written by 09: chrom start end id kind
        regs = [l.rstrip("\n").split("\t") for l in open(bed)]
        sd = path("slices", ref)
        jobs = [(ref, f.split(".")[0], os.path.join(sd, f), bed) for f in sorted(os.listdir(sd))
                if f.endswith(f".{ref}.{CFG['slice_format']}") and os.path.exists(os.path.join(sd, f + ".ok"))]
        log(f"mapq ladder {ref}: {len(jobs)} slices × {len(QS)} thresholds")
        res = {}
        with ProcessPoolExecutor(max_workers=CFG.get("threads", 8) * 2) as ex:
            futs = [ex.submit(one, j) for j in jobs]
            for i, f in enumerate(as_completed(futs), 1):
                s, m, err = f.result()
                if m is None:
                    log(f"mapq ladder {ref} {s}: {err}")
                else:
                    res[s] = m
                if i % 100 == 0:
                    log(f"mapq ladder {ref}: {i}/{len(jobs)}")
        samples = sorted(res)
        np.savez_compressed(out, samples=np.array(samples), qs=np.array(QS), ids=np.array([r[3] for r in regs]),
                            kinds=np.array([r[4] for r in regs]), cn=np.stack([res[s] for s in samples]))
        status(f"S11_mapq_ladder_{ref}", "DONE", f"{len(samples)} slices × {len(QS)} MAPQ thresholds")


if __name__ == "__main__":
    main()
