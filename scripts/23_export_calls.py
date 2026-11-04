#!/usr/bin/env python3
"""S13 — every exact call of the benchmarked methods, as downloadable tables and one browsable page.

Inputs (all already computed): data/tables/calls_long.tsv (integer calls + truth, 13_concordance.py),
data/tables/depth.{ref}.tsv.gz (unrounded naive CN, 09), data/tables/dragen_calls.tsv (11),
data/tables/ctyper_calls.tsv (18), data/meta/…ped_population.txt.

Outputs (docs/ = GitHub Pages):
  docs/calls/calls_wide.tsv     one row per person: population, truth, then <method>.<target> integer calls
  docs/calls/calls_long.tsv     sample, method, target, call, truth (Bayes v0 rows excluded, as in the toolbox)
  docs/calls/naive_cn.tsv       unrounded naive CN per person × reference × MAPQ threshold × gene (+ control depth)
  docs/calls/dragen_calls.tsv   DRAGEN fields as parsed (raw totals, per-site SMN1 values, HBA genotype)
  docs/calls/ctyper_calls.tsv   ctyper gene counts and the chosen pangenome alleles
  explorer/calls.data.json      data for docs/calls.html (built by 19_build_site.py from explorer/calls.html)
"""
import csv
import gzip
import json
import os
from collections import defaultdict

from cfg import REPO_ROOT, path
from common import status

OUT = os.path.join(REPO_ROOT, "docs", "calls")
os.makedirs(OUT, exist_ok=True)
METHODS = ["dragen_4.0.3", "dragen_4.2.7", "dragen_4.4.7", "ctyper_grch38", "ctyper_chm13",
           "naive_grch38_mq0", "naive_chm13_mq0", "naive_grch38_mq20", "naive_chm13_mq20"]
TARGETS = ["SMN1", "SMN2_all", "SMN_total", "HBA_total", "CONTROL_RAD17"]
NAIVE_GENES = ["SMN1", "SMN2", "HBA1", "HBA2", "RAD17"]


def rows(fn):
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main():
    pop = {}
    for line in open(path("meta", "20130606_g1k_3202_samples_ped_population.txt")):
        f = line.split()
        if f and f[0] != "FamilyID":
            pop[f[1]] = (f[5], f[6])

    call = defaultdict(dict)    # sample -> {method.target: int}
    truth = defaultdict(dict)   # sample -> {target: int}
    long_rows = []
    for r in rows(path("tables", "calls_long.tsv")):
        if r["method"] not in METHODS:
            continue
        long_rows.append(r)
        if r["call"] != "":
            call[r["sample"]][f"{r['method']}.{r['target']}"] = int(r["call"])
        if r["truth"] != "":
            truth[r["sample"]][r["target"]] = int(r["truth"])
    samples = sorted(call)

    # unrounded naive CN (and the control median depth it was divided by)
    naive = defaultdict(dict)   # sample -> {ref.mq.gene: float}
    ctrl = {}
    for ref in ["grch38", "chm13"]:
        for r in rows(path("tables", f"depth.{ref}.tsv.gz")):
            if r["kind"] == "gene" and r["id"] in NAIVE_GENES and r["min_mapq"] in ("0", "20"):
                naive[r["sample"]][f"{ref}.{r['min_mapq']}.{r['id']}"] = (float(r["mean_depth"]), float(r["norm_cn"]))
    for r in rows(path("tables", "control_depth.tsv")):
        if r["min_mapq"] in ("0", "20"):
            ctrl[(r["sample"], r["ref"], r["min_mapq"])] = float(r["control_median_depth"])

    dragen = defaultdict(dict)
    for r in rows(path("tables", "dragen_calls.tsv")):
        dragen[r["sample"]][r["dragen_version"]] = r
    ctyper = defaultdict(dict)
    for r in rows(path("tables", "ctyper_calls.tsv")):
        ctyper[r["sample"]][r["ref"]] = r

    # ---- TSVs
    cols = [f"{m}.{t}" for m in METHODS for t in TARGETS if m.startswith("naive") or not t.startswith("CONTROL_")]  # only naive depth reports RAD17
    with open(os.path.join(OUT, "calls_wide.tsv"), "w") as fh:
        fh.write("\t".join(["sample", "population", "superpopulation"] + [f"truth.{t}" for t in TARGETS] + cols) + "\n")
        for s in samples:
            p = pop.get(s, ("", ""))
            fh.write("\t".join([s, *p] + [str(truth[s].get(t, "")) for t in TARGETS] + [str(call[s].get(c, "")) for c in cols]) + "\n")
    with open(os.path.join(OUT, "calls_long.tsv"), "w") as fh:
        fh.write("sample\tmethod\ttarget\tcall\ttruth\n")
        for r in long_rows:
            fh.write("\t".join(r[k] for k in ("sample", "method", "target", "call", "truth")) + "\n")
    with open(os.path.join(OUT, "naive_cn.tsv"), "w") as fh:
        fh.write("sample\tref\tmin_mapq\tcontrol_median_depth\tgene\tmean_depth\tcn_unrounded\tcall\n")
        for s in samples:
            for ref in ["grch38", "chm13"]:
                for q in ("0", "20"):
                    for g in NAIVE_GENES:
                        v = naive[s].get(f"{ref}.{q}.{g}")
                        if v:
                            fh.write(f"{s}\t{ref}\t{q}\t{ctrl.get((s, ref, q), '')}\t{g}\t{v[0]}\t{v[1]}\t{round(v[1])}\n")
    for name in ("dragen_calls.tsv", "ctyper_calls.tsv"):
        src = list(rows(path("tables", name)))
        keep = [k for k in src[0] if k != "source_files"]
        with open(os.path.join(OUT, name), "w") as fh:
            fh.write("\t".join(keep) + "\n")
            for r in src:
                fh.write("\t".join(r[k] for k in keep) + "\n")

    # ---- page data (compact)
    page = dict(methods=METHODS, targets=TARGETS, people=[])
    for s in samples:
        p = pop.get(s, ("", ""))
        d = {"id": s, "pop": p[0], "sp": p[1], "t": [truth[s].get(t) for t in TARGETS],
             "c": [[call[s].get(f"{m}.{t}") for t in TARGETS] for m in METHODS]}
        # hover details: unrounded naive CN, DRAGEN raw totals/genotype, ctyper raw fields
        nv = {}
        for ref in ["grch38", "chm13"]:
            for q in ("0", "20"):
                nv[f"{ref}_{q}"] = {g: round(naive[s][f"{ref}.{q}.{g}"][1], 2) for g in NAIVE_GENES if f"{ref}.{q}.{g}" in naive[s]}
        d["naive"] = nv
        d["dragen"] = {v: {k: r[k] for k in ("SMN1_CN", "SMN2_CN", "SMN2_del78_CN", "SMN_total_CN_raw", "SMN_fulllength_CN_raw",
                                             "SMN1_CN_raw_sites", "HBA_genotype")} for v, r in dragen[s].items()}
        d["ctyper"] = {ref: {k: r[k] for k in ("SMN1", "SMN2", "SMN_partial", "HBA1", "HBA2")} for ref, r in ctyper[s].items()}
        page["people"].append(d)
    with open(os.path.join(REPO_ROOT, "explorer", "calls.data.json"), "w") as fh:
        json.dump(page, fh, separators=(",", ":"))
    status("S13_calls", "DONE", f"docs/calls/*.tsv: {len(samples)} people × {len(METHODS)} methods × {len(TARGETS)} targets")


if __name__ == "__main__":
    main()
