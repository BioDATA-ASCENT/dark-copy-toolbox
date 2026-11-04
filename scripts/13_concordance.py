#!/usr/bin/env python3
"""Stretch 3 — every method vs assembly truth, on the same targets.

Targets (diploid integers):  SMN1 · SMN2_all (= full-length SMN2 + SMN2Δ7-8) · SMN_total · HBA_total
Methods:
  naive_{ref}_mq{0,20}  round(norm_cn) per gene from data/tables/depth.{ref}.tsv.gz (no calibration at all)
  dragen_{4.0.3,4.2.7,4.4.7}   data/tables/dragen_calls.tsv (SMN2_all = SMN2_CN + SMN2_del78_CN; 4.0.3 has no HBA)
  bayes_{tag}           caller/bayes_cn.py outputs in data/tables/bayes_cn.*.tsv
  ctyper_{ref}          ctyper 1.0.5 run on our slices (scripts/18_ctyper_run.py)
  (side table) ctyper annotation copies per assembly in ctyper's PangenomeAlleles table (assembly-based, so a check of
                        OUR truth typing, not a read-based caller); HPRC-R1-era assemblies, joined by sample

Outputs: data/tables/calls_long.tsv (sample, method, target, call, truth)
         data/tables/concordance.tsv (method, target, n, n_correct, accuracy, n_over, n_under)
         data/tables/truth_crosscheck.tsv (every truth vs DRAGEN 4.2.7 disagreement)
"""
import csv
import glob
import gzip
import os
import re
from collections import defaultdict

from cfg import CFG, path
from common import log, status

TB = path("tables")


def tsv(fn):
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ival(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def main():
    truth = {}
    for r in tsv(path("truth", "diploid.tsv")):
        t = {}
        if r["usable_SMN"] == "True":
            t["SMN1"] = int(r["SMN1"])
            t["SMN2_all"] = int(r["SMN2"]) + int(r["SMN_del7-8"])
            t["SMN_total"] = t["SMN1"] + t["SMN2_all"] + int(r["SMN?"])
        if r["usable_HBA"] == "True":
            t["HBA_total"] = int(r["HBA_total"])
        truth[r["sample"]] = t
    # positive control: copies of a unique gene annotated (CAT) on each haplotype, summed when both haplotypes exist
    CTRL = CFG.get("control_gene")
    if CTRL:
        per = defaultdict(list)
        for r in tsv(path("truth", "haplotypes.tsv")):
            fn = path("truth", "cat", f"{r['assembly']}.genes.tsv")
            if os.path.exists(fn):
                per[r["sample"]].append(sum(l.split("\t")[4] == CTRL for l in open(fn)))
        for smp, v in per.items():
            if len(v) == 2:
                truth.setdefault(smp, {})[f"CONTROL_{CTRL}"] = sum(v)
    calls = defaultdict(dict)  # (sample, method) -> {target: call}

    for ref in ["grch38", "chm13"]:
        fn = path("tables", f"depth.{ref}.tsv.gz")
        if not os.path.exists(fn):
            continue
        g = defaultdict(dict)
        for r in tsv(fn):
            if r["kind"] == "gene" and r["id"] in ("SMN1", "SMN2", "HBA1", "HBA2", CTRL) and r["min_mapq"] in ("0", "20"):
                g[(r["sample"], r["min_mapq"])][r["id"]] = float(r["norm_cn"])
        for (s, q), d in g.items():
            c = {k: round(v) for k, v in d.items()}
            out = {}
            if "SMN1" in c and "SMN2" in c:
                out.update(SMN1=c["SMN1"], SMN2_all=c["SMN2"], SMN_total=c["SMN1"] + c["SMN2"])
            if "HBA1" in c and "HBA2" in c:
                out["HBA_total"] = c["HBA1"] + c["HBA2"]
            if CTRL and CTRL in c:
                out[f"CONTROL_{CTRL}"] = c[CTRL]
            calls[(s, f"naive_{ref}_mq{q}")] = out

    if os.path.exists(path("tables", "dragen_calls.tsv")):
        for r in tsv(path("tables", "dragen_calls.tsv")):
            a, b, d = ival(r["SMN1_CN"]), ival(r["SMN2_CN"]), ival(r["SMN2_del78_CN"])
            out = {}
            if a is not None and b is not None:
                out.update(SMN1=a, SMN2_all=b + (d or 0), SMN_total=a + b + (d or 0))
            if ival(r["HBA_total_copies"]) is not None:
                out["HBA_total"] = ival(r["HBA_total_copies"])
            calls[(r["sample"], f"dragen_{r['dragen_version']}")] = out

    for fn in glob.glob(path("tables", "bayes_cn.*.tsv")):
        tag = os.path.basename(fn)[len("bayes_cn."):-4]
        for r in tsv(fn):
            m = calls[(r["sample"], f"bayes_{tag}")]
            if r["family"] == "SMN":
                m.update(SMN1=int(r["map_A"]), SMN2_all=int(r["map_B"]), SMN_total=int(r["map_total"]))
            else:
                m["HBA_total"] = int(r["map_total"])

    # ctyper run on our slices (scripts/18): k-mer genotyping against pangenome alleles
    if os.path.exists(path("tables", "ctyper_calls.tsv")):
        for r in tsv(path("tables", "ctyper_calls.tsv")):
            calls[(r["sample"], f"ctyper_{r['ref']}")] = dict(
                SMN1=int(r["SMN1"]), SMN2_all=int(r["SMN2"]) + int(r["SMN_partial"]), SMN_total=int(r["SMN_alleles"]),
                HBA_total=int(r["HBA_total"]))

    # ctyper pangenome allele annotation: copies per assembly haplotype
    cty = path("precomputed", "ctyper", "PangenomeAlleles_typefix.v1.0.tsv")
    cty_hap = defaultdict(lambda: defaultdict(int))
    if os.path.exists(cty):
        for line in open(cty):
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or not re.match(r"(SMN|HBA)", f[2]):
                continue
            gene = f[2].split(":")[0]
            for mem in f[9].split(","):
                m = re.search(r"_((?:HG|NA)\d{5})_h([12])_", mem)
                if m:
                    cty_hap[(m.group(1), m.group(2))][gene] += 1
        cty_calls = {}
        per = defaultdict(lambda: defaultdict(int))
        nh = defaultdict(set)
        for (s, h), d in cty_hap.items():
            nh[s].add(h)
            for k, v in d.items():
                per[s][k] += v
        for s, d in per.items():
            if len(nh[s]) != 2:
                continue
            out = {"SMN1": d.get("SMN1", 0), "SMN2_all": d.get("SMN2", 0)}
            out["SMN_total"] = out["SMN1"] + out["SMN2_all"]
            out["HBA_total"] = d.get("HBA1", 0) + d.get("HBA2", 0)
            cty_calls[s] = out
        # NOT scored: ctyper's HBA alleles are one bundled HBA2,HBA1,HBAP1 group (no per-gene count) and its SMN
        # alleles cover only the R1 haplotypes that resolved the locus -> 50% "accuracy" is an annotation-coverage
        # artefact, not a caller result. Kept as a side table for reference only.
        with open(os.path.join(TB, "ctyper_annotation_counts.tsv"), "w") as fh:
            fh.write("sample\tSMN1\tSMN2_all\tHBA_group_alleles\n")
            for s, o in sorted(cty_calls.items()):
                fh.write(f"{s}\t{o['SMN1']}\t{o['SMN2_all']}\t{o['HBA_total']}\n")

    long, conc = [], defaultdict(lambda: [0, 0, 0, 0])
    for (s, m), d in calls.items():
        for tgt, c in d.items():
            t = truth.get(s, {}).get(tgt)
            long.append((s, m, tgt, c, "" if t is None else t))
            if t is not None and c is not None:
                k = conc[(m, tgt)]
                k[0] += 1
                k[1] += c == t
                k[2] += c > t
                k[3] += c < t
    with open(os.path.join(TB, "calls_long.tsv"), "w") as fh:
        fh.write("sample\tmethod\ttarget\tcall\ttruth\n")
        for r in sorted(long):
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(TB, "concordance.tsv"), "w") as fh:
        fh.write("method\ttarget\tn\tn_correct\taccuracy\tn_over\tn_under\n")
        for (m, tgt), (n, ok, ov, un) in sorted(conc.items()):
            fh.write(f"{m}\t{tgt}\t{n}\t{ok}\t{ok / n:.3f}\t{ov}\t{un}\n")
    dis = [r for r in long if r[1] == "dragen_4.2.7" and r[4] != "" and r[3] != r[4]]
    with open(os.path.join(TB, "truth_crosscheck.tsv"), "w") as fh:
        fh.write("sample\tmethod\ttarget\tcall\ttruth\n")
        for r in sorted(dis):
            fh.write("\t".join(map(str, r)) + "\n")
    log(f"concordance: {len(conc)} method×target cells; {len(dis)} truth disagreements vs DRAGEN 4.2.7")
    status("S3_concordance", "DONE", f"{len(conc)} method×target cells; {len(dis)} disagreements in truth_crosscheck.tsv")


if __name__ == "__main__":
    main()
