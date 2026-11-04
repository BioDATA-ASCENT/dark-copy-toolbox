#!/usr/bin/env python3
"""Export everything the web explorer needs into one compact JSON (reports/explorer/data.json) and
inline it into reports/explorer/index.html (built from explorer/template.html).

Contents: locus geometry (genes, PSVs) per reference, depth profiles (MAPQ>=0 / >=20, re-binned for size) for
every truth sample + population quantiles over all sliced samples, per-sample truth / DRAGEN / naive calls,
the Bayesian caller's sufficient statistics (so the page can re-run the caller live with other parameters),
and the concordance table.
"""
import csv
import gzip
import importlib.util
import json
import os
import statistics
from collections import defaultdict

import numpy as np

from cfg import REPO_ROOT, path
from common import log, status

OUT = os.path.join(REPO_ROOT, "reports", "explorer")
os.makedirs(OUT, exist_ok=True)
REBIN = {"SMN": 5, "HBA": 2, "CHR21_FALSEDUP": 2}  # merge this many native bins per displayed bin

spec = importlib.util.spec_from_file_location("bayes", os.path.join(REPO_ROOT, "caller", "bayes_cn.py"))
bayes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bayes)


def tsv(fn):
    if not os.path.exists(fn):
        return []
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def r2(x):
    return None if x is None or x != x else round(x, 2)


def main():
    sel = tsv(path("meta", "samples_selected.tsv"))
    truth = {r["sample"]: r for r in tsv(path("truth", "diploid.tsv"))}
    truth_samples = {r["sample"] for r in sel if r["role"] == "truth"}
    data = {"loci": {}, "profiles": {}, "samples": [], "caller": {}, "concordance": tsv(path("tables", "concordance.tsv"))}

    # ---- locus geometry
    for ref in ["grch38", "chm13"]:
        loci = defaultdict(list)
        for l in open(path("regions", f"{ref}.loci.bed")):
            c, s, e, n = l.rstrip("\n").split("\t")[:4]
            loci[n].append((c, int(s), int(e)))
        genes = defaultdict(list)
        for l in open(path("regions", f"{ref}.genes.bed")):
            c, s, e, n, locus, strand = l.rstrip("\n").split("\t")
            genes[locus.replace("_anchor", "")].append(dict(name=n, chrom=c, start=int(s), end=int(e), strand=strand,
                                                             anchor=locus.endswith("_anchor")))
        psv = defaultdict(list)
        for r in tsv(path("regions", f"{ref}.psv.tsv")):
            psv[r["family"]].append(int(r["pos"]))
        fam_of = {"SMN": "SMN", "HBA": "HBA"}
        data["loci"][ref] = {n: dict(windows=w, genes=genes.get(n, []), psv=sorted(psv.get(fam_of.get(n, ""), [])))
                             for n, w in loci.items()}

    # ---- depth profiles
    for ref in ["grch38", "chm13"]:
        prof = defaultdict(lambda: defaultdict(dict))  # locus -> mq -> sample -> [values]
        coords = defaultdict(list)
        rows = tsv(path("tables", f"depth.{ref}.tsv.gz"))
        by = defaultdict(list)
        for r in rows:
            if r["kind"] == "bin" and r["min_mapq"] in ("0", "20"):
                locus = r["id"].split(":")[0]
                by[(locus, r["min_mapq"], r["sample"])].append((int(r["start"]), int(r["end"]), float(r["norm_cn"])))
        for (locus, q, s), v in by.items():
            v.sort()
            k = REBIN.get(locus, 1)
            vals = [sum(x[2] for x in v[i:i + k]) / len(v[i:i + k]) for i in range(0, len(v), k)]
            prof[locus][q][s] = vals
            if not coords[locus]:
                coords[locus] = [[v[i][0], v[min(i + k, len(v)) - 1][1]] for i in range(0, len(v), k)]
        out = {}
        for locus, d in prof.items():
            e = {"bins": coords[locus], "mq": {}}
            for q, per in d.items():
                arr = np.array(list(per.values()))
                e["mq"][q] = {"q10": [r2(x) for x in np.percentile(arr, 10, axis=0)],
                              "q50": [r2(x) for x in np.percentile(arr, 50, axis=0)],
                              "q90": [r2(x) for x in np.percentile(arr, 90, axis=0)],
                              "samples": {s: [r2(x) for x in v] for s, v in per.items() if s in truth_samples}}
            out[locus] = e
        data["profiles"][ref] = out
        log(f"explorer: {ref} profiles for {len(prof)} loci")

    # ---- caller sufficient statistics (mirrors bayes.posterior)
    for ref in ["grch38", "chm13"]:
        if not os.path.exists(path("tables", f"depth.{ref}.tsv.gz")):
            continue
        genes, ctrl, psv = bayes.load_inputs(ref, 0)
        rho = {fam: bayes.site_bias(psv, fam) for fam in bayes.FAMS}
        per = {}
        for s, g in genes.items():
            e = {}
            cv2 = (statistics.pstdev(ctrl[s]) / statistics.mean(ctrl[s])) ** 2 if len(ctrl[s]) > 5 else None
            for fam, (A, B) in bayes.FAMS.items():
                if A not in g or B not in g:
                    continue
                sites = psv.get((s, fam), {})
                use = [(x / rho[fam].get(sid, 1.0), y) for sid, (x, y) in sites.items() if x + y > 0 and (not rho[fam] or sid in rho[fam])]
                e[fam] = dict(sum_norm=round(g[A][0] + g[B][0], 4), nA=round(g[A][0], 3), nB=round(g[B][0], 3),
                              L=(g[A][1] + g[B][1]) / 2, cv2=cv2, X=round(sum(x for x, _ in use), 2),
                              N=round(sum(x + y for x, y in use), 2), nsites=len(use))
            per[s] = e
        data["caller"][ref] = per

    # ---- per-sample table
    dragen = defaultdict(dict)
    for r in tsv(path("tables", "dragen_calls.tsv")):
        dragen[r["sample"]][r["dragen_version"]] = {k: r[k] for k in ["SMN1_CN", "SMN2_CN", "SMN2_del78_CN", "SMN_total_CN_raw",
                                                                      "HBA_genotype", "HBA_total_copies"]}
    cty = defaultdict(dict)
    for r in tsv(path("tables", "ctyper_calls.tsv")):
        cty[r["sample"]][r["ref"]] = {k: r[k] for k in ["SMN1", "SMN2", "SMN_partial", "SMN_alleles", "HBA_total", "HBA_alleles", "CBS", "U2AF1"]}
    for r in sel:
        s = r["sample"]
        t = truth.get(s)
        tr = None
        if t:
            tr = {k: t[k] for k in ["n_haps", "SMN1", "SMN2", "SMN_del7-8", "SMN?", "SMN_converted", "HBA1", "HBA2",
                                    "HBA_hybrid", "HBA_total", "usable_SMN", "usable_HBA", "flags"]}
        data["samples"].append(dict(id=s, role=r["role"], pop=r["population"], superpop=r["superpop"], sex=r["sex"],
                                    truth=tr, dragen=dragen.get(s, {}), ctyper=cty.get(s) or None))
    pg = os.path.join(OUT, "pangraph.json")
    data["pangraph"] = json.load(open(pg)) if os.path.exists(pg) else {}
    data["meta"] = dict(n_selected=len(sel), n_truth=len(truth), refs=list(data["caller"].keys()),
                        native_bin_bp={"SMN": 1000, "HBA": 250, "CHR21_FALSEDUP": 1000}, rebin=REBIN)
    js = json.dumps(data, separators=(",", ":"))
    with open(os.path.join(OUT, "data.json"), "w") as fh:
        fh.write(js)
    tpl = os.path.join(REPO_ROOT, "explorer", "template.html")
    if os.path.exists(tpl):
        html = open(tpl).read().replace("/*__DATA__*/null", js)
        with open(os.path.join(OUT, "index.html"), "w") as fh:
            fh.write(html)
    status("S6_explorer", "DONE", f"reports/explorer/data.json ({len(js) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
