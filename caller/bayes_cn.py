#!/usr/bin/env python3
"""Bayesian paralog copy-number caller — v0, A STARTING POINT for the team (written overnight 4 Nov 2026).

For a paralog pair (A, B) = (SMN1, SMN2) or (HBA1, HBA2), grid over (cA, cB) in {0..5}^2:

  log P(cA, cB | data) = log NB(depth | cA + cB)                  gene-level depth, MAPQ >= 0
                       + log BetaBinom(nA | nA + nB, p(cA, cB))   PSV allele counts, pooled over sites
                       + log prior(cA, cB)                        population prior (optional)

  depth: c_hat = norm_cn(A) + norm_cn(B) at MAPQ >= 0 (reads multi-mapping between the paralogs are
         conserved in the sum), divided by a per-(family, ref) calibration factor = median(c_hat) / 4
         over all samples (no truth used). Converted to pseudo-reads with the gene length; NB variance is
         moment-matched to the sample's own control-window dispersion scaled to gene length.
  PSV:   nA/nB summed over both paralog positions of each site, divided by a per-site allelic bias rho
         (estimated without truth), pooled, shrunk to --eff-sites sites' worth of reads;
         p = (cA(1-e) + cB e)/(cA + cB), beta-binomial concentration --conc.
  prior: frequencies of DRAGEN v4.2.7 calls among the selected 1kGP samples, leaving the test sample out,
         smoothed with a pseudo-count. HBA: DRAGEN reports totals only -> prior on the total, uniform split.

Evaluation (--eval): against data/truth/diploid.tsv (usable rows only). SMN truth = (SMN1, SMN2 + SMN2Δ7-8);
HBA truth = total α-globin copies. Writes posteriors + metrics + a reliability plot.

  python caller/bayes_cn.py --ref grch38 --eval
  python caller/bayes_cn.py --ref chm13 --no-prior --eval
"""
import argparse
import csv
import gzip
import math
import os
import statistics
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import betabinom, nbinom

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from cfg import REPO_ROOT, path  # noqa: E402

GRID = range(6)
FAMS = {"SMN": ("SMN1", "SMN2"), "HBA": ("HBA1", "HBA2")}
READ_LEN = 150


def tsv(fn):
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_inputs(ref, mapq):
    genes, ctrl = defaultdict(dict), defaultdict(list)
    for r in tsv(path("tables", f"depth.{ref}.tsv.gz")):
        if int(r["min_mapq"]) != mapq:
            continue
        if r["kind"] == "gene" and r["id"] in {g for p in FAMS.values() for g in p}:
            genes[r["sample"]][r["id"]] = (float(r["norm_cn"]), int(r["end"]) - int(r["start"]), float(r["mean_depth"]))
        elif r["kind"] == "control":
            ctrl[r["sample"]].append(float(r["norm_cn"]))
    psv = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in tsv(path("tables", f"psv_counts.{ref}.tsv.gz")):
        if r["family"] == "SMN":
            pass  # all 16 SMNCopyNumberCaller sites; 'selected' ones are the most specific
        psv[(r["sample"], r["family"])][r["site_id"]][0] += int(r["n_baseA"])
        psv[(r["sample"], r["family"])][r["site_id"]][1] += int(r["n_baseB"])
    return genes, ctrl, psv


def dragen_prior(exclude=None, pseudo=0.5):
    fn = path("tables", "dragen_calls.tsv")
    smn, hba = np.full((6, 6), pseudo), np.full(11, pseudo)
    if os.path.exists(fn):
        for r in tsv(fn):
            if r["dragen_version"] != "4.2.7" or r["sample"] == exclude:
                continue
            try:
                a, b = int(r["SMN1_CN"]), int(r["SMN2_CN"])
                if a < 6 and b < 6:
                    smn[a, b] += 1
            except ValueError:
                pass
            if r["HBA_total_copies"].isdigit() and int(r["HBA_total_copies"]) <= 10:
                hba[int(r["HBA_total_copies"])] += 1
    lp_smn = np.log(smn / smn.sum())
    t = hba / hba.sum()
    lp_hba = np.array([[math.log(t[a + b]) - math.log(a + b + 1) for b in GRID] for a in GRID])
    return {"SMN": lp_smn, "HBA": lp_hba}


def site_bias(psv, fam):
    """Per-site allelic bias rho (odds of an A-read vs a B-read per copy), estimated without truth:
    SMN: relative to DRAGEN v4.2.7 SMN1/SMN2 calls; HBA: assuming the median sample is balanced (cA = cB)."""
    expect = {}
    if fam == "SMN" and os.path.exists(path("tables", "dragen_calls.tsv")):
        for r in tsv(path("tables", "dragen_calls.tsv")):
            if r["dragen_version"] == "4.2.7" and r["SMN1_CN"].isdigit() and r["SMN2_CN"].isdigit():
                if int(r["SMN1_CN"]) > 0 and int(r["SMN2_CN"]) > 0:
                    expect[r["sample"]] = int(r["SMN1_CN"]) / int(r["SMN2_CN"])
    ratios = defaultdict(list)
    for (s, f), sites in psv.items():
        if f != fam:
            continue
        e = expect.get(s, 1.0) if fam == "SMN" else 1.0
        if fam == "SMN" and s not in expect:
            continue
        for sid, (x, y) in sites.items():
            if x >= 3 and y >= 3:
                ratios[sid].append((x / y) / e)
    return {sid: statistics.median(v) for sid, v in ratios.items() if len(v) >= 10}


def posterior(fam, g, ctrl, sites, calib, prior, eps=0.01, w_depth=1.0, eff_sites=3.0, conc=150.0, rho=None):
    A, B = FAMS[fam]
    (nA, LA, _), (nB, LB, _) = g[A], g[B]
    c_hat = (nA + nB) / calib
    L = (LA + LB) / 2
    lam = 0.5 * L / READ_LEN * statistics.median([30.0])  # pseudo-reads per copy at 30x; only ratios matter
    cv2 = (statistics.pstdev(ctrl) / statistics.mean(ctrl)) ** 2 * (2000 / L) + 0.03 ** 2 if len(ctrl) > 5 else 0.05 ** 2
    k = max(0, round(c_hat * lam))
    # PSV evidence: bias-correct each site (x / rho), pool over sites, then shrink to `eff_sites` independent
    # sites' worth of reads (one 150-bp read covers several PSVs) and use a beta-binomial for overdispersion.
    use = [(x / (rho or {}).get(sid, 1.0), y) for sid, (x, y) in sites.items()
           if x + y > 0 and (not rho or sid in rho)]  # skip uncalibratable sites (e.g. never show B)
    X, N = sum(x for x, _ in use), sum(x + y for x, y in use)
    shrink = min(1.0, eff_sites / len(use)) if use else 0.0
    n_eff, x_eff = int(round(N * shrink)), int(round(X * shrink))
    lp = np.zeros((6, 6))
    for a in GRID:
        for b in GRID:
            mu = max(0.05, a + b) * lam
            var = max(mu * 1.0001, mu + cv2 * mu * mu)
            n_nb, p_nb = mu * mu / (var - mu), mu / var
            ll = w_depth * nbinom.logpmf(k, n_nb, p_nb)
            if a + b > 0 and n_eff > 0:
                p = (a * (1 - eps) + b * eps) / (a + b)
                ll += betabinom.logpmf(x_eff, n_eff, p * conc, (1 - p) * conc)
            lp[a, b] = ll + (prior[fam][a, b] if prior is not None else 0.0)
    lp -= lp.max()
    post = np.exp(lp)
    return post / post.sum(), c_hat


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="grch38", choices=["grch38", "chm13"])
    ap.add_argument("--mapq", type=int, default=0)
    ap.add_argument("--no-prior", action="store_true")
    ap.add_argument("--eps", type=float, default=0.01)
    ap.add_argument("--eff-sites", type=float, default=3.0, help="independent PSV sites' worth of reads (sites are correlated)")
    ap.add_argument("--conc", type=float, default=150.0, help="beta-binomial concentration for PSV allele fraction")
    ap.add_argument("--eval", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    tag = f"{a.ref}.mq{a.mapq}" + (".noprior" if a.no_prior else "")
    out = a.out or path("tables", f"bayes_cn.{tag}.tsv")

    genes, ctrl, psv = load_inputs(a.ref, a.mapq)
    calib = {}
    for fam, (A, B) in FAMS.items():
        vals = [g[A][0] + g[B][0] for g in genes.values() if A in g and B in g]
        calib[fam] = statistics.median(vals) / 4 if vals else 1.0
    rho = {fam: site_bias(psv, fam) for fam in FAMS}
    print("site bias rho:", {f: {k: round(v, 2) for k, v in d.items()} for f, d in rho.items()})
    rows = []
    for s in sorted(genes):
        prior = None if a.no_prior else dragen_prior(exclude=s)  # leave-one-out prior
        for fam, (A, B) in FAMS.items():
            if A not in genes[s] or B not in genes[s]:
                continue
            post, c_hat = posterior(fam, genes[s], ctrl[s], psv.get((s, fam), {}), calib[fam], prior, a.eps, eff_sites=a.eff_sites, conc=a.conc, rho=rho[fam])
            i, j = np.unravel_index(post.argmax(), post.shape)
            tot = np.array([sum(post[x, t - x] for x in GRID if 0 <= t - x < 6) for t in range(11)])
            rows.append(dict(sample=s, family=fam, ref=a.ref, c_hat=round(c_hat, 3), map_A=i, map_B=j, p_map=round(post[i, j], 4),
                             map_total=int(tot.argmax()), p_total=round(tot.max(), 4),
                             post_A=",".join(f"{v:.4f}" for v in post.sum(1)), post_total=",".join(f"{v:.4f}" for v in tot)))
    cols = list(rows[0].keys())
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {out} ({len(rows)} rows); calibration {calib}")
    if a.eval:
        evaluate(rows, tag)


def evaluate(rows, tag):
    fn = path("truth", "diploid.tsv")
    if not os.path.exists(fn):
        print("no truth yet")
        return
    truth = {r["sample"]: r for r in tsv(fn)}
    res, rel = [], []
    for r in rows:
        t = truth.get(r["sample"])
        if not t or t[f"usable_{r['family']}"] != "True":
            continue
        if r["family"] == "SMN":
            tA, tB = int(t["SMN1"]), int(t["SMN2"]) + int(t["SMN_del7-8"])
            pA = [float(x) for x in r["post_A"].split(",")]
            res.append(("SMN1", r["map_A"] == tA, sum((pA[k] - (k == tA)) ** 2 for k in range(6)) if tA < 6 else 2.0))
            res.append(("SMN1+SMN2", (r["map_A"], r["map_B"]) == (tA, tB), None))
            rel.append((r["p_map"], (r["map_A"], r["map_B"]) == (tA, tB)))
        else:
            tt = int(t["HBA_total"])
            pt = [float(x) for x in r["post_total"].split(",")]
            res.append(("HBA_total", r["map_total"] == tt, sum((pt[k] - (k == tt)) ** 2 for k in range(11)) if tt < 11 else 2.0))
            rel.append((r["p_total"], r["map_total"] == tt))
    mfn = path("tables", f"bayes_cn_eval.{tag}.tsv")
    with open(mfn, "w") as fh:
        fh.write("target\tn\taccuracy\tbrier\n")
        for k in sorted({x[0] for x in res}):
            xs = [x for x in res if x[0] == k]
            br = [x[2] for x in xs if x[2] is not None]
            line = f"{k}\t{len(xs)}\t{sum(x[1] for x in xs) / len(xs):.3f}\t{(sum(br) / len(br)) if br else float('nan'):.3f}"
            fh.write(line + "\n")
            print(line)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bins = np.linspace(0, 1, 6)
    xs, ys, ns = [], [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = [ok for p, ok in rel if lo <= p < hi or (hi == 1 and p == 1)]
        if sel:
            xs.append((lo + hi) / 2); ys.append(sum(sel) / len(sel)); ns.append(len(sel))
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot([0, 1], [0, 1], "k:", lw=0.7)
    ax.plot(xs, ys, "o-", color="#1f77b4")
    for x, y, n in zip(xs, ys, ns):
        ax.annotate(str(n), (x, y), fontsize=7, xytext=(3, -10), textcoords="offset points")
    ax.set_xlabel("posterior of MAP call"); ax.set_ylabel("fraction correct"); ax.set_title(f"Reliability ({tag})")
    fig.tight_layout()
    fig.savefig(os.path.join(REPO_ROOT, "reports", f"fig_bayes_reliability_{tag}.png"), dpi=130)


if __name__ == "__main__":
    main()
