#!/usr/bin/env python3
"""Phase 10 — reports/MORNING_REPORT.md + first-look figures. Safe to run at any point; it reports what exists."""
import csv
import gzip
import os
import shutil
import subprocess
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cfg import REPO_ROOT, logs, path  # noqa: E402
from common import status  # noqa: E402

REP = os.path.join(REPO_ROOT, "reports")
os.makedirs(REP, exist_ok=True)


def tsv(fn):
    if not os.path.exists(fn):
        return []
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def du(p):
    try:
        return subprocess.run(["du", "-sh", p], capture_output=True, text=True).stdout.split("\t")[0]
    except Exception:  # noqa: BLE001
        return "?"


def fig_profiles(depth, ref, locus, samples, fn):
    fig, axes = plt.subplots(len(samples), 1, figsize=(10, 2.2 * len(samples)), sharex=True, squeeze=False)
    for ax, s in zip(axes[:, 0], samples):
        for q, colr in [("0", "#7f7f7f"), ("20", "#1f77b4")]:
            rows = [r for r in depth if r["sample"] == s and r["kind"] == "bin" and r["id"].startswith(locus + ":") and r["min_mapq"] == q]
            rows.sort(key=lambda r: int(r["start"]))
            if rows:
                ax.plot([int(r["start"]) / 1e6 for r in rows], [float(r["norm_cn"]) for r in rows], lw=0.7, color=colr,
                        label=f"MAPQ≥{q}")
        ax.axhline(2, color="k", lw=0.5, ls=":")
        ax.set_ylim(0, 8)
        ax.set_ylabel(f"{s}\nnorm CN")
    axes[0, 0].legend(loc="upper right", fontsize=7)
    axes[-1, 0].set_xlabel(f"{ref} position (Mb)")
    fig.suptitle(f"{locus} depth profile on {ref}")
    fig.tight_layout()
    fig.savefig(fn, dpi=130)
    plt.close(fig)


def main():
    lines = ["# Morning report — dark dosage", ""]
    st = os.path.join(REPO_ROOT, "STATUS.md")
    if os.path.exists(st):
        lines += ["## Phase status", "", *[l.rstrip() for l in open(st) if l.startswith("|")], ""]

    sel = tsv(path("meta", "samples_selected.tsv"))
    roles = defaultdict(int)
    for r in sel:
        roles[r["role"]] += 1
    sl = tsv(logs("slices.tsv"))
    ok = defaultdict(set)
    bad = defaultdict(list)
    for r in sl:
        (ok[r["ref"]].add(r["sample"]) if r["state"] == "OK" else bad[r["ref"]].append(f"{r['sample']}:{r['state']}"))
    both = ok["grch38"] & ok["chm13"]
    lines += ["## Data on disk", "",
              f"- selected samples: {len(sel)} ({dict(roles)})",
              f"- slices OK: GRCh38 {len(ok['grch38'])}, CHM13 {len(ok['chm13'])}, both {len(both)}",
              f"- non-OK slice attempts: GRCh38 {len(bad['grch38'])}, CHM13 {len(bad['chm13'])} (see logs/slices.tsv)",
              *[f"- `{d}`: {du(path(d))}" for d in ["refs", "slices", "precomputed", "truth", "tables"] if os.path.exists(path(d))],
              ""]

    dip = tsv(path("truth", "diploid.tsv"))
    if dip:
        usable = [r for r in dip if r["usable_SMN"] == "True"]
        dist = defaultdict(int)
        for r in usable:
            dist[(r["SMN1"], r["SMN2"])] += 1
        hba = defaultdict(int)
        for r in dip:
            if r["usable_HBA"] == "True":
                hba[r["HBA_total"]] += 1
        lines += ["## Truth (HPRC R2 assemblies)", "",
                  f"- samples typed: {len(dip)}; usable SMN {len(usable)}; usable HBA {sum(hba.values())}",
                  "- SMN1,SMN2 diploid copies: " + ", ".join(f"{k[0]}+{k[1]}: {v}" for k, v in sorted(dist.items())),
                  "- α-globin total copies: " + ", ".join(f"{k}: {v}" for k, v in sorted(hba.items())),
                  "- flags worth reading: `data/truth/diploid.tsv` column `flags`", ""]

    figs = []
    for ref in ["grch38", "chm13"]:
        depth = tsv(path("tables", f"depth.{ref}.tsv.gz"))
        if not depth:
            continue
        samples = sorted({r["sample"] for r in depth})[:3]
        for locus in ["SMN", "HBA"]:
            fn = os.path.join(REP, f"fig_{locus}_{ref}_profile.png")
            fig_profiles(depth, ref, locus, samples, fn)
            figs.append(fn)
        # naive gene-level CN (MAPQ>=20) vs truth, GRCh38 SMN1/SMN2 and HBA1/HBA2
        if dip:
            truth = {r["sample"]: r for r in dip}
            for gene in ["SMN1", "SMN2", "HBA1", "HBA2"]:
                pts = []
                for q in ["0", "20"]:
                    for r in depth:
                        if r["kind"] == "gene" and r["id"] == gene and r["min_mapq"] == q and r["sample"] in truth:
                            pts.append((q, float(truth[r["sample"]][gene]), float(r["norm_cn"])))
                if not pts:
                    continue
                fig, ax = plt.subplots(figsize=(4, 4))
                for q, colr in [("0", "#7f7f7f"), ("20", "#1f77b4")]:
                    xs = [p[1] for p in pts if p[0] == q]
                    ys = [p[2] for p in pts if p[0] == q]
                    ax.scatter(xs, ys, s=10, alpha=0.6, color=colr, label=f"MAPQ≥{q}")
                ax.plot([0, 5], [0, 5], "k:", lw=0.6)
                ax.set_xlabel(f"truth {gene} copies (assembly)")
                ax.set_ylabel(f"naive depth CN on {ref}")
                ax.legend(fontsize=7)
                ax.set_title(f"{gene}: naive read depth vs truth")
                fig.tight_layout()
                fn = os.path.join(REP, f"fig_naive_vs_truth_{gene}_{ref}.png")
                fig.savefig(fn, dpi=130)
                plt.close(fig)
                figs.append(fn)
    if figs:
        lines += ["## First-look figures", "", *[f"- `reports/{os.path.basename(f)}`" for f in figs], ""]

    loc = path("precomputed", "locityper", "our_loci.txt")
    if os.path.exists(loc):
        found = open(loc).read().split()
        lines += ["## Precomputed", "", "- Locityper loci covering our genes: " + (", ".join(found) or "NONE"), ""]
    for key in ["v4_0_3", "v4_2_7", "v4_4_7"]:
        wk = path("precomputed", "dragen", key, "wanted_keys.txt")
        if os.path.exists(wk):
            lines.append(f"- DRAGEN {key}: {sum(1 for _ in open(wk))} per-sample files")
    lines.append("")

    dec = os.path.join(REPO_ROOT, "DECISIONS.md")
    if os.path.exists(dec):
        lines += ["## Decisions made overnight", "", open(dec).read(), ""]
    lines += ["## Free disk", "", f"- {shutil.disk_usage(path('.')).free / 1e9:.0f} GB free under data/", ""]
    with open(os.path.join(REP, "MORNING_REPORT.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    status("P10_report", "DONE", "reports/MORNING_REPORT.md")


if __name__ == "__main__":
    main()
