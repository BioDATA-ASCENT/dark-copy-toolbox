#!/usr/bin/env python3
"""S7 — figures for the local pangenome graphs of SMN and HBA.

1. reports/fig_pangraph_{fam}_haplotypes.png: every graph path drawn as a line with its genes (from CAT, re-labelled
   with our sequence-based truth call for SMN1/SMN2 / HBA1/HBA2 copies), aligned on the left anchor.
2. Copies pggb's own odgi renderings (1D `viz` and 2D `draw`) to reports/fig_pangraph_{fam}_odgi_{1d,2d}.png.
3. reports/explorer/pangraph.json: the same layout for the web explorer.
"""
import base64
import csv
import glob
import json
import os
import re
import shutil
import subprocess
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

from cfg import CFG, REPO_ROOT, path
from common import log, status

P = path("pangraph")
PGGB_ENV = CFG.get("pangraph", {}).get("conda_env", "pggb074")
R = os.path.join(REPO_ROOT, "reports")
COL = {"SMN1": "#2a78d6", "SMN2": "#eb6834", "SMN2_del7-8": "#f2a06c", "SMN?": "#999999", "NAIP": "#6b7a8f",
       "SERF1A": "#b9c4d3", "SERF1B": "#b9c4d3", "GTF2H2": "#8fa3bf", "GTF2H2C": "#8fa3bf",
       "HBA1": "#2a78d6", "HBA2": "#eb6834", "HBA_hybrid": "#9b59b6", "HBZ": "#8fa3bf", "HBM": "#b9c4d3", "HBQ1": "#6b7a8f"}
ANCH = "#d0d4da"
NAME = re.compile(r"^(?P<s>[^#]+)#(?P<h>[^#]+)#(?P<c>[^:]+):(?P<lo>\d+)-(?P<hi>\d+)(?P<f>\(-\))?$")


def truth_copies():
    """(sample, hap, contig) -> [(start, end, call)] in contig coordinates (0-based)."""
    d = defaultdict(list)
    for r in csv.DictReader(open(path("truth", "haplotype_copies.tsv")), delimiter="\t"):
        call = r["call"]
        if r["family"] == "SMN" and r.get("exon78") == "absent":
            call = "SMN2_del7-8" if call == "SMN2" else call
        d[(r["sample"], r["hap"], r["contig"])].append((int(r["start"]), int(r["end"]), call, r["family"]))
    return d


def layout(fam, tc):
    paths = defaultdict(list)
    tag = {}
    for r in csv.DictReader(open(os.path.join(P, f"{fam}.genes.tsv")), delimiter="\t"):
        if r["gene"] == "#truth":
            tag[r["path"]] = r["strand"]
            continue
        paths[r["path"]].append([r["gene"], int(r["start"]), int(r["end"]), r["strand"]])
    out = []
    for p, gs in paths.items():
        m = NAME.match(p)
        L = int(m["hi"]) - int(m["lo"]) + 1
        if not m["s"].isupper() or m["s"] in ("GRCH38", "CHM13"):
            pass
        # replace the family's core copies by our truth calls where we have them (assemblies only)
        tcs = tc.get((m["s"], m["h"], m["c"]), [])
        tcs = [t for t in tcs if t[3] == fam]
        if tcs:
            core = {"SMN": {"SMN1", "SMN2"}, "HBA": {"HBA1", "HBA2"}}[fam]
            gs = [g for g in gs if g[0] not in core]
            lo = int(m["lo"]) - 1
            for s, e, call, _ in tcs:
                s, e = s - lo, e - lo
                if m["f"]:
                    s, e = L - e, L - s
                gs.append([call, s, e, "."])
        out.append(dict(path=p, label=f"{m['s']}#{m['h']}" if m["s"] not in ("GRCH38", "CHM13") else m["s"],
                        length=L, truth=tag.get(p, "reference"), genes=sorted(gs, key=lambda g: g[1])))
    refs = [o for o in out if o["truth"] == "reference"]
    rest = sorted([o for o in out if o["truth"] != "reference"], key=lambda o: (o["truth"], o["length"]))
    return refs + rest


def plot(fam, rows):
    h = 0.32 * len(rows) + 1.2
    fig, ax = plt.subplots(figsize=(13, h))
    maxL = max(r["length"] for r in rows)
    for i, r in enumerate(rows):
        y = len(rows) - i
        ax.plot([0, r["length"] / 1e3], [y, y], color="#555", lw=0.8, zorder=1)
        for g, s, e, st in r["genes"]:
            c = COL.get(g, ANCH)
            ax.add_patch(plt.Rectangle((s / 1e3, y - 0.28), max((e - s) / 1e3, maxL / 1e3 / 400), 0.56,
                                       color=c, lw=0, zorder=2))
        ax.text(-maxL / 1e3 * 0.01, y, r["label"], ha="right", va="center", fontsize=7, family="monospace")
        ax.text(r["length"] / 1e3 + maxL / 1e3 * 0.01, y, r["truth"], ha="left", va="center", fontsize=7, color="#333")
    used = sorted({g[0] for r in rows for g in r["genes"]})
    handles = [plt.Rectangle((0, 0), 1, 1, color=COL.get(g, ANCH)) for g in used]
    ax.legend(handles, used, ncol=min(8, len(used)), fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.06 * 12 / h),
              frameon=False)
    ax.set_ylim(0.3, len(rows) + 0.7)
    ax.set_xlim(-maxL / 1e3 * 0.16, maxL / 1e3 * 1.18)
    ax.set_yticks([])
    for s in ["left", "right", "top"]:
        ax.spines[s].set_visible(False)
    ax.set_xlabel("kb from the left unique anchor (paths oriented like GRCh38)")
    ax.set_title(f"{fam} locus: structure of every path in the local pangenome graph (truth call on the right)",
                 fontsize=10, loc="left")
    fig.tight_layout()
    fn = os.path.join(R, f"fig_pangraph_{fam}_haplotypes.png")
    fig.savefig(fn, dpi=150)
    plt.close(fig)
    return fn


def main():
    tc = truth_copies()
    data = {}
    for fam in ["SMN", "HBA"]:
        if not os.path.exists(os.path.join(P, f"{fam}.genes.tsv")):
            continue
        rows = layout(fam, tc)
        log(f"pangraph fig: {plot(fam, rows)}")
        g = {"paths": rows}
        d = os.path.join(P, f"pggb_{fam}")
        for kind, pat in [("1d", "*.final.og.viz_depth_multiqc.png"), ("2d", "*.final.og.lay.draw.png")]:
            hits = sorted(glob.glob(os.path.join(d, pat)))
            if hits:
                dst = os.path.join(R, f"fig_pangraph_{fam}_odgi_{kind}.png")
                shutil.copy(hits[0], dst)
                g[f"odgi_{kind}"] = "data:image/png;base64," + base64.b64encode(open(dst, "rb").read()).decode()
        og = sorted(glob.glob(os.path.join(d, "*.smooth.final.og")))
        if og:
            out = subprocess.run(["conda", "run", "-n", PGGB_ENV, "odgi", "stats", "-i", og[0], "-S"],
                                 capture_output=True, text=True).stdout.strip().splitlines()
            vals = dict(zip(out[0].lstrip("#").split("\t"), map(int, out[1].split("\t")))) if len(out) > 1 else {}
            g["graph"] = dict(graph_bp=vals.get("length"), nodes=vals.get("nodes"), edges=vals.get("edges"),
                              path_bp=sum(r["length"] for r in rows), n_paths=len(rows))
        data[fam] = g
    with open(os.path.join(R, "explorer", "pangraph.json"), "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
    status("S7_pangraph_figs", "DONE", ", ".join(f"{k}: {len(v['paths'])} paths" for k, v in data.items()))


if __name__ == "__main__":
    main()
