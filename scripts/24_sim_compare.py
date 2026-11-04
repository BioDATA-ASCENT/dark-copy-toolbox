#!/usr/bin/env python3
"""S14 — genotype the simulated patients with BOTH alignment-free methods, on identical reads.

  counter : the toolbox's in-browser copy-number counter (31-mer panel from GRCh38; two ratios:
            total = round(SMN reads / reads-per-copy), SMN1 = round(total × SMN1 share at the twin sites))
  ctyper  : ctyper 1.0.5 with the same 4-group matrix as the real runs (scripts/17, 18)

Reads: docs/sim/patient_{A,B,C}.fq.gz (scripts/20). Coverage 10x = all reads (one dataset per patient);
5x and 2x = `eval_reps` random read subsamples each (same subsample for both methods); 30x = `high_reps` fresh
simulations per patient with scripts/20's simulator (standard WGS depth, where ctyper is designed to work; not saved).

ctyper input is FASTA: in this ctyper build the gzipped-FASTQ path inflates k-mer counts ~1000-fold (checked on
HG00097's real reads: FASTQ coefficients in the thousands vs ~1 from the indexed BAM), and an unaligned BAM can't
use the indexed-BAM path. FASTA gives coefficients ~1, as on real data.
ctyper's depth input d (31-mer depth of one copy) is measured from the reads themselves, from the single-copy
control gene MARVELD2, as the counter does: no truth and no nominal coverage is used.

Outputs: data/tables/sim_compare.tsv, reports/explorer/sim_compare.json (merged into the toolbox by 14)
"""
import gzip
import importlib.util
import json
import os
import random
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from cfg import CFG, REPO_ROOT, path
from common import log, revcomp, status

S = CFG["simulation"]
K, RL = 31, S["read_length"]
MATRIX = path("ctyper", "subset_matrix.txt")
TMP = path("ctyper", "simtmp")
os.makedirs(TMP, exist_ok=True)


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(REPO_ROOT, "scripts", fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CT = load("ctyper_run", "18_ctyper_run.py")
SIM = load("simulate_patients", "20_simulate_patients.py")


def read_fastq(fn):
    with gzip.open(fn, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq = fh.readline().strip()
            fh.readline(); fh.readline()
            yield h[1:].strip(), seq


def panel_lut(panel):
    lut = {}
    for k in panel["ctrl"]:
        lut[k] = lut[revcomp(k)] = ("c", None)
    for k in panel["gene"]:
        lut[k] = lut[revcomp(k)] = ("g", None)
    for sid, (a, b) in panel["sites"].items():
        lut[a] = lut[revcomp(a)] = ("A", sid)
        lut[b] = lut[revcomp(b)] = ("B", sid)
    return lut


def hit(seq, lut):
    c = g = False
    st = {}
    for i in range(len(seq) - K + 1):
        x = lut.get(seq[i:i + K])
        if x:
            if x[0] == "c":
                c = True
            elif x[0] == "g":
                g = True
            else:
                st[x[1]] = x[0]
    return c, g, st


def annotate(seq, lut):
    """Marker positions in one read, for the toolbox's read viewer: merged k-mer spans + the twin letters."""
    spans, sites = [], []
    for i in range(len(seq) - K + 1):
        x = lut.get(seq[i:i + K])
        if not x:
            continue
        if x[0] in "cg":
            if spans and spans[-1][2] == x[0] and spans[-1][1] >= i:
                spans[-1][1] = i + K
            else:
                spans.append([i, i + K, x[0]])
        else:
            sites.append([i + K // 2, seq[i + K // 2], x[0], x[1]])
    return spans, sites


def examples(reads, hits, lut, n_each=2):
    """A few reads of each kind: control marker, SMN marker, SMN1 letter, SMN2 letter (first ones in the file)."""
    want = {"c": n_each, "g": n_each, "A": n_each, "B": n_each}
    out = []
    for seq, (c, g, st) in zip(reads, hits):
        kind = "A" if "A" in st.values() else "B" if "B" in st.values() else "g" if g else "c" if c else None
        if kind and want[kind] > 0:
            want[kind] -= 1
            spans, sites = annotate(seq, lut)
            out.append(dict(kind=kind, seq=seq, spans=spans, sites=sites))
        if not any(want.values()):
            break
    order = {"c": 0, "g": 1, "A": 2, "B": 3}
    return sorted(out, key=lambda e: order[e["kind"]])


def counter_call(hits, panel):
    """Exactly the browser's genotypeCounts()."""
    Rc = sum(h[0] for h in hits); Rg = sum(h[1] for h in hits)
    X = sum(v == "A" for h in hits for v in h[2].values()); N = sum(len(h[2]) for h in hits)
    per = Rc / 2 * panel["L"]["gene"] / panel["L"]["ctrl"]
    total = Rg / per if per else float("nan")
    T = max(0, round(total)) if per else 0
    share = X / N if N else float("nan")
    A = min(T, round(T * share)) if N else 0
    return dict(SMN1=A, SMN2=T - A, total=round(total, 3), share=round(share, 3) if N else None, Rc=Rc, Rg=Rg, X=X, N=N)


def kmer_depth(hits, panel):
    """31-mer depth of ONE copy, measured from control-gene reads: reads touching MARVELD2 per copy -> read depth."""
    Rc = sum(h[0] for h in hits)
    L_eff = panel["L"]["ctrl"] + RL - 2 * K          # a read is counted when it carries >= 1 control 31-mer
    read_depth_per_copy = (Rc / 2) * RL / L_eff
    return round(read_depth_per_copy * (RL - K + 1) / RL, 3)


def run_ctyper(seqs, d, genes_of):
    td = tempfile.mkdtemp(dir=TMP)
    try:
        fa, out = os.path.join(td, "r.fa"), os.path.join(td, "r.out")
        with open(fa, "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">r{i}\n{s}\n")
        r = subprocess.run(["conda", "run", "-n", CFG["ctyper"]["conda_env"], "ctyper", "-i", fa, "-m", MATRIX, "-o", out,
                            "-d", str(d), "-n", "1"], capture_output=True, text=True)
        if r.returncode or not os.path.exists(out):
            return None
        cn, res = CT.parse(out, genes_of)
        coefs = {}
        cur = None
        for line in open(out):
            if line.startswith(">"):
                cur = line[2:].split()[0]
            elif line.startswith("regress:") and cur == "SMN_group1":
                for x in line.split(":", 1)[1].strip().strip(",").split(","):
                    if x:
                        a, v = x.rsplit(":", 1)
                        coefs[a] = float(v)
        alleles = res.get("SMN_group1", [])
        return dict(SMN1=cn["SMN1"], SMN2=cn["SMN2"], partial=cn["SMN1_partial"] + cn["SMN2_partial"],
                    alleles=[[a, "+".join(genes_of.get(a, [])) or "?", round(coefs.get(a, 0), 2)] for a in alleles],
                    n_alleles=len(alleles))
    finally:
        shutil.rmtree(td, ignore_errors=True)


def main():
    sim = json.load(open(os.path.join(REPO_ROOT, "reports", "explorer", "sim.json")))
    panel = sim["panel"]
    lut = panel_lut(panel)
    genes_of = CT.allele_genes()
    jobs = []   # (patient, coverage, rep, seqs, hits)
    ex_reads = {}
    for pt in sim["patients"]:
        reads = [s for _, s in read_fastq(os.path.join(REPO_ROOT, "docs", pt["fastq"]))]
        hits = [hit(s, lut) for s in reads]
        ex_reads[pt["id"]] = dict(examples=examples(reads, hits, lut), n_reads=len(reads),
                                  n_marked=sum(1 for h in hits if h[0] or h[1] or h[2]))
        jobs.append((pt, pt["coverage"], 0, reads, hits))
        for cov in S["eval_coverages"]:
            if cov >= pt["coverage"]:
                continue
            frac = cov / pt["coverage"]
            for rep in range(S["eval_reps"]):
                rng = random.Random(1000 + 97 * rep + cov)
                keep = [i for i in range(len(reads)) if rng.random() < frac]
                jobs.append((pt, cov, rep, [reads[i] for i in keep], [hits[i] for i in keep]))
    haps = SIM.read_fa(path("pangraph", "SMN.fa"))
    byid = {k.split(":")[0]: v for k, v in haps.items()}
    hc = S.get("high_coverage")
    for rep in range(S.get("high_reps", 0) if hc else 0):
        for pt in sim["patients"]:
            rng = random.Random(S["seed"] + ord(pt["id"]) + 1000 * (rep + 1))
            reads = [s for _, s in SIM.simulate([byid[h] for h in pt["haplotypes"]], hc, rng)]
            jobs.append((pt, hc, rep, reads, [hit(s, lut) for s in reads]))
    log(f"sim compare: {len(jobs)} read sets")

    def one(job):
        pt, cov, rep, seqs, hs = job
        c = counter_call(hs, panel)
        d = kmer_depth(hs, panel)
        y = run_ctyper(seqs, d, genes_of)
        return pt["id"], cov, rep, len(seqs), c, d, y

    rows = []
    with ThreadPoolExecutor(max_workers=CFG["ctyper"].get("jobs", 16)) as ex:
        for f in as_completed([ex.submit(one, j) for j in jobs]):
            rows.append(f.result())
    rows.sort(key=lambda r: (r[0], -r[1], r[2]))
    truth = {p["id"]: (p["SMN1"], p["SMN2"]) for p in sim["patients"]}

    with open(path("tables", "sim_compare.tsv"), "w") as fh:
        fh.write("patient\tcoverage\trep\tn_reads\tcounter_SMN1\tcounter_SMN2\tcounter_total_raw\tcounter_share\t"
                 "ctyper_d\tctyper_SMN1\tctyper_SMN2\tctyper_partial\tctyper_alleles\ttruth_SMN1\ttruth_SMN2\n")
        for p, cov, rep, n, c, d, y in rows:
            y = y or {}
            fh.write("\t".join(map(str, [p, cov, rep, n, c["SMN1"], c["SMN2"], c["total"], c["share"], d,
                                         y.get("SMN1", ""), y.get("SMN2", ""), y.get("partial", ""),
                                         ";".join(a for a, _, _ in y.get("alleles", [])), *truth[p]])) + "\n")

    # summary for the page
    acc = defaultdict(lambda: Counter())
    runs = defaultdict(lambda: defaultdict(list))
    for p, cov, rep, n, c, d, y in rows:
        t = truth[p]
        acc[cov]["n"] += 1
        acc[cov]["counter"] += (c["SMN1"], c["SMN2"]) == t
        acc[cov]["counter_smn1"] += c["SMN1"] == t[0]
        if y:
            ys2 = y["SMN2"] + y["partial"]
            acc[cov]["ctyper"] += (y["SMN1"], ys2) == t
            acc[cov]["ctyper_smn1"] += y["SMN1"] == t[0]
        runs[p][str(cov)].append(dict(rep=rep, n=n, counter=c, d=d,
                                      ctyper=dict(SMN1=y["SMN1"], SMN2=y["SMN2"] + y["partial"], alleles=y["alleles"]) if y else None))
    out = dict(runs=runs, reads=ex_reads, accuracy={str(k): dict(v) for k, v in sorted(acc.items(), reverse=True)},
               note="30x = 3 fresh simulations per patient; 10x = all published simulated reads (one read set per patient); 5x and 2x = 20 random subsamples each. "
                    "Exact = SMN1 and SMN2 both right (ctyper SMN2 includes partial copies).")
    with open(os.path.join(REPO_ROOT, "reports", "explorer", "sim_compare.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    status("S14_sim_compare", "DONE", "; ".join(
        f"{k}x: counter {v['counter']}/{v['n']}, ctyper {v.get('ctyper', 0)}/{v['n']}" for k, v in sorted(acc.items(), reverse=True)))


if __name__ == "__main__":
    main()
