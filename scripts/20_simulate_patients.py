#!/usr/bin/env python3
"""S10 — simulated SMA patients, genotyped alignment-free (k-mer counting; no read mapping anywhere).

1. Patients: each is two real HPRC R2 SMN-region haplotypes (MARVELD2..BDP1, data/pangraph/SMN.fa, oriented like
   GRCh38) with known SMN1/SMN2 copies (config simulation.patients).
2. Reads: 150 bp pairs, fragment N(400, 50), uniform 0.3% substitutions, both strands, at `coverage` (diploid).
   Only the region is simulated, so k-mers are checked for uniqueness within the region only.
3. K-mer panel, built from GRCh38 alone (never from the patients):
     ctrl  31-mers of MARVELD2 (single-copy anchor) occurring once in the GRCh38 window   -> 2-copy yardstick
     gene  31-mers of the SMN1 gene that also occur in SMN2 and exactly twice in the window -> total SMN copies
     sites SMN1 / SMN2 31-mers centred on each of the 16 SMNCopyNumberCaller sites        -> how copies split
   ctrl/gene k-mers are thinned (every `kmer_step`-th position) to keep the panel small for the browser.
4. Counting (same rules as the browser): a read counts once for ctrl / gene if any of its 31-mers (either strand)
   is in the set; once per site for SMN1 or SMN2 k-mers.
5. Genotype: same model as the toolbox's caller lab — depth term NB(c × reads-per-copy) + pooled beta-binomial over
   the sites (shrunk to `eff_sites` independent sites) — on a 0..5 × 0..5 grid, no population prior.

Outputs: docs/sim/patient_{A,B,C}.fq.gz, reports/explorer/sim.json (panel + per-patient counts + truth),
         data/tables/sim_calls.tsv (calls at several coverages by subsampling reads)
"""
import gzip
import json
import math
import os
import random
from collections import Counter

from cfg import CFG, REPO_ROOT, path
from common import log, read_fasta_region, revcomp, status

S = CFG["simulation"]
K = 31
OUT_FQ = os.path.join(REPO_ROOT, "docs", "sim")
os.makedirs(OUT_FQ, exist_ok=True)


def read_fa(fn):
    seqs, name = {}, None
    for line in open(fn):
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def build_panel():
    genes = {l.split("\t")[3]: l.rstrip("\n").split("\t") for l in open(path("regions", "grch38.genes.bed"))}
    fa = path(CFG["ref_fasta"]["grch38"])
    chrom = genes["SMN1"][0]
    w0, w1 = int(genes["MARVELD2"][1]), int(genes["BDP1"][2]) if "BDP1" in genes else 71567820
    win = read_fasta_region(fa, f"{chrom}:{w0 + 1}-{w1}")[1]
    cnt = Counter(win[i:i + K] for i in range(len(win) - K + 1))
    for i in range(len(win) - K + 1):          # count both strands
        cnt[revcomp(win[i:i + K])] += 0
    rc_cnt = Counter(revcomp(win[i:i + K]) for i in range(len(win) - K + 1))
    tot = lambda k: cnt.get(k, 0) + rc_cnt.get(k, 0)

    def gene_seq(g):
        c, s, e = genes[g][0], int(genes[g][1]), int(genes[g][2])
        return read_fasta_region(fa, f"{c}:{s + 1}-{e}")[1]

    step = S["kmer_step"]
    mar = gene_seq("MARVELD2")
    ctrl = [mar[i:i + K] for i in range(0, len(mar) - K + 1, step) if "N" not in mar[i:i + K] and tot(mar[i:i + K]) == 1]
    smn1, smn2 = gene_seq("SMN1"), gene_seq("SMN2")
    s2 = set(smn2[i:i + K] for i in range(len(smn2) - K + 1))
    gene = [smn1[i:i + K] for i in range(0, len(smn1) - K + 1, step)
            if smn1[i:i + K] in s2 and tot(smn1[i:i + K]) == 2 and "N" not in smn1[i:i + K]]
    sites = {}
    for line in open(path("regions", "grch38.psv.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[0] != "SMN":
            continue
        sid, copy, c, pos = f[1], f[2], f[3], int(f[4])
        seq = read_fasta_region(fa, f"{c}:{pos - K // 2}-{pos + K // 2}")[1]
        if f[9] == "-":
            seq = revcomp(seq)
        sites.setdefault(sid, {})[copy] = seq
    sites = {sid: [d["SMN1"], d["SMN2"]] for sid, d in sorted(sites.items())
             if "SMN1" in d and "SMN2" in d and d["SMN1"] != d["SMN2"] and tot(d["SMN1"]) == 1 and tot(d["SMN2"]) == 1}
    L = {"ctrl": len(mar), "gene": len(smn1)}
    log(f"sim panel: ctrl {len(ctrl)} k-mers, gene {len(gene)}, sites {len(sites)}")
    return dict(k=K, ctrl=ctrl, gene=gene, sites=sites, L=L)


def simulate(haps, cov, rng):
    """Yield (name, seq, qual) single reads (both mates of each pair)."""
    rl, fm, fs, err = S["read_length"], S["frag_mean"], S["frag_sd"], S["error_rate"]
    n = 0
    for h in haps:
        n_pairs = int(len(h) * cov / 2 / (2 * rl))     # cov is diploid: each haplotype gets cov/2
        for _ in range(n_pairs):
            fl = max(rl + 20, int(rng.gauss(fm, fs)))
            st = rng.randrange(0, max(1, len(h) - fl))
            frag = h[st:st + fl]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            for mate in (frag[:rl], revcomp(frag)[:rl]):
                if "N" in mate or len(mate) < rl:
                    continue
                m = list(mate)
                for i in range(rl):
                    if rng.random() < err:
                        m[i] = rng.choice([b for b in "ACGT" if b != m[i]])
                n += 1
                yield f"r{n}", "".join(m)


def count_reads(reads, panel):
    lut = {}
    for k in panel["ctrl"]:
        lut[k] = lut[revcomp(k)] = ("c", None)
    for k in panel["gene"]:
        lut[k] = lut[revcomp(k)] = ("g", None)
    for sid, (a, b) in panel["sites"].items():
        lut[a] = lut[revcomp(a)] = ("A", sid)
        lut[b] = lut[revcomp(b)] = ("B", sid)
    hits = []  # per read: (ctrl, gene, {sid: 'A'/'B'})
    for seq in reads:
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
        if c or g or st:
            hits.append((c, g, st))
    return hits


def summarise(hits, frac, seed=11):
    rng = random.Random(seed)
    R_c = R_g = 0
    sites = Counter()
    for c, g, st in hits:
        if frac < 1 and rng.random() > frac:
            continue
        R_c += c
        R_g += g
        for sid, ab in st.items():
            sites[(sid, ab)] += 1
    return R_c, R_g, sites


def lgamma(x):
    return math.lgamma(x)


def nb_logpmf(k, mu, size):
    p = size / (size + mu)
    return lgamma(k + size) - lgamma(size) - lgamma(k + 1) + size * math.log(p) + k * math.log(1 - p)


def bb_logpmf(x, n, a, b):
    lb = lambda u, v: lgamma(u) + lgamma(v) - lgamma(u + v)
    return lgamma(n + 1) - lgamma(x + 1) - lgamma(n - x + 1) + lb(x + a, n - x + b) - lb(a, b)


def genotype(R_c, R_g, sites, panel, eps=0.01, eff=3, conc=150, size=200):
    per_copy = R_c / 2 * panel["L"]["gene"] / panel["L"]["ctrl"]
    X = sum(v for (sid, ab), v in sites.items() if ab == "A")
    N = sum(sites.values())
    ns = len(panel["sites"])
    xE, nE = round(X * eff / ns), round(N * eff / ns)
    post = {}
    for a in range(6):
        for b in range(6):
            mu = max(0.05, a + b) * per_copy
            ll = nb_logpmf(R_g, mu, size)
            if a + b > 0 and nE > 0:
                p = (a * (1 - eps) + b * eps) / (a + b)
                ll += bb_logpmf(xE, nE, p * conc, (1 - p) * conc)
            post[(a, b)] = ll
    m = max(post.values())
    z = sum(math.exp(v - m) for v in post.values())
    post = {k: math.exp(v - m) / z for k, v in post.items()}
    best = max(post, key=post.get)
    return best, post[best], (R_g / per_copy if per_copy else float("nan")), X, N


def main():
    panel = build_panel()
    haps = read_fa(path("pangraph", "SMN.fa"))
    byid = {k.split(":")[0]: v for k, v in haps.items()}
    out = {"panel": panel, "patients": []}
    rows = []
    for pt in S["patients"]:
        rng = random.Random(S["seed"] + ord(pt["id"]))
        hs = [byid[h] for h in pt["haplotypes"]]
        fq = os.path.join(OUT_FQ, f"patient_{pt['id']}.fq.gz")
        reads = []
        with gzip.open(fq, "wt", compresslevel=6) as fh:
            for name, seq in simulate(hs, S["coverage"], rng):
                fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
                reads.append(seq)
        hits = count_reads(reads, panel)
        R_c, R_g, st = summarise(hits, 1.0)
        call, p, chat, X, N = genotype(R_c, R_g, st, panel)
        log(f"patient {pt['id']}: {len(reads)} reads; ctrl {R_c} gene {R_g} → total≈{chat:.2f}; "
            f"SMN1-site reads {X}/{N}; call {call} p={p:.3f}; truth {pt['SMN1']},{pt['SMN2']}")
        for cov in S["eval_coverages"]:
            frac = min(1.0, cov / S["coverage"])
            for rep in range(S["eval_reps"]):
                rc, rg, s2 = summarise(hits, frac, seed=100 + rep)
                c2, p2, ch2, _, _ = genotype(rc, rg, s2, panel)
                rows.append((pt["id"], cov, rep, c2[0], c2[1], round(p2, 4), pt["SMN1"], pt["SMN2"]))
        out["patients"].append(dict(id=pt["id"], label=pt["label"], haplotypes=pt["haplotypes"], SMN1=pt["SMN1"],
                                    SMN2=pt["SMN2"], coverage=S["coverage"], n_reads=len(reads),
                                    fastq=f"sim/patient_{pt['id']}.fq.gz", fastq_mb=round(os.path.getsize(fq) / 1e6, 1),
                                    hits=[[int(c), int(g), st] for c, g, st in hits]))
    with open(os.path.join(REPO_ROOT, "reports", "explorer", "sim.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    with open(path("tables", "sim_calls.tsv"), "w") as fh:
        fh.write("patient\tcoverage\trep\tcall_SMN1\tcall_SMN2\tposterior\ttruth_SMN1\ttruth_SMN2\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    acc = Counter((r[1], r[3] == r[6] and r[4] == r[7]) for r in rows)
    status("S10_sim", "DONE", "; ".join(f"{c}x: {acc[(c, True)]}/{acc[(c, True)] + acc[(c, False)]} exact"
                                       for c in S["eval_coverages"]))


if __name__ == "__main__":
    main()
