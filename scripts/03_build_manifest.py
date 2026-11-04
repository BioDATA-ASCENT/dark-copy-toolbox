#!/usr/bin/env python3
"""Phase 3 — sample manifest.

Outputs (data/meta/):
  samples_all.tsv       all 1kGP samples: sex, population, superpop, parents, GRCh38 CRAM URLs, HPRC R2 status
  samples_selected.tsv  priority-ordered subset to slice tonight (role: truth | parent | extra)
"""
import csv
import os
import random
import re
from collections import defaultdict

from cfg import CFG, path
from common import POP2SUPER, SAMPLE_RE, log, status

M = path("meta")


def read_table(fn):
    """Whitespace/tab table with a header line (leading '#' allowed)."""
    rows = []
    with open(fn) as fh:
        header = None
        for line in fh:
            if not line.strip():
                continue
            if header is None:
                header = re.split(r"\t" if "\t" in line else r"\s+", line.lstrip("#").strip())
                continue
            if line.startswith("#"):
                continue
            vals = re.split(r"\t" if "\t" in line else r"\s+", line.rstrip("\n"))
            rows.append(dict(zip(header, vals)))
    return header, rows


def col(header, *needles):
    for h in header:
        hl = h.lower()
        if any(n in hl for n in needles):
            return h
    return None


def main():
    samples = {}

    # Pedigree / population
    for fn in ["20130606_g1k_3202_samples_ped_population.txt", "1kGP.3202_samples.pedigree_info.txt"]:
        fp = os.path.join(M, fn)
        if not os.path.exists(fp):
            continue
        hdr, rows = read_table(fp)
        cs, cf, cm, cx = col(hdr, "sample"), col(hdr, "father"), col(hdr, "mother"), col(hdr, "sex")
        cp, csp = col(hdr, "population"), col(hdr, "superpop")
        if cp == csp:
            cp = next((h for h in hdr if h.lower() == "population"), None)
        for r in rows:
            s = r.get(cs, "")
            if not SAMPLE_RE.fullmatch(s):
                continue
            d = samples.setdefault(s, {"sample": s})
            for key, c in [("father", cf), ("mother", cm), ("sex", cx), ("population", cp), ("superpop", csp)]:
                if c and r.get(c) and r[c] not in ("0", "NA") and not d.get(key):
                    d[key] = r[c]
    log(f"pedigree samples: {len(samples)}")

    # Sequence indices: ENA path + population fallback
    for fn in ["1000G_2504_high_coverage.sequence.index", "1000G_698_related_high_coverage.sequence.index"]:
        fp = os.path.join(M, fn)
        if not os.path.exists(fp):
            log(f"missing {fn}")
            continue
        hdr, rows = read_table(fp)
        cpath = col(hdr, "ena_file_path", "file_path") or hdr[0]
        csam = col(hdr, "sample_name") or col(hdr, "sample")
        cpop = next((h for h in hdr if h.lower() == "population"), None)
        for r in rows:
            p = r.get(cpath, "")
            m = SAMPLE_RE.search(r.get(csam, "") or p)
            if not m or not p.endswith(".cram"):
                continue
            d = samples.setdefault(m.group(0), {"sample": m.group(0)})
            d["grch38_cram_ena"] = re.sub(r"^ftp://", "http://", p)
            if cpop and r.get(cpop) and not d.get("population"):
                d["population"] = r[cpop]

    # S3 listing (preferred, anonymous HTTPS)
    s3p = CFG["urls"]["kgp_s3_prefix"]
    bucket = s3p[5:].split("/")[0]
    lst = os.path.join(M, "kgp_s3_crams.txt")
    if os.path.exists(lst):
        for key in open(lst):
            key = key.strip()
            m = SAMPLE_RE.search(os.path.basename(key))
            if not m:
                continue
            d = samples.setdefault(m.group(0), {"sample": m.group(0)})
            url = f"https://{bucket}.s3.amazonaws.com/{key}"
            if key.endswith(".cram"):
                d["grch38_cram_s3"] = url
            elif key.endswith(".crai"):
                d["grch38_crai_s3"] = url

    for d in samples.values():
        mapped = POP2SUPER.get(d.get("population", ""))
        if mapped:
            d["superpop"] = mapped
        elif d.get("superpop") not in {"AFR", "AMR", "EAS", "EUR", "SAS"}:
            d["superpop"] = ""

    # HPRC release 2
    hprc = defaultdict(list)
    with open(os.path.join(M, "hprc_r2_assemblies.csv")) as fh:
        for r in csv.DictReader(fh):
            hprc[r["sample_id"]].append(r["assembly_name"])
    for s, d in samples.items():
        d["in_hprc_r2"] = "yes" if s in hprc else "no"
        d["hprc_haplotypes"] = ",".join(sorted(hprc.get(s, [])))
    in_both = sorted(s for s in hprc if s in samples)
    not_in_kgp = sorted(s for s in hprc if s not in samples)
    log(f"HPRC R2 samples: {len(hprc)}; in 1kGP: {len(in_both)}; not in 1kGP: {len(not_in_kgp)}")

    cols = ["sample", "sex", "population", "superpop", "father", "mother", "in_hprc_r2", "hprc_haplotypes",
            "grch38_cram_s3", "grch38_crai_s3", "grch38_cram_ena"]
    with open(os.path.join(M, "samples_all.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in sorted(samples):
            fh.write("\t".join(samples[s].get(c, "") for c in cols) + "\n")

    # Selection: truth samples round-robin over superpops, then their parents, then random extras.
    sc = CFG["samples"]
    rng = random.Random(sc["seed"])
    by_sp = defaultdict(list)
    for s in in_both:
        if samples[s].get("grch38_cram_s3") or samples[s].get("grch38_cram_ena"):
            by_sp[samples[s].get("superpop") or "UNK"].append(s)
    for v in by_sp.values():
        rng.shuffle(v)
    order, chosen = [], set()
    while any(by_sp.values()):
        for sp in sorted(by_sp):
            if by_sp[sp]:
                s = by_sp[sp].pop()
                order.append((s, "truth"))
                chosen.add(s)
    if sc.get("include_trio_parents"):
        for s, _ in list(order):
            for p in (samples[s].get("father"), samples[s].get("mother")):
                if p and p in samples and p not in chosen:
                    order.append((p, "parent"))
                    chosen.add(p)
    extra_pool = [s for s in samples if s not in chosen and samples[s].get("in_hprc_r2") == "no"]
    rng.shuffle(extra_pool)
    for s in extra_pool[: sc.get("n_extra_random", 0)]:
        order.append((s, "extra"))
    order = order[: sc["max_samples"]]

    sel_cols = ["priority", "sample", "role"] + cols[1:]
    with open(os.path.join(M, "samples_selected.tsv"), "w") as fh:
        fh.write("\t".join(sel_cols) + "\n")
        for i, (s, role) in enumerate(order, 1):
            fh.write("\t".join([str(i), s, role] + [samples[s].get(c, "") for c in cols[1:]]) + "\n")
    n_truth = sum(r == "truth" for _, r in order)
    n_par = sum(r == "parent" for _, r in order)
    no_url = [s for s, _ in order if not samples[s].get("grch38_cram_s3")]
    if no_url:
        log(f"{len(no_url)} selected samples have no S3 CRAM (ENA fallback): {no_url[:10]}")
    status("P3_manifest", "DONE", f"{len(samples)} 1kGP samples; selected {len(order)} ({n_truth} truth, {n_par} parents)")


if __name__ == "__main__":
    main()
