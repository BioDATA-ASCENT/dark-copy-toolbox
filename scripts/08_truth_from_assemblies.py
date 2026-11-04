#!/usr/bin/env python3
"""Phase 8 — truth copy numbers from HPRC release-2 haplotype assemblies (no read mapping).

Per haplotype of each selected truth sample:
  1. stream its CAT GFF3 and keep only SMN-locus / HBA-locus genes + anchors (locates the locus),
  2. pull the locus region from the remote bgzipped assembly (samtools faidx + local .fai/.gzi),
  3. find EVERY copy by sequence — GRCh38 SMN1 (exons 1-8) and HBA2 are projected onto the region —
     so copies CAT mislabels or misses are still counted,
  4. type each copy at paralog-specific sites: SMN1 vs SMN2 at c.840 (+7 more sites), exon 7-8 present?,
     HBA1 vs HBA2 vs hybrid from HBA PSVs.

Outputs (data/truth/): haplotype_copies.tsv (one row per gene copy), haplotypes.tsv, diploid.tsv
Needs: P1 (HPRC tables), P2 (GRCh38 FASTA), P5 (data/regions/grch38.psv.tsv).
"""
import csv
import gzip
import io
import os
import re
import tempfile
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from cfg import CFG, path
from common import log, project_positions, read_fasta_region, revcomp, run, s3_to_https, status

T = path("truth")
os.makedirs(os.path.join(T, "cat"), exist_ok=True)
os.makedirs(os.path.join(T, "idx"), exist_ok=True)
GRCH38 = path(CFG["ref_fasta"]["grch38"])
LOC = CFG["loci"]
WANT = {}
for lname in ["SMN", "HBA"]:
    for g in LOC[lname]["genes"] + LOC[lname].get("anchors", []):
        WANT[g] = lname
ATTR = re.compile(r"(?:^|;)(?:gene_name|Name|source_gene_common_name|gene|gene_id)=([^;]+)")
PAD = {"SMN": 200_000, "HBA": 20_000}


def http_lines(url):
    req = urllib.request.Request(url, headers={"User-Agent": "dark-dosage-hackathon/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r:
        stream = gzip.GzipFile(fileobj=r) if url.endswith(".gz") else r
        for raw in io.TextIOWrapper(stream, encoding="utf-8", errors="replace"):
            yield raw


def cat_genes(asm, url):
    out = os.path.join(T, "cat", f"{asm}.genes.tsv")
    if os.path.exists(out + ".ok"):
        return [l.rstrip("\n").split("\t") for l in open(out)]
    rows = []
    for line in http_lines(url):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        names = set(ATTR.findall(f[8]))
        hit = next((n for n in names if n in WANT), None)
        if hit:
            rows.append([f[0], str(int(f[3]) - 1), f[4], f[6], hit, WANT[hit]])
    with open(out, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    open(out + ".ok", "w").close()
    return rows


def local_idx(url, asm, ext):
    dest = os.path.join(T, "idx", f"{asm}{ext}")
    if not os.path.exists(dest):
        req = urllib.request.Request(url, headers={"User-Agent": "dark-dosage-hackathon/1.0"})
        with urllib.request.urlopen(req, timeout=300) as r, open(dest + ".part", "wb") as fh:
            fh.write(r.read())
        os.replace(dest + ".part", dest)
    return dest


def fai_lengths(fai):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}


def extract(asm_url, fai, gzi, region):
    out = run(["samtools", "faidx", "--fai-idx", fai, "--gzi-idx", gzi, asm_url, region])
    return "".join(out.splitlines()[1:]).upper()


def load_queries():
    """GRCh38 query sequences + site offsets for SMN1 and HBA2."""
    psv = [r for r in csv.DictReader(open(path("regions", "grch38.psv.tsv")), delimiter="\t")]
    smn_start = 70924941
    smn_seq = read_fasta_region(GRCH38, f"chr5:{smn_start}-70954173")[1]
    smn_sites = [(r["site_id"], int(r["pos"]) - smn_start, r["base_SMN1_or_HBA1"], r["base_SMN2_or_HBA2"], r["selected"])
                 for r in psv if r["family"] == "SMN" and r["copy"] == "SMN1"]
    exon78 = (70948273 + 1 - smn_start, 70954173 - smn_start)  # SMN_region_38.bed exon 7-8 block
    genes = {l.split("\t")[3]: l.split("\t") for l in open(path("regions", "grch38.genes.bed"))}
    c, s, e = genes["HBA2"][0], int(genes["HBA2"][1]), int(genes["HBA2"][2])
    pad = 300
    hba_start = s - pad + 1
    hba_seq = read_fasta_region(GRCH38, f"{c}:{hba_start}-{e + pad}")[1]
    hba_sites = [(r["site_id"], int(r["pos"]) - hba_start, r["base_SMN1_or_HBA1"], r["base_SMN2_or_HBA2"], "yes")
                 for r in psv if r["family"] == "HBA" and r["copy"] == "HBA2"]
    return dict(SMN=(smn_seq, smn_sites, exon78), HBA=(hba_seq, hba_sites, None))


def type_copies(family, region_seq, q):
    qseq, sites, exon78 = q
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, dir=T) as fh:
        fh.write(">region\n" + region_seq + "\n")
        tfa = fh.name
    try:
        hits = project_positions(qseq, tfa, [o for _, o, *_ in sites])
    finally:
        for p in (tfa, tfa + ".fai"):
            if os.path.exists(p):
                os.unlink(p)
    min_id = 0.9 if family == "SMN" else 0.93
    hits = [h for h in hits if h["identity"] >= min_id and h["qcov"] >= (0.3 if family == "SMN" else 0.6)]
    hits.sort(key=lambda h: -(h["identity"] * (h["tend"] - h["tstart"])))
    kept = []
    for h in hits:
        if all(min(h["tend"], k["tend"]) - max(h["tstart"], k["tstart"]) < 0.5 * (h["tend"] - h["tstart"]) for k in kept):
            kept.append(h)
    copies = []
    for h in sorted(kept, key=lambda h: h["tstart"]):
        calls = []
        for sid, off, b_a, b_b, sel in sites:
            p = h["sites"].get(off)
            if p is None:
                calls.append((sid, "-", sel))
                continue
            b = p[1] if h["strand"] == "+" else revcomp(p[1])
            calls.append((sid, "A" if b == b_a else "B" if b == b_b else "x", sel))
        pattern = "".join(c for _, c, _ in calls)
        rec = dict(tstart=h["tstart"], tend=h["tend"], strand=h["strand"], identity=round(h["identity"], 4),
                   qcov=round(h["qcov"], 3), pattern=pattern)
        if family == "SMN":
            c840 = next(c for sid, c, _ in calls if "splice" in sid)
            sel = [c for _, c, s in calls if s == "yes" and c in "AB"]
            e78 = max(0, min(h["qend"], exon78[1]) - max(h["qstart"], exon78[0])) / (exon78[1] - exon78[0])
            rec["exon78"] = "present" if e78 >= 0.8 else "absent"
            rec["call"] = ("SMN1" if c840 == "A" else "SMN2" if c840 == "B" else "SMN?")
            if rec["exon78"] == "absent":
                rec["call"] = "SMN_del7-8"
            rec["converted"] = "yes" if sel and 0 < sum(c == "A" for c in sel) < len(sel) else "no"
        else:
            ab = [c for c in pattern if c in "AB"]
            fa = sum(c == "A" for c in ab) / len(ab) if ab else float("nan")
            rec["call"] = "HBA1" if fa >= 0.8 else "HBA2" if fa <= 0.2 else ("HBA_hybrid" if ab else "HBA?")
            rec["exon78"], rec["converted"] = "", ""
        copies.append(rec)
    return copies


def process(sample, hap, asm, asm_url, fai_url, gzi_url, cat_url, queries):
    genes = cat_genes(asm, cat_url)
    fai = local_idx(fai_url, asm, ".fai")
    gzi = local_idx(gzi_url, asm, ".gzi")
    lens = fai_lengths(fai)
    # FASTA contigs are PanSN-named (SAMPLE#HAP#CONTIG); CAT GFF3 uses the bare contig name
    alias = {n.split("#")[-1]: n for n in lens}
    alias.update({n: n for n in lens})
    rows, summary = [], {"sample": sample, "hap": hap, "assembly": asm}
    for fam in ["SMN", "HBA"]:
        core = [g for g in genes if g[5] == fam and g[4] in LOC[fam]["genes"]]
        anchors = {g[4] for g in genes if g[5] == fam and g[4] in LOC[fam].get("anchors", [])}
        by_contig = defaultdict(list)
        for g in core:
            by_contig[g[0]].append(g)
        n_copies = defaultdict(int)
        flags = []
        if not by_contig:
            flags.append("locus_not_found")
        if len(by_contig) > 1:
            flags.append(f"multi_contig({len(by_contig)})")
        for contig, gs in by_contig.items():
            lo = max(0, min(int(g[1]) for g in gs) - PAD[fam])
            fa_contig = alias.get(contig, contig)
            hi = min(lens.get(fa_contig, 10**12), max(int(g[2]) for g in gs) + PAD[fam])
            if lo == 0 or hi == lens.get(fa_contig):
                flags.append(f"near_contig_end:{contig}")
            anchors_here = {g[4] for g in genes if g[0] == contig and g[4] in LOC[fam].get("anchors", [])}
            seq = extract(asm_url, fai, gzi, f"{fa_contig}:{lo + 1}-{hi}")
            for k, c in enumerate(type_copies(fam, seq, queries[fam])):
                n_copies[c["call"]] += 1
                if fam == "SMN" and c["converted"] == "yes":
                    n_copies["SMN_converted"] += 1
                rows.append([sample, hap, asm, fam, contig, lo + c["tstart"] + 1, lo + c["tend"], c["strand"],
                             c["identity"], c["qcov"], c["call"], c["exon78"], c["converted"], c["pattern"],
                             ",".join(sorted(anchors_here))])
        if fam == "SMN":
            for key in ["SMN1", "SMN2", "SMN_del7-8", "SMN?", "SMN_converted"]:
                summary[key] = n_copies.get(key, 0)
        else:
            for key in ["HBA1", "HBA2", "HBA_hybrid", "HBA?"]:
                summary[key] = n_copies.get(key, 0)
            summary["HBA_total"] = sum(n_copies.get(k, 0) for k in ["HBA1", "HBA2", "HBA_hybrid", "HBA?"])
        summary[f"{fam}_anchors"] = ",".join(sorted(anchors))
        summary[f"{fam}_flags"] = ";".join(flags)
    return rows, summary


def main():
    sel = [r for r in csv.DictReader(open(path("meta", "samples_selected.tsv")), delimiter="\t") if r["role"] == "truth"]
    asm = defaultdict(list)
    for r in csv.DictReader(open(path("meta", "hprc_r2_assemblies.csv"))):
        asm[r["sample_id"]].append(r)
    cat = {r["assembly_name"]: r["location"] for r in csv.DictReader(open(path("meta", "hprc_r2_cat_annotation.csv")))}
    queries = load_queries()
    log(f"truth: {len(sel)} samples; SMN sites {len(queries['SMN'][1])}, HBA PSVs {len(queries['HBA'][1])}")

    jobs = []
    for r in sel:
        for a in asm.get(r["sample"], []):
            if a["assembly_name"] not in cat:
                log(f"no CAT annotation for {a['assembly_name']}")
                continue
            jobs.append((r["sample"], a["haplotype"], a["assembly_name"], s3_to_https(a["assembly"]),
                         s3_to_https(a["assembly_fai"]), s3_to_https(a["assembly_gzi"]),
                         s3_to_https(cat[a["assembly_name"]]), queries))
    copy_cols = ["sample", "hap", "assembly", "family", "contig", "start", "end", "strand", "identity", "qcov",
                 "call", "exon78", "converted", "site_pattern", "anchors_on_contig"]
    hap_rows, n_fail = [], 0
    with open(os.path.join(T, "haplotype_copies.tsv"), "w") as fc, ThreadPoolExecutor(max_workers=6) as ex:
        fc.write("\t".join(copy_cols) + "\n")
        futs = {ex.submit(process, *j): j[2] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                rows, summ = fut.result()
                for row in rows:
                    fc.write("\t".join(map(str, row)) + "\n")
                fc.flush()
                hap_rows.append(summ)
            except Exception as e:  # noqa: BLE001
                n_fail += 1
                log(f"truth FAILED {futs[fut]}: {str(e)[:300]}")
            if i % 20 == 0:
                log(f"truth: {i}/{len(jobs)} haplotypes")

    hcols = ["sample", "hap", "assembly", "SMN1", "SMN2", "SMN_del7-8", "SMN?", "SMN_converted", "SMN_anchors",
             "SMN_flags", "HBA1", "HBA2", "HBA_hybrid", "HBA?", "HBA_total", "HBA_anchors", "HBA_flags"]
    with open(os.path.join(T, "haplotypes.tsv"), "w") as fh:
        fh.write("\t".join(hcols) + "\n")
        for r in sorted(hap_rows, key=lambda r: (r["sample"], r["hap"])):
            fh.write("\t".join(str(r.get(c, "")) for c in hcols) + "\n")

    by_s = defaultdict(list)
    for r in hap_rows:
        by_s[r["sample"]].append(r)
    num = ["SMN1", "SMN2", "SMN_del7-8", "SMN?", "SMN_converted", "HBA1", "HBA2", "HBA_hybrid", "HBA?", "HBA_total"]
    with open(os.path.join(T, "diploid.tsv"), "w") as fh:
        fh.write("\t".join(["sample", "n_haps"] + num + ["usable_SMN", "usable_HBA", "flags"]) + "\n")
        for s, rs in sorted(by_s.items()):
            sums = [sum(int(r.get(k, 0) or 0) for r in rs) for k in num]
            ok_smn = len(rs) == 2 and all("locus_not_found" not in r["SMN_flags"] for r in rs)
            ok_hba = len(rs) == 2 and all("locus_not_found" not in r["HBA_flags"] for r in rs)
            flags = ";".join(f"h{r['hap']}:{r['SMN_flags']}|{r['HBA_flags']}" for r in rs if r["SMN_flags"] or r["HBA_flags"])
            fh.write("\t".join(map(str, [s, len(rs)] + sums + [ok_smn, ok_hba, flags])) + "\n")
    state = "DONE" if n_fail == 0 else "PARTIAL"
    status("P8_truth", state, f"{len(hap_rows)}/{len(jobs)} haplotypes typed; {n_fail} failed; data/truth/diploid.tsv")


if __name__ == "__main__":
    main()
