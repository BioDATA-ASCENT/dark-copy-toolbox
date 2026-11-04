#!/usr/bin/env python3
"""Stretch 1 — tidy DRAGEN SMN / HBA calls for the selected samples -> data/tables/dragen_calls.tsv

Formats differ by version:
  v4.0.3  <s>.smn.tsv only (no HBA caller)                              ref: hg38-alt-masked-v2-graph
  v4.2.7  <s>.final.smn.tsv + <s>.final.targeted.json (keys smn, hba)   ref: hg38_alt_masked_graph_v3
  v4.4.7  <s>.targeted.json (keys smn, hba with totalCopyNumber)        ref: hg38-alt_masked.cnv.graph...
HBA total copies come from `totalCopyNumber` when present, else by counting 'a' in the genotype
(aa/aa = 4, -a/aa = 3, --/aa = 2, aaa/aa = 5).
"""
import csv
import glob
import json
import os
import re

from cfg import path
from common import log, status

D = path("precomputed", "dragen")
COLS = ["sample", "dragen_version", "reference", "SMN1_CN", "SMN2_CN", "SMN2_del78_CN", "SMN_total_CN_raw",
        "SMN_fulllength_CN_raw", "SMN1_CN_raw_sites", "isSMA", "isCarrier", "HBA_genotype", "HBA_total_copies",
        "HBA_filter", "HBA_qual", "source_files"]


def hba_count(gt):
    if not gt or "/" not in gt:
        return ""
    return sum(h.count("a") for h in gt.split("/"))


def read_smn_tsv(fn):
    rows = list(csv.reader(open(fn), delimiter="\t"))
    hdr = [h.lstrip("#") for h in rows[0]]
    return dict(zip(hdr, rows[1])) if len(rows) > 1 else {}


def main():
    sel = {r["sample"] for r in csv.DictReader(open(path("meta", "samples_selected.tsv")), delimiter="\t")}
    out = []
    for key in ["v4_0_3", "v4_2_7", "v4_4_7"]:
        per = {}
        for fn in glob.glob(os.path.join(D, key, "files", "data", "individuals", "*", "*", "*")):
            sample = os.path.basename(os.path.dirname(fn))
            if sample not in sel:
                continue
            refname = fn.split(os.sep)[-3]
            rec = per.setdefault(sample, {"sample": sample, "dragen_version": key.replace("_", ".").lstrip("v"),
                                          "reference": refname, "source_files": []})
            if fn.endswith(".smn.tsv"):
                t = read_smn_tsv(fn)
                rec.update({"SMN1_CN": t.get("SMN1_CN", ""), "SMN2_CN": t.get("SMN2_CN", ""),
                            "SMN2_del78_CN": t.get("SMN2delta7-8_CN", ""), "SMN_total_CN_raw": t.get("Total_CN_raw", ""),
                            "SMN_fulllength_CN_raw": t.get("Full_length_CN_raw", ""),
                            "SMN1_CN_raw_sites": t.get("SMN1_CN_raw", ""), "isSMA": t.get("isSMA", ""),
                            "isCarrier": t.get("isCarrier", "")})
                rec["source_files"].append(os.path.basename(fn))
            elif fn.endswith("targeted.json"):
                j = json.load(open(fn))
                smn, hba = j.get("smn") or {}, j.get("hba") or {}
                if smn and not rec.get("SMN1_CN"):
                    rec.update({"SMN1_CN": smn.get("smn1CopyNumber", ""), "SMN2_CN": smn.get("smn2CopyNumber", ""),
                                "SMN2_del78_CN": smn.get("smn2Delta78CopyNumber", ""),
                                "SMN_total_CN_raw": smn.get("totalCopyNumberFloat", ""),
                                "SMN_fulllength_CN_raw": smn.get("fullLengthCopyNumberFloat", "")})
                if hba:
                    gt = hba.get("genotype", "")
                    rec.update({"HBA_genotype": gt,
                                "HBA_total_copies": hba.get("totalCopyNumber", hba_count(gt)),
                                "HBA_filter": hba.get("genotypeFilter", ""),
                                "HBA_qual": hba.get("genotypeQuality", hba.get("genotypeQual", ""))})
                rec["source_files"].append(os.path.basename(fn))
        log(f"DRAGEN {key}: {len(per)} samples with SMN/HBA outputs")
        out += per.values()
    fn = path("tables", "dragen_calls.tsv")
    with open(fn, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in sorted(out, key=lambda r: (r["sample"], r["dragen_version"])):
            r["source_files"] = ",".join(sorted(r["source_files"]))
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")
    status("S1_dragen_table", "DONE", f"{len(out)} sample×version rows -> data/tables/dragen_calls.tsv")


if __name__ == "__main__":
    main()
