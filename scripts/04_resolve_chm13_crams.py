#!/usr/bin/env python3
"""Phase 4 — find CHM13v2.0 CRAM/CRAI for 1kGP samples (AnVIL_T2T_CHRY) via the AnVIL Azul API.

The AWS mirror stores objects under content hashes, so the .crai is NOT at <cram>.crai:
06_slice_crams.sh downloads each CRAI and passes it explicitly with `samtools view -X`.

Output: data/meta/chm13_crams.tsv  (sample, dataset, cram_name, cram_url, crai_url, cram_size)
        data/meta/azul_files_raw.jsonl  (every file hit, for debugging)
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request

from cfg import CFG, path
from common import SAMPLE_RE, log, s3_to_https, status

AZUL = CFG["urls"]["azul_service"].rstrip("/")
M = path("meta")


def get(url, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "dark-dosage-hackathon/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            log(f"azul GET failed ({i + 1}/{tries}): {e} :: {url[:160]}")
            time.sleep(10 * (i + 1))
    raise RuntimeError(f"Azul unreachable: {url}")


def q(endpoint, catalog, filters, size=1, extra=None):
    params = {"catalog": catalog, "filters": json.dumps(filters), "size": size}
    params.update(extra or {})
    return f"{AZUL}/index/{endpoint}?{urllib.parse.urlencode(params)}"


def facet_terms(resp, field):
    return [t["term"] for t in resp.get("termFacets", {}).get(field, {}).get("terms", []) if t.get("term")]


def main():
    cats = get(f"{AZUL}/index/catalogs")
    catalog = cats.get("default_catalog") or sorted(cats.get("catalogs", {}))[-1]
    log(f"Azul catalog: {catalog}")

    titles = facet_terms(get(q("datasets", catalog, {}, size=1)), "datasets.title")
    log(f"{len(titles)} dataset titles; T2T-like: {[t for t in titles if 'T2T' in t.upper()]}")

    def pick(cands):
        low = {t.lower(): t for t in titles}
        for c in cands:
            if c.lower() in low:
                return low[c.lower()]
        return None

    title = pick(CFG["anvil"]["dataset_title_candidates"])
    used_fallback = False
    if not title:
        title = pick(CFG["anvil"]["fallback_dataset_title_candidates"])
        used_fallback = bool(title)
    if not title:
        status("P4_chm13_crams", "FAILED", "no T2T dataset title found in Azul — see logs/run.log")
        raise SystemExit(1)
    log(f"dataset: {title}{' (FALLBACK: CHM13 v1.x alignments)' if used_fallback else ''}")

    base_filter = {"datasets.title": {"is": [title]}}
    formats = facet_terms(get(q("files", catalog, base_filter)), "files.file_format")
    want = [f for f in formats if re.search(r"cra[mi]$", f, re.I)]
    log(f"file formats in dataset: {formats} -> using {want}")
    filters = dict(base_filter, **{"files.file_format": {"is": want}})

    url = q("files", catalog, filters, size=250)
    raw = open(os.path.join(M, "azul_files_raw.jsonl"), "w")
    files = []
    while url:
        resp = get(url)
        for hit in resp.get("hits", []):
            for f in hit.get("files", []):
                raw.write(json.dumps(f) + "\n")
                files.append(f)
        url = resp.get("pagination", {}).get("next")
    raw.close()
    log(f"{len(files)} CRAM/CRAI file records")

    crams, crais = {}, {}
    for f in files:
        name = f.get("file_name") or ""
        mirror = f.get("azul_mirror_uri") or ""
        m = SAMPLE_RE.search(name)
        if not m or not mirror:
            continue
        rec = dict(sample=m.group(0), name=name, url=s3_to_https(mirror), size=f.get("file_size", ""))
        (crais if name.endswith(".crai") else crams).setdefault(m.group(0), []).append(rec)

    out = os.path.join(M, "chm13_crams.tsv")
    n_ok = 0
    with open(out, "w") as fh:
        fh.write("sample\tdataset\tcram_name\tcram_url\tcrai_url\tcram_size\n")
        for s in sorted(crams):
            # prefer names that look CHM13v2-ish, then the largest file
            cands = sorted(crams[s], key=lambda r: (("chm13" in r["name"].lower()) or ("v2" in r["name"].lower()),
                                                   int(r["size"] or 0)), reverse=True)
            c = cands[0]
            crai = next((x for x in crais.get(s, []) if x["name"] in (c["name"] + ".crai", c["name"][:-5] + ".crai")),
                        (crais.get(s) or [None])[0])
            if not crai:
                continue
            fh.write(f"{s}\t{title}\t{c['name']}\t{c['url']}\t{crai['url']}\t{c['size']}\n")
            n_ok += 1
    note = f"{n_ok} samples with CRAM+CRAI in {title}" + (" (FALLBACK v1.x)" if used_fallback else "")
    status("P4_chm13_crams", "DONE" if n_ok > 3000 else "PARTIAL", note)


if __name__ == "__main__":
    main()
