#!/usr/bin/env python3
"""S8a — cut the ctyper k-mer matrix down to the gene groups we study (config ctyper.groups).

Own implementation (we never execute code from data/): stream-decompress HprcCpcHgsvc_cmr_matrix.txt.gz once,
copy the byte ranges listed in its .index, and write a matching .index (same header and columns, new offsets).
Output: data/ctyper/subset_matrix.txt(.index)
"""
import gzip
import os

from cfg import CFG, path
from common import log, status

SRC = path("precomputed", "ctyper", "HprcCpcHgsvc_cmr_matrix.txt")
OUT = path("ctyper")
os.makedirs(OUT, exist_ok=True)


def main():
    dst = os.path.join(OUT, "subset_matrix.txt")
    groups = set(CFG["ctyper"]["groups"])
    if os.path.exists(dst + ".ok") and open(dst + ".ok").read().split() == sorted(groups):
        return
    lines = open(SRC + ".index").read().splitlines()
    header = [l for l in lines if l.startswith("@")]
    rows = [l.split("\t") for l in lines if l.startswith("#") and l.split("\t")[0][1:] in groups]
    rows.sort(key=lambda r: int(r[1]))
    log(f"ctyper subset: {[r[0] for r in rows]}")
    out_rows, pos = [], 0
    with gzip.open(SRC + ".gz", "rb") as fh, open(dst, "wb") as out:
        cur = 0
        for r in rows:
            s, e = int(r[1]), int(r[2])
            while cur < s:                              # skip forward in big chunks
                n = min(64 << 20, s - cur)
                cur += len(fh.read(n))
            block = fh.read(e - s)
            cur += len(block)
            assert block.startswith(b"#" + r[0][1:].encode()), f"offset mismatch at {r[0]}"
            out.write(block)
            out_rows.append([r[0], str(pos), str(pos + len(block))] + r[3:])
            pos += len(block)
    with open(dst + ".index", "w") as fh:
        for h in header:
            fh.write(h + "\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")
    open(dst + ".ok", "w").write("\n".join(sorted(groups)) + "\n")
    status("S8_ctyper_subset", "DONE", f"{len(out_rows)} groups, {pos / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
