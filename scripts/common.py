"""Shared helpers: population map, URL conversion, logging, PSV projection with minimap2."""
import datetime
import os
import re
import subprocess
import sys

from cfg import REPO_ROOT, logs

POP2SUPER = {
    **{p: "AFR" for p in ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"]},
    **{p: "AMR" for p in ["MXL", "PUR", "CLM", "PEL"]},
    **{p: "EAS" for p in ["CHB", "JPT", "CHS", "CDX", "KHV"]},
    **{p: "EUR" for p in ["CEU", "TSI", "FIN", "GBR", "IBS"]},
    **{p: "SAS" for p in ["GIH", "PJL", "BEB", "STU", "ITU"]},
}
SAMPLE_RE = re.compile(r"(?:HG|NA)\d{5}")


def log(msg):
    line = f"[{datetime.datetime.now():%F %T}] {msg}"
    print(line, file=sys.stderr, flush=True)
    os.makedirs(logs(), exist_ok=True)
    with open(logs("run.log"), "a") as fh:
        fh.write(line + "\n")


def status(phase, state, note=""):
    f = os.path.join(REPO_ROOT, "STATUS.md")
    if not os.path.exists(f):
        with open(f, "w") as fh:
            fh.write("# Overnight status\n\n| time | phase | state | note |\n|---|---|---|---|\n")
    with open(f, "a") as fh:
        fh.write(f"| {datetime.datetime.now():%F %H:%M} | {phase} | {state} | {note} |\n")
    log(f"STATUS {phase} {state} {note}")


def s3_to_https(uri):
    if not uri or not uri.startswith("s3://"):
        return uri
    bucket, _, key = uri[5:].partition("/")
    return f"https://{bucket}.s3.amazonaws.com/{key}"


def run(cmd, **kw):
    """Run a command (list), raise with stderr on failure."""
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))}\n{p.stderr[-2000:]}")
    return p.stdout


def read_fasta_region(fasta, region):
    """samtools faidx -> (name, sequence)."""
    out = run(["samtools", "faidx", fasta, region])
    lines = out.strip().splitlines()
    return lines[0][1:], "".join(lines[1:]).upper()


# ---------------------------------------------------------------- PSV projection
def _walk_cs(cs, qstart, tstart, strand, tlen):
    """Yield (query_pos0, target_pos0) for aligned (match/mismatch) columns of a PAF cs:Z (long or short).

    Query coordinates in PAF are on the forward query strand; for strand '-' the target coordinate runs
    backwards. Here qstart/tstart are PAF columns 3/8 (0-based starts) and tlen-based flip is handled
    by the caller passing the right tstart: we walk query forward and target forward on the aligned
    strand, then convert.
    """
    q, t = qstart, 0
    for op in re.findall(r"(:\d+|=[A-Za-z]+|\*[a-z]{2}|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})", cs):
        if op[0] == ":":
            n = int(op[1:])
            for i in range(n):
                yield q + i, t + i
            q += n
            t += n
        elif op[0] == "=":
            n = len(op) - 1
            for i in range(n):
                yield q + i, t + i
            q += n
            t += n
        elif op[0] == "*":
            yield q, t
            q += 1
            t += 1
        elif op[0] == "+":  # insertion in query
            q += len(op) - 1
        elif op[0] == "-":  # deletion from query (extra target bases)
            t += len(op) - 1
        elif op[0] == "~":  # intron (spliced alignment) — not expected here
            n = int(re.findall(r"\d+", op)[0])
            t += n


def project_positions(query_seq, target_fa, query_offsets, preset="asm20", extra=("-N", "20", "-p", "0.3")):
    """Align query_seq (string) to every sequence in target_fa and project 0-based query offsets.

    Returns a list of hits: dict(target, tstart, tend, strand, identity, qcov, sites={offset: (tpos0, tbase)}),
    where tpos0 is 0-based on the target forward strand and tbase is the target base read on the
    forward strand (complement it yourself if strand == '-' and you want query-strand bases).
    """
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(">query\n" + query_seq + "\n")
        qpath = fh.name
    try:
        paf = run(["minimap2", "-x", preset, "-c", "--cs=long", "--secondary=yes", *extra, target_fa, qpath])
    finally:
        os.unlink(qpath)
    want = set(query_offsets)
    hits = []
    tseqs = {}
    for line in paf.splitlines():
        f = line.split("\t")
        qlen, qs, qe, strand = int(f[1]), int(f[2]), int(f[3]), f[4]
        tname, tlen, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
        nmatch, alen = int(f[9]), int(f[10])
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None:
            continue
        sites = {}
        if strand == "+":
            for qpos, toff in _walk_cs(cs, qs, ts, strand, tlen):
                if qpos in want:
                    sites[qpos] = ts + toff
        else:
            # For '-' minimap2 reports cs relative to the reverse-complemented query aligned to the target
            # forward strand: walking target forward corresponds to query positions decreasing from qe-1.
            for qrc, toff in _walk_cs(cs, 0, ts, strand, tlen):
                qpos = qe - 1 - qrc
                if qpos in want:
                    sites[qpos] = ts + toff
        if tname not in tseqs:
            tseqs[tname] = read_fasta_region(target_fa, tname)[1]
        tseq = tseqs[tname]
        hits.append(dict(
            target=tname, tstart=ts, tend=te, strand=strand, qstart=qs, qend=qe,
            identity=nmatch / alen if alen else 0.0, qcov=(qe - qs) / qlen,
            sites={o: (p, tseq[p]) for o, p in sites.items()},
        ))
    return hits


COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]
