#!/usr/bin/env bash
# Shared helpers for the overnight scripts: config access, logging, STATUS.md, retrying downloads.
# Source it:  source "$(dirname "$0")/lib.sh"
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"   # every script runs from the repo root
cfg() { python "$REPO_ROOT/scripts/cfg.py" "$1"; }

DATA="$REPO_ROOT/$(cfg data_dir)"
LOGS="$REPO_ROOT/$(cfg logs_dir)"
THREADS="$(cfg threads)"
mkdir -p "$DATA" "$LOGS"

export AWS_RETRY_MODE=standard AWS_MAX_ATTEMPTS=3   # fail fast on unreachable buckets

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS/run.log" >&2; }

# status PHASE STATE "message"   (STATE: START | DONE | PARTIAL | FAILED | SKIPPED)
status() {
  local f="$REPO_ROOT/STATUS.md"
  [[ -f "$f" ]] || printf '# Overnight status\n\n| time | phase | state | note |\n|---|---|---|---|\n' > "$f"
  printf '| %s | %s | %s | %s |\n' "$(date '+%F %H:%M')" "$1" "$2" "${3:-}" >> "$f"
  log "STATUS $1 $2 ${3:-}"
}

# s3://bucket/key  ->  https://bucket.s3.amazonaws.com/key   (anonymous HTTPS, works with htslib)
s3_to_https() { local u="${1#s3://}"; printf 'https://%s.s3.amazonaws.com/%s\n' "${u%%/*}" "${u#*/}"; }

url_ok() {  # HEAD (or 1-byte range) check
  local url="$1"
  if [[ "$url" == s3://* ]]; then timeout 90 aws s3 ls --no-sign-request --cli-connect-timeout 20 "$url" >/dev/null 2>&1; return; fi
  curl -sfIL --connect-timeout 15 -m 30 "$url" >/dev/null 2>&1 || curl -sfL --connect-timeout 15 -m 30 -r 0-0 -o /dev/null "$url" 2>/dev/null
}

# fetch URL DEST [MD5] — idempotent (DEST.ok marker), resumable, retried with backoff, md5-checked.
fetch() {
  local url="$1" dest="$2" md5="${3:-}" tmp i
  if [[ -s "$dest" && -e "$dest.ok" ]]; then log "have $(basename "$dest")"; return 0; fi
  mkdir -p "$(dirname "$dest")"; tmp="$dest.part"
  for i in 1 2 3 4 5; do
    if [[ "$url" == /* || "$url" == file://* ]]; then   # local file (tests, pre-staged data)
      cp -f "${url#file://}" "$tmp" && break
    elif [[ "$url" == s3://* ]]; then
      aws s3 cp --no-sign-request --only-show-errors --cli-connect-timeout 30 "$url" "$tmp" && break
    else
      if command -v aria2c >/dev/null 2>&1 && \
         aria2c -q -x4 -s4 -c --connect-timeout=30 --allow-overwrite=true --auto-file-renaming=false \
                -d "$(dirname "$tmp")" -o "$(basename "$tmp")" "$url"; then break; fi
      curl -fsSL --connect-timeout 30 --retry 3 -C - -o "$tmp" "$url" && break
    fi
    log "retry $i/5 for $url"; sleep $((i * 20))
  done
  if [[ ! -s "$tmp" ]]; then log "FAILED download $url"; return 1; fi
  if [[ -n "$md5" ]]; then
    local got; got="$(md5sum "$tmp" | cut -d' ' -f1)"
    if [[ "$got" != "$md5" ]]; then log "MD5 mismatch $(basename "$dest"): $got != $md5"; mv -f "$tmp" "$dest.badmd5"; return 1; fi
  fi
  mv -f "$tmp" "$dest"; rm -f "$dest.part.aria2"; touch "$dest.ok"
  log "got $(basename "$dest") ($(du -h "$dest" | cut -f1))"
}

# fetch_first DEST URL1 URL2 ...  — first URL that works wins
fetch_first() {
  local dest="$1"; shift
  local u
  for u in "$@"; do fetch "$u" "$dest" && return 0; done
  return 1
}

free_gb() { df -Pk "$1" | awk 'NR==2{printf "%d", $4/1024/1024}'; }
