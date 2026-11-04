#!/usr/bin/env bash
# Phase 6 — stream ONLY the locus + control regions out of remote 1kGP CRAMs, for both references.
#
#   06_slice_crams.sh --check-ref        # verify local FASTAs match the CRAM @SQ M5 (run before bulk)
#   06_slice_crams.sh [--ref grch38|chm13|both] [--limit N]
#   06_slice_crams.sh --one SAMPLE REF   # internal worker (used via xargs)
#
# Output: data/slices/{ref}/{sample}.{ref}.{cram|bam} (+ index, + .ok); one line per attempt in logs/slices.tsv
source "$(dirname "$0")/lib.sh"

M="$DATA/meta"; G="$DATA/regions"; S="$DATA/slices"
declare -A REF=( [grch38]="$DATA/$(cfg ref_fasta.grch38)" [chm13]="$DATA/$(cfg ref_fasta.chm13)" )
FMT="$(cfg slice_format)"; EXT="$FMT"; [[ "$FMT" == cram ]] && FMTFLAG="-C" || FMTFLAG="-b"

# cram_url / crai_url for SAMPLE REF -> prints "cram<TAB>crai" or nothing
urls_for() {
  local s="$1" ref="$2"
  if [[ "$ref" == grch38 ]]; then
    awk -F'\t' -v s="$s" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
      $h["sample"]==s{c=$h["grch38_cram_s3"]; x=$h["grch38_crai_s3"]; if(c==""){c=$h["grch38_cram_ena"]; x=c".crai"} if(x=="")x=c".crai"; print c"\t"x}' \
      "$M/samples_all.tsv"
  else
    [[ -s "$M/chm13_crams.tsv" ]] && awk -F'\t' -v s="$s" '$1==s{print $4"\t"$5; exit}' "$M/chm13_crams.tsv"
  fi
  return 0
}

# first selected sample that has URLs for REF (empty if none) — never fails under set -e
first_with_urls() {
  local ref="$1" x
  while read -r x; do [[ -n "$(urls_for "$x" "$ref")" ]] && { echo "$x"; return 0; }; done \
    < <(awk -F'\t' 'NR>1{print $2}' "$M/samples_selected.tsv")
  return 0
}

seq_md5_local() { samtools faidx "$1" "$2" | grep -v '^>' | tr -d '\n' | tr 'a-z' 'A-Z' | md5sum | cut -d' ' -f1; }

check_ref() {
  local ref s line cram crai rc=0 out="$LOGS/ref_check.tsv"
  echo -e "ref\tsample\tcontig\tcram_M5\tlocal_md5\tmatch" > "$out"
  for ref in grch38 chm13; do
    s="$(first_with_urls "$ref")"
    [[ -z "$s" ]] && { log "check-ref: no $ref CRAM URL for any selected sample"; rc=1; continue; }
    line="$(urls_for "$s" "$ref")"; cram="${line%%$'\t'*}"
    hdr="$(timeout 300 samtools view -H "$cram" 2>>"$LOGS/run.log" || true)"
    [[ -z "$hdr" ]] && { log "check-ref: cannot read header of $cram"; rc=1; continue; }
    for contig in chr5 chr16 chr21 chr1; do
      m5="$(grep -P "^@SQ\tSN:${contig}\t" <<<"$hdr" | grep -oP 'M5:\K[0-9a-f]+' || true)"
      loc="$(seq_md5_local "${REF[$ref]}" "$contig")"
      if [[ -z "$m5" ]]; then match=no_M5_tag; elif [[ "$m5" == "$loc" ]]; then match=yes; else match=NO; rc=1; fi
      echo -e "$ref\t$s\t$contig\t${m5:-NA}\t$loc\t$match" >> "$out"
    done
  done
  cat "$out" >&2
  return $rc
}

one() {
  local s="$1" ref="$2" line cram crai dir out tmp t0 n
  dir="$S/$ref"; out="$dir/$s.$ref.$EXT"; mkdir -p "$dir" "$DATA/crai/$ref"
  [[ -e "$out.ok" ]] && return 0
  [[ -s "$LOGS/slices.tsv" ]] || echo -e "sample\tref\tstate\tseconds\treads\tbytes" > "$LOGS/slices.tsv"
  line="$(urls_for "$s" "$ref")"
  if [[ -z "$line" ]]; then echo -e "$s\t$ref\tNO_URL\t0\t0\t0" >> "$LOGS/slices.tsv"; return 0; fi
  cram="${line%%$'\t'*}"; crai="${line#*$'\t'}"
  fetch "$crai" "$DATA/crai/$ref/$s.crai" || { echo -e "$s\t$ref\tNO_CRAI\t0\t0\t0" >> "$LOGS/slices.tsv"; return 0; }
  tmp="$dir/.$s.$ref.tmp.$EXT"
  for attempt in 1 2 3; do
    t0=$SECONDS
    if timeout 5400 samtools view -@ 2 -T "${REF[$ref]}" -M -L "$G/$ref.slice.bed" $FMTFLAG \
         -o "$tmp" -X "$cram" "$DATA/crai/$ref/$s.crai" 2>>"$LOGS/slice_stderr.$ref.log" \
       && samtools index "$tmp" && samtools quickcheck "$tmp"; then
      n="$(samtools view -c -T "${REF[$ref]}" "$tmp")"
      mv -f "$tmp" "$out"; mv -f "$tmp".*i "$out.$( [[ $EXT == cram ]] && echo crai || echo bai )" 2>/dev/null || samtools index "$out"
      touch "$out.ok"
      echo -e "$s\t$ref\tOK\t$((SECONDS - t0))\t$n\t$(stat -c %s "$out")" >> "$LOGS/slices.tsv"
      return 0
    fi
    rm -f "$tmp" "$tmp".*i; sleep $((attempt * 60))
  done
  echo -e "$s\t$ref\tFAILED\t$((SECONDS - t0))\t0\t0" >> "$LOGS/slices.tsv"
  return 0
}

case "${1:-}" in
  --one) one "$2" "$3"; exit 0 ;;
  --check-ref)
    if check_ref; then status P6_ref_check DONE "local FASTAs match CRAM M5 tags"; else status P6_ref_check FAILED "see logs/ref_check.tsv"; exit 1; fi
    exit 0 ;;
esac

which_ref=both; limit=0
while [[ $# -gt 0 ]]; do case "$1" in --ref) which_ref="$2"; shift 2 ;; --limit) limit="$2"; shift 2 ;; *) echo "unknown arg $1"; exit 2 ;; esac; done
refs=(grch38 chm13); [[ "$which_ref" != both ]] && refs=("$which_ref")

# Reachability gate: if the first sample's index can't be fetched for a reference, skip that reference
# instead of burning ~5 min of retries per sample.
gated=()
for r in "${refs[@]}"; do
  first="$(first_with_urls "$r")"
  line="$( [[ -n "$first" ]] && urls_for "$first" "$r" || true)"
  if [[ -z "$line" ]]; then log "P6: no URLs for $r — skipping"; status P6_slices SKIPPED "$r: no CRAM URLs"; continue; fi
  if [[ "${line#*$'\t'}" == /* ]] || url_ok "${line#*$'\t'}"; then gated+=("$r")
  else log "P6: $r index unreachable (${line#*$'\t'}) — skipping $r"; status P6_slices SKIPPED "$r: CRAM host unreachable"; fi
done
refs=("${gated[@]}")
(( ${#refs[@]} )) || { status P6_slices FAILED "no reachable reference"; exit 1; }

status P6_slices START "refs=${refs[*]} limit=$limit jobs=$(cfg slice_jobs)"
[[ -f "$LOGS/slices.tsv" ]] || echo -e "sample\tref\tstate\tseconds\treads\tbytes" > "$LOGS/slices.tsv"
jobs="$LOGS/slice_jobs.txt"; : > "$jobs"
awk -F'\t' 'NR>1{print $2}' "$M/samples_selected.tsv" | { if (( limit > 0 )); then head -n "$limit"; else cat; fi; } \
  | while read -r s; do for r in "${refs[@]}"; do echo "$s $r"; done; done > "$jobs"
xargs -P "$(cfg slice_jobs)" -n 2 bash "$REPO_ROOT/scripts/06_slice_crams.sh" --one < "$jobs"

for r in "${refs[@]}"; do
  ok=$(ls "$S/$r"/*.ok 2>/dev/null | wc -l); tot=$(grep -c " $r$" "$jobs" || true)
  log "$r slices: $ok/$tot OK"
done
nfail=$(awk -F'\t' '$3!="OK" && NR>1' "$LOGS/slices.tsv" | wc -l)
status P6_slices DONE "grch38 $(ls "$S/grch38"/*.ok 2>/dev/null | wc -l), chm13 $(ls "$S/chm13"/*.ok 2>/dev/null | wc -l) slices; $nfail non-OK attempts (logs/slices.tsv)"
