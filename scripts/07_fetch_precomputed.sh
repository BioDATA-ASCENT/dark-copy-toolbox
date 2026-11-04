#!/usr/bin/env bash
# Phase 7 — precomputed results: Locityper + ctyper (Zenodo, md5-checked) and DRAGEN per-sample
# targeted-caller outputs (SMN, HBA, ...) for the selected samples. Independent of slicing.
source "$(dirname "$0")/lib.sh"
status P7_precomputed START
P="$DATA/precomputed"; mkdir -p "$P/locityper" "$P/ctyper" "$P/dragen"
rc=0

zen() {  # zen TOOL KEY(files|optional_files)
  local tool="$1" key="$2" base name md5
  base="$(cfg zenodo.$tool.base)"
  while read -r name md5; do
    fetch "$base/$name?download=1" "$P/$tool/$name" "$md5" || rc=1
  done < <(python -c "
import sys; sys.path.insert(0, 'scripts')
from cfg import CFG
for k, v in (CFG['zenodo']['$tool'].get('$key') or {}).items(): print(k, v)")
}
zen locityper files
zen ctyper files
[[ "$(cfg fetch_ctyper_matrix)" == "True" ]] && zen ctyper optional_files

# Unpack Locityper benchmark + database (db is ~0.8 GB compressed)
L="$P/locityper"
for t in locityper-benchmarking db; do
  if [[ -s "$L/$t.tar.gz" && ! -e "$L/$t.unpacked" ]]; then
    mkdir -p "$L/$t" && tar -xzf "$L/$t.tar.gz" -C "$L/$t" && touch "$L/$t.unpacked"
  fi
done
# Which Locityper loci cover our genes?
if [[ -d "$L/db" ]]; then
  find "$L/db" -maxdepth 3 -type d | sed "s#$L/db/##" | sort > "$L/db_dirs.txt"
  grep -iE 'SMN|NAIP|SERF1|GTF2H2|HBA|HBZ|HBQ|HBM|CBS|U2AF1' "$L/db_dirs.txt" > "$L/our_loci.txt" || true
  log "Locityper loci matching our genes: $(tr '\n' ' ' < "$L/our_loci.txt")"
fi

# DRAGEN: list each bucket once, then pull small per-sample files for selected samples.
sel="$DATA/meta/samples_selected.tsv"
pat="$(python -c "import sys; sys.path.insert(0,'scripts'); from cfg import CFG; print('|'.join(CFG['dragen_file_patterns']))")"
maxb=$(( $(cfg dragen_max_file_mb) * 1024 * 1024 ))
for key in v4_0_3 v4_2_7 v4_4_7; do
  b="$(cfg dragen_buckets.$key)"; d="$P/dragen/$key"; mkdir -p "$d"
  if [[ ! -s "$d/listing.txt.gz" ]]; then
    log "listing $b (one-off, can take a few minutes)"
    aws s3 ls --no-sign-request --recursive "$b" | gzip > "$d/listing.txt.gz.tmp" && mv "$d/listing.txt.gz.tmp" "$d/listing.txt.gz" || { log "cannot list $b"; rc=1; continue; }
  fi
  if [[ -s "$sel" ]]; then
    awk -F'\t' 'NR>1{print $2}' "$sel" > "$d/samples.txt"
    zcat "$d/listing.txt.gz" | grep -F -f "$d/samples.txt" | grep -iE "$pat" \
      | grep -vE "$(cfg dragen_exclude_pattern)" \
      | awk -v m="$maxb" '$3<=m{print $4}' > "$d/wanted_keys.txt" || true
    log "DRAGEN $key: $(wc -l < "$d/wanted_keys.txt") files for selected samples"
    # parallel anonymous HTTPS fetches (aws-cli start-up made serial s3 cp ~13 s/file); skips files on disk
    # path-style URL in the bucket's own region (bucket names with dots break virtual-host TLS)
    bk="${b#s3://}"; bk="${bk%/}"
    rg="$(curl -sI "https://s3.amazonaws.com/$bk/" | tr -d '\r' | awk -F': ' 'tolower($1)=="x-amz-bucket-region"{print $2}')"
    bh="https://s3.${rg:-us-east-1}.amazonaws.com/$bk/"; mkdir -p "$d/files"
    grep -v -x -F -f <(cd "$d/files" && find . -type f -size +0 | sed 's#^\./##'; echo __none__) "$d/wanted_keys.txt" \
      | xargs -r -P 16 -I{} sh -c 'mkdir -p "$(dirname "$2/{}")" && curl -fsS --retry 3 -o "$2/{}" "$1{}" || echo "FAILED {}"' _ "$bh" "$d/files" || true
  fi
done

if (( rc )); then status P7_precomputed PARTIAL "see logs/run.log"; exit 1; fi
status P7_precomputed DONE "Locityper, ctyper tables, DRAGEN per-sample outputs"
