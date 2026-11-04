#!/usr/bin/env bash
# Phase 1 — small metadata: 1kGP sample tables + CRAM listings, HPRC release-2 tables, SMN caller files.
source "$(dirname "$0")/lib.sh"
status P1_metadata START
M="$DATA/meta"; mkdir -p "$M"
rc=0

fetch "$(cfg urls.kgp_pedigree)" "$M/1kGP.3202_samples.pedigree_info.txt" || rc=1
fetch "$(cfg urls.kgp_population_ped)" "$M/20130606_g1k_3202_samples_ped_population.txt" \
  || log "population ped not found at the standard path; 03_build_manifest.py falls back to the sequence index POPULATION column"

s3p="$(cfg urls.kgp_s3_prefix)"
fetch_first "$M/1000G_2504_high_coverage.sequence.index" \
  "$(cfg urls.kgp_seq_index_2504)" "${s3p}1000G_2504_high_coverage.sequence.index" || rc=1
fetch_first "$M/1000G_698_related_high_coverage.sequence.index" \
  "$(cfg urls.kgp_seq_index_698)" "${s3p}additional_698_related/1000G_698_related_high_coverage.sequence.index" || rc=1

# Authoritative CRAM locations on S3 (anonymous). Narrow prefixes first, whole prefix as fallback.
lst="$M/kgp_s3_crams.txt"
if [[ ! -s "$lst" ]]; then
  { aws s3 ls --no-sign-request --recursive "${s3p}data/" || true
    aws s3 ls --no-sign-request --recursive "${s3p}additional_698_related/data/" || true; } \
    | awk '{print $4}' | grep -E '\.cram(\.crai)?$' > "$lst.tmp" || true
  if [[ "$(grep -c '\.cram$' "$lst.tmp" || true)" -lt 3000 ]]; then
    log "narrow S3 listing found $(grep -c '\.cram$' "$lst.tmp" || true) CRAMs; listing the whole prefix"
    aws s3 ls --no-sign-request --recursive "$s3p" | awk '{print $4}' | grep -E '\.cram(\.crai)?$' > "$lst.tmp" || true
  fi
  mv "$lst.tmp" "$lst"
fi
log "S3 CRAMs listed: $(grep -c '\.cram$' "$lst" || true)"

fetch "$(cfg urls.hprc_r2_assemblies)"     "$M/hprc_r2_assemblies.csv" || rc=1
fetch "$(cfg urls.hprc_r2_cat_annotation)" "$M/hprc_r2_cat_annotation.csv" || rc=1
fetch "$(cfg urls.hprc_r2_illumina)"       "$M/hprc_r2_illumina.csv" || true
fetch "$(cfg urls.smncaller_region_38)"    "$M/SMN_region_38.bed" || rc=1
fetch "$(cfg urls.smncaller_snp_38)"       "$M/SMN_SNP_38.txt" || rc=1

if (( rc )); then status P1_metadata PARTIAL "some metadata missing — see logs/run.log"; exit 1; fi
status P1_metadata DONE "$(grep -c '\.cram$' "$lst" || true) S3 CRAMs; HPRC R2 tables ok"
