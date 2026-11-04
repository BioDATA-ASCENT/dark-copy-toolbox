#!/usr/bin/env bash
# Phase 2 — reference FASTAs, annotations, chain, SD track. ~7 GB. Safe to run in the background.
source "$(dirname "$0")/lib.sh"
status P2_references START
R="$DATA/refs"; mkdir -p "$R"
rc=0

# GRCh38 exactly as used by the NYGC 1kGP CRAMs (needed to decode CRAM slices).
fetch_first "$R/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  "$(cfg urls.grch38_fasta)" "$(cfg urls.grch38_fasta_s3)" || rc=1
[[ -s "$R/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" ]] || samtools faidx "$R/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# CHM13v2.0. chr1-22/X are identical across the analysis-set variants, so this FASTA decodes
# autosomal/chrX slices of the AnVIL_T2T_CHRY CRAMs; 06_slice_crams.sh --check-ref confirms by MD5.
fetch "$(cfg urls.chm13_fasta_gz)" "$R/chm13v2.0.fa.gz" || rc=1
if [[ -s "$R/chm13v2.0.fa.gz" && ! -s "$R/chm13v2.0.fa.fai" ]]; then
  pigz -dc "$R/chm13v2.0.fa.gz" > "$R/chm13v2.0.fa.tmp" && mv "$R/chm13v2.0.fa.tmp" "$R/chm13v2.0.fa"
  samtools faidx "$R/chm13v2.0.fa"
fi
aws s3 ls --no-sign-request "$(cfg urls.chm13_analysis_set_s3_prefix)" > "$R/chm13_analysis_set.listing.txt" 2>/dev/null || true

# CHM13 gene annotation: newest RefSeq Liftoff GFF3 in the annotation folder.
aws s3 ls --no-sign-request "$(cfg urls.chm13_annotation_s3_prefix)" > "$R/chm13_annotation.listing.txt" 2>/dev/null || true
gff="$(awk '{print $4}' "$R/chm13_annotation.listing.txt" | grep -iE 'refseq.*liftoff.*\.gff3(\.gz)?$' | sort -V | tail -1 || true)"
[[ -n "$gff" ]] || gff="$(awk '{print $4}' "$R/chm13_annotation.listing.txt" | grep -iE '\.gff3(\.gz)?$' | grep -viE 'repeat|rm|trf|sat' | sort -V | tail -1 || true)"
if [[ -n "$gff" ]]; then
  if [[ "$gff" == *.gz ]]; then dest="$R/chm13v2.0.genes.gff3.gz"; else dest="$R/chm13v2.0.genes.gff3"; fi
  fetch "$(cfg urls.chm13_annotation_s3_prefix)$gff" "$dest" || rc=1
  echo "$gff" > "$R/chm13_annotation.chosen.txt"
else
  log "no CHM13 gene GFF3 found in the annotation listing"; rc=1
fi
# CHM13 segmental duplications (optional, for QC of control blocks)
sd="$(awk '{print $4}' "$R/chm13_annotation.listing.txt" | grep -iE '(sedef|segdup|_SD\.|SDs)' | grep -iE '\.bed(\.gz)?$' | head -1 || true)"
[[ -n "$sd" ]] && fetch "$(cfg urls.chm13_annotation_s3_prefix)$sd" "$R/chm13v2.0.segdups.bed${sd##*.bed}" || true

# GRCh38 annotation + SDs
gtf_url="$(cfg urls.gencode_gtf)"
fetch_first "$R/gencode.basic.annotation.gtf.gz" "$gtf_url" \
  "$(echo "$gtf_url" | sed 's/release_49/release_48/; s/v49/v48/')" \
  "$(echo "$gtf_url" | sed 's/release_49/release_47/; s/v49/v47/')" || rc=1
fetch "$(cfg urls.grch38_superdups)" "$R/hg38.genomicSuperDups.txt.gz" || true

fetch "$(cfg urls.chain_grch38_to_chm13)" "$R/grch38-chm13v2.chain" || rc=1

( cd "$R" && md5sum ./*.fa ./*.gz ./*.chain 2>/dev/null > md5sums.txt ) || true
if (( rc )); then status P2_references PARTIAL "see logs/run.log"; exit 1; fi
status P2_references DONE "GRCh38 + CHM13 + annotations + chain"
