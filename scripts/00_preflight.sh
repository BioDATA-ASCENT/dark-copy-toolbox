#!/usr/bin/env bash
# Phase 0 — tools, disk, connectivity and bandwidth. Writes logs/preflight.md.
# Exit 1 only for hard blockers (missing tools, too little disk). Unreachable hosts are reported, not fatal.
source "$(dirname "$0")/lib.sh"
status P0_preflight START
out="$LOGS/preflight.md"
: > "$out"

echo "## Tools" >> "$out"
missing=0
for t in samtools minimap2 CrossMap aws curl pigz python; do          # required tonight
  if command -v "$t" >/dev/null 2>&1; then
    v="$("$t" --version 2>&1 | head -1 || true)"; echo "- $t: ${v:-present}" >> "$out"
  else
    echo "- $t: **MISSING (required)**" >> "$out"; missing=1
  fi
done
for t in aria2c bcftools bedtools seqkit locityper ctyper parascopy smn_caller.py; do   # nice to have
  command -v "$t" >/dev/null 2>&1 && echo "- $t: present" >> "$out" || echo "- $t: missing (optional)" >> "$out"
done
python -c "import pysam, numpy, scipy, pandas, matplotlib, yaml" 2>>"$out" && echo "- python libs: ok" >> "$out" || missing=1

echo -e "\n## Disk" >> "$out"
gb="$(free_gb "$DATA")"; need="$(cfg min_free_gb)"
echo "- free under $DATA: ${gb} GB (need ${need})" >> "$out"

echo -e "\n## Connectivity" >> "$out"
check() { if url_ok "$2"; then echo "- OK   $1  $2" >> "$out"; else echo "- FAIL $1  $2" >> "$out"; fi; }
check kgp_pedigree        "$(cfg urls.kgp_pedigree)"
check kgp_seq_index_2504  "$(cfg urls.kgp_seq_index_2504)"
check kgp_seq_index_698   "$(cfg urls.kgp_seq_index_698)"
check kgp_population_ped  "$(cfg urls.kgp_population_ped)"
check kgp_s3              "$(cfg urls.kgp_s3_prefix)"
check grch38_fasta        "$(cfg urls.grch38_fasta)"
check grch38_fasta_s3     "$(cfg urls.grch38_fasta_s3)"
check chm13_fasta         "$(cfg urls.chm13_fasta_gz)"
check chm13_annotation    "$(cfg urls.chm13_annotation_s3_prefix)"
check chain               "$(cfg urls.chain_grch38_to_chm13)"
check gencode             "$(cfg urls.gencode_gtf)"
check ucsc_superdups      "$(cfg urls.grch38_superdups)"
check hprc_tables         "$(cfg urls.hprc_r2_assemblies)"
check zenodo_locityper    "$(cfg zenodo.locityper.base)/README.md?download=1"
check zenodo_ctyper       "https://zenodo.org/records/16340156"
check azul                "$(cfg urls.azul_service)/index/catalogs"
check dragen_v4_2_7       "$(cfg dragen_buckets.v4_2_7)"
check dragen_v4_4_7       "$(cfg dragen_buckets.v4_4_7)"
check ena_ftp_http        "http://ftp.sra.ebi.ac.uk/vol1/"

echo -e "\n## Bandwidth (first 200 MB of the GRCh38 FASTA)" >> "$out"
speed="$(curl -s -r 0-200000000 -o /dev/null -w '%{speed_download}' -m 120 "$(cfg urls.grch38_fasta)" || echo 0)"
echo "- EBI: $(python -c "print(round(${speed:-0}/1e6,1))") MB/s" >> "$out"
speed2="$(curl -s -r 0-200000000 -o /dev/null -w '%{speed_download}' -m 120 "$(s3_to_https "$(cfg urls.grch38_fasta_s3)")" || echo 0)"
echo "- S3 : $(python -c "print(round(${speed2:-0}/1e6,1))") MB/s" >> "$out"

cat "$out" >&2
if (( missing )); then status P0_preflight FAILED "missing tools — see logs/preflight.md"; exit 1; fi
if (( gb < need )); then status P0_preflight FAILED "only ${gb} GB free"; exit 1; fi
status P0_preflight DONE "see logs/preflight.md"
