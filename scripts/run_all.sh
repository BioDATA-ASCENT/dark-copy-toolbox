#!/usr/bin/env bash
# Convenience orchestrator: runs every phase in dependency order, with the slow independent phases
# (slicing, truth, precomputed downloads) in parallel. Every phase is idempotent, so re-running resumes.
#   bash scripts/run_all.sh            # everything
#   bash scripts/run_all.sh --limit 20 # smoke run on the top-20 priority samples
source "$(dirname "$0")/lib.sh"
limit_args=("$@")
py() { python "scripts/$1"; }

bash scripts/00_preflight.sh
bash scripts/01_fetch_metadata.sh || log "P1 partial — continuing"
bash scripts/02_fetch_references.sh > "$LOGS/p2.log" 2>&1 & p2=$!
py 03_build_manifest.py
py 04_resolve_chm13_crams.py || log "P4 failed — CHM13 arm will be skipped by 06"
bash scripts/07_fetch_precomputed.sh > "$LOGS/p7.log" 2>&1 & p7=$!
wait "$p2" || log "P2 partial — check logs/p2.log"
py 05_make_regions.py
bash scripts/06_slice_crams.sh --check-ref || log "ref check failed — see logs/ref_check.tsv (CHM13 slices may fail)"
bash scripts/06_slice_crams.sh "${limit_args[@]}" > "$LOGS/p6.log" 2>&1 & p6=$!
python scripts/08_truth_from_assemblies.py > "$LOGS/p8.log" 2>&1 & p8=$!
wait "$p6"; wait "$p8" || log "P8 partial — check logs/p8.log"; wait "$p7" || log "P7 partial — check logs/p7.log"
py 09_depth_and_psv.py
py 10_morning_report.py
log "run_all finished"
