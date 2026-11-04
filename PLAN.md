# Overnight plan — "Dark dosage" (BioDATA ASCENT code-a-thon, Wed 4 Nov 2026)

**You are Claude Code, running unattended overnight. Nobody will answer questions until ~08:00 ET.**
Read `CLAUDE.md` first (ground rules), then execute this plan top to bottom.

## Mission

Tomorrow (09:00–17:00 ET) the team asks: *how well do read-depth copy-number callers resolve genes in
segmental duplications (SMN1/SMN2, HBA1/HBA2), and does a T2T or pangenome reference beat linear GRCh38?*
They will **not** map reads. Tonight you prepare everything they need:

1. references + annotations,
2. a priority-ordered sample manifest (HPRC release-2 samples that are also in 1kGP, then their parents),
3. **region slices** of the existing 1kGP 30x CRAMs on GRCh38 *and* CHM13v2.0 (SMN locus, HBA cluster,
   a chr21 false-duplication demo locus, 30 contiguous normalisation blocks) — streamed, never whole CRAMs,
4. **truth copy numbers** from HPRC release-2 haplotype assemblies (sequence-based, per gene copy),
5. **precomputed pangenome/paralog-aware results** (Locityper, ctyper tables, DRAGEN v4.0.3/v4.2.7/v4.4.7),
6. depth + paralog-specific-site (PSV) allele-count tables, first-look figures, and `reports/MORNING_REPORT.md`.

## Arms the team will compare tomorrow (so you know why each download matters)

| Arm | Source | Prepared tonight by |
|---|---|---|
| Naive depth, GRCh38 (MAPQ≥0 vs ≥20) | GRCh38 slices | P6, P9 |
| Naive depth, CHM13 | CHM13 slices | P6, P9 |
| Paralog-aware linear (DRAGEN SMN/HBA callers, graph-hg38) | DRAGEN v4.0.3 / v4.2.7 per-sample outputs | P7 |
| Pangenome-reference mapping (DRAGEN multigenome) | DRAGEN v4.4.7 per-sample outputs | P7 |
| Pangenome haplotype genotyping | Locityper 1kGP genotypes + leave-one-out benchmark; ctyper allele tables | P7 |
| Team's Bayesian caller (depth + PSV allele fraction + population prior) | depth + PSV tables | P9 (inputs only) |
| **Truth** | HPRC R2 assemblies | P8 |

## Phases

All scripts are idempotent: re-running skips finished work (`*.ok` markers). Run from the repo root with
the `darkdosage` env active. Each script appends to `STATUS.md`.

| Phase | Command | Needs | Rough time | Done when |
|---|---|---|---|---|
| P0 | `bash scripts/00_preflight.sh` | – | 5 min | `logs/preflight.md` exists, exit 0 |
| P1 | `bash scripts/01_fetch_metadata.sh` | P0 | 5–15 min | `data/meta/kgp_s3_crams.txt` ≥ 3,000 CRAMs |
| P2 | `bash scripts/02_fetch_references.sh` **(background)** | P0 | 15–60 min | both FASTAs + `.fai`, GTF, CHM13 GFF3, chain |
| P3 | `python scripts/03_build_manifest.py` | P1 | 1 min | `data/meta/samples_selected.tsv` |
| P4 | `python scripts/04_resolve_chm13_crams.py` | P1 | 5 min | `data/meta/chm13_crams.tsv` ≥ 3,000 rows |
| P7 | `bash scripts/07_fetch_precomputed.sh` **(background)** | P3 | 30–90 min | Zenodo files md5-OK, DRAGEN keys listed |
| P5 | `python scripts/05_make_regions.py` | P2 | 5–10 min | `data/regions/*.slice.bed`, `*.psv.tsv` |
| P6a | `bash scripts/06_slice_crams.sh --check-ref` | P4, P5 | 2 min | all `match=yes` in `logs/ref_check.tsv` |
| P6 | `bash scripts/06_slice_crams.sh` **(background)** | P6a | 2–5 h | slices for both refs |
| P8 | `python scripts/08_truth_from_assemblies.py` **(background, parallel to P6)** | P5 | 1–3 h | `data/truth/diploid.tsv` |
| P9 | `python scripts/09_depth_and_psv.py` | P6 | 15–30 min | `data/tables/*.tsv.gz` |
| P10 | `python scripts/10_morning_report.py` | any | 2 min | `reports/MORNING_REPORT.md` |

`scripts/run_all.sh` chains these in the right order — useful for re-runs, but tonight **run phases one by one**
so you can check each "done when" and react.

### Running long phases

Start P2, P6, P7, P8 as background commands with their output redirected to `logs/p<N>.log`. Check in every
10–20 minutes (`tail` the log, `wc -l logs/slices.tsv`, `cat STATUS.md`). Don't sit in one blocking foreground
command for more than ~10 minutes. Background jobs survive between your turns; keep working on other phases
while they run.

### Smoke test before the bulk slice run

Before launching P6 in full, run `bash scripts/06_slice_crams.sh --limit 2` and check: two `.ok` slices per
reference, `samtools view -c` > 100k reads each, `logs/slices.tsv` seconds per slice. Then estimate the total
time: (seconds per slice × 2 refs × N samples) / `slice_jobs`. If it won't finish by ~06:00, lower
`samples.max_samples` in `config/config.yaml` (truth samples come first, so a shorter list stays balanced
across superpopulations), or raise `slice_jobs` if bandwidth is the limit and nothing is failing. Log the
change in `DECISIONS.md`.

## Acceptance checks (do these, don't just trust exit codes)

- **P1**: `head data/meta/kgp_s3_crams.txt` looks like `1000G_2504_high_coverage/.../HG00096.final.cram`.
- **P3**: about 200 HPRC-R2 samples in 1kGP; `role=truth` rows come first and alternate superpopulations.
- **P4**: spot-check one row: `samtools view -H "<cram_url>" | grep -m3 '^@SQ'` shows `chr1` with a CHM13 length
  (248,387,328). `dataset` must be `ANVIL_T2T_CHRY` (CHM13 v2.0). If the script fell back to `ANVIL_T2T`
  (CHM13 v1.x), say so prominently in the report.
- **P5**: `data/regions/projection_check.tsv` — the self-test must recover ≥ 90% of SMN_SNP_38 coordinates.
  `chm13.psv.tsv` should contain ≥ 2 SMN-like copies (CHM13 has 3 SMN-region repeat copies), and both
  references should have tens of HBA PSVs. GRCh38 `genes.bed` must contain CBSL/U2AF1L5/KCNE1B etc. (the false
  duplication); if GENCODE names differ, find them with `zgrep -w gene_name data/refs/gencode*.gtf.gz | grep chr21`.
- **P6a**: if CHM13 M5 mismatches, list `data/refs/chm13_analysis_set.listing.txt`, download the matching
  analysis-set FASTA (try `chm13v2.0_maskedY_rCRS.fa.gz` first), recheck, then update `ref_fasta.chm13` in the
  config. chr1-22 should be identical across variants, so a mismatch there means a wrong reference version.
- **P8**: sanity: per haplotype SMN1 + SMN2 is usually 2 (range 1–4); α-globin per haplotype usually 2
  (1 = −α3.7/−α4.2, 3 = triplication, 0 = −−SEA etc., mostly in EAS/SEA samples). `locus_not_found` haplotypes are
  **missing data, not CN 0**. Inspect 3 random haplotypes in `haplotype_copies.tsv` by hand.
- **P9**: control median depth ~25–40x per sample; on GRCh38 the SMN genes should show a MAPQ≥20 depth drop
  relative to MAPQ≥0 (that is the effect the team will measure).

## Fallbacks (decide and log in `DECISIONS.md`; don't stop)

- **EBI unreachable** → scripts already fall back to S3 mirrors; GRCh38 CRAMs come from S3 anyway.
- **S3 listing of 1kGP CRAMs fails** → 03 uses ENA paths from the sequence index (`http://ftp.sra.ebi.ac.uk/...`).
- **Azul API fails or no CHM13v2 CRAMs found** → retry later in the night (it may be transient). If it still
  fails: run P6 with `--ref grch38`, keep everything else, and write in the report exactly what the team must do
  in the morning (export a file manifest for AnVIL_T2T_CHRY from the AnVIL Data Explorer, whose
  `files.azul_mirror_uri` column gives S3 URIs, and save it as `data/meta/chm13_crams.tsv`). The CHM13 arm is then
  delayed, not lost.
- **A slice keeps failing** → it is logged and skipped. If > 10% fail for one reference, read
  `logs/slice_stderr.<ref>.log`, fix the cause (throttling → fewer `slice_jobs`; wrong URL scheme → convert;
  reference mismatch → P6a), then re-run (finished slices are skipped).
- **Truth haplotype fails** (CAT stream or remote faidx) → re-run P8 later; it skips finished CAT files. If a
  whole assembly region can't be read remotely, leave it as missing. **Don't** download whole assemblies.
- **Disk < 20 GB free** → stop starting new downloads, finish tables and the report.
- **A script has a bug** → make the minimal fix, re-run, note it in `DECISIONS.md` (file, what, why). Commit.

## Stretch goals (only once P0–P10 are DONE and it's before ~05:30)

Do them in this order; each must leave a tidy TSV in `data/tables/` and a line in the report.

1. **DRAGEN tidy table** → `data/tables/dragen_calls.tsv`: sample, dragen_version, SMN1_CN, SMN2_CN,
   SMN2_del78_CN, SMN total/full-length raw CN, HBA genotype + total copies, filters. Inspect a few files first;
   formats differ by version (`smn-tsv` in 4.0.3; targeted JSON in 4.2.7/4.4.7).
2. **Locityper tidy tables**: (a) which db loci contain SMN1/SMN2/HBA1/HBA2 (`our_loci.txt`); (b) for those
   loci, the leave-one-out benchmark rows from `locityper-benchmarking` for our truth samples; (c) the
   1kGP genotype (haplotype pair) for every selected sample. Then the key derived table: **Locityper-implied copy
   numbers**. Type the two predicted haplotypes' sequences (from the db's `haplotypes.fa.gz`) with
   `type_copies()` from `scripts/08_truth_from_assemblies.py` → `data/tables/locityper_cn.tsv`. Note that
   non-LOO 1kGP genotypes for HPRC samples are circular (their own haplotypes are in the panel): flag them.
3. **Cross-checks** → `data/tables/truth_crosscheck.tsv`: truth vs DRAGEN v4.2.7 SMN1/SMN2 on HPRC samples,
   truth vs ctyper `PangenomeAlleles_*` where the assembly names can be joined. List every disagreement. A few
   are expected (assembly errors, caller errors); many means a truth-typing bug, so investigate.
4. **Snakefile** wrapping P3→P10 (rules call the existing scripts) so the workflow is reproducible tomorrow.
5. **Bayesian caller v0** (`caller/bayes_cn.py`, clearly marked as a starting point): grid over (CN_A, CN_B)
   ∈ {0..5}², negative-binomial likelihood of gene-level depth (MAPQ≥0, normalised by controls),
   binomial likelihood of PSV allele counts with a small error rate, and an optional population prior from
   Locityper/DRAGEN 1kGP calls. Output posteriors per sample; evaluate against truth with leave-one-out
   (accuracy, Brier score, reliability plot). Keep it under ~200 lines; the team will refine it.

## Hard deadline

By **07:30 ET** `reports/MORNING_REPORT.md` must exist and be current. Re-run P10 at the end. At the top of
the report, add a 5-line human summary: what's ready, what isn't, the single most important caveat, and the
first thing the team should do at 09:00. Then commit and stop.

## What was verified while writing this plan (3 Nov, 21:00–22:00 ET)

- Env solves and every tool runs (Linux); pipeline tested end to end on synthetic data (slicing via explicit CRAI,
  idempotent re-runs, depth/PSV tables, report; PSV projection exact on both strands; copy typing of
  SMN1/SMN2/Δ7-8/converted copies).
- HPRC R2 tables (GitHub raw): 234 samples / 466 haplotypes with S3 FASTA + `.fai` + `.gzi`; CAT GFF3 for 462 haplotypes.
- SMNCopyNumberCaller GRCh38 files: 16 SMN1/SMN2 differentiating sites, 3,000 normalisation windows.
- Zenodo: Locityper record 15879519 and ctyper record 16340156 file lists + md5.
- DRAGEN 1kGP buckets and AnVIL Azul service exist; CHM13 1kGP alignments are in AnVIL_T2T_CHRY (AWS mirror).
- **Not reachable from the planning sandbox, so preflight checks them**: EBI FTP paths for the sequence
  indices / population file / GRCh38 FASTA, GENCODE v49 URL, the S3 key layout of the 1kGP CRAMs, the Azul field
  names (`datasets.title`, `files.file_format`, `azul_mirror_uri`), the CHM13 annotation file names.
