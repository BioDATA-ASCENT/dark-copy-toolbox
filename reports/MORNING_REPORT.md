# Morning report — dark dosage

## Phase status

| time | phase | state | note |
|---|---|---|---|
| 2026-11-04 03:38 | P0_preflight | START |  |
| 2026-11-04 03:40 | P0_preflight | DONE | see logs/preflight.md |
| 2026-11-04 03:40 | P1_metadata | START |  |
| 2026-11-04 03:40 | P2_references | START |  |
| 2026-11-04 03:43 | P2_references | DONE | GRCh38 + CHM13 + annotations + chain |
| 2026-11-04 03:46 | P1_metadata | DONE | 3202 S3 CRAMs; HPRC R2 tables ok |
| 2026-11-04 03:46 | P3_manifest | DONE | 3202 1kGP samples; selected 150 (150 truth, 0 parents) |
| 2026-11-04 03:48 | P4_chm13_crams | DONE | 3202 samples with CRAM+CRAI in ANVIL_T2T_CHRY |
| 2026-11-04 03:48 | P7_precomputed | START |  |
| 2026-11-04 03:48 | P5_regions | DONE | genes/loci/bins/controls/slice/psv for grch38 + chm13 |
| 2026-11-04 03:50 | P5_regions | DONE | genes/loci/bins/controls/slice/psv for grch38 + chm13 |
| 2026-11-04 03:50 | P6_ref_check | DONE | local FASTAs match CRAM M5 tags |
| 2026-11-04 03:51 | P6_slices | START | refs=grch38 chm13 limit=2 jobs=6 |
| 2026-11-04 03:51 | P6_slices | DONE | grch38 2, chm13 2 slices; 0 non-OK attempts (logs/slices.tsv) |
| 2026-11-04 03:51 | P3_manifest | DONE | 3202 1kGP samples; selected 650 (223 truth, 242 parents) |
| 2026-11-04 03:52 | P6_slices | START | refs=grch38 chm13 limit=0 jobs=8 |
| 2026-11-04 03:55 | P7_precomputed | START |  |
| 2026-11-04 03:56 | P7_precomputed | START |  |
| 2026-11-04 04:16 | P9_tables | DONE | 293 sample×ref depth tables; median control depth 35.819x |
| 2026-11-04 04:17 | P7_precomputed | DONE | Locityper, ctyper tables, DRAGEN per-sample outputs |
| 2026-11-04 04:19 | S1_dragen_table | DONE | 1950 sample×version rows -> data/tables/dragen_calls.tsv |
| 2026-11-04 04:20 | S2_locityper | DONE | loci ['CBS', 'HBM', 'MARVELD2']; LOO rows 360; 1kGP gts 1950; CN typed 0 (no SMN/HBA locus in db) |
| 2026-11-04 04:25 | S6_explorer | DONE | reports/explorer/data.json (2.3 MB) |
| 2026-11-04 04:29 | S6_explorer | DONE | reports/explorer/data.json (2.3 MB) |

## Data on disk

- selected samples: 650 ({'truth': 223, 'parent': 242, 'extra': 185})
- slices OK: GRCh38 276, CHM13 276, both 276
- non-OK slice attempts: GRCh38 0, CHM13 0 (see logs/slices.tsv)
- `refs`: 7.1G
- `slices`: 14G
- `precomputed`: 3.1G
- `truth`: 95M
- `tables`: 45M

## First-look figures

- `reports/fig_SMN_grch38_profile.png`
- `reports/fig_HBA_grch38_profile.png`
- `reports/fig_SMN_chm13_profile.png`
- `reports/fig_HBA_chm13_profile.png`


## Web explorer

- `reports/explorer/index.html` — open in a browser (self-contained, works offline): locus browser (GRCh38 vs CHM13, MAPQ≥0 vs ≥20, population band), per-sample card with truth / DRAGEN / posterior, live Bayes caller lab with sliders re-scored against truth, and the scoreboard.

## Locityper

- db loci overlapping our windows: CBS (CBS,U2AF1,CRYAA), HBM (HBZ,HBM), MARVELD2 (MARVELD2)
- **No Locityper locus contains SMN1/SMN2 or HBA1/HBA2** → no Locityper-implied CN for those genes (`data/tables/locityper_cn.tsv` says so). LOO benchmark rows: `data/tables/locityper_loo.tsv`; 1kGP genotypes (with `circular` flag): `data/tables/locityper_1kgp.tsv`.

## Precomputed

- Locityper loci covering our genes: db/loci/CBS, db/loci/HBM

- DRAGEN v4_0_3: 2604 per-sample files
- DRAGEN v4_2_7: 5438 per-sample files
- DRAGEN v4_4_7: 4788 per-sample files

## Decisions made overnight

# Decisions (overnight, times GMT; ET = GMT-5)

- 03:40 P0: EBI 4.9 MB/s vs S3 28 MB/s → 02_fetch_references.sh now tries the S3 GRCh38 FASTA first (EBI as fallback). GENCODE v49 preflight FAIL was transient (HTTP 200 on retry).
- 03:52 P6 smoke: 30 s/slice, ~1.05M reads, 15–33 MB per slice, 0 failures → raised samples.max_samples 150→650 and n_extra_random 0→200 (all 223 HPRC∩1kGP truth samples + trio parents + 200 random non-HPRC 1kGP samples for a population prior), slice_jobs 6→8. Estimated ≈ 2×~520×30 s/8 ≈ 65 min, ~25 GB.
- 03:57 P8 bug: HPRC R2 FASTAs use PanSN contig names (SAMPLE#HAP#CONTIG) while CAT GFF3 uses bare names → every faidx failed. Fix in 08_truth_from_assemblies.py: alias bare→PanSN names from the .fai. Re-ran.
- 03:50 P5: GENCODE v49 has no names for the chr21 false-dup copies (CBSL/U2AF1L5/CRYAA2) → located by aligning CBS/U2AF1/CRYAA to chr21 (config 'paralog_of'); KCNE1/KCNE1B dropped (KCNE1 is 9 Mb away, made a 9-Mb window). CHM13 RefSeq Liftoff v5.2 lacks CBS/U2AF1/CRYAA → lifted GRCh38 records via chain. Loci now get one window per gene cluster (slice BED 13.0 → 4.6 Mb).
- 04:10 P9: samtools bedcov on CRAM took ~3 min/slice (container re-decoded per region) → slice decoded once to a temp uncompressed BAM, then bedcov ×3 + PSV counts on it (~1 s/slice with 8 procs).
- 04:00 P7: DRAGEN fetch pulled HLA-folder metrics serially via aws-cli (~13 s/file) → patterns tightened (smn.tsv, hba, targeted, cnv_metrics, ploidy, wgs_coverage), /hla/ excluded, parallel curl (16) on path-style HTTPS (bucket 1000genomes-dragen-v4.0.3 has dots → virtual-host TLS fails).
- 04:00 Locityper DB (326 loci) has NO SMN1/SMN2 or HBA1/HBA2 locus (HBM locus ends at chr16:166,764, before HBA2). Locityper arm therefore only covers the chr21 CBS false-dup demo (locus CBS). our_loci.txt check in 07 looked one directory too shallow (archive unpacks to db/db/loci).
- 04:25 P8: fragmented assemblies produce fake copies/deletions (HG00232 hap2: contig ends inside SMN2 → looked like SMN2Δ7-8; hap1: locus split over 2 contigs). Added completeness rule: haplotype usable only if one contig holds all core genes plus ≥1 unique anchor gene on each side (config anchor_sides), and no copy within 5 kb of a contig end. Flags incomplete_locus / copy_at_contig_end → usable_SMN/HBA False. Restarted P8 (CAT/idx caches kept).


## Free disk

- 242 GB free under data/

