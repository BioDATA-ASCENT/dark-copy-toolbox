# Overnight status

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
| 2026-11-04 04:30 | P10_report | DONE | reports/MORNING_REPORT.md |
| 2026-11-04 05:20 | P6_slices | DONE | grch38 648, chm13 650 slices; 2 non-OK attempts (logs/slices.tsv) |
| 2026-11-04 05:33 | P8_truth | DONE | 446/446 haplotypes typed; 0 failed; data/truth/diploid.tsv |
| 2026-11-04 05:55 | P9_tables | DONE | 1298 sample×ref depth tables; median control depth 34.733x |
| 2026-11-04 05:57 | S3_concordance | DONE | 47 method×target cells; 91 disagreements in truth_crosscheck.tsv |
| 2026-11-04 05:58 | S3_concordance | DONE | 43 method×target cells; 11 disagreements in truth_crosscheck.tsv |
| 2026-11-04 05:59 | S6_explorer | DONE | reports/explorer/data.json (3.6 MB) |
| 2026-11-04 06:00 | S6_explorer | DONE | reports/explorer/data.json (3.6 MB) |
| 2026-11-04 06:01 | S6_explorer | DONE | reports/explorer/data.json (3.6 MB) |
| 2026-11-04 06:02 | P10_report | DONE | reports/MORNING_REPORT.md |
| 2026-11-04 06:03 | P10_report | DONE | reports/MORNING_REPORT.md |
| 2026-11-04 12:28 | S7_pangraph_fa_SMN | DONE | 18 paths (35.7 Mb) |
| 2026-11-04 12:30 | S7_pangraph_fa_HBA | DONE | 50 paths (7.5 Mb) |
| 2026-11-04 12:31 | S7_pangraph_figs | DONE | SMN: 18 paths, HBA: 50 paths |
| 2026-11-04 12:37 | S7_pangraph_figs | DONE | SMN: 18 paths, HBA: 50 paths |
| 2026-11-04 12:39 | S6_explorer | DONE | reports/explorer/data.json (3.8 MB) |
| 2026-11-04 12:51 | S8_ctyper_subset | DONE | 2 groups, 1.4 MB |
| 2026-11-04 12:51 | S8_ctyper_subset | DONE | 4 groups, 2.9 MB |
| 2026-11-04 13:11 | S8_ctyper | DONE | 1298/1298 sample×ref runs; data/tables/ctyper_calls.tsv |
| 2026-11-04 13:12 | S3_concordance | DONE | 51 method×target cells; 11 disagreements in truth_crosscheck.tsv |
| 2026-11-04 13:13 | S6_explorer | DONE | reports/explorer/data.json (3.9 MB) |
| 2026-11-04 13:14 | P10_report | DONE | reports/MORNING_REPORT.md |
| 2026-11-04 14:07 | S6_explorer | DONE | reports/explorer/data.json (4.1 MB) |
| 2026-11-04 14:08 | S6_explorer | DONE | reports/explorer/data.json (4.1 MB) |
| 2026-11-04 14:19 | S3_concordance | DONE | 55 method×target cells; 11 disagreements in truth_crosscheck.tsv |
| 2026-11-04 14:20 | S6_explorer | DONE | reports/explorer/data.json (4.1 MB) |
| 2026-11-04 14:20 | S9_site | DONE | docs/index.html + docs/outline.html |
