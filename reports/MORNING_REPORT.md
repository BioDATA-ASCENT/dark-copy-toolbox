# Morning report — dark dosage

## Summary for the 09:00 team (written 2026-11-04 06:05 GMT)

1. **Ready:** P0–P10 and all 5 stretch goals ran end to end. That covers 650 1kGP samples sliced on GRCh38 (648) and CHM13 (650), and HPRC R2 truth for 223 samples (206 usable for SMN, 221 for α-globin). The interactive **Dark Dosage Explorer** is at `reports/explorer/index.html` (also published privately at https://claude.ai/artifact/GSTPPWSrGeTLuTCjyX4Kq1; share it from its Share menu before the team can open it).
2. **Results:** DRAGEN 4.x gets SMN1 right in 205/206 and α-globin total in 220/221. The Bayes v0 caller gets 94.7% SMN1 and 95.0% α-globin, or 94.2% / 86.9% without the prior. Naive MAPQ≥20 depth gets ≤0.5% of SMN1 calls right on both references. CHM13 fixes the chr21 CBS/U2AF1 false duplication: median CN goes from 0.27 at MAPQ≥20 on GRCh38 to 2.27 on CHM13.
3. **Not ready:** Locityper has no SMN or HBA locus, so it only covers the chr21 CBS demo. ctyper is not scored (bundled HBA allele, incomplete SMN coverage). The α-globin HBA1/HBA2 split is not scored, only the total, because −α3.7 hybrids are typed as a single copy.
4. **Most important caveat:** Bayes v0 is not independent of DRAGEN. Its prior (leave-one-out) and its SMN per-site PSV bias are both calibrated on DRAGEN 4.2.7 calls. Quote the `noprior` rows when claiming what read depth alone achieves. Truth itself agrees with DRAGEN except in 6 samples (`data/tables/truth_crosscheck.tsv`).
5. **First thing at 09:00:** open the explorer. Then look by hand at the 6 truth-vs-DRAGEN disagreements (HG00253, HG00597, HG00658, HG02668, NA18959, NA20905) in the sample card to decide whether truth or DRAGEN is wrong, before any accuracy number goes on a slide.

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

## Data on disk

- selected samples: 650 ({'truth': 223, 'parent': 242, 'extra': 185})
- slices OK: GRCh38 648, CHM13 650, both 648
- non-OK slice attempts: GRCh38 2, CHM13 0 (see logs/slices.tsv)
- `refs`: 7.1G
- `slices`: 31G
- `precomputed`: 3.1G
- `truth`: 325M
- `tables`: 196M

## Truth (HPRC R2 assemblies)

- samples typed: 223; usable SMN 206; usable HBA 221
- SMN1,SMN2 diploid copies: 1+1: 1, 1+4: 1, 2+0: 5, 2+1: 55, 2+2: 96, 2+3: 3, 3+0: 5, 3+1: 28, 3+2: 3, 4+0: 8, 4+1: 1
- α-globin total copies: 2: 2, 3: 23, 4: 189, 5: 6, 6: 1
- flags worth reading: `data/truth/diploid.tsv` column `flags`

## First-look figures

- `reports/fig_SMN_grch38_profile.png`
- `reports/fig_HBA_grch38_profile.png`
- `reports/fig_naive_vs_truth_SMN1_grch38.png`
- `reports/fig_naive_vs_truth_SMN2_grch38.png`
- `reports/fig_naive_vs_truth_HBA1_grch38.png`
- `reports/fig_naive_vs_truth_HBA2_grch38.png`
- `reports/fig_SMN_chm13_profile.png`
- `reports/fig_HBA_chm13_profile.png`
- `reports/fig_naive_vs_truth_SMN1_chm13.png`
- `reports/fig_naive_vs_truth_SMN2_chm13.png`
- `reports/fig_naive_vs_truth_HBA1_chm13.png`
- `reports/fig_naive_vs_truth_HBA2_chm13.png`

## Method concordance with assembly truth (exact match, usable truth samples)

| method | SMN1 | SMN2 incl. Δ7-8 | SMN total | α-globin total |
|---|---|---|---|---|
| bayes_chm13.mq0 | 95.1% (196/206) | 93.7% (193/206) | 97.6% (201/206) | 94.6% (209/221) |
| bayes_chm13.mq0.noprior | 93.7% (193/206) | 92.7% (191/206) | 98.1% (202/206) | 87.3% (193/221) |
| bayes_grch38.mq0 | 94.7% (195/206) | 93.2% (192/206) | 97.6% (201/206) | 95.0% (210/221) |
| bayes_grch38.mq0.noprior | 94.2% (194/206) | 92.7% (191/206) | 97.6% (201/206) | 86.9% (192/221) |
| dragen_4.0.3 | 99.5% (205/206) | 98.1% (202/206) | 97.6% (201/206) | – |
| dragen_4.2.7 | 99.5% (205/206) | 98.1% (202/206) | 97.6% (201/206) | 99.5% (220/221) |
| dragen_4.4.7 | 99.5% (205/206) | 98.1% (202/206) | 97.6% (201/206) | 99.5% (220/221) |
| naive_chm13_mq0 | 75.2% (155/206) | 71.8% (148/206) | 92.7% (191/206) | 16.7% (37/221) |
| naive_chm13_mq20 | 0.5% (1/206) | 9.7% (20/206) | 0.0% (0/206) | 64.7% (143/221) |
| naive_grch38_mq0 | 75.7% (156/206) | 57.3% (118/206) | 80.1% (165/206) | 5.9% (13/221) |
| naive_grch38_mq20 | 0.0% (0/206) | 8.3% (17/206) | 0.0% (0/206) | 56.1% (124/221) |

Full table: `data/tables/concordance.tsv`; per-call rows `data/tables/calls_long.tsv`; truth vs DRAGEN 4.2.7 / ctyper disagreements `data/tables/truth_crosscheck.tsv`.

- Bayes v0 `hm13.mq0.noprior`: HBA_total acc 0.873 brier 0.218 (n=221); SMN1 acc 0.937 brier 0.106 (n=206); SMN1+SMN2 acc 0.922 brier nan (n=206)
- Bayes v0 `hm13.mq0`: HBA_total acc 0.946 brier 0.072 (n=221); SMN1 acc 0.951 brier 0.085 (n=206); SMN1+SMN2 acc 0.932 brier nan (n=206)
- Bayes v0 `rch38.mq0.noprior`: HBA_total acc 0.869 brier 0.219 (n=221); SMN1 acc 0.942 brier 0.099 (n=206); SMN1+SMN2 acc 0.922 brier nan (n=206)
- Bayes v0 `rch38.mq0`: HBA_total acc 0.950 brier 0.072 (n=221); SMN1 acc 0.947 brier 0.082 (n=206); SMN1+SMN2 acc 0.927 brier nan (n=206)

## Web explorer

- `reports/explorer/index.html` — open in a browser (self-contained, works offline): locus browser (GRCh38 vs CHM13, MAPQ≥0 vs ≥20, population band), per-sample card with truth / DRAGEN / posterior, live Bayes caller lab with sliders re-scored against truth, and the scoreboard.

## Truth sanity (per complete haplotype)

- haplotypes typed 446; complete SMN locus 419; complete HBA locus 443
- SMN copies (SMN1+SMN2+Δ7-8) per haplotype: 1: 69, 2: 340, 3: 10
- α-globin copies per haplotype: 1: 27, 2: 409, 3: 6, 4: 1

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
- 2026-11-04 05:35 | P6 | HG00619/HG00620 GRCh38 have no CRAI in the S3 index → 648/650 GRCh38 slices; accepted (<1% loss), no fallback attempted.
- 2026-11-04 05:35 | P8 | Truth sanity: 206/223 samples usable for SMN, 221/223 for HBA; HBA totals 4 (189), 3 (23), 5 (6), 2 (2), 6 (1); −α3.7 haplotypes appear as a single HBA1-like copy (hybrid), so HBA1/HBA2 split is not scored, only HBA_total.
- 2026-11-04 06:00 | S3 | ctyper PangenomeAlleles annotation dropped from scoring: HBA alleles are a bundled HBA2,HBA1,HBAP1 group and SMN alleles cover only R1 haplotypes that resolved the locus (50% 'accuracy' = coverage artefact). Counts kept in data/tables/ctyper_annotation_counts.tsv.


## Free disk

- 223 GB free under data/

