# Decisions (overnight, times GMT; ET = GMT-5)

- 03:40 P0: EBI 4.9 MB/s vs S3 28 MB/s → 02_fetch_references.sh now tries the S3 GRCh38 FASTA first (EBI as fallback). GENCODE v49 preflight FAIL was transient (HTTP 200 on retry).
- 03:52 P6 smoke: 30 s/slice, ~1.05M reads, 15–33 MB per slice, 0 failures → raised samples.max_samples 150→650 and n_extra_random 0→200 (all 223 HPRC∩1kGP truth samples + trio parents + 200 random non-HPRC 1kGP samples for a population prior), slice_jobs 6→8. Estimated ≈ 2×~520×30 s/8 ≈ 65 min, ~25 GB.
- 03:57 P8 bug: HPRC R2 FASTAs use PanSN contig names (SAMPLE#HAP#CONTIG) while CAT GFF3 uses bare names → every faidx failed. Fix in 08_truth_from_assemblies.py: alias bare→PanSN names from the .fai. Re-ran.
- 03:50 P5: GENCODE v49 has no names for the chr21 false-dup copies (CBSL/U2AF1L5/CRYAA2) → located by aligning CBS/U2AF1/CRYAA to chr21 (config 'paralog_of'); KCNE1/KCNE1B dropped (KCNE1 is 9 Mb away, made a 9-Mb window). CHM13 RefSeq Liftoff v5.2 lacks CBS/U2AF1/CRYAA → lifted GRCh38 records via chain. Loci now get one window per gene cluster (slice BED 13.0 → 4.6 Mb).
