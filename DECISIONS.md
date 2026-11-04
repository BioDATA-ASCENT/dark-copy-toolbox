# Decisions (overnight, times GMT; ET = GMT-5)

- 03:40 P0: EBI 4.9 MB/s vs S3 28 MB/s → 02_fetch_references.sh now tries the S3 GRCh38 FASTA first (EBI as fallback). GENCODE v49 preflight FAIL was transient (HTTP 200 on retry).
- 03:52 P6 smoke: 30 s/slice, ~1.05M reads, 15–33 MB per slice, 0 failures → raised samples.max_samples 150→650 and n_extra_random 0→200 (all 223 HPRC∩1kGP truth samples + trio parents + 200 random non-HPRC 1kGP samples for a population prior), slice_jobs 6→8. Estimated ≈ 2×~520×30 s/8 ≈ 65 min, ~25 GB.
