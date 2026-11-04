# Decisions (overnight, times GMT; ET = GMT-5)

- 03:40 P0: EBI 4.9 MB/s vs S3 28 MB/s → 02_fetch_references.sh now tries the S3 GRCh38 FASTA first (EBI as fallback). GENCODE v49 preflight FAIL was transient (HTTP 200 on retry).
