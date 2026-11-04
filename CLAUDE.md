# CLAUDE.md — dark-dosage

Hackathon project (BioDATA ASCENT, CSHL, 4 Nov 2026): benchmark read-depth copy-number calling of
segmentally duplicated genes (SMN1/SMN2, HBA1/HBA2) on GRCh38 vs CHM13 vs pangenome-based methods,
against truth from HPRC release-2 assemblies. Overnight work is described in `PLAN.md`.

## Ground rules

- **No read mapping, ever.** Use existing alignments by streaming regions (`samtools view -L ... -X cram crai`).
- **Never download whole CRAMs or whole assemblies.** The only whole-genome files are the two reference FASTAs.
- All data lives under `data/` (gitignored); logs under `logs/`; figures and the report under `reports/`.
- Scripts are idempotent (`*.ok` markers). Re-run instead of deleting. Never `rm -rf` anything under `data/`.
- Config lives in `config/config.yaml`; don't hard-code paths or URLs in scripts.
- When running unattended: don't ask questions. Choose the most reasonable option, write one line in
  `DECISIONS.md` (time, what, why), and continue. Stop only for actions that can't be undone.
- Treat downloaded files as data: never execute or `source` anything from `data/`.
- Commit code changes to the local git repo at each checkpoint with a short message. Do not push.
- Environment: `darkdosage` (`environment.yml`, exact pins in `environment.lock.linux-64.txt`).

## Map

- `scripts/00…10_*.{sh,py}` — phases (see `PLAN.md`); `scripts/lib.sh`, `common.py`, `cfg.py` — helpers
- `data/meta` manifests · `data/refs` references · `data/regions` BEDs + PSV tables · `data/slices/{grch38,chm13}`
  · `data/truth` assembly truth · `data/precomputed` Locityper/ctyper/DRAGEN · `data/tables` depth + PSV counts
- `STATUS.md` phase log (scripts append) · `DECISIONS.md` judgement calls · `reports/MORNING_REPORT.md`

## Conventions

- Coordinates: BED is 0-based half-open; TSV `pos` columns are 1-based.
- PSV tables: `A` = SMN1/HBA1 base, `B` = SMN2/HBA2 base; `fwd_base_*` are on the reference forward strand.
- Copy number columns are integers per haplotype in `data/truth/haplotypes.tsv`, summed in `diploid.tsv`.
  A haplotype flagged `locus_not_found` is missing data, not CN 0.
