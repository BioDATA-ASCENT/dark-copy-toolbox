# Dark Copy Toolbox (dark-dosage)

**Live site:** [the toolbox](https://jacgonisa.github.io/dark-copy-toolbox/) (interactive, step by step) ·
[project outline](https://jacgonisa.github.io/dark-copy-toolbox/outline.html) (what we downloaded, how we processed it, results, problems).
Both pages are static files in [`docs/`](docs/), built by `scripts/14_export_explorer.py` + `scripts/19_build_site.py`.

**Does a T2T or pangenome reference fix copy-number calling in segmental duplications?**
BioDATA ASCENT code-a-thon, CSHL, 4 Nov 2026.

We benchmark read-depth copy-number calls for paralogous genes (SMN1/SMN2 for spinal muscular atrophy,
HBA1/HBA2 for α-thalassaemia) from 1000 Genomes 30x short reads aligned to GRCh38 vs T2T-CHM13, plus
pangenome-based genotypers, against truth from HPRC release-2 assemblies. **No reads are mapped**: we stream
small regions out of existing public alignments.

## 1. Set up the environment

Full environment (Linux or Intel Mac; used by the overnight run):

```bash
mamba env create -f environment.yml        # ~600 MB download
mamba activate darkdosage
# exact versions (Linux x86-64 only):  mamba create -n darkdosage --file environment.lock.linux-64.txt
```

Apple Silicon Mac: either build the Intel version under Rosetta
(`CONDA_SUBDIR=osx-64 mamba env create -f environment.yml`), or, if you only need to analyse the outputs
tomorrow, use the light environment:

```bash
mamba env create -f environment.laptop.yml && mamba activate darkdosage-lite
```

`environment.laptop.yml` solves on osx-arm64, osx-64 and linux-64. The full env solves on linux-64 and
osx-64; ctyper has no osx-arm64 build.

## 2. Overnight run (Claude Code)

Use a Linux machine or VM that stays awake, with ≥ 60 GB free disk. An AWS VM in us-east-1 is fastest
because most of the data is on S3 there, with free egress.

```bash
cd dark-dosage
mamba activate darkdosage
tmux new -s overnight
claude --permission-mode auto
```

Then tell it:

> Read PLAN.md and execute it end to end. I'm asleep — don't wait for me; log decisions in DECISIONS.md.

Watch the first two phases (about 10 minutes) before leaving it: preflight and metadata. They show whether the
data hosts are reachable from that machine. Headless alternative (Claude Code v2.1.259 or later):

```bash
claude -p "Read PLAN.md and execute it end to end, unattended." \
  --permission-mode auto --permission-prompts none > logs/claude_overnight.log 2>&1
```

If auto mode isn't available on your plan, run it interactively in `tmux` with `--permission-mode acceptEdits`
and approve Bash once per command pattern before bed. Use `--dangerously-skip-permissions` only inside a
throwaway container or VM.

No Claude Code? `bash scripts/run_all.sh` runs every phase in order (idempotent, re-runnable).

## 3. What you get in the morning

| Path | Contents |
|---|---|
| `reports/MORNING_REPORT.md` | what's ready, what failed, first-look figures |
| `data/meta/samples_selected.tsv` | priority-ordered samples (HPRC R2 ∩ 1kGP first, then their parents) |
| `data/slices/{grch38,chm13}/` | per-sample CRAM slices: SMN locus, HBA cluster, chr21 demo, 30 control blocks |
| `data/truth/diploid.tsv` | assembly-derived SMN1/SMN2/SMN2Δ7-8/HBA1/HBA2 copies per sample |
| `data/tables/depth.*.tsv.gz` | depth per bin/gene/control at MAPQ ≥ 0/1/20, normalised CN |
| `data/tables/psv_counts.*.tsv.gz` | allele counts at SMN1/SMN2 and HBA1/HBA2 paralog-specific sites |
| `data/precomputed/` | Locityper (1kGP genotypes, LOO benchmark, db), ctyper allele tables, DRAGEN per-sample outputs |

## 4. Plan for the day (4 h of hacking)

1. **09:00** Read the morning report; split into pairs: *linear vs T2T depth*, *precomputed callers vs truth*,
   *Bayesian caller*.
2. **Figure 1**: depth profile across SMN and HBA, GRCh38 vs CHM13, MAPQ ≥ 0 vs ≥ 20.
3. **Figure 2**: concordance with truth per method (naive GRCh38, naive CHM13, DRAGEN 4.2.7 graph-hg38,
   DRAGEN 4.4.7 multigenome, Locityper LOO).
4. **Figure 3**: Bayesian caller (depth + PSV allele fraction + population prior), calibration and reliability plot.

## 5. The tool: Dark Copy Toolbox

**Benchmark scope:** naive read depth (GRCh38 / CHM13 × MAPQ ≥ 0 / ≥ 20) and ctyper (pangenome k-mers, on GRCh38-
and CHM13-aligned reads), against HPRC R2 assembly truth, with RAD17 as a single-copy positive control.
DRAGEN is **not benchmarked** (its SMN/HBA callers are purpose-built for these exact loci); its public calls were
used only to cross-check the truth (outline §6a). Locityper is labelled "not run" (it maps reads).

- **`docs/index.html`** (built from `explorer/template.html` by `scripts/14_export_explorer.py` + `19_build_site.py`):
  six steps — the problem, reads on a reference (MAPQ slider, zoom, exon models), the pangenome, one person
  (truth vs ctyper), compare methods (scoreboard with 95% CIs), and alignment-free genotyping of simulated patients.
- **`docs/calls.html`** + **`docs/calls/*.tsv`** (`scripts/23_export_calls.py`): every call of every benchmarked
  method for all 650 people, unrounded naive CN, ctyper gene counts.
- **`docs/outline.html`** (`explorer/outline.html`): data sources with exact download commands, pipeline, results,
  problems, and the internal experiments (Bayes v0, DRAGEN).
- Internal, not in the deliverable: **`caller/bayes_cn.py`**, the Bayes v0 caller on the command line:
  `python caller/bayes_cn.py --ref grch38 --eval` (also `--ref chm13`, `--no-prior`, `--eff-sites`, `--conc`, `--eps`).
- **`scripts/13_concordance.py`** — every method vs truth on the same targets → `data/tables/concordance.tsv`.
- **`snakemake -j 4`** reruns everything from P3 on (all steps idempotent).

## Data sources

1000 Genomes 30x (NYGC, GRCh38) · AnVIL_T2T_CHRY (CHM13v2.0 alignments, AnVIL/AWS mirror) · HPRC release 2
assemblies + CAT annotation · T2T-CHM13v2.0 · GENCODE · Illumina DRAGEN 1kGP reanalysis (AWS Open Data) ·
Locityper (Zenodo 15879519) · ctyper (Zenodo 16340156) · SMNCopyNumberCaller site/region files.
