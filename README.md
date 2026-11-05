# Dark Copy Toolbox

**Can a better reference genome count copies of near-identical genes?**
BioDATA ASCENT code-a-thon, Cold Spring Harbor Laboratory, November 2026.

👉 **[Open the toolbox](https://biodata-ascent.github.io/dark-copy-toolbox/)** ·
[Project outline](https://biodata-ascent.github.io/dark-copy-toolbox/outline.html) ·
[Every call, every person](https://biodata-ascent.github.io/dark-copy-toolbox/calls.html)

## Goal

Some disease genes come in near-identical pairs ("twins"). Short reads can't tell the twins apart, so standard
pipelines miscount them. Two examples:

- **SMN1 / SMN2**: having no working SMN1 causes spinal muscular atrophy (SMA); having one makes you a carrier.
- **HBA1 / HBA2**: losing α-globin copies causes α-thalassaemia.

We asked whether newer references fix this:

1. **A better linear reference.** Does the complete T2T-CHM13 genome give better copy-number calls than GRCh38?
2. **A pangenome.** Does ctyper, which genotypes against hundreds of real assembled haplotypes, do better than counting read depth?

The truth comes from HPRC release 2: the copies in each person's own assembled genome. The reads are public
1000 Genomes 30x short reads. We map no reads ourselves; we only stream the small regions we need out of
existing alignments.

## Deliverables

| | What it is |
|---|---|
| **[Toolbox](https://biodata-ascent.github.io/dark-copy-toolbox/)** | Interactive walk-through in six steps: the problem, reads on a reference, the pangenome, one person, the method scoreboard, and genotyping simulated patients. |
| **[Outline](https://biodata-ascent.github.io/dark-copy-toolbox/outline.html)** | Data sources, pipeline, results and caveats on one page. |
| **[All calls](https://biodata-ascent.github.io/dark-copy-toolbox/calls.html)** | Every call by every method for every person, with downloadable tables ([`docs/calls/`](docs/calls/)). |
| **Pipeline** | [`scripts/`](scripts/) and a [`Snakefile`](Snakefile) that rebuild every table and the site. |

## Key results

The benchmark covers 206 people for SMN and 221 for α-globin. Each cell is the share whose copy number exactly matches their assembly.

| Method | SMN1 | SMN total | α-globin total |
|---|---|---|---|
| **ctyper** (pangenome k-mers) | **91%** | 92% | **77%** |
| Naive depth, all reads (GRCh38 or CHM13) | 76% | **98%** | 3% |
| Naive depth, uniquely mapped reads only | 0% | – | 34–88% |
| Always guess "2 copies" | 77% | – | – |

- **Switching to T2T-CHM13 alone doesn't help.** On both references the twins share the same reads, so depth gets the twins' combined count right but not the split between them.
- **The pangenome does help.** ctyper is the only method that gets SMN1, the gene that decides SMA status, right for about 9 in 10 people.
- **The pipeline itself is sound.** A single-copy control gene (RAD17) comes out right for 223 of 223 people.
- **Simulated patients can be genotyped without aligning reads.** For an SMA patient, a carrier and an unaffected person, a small marker-counting method gets SMN1 right at all coverages from 30x down to 5x (and 50 of 60 at 2x). ctyper is right at 10x and above but fails increasingly below that.

## Team

- Jacob González Isa ([@jacgonisa](https://github.com/jacgonisa))
- Qingyao Sun ([@nalzok](https://github.com/nalzok))
- [@roberimun](https://github.com/roberimun)
- [@xzengComBio](https://github.com/xzengComBio)

## Run it

```bash
mamba env create -f environment.yml && mamba activate darkdosage
snakemake -j 4          # or: bash scripts/run_all.sh
```

## Data

The data comes from public sources:
- 1000 Genomes 30x reads aligned to GRCh38 (NYGC) and to T2T-CHM13 (AnVIL)
- HPRC release 2 assemblies
- the T2T-CHM13v2.0 reference
- GENCODE gene annotation
- the precomputed ctyper database (Zenodo 16340156)
