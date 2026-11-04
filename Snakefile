# Reproducible wrapper for P3 -> P10 + stretch steps. Every rule calls the existing idempotent script.
#   snakemake -j 4                 # everything (slicing is the slow part; it resumes)
#   snakemake -j 4 explorer        # rebuild the web explorer from existing tables
# Downloads (P1 metadata, P2 references, P7 precomputed) are run by the shell scripts and treated as inputs.
D = "data"

rule all:
    input: "reports/MORNING_REPORT.md", "reports/explorer/index.html"

rule manifest:
    input: f"{D}/meta/kgp_s3_crams.txt", f"{D}/meta/hprc_r2_assemblies.csv", "config/config.yaml"
    output: f"{D}/meta/samples_selected.tsv", f"{D}/meta/samples_all.tsv"
    shell: "python scripts/03_build_manifest.py"

rule chm13_crams:
    output: f"{D}/meta/chm13_crams.tsv"
    shell: "python scripts/04_resolve_chm13_crams.py"

rule regions:
    input: f"{D}/refs/grch38-chm13v2.chain", f"{D}/meta/SMN_SNP_38.txt", "config/config.yaml"
    output: f"{D}/regions/grch38.slice.bed", f"{D}/regions/chm13.slice.bed", f"{D}/regions/grch38.psv.tsv", f"{D}/regions/chm13.psv.tsv"
    shell: "python scripts/05_make_regions.py"

rule slices:
    input: rules.manifest.output, rules.chm13_crams.output, rules.regions.output
    output: touch("logs/slices.done")
    shell: "bash scripts/06_slice_crams.sh --check-ref && bash scripts/06_slice_crams.sh"

rule truth:
    input: rules.manifest.output, rules.regions.output
    output: f"{D}/truth/diploid.tsv", f"{D}/truth/haplotypes.tsv"
    shell: "python scripts/08_truth_from_assemblies.py"

rule tables:
    input: "logs/slices.done"
    output: f"{D}/tables/depth.grch38.tsv.gz", f"{D}/tables/psv_counts.grch38.tsv.gz", f"{D}/tables/control_depth.tsv"
    shell: "python scripts/09_depth_and_psv.py"

rule dragen:
    input: rules.manifest.output
    output: f"{D}/tables/dragen_calls.tsv"
    shell: "python scripts/11_dragen_table.py"

rule locityper:
    input: rules.manifest.output, rules.regions.output
    output: f"{D}/tables/locityper_our_loci.tsv", f"{D}/tables/locityper_loo.tsv", f"{D}/tables/locityper_1kgp.tsv"
    shell: "python scripts/12_locityper_tables.py"

rule bayes:
    input: rules.tables.output, rules.dragen.output, rules.truth.output
    output: f"{D}/tables/bayes_cn.grch38.mq0.tsv", f"{D}/tables/bayes_cn.chm13.mq0.tsv"
    shell: "python caller/bayes_cn.py --ref grch38 --eval && python caller/bayes_cn.py --ref chm13 --eval"

rule concordance:
    input: rules.bayes.output, rules.truth.output
    output: f"{D}/tables/concordance.tsv", f"{D}/tables/truth_crosscheck.tsv"
    shell: "python scripts/13_concordance.py"

rule explorer:
    input: rules.concordance.output, "explorer/template.html"
    output: "reports/explorer/index.html"
    shell: "python scripts/14_export_explorer.py"

rule report:
    input: rules.concordance.output, rules.locityper.output
    output: "reports/MORNING_REPORT.md"
    shell: "python scripts/10_morning_report.py"
