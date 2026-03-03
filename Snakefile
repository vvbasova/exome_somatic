# Snakefile
wildcard_constraints:
    sample = r"[A-Za-z0-9_.-]+"

configfile: "config/config_01.yaml"


include: "rules/00_common.smk"
include: "rules/01_qc.smk"
include: "rules/02_align.smk"
include: "rules/03_dedup_and_metrics.smk"
include: "rules/04_bqsr.smk"
include: "rules/05_mutect2_calling.smk"
include: "rules/06_mutect2_filtering.smk"
include: "rules/07_manta.smk"
include: "rules/08_strelka.smk"
include: "rules/09_deepvariant.smk"
include: "rules/10_normalize_vcf.smk"
include: "rules/11_merge_callers.smk"
include: "rules/12_vep.smk"

rule all:
    input:
        [
            expand("results/16_vep/{sample}/{sample}.merged.vep.vcf",sample=SAMPLES),

            expand("results/02_fastp/{sample}/{sample}.fastp.html", sample=SAMPLES),

            expand("results/03_align_stats/{sample}/{sample}.flagstat.txt", sample=SAMPLES),
            expand("results/03_align_stats/{sample}/{sample}.idxstats.txt", sample=SAMPLES),

            expand("results/04_dedup/{sample}/{sample}.markdup.metrics.txt", sample=SAMPLES),
            expand("results/04_dedup_stats/{sample}/{sample}.flagstat.txt", sample=SAMPLES),
            expand("results/04_dedup_stats/{sample}/{sample}.idxstats.txt", sample=SAMPLES),

            expand("results/05_alignment_metrics/{sample}/{sample}.alignment_metrics.txt", sample=SAMPLES),
            expand("results/05_hs_metrics/{sample}/{sample}.hs_metrics.txt",sample=SAMPLES),
        ]