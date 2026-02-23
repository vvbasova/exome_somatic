# Snakefile

configfile: "config/config_01.yaml"


include: "rules/00_common.smk"
include: "rules/01_qc.smk"
include: "rules/02_align.smk"
include: "rules/03_dedup_and_metrics.smk"
include: "rules/04_bqsr.smk"
include: "rules/05_mutect2_calling.smk"
include: "rules/06_mutect2_filtering.smk"
include: "rules/07_vep.smk"


rule all:
    input:
        [
            expand("results/11_vep/{sample}/{sample}.vep.vcf", sample=SAMPLES),
            expand("results/11_vep/{sample}/{sample}.vep.stats.html", sample=SAMPLES),

            expand("results/01_fastqc_raw/{sample}/{sample}_R1_fastqc.html", sample=SAMPLES),
            expand("results/01_fastqc_raw/{sample}/{sample}_R2_fastqc.html", sample=SAMPLES),
            expand("results/02_fastp/{sample}/{sample}.fastp.html", sample=SAMPLES),

            expand("results/03_align_stats/{sample}/{sample}.flagstat.txt", sample=SAMPLES),
            expand("results/03_align_stats/{sample}/{sample}.idxstats.txt", sample=SAMPLES),

            expand("results/04_dedup/{sample}/{sample}.markdup.metrics.txt", sample=SAMPLES),
            expand("results/04_dedup_stats/{sample}/{sample}.flagstat.txt", sample=SAMPLES),
            expand("results/04_dedup_stats/{sample}/{sample}.idxstats.txt", sample=SAMPLES),

            expand("results/05_alignment_metrics/{sample}/{sample}.alignment_metrics.txt", sample=SAMPLES),
            expand("results/05_hs_metrics/{sample}/{sample}.hs_metrics.txt",sample=SAMPLES),
        ]