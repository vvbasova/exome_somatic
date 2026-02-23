# rules/04_bqsr.smk

rule base_recalibrator:
    input:
        ref=REF_FASTA,
        bam="results/04_dedup/{sample}/{sample}.dedup.bam"
    output:
        table="results/06_bqsr_tables/{sample}/{sample}.recal.table"
    log:
        "logs/06_bqsr_tables/{sample}.baserecal.log"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=480
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/06_bqsr_tables/{wildcards.sample} logs/06_bqsr_tables

        gatk BaseRecalibrator \
          -R {input.ref} \
          -I {input.bam} \
          --known-sites {DBSNP} \
          --known-sites {MILLS} \
          -O {output.table} \
          > {log} 2>&1
        """


rule apply_bqsr:
    input:
        ref=REF_FASTA,
        bam="results/04_dedup/{sample}/{sample}.dedup.bam",
        table="results/06_bqsr_tables/{sample}/{sample}.recal.table"
    output:
        bam="results/07_bqsr_bam/{sample}/{sample}.recal.bam",
        bai="results/07_bqsr_bam/{sample}/{sample}.recal.bam.bai"
    log:
        "logs/07_bqsr_bam/{sample}.applybqsr.log"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=480
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/07_bqsr_bam/{wildcards.sample} logs/07_bqsr_bam

        gatk ApplyBQSR \
          -R {input.ref} \
          -I {input.bam} \
          --bqsr-recal-file {input.table} \
          -O {output.bam} \
          > {log} 2>&1

        samtools index {output.bam} 2>> {log}
        """