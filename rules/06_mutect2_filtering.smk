# rules/06_mutect2_filtering.smk

rule get_pileup_summaries:
    input:
        bam="results/07_bqsr_bam/{sample}/{sample}.recal.bam",
        common=COMMON,
        lims=TARGET_INTERVALS
    output:
        pileups="results/09_mutect2_artifacts/{sample}/{sample}.pileups.table"
    log:
        "logs/09_mutect2_artifacts/{sample}.pileups.log"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=360
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/09_mutect2_artifacts/{wildcards.sample} logs/09_mutect2_artifacts

        export GATK_JAVA_OPTS="-Xmx12g"

        gatk GetPileupSummaries \
          -I {input.bam} \
          -V {input.common} \
          -L {input.lims} \
          -O {output.pileups} \
          > {log} 2>&1
        """


rule calculate_contamination:
    input:
        pileups="results/09_mutect2_artifacts/{sample}/{sample}.pileups.table"
    output:
        contam="results/09_mutect2_artifacts/{sample}/{sample}.contamination.table"
    log:
        "logs/09_mutect2_artifacts/{sample}.contam.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=360
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail

        export GATK_JAVA_OPTS="-Xmx3g"

        gatk CalculateContamination \
          -I {input.pileups} \
          -O {output.contam} \
          > {log} 2>&1
        """


rule filter_mutect_calls:
    input:
        ref=REF_FASTA,
        vcf="results/08_mutect2_unfiltered/{sample}/{sample}.unfiltered.vcf.gz",
        contam="results/09_mutect2_artifacts/{sample}/{sample}.contamination.table",
        rom="results/09_mutect2_artifacts/{sample}/{sample}.read-orientation-model.tar.gz"
    output:
        vcf="results/10_mutect2_filtered/{sample}/{sample}.filtered.vcf.gz"
    log:
        "logs/10_mutect2_filtered/{sample}.filtermutectcalls.log"
    threads: 1
    resources:
        mem_mb=12000,
        runtime=360
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/10_mutect2_filtered/{wildcards.sample} logs/10_mutect2_filtered

        export GATK_JAVA_OPTS="-Xmx10g"

        gatk FilterMutectCalls \
          -V {input.vcf} \
          -R {input.ref} \
          --contamination-table {input.contam} \
          --ob-priors {input.rom} \
          -O {output.vcf} \
          > {log} 2>&1
        """
