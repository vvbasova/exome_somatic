# rules/05_mutect2_calling.smk

rule mutect2:
    input:
        ref=REF_FASTA,
        bam="results/07_bqsr_bam/{sample}/{sample}.recal.bam"
    output:
        vcf="results/08_mutect2_unfiltered/{sample}/{sample}.unfiltered.vcf.gz",
        f1r2="results/08_mutect2_unfiltered/{sample}/{sample}.f1r2.tar.gz"
    log:
        "logs/08_mutect2_unfiltered/{sample}.mutect2.log"
    threads: 4
    resources:
        mem_mb=24000,
        runtime=1440
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/08_mutect2_unfiltered/{wildcards.sample} logs/08_mutect2_unfiltered

        export GATK_JAVA_OPTS="-Xmx20g"

        gatk Mutect2 \
          -R {input.ref} \
          -I {input.bam} \
          --tumor-sample {wildcards.sample} \
          --germline-resource {GNOMAD} \
          --panel-of-normals {PON} \
          --f1r2-tar-gz {output.f1r2} \
          --native-pair-hmm-threads {threads} \
          -O {output.vcf} \
          > {log} 2>&1
        """


rule learn_read_orientation_model:
    input:
        f1r2="results/08_mutect2_unfiltered/{sample}/{sample}.f1r2.tar.gz"
    output:
        rom="results/09_mutect2_artifacts/{sample}/{sample}.read-orientation-model.tar.gz"
    log:
        "logs/09_mutect2_artifacts/{sample}.learnrom.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/09_mutect2_artifacts/{wildcards.sample} logs/09_mutect2_artifacts

        export GATK_JAVA_OPTS="-Xmx6g"

        gatk LearnReadOrientationModel \
          -I {input.f1r2} \
          -O {output.rom} \
          > {log} 2>&1
        """