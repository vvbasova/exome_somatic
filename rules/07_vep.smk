# rules/07_vep.smk

rule vep_annotate:
    input:
        vcf="results/10_mutect2_filtered/{sample}/{sample}.filtered.vcf.gz",
        fasta=REF_FASTA
    output:
        vcf="results/11_vep/{sample}/{sample}.vep.vcf",
        stats="results/11_vep/{sample}/{sample}.vep.stats.html"
    log:
        "logs/11_vep/{sample}.vep.log"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600
    conda:
        "../envs/vep.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/11_vep/{wildcards.sample} logs/11_vep

        vep \
          --offline \
          --cache \
          --dir_cache {VEP_CACHE} \
          --species {VEP_SPECIES} \
          --assembly {VEP_ASSEMBLY} \
          --fasta {input.fasta} \
          --vcf \
          -i {input.vcf} \
          -o {output.vcf} \
          --force_overwrite \
          --fork {threads} \
          --everything \
          --stats_file {output.stats} \
          > {log} 2>&1
        """


rule vep_annotate:
    input:
        vcf="results/10_mutect2_filtered/{sample}/{sample}.filtered.vcf.gz",
        fasta=REF_FASTA
    output:
        vcf="results/11_vep/{sample}/{sample}.vep.vcf",
        stats="results/11_vep/{sample}/{sample}.vep.stats.html"
    log:
        "logs/11_vep/{sample}.vep.log"
    threads: 16
    resources:
        mem_mb=64000,
        time="10:00:00"
    conda:
        "../envs/vep.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/11_vep/{wildcards.sample} logs/11_vep

        vep \
          --offline \
          --cache \
          --dir_cache {VEP_CACHE} \
          --species {VEP_SPECIES} \
          --assembly {VEP_ASSEMBLY} \
          --fasta {input.fasta} \
          --vcf \
          -i {input.vcf} \
          -o {output.vcf} \
          --force_overwrite \
          --fork {threads} \
          --everything \
          --stats_file {output.stats} \
          > {log} 2>&1
        """