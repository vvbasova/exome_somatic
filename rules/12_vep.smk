rule run_vep:
    input:
        vcf="results/15_merged/{sample}/{sample}.merged.vcf.gz",
        fasta=REF_FASTA
    output:
        vcf="results/16_vep/{sample}/{sample}.merged.vep.vcf",
        stats="results/16_vep/{sample}/{sample}.merged.vep.stats.html"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240
    conda:
        "../envs/vep.yaml"
    log:
        "logs/16_vep/{sample}.vep.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/16_vep/{wildcards.sample} logs/16_vep

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
