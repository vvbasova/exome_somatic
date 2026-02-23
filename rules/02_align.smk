rule bwa_mem:
    input:
        r1="results/02_fastp/{sample}/{sample}_R1.trim.fq.gz",
        r2="results/02_fastp/{sample}/{sample}_R2.trim.fq.gz",
        ref=REF_FASTA
    output:
        bam="results/03_align/{sample}/{sample}.sorted.bam",
        bai="results/03_align/{sample}/{sample}.sorted.bam.bai"
    log:
        "logs/03_align/{sample}.log"
    threads: 8
    resources:
        mem_mb=24000,
        runtime=360
    conda:
        "../envs/bwa.yaml"
    params:
        rg=rg_line
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/03_align/{wildcards.sample} logs/03_align

        bwa mem -t {threads} -R '{params.rg}' \
          {input.ref} {input.r1} {input.r2} 2>> {log} \
        | samtools sort -@ {threads} -o {output.bam} - 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}
        """