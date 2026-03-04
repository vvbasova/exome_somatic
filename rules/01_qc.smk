# rules/01_qc.smk
rule fastp_trim:
    input:
        r1=get_r1,
        r2=get_r2
    output:
        r1="results/02_fastp/{sample}/{sample}_R1.trim.fq.gz",
        r2="results/02_fastp/{sample}/{sample}_R2.trim.fq.gz",
        html="results/02_fastp/{sample}/{sample}.fastp.html",
        json="results/02_fastp/{sample}/{sample}.fastp.json"
    log:
        "logs/02_fastp/{sample}.log"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    conda:
        "../envs/fastp.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/02_fastp/{wildcards.sample} logs/02_fastp

        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --thread {threads} \
          --html {output.html} \
          --json {output.json} \
          > {log} 2>&1
        """