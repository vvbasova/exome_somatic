# rules/01_qc.smk

rule fastqc_raw:
    input:
        r1=get_r1,
        r2=get_r2
    output:
        html1="results/01_fastqc_raw/{sample}/{sample}_R1_fastqc.html",
        zip1 ="results/01_fastqc_raw/{sample}/{sample}_R1_fastqc.zip",
        html2="results/01_fastqc_raw/{sample}/{sample}_R2_fastqc.html",
        zip2 ="results/01_fastqc_raw/{sample}/{sample}_R2_fastqc.zip"
    log:
        "logs/01_fastqc_raw/{sample}.log"
    threads: 2
    resources:
        mem_mb=4000,
        runtime=30
    conda:
        "../envs/fastqc.yaml"
    shell:
        r"""
        set -euo pipefail

        outdir="results/01_fastqc_raw/{wildcards.sample}"
        mkdir -p "$outdir" logs/01_fastqc_raw

        fastqc -t {threads} -o "$outdir" "{input.r1}" "{input.r2}" >> "{log}" 2>&1

        r1base=$(basename "{input.r1}")
        r1base=${{r1base%.gz}}; r1base=${{r1base%.fastq}}; r1base=${{r1base%.fq}}

        r2base=$(basename "{input.r2}")
        r2base=${{r2base%.gz}}; r2base=${{r2base%.fastq}}; r2base=${{r2base%.fq}}

        mv "$outdir/${{r1base}}_fastqc.html" "{output.html1}"
        mv "$outdir/${{r1base}}_fastqc.zip"  "{output.zip1}"
        mv "$outdir/${{r2base}}_fastqc.html" "{output.html2}"
        mv "$outdir/${{r2base}}_fastqc.zip"  "{output.zip2}"
        """


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