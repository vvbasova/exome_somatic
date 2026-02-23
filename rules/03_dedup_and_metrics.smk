# rules/03_dedup_metrics.smk

rule flagstat_align:
    input:
        bam="results/03_align/{sample}/{sample}.sorted.bam"
    output:
        txt="results/03_align_stats/{sample}/{sample}.flagstat.txt"
    log:
        "logs/03_align_stats/{sample}.flagstat.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/03_align_stats/{wildcards.sample} logs/03_align_stats
        samtools flagstat {input.bam} > {output.txt} 2> {log}
        """


rule idxstats_align:
    input:
        bam="results/03_align/{sample}/{sample}.sorted.bam"
    output:
        txt="results/03_align_stats/{sample}/{sample}.idxstats.txt"
    log:
        "logs/03_align_stats/{sample}.idxstats.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/03_align_stats/{wildcards.sample} logs/03_align_stats
        samtools idxstats {input.bam} > {output.txt} 2> {log}
        """


rule markduplicates:
    input:
        bam="results/03_align/{sample}/{sample}.sorted.bam"
    output:
        bam="results/04_dedup/{sample}/{sample}.dedup.bam",
        bai="results/04_dedup/{sample}/{sample}.dedup.bai",
        metrics="results/04_dedup/{sample}/{sample}.markdup.metrics.txt"
    log:
        "logs/04_dedup/{sample}.markduplicates.log"
    threads: 1
    resources:
        mem_mb=24000,
        runtime=480
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/04_dedup/{wildcards.sample} logs/04_dedup

        picard MarkDuplicates \
          I={input.bam} \
          O={output.bam} \
          M={output.metrics} \
          CREATE_INDEX=true \
          VALIDATION_STRINGENCY=SILENT \
          > {log} 2>&1
        """


rule flagstat_dedup:
    input:
        bam="results/04_dedup/{sample}/{sample}.dedup.bam"
    output:
        txt="results/04_dedup_stats/{sample}/{sample}.flagstat.txt"
    log:
        "logs/04_dedup_stats/{sample}.flagstat.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/04_dedup_stats/{wildcards.sample} logs/04_dedup_stats
        samtools flagstat {input.bam} > {output.txt} 2> {log}
        """


rule idxstats_dedup:
    input:
        bam="results/04_dedup/{sample}/{sample}.dedup.bam"
    output:
        txt="results/04_dedup_stats/{sample}/{sample}.idxstats.txt"
    log:
        "logs/04_dedup_stats/{sample}.idxstats.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/04_dedup_stats/{wildcards.sample} logs/04_dedup_stats
        samtools idxstats {input.bam} > {output.txt} 2> {log}
        """


rule collect_alignment_summary_metrics:
    input:
        ref=REF_FASTA,
        bam="results/04_dedup/{sample}/{sample}.dedup.bam"
    output:
        metrics="results/05_alignment_metrics/{sample}/{sample}.alignment_metrics.txt"
    log:
        "logs/05_alignment_metrics/{sample}.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/05_alignment_metrics/{wildcards.sample} logs/05_alignment_metrics

        picard CollectAlignmentSummaryMetrics \
          R={input.ref} \
          I={input.bam} \
          O={output.metrics} \
          VALIDATION_STRINGENCY=SILENT \
          > {log} 2>&1
        """


rule collect_hs_metrics:
    input:
        bam="results/04_dedup/{sample}/{sample}.dedup.bam",
        ref=REF_FASTA,
        baits=BAIT_INTERVALS,
        targets=TARGET_INTERVALS
    output:
        metrics="results/05_hs_metrics/{sample}/{sample}.hs_metrics.txt"
    log:
        "logs/05_hs_metrics/{sample}.log"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/05_hs_metrics/{wildcards.sample} logs/05_hs_metrics

        picard CollectHsMetrics \
          I={input.bam} \
          O={output.metrics} \
          R={input.ref}  \
          BAIT_INTERVALS={input.baits} \
          TARGET_INTERVALS={input.targets} \
          > {log} 2>&1
        """