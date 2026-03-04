rule manta:
    input:
        bam="results/07_bqsr_bam/{sample}/{sample}.recal.bam",
        ref=REF_FASTA
    output:
        indels="results/11_manta/{sample}/results/variants/candidateSmallIndels.vcf.gz",
        indels_tbi="results/11_manta/{sample}/results/variants/candidateSmallIndels.vcf.gz.tbi",
        sv="results/11_manta/{sample}/results/variants/diploidSV.vcf.gz",
        sv_tbi="results/11_manta/{sample}/results/variants/diploidSV.vcf.gz.tbi"
    log:
        "logs/11_manta/{sample}.log"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=480
    conda:
        "../envs/strelka_manta.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=results/11_manta/{wildcards.sample}
        mkdir -p $OUTDIR

        configManta.py \
            --bam {input.bam} \
            --referenceFasta {input.ref} \
            --runDir $OUTDIR

        $OUTDIR/runWorkflow.py -m local -j {threads} > {log} 2>&1
        """