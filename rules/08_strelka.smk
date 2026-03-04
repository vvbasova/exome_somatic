rule strelka:
    input:
        bam="results/07_bqsr_bam/{sample}/{sample}.recal.bam",
        ref=REF_FASTA,
        manta_indels="results/11_manta/{sample}/results/variants/candidateSmallIndels.vcf.gz"
    output:
        vcf="results/12_strelka/{sample}/run/results/variants/variants.vcf.gz",
        tbi="results/12_strelka/{sample}/run/results/variants/variants.vcf.gz.tbi"
    log:
        "logs/12_strelka/{sample}.log"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=480
    conda:
        "../envs/strelka_manta.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=results/12_strelka/{wildcards.sample}
        mkdir -p $OUTDIR

        configureStrelkaGermlineWorkflow.py \
            --bam {input.bam} \
            --referenceFasta {input.ref} \
            --indelCandidates {input.manta_indels} \
            --runDir $OUTDIR/run

        $OUTDIR/run/runWorkflow.py -m local -j {threads} > {log} 2>&1
        """