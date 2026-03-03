rule deepvariant:
    input:
        bam="results/07_bqsr_bam/{sample}/{sample}.recal.bam",
        bai="results/07_bqsr_bam/{sample}/{sample}.recal.bam.bai",
        ref=REF_FASTA
    output:
        vcf="results/13_deepvariant/{sample}/{sample}.DV.vcf.gz",
        gvcf="results/13_deepvariant/{sample}/{sample}.DV.g.vcf.gz"
    log:
        "logs/13_deepvariant/{sample}.log"
    threads: 16
    resources:
        mem_mb=32000,
        runtime=720
    shell:
        r"""
        set -euo pipefail

        OUTDIR=results/13_deepvariant/{wildcards.sample}
        mkdir -p $OUTDIR
        mkdir -p $OUTDIR/tmp
        mkdir -p $OUTDIR/logs

        singularity exec \
            --bind $PWD:/work \
            --bind /mnt/resources:/mnt/resources \
            --bind $OUTDIR/tmp:/tmp \
            --pwd /work \
            --env TMPDIR=/tmp \
            --env TEMP=/tmp \
            --env TMP=/tmp \
            {DEEPVARIANT_CONTAINER} \
            /opt/deepvariant/bin/run_deepvariant \
                --model_type=WES \
                --ref={input.ref} \
                --reads=/work/{input.bam} \
                --output_vcf=/work/{output.vcf} \
                --output_gvcf=/work/{output.gvcf} \
                --num_shards={threads} \
                --intermediate_results_dir=/tmp \
                --logging_dir=/work/$OUTDIR/logs \
                > {log} 2>&1
        """