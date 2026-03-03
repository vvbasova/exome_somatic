rule merge_callers:
    input:
        m2="results/14_normalized/{sample}/{sample}.M2.normalized.vcf.gz",
        s2="results/14_normalized/{sample}/{sample}.S2.normalized.vcf.gz",
        dv="results/14_normalized/{sample}/{sample}.DV.normalized.vcf.gz"
    output:
        vcf="results/15_merged/{sample}/{sample}.merged.vcf.gz",
        tbi="results/15_merged/{sample}/{sample}.merged.vcf.gz.tbi"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/merge/{sample}.merge.log"
    shell:
        r"""
        mkdir -p results/15_merged/{wildcards.sample}

        python scripts/merge_M2_S2_DV_vcf.py \
            {input.m2} \
            {input.s2} \
            {input.dv} \
            {output.vcf} \
            {wildcards.sample} \
            --log {log}

        """