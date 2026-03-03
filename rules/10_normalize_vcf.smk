rule normalize_mutect2:
    input:
        "results/10_mutect2_filtered/{sample}/{sample}.filtered.vcf.gz"
    output:
        "results/14_normalized/{sample}/{sample}.M2.normalized.vcf.gz",
        "results/14_normalized/{sample}/{sample}.M2.normalized.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        TMP1=$(mktemp -u).split.vcf.gz
        bcftools norm -m -both {input} -Oz -o $TMP1
        bcftools norm -f {REF_FASTA} $TMP1 -Oz -o {output[0]}
        tabix -p vcf {output[0]}
        rm -f $TMP1
        """


rule normalize_strelka2:
    input:
        "results/12_strelka/{sample}/run/results/variants/variants.vcf.gz"
    output:
        "results/14_normalized/{sample}/{sample}.S2.normalized.vcf.gz",
        "results/14_normalized/{sample}/{sample}.S2.normalized.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        TMP1=$(mktemp -u).split.vcf.gz
        bcftools norm -m -both {input} -Oz -o $TMP1
        bcftools norm -f {REF_FASTA} $TMP1 -Oz -o {output[0]}
        tabix -p vcf {output[0]}
        rm -f $TMP1
        """


rule normalize_deepvariant:
    input:
        "results/13_deepvariant/{sample}/{sample}.DV.vcf.gz"
    output:
        "results/14_normalized/{sample}/{sample}.DV.normalized.vcf.gz",
        "results/14_normalized/{sample}/{sample}.DV.normalized.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        TMP1=$(mktemp -u).split.vcf.gz
        bcftools norm -m -both {input} -Oz -o $TMP1
        bcftools norm -f {REF_FASTA} $TMP1 -Oz -o {output[0]}
        tabix -p vcf {output[0]}
        rm -f $TMP1
        """