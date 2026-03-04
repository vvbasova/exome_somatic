rule make_metrics_report:
    input:
        fastp_json="results/02_fastp/{sample}/{sample}.fastp.json",
        alignment_metrics="results/05_alignment_metrics/{sample}/{sample}.alignment_metrics.txt",
        hs_metrics="results/05_hs_metrics/{sample}/{sample}.hs_metrics.txt"
    output:
        txt="results/17_final_tables/{sample}/{sample}.qc_report.txt"
    log:
        "logs/17_final_tables/{sample}.qc_report.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    conda:
        "../envs/parcer.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/17_final_tables/{wildcards.sample} logs/17_final_tables

        python scripts/make_metrics_report.py \
          --sample {wildcards.sample} \
          --fastp-json {input.fastp_json} \
          --alignment-metrics {input.alignment_metrics} \
          --hs-metrics {input.hs_metrics} \
          --outdir results/17_final_tables/{wildcards.sample} \
          --prefix {wildcards.sample} \
          > {log} 2>&1
        """


rule make_variant_tables:
    input:
        vcf="results/16_vep/{sample}/{sample}.merged.vep.vcf.gz",
        metrics_txt="results/17_final_tables/{sample}/{sample}.qc_report.txt",
        info_xlsx="data/info_list.xlsx",
        hotspots="data/hotspots_v2.xls",
    output:
        raw_csv="results/17_final_tables/{sample}/{sample}_raw.csv",
        filtered_csv="results/17_final_tables/{sample}/{sample}_filtered.csv",
        xlsx="results/17_final_tables/{sample}/{sample}_vep_annotated.xlsx",
    log:
        "logs/17_final_tables/{sample}.make_variant_tables.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/parcer.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/17_final_tables/{wildcards.sample} logs/17_final_tables

        python scripts/make_final_tables.py \
          {input.vcf} \
          --sample {wildcards.sample} \
          --outdir results/17_final_tables/{wildcards.sample} \
          --hotspots {input.hotspots} \
          --info-xlsx {input.info_xlsx} \
          --metrics-txt {input.metrics_txt} \
          > {log} 2>&1

        test -s {output.raw_csv}
        test -s {output.filtered_csv}
        test -s {output.xlsx}
        """