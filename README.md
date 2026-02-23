# Exome Processing Pipeline (Snakemake + Slurm)

Pipeline for processing whole-exome sequencing (WES) data:

**FASTQ → QC → Alignment → Deduplication → BQSR → Mutect2 → VEP annotation**

---

# 1. Installation

---

## Install Snakemake + Slurm executor 

```bash
conda install bioconda::snakemake
conda install bioconda::snakemake-executor-plugin-slurm
```



# 2. Project Structure

---

```
exome_processing/
│
├── config/
│   ├── config.yaml
│   └── samplesheet.csv
│
├── envs/
│   ├── bwa.yaml
│   ├── fastp.yaml
│   ├── fastqc.yaml
│   ├── gatk.yaml
│   └── picard.yaml
│
├── logs/
│   ├── 01_fastqc_raw/
│   ├── 02_fastp/
│   ├── 03_align/
│   ├── 03_align_stats/
│   ├── 04_dedup/
│   ├── 04_dedup_stats/
│   ├── 05_alignment_metrics/
│   ├── 05_hs_metrics/
│   ├── 06_bqsr_tables/
│   ├── 07_bqsr_bam/
│   ├── 08_mutect2_unfiltered/
│   ├── 09_mutect2_artifacts/
│   ├── 10_mutect2_filtered/
│   └── 11_vep/
│
├── results/
│   ├── 01_fastqc_raw/
│   ├── 02_fastp/
│   ├── 03_align/
│   ├── 03_align_stats/
│   ├── 04_dedup/
│   ├── 04_dedup_stats/
│   ├── 05_alignment_metrics/
│   ├── 05_hs_metrics/
│   ├── 06_bqsr_tables/
│   ├── 07_bqsr_bam/
│   ├── 08_mutect2_unfiltered/
│   ├── 09_mutect2_artifacts/
│   ├── 10_mutect2_filtered/
│   └── 11_vep/
│
├── rules/
│   ├── 00_common.smk
│   ├── 01_qc.smk
│   ├── 02_align.smk
│   ├── 03_dedup_and_metrics.smk
│   ├── 04_bqsr.smk
│   ├── 05_mutect2_calling.smk
│   ├── 06_mutect2_filtering.smk
│   └── 07_vep.smk
│
├── Snakefile
└── README.md
```

Each subdirectory in `results` and `logs` directories will have a subdirectory for each sample.

# 3. Configuration Files

---

## 3.1 samplesheet.csv

CSV format:

```
sample_id,fastq_R1,fastq_R2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

- Header must be exactly:
  - `sample_id`
  - `fastq_R1`
  - `fastq_R2`
- No duplicate sample IDs
- Full paths to FASTQ files
- All programs will use sample IDs in output files names

## 3.2 config.yaml 

The file should contain full paths to reference genome and
index files, VEP cache directory and GATK files.

Files for GATK utilities can be downloaded from GATK cloud storage:

https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38
https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0


# 4. Rules

---

### rules/00_common.smk

- Parses `samplesheet.csv` 
- Defines global variables

---

### rules/01_qc.smk

- `fastp_trim`: trimming + HTML/JSON report

---

### rules/02_align.smk

- `bwa_mem`: `bwa mem` + `samtools sort` + `samtools index`

---

### rules/03_dedup_and_metrics.smk

- `flagstat_align`, `idxstats_align`: samtools stats after alignment
- `markduplicates`: Picard `MarkDuplicates` + metrics + BAM index
- `flagstat_dedup`, `idxstats_dedup`: samtools stats after dedup
- `collect_alignment_summary_metrics`: Picard `CollectAlignmentSummaryMetrics` (post-dedup)
- `collect_hs_metrics`: Picard `CollectHsMetrics` (post-dedup)

---

### rules/04_bqsr.smk

- `base_recalibrator`: GATK `BaseRecalibrator` → `recal.table`
- `apply_bqsr`: GATK `ApplyBQSR` → recalibrated BAM + index

---

### rules/05_mutect2_calling.smk

- `mutect2`: GATK `Mutect2` → unfiltered VCF + `f1r2.tar.gz`
- `learn_read_orientation_model`: `LearnReadOrientationModel` → orientation model tar

---

### rules/06_mutect2_filtering.smk

- `get_pileup_summaries`: `GetPileupSummaries` → pileups table
- `calculate_contamination`: `CalculateContamination` → contamination table
- `filter_mutect_calls`: `FilterMutectCalls` → filtered VCF

---

### rules/07_vep.smk

- `vep_annotate`: VEP offline annotation 


# 5. Example run command 

---

```bash
snakemake \
  --executor slurm \
  --use-conda \
  --cores 32 \
  --jobs 20
```