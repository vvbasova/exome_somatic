# rules/00_common.smk

import os
import csv
import re

###############################################################################
# Samplesheet parsing
###############################################################################

REQUIRED_SAMPLE_COLUMNS = {"sample_id", "fastq_R1", "fastq_R2"}

def load_samples(path):
    """
    Returns:
      - samples: list of sample IDs
      - samples_dict: { sample_id: {"r1": path, "r2": path} }
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"No samplesheet: {path}")

    samples = []
    samples_dict = {}

    with open(path, newline="") as f:
        reader = csv.DictReader(f)

        if not reader.fieldnames:
            raise ValueError(f"Empty header in samplesheet: {path}")

        fields = set([c.strip() for c in reader.fieldnames])
        if fields != REQUIRED_SAMPLE_COLUMNS:
            raise ValueError(
                f"{path} must contain columns {sorted(REQUIRED_SAMPLE_COLUMNS)}, "
                f"but has {sorted(fields)}"
            )

        for row in reader:
            sample = row["sample_id"].strip()
            r1 = row["fastq_R1"].strip()
            r2 = row["fastq_R2"].strip()

            if not sample or not r1 or not r2:
                raise ValueError(f"Empty value in row: {row}")

            if sample in samples_dict:
                raise ValueError(f"Duplicate sample_id in samplesheet: {sample}")

            samples.append(sample)
            samples_dict[sample] = {"r1": r1, "r2": r2}

    if not samples:
        raise ValueError(f"No samples found in samplesheet: {path}")

    return samples, samples_dict


SAMPLESHEET = config["samplesheet"]
SAMPLES, SAMPLES_DICT = load_samples(SAMPLESHEET)

for s in SAMPLES:
    if not re.match(r"^[A-Za-z0-9_.-]+$", s):
        raise ValueError(f"Invalid sample name (non-ASCII or forbidden chars): {s}")

###############################################################################
# Config parsing
###############################################################################

REF_FASTA = config["reference"]["fasta"]
REF_FAI   = config["reference"].get("fai", "")
REF_DICT  = config["reference"].get("dict", "")

# GATK known sites (BQSR)
DBSNP = config["gatk"]["known_sites"]["dbsnp"]
MILLS = config["gatk"]["known_sites"]["mills"]

# Mutect2 resources
GNOMAD = config["mutect2"]["germline_resource"]
PON    = config["mutect2"]["panel_of_normals"]
COMMON = config["mutect2"]["common_variants"]

# VEP config
VEP_CACHE   = config.get("vep", {}).get("cache_dir", "")
VEP_SPECIES = config.get("vep", {}).get("species", "homo_sapiens")
VEP_ASSEMBLY = config.get("vep", {}).get("assembly", "GRCh38")

# Intervals
INTERVALS = config["intervals"]
TARGET_INTERVALS = INTERVALS["target"]
BAIT_INTERVALS = INTERVALS["bait"]

DEEPVARIANT_CONTAINER = config["deepvariant"]["container"]
###############################################################################
# Helper functions
###############################################################################

def get_r1(wc):
    return SAMPLES_DICT[wc.sample]["r1"]

def get_r2(wc):
    return SAMPLES_DICT[wc.sample]["r2"]


def rg_line(wc):
    """
    Build read group line for bwa mem -R.
    SM/ID are derived from sample; other fields come from config.readgroup.
    """
    rg = config.get("readgroup", {})
    pl = rg.get("PL", "ILLUMINA")
    lb = rg.get("LB", "lib1")
    pu = rg.get("PU", "unit1")
    cn = rg.get("CN", "YourCenter")

    return (
        f"@RG\\tID:{wc.sample}"
        f"\\tSM:{wc.sample}"
        f"\\tPL:{pl}"
        f"\\tLB:{lb}"
        f"\\tPU:{pu}"
        f"\\tCN:{cn}"
    )