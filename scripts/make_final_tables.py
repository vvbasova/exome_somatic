#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam


###############################################################################
# config
###############################################################################

DEFAULT_MIN_DP = 30
DEFAULT_MIN_VAF = 0.035
DEFAULT_MAX_POP_AF = 0.05

DEFAULT_M2_BAD_FILTERS = {"contamination", "map_qual", "base_qual"}
DEFAULT_S2_BAD_FILTERS = {"LowDepth", "LowGQX", "HighDepth", "HighDPFRatio", "SiteConflict"}
DEFAULT_DV_BAD_FILTERS = {"RefCall"}

DEFAULT_KEEP_IF_CLINSIG_HAS = ("pathogenic", "likely_pathogenic")
DEFAULT_HOTSPOT_COL = "is_in_hotspot_position"

DEFAULT_KEEP_IMPACT = {"HIGH", "MODERATE"}
DEFAULT_DROP_CONSEQUENCE = {
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intergenic_variant",
}

EXCEL_MAX_ROWS = 1_000_000


###############################################################################
# logger
###############################################################################

def setup_logger(verbosity: int) -> logging.Logger:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    return logging.getLogger("make_reports")


###############################################################################
# vcf parsing
###############################################################################

def flatten_pysam_value(v: Any) -> Optional[str]:
    if v is None:
        return None
    if isinstance(v, (str, int, float)):
        return str(v)
    if v is True:
        return "True"
    if v is False:
        return "False"
    if isinstance(v, (tuple, list)):
        return ",".join(map(str, v))
    return str(v)


def csq_fields_from_header(vf: pysam.VariantFile, csq_tag: str = "CSQ") -> List[str]:
    if csq_tag not in vf.header.info:
        return []
    desc = vf.header.info[csq_tag].description or ""
    if "Format:" not in desc:
        return []
    fmt = desc.split("Format:", 1)[1].strip()
    return fmt.split("|") if fmt else []


def parse_vep_union_vcf(
    vcf_path: Path,
    sample_name: Optional[str] = None,
    sample_index: int = 0,
    csq_tag: str = "CSQ",
    explode_csq: bool = True,
    keep_all_info: bool = True,
    keep_all_format: bool = True,
) -> pd.DataFrame:
    vf = pysam.VariantFile(str(vcf_path))

    samples = list(vf.header.samples)
    chosen_sample = None
    if samples:
        if sample_name is not None:
            if sample_name not in samples:
                raise ValueError(f"Sample '{sample_name}' not found in VCF samples: {samples}")
            chosen_sample = sample_name
        else:
            if sample_index >= len(samples):
                raise ValueError(f"sample_index={sample_index} out of range; samples={samples}")
            chosen_sample = samples[sample_index]

    csq_fields = csq_fields_from_header(vf, csq_tag=csq_tag)
    n_csq = len(csq_fields)

    rows: List[Dict[str, Any]] = []
    append = rows.append
    flat = flatten_pysam_value

    for rec in vf.fetch():
        contig = rec.contig
        pos = rec.pos
        rid = rec.id if rec.id is not None else "."
        ref = rec.ref

        alts = list(rec.alts) if rec.alts else []
        alt0 = alts[0] if alts else None
        alt_list = ",".join(alts) if alts else None

        qual = rec.qual
        filt = "PASS" if len(rec.filter.keys()) == 0 else ";".join(list(rec.filter.keys()))

        base: Dict[str, Any] = {
            "CHROM": contig,
            "POS": pos,
            "ID": rid,
            "REF": ref,
            "ALT": alt0,
            "ALT_LIST": alt_list,
            "QUAL": qual,
            "FILTER": filt,
            "SAMPLE_NAME": chosen_sample,
        }

        if keep_all_info:
            for k, v in rec.info.items():
                if k == csq_tag:
                    continue
                base[f"INFO_{k}"] = flat(v)

        if chosen_sample is not None and chosen_sample in rec.samples and keep_all_format:
            s = rec.samples[chosen_sample]
            for k, v in s.items():
                base[f"FORMAT_{k}"] = flat(v)

        csq_val = rec.info.get(csq_tag, None)

        if not csq_val:
            if explode_csq:
                append(dict(base))
            else:
                out = dict(base)
                out[f"INFO_{csq_tag}"] = None
                append(out)
            continue

        if not explode_csq:
            out = dict(base)
            out[f"INFO_{csq_tag}"] = flat(csq_val)
            append(out)
            continue

        for csq in csq_val:
            parts = csq.split("|")
            out = dict(base)

            if n_csq:
                lim = min(len(parts), n_csq)
                for i in range(lim):
                    out[csq_fields[i]] = parts[i]

                allele = out.get("Allele") or out.get("ALLELE")
                if allele:
                    out["ALT"] = allele
            else:
                out[csq_tag] = csq

            append(out)

    df = pd.DataFrame(rows)

    for col in ("POS", "QUAL"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "INFO_CALLERS" in df.columns:
        s = df["INFO_CALLERS"].fillna("")
        df["CALLER_M2"] = s.str.contains(r"(?:^|,)M2(?:,|$)", regex=True)
        df["CALLER_S2"] = s.str.contains(r"(?:^|,)S2(?:,|$)", regex=True)
        df["CALLER_DV"] = s.str.contains(r"(?:^|,)DV(?:,|$)", regex=True)
        df["N_CALLERS"] = (
            df["CALLER_M2"].astype("int8")
            + df["CALLER_S2"].astype("int8")
            + df["CALLER_DV"].astype("int8")
        )

    return df


###############################################################################
# hotspot
###############################################################################

def load_hotspots_df(hotspot_path: Path) -> pd.DataFrame:
    hs = pd.read_excel(
        hotspot_path,
        engine="xlrd" if hotspot_path.suffix.lower() == ".xls" else None,
        usecols=["Hugo_Symbol", "Amino_Acid_Position"],
    )

    hs = hs.rename(columns={
        "Hugo_Symbol": "SYMBOL",
        "Amino_Acid_Position": "Protein_position",
    })

    hs["SYMBOL"] = hs["SYMBOL"].astype("string")
    hs["Protein_position"] = pd.to_numeric(hs["Protein_position"], errors="coerce").astype("Int64")
    hs = hs.dropna(subset=["SYMBOL", "Protein_position"])
    hs = hs.drop_duplicates(subset=["SYMBOL", "Protein_position"]).reset_index(drop=True)
    return hs


def annotate_hotspots(df: pd.DataFrame, hs: pd.DataFrame, hotspot_col: str = DEFAULT_HOTSPOT_COL) -> pd.DataFrame:
    out = df

    out["SYMBOL"] = out.get("SYMBOL", pd.Series(pd.NA, index=out.index)).astype("string")
    out["Protein_position"] = pd.to_numeric(out.get("Protein_position", pd.Series(pd.NA, index=out.index)),
                                           errors="coerce").astype("Int64")

    hs_index = pd.MultiIndex.from_frame(hs[["SYMBOL", "Protein_position"]])
    df_index = pd.MultiIndex.from_frame(out[["SYMBOL", "Protein_position"]])

    out[hotspot_col] = df_index.isin(hs_index)
    return out


###############################################################################
# csv
###############################################################################

def write_csv(
    df: pd.DataFrame,
    out_path: Path,
    compress: bool = False,
    chunksize: Optional[int] = None,
) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if compress else None

    if chunksize is None:
        df.to_csv(
            out_path,
            index=False,
            sep=",",
            na_rep="",
            lineterminator="\n",
            compression=compression,
        )
        return out_path

    first = True
    for start in range(0, len(df), chunksize):
        chunk = df.iloc[start:start + chunksize]
        chunk.to_csv(
            out_path,
            index=False,
            sep=",",
            na_rep="",
            lineterminator="\n",
            compression=compression,
            mode="w" if first else "a",
            header=first,
        )
        first = False

    return out_path


###############################################################################
# filter noise
###############################################################################

def split_filter_cell(x: Any) -> set[str]:
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return set()
    s = str(x).strip()
    if not s or s in {".", "PASS"}:
        return set()
    s = s.replace(";", ",")
    return {t.strip() for t in s.split(",") if t.strip()}


def has_any_bad(filter_series: pd.Series, bad_set: set[str]) -> pd.Series:
    return filter_series.apply(lambda v: len(split_filter_cell(v) & bad_set) > 0)


def clinsig_keep(series: pd.Series, keep_terms: Tuple[str, ...]) -> pd.Series:
    s = series.fillna("").astype(str).str.lower()
    pat = "|".join(map(re.escape, keep_terms))
    return s.str.contains(pat, regex=True)


def safe_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def alt_count_from_ad_biallelic(ad_cell: Any) -> float:
    if ad_cell is None or (isinstance(ad_cell, float) and pd.isna(ad_cell)):
        return np.nan
    s = str(ad_cell).strip()
    if not s or s == ".":
        return np.nan
    s = s.split(";")[0]
    parts = [p for p in s.split(",") if p != ""]
    if len(parts) < 2:
        return np.nan
    try:
        return float(parts[1])
    except Exception:
        return np.nan


def add_aggregate_metrics(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()

    dp_cols = [c for c in ("FORMAT_M2_DP", "FORMAT_S2_DP", "FORMAT_DV_DP") if c in d.columns]
    for c in dp_cols:
        d[c] = safe_num(d[c])
    d["DP_AGG"] = d[dp_cols].max(axis=1) if dp_cols else np.nan

    if "FORMAT_S2_AD" in d.columns and "FORMAT_S2_DP" in d.columns:
        d["S2_ALT_i"] = d["FORMAT_S2_AD"].map(alt_count_from_ad_biallelic)
        d["S2_VAF_f"] = d["S2_ALT_i"] / safe_num(d["FORMAT_S2_DP"])

    vaf_parts = []
    if "FORMAT_M2_AF" in d.columns:
        vaf_parts.append(safe_num(d["FORMAT_M2_AF"]))
    if "FORMAT_DV_VAF" in d.columns:
        vaf_parts.append(safe_num(d["FORMAT_DV_VAF"]))
    if "S2_VAF_f" in d.columns:
        vaf_parts.append(safe_num(d["S2_VAF_f"]))

    d["VAF_AGG"] = pd.concat(vaf_parts, axis=1).max(axis=1) if vaf_parts else np.nan

    pop = safe_num(d["MAX_AF"]) if "MAX_AF" in d.columns else None
    if pop is None or pop.isna().all():
        pop = safe_num(d["AF"]) if "AF" in d.columns else None
    d["POP_AF"] = pop if pop is not None else np.nan

    return d


def filter_noise(
    df: pd.DataFrame,
    min_dp: int,
    min_vaf: float,
    max_pop_af: float,
    hotspot_col: str,
    clinsig_terms: Tuple[str, ...],
    m2_bad_filters: set[str],
    s2_bad_filters: set[str],
    dv_bad_filters: set[str],
) -> pd.DataFrame:
    d = add_aggregate_metrics(df)

    dp = safe_num(d["DP_AGG"]).fillna(0)
    vaf = safe_num(d["VAF_AGG"]).fillna(0)
    popaf = safe_num(d["POP_AF"])

    base_ok = (dp >= min_dp) & (vaf >= min_vaf) & (popaf.isna() | (popaf <= max_pop_af))

    hotspot_keep_mask = (
        d[hotspot_col].astype(bool)
        if hotspot_col in d.columns else pd.Series(False, index=d.index)
    )

    clinsig_keep_mask = (
        clinsig_keep(d["CLIN_SIG"], clinsig_terms)
        if "CLIN_SIG" in d.columns else pd.Series(False, index=d.index)
    )

    m2_bad = has_any_bad(d.get("INFO_M2_FILTER", pd.Series("", index=d.index)), m2_bad_filters)
    s2_bad = has_any_bad(d.get("INFO_S2_FILTER", pd.Series("", index=d.index)), s2_bad_filters)
    dv_bad = has_any_bad(d.get("INFO_DV_FILTER", pd.Series("", index=d.index)), dv_bad_filters)

    has_m2 = d.get("CALLER_M2", pd.Series(False, index=d.index)).astype(bool)
    has_s2 = d.get("CALLER_S2", pd.Series(False, index=d.index)).astype(bool)
    has_dv = d.get("CALLER_DV", pd.Series(False, index=d.index)).astype(bool)

    all_bad_supported = ((~has_m2) | m2_bad) & ((~has_s2) | s2_bad) & ((~has_dv) | dv_bad)

    keep = (base_ok & ~all_bad_supported) | hotspot_keep_mask | clinsig_keep_mask
    return d.loc[keep].copy()


###############################################################################
# xlsx table
###############################################################################

def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def first_nonnull(*series_list: pd.Series) -> pd.Series:
    out = series_list[0].copy()
    for s in series_list[1:]:
        out = out.where(~out.isna(), s)
    return out


def parse_ad_alt_sum(ad_str: Any) -> float:
    if ad_str is None or (isinstance(ad_str, float) and np.isnan(ad_str)):
        return np.nan
    s = str(ad_str).strip()
    if not s or s in {".", "None"}:
        return np.nan
    try:
        parts = [float(x) for x in s.replace(";", ",").split(",") if x != ""]
    except Exception:
        return np.nan
    if len(parts) < 2:
        return np.nan
    return float(np.nansum(parts[1:]))


def add_final_metrics(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()

    if "DP_AGG" in d.columns:
        d["DP_FINAL"] = to_num(d["DP_AGG"])
    else:
        dp_cols = [c for c in ("FORMAT_M2_DP", "FORMAT_S2_DP", "FORMAT_DV_DP") if c in d.columns]
        d["DP_FINAL"] = d[dp_cols].apply(to_num).max(axis=1) if dp_cols else np.nan

    if "VAF_AGG" in d.columns:
        d["VAF_FINAL"] = to_num(d["VAF_AGG"])
    else:
        vaf_cols = []
        if "FORMAT_M2_AF" in d.columns:
            vaf_cols.append(to_num(d["FORMAT_M2_AF"]))
        if "FORMAT_DV_VAF" in d.columns:
            vaf_cols.append(to_num(d["FORMAT_DV_VAF"]))
        d["VAF_FINAL"] = pd.concat(vaf_cols, axis=1).max(axis=1) if vaf_cols else np.nan

    m2_ad = d.get("FORMAT_M2_AD", pd.Series(np.nan, index=d.index, dtype="object"))
    s2_ad = d.get("FORMAT_S2_AD", pd.Series(np.nan, index=d.index, dtype="object"))
    dv_ad = d.get("FORMAT_DV_AD", pd.Series(np.nan, index=d.index, dtype="object"))
    d["AD_FINAL"] = first_nonnull(m2_ad, s2_ad, dv_ad)

    d["ALT_COUNT_FINAL"] = d["AD_FINAL"].map(parse_ad_alt_sum)

    max_af = to_num(d.get("MAX_AF", pd.Series(np.nan, index=d.index)))
    af = to_num(d.get("AF", pd.Series(np.nan, index=d.index)))
    d["POP_AF_FINAL"] = max_af.where(~max_af.isna(), af)

    return d


def split_interesting(
    df_clean: pd.DataFrame,
    keep_impact: set[str],
    drop_consequence: set[str],
    hotspot_col: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = df_clean.copy()

    hotspot_keep_mask = d.get(hotspot_col, pd.Series(False, index=d.index)).astype(bool)
    impact_ok = d.get("IMPACT", "").astype(str).isin(keep_impact)

    cons = d.get("Consequence", "").fillna("").astype(str)
    has_drop = cons.str.contains("|".join(drop_consequence), regex=True)
    consequence_ok = ~has_drop

    interesting = (impact_ok & consequence_ok) | hotspot_keep_mask
    return d.loc[interesting].copy(), d.loc[~interesting].copy()


def prepare_interpretation_tables(
    df_clean: pd.DataFrame,
    hotspot_col: str,
    keep_impact: set[str],
    drop_consequence: set[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = add_final_metrics(df_clean)

    rename_map = {
        "INFO_M2_FILTER": "M2_FILTER",
        "INFO_S2_FILTER": "S2_FILTER",
        "INFO_DV_FILTER": "DV_FILTER",
        "INFO_CALLERS": "CALLERS",
    }
    d = d.rename(columns={k: v for k, v in rename_map.items() if k in d.columns})

    final_rename = {
        "DP_FINAL": "DP",
        "VAF_FINAL": "VAF",
        "AD_FINAL": "AD",
        "POP_AF_FINAL": "POP_AF",
        "ALT_COUNT_FINAL": "ALT_COUNT",
    }
    d = d.rename(columns={k: v for k, v in final_rename.items() if k in d.columns})

    if "CALLERS" not in d.columns and "INFO_CALLERS" in d.columns:
        d["CALLERS"] = d["INFO_CALLERS"]

    df_main, df_other = split_interesting(d, keep_impact, drop_consequence, hotspot_col)

    desired = [
        "SYMBOL", "CHROM", "POS", "REF", "ALT",
        "HGVSc", "HGVSp", "Consequence", "IMPACT",
        "CANONICAL", "MANE_SELECT", "ENSP",
        "DP", "VAF", "AD",
        "MAX_AF", "AF",
        "CLIN_SIG", "SIFT", "PolyPhen", "PUBMED", "Existing_variation",
        hotspot_col, "PHENO", "GENE_PHENO", "CALLERS", "CALLER_M2",
        "CALLER_S2", "CALLER_DV", "N_CALLERS",
        "M2_FILTER", "S2_FILTER", "DV_FILTER",
    ]
    desired = [c for c in desired if c in d.columns]

    df_main = df_main.loc[:, desired].copy()
    df_other = df_other.loc[:, desired].copy()
    return df_main, df_other


###############################################################################
# METRICS sheet in xlsx
###############################################################################

def read_metrics_report_txt(path: Path) -> pd.DataFrame:

    text = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    kv: Dict[str, str] = {}

    for ln in text:
        ln = ln.strip()
        if not ln or ":" not in ln:
            continue
        k, v = ln.split(":", 1)
        kv[k.strip()] = v.strip()

    sample = kv.get("SAMPLE", "")
    if not sample and "SAMPLE" not in kv:
        sample = kv.get("SAMPLE", "")

    ordered = [
        ("SAMPLE", "SAMPLE"),
        ("Total reads", "Total reads"),
        ("Mean read length", "Mean read length"),
        ("Mean read quality", "Mean read quality"),
        ("% Q30 bases", "% Q30 bases"),
        ("% mapped reads", "% mapped reads"),
        ("% properly paired", "% properly paired"),
        ("Duplication rate (%)", "Duplication rate (%)"),
        ("Mean target coverage", "Mean target coverage"),
        ("Median target coverage", "Median target coverage"),
        ("% target bases 20x", "% target bases 20x"),
        ("% target bases 100x", "% target bases 100x"),
    ]

    row: Dict[str, Any] = {}
    for src, dst in ordered:
        if src == "SAMPLE":
            row[dst] = kv.get("SAMPLE", kv.get("SAMPLE", ""))
        else:
            row[dst] = kv.get(src)

    def _num(x: Any) -> Any:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        try:
            return float(s) if any(ch in s for ch in ".eE") else int(s)
        except Exception:
            try:
                return float(s)
            except Exception:
                return s

    for k in list(row.keys()):
        if k == "SAMPLE":
            continue
        row[k] = _num(row[k])

    if not row.get("SAMPLE"):
        m = re.search(r"^SAMPLE:\s*(\S+)\s*$", "\n".join(text), flags=re.MULTILINE)
        if m:
            row["SAMPLE"] = m.group(1)

    return pd.DataFrame([row])


###############################################################################
# xlsx writing
###############################################################################

def write_excel_report(
    df_main: pd.DataFrame,
    df_other: pd.DataFrame,
    info_xlsx_path: Path,
    out_xlsx_path: Path,
    metrics_df: Optional[pd.DataFrame] = None,
    max_rows_per_sheet: int = EXCEL_MAX_ROWS,
    engine: str = "xlsxwriter",
) -> None:
    out_xlsx_path.parent.mkdir(parents=True, exist_ok=True)
    info_book = pd.read_excel(info_xlsx_path, sheet_name=None)

    with pd.ExcelWriter(out_xlsx_path, engine=engine) as w:
        df_main.to_excel(w, sheet_name="MAIN", index=False)
        ws = w.sheets["MAIN"]
        ws.freeze_panes(1, 0)
        if len(df_main) > 0 and len(df_main.columns) > 0:
            ws.autofilter(0, 0, min(len(df_main), max_rows_per_sheet), len(df_main.columns) - 1)

        total = len(df_other)
        if total > 0:
            n_sheets = (total - 1) // max_rows_per_sheet + 1
            for i in range(n_sheets):
                start = i * max_rows_per_sheet
                end = min(start + max_rows_per_sheet, total)
                chunk = df_other.iloc[start:end]
                sheet_name = "OTHER" if i == 0 else f"OTHER{i+1}"
                chunk.to_excel(w, sheet_name=sheet_name, index=False)
                ws2 = w.sheets[sheet_name]
                ws2.freeze_panes(1, 0)
                if len(chunk) > 0 and len(chunk.columns) > 0:
                    ws2.autofilter(0, 0, len(chunk), len(chunk.columns) - 1)

        if metrics_df is not None:
            metrics_df.to_excel(w, sheet_name="METRICS", index=False)
            ws_m = w.sheets["METRICS"]
            ws_m.freeze_panes(1, 0)
            if len(metrics_df) > 0 and len(metrics_df.columns) > 0:
                ws_m.autofilter(0, 0, len(metrics_df), len(metrics_df.columns) - 1)

        if len(info_book) == 1:
            only_name = next(iter(info_book))
            info_book[only_name].to_excel(w, sheet_name="INFO", index=False)
            ws3 = w.sheets["INFO"]
            ws3.freeze_panes(1, 0)
            if len(info_book[only_name]) > 0 and len(info_book[only_name].columns) > 0:
                ws3.autofilter(0, 0, len(info_book[only_name]), len(info_book[only_name].columns) - 1)
        else:
            first = True
            for name, df_info in info_book.items():
                sheet = "INFO" if first else f"INFO_{name}"[:31]
                first = False
                df_info.to_excel(w, sheet_name=sheet, index=False)
                ws3 = w.sheets[sheet]
                ws3.freeze_panes(1, 0)
                if len(df_info) > 0 and len(df_info.columns) > 0:
                    ws3.autofilter(0, 0, len(df_info), len(df_info.columns) - 1)


###############################################################################
# CLI
###############################################################################

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Parse VEP-annotated union VCF (M2/S2/DV), annotate hotspots, filter noise, export CSV + Excel."
    )
    ap.add_argument("vcf", type=Path, help="VEP-annotated VCF (plain .vcf or .vcf.gz)")
    ap.add_argument("--sample", required=True, help="Sample name in VCF")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory")

    ap.add_argument("--hotspots", type=Path, required=True, help="Cancer hotspots Excel (.xls/.xlsx)")
    ap.add_argument("--info-xlsx", type=Path, required=True, help="Info list Excel for INFO sheet")
    ap.add_argument("--metrics-txt", type=Path, default=None, help="QC metrics text report to add METRICS sheet")

    ap.add_argument("--min-dp", type=int, default=DEFAULT_MIN_DP)
    ap.add_argument("--min-vaf", type=float, default=DEFAULT_MIN_VAF)
    ap.add_argument("--max-pop-af", type=float, default=DEFAULT_MAX_POP_AF)

    ap.add_argument("--csv-chunksize", type=int, default=500_000)
    ap.add_argument("--compress-csv", action="store_true")

    ap.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity (-v, -vv)")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    log = setup_logger(args.verbose)

    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Parsing VCF: %s", args.vcf)
    df_raw = parse_vep_union_vcf(
        vcf_path=args.vcf,
        sample_name=args.sample,
        explode_csq=True,
        keep_all_info=True,
        keep_all_format=True,
    )
    log.info("Parsed rows=%d cols=%d", len(df_raw), len(df_raw.columns))

    log.info("Loading hotspots: %s", args.hotspots)
    hs = load_hotspots_df(args.hotspots)

    log.info("Annotating hotspots")
    df_raw = annotate_hotspots(df_raw, hs, hotspot_col=DEFAULT_HOTSPOT_COL)

    raw_csv = outdir / f"{args.sample}_raw.csv"
    if args.compress_csv:
        raw_csv = raw_csv.with_suffix(".csv.gz")

    log.info("Writing raw CSV: %s", raw_csv)
    write_csv(df_raw, raw_csv, compress=args.compress_csv, chunksize=args.csv_chunksize)

    log.info("Filtering noise")
    df_clean = filter_noise(
        df=df_raw,
        min_dp=args.min_dp,
        min_vaf=args.min_vaf,
        max_pop_af=args.max_pop_af,
        hotspot_col=DEFAULT_HOTSPOT_COL,
        clinsig_terms=DEFAULT_KEEP_IF_CLINSIG_HAS,
        m2_bad_filters=DEFAULT_M2_BAD_FILTERS,
        s2_bad_filters=DEFAULT_S2_BAD_FILTERS,
        dv_bad_filters=DEFAULT_DV_BAD_FILTERS,
    )
    log.info("Clean rows=%d", len(df_clean))

    clean_csv = outdir / f"{args.sample}_filtered.csv"
    if args.compress_csv:
        clean_csv = outdir / f"{args.sample}_filtered.csv.gz"
    log.info("Writing clean CSV: %s", clean_csv)
    write_csv(df_clean, clean_csv, compress=args.compress_csv, chunksize=args.csv_chunksize)

    log.info("Preparing MAIN/OTHER for interpretation")
    df_main, df_other = prepare_interpretation_tables(
        df_clean=df_clean,
        hotspot_col=DEFAULT_HOTSPOT_COL,
        keep_impact=DEFAULT_KEEP_IMPACT,
        drop_consequence=DEFAULT_DROP_CONSEQUENCE,
    )

    metrics_df = None
    if args.metrics_txt is not None:
        log.info("Reading metrics report: %s", args.metrics_txt)
        metrics_df = read_metrics_report_txt(args.metrics_txt)

    out_xlsx = outdir / f"{args.sample}_vep_annotated.xlsx"
    log.info("Writing Excel: %s", out_xlsx)
    write_excel_report(
        df_main=df_main,
        df_other=df_other,
        info_xlsx_path=args.info_xlsx,
        out_xlsx_path=out_xlsx,
        metrics_df=metrics_df,
        max_rows_per_sheet=EXCEL_MAX_ROWS,
        engine="xlsxwriter",
    )

    log.warning("Done. Outputs:\n  %s\n  %s\n  %s", raw_csv, clean_csv, out_xlsx)


if __name__ == "__main__":
    main()