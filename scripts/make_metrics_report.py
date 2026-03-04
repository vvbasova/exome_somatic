#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple
import pandas as pd


def _to_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        if isinstance(x, str) and x.strip() == "":
            return None
        return float(x)
    except Exception:
        return None


def _to_int(x: Any) -> Optional[int]:
    try:
        if x is None:
            return None
        if isinstance(x, str) and x.strip() == "":
            return None
        return int(float(x))
    except Exception:
        return None


def _pct(x: Optional[float]) -> Optional[float]:
    if x is None:
        return None
    return x * 100.0


def _fmt_num(x: Optional[float], ndigits: int = 2) -> str:
    if x is None:
        return "NA"
    return f"{x:.{ndigits}f}"


def _fmt_int(x: Optional[int]) -> str:
    if x is None:
        return "NA"
    return str(x)


def _mean_q_from_quality_curves(report: dict, key: str) -> Optional[float]:
    blk = report.get(key, {}) or {}
    qc = blk.get("quality_curves", {}) or {}
    arr = qc.get("mean")
    if not isinstance(arr, list) or len(arr) == 0:
        return None

    vals = []
    for x in arr:
        fx = _to_float(x)
        if fx is not None:
            vals.append(fx)

    if not vals:
        return None
    return sum(vals) / len(vals)


def _read_picard_table(path: Path, metrics_class_contains: str) -> pd.DataFrame:
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()

    start_idx = None
    for i, ln in enumerate(lines):
        if ln.startswith("## METRICS CLASS") and metrics_class_contains in ln:
            start_idx = i
            break
    if start_idx is None:
        raise ValueError(f"Could not find METRICS CLASS containing '{metrics_class_contains}' in {path}")

    j = start_idx + 1
    while j < len(lines) and (not lines[j].strip() or lines[j].startswith("##")):
        j += 1
    if j >= len(lines):
        raise ValueError(f"Malformed picard file {path}: missing header row after metrics class section")

    header = lines[j].rstrip("\n").split("\t")
    ncol = len(header)
    j += 1

    data_rows = []
    while j < len(lines):
        ln = lines[j]
        if not ln.strip():
            break
        if ln.startswith("##"):
            break

        parts = ln.rstrip("\n").split("\t")
        if len(parts) < ncol:
            parts = parts + [""] * (ncol - len(parts))
        elif len(parts) > ncol:
            parts = parts[: ncol - 1] + ["\t".join(parts[ncol - 1 :])]

        data_rows.append(parts)
        j += 1

    if not data_rows:
        raise ValueError(f"No data rows found for metrics class '{metrics_class_contains}' in {path}")

    return pd.DataFrame(data_rows, columns=header)


def parse_fastp_report(fastp_json: Path, prefer_after_filtering: bool = True) -> Dict[str, Any]:
    report = json.loads(Path(fastp_json).read_text(encoding="utf-8", errors="replace"))

    summary = report.get("summary", {})
    before = summary.get("before_filtering", {}) or {}
    after = summary.get("after_filtering", {}) or {}
    chosen = after if prefer_after_filtering else before

    total_reads = _to_int(chosen.get("total_reads"))
    q30_pct = _pct(_to_float(chosen.get("q30_rate")))

    r1 = _mean_q_from_quality_curves(report, "read1_after_filtering")
    r2 = _mean_q_from_quality_curves(report, "read2_after_filtering")
    if r1 is None and r2 is None:
        r1 = _mean_q_from_quality_curves(report, "read1_before_filtering")
        r2 = _mean_q_from_quality_curves(report, "read2_before_filtering")

    if r1 is None and r2 is None:
        mean_read_quality = None
    elif r1 is None:
        mean_read_quality = r2
    elif r2 is None:
        mean_read_quality = r1
    else:
        mean_read_quality = (r1 + r2) / 2.0

    dup = report.get("duplication", {}) or {}
    dup_pct = _pct(_to_float(dup.get("rate")))

    return {
        "FASTP_TOTAL_READS": total_reads,
        "FASTP_Q30_PCT": q30_pct,
        "FASTP_MEAN_READ_QUALITY": mean_read_quality,
        "FASTP_DUPLICATION_PCT": dup_pct,
        "FASTP_USED_SECTION": "after_filtering" if prefer_after_filtering else "before_filtering",
    }


def parse_alignment_summary_metrics(alignment_metrics_txt: Path) -> Dict[str, Any]:
    df = _read_picard_table(alignment_metrics_txt, "picard.analysis.AlignmentSummaryMetrics")
    if "CATEGORY" not in df.columns:
        raise ValueError(f"AlignmentSummaryMetrics table missing CATEGORY column in {alignment_metrics_txt}")

    pair = df[df["CATEGORY"] == "PAIR"]
    row = pair.iloc[0] if not pair.empty else df.iloc[-1]

    return {
        "ALIGN_TOTAL_READS": _to_int(row.get("TOTAL_READS")),
        "ALIGN_MAPPED_PCT": _pct(_to_float(row.get("PCT_PF_READS_ALIGNED"))),
        "ALIGN_PROPERLY_PAIRED_PCT": _pct(_to_float(row.get("PCT_READS_ALIGNED_IN_PAIRS"))),
        "ALIGN_MEAN_READ_LENGTH": _to_float(row.get("MEAN_READ_LENGTH")),
    }


def parse_hs_metrics(hs_metrics_txt: Path) -> Dict[str, Any]:
    df = _read_picard_table(hs_metrics_txt, "picard.analysis.directed.HsMetrics")
    row = df.iloc[0]

    return {
        "HS_MEAN_TARGET_COVERAGE": _to_float(row.get("MEAN_TARGET_COVERAGE")),
        "HS_MEDIAN_TARGET_COVERAGE": _to_float(row.get("MEDIAN_TARGET_COVERAGE")),
        "HS_PCT_TARGET_BASES_20X": _pct(_to_float(row.get("PCT_TARGET_BASES_20X"))),
        "HS_PCT_TARGET_BASES_100X": _pct(_to_float(row.get("PCT_TARGET_BASES_100X"))),
    }


def build_sample_qc(
    sample: str,
    fastp_json: Path,
    alignment_metrics_txt: Path,
    hs_metrics_txt: Path,
) -> Tuple[pd.DataFrame, str]:
    fp = parse_fastp_report(fastp_json, prefer_after_filtering=True)
    am = parse_alignment_summary_metrics(alignment_metrics_txt)
    hs = parse_hs_metrics(hs_metrics_txt)

    total_reads = am.get("ALIGN_TOTAL_READS") or fp.get("FASTP_TOTAL_READS")

    metrics_row = {
        "SAMPLE": sample,
        "Total reads": total_reads,
        "Mean read length": am.get("ALIGN_MEAN_READ_LENGTH"),
        "Mean read quality": fp.get("FASTP_MEAN_READ_QUALITY"),
        "% Q30 bases": fp.get("FASTP_Q30_PCT"),
        "% mapped reads": am.get("ALIGN_MAPPED_PCT"),
        "% properly paired": am.get("ALIGN_PROPERLY_PAIRED_PCT"),
        "Duplication rate (%)": fp.get("FASTP_DUPLICATION_PCT"),
        "Mean target coverage": hs.get("HS_MEAN_TARGET_COVERAGE"),
        "Median target coverage": hs.get("HS_MEDIAN_TARGET_COVERAGE"),
        "% target bases 20x": hs.get("HS_PCT_TARGET_BASES_20X"),
        "% target bases 100x": hs.get("HS_PCT_TARGET_BASES_100X"),
    }
    metrics_df = pd.DataFrame([metrics_row])

    txt_lines = [
        f"SAMPLE: {sample}",
        "",
        f"Total reads: {_fmt_int(total_reads)}",
        f"Mean read length: {_fmt_num(am.get('ALIGN_MEAN_READ_LENGTH'))}",
        f"Mean read quality: {_fmt_num(fp.get('FASTP_MEAN_READ_QUALITY'))}",
        f"% Q30 bases: {_fmt_num(fp.get('FASTP_Q30_PCT'))}",
        f"% mapped reads: {_fmt_num(am.get('ALIGN_MAPPED_PCT'))}",
        f"% properly paired: {_fmt_num(am.get('ALIGN_PROPERLY_PAIRED_PCT'))}",
        f"Duplication rate (%): {_fmt_num(fp.get('FASTP_DUPLICATION_PCT'))}",
        f"Mean target coverage: {_fmt_num(hs.get('HS_MEAN_TARGET_COVERAGE'))}",
        f"Median target coverage: {_fmt_num(hs.get('HS_MEDIAN_TARGET_COVERAGE'))}",
        f"% target bases 20x: {_fmt_num(hs.get('HS_PCT_TARGET_BASES_20X'))}",
        f"% target bases 100x: {_fmt_num(hs.get('HS_PCT_TARGET_BASES_100X'))}",
    ]
    qc_txt = "\n".join(txt_lines)
    return metrics_df, qc_txt


def _cli() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Parse fastp + Picard metrics and produce QC metrics sheet + txt report.")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastp-json", type=Path, required=True)
    ap.add_argument("--alignment-metrics", type=Path, required=True)
    ap.add_argument("--hs-metrics", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--prefix", default=None, help="Output prefix (default: sample)")
    return ap.parse_args()


def main() -> None:
    args = _cli()
    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = args.prefix or args.sample

    metrics_df, qc_txt = build_sample_qc(
        sample=args.sample,
        fastp_json=args.fastp_json,
        alignment_metrics_txt=args.alignment_metrics,
        hs_metrics_txt=args.hs_metrics,
    )

    out_txt = outdir / f"{prefix}.qc_report.txt"

    out_txt.write_text(qc_txt, encoding="utf-8")

    print(f"[OK] Wrote: {out_txt}")

if __name__ == "__main__":
    main()