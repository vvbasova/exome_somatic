#!/usr/bin/env python3
import argparse
import logging
import os
import sys
import tempfile
import gzip
import pysam


#header normalization

def norm_number(n):
    """Normalize VCF header Number field to a valid VCF token."""
    if n is None:
        return "."
    if isinstance(n, int):
        return "." if n < 0 else str(n)
    return str(n)


def norm_type(t):
    """Normalize VCF header Type field to allowed VCF types."""
    if t in ("Integer", "Float", "String", "Flag"):
        return t
    return "String"


def safe_desc(desc, prefix):
    """Prefix descriptions to keep provenance of copied tags."""
    return f"[{prefix}] {(desc or '').strip()}".strip()


def record_key(rec):
    """Key for biallelic records only."""
    if not rec.alts or len(rec.alts) != 1:
        return None
    return (rec.contig, rec.pos, rec.ref, rec.alts[0])


def to_python(v):
    """Convert pysam/numpy-like values to plain Python scalars/tuples."""
    try:
        import numpy as np
        if isinstance(v, np.generic):
            v = v.item()
    except Exception:
        pass

    if isinstance(v, list):
        return tuple(to_python(x) for x in v)
    if isinstance(v, tuple):
        return tuple(to_python(x) for x in v)
    return v


def gt_tuple_to_str(gt, phased=False):
    """Convert pysam GT tuple into VCF-like string (0/1)."""
    if gt is None:
        return "./."
    if not isinstance(gt, (tuple, list)):
        return str(gt)

    sep = "|" if phased else "/"
    parts = ["." if a is None else str(a) for a in gt]

    if len(parts) == 0:
        return "./."
    if len(parts) == 1:
        return parts[0]
    return sep.join(parts)


def pl_to_str(pl):
    """Store PL as string."""
    if pl is None:
        return "."
    pl = to_python(pl)
    if isinstance(pl, tuple):
        return ",".join(str(x) for x in pl)
    return str(pl)


def load_records(path, log):
    """Load biallelic records into a dict keyed by (contig,pos,ref,alt)."""
    vf = pysam.VariantFile(path)
    m = {}
    skipped_multi = 0

    for rec in vf.fetch():
        k = record_key(rec)
        if k is None:
            skipped_multi += 1
            continue
        m[k] = rec

    if skipped_multi:
        log.warning(
            "File %s: skipped %d non-biallelic records (expected biallelic inputs).",
            path, skipped_multi
        )
    return vf, m


#trimming vectors

def expected_len(number_code, n_alleles):
    """
    Determine expected vector length for INFO/FORMAT given VCF Number code.
    """
    if number_code is None:
        return None

    if isinstance(number_code, int):
        if number_code < 0:
            return None
        return number_code

    s = str(number_code)
    if s == ".":
        return None
    if s == "A":
        return max(n_alleles - 1, 0)   # biallelic -> 1
    if s == "R":
        return n_alleles              # biallelic -> 2
    if s == "G":
        return None                   # do not enforce G here
    try:
        return int(s)
    except Exception:
        return None


def trim_to_expected(v, exp_len):
    """Trim tuple vectors to expected length"""
    if exp_len is None:
        return v
    if isinstance(v, tuple):
        return v[:exp_len] if len(v) > exp_len else v
    return v


def rewrite_filter_column_to_dot(in_vcf_gz, out_vcf_gz, log):
    """ FILTER='.' in out file"""
    log.info("Rewriting FILTER column to '.' in %s -> %s", in_vcf_gz, out_vcf_gz)

    with gzip.open(in_vcf_gz, "rt") as fin, pysam.BGZFile(out_vcf_gz, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line.encode("utf-8"))
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 7:
                cols[6] = "."  # FILTER column
            fout.write(("\t".join(cols) + "\n").encode("utf-8"))


#logging


def setup_logger(log_path, level):
    log = logging.getLogger("vcfmerge")
    log.setLevel(level)
    log.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")
    handler = logging.StreamHandler(sys.stderr) if not log_path else logging.FileHandler(log_path)
    handler.setFormatter(fmt)
    log.addHandler(handler)
    return log


#main

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Union-merge 3 normalized biallelic VCFs into one VCF with prefixed INFO/FORMAT. "
            "Caller GT/PL stored as strings. Allele-specific vectors are trimmed to match output header. "
            "Output FILTER is set to '.' for all records."
        )
    )
    ap.add_argument("m2")
    ap.add_argument("s2")
    ap.add_argument("dv")
    ap.add_argument("out")          # final .vcf.gz
    ap.add_argument("sample_name")

    # Strelka AD number fix
    ap.add_argument("--no-fix-strelka-ad-number", action="store_true",
                    help="Disable: force S2_AD/S2_ADF/S2_ADR Number=R in output header (default: enabled).")

    # Logging
    ap.add_argument("--log", default=None,
                    help="Log file path. Default: stderr.")
    ap.add_argument("--log-level", default="INFO",
                    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                    help="Logging verbosity (default: INFO).")

    args = ap.parse_args()
    log = setup_logger(args.log, getattr(logging, args.log_level))

    fix_strelka_ad = not args.no_fix_strelka_ad_number
    log.info("fix_strelka_ad_number=%s", fix_strelka_ad)

    # Load
    m2_vf, m2_map = load_records(args.m2, log)
    s2_vf, s2_map = load_records(args.s2, log)
    dv_vf, dv_map = load_records(args.dv, log)

    keys = set(m2_map) | set(s2_map) | set(dv_map)
    log.info("Union variants: %d", len(keys))

    # Build header
    out_header = pysam.VariantHeader()

    # contigs
    contigs = set()
    for vf in (m2_vf, s2_vf, dv_vf):
        contigs.update(list(vf.header.contigs))
    for ctg in sorted(contigs):
        length = None
        for vf in (m2_vf, s2_vf, dv_vf):
            if ctg in vf.header.contigs and vf.header.contigs[ctg].length:
                length = vf.header.contigs[ctg].length
                break
        if length:
            out_header.contigs.add(ctg, length=length)
        else:
            out_header.contigs.add(ctg)

    out_header.add_sample(args.sample_name)

    out_header.add_line(
        '##INFO=<ID=CALLERS,Number=1,Type=String,Description="Callers supporting this variant (comma-separated from {M2,S2,DV})">'
    )

    def add_prefixed_defs(vf, prefix):
        # INFO
        for info_id, info in vf.header.info.items():
            pid = f"{prefix}_{info_id}"
            if pid in out_header.info:
                continue
            out_header.add_line(
                f'##INFO=<ID={pid},Number={norm_number(info.number)},Type={norm_type(info.type)},Description="{safe_desc(info.description, prefix)}">'
            )

        # FORMAT
        for fmt_id, fmt in vf.header.formats.items():
            pid = f"{prefix}_{fmt_id}"
            if pid in out_header.formats:
                continue

            if fmt_id in ("GT", "PL"):
                out_header.add_line(
                    f'##FORMAT=<ID={pid},Number=1,Type=String,Description="[{prefix}] {fmt_id} stored as string">'
                )
                continue

            num = norm_number(fmt.number)
            typ = norm_type(fmt.type)
            desc = safe_desc(fmt.description, prefix)

            if fix_strelka_ad and prefix == "S2" and fmt_id in ("AD", "ADF", "ADR") and num == ".":
                num = "R"

            out_header.add_line(
                f'##FORMAT=<ID={pid},Number={num},Type={typ},Description="{desc}">'
            )

        # QUAL/FILTER into INFO
        for extra in ("QUAL", "FILTER"):
            pid = f"{prefix}_{extra}"
            if pid in out_header.info:
                continue
            out_header.add_line(
                f'##INFO=<ID={pid},Number=1,Type=String,Description="[{prefix}] original VCF {extra} column">'
            )

    add_prefixed_defs(m2_vf, "M2")
    add_prefixed_defs(s2_vf, "S2")
    add_prefixed_defs(dv_vf, "DV")

    out_dir = os.path.dirname(os.path.abspath(args.out)) or "."
    os.makedirs(out_dir, exist_ok=True)

    with tempfile.NamedTemporaryFile(prefix="merge_tmp_", suffix=".vcf.gz", dir=out_dir, delete=False) as tmp:
        tmp_path = tmp.name

    log.info("Writing intermediate merged VCF: %s", tmp_path)
    out_vf = pysam.VariantFile(tmp_path, "wz", header=out_header)

    def set_info_safe(rec_out, tag, value):
        value = to_python(value)
        meta = rec_out.header.info.get(tag)
        exp_len = expected_len(getattr(meta, "number", None), n_alleles=len(rec_out.alleles))
        value2 = trim_to_expected(value, exp_len)

        if isinstance(value2, bool):
            if value2:
                rec_out.info[tag] = True
            return

        rec_out.info[tag] = value2

    def copy_from(rec_src, prefix, rec_out):
        # QUAL/FILTER into INFO
        rec_out.info[f"{prefix}_QUAL"] = "." if rec_src.qual is None else str(rec_src.qual)

        filt = list(rec_src.filter.keys())
        rec_out.info[f"{prefix}_FILTER"] = "PASS" if len(filt) == 0 else ",".join(filt)

        # INFO
        for k in rec_src.info.keys():
            pk = f"{prefix}_{k}"
            set_info_safe(rec_out, pk, rec_src.info.get(k))

        # FORMAT
        if len(rec_src.samples) > 0:
            src_sname = list(rec_src.samples.keys())[0]
            sdata = rec_src.samples[src_sname]
            out_s = rec_out.samples[args.sample_name]

            phased = getattr(sdata, "phased", False)
            out_s[f"{prefix}_GT"] = gt_tuple_to_str(sdata.get("GT"), phased=phased)
            out_s[f"{prefix}_PL"] = pl_to_str(sdata.get("PL"))

            for fk in sdata.keys():
                if fk in ("GT", "PL"):
                    continue
                out_tag = f"{prefix}_{fk}"
                val = to_python(sdata.get(fk))

                meta = rec_out.header.formats.get(out_tag)
                exp_len = expected_len(getattr(meta, "number", None), n_alleles=len(rec_out.alleles))
                val2 = trim_to_expected(val, exp_len)

                out_s[out_tag] = val2

    for k in sorted(keys, key=lambda x: (x[0], x[1], x[2], x[3])):
        ctg, pos, ref, alt = k
        rec_out = out_vf.new_record(contig=ctg, start=pos - 1, alleles=(ref, alt))
        rec_out.id = "."
        rec_out.qual = None

        callers = []
        m2_rec = m2_map.get(k)
        s2_rec = s2_map.get(k)
        dv_rec = dv_map.get(k)

        if m2_rec is not None:
            callers.append("M2")
            copy_from(m2_rec, "M2", rec_out)
        if s2_rec is not None:
            callers.append("S2")
            copy_from(s2_rec, "S2", rec_out)
        if dv_rec is not None:
            callers.append("DV")
            copy_from(dv_rec, "DV", rec_out)

        rec_out.info["CALLERS"] = ",".join(callers)
        out_vf.write(rec_out)

    out_vf.close()

    rewrite_filter_column_to_dot(tmp_path, args.out, log)

    try:
        pysam.tabix_index(args.out, preset="vcf", force=True)
        log.info("Tabix index written: %s.tbi", args.out)
    except Exception as e:
        log.warning("Could not tabix-index output (%s). You can run: tabix -p vcf %s", e, args.out)

    try:
        os.unlink(tmp_path)
        tbi_tmp = tmp_path + ".tbi"
        if os.path.exists(tbi_tmp):
            os.unlink(tbi_tmp)
    except Exception:
        pass

    log.info("DONE: %s", args.out)


if __name__ == "__main__":
    main()