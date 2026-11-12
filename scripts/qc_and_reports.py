#!/usr/bin/env python3
"""
qc_and_reports.py -- MASPHL viral NGS, Phase 3: QC and result files.

Reads the loader's sample rows and the TheiaCoV results for one run and writes:

  qc_rows.tsv              assembly_status, qc_fail_reasons (+ masphl_id for MPXV)
                           for the data table
  reports/                 MPXV/WNV: <prefix>_results_<YYYY-MM-DD>.tsv
                           SC2:      <btb>_<YYYY-MM-DD>.run_results.csv
                                     <btb>_<YYYY-MM-DD>.lims_file.csv
  <prefix>_<YYYY-MM-DD>.run_summary.txt
  CONTAMINATION            "true"/"false"

QC (same rule for every sample, controls included):
  PASS when percent_reference_coverage >= min coverage (SC2 90, MPXV/WNV 80)
        and meanmapq_trim >= 20 and meanbaseq_trim >= 20
  (>= matches the validated SPHL LIMS prep code.)
Negative controls: contaminated when assembly_length_unambiguous >= 600.

SC2 files reproduce SPHL LIMS prep + gen v1.1.0 column for column, including
its defaults for missing values (numbers 0, text NA), INVALID lineage for FAIL,
and the run_results "nextclade_lineage" column holding nextclade_clade.
Differences from v1.1.0: organism is "SARS-CoV-2" (was "SARS-CoV 2"), no
trailing comma on run_results rows, and the date follows US Eastern time.

Standard library only. Embedded in tasks/task_qc_reports.wdl by
tools/embed_scripts.py. Do not use the WDL placeholder characters (tilde-brace)
or three closing angle brackets in this file.
"""

import argparse
import csv
import datetime
import os
import sys

SCRIPT_VERSION = "0.3.0"

MPXV_COLUMNS = """batchid assembly_status assembly_length_unambiguous assembly_mean_coverage
assembly_method consensus_n_variant_min_depth fastq_scan_num_reads_clean_pairs
fastq_scan_num_reads_raw_pairs kraken_human kraken_human_dehosted kraken_sc2 kraken_sc2_dehosted
kraken_target_organism kraken_target_organism_dehosted kraken_target_organism_name masphl_id
meanbaseq_trim meanmapq_trim nextclade_aa_dels nextclade_aa_subs nextclade_clade nextclade_docker
nextclade_ds_tag nextclade_lineage nextclade_qc nextclade_version number_Degenerate number_N
number_Total percent_reference_coverage primer_bed_name primer_trimmed_read_percent
read_screen_clean read_screen_raw sample_sheet_id theiacov_illumina_pe_analysis_date
theiacov_illumina_pe_version""".split()

WNV_COLUMNS = """batchid assembly_status assembly_length_unambiguous assembly_mean_coverage
assembly_method consensus_n_variant_min_depth fastq_scan_num_reads_clean_pairs
fastq_scan_num_reads_raw_pairs kraken_human kraken_human_dehosted kraken_target_organism
kraken_target_organism_dehosted kraken_target_organism_name meanbaseq_trim meanmapq_trim
number_Degenerate number_N number_Total percent_reference_coverage primer_bed_name
primer_trimmed_read_percent read_screen_clean read_screen_raw theiacov_illumina_pe_analysis_date
theiacov_illumina_pe_version""".split()

# added after the SOP columns in the MPXV/WNV reports
EXTRA_REPORT_COLUMNS = {"MPXV": ["analysis_route", "nextclade_clade_all_clades",
                                 "negative_control", "qc_fail_reasons"],
                        "WNV": ["negative_control", "qc_fail_reasons"]}

RUN_RESULTS_HEADER = ("sample_id,batch_id,seq_date,assembly_status,percent_reference_coverage,"
                      "mean_depth,meanbaseq_trim,meanmapq_trim,nextclade_lineage,pangolin_lineage,"
                      "qc_reads_raw,qc_reads_clean,%_human_reads,%_SARS-COV-2_reads,dehosted_%human,"
                      "dehosted_%SC2,num_N,num_degenerate,num_ACTG,num_total,pangolin_version,"
                      "AA_substitutions,AA_deletions")
LIMS_HEADER = "sample_id,assembly_status,tool_lineage,lineage_to_maven,pango_version,organism,test"


def read_tsv(path):
    """Rows keyed by sample ID (the entity:<table>_id column)."""
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(rd, None)
        if not header:
            return {}, []
        out = {}
        for r in rd:
            if r:
                out[r[0]] = dict(zip(header[1:], r[1:]))
        return out, header


def num(v):
    if v is None:
        return None
    v = v.strip()
    if v in ("", "NA", "N/A", "nan", "NaN"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def fmt(x):
    return f"{x:g}" if x is not None else "missing"


# --------------------------------------------------------------------------
# QC
# --------------------------------------------------------------------------

def qc_sample(row, opts):
    """Return (status, reasons, contaminated or None)."""
    reasons = []
    if row.get("load_status") == "no_reads":
        reasons.append("no reads in BaseSpace")
    elif not row.get("_analyzed"):
        reasons.append("no TheiaCoV results")
    cov, mapq, baseq = (num(row.get("percent_reference_coverage")), num(row.get("meanmapq_trim")),
                        num(row.get("meanbaseq_trim")))
    if not reasons:
        if cov is None or cov < opts["min_cov"]:
            reasons.append(f"percent_reference_coverage {fmt(cov)} < {fmt(opts['min_cov'])}")
        if mapq is None or mapq < opts["min_mapq"]:
            reasons.append(f"meanmapq_trim {fmt(mapq)} < {fmt(opts['min_mapq'])}")
        if baseq is None or baseq < opts["min_baseq"]:
            reasons.append(f"meanbaseq_trim {fmt(baseq)} < {fmt(opts['min_baseq'])}")
    contaminated = None
    if row.get("negative_control") == "yes":
        acgt = num(row.get("assembly_length_unambiguous")) or 0
        contaminated = acgt >= opts["control_max"]
    return ("FAIL" if reasons else "PASS"), reasons, contaminated


def masphl_id(sample_id, row, prefix):
    year = (row.get("basespace_run_date") or "")[:4]
    if not year.isdigit():
        year = str(datetime.date.today().year)
    return f"{prefix}{year}_{sample_id}"


# --------------------------------------------------------------------------
# SC2 LIMS files (reproduce SPHL LIMS prep + gen v1.1.0)
# --------------------------------------------------------------------------

def _text(v, default="NA"):
    return v if v not in (None, "") else default


def _float_str(v):
    return v if v not in (None, "") else "0.0"


def _int_str(v):
    return v if v not in (None, "") else "0"


def sc2_rows(samples, order):
    run_results, lims = [], []
    for sid in order:
        r = samples[sid]
        status = r["assembly_status"]
        pango = _text(r.get("pango_lineage"))
        lineage = pango if status == "PASS" else "INVALID"
        pango_version = _text(r.get("pangolin_assignment_version"))
        run_results.append([
            sid, r.get("batchid", ""), r.get("theiacov_illumina_pe_analysis_date", ""), status,
            _float_str(r.get("percent_reference_coverage")),
            _float_str(r.get("assembly_mean_coverage")),
            _float_str(r.get("meanbaseq_trim")), _float_str(r.get("meanmapq_trim")),
            _text(r.get("nextclade_clade")).replace(",", ""),       # header says nextclade_lineage
            pango,
            _int_str(r.get("fastq_scan_num_reads_raw_pairs")),
            _int_str(r.get("fastq_scan_num_reads_clean_pairs")),
            _float_str(r.get("kraken_human")), _float_str(r.get("kraken_sc2")),
            _float_str(r.get("kraken_human_dehosted")), _float_str(r.get("kraken_sc2_dehosted")),
            _int_str(r.get("number_N")), _int_str(r.get("number_Degenerate")),
            _int_str(r.get("assembly_length_unambiguous")), _int_str(r.get("number_Total")),
            pango_version,
            _text(r.get("nextclade_aa_subs")).replace(",", "|"),
            _text(r.get("nextclade_aa_dels")).replace(",", "|"),
        ])
        lims.append([sid, status, lineage, lineage, pango_version, "SARS-CoV-2", "SARS-CoV-2 Sequencing"])
    return run_results, lims


def write_plain_csv(path, header, rows):
    # values never contain commas after the replacements above; write them as-is
    # (like LIMS gen did) so BtB sees exactly the same text
    with open(path, "w", newline="") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write(",".join(str(v) for v in r) + "\n")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def build(opts):
    loader, _h = read_tsv(opts["loader_rows"])
    results, _h2 = read_tsv(opts["theiacov_results"]) if opts["theiacov_results"] else ({}, [])
    order = list(loader.keys())
    samples = {}
    for sid in order:
        r = dict(loader[sid])
        if sid in results:
            r.update(results[sid])
            r["_analyzed"] = True
        samples[sid] = r

    pathogen = opts["pathogen"]
    qc_out, alerts, controls, fails = [], [], [], []
    for sid in order:
        r = samples[sid]
        status, reasons, contaminated = qc_sample(r, opts)
        r["assembly_status"] = status
        r["qc_fail_reasons"] = "; ".join(reasons)
        if pathogen == "MPXV":
            r["masphl_id"] = masphl_id(sid, r, opts["masphl_id_prefix"])
        qc_out.append(r)
        if contaminated is not None:
            controls.append((sid, r.get("assembly_length_unambiguous") or "0",
                             r.get("fastq_scan_num_reads_raw_pairs") or r.get("basespace_read_pairs_fetched") or "",
                             contaminated))
            if contaminated:
                alerts.append(sid)
        elif status == "FAIL":
            fails.append((sid, r["qc_fail_reasons"]))
    return samples, order, qc_out, alerts, controls, fails


def write_all(opts, samples, order, alerts, controls, fails, check_summary_text, replaced):
    out = opts["out_dir"]
    os.makedirs(os.path.join(out, "reports"), exist_ok=True)
    pathogen, prefix, date = opts["pathogen"], opts["file_prefix"], opts["date"]
    table = opts["table_name"]

    qc_cols = ["assembly_status", "qc_fail_reasons"] + (["masphl_id"] if pathogen == "MPXV" else []) \
        + ["masphl_qc_version"]
    with open(os.path.join(out, "qc_rows.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([f"entity:{table}_id"] + qc_cols)
        for sid in order:
            r = samples[sid]
            w.writerow([sid] + [r.get(c, "") for c in qc_cols[:-1]] + [opts["version"]])

    report_files = []
    if pathogen == "SC2":
        rr, lims = sc2_rows(samples, order)
        p1 = os.path.join(out, "reports", f"{prefix}_{date}.run_results.csv")
        p2 = os.path.join(out, "reports", f"{prefix}_{date}.lims_file.csv")
        write_plain_csv(p1, RUN_RESULTS_HEADER, rr)
        write_plain_csv(p2, LIMS_HEADER, lims)
        report_files += [p1, p2]
    else:
        cols = (MPXV_COLUMNS if pathogen == "MPXV" else WNV_COLUMNS) + EXTRA_REPORT_COLUMNS[pathogen]
        p = os.path.join(out, "reports", f"{prefix}_results_{date}.tsv")
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sample_id"] + cols)
            for sid in order:
                r = samples[sid]
                w.writerow([sid] + [r.get(c, "") for c in cols])
        report_files.append(p)

    n_pass = sum(1 for s in order if samples[s]["assembly_status"] == "PASS")
    lines = [
        f"MASPHL viral NGS run summary -- {pathogen} -- {prefix} -- {date}",
        f"QC version {SCRIPT_VERSION}; workflow {opts['version']}",
        "",
        f"QC rule: percent_reference_coverage >= {fmt(opts['min_cov'])}, meanmapq_trim >= {fmt(opts['min_mapq'])}, "
        f"meanbaseq_trim >= {fmt(opts['min_baseq'])}",
        f"Samples: {len(order)}   PASS: {n_pass}   FAIL: {len(order) - n_pass}",
        "",
    ]
    if alerts:
        lines += ["!" * 70,
                  "CONTAMINATION ALERT -- notify the Laboratory Supervisor",
                  f"Negative control(s) with assembly_length_unambiguous >= {opts['control_max']}: "
                  + ", ".join(alerts),
                  "!" * 70, ""]
    lines.append("NEGATIVE CONTROLS (contaminated when unambiguous bases >= "
                 f"{opts['control_max']})")
    if controls:
        for sid, acgt, raw, bad in controls:
            lines.append(f"  {sid:20} unambiguous bases {acgt:>8}   raw read pairs {raw:>10}   "
                         f"{'CONTAMINATED' if bad else 'clean'}")
    else:
        lines.append("  none identified in this run")
    lines.append("")
    lines.append("FAILED SAMPLES (controls excluded)")
    if fails:
        if pathogen == "SC2":
            lines.append("  Check Ct values in BtB before reporting INVALID (Ct < 30 -> resequence).")
        for sid, why in fails:
            lines.append(f"  {sid:20} {why}")
    else:
        lines.append("  none")
    lines.append("")
    if pathogen == "MPXV":
        rerun = [s for s in order if samples[s].get("analysis_route") == "mpxv_clade_i_rerun"]
        lines.append("MPXV CLADE I RERUNS: " + (", ".join(
            f"{s} ({samples[s].get('nextclade_clade_all_clades', '')})" for s in rerun) or "none"))
        lines.append("")
    if replaced:
        lines.append("ROWS REPLACED FROM AN EARLIER BASESPACE RUN")
        lines += [f"  {sid:20} previous run {old}, previous batchid {b}" for sid, old, b in replaced]
        lines.append("")
    lines.append("FILES")
    lines += ["  " + os.path.basename(p) for p in report_files]
    lines += ["", "-" * 70, "LOADER CHECK", check_summary_text.rstrip()]
    summary_path = os.path.join(out, f"{prefix}_{date}.run_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(out, "CONTAMINATION"), "w") as fh:
        fh.write("true" if alerts else "false")
    return summary_path, report_files


def read_replaced(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return [(r["sample_id"], r["previous_basespace_run_id"], r.get("previous_batchid", ""))
                for r in rd]


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--loader-rows", required=True)
    p.add_argument("--theiacov-results", default="")
    p.add_argument("--check-summary", default="")
    p.add_argument("--replaced-rows", default="")
    p.add_argument("--pathogen", required=True)
    p.add_argument("--file-prefix", required=True)
    p.add_argument("--min-percent-reference-coverage", default="")
    p.add_argument("--min-meanmapq", default="20")
    p.add_argument("--min-meanbaseq", default="20")
    p.add_argument("--control-max-unambiguous-bases", default="600")
    p.add_argument("--masphl-id-prefix", default="MPXV_USA_")
    p.add_argument("--date", default="", help="YYYY-MM-DD for file names (default: today)")
    p.add_argument("--table-name", default="sample")
    p.add_argument("--version", default="dev")
    p.add_argument("--out-dir", required=True)
    a = p.parse_args(argv)

    pathogen = a.pathogen.strip().upper()
    if pathogen not in ("SC2", "MPXV", "WNV"):
        sys.stderr.write(f"ERROR: unknown pathogen {a.pathogen}\n")
        return 1
    min_cov = float(a.min_percent_reference_coverage) if a.min_percent_reference_coverage.strip() \
        else (90.0 if pathogen == "SC2" else 80.0)
    opts = {
        "loader_rows": a.loader_rows, "theiacov_results": a.theiacov_results,
        "pathogen": pathogen, "file_prefix": a.file_prefix, "min_cov": min_cov,
        "min_mapq": float(a.min_meanmapq), "min_baseq": float(a.min_meanbaseq),
        "control_max": int(float(a.control_max_unambiguous_bases)),
        "masphl_id_prefix": a.masphl_id_prefix, "date": a.date or datetime.date.today().isoformat(),
        "table_name": a.table_name, "version": a.version, "out_dir": a.out_dir,
    }
    samples, order, _qc, alerts, controls, fails = build(opts)
    if not order:
        sys.stderr.write("ERROR: no samples in the loader rows\n")
        return 1
    check_text = open(a.check_summary).read() if a.check_summary and os.path.exists(a.check_summary) else ""
    summary, _files = write_all(opts, samples, order, alerts, controls, fails, check_text,
                                read_replaced(a.replaced_rows))
    with open(summary) as fh:
        sys.stdout.write(fh.read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
