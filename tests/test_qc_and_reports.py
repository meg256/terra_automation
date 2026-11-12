"""
Phase 3: QC PASS/FAIL, MPXV/WNV reports, SC2 run_results + lims_file.

MPXV: real run MA26-32-Mpox-20260911 (export from sample_test).
SC2:  the anonymised run_results / lims_file examples produced by SPHL LIMS
      prep + gen v1.1.0 are rebuilt from TheiaCoV-style inputs and must match
      line for line (apart from the documented changes).
"""
import csv
import os
import sys

HERE = os.path.dirname(__file__)
FIX = os.path.join(HERE, "fixtures")
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import qc_and_reports as q  # noqa: E402


def run(tmp_path, pathogen, loader, results, prefix, **kw):
    out = tmp_path / "out"
    args = ["--loader-rows", str(loader), "--theiacov-results", str(results), "--pathogen", pathogen,
            "--file-prefix", prefix, "--date", "2026-11-12", "--table-name", "sample_test",
            "--version", "test", "--out-dir", str(out)]
    for k, v in kw.items():
        args += ["--" + k.replace("_", "-"), v]
    assert q.main(args) == 0
    return out


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------ MPXV ---

def test_mpxv_real_run(tmp_path):
    out = run(tmp_path, "MPXV", os.path.join(FIX, "mpxv_MA26-32_loader_rows.tsv"),
              os.path.join(FIX, "mpxv_MA26-32_theiacov_results.tsv"), "MA26-32-Mpox-20260911")
    qc = {r["entity:sample_test_id"]: r for r in tsv(out / "qc_rows.tsv")}
    assert {k: v["assembly_status"] for k, v in qc.items()} == {
        "MA0271": "PASS", "MA0272": "PASS", "MA0273": "FAIL", "NTCMA2632": "FAIL"}
    assert qc["MA0273"]["qc_fail_reasons"] == "percent_reference_coverage 77.38 < 80"
    assert qc["MA0271"]["masphl_id"] == "MPXV_USA_2026_MA0271"        # year of the BaseSpace run

    report = tsv(out / "reports" / "MA26-32-Mpox-20260911_results_2026-11-12.tsv")
    assert list(report[0].keys())[:3] == ["sample_id", "batchid", "assembly_status"]
    assert [c for c in q.MPXV_COLUMNS if c not in report[0]] == []
    assert report[1]["analysis_route"] == "mpxv_clade_i_rerun"

    summary = (out / "MA26-32-Mpox-20260911_2026-11-12.run_summary.txt").read_text()
    assert "CONTAMINATION ALERT" not in summary
    assert "NTCMA2632" in summary and "clean" in summary and "3423222" in summary
    assert "MPXV CLADE I RERUNS: MA0272 (Ib), MA0273 (Ib)" in summary
    assert (out / "CONTAMINATION").read_text() == "false"


def test_mpxv_contaminated_control(tmp_path):
    rows = list(csv.reader(open(os.path.join(FIX, "mpxv_MA26-32_theiacov_results.tsv"), newline=""),
                           delimiter="\t"))
    col = rows[0].index("assembly_length_unambiguous")
    for r in rows[1:]:
        if r[0] == "NTCMA2632":
            r[col] = "600"
    res = tmp_path / "res.tsv"
    with open(res, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    out = run(tmp_path, "MPXV", os.path.join(FIX, "mpxv_MA26-32_loader_rows.tsv"), res, "X")
    s = (out / "X_2026-11-12.run_summary.txt").read_text()
    assert "CONTAMINATION ALERT" in s and "NTCMA2632" in s          # exactly 600 = contaminated
    assert (out / "CONTAMINATION").read_text() == "true"


def test_threshold_is_greater_or_equal(tmp_path):
    loader = tmp_path / "l.tsv"
    loader.write_text("entity:sample_test_id\tbatchid\tnegative_control\tload_status\n"
                      "A\tB1\tno\tloaded\nB\tB1\tno\tloaded\n")
    res = tmp_path / "r.tsv"
    res.write_text("entity:sample_test_id\tpercent_reference_coverage\tmeanmapq_trim\tmeanbaseq_trim\n"
                   "A\t80\t20\t20\nB\t79.99\t20\t19.99\n")
    out = run(tmp_path, "WNV", loader, res, "B1")
    qc = {r["entity:sample_test_id"]: r for r in tsv(out / "qc_rows.tsv")}
    assert qc["A"]["assembly_status"] == "PASS"
    assert qc["B"]["qc_fail_reasons"] == ("percent_reference_coverage 79.99 < 80; "
                                          "meanbaseq_trim 19.99 < 20")
    assert "masphl_id" not in qc["A"]                               # MPXV only


def test_no_reads_sample_fails_and_is_reported(tmp_path):
    loader = tmp_path / "l.tsv"
    loader.write_text("entity:sample_test_id\tbatchid\tnegative_control\tload_status\n"
                      "NTC1\tB1\tyes\tno_reads\nA\tB1\tno\tloaded\n")
    res = tmp_path / "r.tsv"
    res.write_text("entity:sample_test_id\tpercent_reference_coverage\tmeanmapq_trim\tmeanbaseq_trim\t"
                   "assembly_length_unambiguous\nA\t95\t50\t35\t10000\n")
    out = run(tmp_path, "WNV", loader, res, "B1")
    qc = {r["entity:sample_test_id"]: r for r in tsv(out / "qc_rows.tsv")}
    assert qc["NTC1"]["assembly_status"] == "FAIL"
    assert qc["NTC1"]["qc_fail_reasons"] == "no reads in BaseSpace"
    s = (out / "B1_2026-11-12.run_summary.txt").read_text()
    assert "NTC1" in s and "clean" in s


def test_wnv_report_columns(tmp_path):
    out = run(tmp_path, "WNV", os.path.join(FIX, "mpxv_MA26-32_loader_rows.tsv"),
              os.path.join(FIX, "mpxv_MA26-32_theiacov_results.tsv"), "WNV4")
    rep = tsv(out / "reports" / "WNV4_results_2026-11-12.tsv")
    assert list(rep[0].keys()) == ["sample_id"] + q.WNV_COLUMNS + ["negative_control", "qc_fail_reasons"]


# ------------------------------------------------------------------- SC2 ---

def sc2_inputs_from_example(tmp_path):
    """Turn the example run_results.csv back into loader rows + TheiaCoV results."""
    ex = list(csv.DictReader(open(os.path.join(FIX, "sc2_example_run_results.csv"), newline="")))
    loader = tmp_path / "loader.tsv"
    res = tmp_path / "results.tsv"
    with open(loader, "w") as fh:
        fh.write("entity:sample_test_id\tbatchid\tnegative_control\tload_status\n")
        for r in ex:
            ctrl = "yes" if r["sample_id"].startswith(("NEG", "NTC")) else "no"
            fh.write(f"{r['sample_id']}\t{r['batch_id']}\t{ctrl}\tloaded\n")
    cols = {"theiacov_illumina_pe_analysis_date": "seq_date",
            "percent_reference_coverage": "percent_reference_coverage",
            "assembly_mean_coverage": "mean_depth", "meanbaseq_trim": "meanbaseq_trim",
            "meanmapq_trim": "meanmapq_trim", "nextclade_clade": "nextclade_lineage",
            "pango_lineage": "pangolin_lineage", "fastq_scan_num_reads_raw_pairs": "qc_reads_raw",
            "fastq_scan_num_reads_clean_pairs": "qc_reads_clean", "kraken_human": "%_human_reads",
            "kraken_sc2": "%_SARS-COV-2_reads", "kraken_human_dehosted": "dehosted_%human",
            "kraken_sc2_dehosted": "dehosted_%SC2", "number_N": "num_N",
            "number_Degenerate": "num_degenerate", "assembly_length_unambiguous": "num_ACTG",
            "number_Total": "num_total", "pangolin_assignment_version": "pangolin_version",
            "nextclade_aa_subs": "AA_substitutions", "nextclade_aa_dels": "AA_deletions"}
    with open(res, "w") as fh:
        fh.write("entity:sample_test_id\t" + "\t".join(cols) + "\n")
        for r in ex:
            vals = []
            for tc, rr in cols.items():
                v = r[rr]
                if tc in ("nextclade_aa_subs", "nextclade_aa_dels"):
                    v = "" if v == "NA" else v.replace("|", ",")   # TheiaCoV uses commas
                if tc == "nextclade_clade" and v == "NA":
                    v = ""
                vals.append(v)
            fh.write(r["sample_id"] + "\t" + "\t".join(vals) + "\n")
    return loader, res


def test_sc2_files_match_lims_gen_v110(tmp_path):
    loader, res = sc2_inputs_from_example(tmp_path)
    out = run(tmp_path, "SC2", loader, res, "28591")
    got_rr = (out / "reports" / "28591_2026-11-12.run_results.csv").read_text().splitlines()
    want_rr = open(os.path.join(FIX, "sc2_example_run_results.csv")).read().splitlines()
    assert got_rr == want_rr          # example has no trailing comma; v1.1.0 added one (dropped on purpose)

    got_lims = (out / "reports" / "28591_2026-11-12.lims_file.csv").read_text().splitlines()
    want_lims = [l.replace("SARS-CoV 2,", "SARS-CoV-2,")       # organism string corrected on purpose
                 for l in open(os.path.join(FIX, "sc2_example_lims_file.csv")).read().splitlines()]
    assert got_lims == want_lims


def test_sc2_summary_lists_ct_review(tmp_path):
    loader, res = sc2_inputs_from_example(tmp_path)
    out = run(tmp_path, "SC2", loader, res, "28591")
    s = (out / "28591_2026-11-12.run_summary.txt").read_text()
    assert "Check Ct values in BtB" in s
    fails_block = s.split("FAILED SAMPLES")[1].split("\n\n")[0]
    assert "25MA1003" in fails_block and "NTC" not in fails_block   # controls listed separately


def test_sc2_missing_values_use_lims_prep_defaults(tmp_path):
    loader = tmp_path / "l.tsv"
    loader.write_text("entity:sample_test_id\tbatchid\tnegative_control\tload_status\n"
                      "NTC9\t28591\tyes\tno_reads\n")
    res = tmp_path / "r.tsv"
    res.write_text("entity:sample_test_id\tpercent_reference_coverage\n")
    out = run(tmp_path, "SC2", loader, res, "28591")
    rr = (out / "reports" / "28591_2026-11-12.run_results.csv").read_text().splitlines()[1]
    assert rr == "NTC9,28591,,FAIL,0.0,0.0,0.0,0.0,NA,NA,0,0,0.0,0.0,0.0,0.0,0,0,0,0,NA,NA,NA"
    lims = (out / "reports" / "28591_2026-11-12.lims_file.csv").read_text().splitlines()[1]
    assert lims == "NTC9,FAIL,INVALID,INVALID,NA,SARS-CoV-2,SARS-CoV-2 Sequencing"
