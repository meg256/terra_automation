"""
Simulated BaseSpace runs for check_run_and_samples.py.

Each test builds a fake snapshot (what get_basespace_run_info would download)
and runs the real script on it. These cover the cases we can't test on real
MASPHL data yet: NextSeq lane splits, multi-project runs, duplicate datasets,
and BaseSpace's '_' -> '-' renaming.

Run:  python -m pytest tests/ -q
"""
import csv
import json
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import check_run_and_samples as crs  # noqa: E402

RUN_ID = "332775445"
FASTQ = "illumina.fastq.v1.8"


# ---------------------------------------------------------------- helpers ---

def ds(name, project, ds_id, reads=50000, dtype=FASTQ, date="2026-11-10T15:00:00Z"):
    return {"Id": ds_id, "Name": name, "Project.Name": project, "Project.Id": "p_" + project,
            "DataSetType.Id": dtype, "DateCreated": date, "_reads": reads}


def v1_sheet(rows):
    """rows: (Sample_ID, Sample_Name, Sample_Project)"""
    head = ("[Header],,,\nLocal Run Manager Analysis Id,1,,\nExperiment Name,RUN,,\n,,,\n"
            "[Reads],,,\n151,,,\n151,,,\n,,,\n[Data],,,\n"
            "Sample_ID,Sample_Name,Description,Sample_Project\n")
    return head + "".join(f"{a},{b},,{c}\n" for a, b, c in rows)


def v2_sheet(rows):
    """MiSeq i100 / NextSeq 2000 style sample sheet (v2)."""
    out = ("[Header],\nFileFormatVersion,2\nRunName,RUN\n\n[BCLConvert_Settings],\n"
           "SoftwareVersion,4.2\n\n[BCLConvert_Data],,,\nSample_ID,Index,Index2\n")
    out += "".join(f"{a},AAAA,CCCC\n" for a, _b, _c in rows)
    out += "\n[Cloud_Data],,\nSample_ID,ProjectName,LibraryName\n"
    out += "".join(f"{a},{c},{b}\n" for a, b, c in rows)
    return out


def make_snapshot(tmp_path, datasets, sheet=None, status="Complete", attrs=True):
    snap = tmp_path / "snap"
    (snap / "attributes").mkdir(parents=True)
    (snap / "run.json").write_text(json.dumps({
        "Id": int(RUN_ID), "Name": "MA26-12-WNV-20261110", "ExperimentName": "MA26-12-WNV-20261110",
        "Status": status, "DateCreated": "2026-11-10T14:00:00.0000000Z"}))
    with open(snap / "datasets.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Id", "Name", "Project.Name", "Project.Id", "DataSetType.Id", "DateCreated"])
        for d in datasets:
            w.writerow([d["Id"], d["Name"], d["Project.Name"], d["Project.Id"],
                        d["DataSetType.Id"], d["DateCreated"]])
            if attrs and d["_reads"] is not None:
                (snap / "attributes" / f"{d['Id']}.csv").write_text(
                    f"Name,Value\nTotalReadsPF,{d['_reads']}\nIsPairedEnd,true\n")
    if sheet is not None:
        (snap / "run_files" / "Files").mkdir(parents=True)
        (snap / "run_files" / "Files" / "SampleSheet.csv").write_text(sheet)
    return snap


def run(tmp_path, snap, **kw):
    out = tmp_path / "out"
    args = ["--snapshot-dir", str(snap), "--out-dir", str(out), "--run-id", kw.pop("run_id", RUN_ID),
            "--pathogen", kw.pop("pathogen", "WNV"), "--loader-version", "test"]
    for k, v in kw.items():
        args += ["--" + k.replace("_", "-"), v]
    code = crs.main(args)
    summary = (out / "check_summary.txt").read_text()
    return code, out, summary


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fetch(out):
    with open(out / "to_fetch.tsv") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


STD = [ds("NTC-01", "WNV-PROJ", "ds.1", reads=120),
       ds("MA0101", "WNV-PROJ", "ds.2"),
       ds("MA0102", "WNV-PROJ", "ds.3")]
STD_SHEET = v1_sheet([("NTC-01", "NTC-01", "MA26-12-WNV-20261110"),
                      ("MA0101", "MA0101", "MA26-12-WNV-20261110"),
                      ("MA0102", "MA0102", "MA26-12-WNV-20261110")])


# ------------------------------------------------------------------ tests ---

def test_normal_single_project_run(tmp_path):
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, STD, STD_SHEET), project="WNV-PROJ")
    assert code == 0
    r = rows(out / "sample_rows.tsv")
    assert [x["entity:sample_id"] for x in r] == ["NTC-01", "MA0101", "MA0102"]
    assert {x["batchid"] for x in r} == {"MA26-12-WNV-20261110"}      # from Sample_Project
    assert [x["negative_control"] for x in r] == ["yes", "no", "no"]
    assert (out / "sample_set_name.txt").read_text() == "MA26-12-WNV-20261110"
    assert fetch(out) == [["NTC-01", "ds.1"], ["MA0101", "ds.2"], ["MA0102", "ds.3"]]
    assert "WARNINGS" not in summary


def test_sc2_requires_five_digit_btb_id(tmp_path):
    snap = make_snapshot(tmp_path, STD, STD_SHEET)
    code, _o, summary = run(tmp_path, snap, pathogen="SC2", project="WNV-PROJ")
    assert code == 1 and "btb_batch_id is required" in summary
    code, _o, summary = run(tmp_path, snap, pathogen="SC2", project="WNV-PROJ", btb_batch_id="2859")
    assert code == 1 and "exactly 5 digits" in summary
    code, out, _s = run(tmp_path, snap, pathogen="SC2", project="WNV-PROJ", btb_batch_id="28591")
    assert code == 0
    assert {x["batchid"] for x in rows(out / "sample_rows.tsv")} == {"28591"}
    assert (out / "sample_set_name.txt").read_text() == "28591"


def test_btb_id_ignored_for_wnv(tmp_path):
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, STD, STD_SHEET),
                        project="WNV-PROJ", btb_batch_id="12345")
    assert code == 0
    assert {x["batchid"] for x in rows(out / "sample_rows.tsv")} == {"MA26-12-WNV-20261110"}


def test_bad_run_id_and_missing_project(tmp_path):
    snap = make_snapshot(tmp_path, STD, STD_SHEET)
    code, _o, summary = run(tmp_path, snap, run_id="MA26-12-WNV", project="WNV-PROJ")
    assert code == 1 and "number from the run's BaseSpace URL" in summary
    code, _o, summary = run(tmp_path, snap)
    assert code == 1 and "basespace_project is empty" in summary


def test_project_not_in_run_lists_projects(tmp_path):
    code, _o, summary = run(tmp_path, make_snapshot(tmp_path, STD, STD_SHEET), project="Mpox-PROJ")
    assert code == 1
    assert "Project 'Mpox-PROJ' has no FASTQ datasets" in summary and "WNV-PROJ" in summary


def test_multi_project_run_filters_to_target(tmp_path):
    data = STD + [ds("WW-001", "Wastewater", "ds.9"), ds("NTC-WW", "Wastewater", "ds.10", reads=40)]
    sheet = v1_sheet([("NTC-01", "NTC-01", "WNV-PROJ"), ("MA0101", "MA0101", "WNV-PROJ"),
                      ("MA0102", "MA0102", "WNV-PROJ"), ("WW-001", "WW-001", "Wastewater"),
                      ("NTC-WW", "NTC-WW", "Wastewater")])
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, data, sheet), project="WNV-PROJ")
    assert code == 0
    assert [x["entity:sample_id"] for x in rows(out / "sample_rows.tsv")] == ["NTC-01", "MA0101", "MA0102"]
    report = {x["basespace_name"]: x["status"] for x in rows(out / "check_report.tsv")}
    assert report["WW-001"] == "other_project" and report["NTC-WW"] == "other_project"
    assert "WARNINGS" not in summary       # wastewater rows on the sheet are not "missing"


def test_nextseq_lanes_are_joined_in_lane_order(tmp_path):
    data = [ds("MA0201_L002", "P", "ds.b", reads=10), ds("MA0201_L001", "P", "ds.a", reads=12),
            ds("MA0201_L004", "P", "ds.d", reads=9), ds("MA0201_L003", "P", "ds.c", reads=11),
            ds("MA02010_L001", "P", "ds.z", reads=5)]          # must NOT join MA0201
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, data), project="P")
    assert code == 0
    assert fetch(out) == [["MA0201", "ds.a,ds.b,ds.c,ds.d"], ["MA02010", "ds.z"]]
    r = {x["entity:sample_id"]: x for x in rows(out / "sample_rows.tsv")}
    assert r["MA0201"]["basespace_reads_pf"] == "42"


def test_underscore_renamed_to_hyphen_keeps_original_id(tmp_path):
    data = [ds("MA-0301", "P", "ds.1"), ds("NTC-1", "P", "ds.2", reads=3)]
    sheet = v1_sheet([("MA_0301", "MA_0301", "P"), ("NTC_1", "NTC_1", "P")])
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, data, sheet), project="P")
    assert code == 0
    r = rows(out / "sample_rows.tsv")
    assert [x["entity:sample_id"] for x in r] == ["MA_0301", "NTC_1"]
    assert r[1]["negative_control"] == "yes"
    assert "WARNINGS" not in summary


def test_sample_id_differs_from_sample_name(tmp_path):
    # the SAMIDI2 / SAMID2 case from the example sheet: dataset is named after Sample_ID
    data = [ds("SAMIDI2", "P", "ds.1")]
    sheet = v1_sheet([("SAMIDI2", "SAMID2", "P")])
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, data, sheet), project="P")
    assert code == 0
    r = rows(out / "sample_rows.tsv")[0]
    assert r["entity:sample_id"] == "SAMIDI2" and r["sample_sheet_id"] == "SAMID2"


def test_duplicate_dataset_is_reported_not_loaded(tmp_path):
    data = STD + [ds("MA0101", "WNV-PROJ", "ds.2b", date="2026-11-11T09:00:00Z")]
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, data, STD_SHEET), project="WNV-PROJ")
    assert code == 0
    assert "MA0101" not in [x["entity:sample_id"] for x in rows(out / "sample_rows.tsv")]
    assert "requeued sample" in summary
    report = {x["dataset_ids"]: x["status"] for x in rows(out / "check_report.tsv")}
    assert report["ds.2,ds.2b"] == "duplicate_datasets"


def test_name_collision_after_normalising(tmp_path):
    data = [ds("MA_01", "P", "ds.1"), ds("MA-01", "P", "ds.2")]
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, data), project="P")
    # neither is loaded -> nothing left -> fails
    assert code == 1 and "No samples to load" in summary


def test_undetermined_and_non_fastq_are_ignored(tmp_path):
    data = STD + [ds("Undetermined", "WNV-PROJ", "ds.u"),
                  ds("some-analysis", "WNV-PROJ", "ds.x", dtype="common.files")]
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, data, STD_SHEET), project="WNV-PROJ")
    assert code == 0
    assert len(rows(out / "sample_rows.tsv")) == 3


def test_zero_read_control_added_but_not_fetched(tmp_path):
    data = [ds("NTC-01", "P", "ds.1", reads=0), ds("MA0101", "P", "ds.2")]
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, data), project="P")
    assert code == 0
    r = {x["entity:sample_id"]: x for x in rows(out / "sample_rows.tsv")}
    assert r["NTC-01"]["load_status"] == "no_reads"
    assert fetch(out) == [["MA0101", "ds.2"]]
    assert "WARNINGS" not in summary         # empty control is expected -> note, not warning


def test_zero_read_patient_sample_is_a_warning(tmp_path):
    data = [ds("NTC-01", "P", "ds.1"), ds("MA0101", "P", "ds.2", reads=0)]
    code, _o, summary = run(tmp_path, make_snapshot(tmp_path, data), project="P")
    assert code == 0 and "WARNINGS" in summary and "MA0101: 0 reads" in summary


def test_unknown_read_count_still_fetched(tmp_path):
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, STD, STD_SHEET, attrs=False), project="WNV-PROJ")
    assert code == 0 and len(fetch(out)) == 3


def test_samplesheet_mismatches_are_warnings(tmp_path):
    sheet = v1_sheet([("MA0101", "MA0101", "WNV-PROJ"), ("MA0102", "MA0102", "WNV-PROJ"),
                      ("MA0199", "MA0199", "WNV-PROJ")])
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, STD, sheet), project="WNV-PROJ")
    assert code == 0
    assert "MA0199: on the sample sheet but has no FASTQ dataset" in summary
    assert "NTC-01: has FASTQs in BaseSpace but is not on the sample sheet" in summary


def test_no_samplesheet_uses_project_as_batchid(tmp_path):
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, STD), project="WNV-PROJ")
    assert code == 0
    assert {x["batchid"] for x in rows(out / "sample_rows.tsv")} == {"WNV-PROJ"}
    assert "No sample sheet found" in summary


def test_v2_samplesheet_miseq_i100(tmp_path):
    sheet = v2_sheet([("NTC-01", "NTC-01", "MA26-15-Mpox-20261111"),
                      ("MA0101", "MA0101", "MA26-15-Mpox-20261111"),
                      ("MA0102", "MA0102", "MA26-15-Mpox-20261111")])
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, STD, sheet), project="WNV-PROJ", pathogen="MPXV")
    assert code == 0
    assert {x["batchid"] for x in rows(out / "sample_rows.tsv")} == {"MA26-15-Mpox-20261111"}


def test_only_these_samples_any_project_no_auto_controls(tmp_path):
    data = STD + [ds("WW-001", "Wastewater", "ds.9")]
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, data, STD_SHEET),
                             project="WNV-PROJ", only_these_samples=" MA0102,WW-001 ")
    assert code == 0
    assert [x["entity:sample_id"] for x in rows(out / "sample_rows.tsv")] == ["MA0102", "WW-001"]
    assert "No negative controls" in summary


def test_only_these_samples_missing_id_fails(tmp_path):
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, STD, STD_SHEET),
                             only_these_samples="MA0101, MA9999")
    assert code == 1 and "MA9999" in summary and "not in run" in summary
    assert not (out / "sample_rows.tsv").exists()


def test_only_these_samples_space_separated_fails(tmp_path):
    code, _o, summary = run(tmp_path, make_snapshot(tmp_path, STD, STD_SHEET),
                            only_these_samples="MA0101 MA0102")
    assert code == 1 and "Separate sample IDs with commas" in summary


def test_only_these_samples_in_two_projects_fails(tmp_path):
    data = STD + [ds("MA0101", "Other", "ds.7")]
    code, _o, summary = run(tmp_path, make_snapshot(tmp_path, data), only_these_samples="MA0101")
    assert code == 1 and "more than one project" in summary


def test_only_these_samples_with_underscore_typed(tmp_path):
    data = [ds("MA-0301", "P", "ds.1")]
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, data), only_these_samples="MA_0301")
    assert code == 0
    assert rows(out / "sample_rows.tsv")[0]["entity:sample_id"] == "MA_0301"


def test_extra_negative_control(tmp_path):
    data = STD + [ds("BLANK7", "WNV-PROJ", "ds.8", reads=10)]
    code, out, summary = run(tmp_path, make_snapshot(tmp_path, data), project="WNV-PROJ",
                             extra_negative_control_ids="BLANK7, NOTHERE")
    assert code == 0
    r = {x["entity:sample_id"]: x["negative_control"] for x in rows(out / "sample_rows.tsv")}
    assert r["BLANK7"] == "yes" and r["MA0101"] == "no"
    assert "lists 'NOTHERE', but no such sample was loaded" in summary


def test_control_prefix_is_case_insensitive(tmp_path):
    data = [ds("ntc-a", "P", "ds.1", reads=5), ds("NegMA01", "P", "ds.2", reads=5), ds("MA1", "P", "ds.3")]
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, data), project="P")
    r = {x["entity:sample_id"]: x["negative_control"] for x in rows(out / "sample_rows.tsv")}
    assert r == {"ntc-a": "yes", "NegMA01": "yes", "MA1": "no"}


def test_incomplete_run(tmp_path):
    snap = make_snapshot(tmp_path, STD, STD_SHEET, status="Running")
    code, _o, summary = run(tmp_path, snap, project="WNV-PROJ")
    assert code == 1 and "not 'Complete'" in summary
    code, _o, summary = run(tmp_path, snap, project="WNV-PROJ", allow_incomplete_run="true")
    assert code == 0 and "allow_incomplete_run = true" in summary


def test_invalid_terra_name_fails(tmp_path):
    data = [ds("MA 01", "P", "ds.1")]
    code, _o, summary = run(tmp_path, make_snapshot(tmp_path, data), project="P")
    assert code == 1 and "characters Terra does not allow" in summary


def test_custom_table_name_headers(tmp_path):
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, STD, STD_SHEET),
                        project="WNV-PROJ", table_name="sample_test")
    assert code == 0
    assert (out / "sample_rows.tsv").read_text().startswith("entity:sample_test_id\t")
    assert (out / "membership.tsv").read_text().startswith("membership:sample_test_set_id\tsample_test\n")


def test_real_example_samplesheet_parses():
    """The anonymised SampleSheet.csv from the conversation (v1, Local Run Manager)."""
    path = os.path.join(HERE, "fixtures", "SampleSheet_example_v1.csv")
    got = crs.parse_samplesheet(path)
    assert got == [
        {"sample_id": "NEGMAXXXXCoV", "sample_name": "NEGMAXXXXCoV", "project": "RUN_NAME"},
        {"sample_id": "SAMID1", "sample_name": "SAMID1", "project": "RUN_NAME"},
        {"sample_id": "SAMIDI2", "sample_name": "SAMID2", "project": "RUN_NAME"},
    ]


def test_embedded_copy_in_wdl_matches_script():
    """tasks/task_basespace.wdl must carry the current version of the script."""
    sys.path.insert(0, os.path.join(HERE, "..", "tools"))
    import embed_scripts
    assert embed_scripts.check(), "Embedded script is stale: run python tools/embed_scripts.py"


def test_bcl_convert_common_fastq_type(tmp_path):
    """Run 335865530: BCL Convert writes FASTQs as common.fastq, plus common.files folders."""
    data = [ds("MA0278", "MA26-38-Mpox", "ds.1", dtype="common.fastq"),
            ds("MA2638MpoxNTC", "MA26-38-Mpox", "ds.2", reads=50, dtype="common.fastq"),
            ds("CM26-1123", "MA26-38-WNV", "ds.3", dtype="common.fastq"),
            ds("MA0278 (BCL Convert)", "MA26-38-WNV-Mpox-20261113_335865530", "ds.4", dtype="common.files"),
            ds("BCL Convert Logs", "MA26-38-WNV-Mpox-20261113_335865530", "ds.5", dtype="common.files")]
    code, out, _s = run(tmp_path, make_snapshot(tmp_path, data), project="MA26-38-Mpox", pathogen="MPXV",
                        extra_negative_control_ids="MA2638MpoxNTC")
    assert code == 0
    r = {x["entity:sample_id"]: x for x in rows(out / "sample_rows.tsv")}
    assert set(r) == {"MA0278", "MA2638MpoxNTC"}
    assert r["MA2638MpoxNTC"]["negative_control"] == "yes"
    t2 = tmp_path / "second"
    t2.mkdir()
    code, out, _s = run(t2, make_snapshot(t2, data), pathogen="MPXV",
                        only_these_samples="MA0278, MA2638MpoxNTC", extra_negative_control_ids="MA2638MpoxNTC")
    assert code == 0 and fetch(out) == [["MA0278", "ds.1"], ["MA2638MpoxNTC", "ds.2"]]


def test_no_fastq_message_lists_types(tmp_path):
    data = [ds("X (BCL Convert)", "P", "ds.1", dtype="common.files")]
    code, _o, summary = run(tmp_path, make_snapshot(tmp_path, data), project="P")
    assert code == 1 and "common.files x1" in summary
