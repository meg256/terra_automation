"""Phase 2: quick input check, MPXV clade I trigger, generated TheiaCoV columns."""
import os
import re
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import check_run_and_samples as crs  # noqa: E402
import gen_theiacov_outputs as gen  # noqa: E402


def validate(tmp_path, problems=None, **kw):
    out = tmp_path / "out"
    args = ["--validate-only", "--out-dir", str(out), "--run-id", kw.pop("run_id", "332775445"),
            "--pathogen", kw.pop("pathogen", "WNV")]
    if problems is not None:
        pf = tmp_path / "problems.txt"
        pf.write_text("\n".join(problems) + "\n")
        args += ["--problems-file", str(pf)]
    for k, v in kw.items():
        args += ["--" + k.replace("_", "-"), v]
    code = crs.main(args)
    return code, (out / "check_summary.txt").read_text(), out


def test_validate_only_ok_needs_no_snapshot(tmp_path):
    code, summary, out = validate(tmp_path, project="P")
    assert code == 0 and summary.startswith("Inputs OK")
    # the WDL task reads these files even in validate-only mode
    for f in ("sample_rows.tsv", "membership.tsv", "to_fetch.tsv", "check_report.tsv", "sample_set_name.txt"):
        assert (out / f).exists()


def test_validate_only_catches_missing_btb_fast(tmp_path):
    code, summary, _ = validate(tmp_path, pathogen="SC2", project="P")
    assert code == 1 and "btb_batch_id is required" in summary


def test_workflow_problems_are_reported_with_input_problems(tmp_path):
    code, summary, _ = validate(tmp_path, pathogen="SC2", project="P",
                                problems=["mpxv_clade_i_reference_genome must be set"])
    assert code == 1
    assert "mpxv_clade_i_reference_genome must be set" in summary
    assert "btb_batch_id is required" in summary      # all problems at once


def test_empty_problems_file_is_fine(tmp_path):
    code, _s, _o = validate(tmp_path, project="P", problems=[])
    assert code == 0


# Same pattern as the workflow: sub(first_clade, "^I[ab]?(/.*)?$", "CLADE_I") == "CLADE_I"
CLADE_I_PATTERN = r"^I[ab]?(/.*)?$"


@pytest.mark.parametrize("clade,rerun", [
    ("I", True), ("Ia", True), ("Ib", True), ("Ib/sh2023", True),
    ("IIa", False), ("IIb", False), ("II", False), ("NA", False), ("", False),
    ("Ic", False), ("lineage-b.1", False),
])
def test_clade_i_trigger(clade, rerun):
    assert (re.sub(CLADE_I_PATTERN, "CLADE_I", clade) == "CLADE_I") is rerun


def test_workflow_uses_the_tested_clade_pattern():
    with open(os.path.join(HERE, "..", "workflows", "wf_masphl_viral_ngs.wdl")) as fh:
        assert f'sub(first_clade, "{CLADE_I_PATTERN}", "CLADE_I")' in fh.read()


def test_generated_theiacov_columns_up_to_date():
    assert gen.check(), "run: python tools/gen_theiacov_outputs.py"


def test_report_and_lims_columns_are_theiacov_outputs():
    outs = set(gen.load_outputs())
    needed = """assembly_length_unambiguous assembly_mean_coverage assembly_method
      consensus_n_variant_min_depth fastq_scan_num_reads_clean_pairs fastq_scan_num_reads_raw_pairs
      kraken_human kraken_human_dehosted kraken_sc2 kraken_sc2_dehosted kraken_target_organism
      kraken_target_organism_dehosted kraken_target_organism_name meanbaseq_trim meanmapq_trim
      nextclade_aa_dels nextclade_aa_subs nextclade_clade nextclade_docker nextclade_ds_tag
      nextclade_lineage nextclade_qc nextclade_version number_Degenerate number_N number_Total
      percent_reference_coverage primer_bed_name primer_trimmed_read_percent read_screen_clean
      read_screen_raw theiacov_illumina_pe_analysis_date theiacov_illumina_pe_version
      pango_lineage pangolin_assignment_version assembly_fasta read1_dehosted read2_dehosted""".split()
    assert [c for c in needed if c not in outs] == []
