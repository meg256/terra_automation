version 1.0

# =============================================================================
# MASPHL BaseSpace Load  (Phase 1 of the MASPHL viral NGS end-to-end workflow)
# -----------------------------------------------------------------------------
# Replaces the manual sample sheet edit + TSV upload + BaseSpace_Fetch steps for
# SC2, MPXV and WNV. Launch in Terra with
#     "Run workflow with inputs defined by file paths"   (no data table rows)
#
# What it does
#   1. get_basespace_run_info  snapshot of the run (datasets, read counts, sheet)
#   2. check_run_and_samples   decide what to load; stop early on bad inputs
#   3. register_samples        write sample rows + sample set to the data table
#   4. fetch_basespace_reads   download FASTQs per sample, join lanes
#   5. write_reads             write read1/read2 paths back to the rows
#
# Inputs users fill in per run:  basespace_run_id, btb_batch_id (SC2 only),
#   and optionally only_these_samples / extra_negative_control_ids.
# Everything else is a saved default in each pathogen workspace.
# =============================================================================

import "../tasks/task_basespace.wdl" as basespace
import "../tasks/task_terra_table.wdl" as terra_table

workflow basespace_load {
  meta {
    description: "Load one BaseSpace run into a Terra workspace: check the run, register samples, fetch reads."
  }
  parameter_meta {
    basespace_run_id: "Number from the run's BaseSpace URL: basespace.illumina.com/run/<number>/details"
    btb_batch_id: "SC2 only: 5-digit BtB batch ID. Becomes batchid and the sample set name."
    only_these_samples: "Optional. Comma-separated sample IDs to load, e.g. 'MA01, MA02, NTC01'. Replaces the project filter; samples may come from any project in the run; controls are NOT added automatically; any ID not found fails the run."
    extra_negative_control_ids: "Optional. Comma-separated IDs of negative controls whose names don't start with NEG or NTC."
    pathogen: "Workspace default: SC2, MPXV or WNV."
    basespace_project: "Workspace default (or per run): BaseSpace project holding this pathogen's samples. Ignored when only_these_samples is set."
    table_name: "Data table to write to. Use e.g. 'sample_test' while testing."
    write_to_table: "false = check and fetch only, write nothing to the data table."
  }

  input {
    # ---- per run ---------------------------------------------------------
    String basespace_run_id
    String? btb_batch_id
    String? only_these_samples
    String? extra_negative_control_ids

    # ---- workspace defaults -----------------------------------------------
    String pathogen
    String? basespace_project
    String access_token
    String api_server = "https://api.basespace.illumina.com"
    String terra_project
    String terra_workspace

    # ---- rarely changed ----------------------------------------------------
    String table_name = "sample"
    Boolean write_to_table = true
    String? sample_set_name
    String control_id_regex = "^(NEG|NTC)"
    Boolean allow_incomplete_run = false
    Boolean clear_old_columns_on_rerun = true
    String loader_version = "phase2-dev"
    # Problems found by a calling workflow (e.g. missing MPXV clade I settings);
    # reported together with the loader's own input checks.
    Array[String] other_input_problems = []
  }

  # Seconds-long input check, so a missing BtB ID etc. fails before BaseSpace is touched
  call basespace.check_run_and_samples as check_inputs {
    input:
      validate_only = true,
      other_input_problems = other_input_problems,
      pathogen = pathogen,
      basespace_run_id = basespace_run_id,
      btb_batch_id = btb_batch_id,
      basespace_project = basespace_project,
      only_these_samples = only_these_samples,
      extra_negative_control_ids = extra_negative_control_ids,
      control_id_regex = control_id_regex,
      sample_set_name = sample_set_name,
      table_name = table_name,
      loader_version = loader_version
  }

  call basespace.get_basespace_run_info {
    input:
      wait_for = check_inputs.check_summary,
      basespace_run_id = basespace_run_id,
      api_server = api_server,
      access_token = access_token
  }

  call basespace.check_run_and_samples {
    input:
      snapshot = get_basespace_run_info.snapshot,
      pathogen = pathogen,
      basespace_run_id = basespace_run_id,
      btb_batch_id = btb_batch_id,
      basespace_project = basespace_project,
      only_these_samples = only_these_samples,
      extra_negative_control_ids = extra_negative_control_ids,
      control_id_regex = control_id_regex,
      allow_incomplete_run = allow_incomplete_run,
      sample_set_name = sample_set_name,
      table_name = table_name,
      loader_version = loader_version
  }

  if (write_to_table) {
    call terra_table.update_terra_table as register_samples {
      input:
        terra_project = terra_project,
        terra_workspace = terra_workspace,
        table_name = table_name,
        rows_tsv = check_run_and_samples.sample_rows_tsv,
        membership_tsv = check_run_and_samples.membership_tsv,
        check_previous_run = true,
        clear_old_columns_on_rerun = clear_old_columns_on_rerun
    }
  }

  scatter (row in check_run_and_samples.fetch_table) {
    call basespace.fetch_basespace_reads {
      input:
        sample_id = row[0],
        dataset_ids = row[1],
        api_server = api_server,
        access_token = access_token
    }
    String fetched_id = row[0]
    String read1_path = fetch_basespace_reads.read1
    String read2_path = fetch_basespace_reads.read2
  }

  if (write_to_table) {
    call terra_table.make_read_rows {
      input:
        table_name = table_name,
        sample_ids = fetched_id,
        read1 = read1_path,
        read2 = read2_path,
        read_pairs = fetch_basespace_reads.read_pairs
    }
    # wait_for: must run after register_samples, which may clear old columns
    call terra_table.update_terra_table as write_reads {
      input:
        terra_project = terra_project,
        terra_workspace = terra_workspace,
        table_name = table_name,
        rows_tsv = make_read_rows.read_rows_tsv,
        check_previous_run = false,
        wait_for = register_samples.update_log
    }
  }

  output {
    File check_summary = check_run_and_samples.check_summary
    File check_report = check_run_and_samples.check_report
    String sample_set = check_run_and_samples.sample_set
    File basespace_snapshot = get_basespace_run_info.snapshot
    String basespace_cli_version = get_basespace_run_info.bs_cli_version
    Array[String] loaded_sample_ids = fetched_id
    Array[File] read1 = fetch_basespace_reads.read1
    Array[File] read2 = fetch_basespace_reads.read2
    Array[Int] read_pairs = fetch_basespace_reads.read_pairs
    File? replaced_rows = register_samples.replaced_rows
    File? table_update_log = register_samples.update_log
    File? read_paths_update_log = write_reads.update_log
  }
}
