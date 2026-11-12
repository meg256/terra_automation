# Phase 3 testing: QC and result files

Same workflow (`MASPHL_Viral_NGS`), same setup as Phase 2, still writing to
`sample_test` in Workflow DEV.

## 0. Config changes from Phase 2

- [ ] **`organism` is no longer an input.** It is set from `pathogen`
      (SC2 -> `sars-cov-2`, MPXV, WNV). Terra may show it as a removed input; just save.
- [ ] **MPXV:** `primer_bed` = `workspace.MPXV_yale_primer_bed` (the all-clades /
      clade II bed, matching `reference_genome` = MT903345). Only the
      `mpxv_clade_i_*` inputs use the NC_003310 files. The workflow now refuses
      to start if the main and clade I primer bed (or reference) are the same file.
- New optional inputs, all defaulting to the SOP values:
  `min_percent_reference_coverage` (SC2 90, MPXV/WNV 80), `min_meanmapq` (20),
  `min_meanbaseq` (20), `control_max_unambiguous_bases` (600),
  `masphl_id_prefix` (`MPXV_USA_`).

## 1. MPXV: rerun MA26-32-Mpox-20260911 with the corrected primer bed

- [ ] `primer_bed_name` is the clade II bed for MA0271 and the NTC, and the
      clade I bed for MA0272/MA0273.
- [ ] MA0271's TheiaCoV values now match production.
- [ ] `assembly_status`: PASS/FAIL as your reviewer called it in production.
      (With the earlier, wrong primer bed: MA0271 PASS, MA0272 PASS, MA0273 FAIL
      on coverage 77.38, NTC FAIL.)
- [ ] `masphl_id` = `MPXV_USA_2026_<sample>`.
- [ ] Output `result_files`: `MA26-32-Mpox-20260911_results_<date>.tsv` has the
      SOP columns in SOP order, plus `analysis_route`,
      `nextclade_clade_all_clades`, `negative_control`, `qc_fail_reasons` at the end.
- [ ] Output `run_summary`: control section shows the NTC as clean; clade I
      reruns listed.

## 2. WNV and SC2: one past run each

- [ ] PASS/FAIL matches the production call for every sample.
- [ ] WNV: `<batchid>_results_<date>.tsv` has the WNV SOP columns.
- [ ] SC2: `<BtB ID>_<date>.run_results.csv` and `.lims_file.csv`. Compare with
      the production files for that run: identical except
      - organism `SARS-CoV-2` (was `SARS-CoV 2`)
      - no extra comma at the end of each run_results row
      - the date in the file name
      Values can also differ if the analysis date differs (`seq_date` is the
      TheiaCoV analysis date, as before).
- [ ] SC2 run summary lists FAIL samples under "Check Ct values in BtB".

## 3. Contamination alert

Without a real contaminated control, trust the simulated test
(`tests/test_qc_and_reports.py`). In a real run with one:

- [ ] `run_summary` starts with a `CONTAMINATION ALERT -- notify the Laboratory
      Supervisor` block naming the control.
- [ ] Output `contamination_alert` = true; sample set column `contamination_alert` = yes.

## 4. Sample set row

- [ ] The run's sample set (`sample_test_set`) has `run_summary`,
      `contamination_alert`, and `run_report` (MPXV/WNV) or `run_results_file` +
      `btb_lims_file` (SC2, same names as the old LIMS gen outputs), each a
      clickable file link.

## 5. Downstream

- [ ] SC2: the validated Mercury_PE_Prep still runs on PASS samples from the
      table (`assembly_status` is written exactly as before).
