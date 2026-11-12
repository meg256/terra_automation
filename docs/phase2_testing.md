# Phase 2 testing: MASPHL_Viral_NGS (load + TheiaCoV)

Goal: the end-to-end workflow gives **the same TheiaCoV results** as today's
manual process, for all three pathogens, and the MPXV clade I rerun happens
on its own. Keep writing to `sample_test` in **Workflow DEV**.

## 0. Setup

1. Push, then in Dockstore check that `MASPHL_Viral_NGS` appears (the GitHub
   App Logs should show a successful push). Publish it and import it into
   Workflow DEV. `MASPHL_BaseSpace_Load` stays registered; you don't need it
   in the workspace anymore.
2. Choose **Run workflow with inputs defined by file paths**.
3. Inputs, then **Save**:
   - **Same as the Phase 1 config:** `access_token`, `api_server`,
     `terra_project` (`"masphl-sequencing"`), `terra_workspace`
     (`"Workflow DEV"`), `table_name` (`"sample_test"`), `pathogen`,
     `basespace_project`.
   - **Copied from the pathogen workspace's current TheiaCoV config**:
     `organism`, `reference_genome`, `primer_bed`, `nextclade_dataset_name`,
     `nextclade_dataset_tag`, and `reference_gff` /
     `reference_gene_locations_bed` / `genome_length` if you set them today.
     Workflow DEV needs the same `workspace.*` reference files as the
     production workspace (copy them under Data > Workspace Data if missing).
   - **Any other TheiaCoV setting you change today** (e.g. `min_depth`): set
     it in the Inputs tab under the task name `theiacov_illumina_pe`, as now.
     For MPXV, set it again under `theiacov_illumina_pe_clade_i`.
   - **MPXV only**: `mpxv_clade_i_reference_genome`,
     `mpxv_clade_i_primer_bed`, `mpxv_clade_i_nextclade_dataset_tag`
     (and the clade I gff / gene bed if your clade I config sets them).
     `mpxv_clade_i_nextclade_dataset_name` defaults to `nextstrain/mpox/clade-i`.

- [ ] In the Inputs tab, TheiaCoV's own settings (e.g.
      `theiacov_illumina_pe.min_depth`) are listed. If they are missing, tell
      me: it means Terra isn't showing nested settings, and I'll add the ones
      you use as top-level inputs.

## 1. Quick failure (seconds, not minutes)

- [ ] SC2 with no `btb_batch_id`: fails at `check_inputs` within about a
      minute, before BaseSpace is contacted.
- [ ] MPXV with `mpxv_clade_i_reference_genome` blank: fails the same way,
      naming the missing setting.

## 2. One past run per pathogen: same results as production

Pick runs already analyzed the current way, so production rows exist to
compare against (production `sample` table vs `sample_test`).

For a few samples per run, including a control and a FAIL-looking sample:

- [ ] `percent_reference_coverage`, `meanmapq_trim`, `meanbaseq_trim`
- [ ] `assembly_length_unambiguous`, `number_N`, `number_Total`
- [ ] `nextclade_clade`, `nextclade_lineage`, `nextclade_aa_subs`
- [ ] SC2: `pango_lineage`, `pangolin_assignment_version`
- [ ] `theiacov_illumina_pe_version` = PHB v4.0.0

Identical values mean the configuration matches. A difference almost always
means a TheiaCoV setting that differs from the production config. Send me the
column and the two values.

Also check:
- [ ] `analysis_route` is `standard` (SC2/WNV) or `mpxv_all_clades` (MPXV)
- [ ] `assembly_fasta`, `read1_dehosted`, `read2_dehosted` are filled (Mercury needs them)
- [ ] Controls with 0 reads (if any) have `load_status = no_reads` and no TheiaCoV columns

## 3. MPXV clade I

If a past clade I sample exists, load it with
`only_these_samples = "<that sample>, <a clade II sample>"`:

- [ ] clade I sample: `analysis_route = mpxv_clade_i_rerun`,
      `nextclade_clade_all_clades` = I/Ia/Ib, and the remaining columns match
      your manual clade I rerun (reference, `primer_bed_name`, Nextclade dataset).
- [ ] clade II sample: `analysis_route = mpxv_all_clades`.

With no clade I sample available, the trigger rule is covered by the simulated
tests (I, Ia, Ib and `Ib/...` trigger; IIa, IIb, II and NA don't), and the
rerun step itself waits for a real case.

## 4. Downstream workflows still work

The new workflow writes the same column names Terra's default output mapping
does, so existing workflows that read the table should run unchanged.

- [ ] SC2: run the current LIMS prep and gen on the new sample set in
      `sample_test` (choose `sample_test` / `sample_test_set` as the data
      table). The lims_file and run_results match production for that run.
- [ ] SC2: Mercury_PE_Prep runs on PASS samples from `sample_test` (stop
      before uploading anything).

## 5. Relaunch

- [ ] Relaunching the same run reuses earlier results (TheiaCoV shows
      "Cache hit") and finishes quickly, with no duplicate rows.

## Known behavior

- If TheiaCoV itself crashes on one sample, the whole run stops before the
  results are written to the table (each task already retries once). Relaunch
  after fixing; finished samples are reused from cache. Tell me if this
  happens often in practice.
- Columns that are empty for every sample in the run (flu-, HIV- or
  measles-only TheiaCoV outputs) are not written to the table.
