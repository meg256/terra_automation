# Phase 1 testing: MASPHL_BaseSpace_Load

Test in the scratch workspace, writing to a test table (`sample_test`) so no
production rows are touched. Keep each launch's `check_summary.txt` - those are
the records to compare and send back.

## 0. One-time setup

1. Push this repo to the lab's GitHub org and register it on Dockstore (the
   `.dockstore.yml` registers `MASPHL_BaseSpace_Load`). Import it into the
   scratch workspace from Dockstore, using the development branch as the version.
2. In the scratch workspace's **Data > Workspace Data**, add the same two
   BaseSpace attributes the production workspaces use (`access_token`,
   `api_server`).
3. In the workflow's Inputs tab, set and **save** these defaults:

   | Input | Value |
   |---|---|
   | access_token | `workspace.access_token` |
   | api_server | `"https://api.basespace.illumina.com"` |
   | terra_project | first part of the workspace URL, e.g. `"masphl-sequencing"` |
   | terra_workspace | second part of the URL, with real spaces, e.g. `"Workflow DEV"` |
   | table_name | `"sample_test"` |
   | pathogen | change per test: `"SC2"`, `"MPXV"` or `"WNV"` |
   | basespace_project | per test (see note below) |

4. Choose **Run workflow with inputs defined by file paths**.

**basespace_project for now:** your current runs put the run name in
`Sample_Project`, so the BaseSpace project changes every run. Until the
pathogen-specific projects exist, enter it per run (it is the project name
shown on the run's Biosamples tab in BaseSpace). If you enter the wrong one,
the workflow stops and lists the projects that are in the run.

## 1. Smoke test (checks the Docker images and BaseSpace access)

Run one WNV run with `write_to_table = false`.

- [ ] Workflow succeeds.
- [ ] `basespace_cli_version` output shows a version (the BaseSpace image works).
- [ ] `check_summary` lists the expected samples and controls.
- [ ] `read_pairs` roughly matches what BaseSpace shows for each sample.

Three things I could only build from Illumina's documentation, not a live
run. If anything breaks, it is most likely one of these - send me the error:

- the run status text (the check expects `Complete`; if your finished runs
  say something else, the check stops and prints the actual value)
- the column names in the dataset listing (`Project.Name`, `DataSetType.Id`)
- the read-count attribute name (`TotalReadsPF`; if it's missing, samples are
  still fetched and the summary says "read count not available")

## 2. One run per pathogen, written to `sample_test`

For one past SC2, MPXV and WNV run each (set `btb_batch_id` for SC2):

- [ ] Same sample list as the hand-edited sample sheet used for that run.
- [ ] `batchid`: SC2 = the 5-digit BtB ID; MPXV/WNV = the `Sample_Project` value.
- [ ] `negative_control` = yes for the NEG/NTC samples only.
- [ ] `sample_sheet_id` matches what BaseSpace_Fetch used before.
- [ ] A sample set exists (SC2: the BtB ID; MPXV/WNV: the `Sample_Project` value)
      and contains every loaded sample.
- [ ] `read1` / `read2` are filled for every sample with reads.
- [ ] Warnings in `check_summary` make sense (e.g. a sample on the sheet that
      failed demultiplexing).

## 3. Same reads as BaseSpace_Fetch v1.4.1

For 2-3 samples already fetched in production with BaseSpace_Fetch v1.4.1,
compare file sizes of `read1`/`read2` (the data table's file links show size).
MiSeq samples have one file per read, so the files should be identical.

- [ ] Sizes match.

## 4. Re-launch the same run

Launch the same run again into `sample_test`.

- [ ] No duplicate rows; `replaced_rows` output is empty (same run = plain update).
- [ ] Faster second time (Terra reuses earlier results).

## 5. Repeat-sample overwrite

Without a real resequenced sample, fake one:

1. In `sample_test`, edit one sample's `basespace_run_id` to `999` and add a
   dummy column value to that row (e.g. `assembly_status = PASS`).
2. Re-launch the run.

- [ ] `replaced_rows` lists that sample with previous run `999`.
- [ ] The dummy `assembly_status` is gone; the row has the new load's values.

## 6. only_these_samples

- [ ] `"S1, S2, NTC1"` (three real IDs) loads exactly those three, nothing else.
- [ ] Adding a made-up ID (`"S1, FAKE99"`) fails and names `FAKE99`; nothing is written.
- [ ] Spaces instead of commas (`"S1 S2"`) fails with a clear message.

## 7. Inputs that should stop the run early

- [ ] Run name or a typo instead of the run number.
- [ ] SC2 with no `btb_batch_id`, and with a 4-digit one.
- [ ] Wrong `basespace_project` (the error lists the run's projects).
- [ ] A copy of the token with one character changed (clear authentication error).

## 8. extra_negative_control_ids

- [ ] Listing a regular sample as an extra control marks it `negative_control = yes`.

## Covered only by simulated tests (no real data yet)

NextSeq lane-split samples, runs with several projects (e.g. wastewater on the
same run), duplicate datasets from a requeued sample, and BaseSpace's `_` to
`-` renaming. See `tests/test_check_run_and_samples.py`. Run these against a
real example as soon as one exists.
