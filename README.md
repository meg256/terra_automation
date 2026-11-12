# MASPHL viral NGS workflows (SC2 / MPXV / WNV)

End-to-end Terra workflow for MASPHL's viral amplicon runs, built in phases.
Labs import one workflow, `MASPHL_Viral_NGS`; each phase adds to it.

| Phase | What it adds | Status |
|---|---|---|
| 1 | Load: check run, register samples, fetch reads (`wf_basespace_load.wdl`) | tested on SC2/MPXV/WNV runs |
| 2 | TheiaCoV v4.0.0 per sample (imported unmodified, pinned by tag) + MPXV clade I rerun; results to the data table | **ready for testing** |
| 3 | QC PASS/FAIL + reports (replaces SPHL LIMS prep/gen for SC2; reports for MPXV/WNV) | planned |

Mercury (GISAID/GenBank prep) is out of scope; run the validated standalone
Mercury workflows afterwards as today.

## Repository layout

```
.dockstore.yml                    registers each workflow with Dockstore
workflows/wf_masphl_viral_ngs.wdl  the end-to-end workflow (import this one)
workflows/wf_basespace_load.wdl    loading step (sub-workflow; also launchable alone)
tasks/task_basespace.wdl           get_basespace_run_info, check_run_and_samples, fetch_basespace_reads
tasks/task_terra_table.wdl         update_terra_table, make_read_rows, make_result_rows
scripts/check_run_and_samples.py   the loader's decision logic (source of truth)
config/theiacov_v4.0.0_outputs.tsv TheiaCoV v4.0.0 output list (drives the table columns)
tools/embed_scripts.py             copies scripts/ into the WDL tasks
tools/gen_theiacov_outputs.py      writes the TheiaCoV-output -> table-column WDL
tests/                             simulated runs and checks (pytest) + example inputs
docs/phase1_testing.md, docs/phase2_testing.md   checklists for testing on real runs
```

## Editing the Python logic

The tasks use plain public images, so the script travels inside the WDL.
Edit `scripts/check_run_and_samples.py`, then:

```
python tools/embed_scripts.py          # refresh the copy inside tasks/task_basespace.wdl
python tools/gen_theiacov_outputs.py   # refresh the TheiaCoV column wiring
python -m pytest tests/ -q             # fails if either is stale
miniwdl check workflows/wf_masphl_viral_ngs.wdl
```

## Updating TheiaCoV later

1. Change the version in the `import` URL in `wf_masphl_viral_ngs.wdl`.
2. Regenerate the output list from that tag and the column wiring:
   `python tools/gen_theiacov_outputs.py --from-wdl <checkout>/workflows/theiacov/wf_theiacov_illumina_pe.wdl --version 4.x.y`,
   set `VERSION` in that tool, then run it again without arguments.
3. Revalidate (docs/phase2_testing.md, section 2).

## Docker images (all public; pin by digest before validation)

| Task | Image | Needs |
|---|---|---|
| get_basespace_run_info, fetch_basespace_reads | `us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1` | `bs`, bash, gzip, tar |
| check_run_and_samples, make_read_rows, make_result_rows | `python:3.12-slim` | python3 (standard library only) |
| update_terra_table | `python:3.12-slim` + `firecloud==0.16.41` installed at run time | python3, internet access to PyPI |

Each image is a task input (`docker`), so it can be swapped without code changes.
`update_terra_table` skips the install if the image already has FISS.

## Attribution

Phase 1 replaces Theiagen's BaseSpace_Fetch v1.4.1 (terra_utilities) for this
lab's use; the fetch step follows the same approach (BaseSpace CLI download,
join lane files). TheiaCoV is Theiagen's Public Health Bioinformatics (PHB)
and is imported unmodified.

TheiaCoV_Illumina_PE_PHB v4.0.0 is called directly from Theiagen's repository by
version tag; no copy of it is kept here.
