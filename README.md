# MASPHL viral NGS workflows (SC2 / MPXV / WNV)

End-to-end Terra workflow for MASPHL's viral amplicon runs, with each part testable as separate dockstore workflow

| Phase | Workflow | Status |
|---|---|---|
| 1 | `MASPHL_BaseSpace_Load` - check run, register samples, fetch reads | **ready for testing** |
| 2 | Analysis - TheiaCoV v4.0.0 (untouched, imported by tag) + MPXV clade I rerun | planned |
| 3 | QC + reports (replaces SPHL LIMS prep/gen for SC2; reports for MPXV/WNV) | planned |
| - | Master workflow chaining 1-3 | planned |

Mercury (GISAID/GenBank prep) is excluded from this end-to-end workflow; run the validated standalone Mercury workflows afterwards as today.

## Repository layout

```
.dockstore.yml                    registers each workflow with Dockstore
workflows/wf_basespace_load.wdl   initial workflow
tasks/task_basespace.wdl          get_basespace_run_info, check_run_and_samples, fetch_basespace_reads
tasks/task_terra_table.wdl        update_terra_table, make_read_rows
scripts/check_run_and_samples.py  the loader's decision logic (source of truth)
tools/embed_scripts.py            copies scripts/ into the WDL tasks
tests/                            simulated BaseSpace runs (pytest) + example inputs
docs/phase1_testing.md            checklist for testing on real runs
```

## Editing the Python logic

The tasks use plain public images, so the script travels inside the WDL.
Edit `scripts/check_run_and_samples.py`, then:

```
python tools/embed_scripts.py        # refresh the copy inside tasks/task_basespace.wdl
python -m pytest tests/ -q           # fails if the embedded copy is stale
miniwdl check workflows/wf_basespace_load.wdl
```

## Docker images (all public; pin by digest before validation)

| Task | Image | Needs |
|---|---|---|
| get_basespace_run_info, fetch_basespace_reads | `us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1` | `bs`, bash, gzip, tar |
| check_run_and_samples, make_read_rows | `python:3.12-slim` | python3 (standard library only) |
| update_terra_table | `python:3.12-slim` + `firecloud==0.16.41` installed at run time | python3, internet access to PyPI |

Each image is a task input (`docker`), so it can be swapped without code changes.
`update_terra_table` skips the install if the image already has FISS.

## Attribution

Phase 1 replaces Theiagen's BaseSpace_Fetch v1.4.1 (terra_utilities) for this
lab's use; the fetch step follows the same approach (BaseSpace CLI download,
join lane files). TheiaCoV is Theiagen's Public Health Bioinformatics (PHB)
and is imported unmodified.
