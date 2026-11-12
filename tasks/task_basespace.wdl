version 1.0

# BaseSpace tasks for the MASPHL viral NGS loader.
#
#   get_basespace_run_info  - bash + BaseSpace CLI only. Takes a snapshot of the
#                             run: run details, dataset list, read counts and
#                             sample sheet. No decisions are made here.
#   check_run_and_samples   - python3 standard library only. Decides what to load
#                             from the snapshot (logic in scripts/check_run_and_samples.py,
#                             embedded below by tools/embed_scripts.py).
#   fetch_basespace_reads   - bash + BaseSpace CLI. Downloads one sample's FASTQ
#                             dataset(s) by exact dataset ID and joins lanes.
#
# The BaseSpace access token is written to a config file and never echoed to
# stdout/stderr. Note that Cromwell still stores each task's command script
# (which contains the token) in the workspace bucket, as BaseSpace_Fetch does
# today, so keep workspace sharing limited to lab staff.

task get_basespace_run_info {
  input {
    File? wait_for                  # only used to order calls (runs after the input check)
    String basespace_run_id
    String api_server
    String access_token
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"
    Int disk_size = 20
  }
  command <<<
    set -euo pipefail
    if ! command -v bs >/dev/null 2>&1; then
      echo "ERROR: the docker image '~{docker}' has no BaseSpace CLI ('bs'). Use an image that does." >&2
      exit 1
    fi

    # ---- authenticate (token goes to a config file, never to the log) ------
    set +x
    export HOME="$PWD/bs_home"
    mkdir -p "$HOME/.basespace"
    printf 'apiServer = %s\naccessToken = %s\n' '~{api_server}' '~{access_token}' > "$HOME/.basespace/default.cfg"
    chmod 600 "$HOME/.basespace/default.cfg"

    mkdir -p snapshot/attributes snapshot/run_files
    { bs --version 2>&1 || true; } > snapshot/bs_version.txt
    echo "BaseSpace CLI: $(cat snapshot/bs_version.txt)"

    if ! bs whoami -F Name -F Host -f csv > snapshot/whoami.csv 2> whoami.err; then
      echo "ERROR: BaseSpace did not accept the access token (check workspace attribute / input access_token and api_server)." >&2
      cat whoami.err >&2
      exit 1
    fi
    echo "Authenticated to BaseSpace as: $(tail -n +2 snapshot/whoami.csv)"

    # ---- run details ------------------------------------------------------
    if ! bs get run -i '~{basespace_run_id}' -f json > snapshot/run.json 2> run.err; then
      echo "ERROR: could not open BaseSpace run '~{basespace_run_id}'. Use the number from the run URL" >&2
      echo "       (basespace.illumina.com/run/<number>/details) and check the token's account can see the run." >&2
      cat run.err >&2
      exit 1
    fi

    # ---- every dataset in the run (all projects) ---------------------------
    bs list dataset --input-run='~{basespace_run_id}' \
      -F Id -F Name -F Project.Name -F Project.Id -F DataSetType.Id -F DateCreated \
      -f csv > snapshot/datasets.csv
    echo "Datasets in run: $(( $(wc -l < snapshot/datasets.csv) - 1 ))"

    # ---- read counts per dataset (used to spot empty samples) --------------
    tail -n +2 snapshot/datasets.csv | cut -d, -f1 | tr -d '"\r' | while read -r ds_id; do
      [ -n "$ds_id" ] || continue
      bs list attributes dataset -i "$ds_id" -f csv > "snapshot/attributes/${ds_id}.csv" 2>/dev/null \
        || rm -f "snapshot/attributes/${ds_id}.csv"
    done

    # ---- sample sheet (optional cross-check) -------------------------------
    # --extension=csv works on old CLI versions (Theiagen's image has bs 1.2.1,
    # which has no --include/--exclude). Only the run's small .csv files come down.
    # The log and the CLI's own help are kept in the snapshot for troubleshooting.
    mkdir -p snapshot/logs
    bs download run --help > snapshot/logs/bs_download_run_help.txt 2>&1 || true
    if bs download run -i '~{basespace_run_id}' -o snapshot/run_files --extension=csv \
         > snapshot/logs/samplesheet_download.log 2>&1; then
      echo "Run .csv files downloaded: $(find snapshot/run_files -type f -name '*.csv' | wc -l)"
      find snapshot/run_files -type f -iname 'samplesheet*.csv' | sed 's/^/  sample sheet: /'
    else
      echo "WARNING: could not download the run's sample sheet; the cross-check will be skipped."
      cat snapshot/logs/samplesheet_download.log
    fi
    # keep only sample sheets (other run csvs are not needed downstream)
    find snapshot/run_files -type f ! -iname 'samplesheet*.csv' -delete 2>/dev/null || true

    tar czf basespace_snapshot.tar.gz snapshot
  >>>
  output {
    File snapshot = "basespace_snapshot.tar.gz"
    String bs_cli_version = read_string("snapshot/bs_version.txt")
  }
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 1
    maxRetries: 1
  }
}

task check_run_and_samples {
  input {
    File? snapshot                  # not needed when validate_only = true
    Boolean validate_only = false   # quick input check before touching BaseSpace
    Array[String] other_input_problems = []
    String pathogen
    String basespace_run_id
    String? btb_batch_id
    String? basespace_project
    String? only_these_samples
    String? extra_negative_control_ids
    String control_id_regex = "^(NEG|NTC)"
    Boolean allow_incomplete_run = false
    String? sample_set_name
    String table_name = "sample"
    String loader_version
    String docker = "python:3.12-slim"
  }
  command <<<
    set -euo pipefail
    if [ '~{validate_only}' = 'false' ]; then
      tar xzf '~{snapshot}'
    fi

    cat > check_run_and_samples.py <<'PYEOF'
    # BEGIN EMBEDDED scripts/check_run_and_samples.py (generated by tools/embed_scripts.py; do not edit here)
    #!/usr/bin/env python3
    """
    check_run_and_samples.py -- MASPHL viral NGS loader, step 2 of 4.
    
    Reads a snapshot of one BaseSpace run (taken by the get_basespace_run_info
    task) and decides which samples to load into Terra.
    
      * checks the inputs (pathogen, BtB batch ID, project / sample list)
      * confirms the run finished and the target project is in it
      * matches each sample to its FASTQ dataset(s), joining lane-split datasets
      * matches BaseSpace names back to the sample sheet, allowing for BaseSpace
        turning "_" into "-"
      * marks negative controls (NEG/NTC prefix, or listed by the user)
      * flags problems: duplicate datasets, empty samples, sample-sheet mismatches,
        requested IDs that are not in the run
    
    Standard library only, so it runs in any python3 image and can be tested
    offline (see tests/). A copy of this file is embedded in
    tasks/task_basespace.wdl by tools/embed_scripts.py -- edit THIS file, then
    re-run the embed tool. Do not use the WDL placeholder characters (tilde-brace)
    or three closing angle brackets anywhere in this file; the embed tool refuses.
    
    Snapshot directory layout (written by get_basespace_run_info):
      run.json                 output of: bs get run -i RUN -f json
      datasets.csv             output of: bs list dataset --input-run=RUN -F ... -f csv
      attributes/<id>.csv      output of: bs list attributes dataset -i <id> -f csv
      run_files/**             SampleSheet*.csv downloaded from the run (optional)
    
    Exit code 0 = OK (warnings may exist), 1 = check failed (reasons printed).
    """
    
    import argparse
    import csv
    import json
    import os
    import re
    import sys
    from collections import OrderedDict
    
    SCRIPT_VERSION = "0.2.0"
    
    PATHOGENS = ("SC2", "MPXV", "WNV")
    LANE_RE = re.compile(r"_L(\d{1,3})$")
    ENTITY_NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
    BTB_BATCH_RE = re.compile(r"^\d{5}$")
    RUN_ID_RE = re.compile(r"^\d+$")
    DEFAULT_CONTROL_REGEX = r"^(NEG|NTC)"
    
    
    class CheckFailed(Exception):
        """Raised with a list of human-readable reasons."""
    
        def __init__(self, reasons):
            super().__init__("; ".join(reasons))
            self.reasons = reasons
    
    
    # --------------------------------------------------------------------------
    # Small helpers
    # --------------------------------------------------------------------------
    
    def norm(name):
        """Key used to compare names: BaseSpace turns '_' into '-' in sample names."""
        return name.strip().replace("_", "-")
    
    
    def parse_id_list(text, label):
        """'MA01,MA02, NTC01' -> ['MA01', 'MA02', 'NTC01'] (order kept, dups dropped)."""
        if text is None or not text.strip():
            return []
        out = []
        for part in text.split(","):
            item = part.strip()
            if not item:
                continue
            if re.search(r"\s", item):
                raise CheckFailed([
                    f"{label}: '{item}' contains a space. Separate sample IDs with commas, "
                    f"e.g. 'MA01, MA02, NTC01'."
                ])
            if item not in out:
                out.append(item)
        return out
    
    
    def first_present(record, *keys, default=""):
        for k in keys:
            v = record.get(k)
            if v not in (None, ""):
                return v
        return default
    
    
    # --------------------------------------------------------------------------
    # Snapshot readers
    # --------------------------------------------------------------------------
    
    def load_run(path):
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        # `bs get run -f json` may return the object directly or wrapped in a list
        if isinstance(data, list):
            data = data[0] if data else {}
        return {
            "id": str(first_present(data, "Id", "id")),
            "name": first_present(data, "Name", "name"),
            "experiment_name": first_present(data, "ExperimentName", "experimentName"),
            "status": first_present(data, "Status", "status"),
            "date_created": first_present(data, "DateCreated", "dateCreated"),
        }
    
    
    def load_datasets(path):
        datasets = []
        with open(path, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                datasets.append({
                    "id": first_present(row, "Id"),
                    "name": first_present(row, "Name"),
                    "project": first_present(row, "Project.Name", "ProjectName"),
                    "project_id": first_present(row, "Project.Id", "ProjectId"),
                    "type": first_present(row, "DataSetType.Id", "DataSetTypeId"),
                    "date_created": first_present(row, "DateCreated"),
                })
        return datasets
    
    
    def load_read_counts(attr_dir, dataset_ids):
        """dataset id -> TotalReadsPF (int) or None when unknown."""
        counts = {}
        for ds_id in dataset_ids:
            counts[ds_id] = None
            path = os.path.join(attr_dir, f"{ds_id}.csv")
            if not os.path.exists(path):
                continue
            with open(path, newline="", encoding="utf-8-sig") as fh:
                for row in csv.DictReader(fh):
                    if first_present(row, "Name") == "TotalReadsPF":
                        try:
                            counts[ds_id] = int(float(first_present(row, "Value")))
                        except ValueError:
                            pass
        return counts
    
    
    def find_samplesheet(run_files_dir):
        """Pick the run's sample sheet: prefer a file named exactly SampleSheet.csv,
        closest to the top of the run folder."""
        found = []
        if run_files_dir and os.path.isdir(run_files_dir):
            for root, _dirs, files in os.walk(run_files_dir):
                for f in files:
                    if f.lower().startswith("samplesheet") and f.lower().endswith(".csv"):
                        p = os.path.join(root, f)
                        exact = f == "SampleSheet.csv"
                        found.append((0 if exact else 1, p.count(os.sep), p))
        found.sort()
        return [p for _a, _b, p in found]
    
    
    def parse_samplesheet(path):
        """Parse a v1 ([Data]) or v2 ([BCLConvert_Data] / [Cloud_Data]) sample sheet.
    
        Returns a list of {'sample_id', 'sample_name', 'project'} in sheet order.
        Columns for the same Sample_ID are merged across sections.
        """
        sections = OrderedDict()
        current = None
        with open(path, newline="", encoding="utf-8-sig") as fh:
            for row in csv.reader(fh):
                cells = [c.strip() for c in row]
                if not any(cells):
                    continue
                first = cells[0]
                if first.startswith("[") and first.endswith("]"):
                    current = first[1:-1]
                    sections[current] = []
                    continue
                if current is not None:
                    sections[current].append(cells)
    
        merged = OrderedDict()
        for lines in sections.values():
            if not lines or "Sample_ID" not in lines[0]:
                continue
            header = lines[0]
            sid_idx = header.index("Sample_ID")
            for cells in lines[1:]:
                if sid_idx >= len(cells) or not cells[sid_idx]:
                    continue
                rec = merged.setdefault(cells[sid_idx], {})
                for h, v in zip(header, cells):
                    if h and v and h not in rec:
                        rec[h] = v
    
        out = []
        for sid, rec in merged.items():
            out.append({
                "sample_id": sid,
                "sample_name": first_present(rec, "Sample_Name", "LibraryName", default=sid),
                "project": first_present(rec, "Sample_Project", "ProjectName"),
            })
        return out
    
    
    # --------------------------------------------------------------------------
    # Core logic (pure; unit-tested)
    # --------------------------------------------------------------------------
    
    def validate_inputs(opts):
        reasons = list(opts.get("other_input_problems", []))
        if opts["pathogen"] not in PATHOGENS:
            reasons.append(f"pathogen must be one of {', '.join(PATHOGENS)} (got '{opts['pathogen']}').")
        if not RUN_ID_RE.match(opts["run_id"] or ""):
            reasons.append(
                f"basespace_run_id must be the number from the run's BaseSpace URL, e.g. "
                f"basespace.illumina.com/run/332775445/details -> 332775445 (got '{opts['run_id']}')."
            )
        if opts["pathogen"] == "SC2":
            if not opts.get("btb_batch_id"):
                reasons.append("btb_batch_id is required for SC2 runs (the 5-digit BtB batch ID).")
            elif not BTB_BATCH_RE.match(opts["btb_batch_id"]):
                reasons.append(f"btb_batch_id must be exactly 5 digits (got '{opts['btb_batch_id']}').")
        if not opts["only_these_samples"] and not opts.get("project"):
            reasons.append(
                "basespace_project is empty. Set the workspace default, or list the samples "
                "to load in only_these_samples."
            )
        try:
            re.compile(opts["control_regex"])
        except re.error as e:
            reasons.append(f"control_id_regex is not a valid pattern: {e}")
        if opts.get("sample_set_name") and not ENTITY_NAME_RE.match(opts["sample_set_name"]):
            reasons.append(
                f"sample_set_name '{opts['sample_set_name']}' may only contain letters, numbers, "
                f"'.', '_' and '-'."
            )
        if reasons:
            raise CheckFailed(reasons)
    
    
    def group_datasets(datasets):
        """Group FASTQ datasets into samples: same project + same name once lane
        suffixes are removed and '_'/'-' are treated alike."""
        groups = OrderedDict()
        skipped = []
        for d in datasets:
            if not d["type"].lower().startswith("illumina.fastq"):
                continue
            if d["name"].lower().startswith("undetermined"):
                skipped.append(dict(d, status="undetermined", note="Undetermined reads; never loaded."))
                continue
            m = LANE_RE.search(d["name"])
            lane = int(m.group(1)) if m else None
            base = d["name"][: m.start()] if m else d["name"]
            key = (d["project"], norm(base))
            g = groups.setdefault(key, {"project": d["project"], "key": norm(base),
                                        "base_names": [], "datasets": []})
            if base not in g["base_names"]:
                g["base_names"].append(base)
            g["datasets"].append(dict(d, lane=lane, base=base))
    
        for g in groups.values():
            names = [x["name"] for x in g["datasets"]]
            lanes = [x["lane"] for x in g["datasets"]]
            g["problem"] = None
            if len(names) != len(set(names)):
                g["problem"] = ("duplicate_datasets",
                                "More than one dataset with the same name (requeued sample?). "
                                "Remove the extra dataset in BaseSpace, or load it with only_these_samples "
                                "after cleaning up.")
            elif None in lanes and any(x is not None for x in lanes):
                g["problem"] = ("duplicate_datasets",
                                "Both a lane-split and an un-split dataset exist for this sample.")
            elif len(g["base_names"]) > 1:
                g["problem"] = ("name_collision",
                                f"Different names become identical once '_' and '-' are treated alike: "
                                f"{', '.join(g['base_names'])}.")
            g["datasets"].sort(key=lambda x: (x["lane"] or 0, x["date_created"], x["name"]))
        return groups, skipped
    
    
    def run_check(run, datasets, sheet_rows, read_counts, opts):
        """Return a result dict. Raises CheckFailed with all fatal reasons."""
        validate_inputs(opts)
        failures, warnings, notes = [], [], []
    
        # ---- run status -------------------------------------------------------
        if run["id"] and run["id"] != opts["run_id"]:
            failures.append(f"BaseSpace returned run {run['id']} for requested run {opts['run_id']}.")
        if run["status"] != "Complete":
            msg = (f"Run status in BaseSpace is '{run['status']}', not 'Complete'. FASTQs may be "
                   f"missing or incomplete.")
            if opts["allow_incomplete_run"]:
                warnings.append(msg + " Continuing because allow_incomplete_run = true.")
            else:
                failures.append(msg + " Wait for the run to finish, or set allow_incomplete_run = true "
                                      "if you are sure the FASTQs are final.")
    
        groups, skipped = group_datasets(datasets)
        if not groups:
            raise CheckFailed(failures + [
                f"No FASTQ datasets found in run {opts['run_id']}. Has FASTQ generation finished, and "
                f"can this BaseSpace token see the run's projects?"
            ])
        projects_in_run = sorted({g["project"] for g in groups.values()})
    
        # ---- sample sheet index ----------------------------------------------
        sheet_by_key = {}
        for rec in sheet_rows:
            for k in {norm(rec["sample_id"]), norm(rec["sample_name"])}:
                sheet_by_key.setdefault(k, rec)
    
        # ---- choose samples -----------------------------------------------------
        only = opts["only_these_samples"]
        chosen = []          # list of (group, typed_id or None)
        if only:
            notes.append("only_these_samples was given: basespace_project is ignored and only the "
                         "listed samples are loaded (controls are not added automatically).")
            missing, ambiguous, broken = [], [], []
            for typed in only:
                keys = {norm(typed)}
                rec = sheet_by_key.get(norm(typed))
                if rec:
                    keys.add(norm(rec["sample_id"]))
                cands = [g for g in groups.values() if g["key"] in keys]
                if not cands:
                    missing.append(typed)
                elif len(cands) > 1:
                    ambiguous.append(f"{typed} (projects: {', '.join(sorted(g['project'] for g in cands))})")
                elif cands[0]["problem"]:
                    broken.append(f"{typed} ({cands[0]['problem'][1]})")
                elif any(c is cands[0] for c, _t in chosen):
                    continue
                else:
                    chosen.append((cands[0], typed))
            if missing:
                failures.append(f"These samples from only_these_samples are not in run {opts['run_id']}: "
                                f"{', '.join(missing)}.")
            if ambiguous:
                failures.append("These samples appear in more than one project of this run, so the "
                                "workflow can't tell which to load: " + "; ".join(ambiguous) + ".")
            if broken:
                failures.append("These listed samples can't be loaded as-is: " + "; ".join(broken))
        else:
            project = opts["project"]
            in_project = [g for g in groups.values() if g["project"] == project]
            if not in_project:
                failures.append(f"Project '{project}' has no FASTQ datasets in run {opts['run_id']}. "
                                f"Projects in this run: {', '.join(projects_in_run) or 'none'}.")
            for g in in_project:
                chosen.append((g, None))
    
        chosen_ids = {id(g) for g, _t in chosen}
    
        # ---- build sample records ---------------------------------------------
        samples, report = [], []
        for d in skipped:
            report.append({"basespace_name": d["name"], "project": d["project"], "dataset_ids": d["id"],
                           "lanes": "", "reads_pf": "", "sample_id": "", "status": d["status"],
                           "note": d["note"]})
    
        batchid_fallbacks = 0
        for g in groups.values():
            ds = g["datasets"]
            ds_ids = [x["id"] for x in ds]
            lanes = ",".join(str(x["lane"]) for x in ds if x["lane"] is not None)
            counts = [read_counts.get(i) for i in ds_ids]
            reads = None if any(c is None for c in counts) else sum(counts)
            base_row = {"basespace_name": "|".join(x["name"] for x in ds), "project": g["project"],
                        "dataset_ids": ",".join(ds_ids), "lanes": lanes,
                        "reads_pf": "" if reads is None else str(reads)}
    
            if id(g) not in chosen_ids:
                status = "not_requested" if only else "other_project"
                note = "Not in only_these_samples." if only else f"In project '{g['project']}', not the target project."
                report.append(dict(base_row, sample_id="", status=status, note=note))
                continue
    
            if g["problem"]:
                status, note = g["problem"]
                warnings.append(f"{'/'.join(g['base_names'])}: {note} Not loaded.")
                report.append(dict(base_row, sample_id="", status=status, note=note))
                continue
    
            typed = next((t for gg, t in chosen if gg is g), None)
            rec = sheet_by_key.get(g["key"])
            sample_id = rec["sample_id"] if rec else (typed or g["base_names"][0])
            sample_sheet_id = rec["sample_name"] if rec else g["base_names"][0]
    
            if opts["pathogen"] == "SC2":
                batchid = opts["btb_batch_id"]
            elif rec and rec["project"]:
                batchid = rec["project"]
            else:
                batchid = g["project"]
                batchid_fallbacks += 1
    
            is_control = bool(re.match(opts["control_regex"], sample_id, re.IGNORECASE)) or \
                norm(sample_id) in {norm(x) for x in opts["extra_negative_control_ids"]}
    
            if reads == 0:
                load_status = "no_reads"
                note = "BaseSpace reports 0 reads; added to the table but not fetched."
                (notes if is_control else warnings).append(f"{sample_id}: 0 reads in BaseSpace; not fetched.")
            else:
                load_status = "loaded"
                note = "" if reads is not None else "Read count not available from BaseSpace."
            if sheet_rows and not rec:
                warnings.append(f"{sample_id}: has FASTQs in BaseSpace but is not on the sample sheet.")
                note = (note + " Not on sample sheet.").strip()
    
            samples.append({
                "sample_id": sample_id, "sample_sheet_id": sample_sheet_id, "batchid": batchid,
                "project": g["project"], "dataset_ids": ds_ids,
                "dataset_names": [x["name"] for x in ds], "reads_pf": reads,
                "negative_control": is_control, "load_status": load_status,
            })
            report.append(dict(base_row, sample_id=sample_id, status=load_status, note=note))
    
        # ---- sample sheet rows with no data ------------------------------------
        if sheet_rows and not only:
            has_projects = any(r["project"] for r in sheet_rows)
            group_keys = {g["key"] for g in groups.values() if g["project"] == opts.get("project")}
            for r in sheet_rows:
                if has_projects and r["project"] != opts.get("project"):
                    continue
                if norm(r["sample_id"]) not in group_keys and norm(r["sample_name"]) not in group_keys:
                    warnings.append(f"{r['sample_id']}: on the sample sheet but has no FASTQ dataset in BaseSpace.")
                    report.append({"basespace_name": "", "project": r["project"], "dataset_ids": "",
                                   "lanes": "", "reads_pf": "", "sample_id": r["sample_id"],
                                   "status": "samplesheet_only",
                                   "note": "On the sample sheet but no FASTQ dataset in BaseSpace."})
        if not sheet_rows:
            notes.append("No sample sheet found for this run; skipped the sample sheet cross-check.")
        if batchid_fallbacks:
            notes.append(f"batchid for {batchid_fallbacks} sample(s) came from the BaseSpace project "
                         f"name because the sample sheet had no Sample_Project for them.")
    
        # ---- whole-run checks ---------------------------------------------------
        seen = {}
        for s in samples:
            if not ENTITY_NAME_RE.match(s["sample_id"]):
                failures.append(f"Sample ID '{s['sample_id']}' has characters Terra does not allow "
                                f"(use letters, numbers, '.', '_' and '-').")
            if s["sample_id"] in seen:
                failures.append(f"Two different datasets would both load as sample '{s['sample_id']}'.")
            seen[s["sample_id"]] = True
    
        loaded_ids = {norm(s["sample_id"]) for s in samples}
        for x in opts["extra_negative_control_ids"]:
            if norm(x) not in loaded_ids:
                warnings.append(f"extra_negative_control_ids lists '{x}', but no such sample was loaded.")
    
        to_fetch = [s for s in samples if s["load_status"] == "loaded"]
        if not samples and not failures:
            failures.append("No samples to load after applying the checks above.")
        elif samples and not to_fetch and not failures:
            failures.append("Every selected sample has 0 reads in BaseSpace; nothing to fetch.")
    
        controls = [s["sample_id"] for s in samples if s["negative_control"]]
        if samples and not controls:
            notes.append("No negative controls were identified in this load.")
    
        # ---- sample set name ------------------------------------------------------
        batchids = sorted({s["batchid"] for s in samples})
        if opts.get("sample_set_name"):
            set_name = opts["sample_set_name"]
        elif opts["pathogen"] == "SC2":
            set_name = opts["btb_batch_id"]
        elif len(batchids) == 1:
            set_name = batchids[0]
        else:
            set_name = f"run_{opts['run_id']}"
        if samples and not ENTITY_NAME_RE.match(set_name):
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", set_name)
            notes.append(f"Sample set name '{set_name}' changed to '{safe}' (Terra naming rules).")
            set_name = safe
    
        if failures:
            raise CheckFailed(failures)
    
        return {"samples": samples, "to_fetch": to_fetch, "report": report, "warnings": warnings,
                "notes": notes, "set_name": set_name, "controls": controls,
                "projects_in_run": projects_in_run}
    
    
    # --------------------------------------------------------------------------
    # Output writers
    # --------------------------------------------------------------------------
    
    ROW_COLUMNS = ["sample_sheet_id", "batchid", "pathogen", "basespace_run_id", "basespace_run_name",
                   "basespace_run_date", "basespace_project", "basespace_dataset_names",
                   "basespace_reads_pf", "negative_control", "load_status", "loader_version"]
    
    REPORT_COLUMNS = ["sample_id", "status", "basespace_name", "project", "dataset_ids", "lanes",
                      "reads_pf", "note"]
    
    
    def run_date(run):
        return (run.get("date_created") or "")[:10]
    
    
    def write_outputs(result, run, opts, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        table = opts["table_name"]
    
        with open(os.path.join(out_dir, "sample_rows.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([f"entity:{table}_id"] + ROW_COLUMNS)
            for s in result["samples"]:
                w.writerow([s["sample_id"], s["sample_sheet_id"], s["batchid"], opts["pathogen"],
                            opts["run_id"], run["experiment_name"] or run["name"], run_date(run),
                            s["project"], "|".join(s["dataset_names"]),
                            "" if s["reads_pf"] is None else s["reads_pf"],
                            "yes" if s["negative_control"] else "no", s["load_status"],
                            opts["loader_version"]])
    
        with open(os.path.join(out_dir, "membership.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([f"membership:{table}_set_id", table])
            for s in result["samples"]:
                w.writerow([result["set_name"], s["sample_id"]])
    
        # no header: read by WDL read_tsv() -> [sample_id, comma-separated dataset ids]
        with open(os.path.join(out_dir, "to_fetch.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            for s in result["to_fetch"]:
                w.writerow([s["sample_id"], ",".join(s["dataset_ids"])])
    
        with open(os.path.join(out_dir, "check_report.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            for r in result["report"]:
                w.writerow(r)
    
        with open(os.path.join(out_dir, "sample_set_name.txt"), "w") as fh:
            fh.write(result["set_name"])
    
        with open(os.path.join(out_dir, "check_summary.txt"), "w") as fh:
            fh.write(summary_text(result, run, opts))
    
    
    def summary_text(result, run, opts):
        s = result["samples"]
        lines = [
            "MASPHL viral NGS loader -- check_run_and_samples",
            f"check version     : {SCRIPT_VERSION}",
            f"pathogen          : {opts['pathogen']}",
            f"BaseSpace run     : {opts['run_id']} ({run['experiment_name'] or run['name']}), "
            f"status {run['status']}, created {run_date(run)}",
            f"target project    : {'(ignored: only_these_samples given)' if opts['only_these_samples'] else opts.get('project')}",
            f"projects in run   : {', '.join(result['projects_in_run'])}",
        ]
        if opts["pathogen"] == "SC2":
            lines.append(f"BtB batch ID      : {opts['btb_batch_id']}")
        lines += [
            f"Terra table / set : {opts['table_name']} / {result['set_name']}",
            "",
            f"Samples loaded    : {len(s)}  (to fetch: {len(result['to_fetch'])}, "
            f"0 reads: {sum(1 for x in s if x['load_status'] == 'no_reads')})",
            f"Negative controls : {', '.join(result['controls']) or 'none identified'}",
            "",
        ]
        if result["warnings"]:
            lines.append("WARNINGS")
            lines += [f"  - {w}" for w in result["warnings"]]
            lines.append("")
        if result["notes"]:
            lines.append("NOTES")
            lines += [f"  - {n}" for n in result["notes"]]
            lines.append("")
        lines.append("Per-dataset detail: check_report.tsv")
        return "\n".join(lines) + "\n"
    
    
    # --------------------------------------------------------------------------
    # CLI
    # --------------------------------------------------------------------------
    
    def build_opts(a):
        return {
            "pathogen": (a.pathogen or "").strip().upper(),
            "run_id": (a.run_id or "").strip(),
            "btb_batch_id": (a.btb_batch_id or "").strip(),
            "project": (a.project or "").strip(),
            "only_these_samples": parse_id_list(a.only_these_samples, "only_these_samples"),
            "extra_negative_control_ids": parse_id_list(a.extra_negative_control_ids,
                                                        "extra_negative_control_ids"),
            "control_regex": a.control_regex or DEFAULT_CONTROL_REGEX,
            "allow_incomplete_run": a.allow_incomplete_run.lower() == "true",
            "sample_set_name": (a.sample_set_name or "").strip(),
            "table_name": a.table_name,
            "loader_version": a.loader_version,
            "other_input_problems": read_problems(a.problems_file),
        }
    
    
    def read_problems(path):
        """Extra input problems found by the calling workflow (one per line)."""
        if not path or not os.path.exists(path):
            return []
        with open(path) as fh:
            return [line.strip() for line in fh if line.strip()]
    
    
    def write_empty_outputs(out_dir):
        """validate-only mode: create the files the WDL task expects, empty."""
        for name in ("sample_rows.tsv", "membership.tsv", "to_fetch.tsv", "check_report.tsv",
                     "sample_set_name.txt"):
            open(os.path.join(out_dir, name), "w").close()
    
    
    def main(argv=None):
        p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
        p.add_argument("--snapshot-dir", default="")
        p.add_argument("--out-dir", required=True)
        p.add_argument("--validate-only", action="store_true",
                       help="only check the inputs (no BaseSpace snapshot needed)")
        p.add_argument("--problems-file", default="",
                       help="file of extra input problems found by the workflow, one per line")
        p.add_argument("--pathogen", required=True)
        p.add_argument("--run-id", required=True)
        p.add_argument("--btb-batch-id", default="")
        p.add_argument("--project", default="")
        p.add_argument("--only-these-samples", default="")
        p.add_argument("--extra-negative-control-ids", default="")
        p.add_argument("--control-regex", default=DEFAULT_CONTROL_REGEX)
        p.add_argument("--allow-incomplete-run", default="false")
        p.add_argument("--sample-set-name", default="")
        p.add_argument("--table-name", default="sample")
        p.add_argument("--loader-version", default="dev")
        a = p.parse_args(argv)
    
        try:
            opts = build_opts(a)
            validate_inputs(opts)
            if a.validate_only:
                os.makedirs(a.out_dir, exist_ok=True)
                write_empty_outputs(a.out_dir)
                with open(os.path.join(a.out_dir, "check_summary.txt"), "w") as fh:
                    fh.write("Inputs OK.\n")
                sys.stdout.write("Inputs OK.\n")
                return 0
            run = load_run(os.path.join(a.snapshot_dir, "run.json"))
            datasets = load_datasets(os.path.join(a.snapshot_dir, "datasets.csv"))
            counts = load_read_counts(os.path.join(a.snapshot_dir, "attributes"),
                                      [d["id"] for d in datasets])
            sheets = find_samplesheet(os.path.join(a.snapshot_dir, "run_files"))
            sheet_rows = []
            if sheets:
                sheet_rows = parse_samplesheet(sheets[0])
            result = run_check(run, datasets, sheet_rows, counts, opts)
            if sheets:
                result["notes"].insert(0, f"Sample sheet used: {os.path.relpath(sheets[0], a.snapshot_dir)}")
            write_outputs(result, run, opts, a.out_dir)
        except CheckFailed as e:
            os.makedirs(a.out_dir, exist_ok=True)
            msg = "CHECK FAILED -- nothing was loaded.\n" + "\n".join(f"  - {r}" for r in e.reasons) + "\n"
            with open(os.path.join(a.out_dir, "check_summary.txt"), "w") as fh:
                fh.write(msg)
            sys.stderr.write(msg)
            return 1
    
        with open(os.path.join(a.out_dir, "check_summary.txt")) as fh:
            sys.stdout.write(fh.read())
        return 0
    
    
    if __name__ == "__main__":
        sys.exit(main())
    # END EMBEDDED scripts/check_run_and_samples.py
    PYEOF

    python3 check_run_and_samples.py \
      ~{true="--validate-only" false="" validate_only} \
      --problems-file '~{write_lines(other_input_problems)}' \
      --snapshot-dir snapshot \
      --out-dir out \
      --pathogen '~{pathogen}' \
      --run-id '~{basespace_run_id}' \
      --btb-batch-id '~{default="" btb_batch_id}' \
      --project '~{default="" basespace_project}' \
      --only-these-samples '~{default="" only_these_samples}' \
      --extra-negative-control-ids '~{default="" extra_negative_control_ids}' \
      --control-regex '~{control_id_regex}' \
      --allow-incomplete-run '~{allow_incomplete_run}' \
      --sample-set-name '~{default="" sample_set_name}' \
      --table-name '~{table_name}' \
      --loader-version '~{loader_version}'
  >>>
  output {
    File check_summary = "out/check_summary.txt"
    File check_report = "out/check_report.tsv"
    File sample_rows_tsv = "out/sample_rows.tsv"
    File membership_tsv = "out/membership.tsv"
    Array[Array[String]] fetch_table = read_tsv("out/to_fetch.tsv")
    String sample_set = read_string("out/sample_set_name.txt")
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 10 SSD"
    preemptible: 1
  }
}

task fetch_basespace_reads {
  input {
    String sample_id
    String dataset_ids            # comma-separated, already in lane order
    String api_server
    String access_token
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 4
  }
  command <<<
    set -euo pipefail

    set +x
    export HOME="$PWD/bs_home"
    mkdir -p "$HOME/.basespace"
    printf 'apiServer = %s\naccessToken = %s\n' '~{api_server}' '~{access_token}' > "$HOME/.basespace/default.cfg"
    chmod 600 "$HOME/.basespace/default.cfg"

    IFS=',' read -ra IDS <<< '~{dataset_ids}'
    n=0
    for ds_id in "${IDS[@]}"; do
      n=$((n + 1))
      dir=$(printf 'ds_%03d' "$n")
      mkdir -p "$dir"
      echo "Downloading dataset ${ds_id} (${n} of ${#IDS[@]})"
      if ! bs download dataset -i "$ds_id" -o "$dir" --extension=fastq.gz --retry > "${dir}.log" 2>&1; then
        cat "${dir}.log" >&2
        echo "ERROR: download of BaseSpace dataset ${ds_id} failed." >&2
        exit 1
      fi
    done

    # Collect R1/R2 files dataset by dataset (lane order), then file name order.
    : > r1_files.txt
    : > r2_files.txt
    for dir in ds_[0-9][0-9][0-9]; do
      find "$dir" -type f \( -name '*_R1_*.fastq.gz' -o -name '*_R1.fastq.gz' \) | sort >> r1_files.txt
      find "$dir" -type f \( -name '*_R2_*.fastq.gz' -o -name '*_R2.fastq.gz' \) | sort >> r2_files.txt
    done
    echo "R1 files:"; cat r1_files.txt
    echo "R2 files:"; cat r2_files.txt

    n1=$(wc -l < r1_files.txt)
    n2=$(wc -l < r2_files.txt)
    if [ "$n1" -lt 1 ]; then
      echo "ERROR: no R1 FASTQ files found for ~{sample_id} in datasets ~{dataset_ids}." >&2
      exit 1
    fi
    if [ "$n1" -ne "$n2" ]; then
      echo "ERROR: ~{sample_id} has ${n1} R1 files but ${n2} R2 files; expected paired-end data." >&2
      exit 1
    fi

    # gzip files can be joined end to end
    while read -r f; do cat "$f"; done < r1_files.txt > '~{sample_id}_R1.fastq.gz'
    while read -r f; do cat "$f"; done < r2_files.txt > '~{sample_id}_R2.fastq.gz'
    gzip -t '~{sample_id}_R1.fastq.gz'
    gzip -t '~{sample_id}_R2.fastq.gz'

    lines1=$(zcat '~{sample_id}_R1.fastq.gz' | wc -l)
    lines2=$(zcat '~{sample_id}_R2.fastq.gz' | wc -l)
    if [ "$lines1" -ne "$lines2" ]; then
      echo "ERROR: R1 and R2 have different read counts for ~{sample_id}." >&2
      exit 1
    fi
    echo $(( lines1 / 4 )) > READ_PAIRS
    echo "~{sample_id}: $(cat READ_PAIRS) read pairs from ${n1} file pair(s)"
  >>>
  output {
    File read1 = "~{sample_id}_R1.fastq.gz"
    File read2 = "~{sample_id}_R2.fastq.gz"
    Int read_pairs = read_int("READ_PAIRS")
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 1
    maxRetries: 1
  }
}
