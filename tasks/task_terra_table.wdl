version 1.0

# Writing to the Terra data table from inside a workflow.
#
# Runs as the launching user's pet service account (the same access the user
# has), through the Terra API via FISS (the firecloud Python package).

task update_terra_table {
  input {
    String terra_project              # workspace namespace (billing project)
    String terra_workspace            # workspace name
    String table_name = "sample"
    File rows_tsv                     # entity:<table>_id + columns
    File? membership_tsv              # membership:<table>_set_id + <table>
    # Before writing, look for existing rows that came from a DIFFERENT BaseSpace
    # run (a repeat/resequenced sample). Those rows are overwritten, and listed.
    Boolean check_previous_run = false
    # For those replaced rows, also delete columns this write doesn't set (old
    # reads, old TheiaCoV results, old QC), so old and new results never mix.
    Boolean clear_old_columns_on_rerun = true
    File? wait_for                    # only used to order calls in the workflow
    String firecloud_version = "0.16.41"
    String docker = "python:3.12-slim"
  }
  command <<<
    set -euo pipefail
    if ! python3 -c "import firecloud" 2>/dev/null; then
      echo "Installing FISS (firecloud==~{firecloud_version})"
      pip install --quiet --no-cache-dir "firecloud==~{firecloud_version}"
    fi

    python3 <<'PYEOF'
    import csv, sys
    from firecloud import api as fapi

    ns, ws, table = "~{terra_project}", "~{terra_workspace}", "~{table_name}"
    rows_path, membership_path = "~{rows_tsv}", "~{default='' membership_tsv}"
    check_previous = "~{check_previous_run}" == "true"
    clear_old = "~{clear_old_columns_on_rerun}" == "true"

    def die(msg, resp=None):
        sys.stderr.write("ERROR: " + msg + "\n")
        if resp is not None:
            sys.stderr.write(f"HTTP {resp.status_code}: {resp.text[:2000]}\n")
        sys.exit(1)

    with open(rows_path, newline="") as fh:
        data = list(csv.reader(fh, delimiter="\t"))
    header, body = data[0], data[1:]
    if header[0] != f"entity:{table}_id":
        die(f"rows file is for '{header[0]}', expected 'entity:{table}_id'")
    columns = header[1:]
    log = [f"Workspace {ns}/{ws}, table '{table}': {len(body)} row(s), columns: {', '.join(columns)}"]

    replaced = []
    if check_previous:
        for row in body:
            sid, new = row[0], dict(zip(columns, row[1:]))
            resp = fapi.get_entity(ns, ws, table, sid)
            if resp.status_code == 404:
                continue
            if resp.status_code != 200:
                die(f"could not read existing row '{sid}'", resp)
            attrs = resp.json().get("attributes", {})
            old_run = str(attrs.get("basespace_run_id", "")).strip()
            if old_run == new.get("basespace_run_id", ""):
                continue                      # same run reprocessed: plain update
            stale = sorted(a for a in attrs if a not in columns)
            replaced.append([sid, old_run or "unknown (loaded before this workflow)",
                             str(attrs.get("batchid", "")), str(len(stale)) if clear_old else "0"])
            if clear_old and stale:
                ops = [{"op": "RemoveAttribute", "attributeName": a} for a in stale]
                resp = fapi.update_entity(ns, ws, table, sid, ops)
                if not resp.ok:
                    die(f"could not clear old columns on '{sid}'", resp)
        log.append(f"Rows replaced from a different BaseSpace run: {len(replaced)}")

    with open(rows_path) as fh:
        resp = fapi.upload_entities(ns, ws, fh.read(), model="flexible")
    if not resp.ok:
        die("upload of rows failed", resp)
    log.append("Rows written.")

    if membership_path:
        with open(membership_path) as fh:
            resp = fapi.upload_entities(ns, ws, fh.read(), model="flexible")
        if not resp.ok:
            die("upload of sample set membership failed", resp)
        log.append("Sample set written.")

    with open("replaced_rows.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "previous_basespace_run_id", "previous_batchid", "old_columns_removed"])
        w.writerows(replaced)
    with open("table_update_log.txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n".join(log))
    PYEOF
  >>>
  output {
    File update_log = "table_update_log.txt"
    File replaced_rows = "replaced_rows.tsv"
  }
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 10 SSD"
    preemptible: 1
    maxRetries: 1
  }
}

task make_read_rows {
  input {
    String table_name = "sample"
    Array[String] sample_ids
    Array[String] read1
    Array[String] read2
    Array[Int] read_pairs
    String timezone = "America/New_York"
    String docker = "python:3.12-slim"
  }
  command <<<
    set -euo pipefail
    today=$(TZ='~{timezone}' date +%Y-%m-%d)
    printf 'entity:%s_id\tread1\tread2\tbasespace_read_pairs_fetched\tbasespace_fetch_date\n' '~{table_name}' > read_rows.tsv
    paste '~{write_lines(sample_ids)}' '~{write_lines(read1)}' '~{write_lines(read2)}' '~{write_lines(read_pairs)}' \
      | awk -v d="$today" 'BEGIN{OFS="\t"} {print $0, d}' >> read_rows.tsv
    cat read_rows.tsv
  >>>
  output {
    File read_rows_tsv = "read_rows.tsv"
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 10 SSD"
    preemptible: 1
  }
}
