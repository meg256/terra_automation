#!/usr/bin/env python3
"""
Copy Python scripts from scripts/ into the WDL task files that run them.

Why: the tasks use plain public images (python:3.12-slim), so the script has to
travel inside the WDL. scripts/ stays the single source of truth and is what
the tests exercise; this tool keeps the WDL copy identical.

  python tools/embed_scripts.py          # update the WDL files
  python tools/embed_scripts.py --check  # exit 1 if any embedded copy is stale

Each embed point in a WDL looks like:
    # BEGIN EMBEDDED scripts/<name>.py ...
    # END EMBEDDED scripts/<name>.py
Everything between the two marker lines is replaced, indented to match.
"""
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
WDL_FILES = ["tasks/task_basespace.wdl", "tasks/task_qc_reports.wdl"]
BEGIN_RE = re.compile(r"^(?P<indent>[ \t]*)# BEGIN EMBEDDED (?P<path>scripts/\S+\.py)\b.*$")
FORBIDDEN = ["~{", ">>>"]   # would break a WDL heredoc command block


def render(wdl_text):
    out, lines, i = [], wdl_text.split("\n"), 0
    while i < len(lines):
        m = BEGIN_RE.match(lines[i])
        out.append(lines[i])
        i += 1
        if not m:
            continue
        indent, rel = m.group("indent"), m.group("path")
        end_marker = f"# END EMBEDDED {rel}"
        while i < len(lines) and lines[i].strip() != end_marker:
            i += 1
        if i == len(lines):
            raise SystemExit(f"missing '{end_marker}'")
        with open(os.path.join(ROOT, rel)) as fh:
            script = fh.read().rstrip("\n")
        for bad in FORBIDDEN:
            if bad in script:
                raise SystemExit(f"{rel} contains '{bad}', which breaks the WDL command block")
        if "PYEOF" in script:
            raise SystemExit(f"{rel} contains the heredoc terminator PYEOF")
        # indent every line (blank ones too) so WDL's common-indent stripping is uniform
        out.extend(indent + line for line in script.split("\n"))
        out.append(lines[i])          # END marker
        i += 1
    return "\n".join(out)


def check():
    for rel in WDL_FILES:
        path = os.path.join(ROOT, rel)
        with open(path) as fh:
            text = fh.read()
        if render(text) != text:
            return False
    return True


def main():
    if "--check" in sys.argv:
        ok = check()
        print("embedded scripts up to date" if ok else "embedded scripts are STALE")
        return 0 if ok else 1
    for rel in WDL_FILES:
        path = os.path.join(ROOT, rel)
        with open(path) as fh:
            text = fh.read()
        new = render(text)
        if new != text:
            with open(path, "w") as fh:
                fh.write(new)
            print(f"updated {rel}")
        else:
            print(f"{rel} already up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
