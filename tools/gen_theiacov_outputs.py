#!/usr/bin/env python3
"""
Generate the WDL that copies TheiaCoV's outputs into the data table.

Why: the end-to-end workflow is launched without data-table rows, so Terra's
usual "Outputs -> this.<name>" mapping does not happen. Instead the workflow
gathers every TheiaCoV output into a TSV row and writes it with FISS. Listing
~195 outputs by hand is error-prone, so this tool writes those WDL lines from
config/theiacov_v<version>_outputs.tsv.

  python tools/gen_theiacov_outputs.py           # update the workflow
  python tools/gen_theiacov_outputs.py --check   # exit 1 if out of date

To move to a new TheiaCoV version, regenerate the outputs list from that tag:
  python tools/gen_theiacov_outputs.py --from-wdl <checkout>/workflows/theiacov/wf_theiacov_illumina_pe.wdl --version 4.x.y
(needs `pip install miniwdl`), then update the import URL in the workflow.

Rules
  * Every output becomes a column with the same name, as Terra's defaults do.
  * For MPXV clade I reruns, the rerun's value is used for every column.
  * Array outputs (flu-only irma_bams) are skipped.
  * File outputs are written as their gs:// path; nothing is downloaded.
"""
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
WORKFLOW = os.path.join(ROOT, "workflows", "wf_masphl_viral_ngs.wdl")
VERSION = "4.0.0"
OUTPUTS = os.path.join(ROOT, "config", f"theiacov_v{VERSION}_outputs.tsv")

PRIMARY_CALL = "theiacov_illumina_pe"
RERUN_CALL = "theiacov_illumina_pe_clade_i"


def load_outputs(path=OUTPUTS):
    outs = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            typ, name = line.rstrip("\n").split("\t")
            if typ.startswith("Array"):
                continue
            outs.append(name)
    return outs


def blocks(outs):
    cols = ["Array[String] theiacov_columns = ["]
    cols += [f'  "{n}"{"," if i < len(outs) - 1 else ""}' for i, n in enumerate(outs)]
    cols.append("]")
    vals = ["Array[String] theiacov_values = ["]
    for i, n in enumerate(outs):
        comma = "," if i < len(outs) - 1 else ""
        vals.append(f'  if (rerun_clade_i) then "~{{{RERUN_CALL}.{n}}}" else "~{{{PRIMARY_CALL}.{n}}}"{comma}')
    vals.append("]")
    return {"theiacov columns": cols, "theiacov values": vals}


def render(text, outs):
    gen = blocks(outs)
    for key, lines in gen.items():
        pat = re.compile(
            rf"(?P<indent>[ \t]*)# BEGIN GENERATED {key}[^\n]*\n.*?[ \t]*# END GENERATED {key}",
            re.S)
        m = pat.search(text)
        if not m:
            raise SystemExit(f"markers for '{key}' not found in workflow")
        ind = m.group("indent")
        body = "\n".join(ind + ln for ln in lines)
        new = (f"{ind}# BEGIN GENERATED {key} (tools/gen_theiacov_outputs.py; do not edit by hand)\n"
               f"{body}\n{ind}# END GENERATED {key}")
        text = text[:m.start()] + new + text[m.end():]
    return text


def from_wdl(path, version):
    import WDL  # miniwdl
    wf = WDL.load(path).workflow
    out = os.path.join(ROOT, "config", f"theiacov_v{version}_outputs.tsv")
    with open(out, "w") as fh:
        fh.write(f"# TheiaCoV_Illumina_PE_PHB v{version} workflow outputs (type<TAB>name)\n")
        for o in wf.outputs:
            fh.write(f"{o.type}\t{o.name}\n")
    print(f"wrote {out}")


def check():
    with open(WORKFLOW) as fh:
        text = fh.read()
    return render(text, load_outputs()) == text


def main():
    if "--from-wdl" in sys.argv:
        path = sys.argv[sys.argv.index("--from-wdl") + 1]
        version = sys.argv[sys.argv.index("--version") + 1]
        from_wdl(path, version)
        return 0
    if "--check" in sys.argv:
        ok = check()
        print("generated TheiaCoV columns up to date" if ok else "generated TheiaCoV columns are STALE")
        return 0 if ok else 1
    with open(WORKFLOW) as fh:
        text = fh.read()
    new = render(text, load_outputs())
    if new != text:
        with open(WORKFLOW, "w") as fh:
            fh.write(new)
        print("updated workflows/wf_masphl_viral_ngs.wdl")
    else:
        print("workflow already up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
