#!/usr/bin/env python3
"""Generate Markdown reference docs from INI files, using INI comments as the single source of truth.

Usage:
  python tools/gen_ini_docs.py \
    --ini call_atlas.ini --out docs/config/initialization.md --title "Initialization" \
    --ini config_file.ini --out docs/config/configuration.md --title "Configuration" \
    --ini settings_file.ini --out docs/config/settings.md --title "Settings"

Conventions:
- Full-line comments starting with '#' or ';' immediately above a section header become the section description.
- Full-line comments immediately above a key become that key's help text.
- Inline comments after a value (e.g. 'x = 1  # help') are appended to the help text.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import pathlib
import re

SECTION_RE = re.compile(r'^\s*\[(.+?)\]\s*$')
KEY_RE = re.compile(r'^\s*([^=:#;]+?)\s*=\s*(.*?)\s*$')


def parse_ini_with_comments(path: pathlib.Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    data = {"_sections": []}
    current = None
    pending_comments: list[str] = []
    in_preamble = True

    def flush_section():
        nonlocal current
        if current is not None:
            data["_sections"].append(current)
        current = None

    for raw in lines:
        line = raw.rstrip("\n")
        stripped = line.strip()

        if stripped == "":
            if pending_comments:
                pending_comments.append("")
            continue

        if stripped.startswith("#") or stripped.startswith(";"):
            pending_comments.append(stripped[1:].lstrip())
            continue

        msec = SECTION_RE.match(line)
        if msec:
            in_preamble = False
            flush_section()
            sec_name = msec.group(1).strip()
            sec_desc = "\n".join([c for c in pending_comments if c != ""]).strip()
            current = {"name": sec_name, "description": sec_desc, "keys": []}
            pending_comments = []
            continue

        mkey = KEY_RE.match(line)
        if mkey and current is not None:
            in_preamble = False
            key = mkey.group(1).strip()
            rest = mkey.group(2).strip()
            value = rest
            inline = ""

            # Split inline comment on whitespace + (# or ;) marker
            mm = re.search(r'(?<!\\)\s([#;])\s*', rest)
            if mm:
                idx = mm.start(1) - 1
                value = rest[:idx].rstrip()
                inline = rest[mm.end():].strip()

            desc = "\n".join([c for c in pending_comments if c != ""]).strip()
            if inline:
                desc = (desc + ("\n" if desc else "") + inline).strip()

            current["keys"].append({"key": key, "value": value, "description": desc})
            pending_comments = []
            continue

        # Anything else is kept as a note within the section (useful if your INI has special syntax)
        if not in_preamble and current is not None:
            current.setdefault("notes", []).append(line)

    flush_section()
    return data


def md_escape(s: str) -> str:
    return s.replace("|", "\\|")


def render_markdown(title: str, parsed: dict, source_path: pathlib.Path) -> str:
    now = _dt.datetime.now().strftime("%Y-%m-%d")
    out: list[str] = []
    out.append(f"# {title} INI reference (WIP)")
    out.append("")
    out.append(f"*Generated from* `{source_path.as_posix()}` on {now}.")
    out.append("")
    out.append("> This page is generated from the INI file comments so it stays in sync. Add/edit comments in the INI file and re-run the generator.")
    out.append("")
    out.append("## File purpose")
    out.append("")
    out.append("_TODO: Add a short narrative description here._")
    out.append("")
    out.append("## Sections and options")
    out.append("")

    for sec in parsed["_sections"]:
        out.append(f"### [{sec['name']}]")
        out.append("")
        if sec.get("description"):
            out.append(sec["description"])
            out.append("")
        if sec.get("keys"):
            out.append("| Key | Default / Example | Help |")
            out.append("|---|---|---|")
            for k in sec["keys"]:
                key = md_escape(k["key"])
                val = md_escape(k["value"])
                desc = md_escape(k["description"].replace("\n", "<br>")) if k.get("description") else ""
                out.append(f"| `{key}` | `{val}` | {desc} |")
            out.append("")
        if sec.get("notes"):
            out.append("**Notes (unparsed lines):**")
            out.append("")
            out.append("```ini")
            out.extend(sec["notes"])
            out.append("```")
            out.append("")

    out.append("## Full file (verbatim)")
    out.append("")
    out.append("```ini")
    out.append(source_path.read_text(encoding='utf-8', errors='replace').rstrip())
    out.append("```")
    out.append("")
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ini", action="append", required=True, help="Path to an INI file (repeatable).")
    ap.add_argument("--out", action="append", required=True, help="Output markdown path (repeatable, same count as --ini).")
    ap.add_argument("--title", action="append", required=True, help="Title used in the markdown H1 (repeatable).")
    args = ap.parse_args()

    if not (len(args.ini) == len(args.out) == len(args.title)):
        raise SystemExit("ERROR: --ini, --out, and --title must be provided the same number of times.")

    for ini_path, out_path, title in zip(args.ini, args.out, args.title):
        ini_p = pathlib.Path(ini_path)
        out_p = pathlib.Path(out_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        parsed = parse_ini_with_comments(ini_p)
        md = render_markdown(title, parsed, ini_p)
        out_p.write_text(md, encoding="utf-8")
        print(f"Wrote {out_p} from {ini_p}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
