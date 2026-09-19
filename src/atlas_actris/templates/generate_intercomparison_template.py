#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate ATLAS intercomparison INI templates and schema reference docs."""

from __future__ import annotations

import argparse
import html
import textwrap
from pathlib import Path
from typing import Any, Mapping

from atlas_actris.utils.parse_intercomparison_file import (
    GENERAL_SCHEMA,
    PLOTTING_SCHEMA,
    DATASET_SCHEMA,
    CHANNEL_GROUP_SCHEMA,
    PAIR_GROUP_SCHEMA,
)
from atlas_actris.templates.intercomparison_template_flavor import (
    GENERAL_FLAVOR,
    PLOTTING_FLAVOR,
    DATASET_FLAVOR,
    CHANNEL_GROUP_FLAVOR,
    PAIR_GROUP_FLAVOR,
    SPECIAL_FLAVOR,
)
from atlas_actris.templates.template_profiles import (
    GENERATED_TEMPLATES_DIRECTORY,
    INTERCOMPARISON_TEMPLATE_KEYS,
)

GENERAL_TEMPLATE_KEYS = INTERCOMPARISON_TEMPLATE_KEYS["general"]
PLOTTING_TEMPLATE_KEYS = INTERCOMPARISON_TEMPLATE_KEYS["plotting"]
DATASET_TEMPLATE_KEYS = INTERCOMPARISON_TEMPLATE_KEYS["dataset"]
CHANNEL_GROUP_TEMPLATE_KEYS = INTERCOMPARISON_TEMPLATE_KEYS["channel_group"]
PAIR_GROUP_TEMPLATE_KEYS = INTERCOMPARISON_TEMPLATE_KEYS["pair_group"]


def validate_template_layout() -> None:
    """Fail early when the explicit template whitelist drifts from a schema."""

    blocks = {
        "general": (GENERAL_SCHEMA, GENERAL_FLAVOR),
        "plotting": (PLOTTING_SCHEMA, PLOTTING_FLAVOR),
        "dataset": (DATASET_SCHEMA, DATASET_FLAVOR),
        "channel_group": (CHANNEL_GROUP_SCHEMA, CHANNEL_GROUP_FLAVOR),
        "pair_group": (PAIR_GROUP_SCHEMA, PAIR_GROUP_FLAVOR),
    }
    errors: list[str] = []
    for block, keys in INTERCOMPARISON_TEMPLATE_KEYS.items():
        if block not in blocks:
            errors.append(f"intercomparison: unknown template block {block!r}")
            continue
        schema, flavor = blocks[block]
        duplicates = sorted({key for key in keys if keys.count(key) > 1})
        unknown = sorted(set(keys) - set(schema))
        missing_flavor = sorted(set(keys) - set(flavor))
        if not keys:
            errors.append(f"intercomparison/{block}: template contains no keys")
        if duplicates:
            errors.append(f"intercomparison/{block}: duplicate keys: {duplicates}")
        if unknown:
            errors.append(f"intercomparison/{block}: keys not present in schema: {unknown}")
        if missing_flavor:
            errors.append(f"intercomparison/{block}: selected keys missing flavor: {missing_flavor}")
    if errors:
        raise ValueError("\n".join(errors))


def _format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item) for item in value)
    return str(value)


def _default_display(meta: Mapping[str, Any]) -> str:
    value = meta.get("default")
    if value is None:
        if meta.get("category") == "mandatory":
            return "none (must be provided)"
        return "empty (None)"
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item) for item in value) if value else "empty list"
    if isinstance(value, Path):
        display = str(value)
        if not value.is_absolute() and not display.startswith("."):
            display = f"./{display}"
        return display
    return str(value)


def _allowed_display(meta: Mapping[str, Any]) -> str:
    parts: list[str] = []
    if meta.get("allowed") is not None:
        parts.append(", ".join(str(item) for item in meta["allowed"]))
    if "min" in meta:
        parts.append(f"> {meta['min']}")
    if "max" in meta:
        parts.append(f"<= {meta['max']}")
    if meta.get("size") is not None:
        parts.append(f"size={meta['size']}")
    return "; ".join(parts)


def _wrap(text: str, width: int = 88) -> list[str]:
    return [f"# {line}" if line else "#" for line in textwrap.wrap(text, width=width)]


def _render_entry(
    key: str,
    meta: Mapping[str, Any],
    flavor: Mapping[str, Any],
    *,
    include_flavor: bool,
) -> str:
    lines: list[str] = []
    if include_flavor:
        lines.extend(_wrap(str(flavor["description"])))
        if meta.get("category") == "mandatory":
            lines.append("# Required: yes")
        else:
            lines.append("# Required: no")
        lines.append(f"# Default: {_default_display(meta)}")
        allowed = _allowed_display(meta)
        if allowed:
            lines.append(f"# Allowed / limits: {allowed}")
        if flavor.get("example"):
            lines.append(f"# Example: {flavor['example']}")
    lines.append(f"{key} =")
    return "\n".join(lines)


def _header() -> str:
    return "\n".join([
        "# " + "=" * 76,
        "# ATLAS intercomparison initialization file",
        "#",
        "# Generated from parser schemas, flavor metadata, and template_profiles.py.",
        "# Values are intentionally left empty. Parser defaults are applied at runtime.",
        "# Regenerate with: atlas-generate-templates",
        "#",
        "# Dataset sections define where/how data are read.",
        "# Channel/pair group sections define which corresponding products are compared.",
        "#",
        "# " + "=" * 76,
        "",
    ])


def _render_schema_block(
    lines: list[str],
    schema: Mapping[str, Mapping[str, Any]],
    flavor: Mapping[str, Mapping[str, Any]],
    keys: list[str],
    *,
    include_flavor: bool,
) -> None:
    for key in keys:
        lines.extend([
            _render_entry(key, schema[key], flavor[key], include_flavor=include_flavor),
            "",
        ])


def render_ini(*, include_flavor: bool = True) -> str:
    validate_template_layout()
    lines = [_header(), "[general]", ""]
    _render_schema_block(
        lines, GENERAL_SCHEMA, GENERAL_FLAVOR, GENERAL_TEMPLATE_KEYS,
        include_flavor=include_flavor,
    )

    lines.extend([
        "",
        "# Global plotting defaults. Group-level plotting values override these.",
        "[plotting]", "",
    ])
    _render_schema_block(
        lines, PLOTTING_SCHEMA, PLOTTING_FLAVOR, PLOTTING_TEMPLATE_KEYS,
        include_flavor=include_flavor,
    )

    lines.extend([
        "",
        "# Repeat one [dataset:<dataset_id>] section per comparison dataset.",
        "# Dataset IDs are arbitrary; reference status is defined only by reference = True.",
        "[dataset:dataset_a_reference]", "",
    ])
    _render_schema_block(
        lines, DATASET_SCHEMA, DATASET_FLAVOR, DATASET_TEMPLATE_KEYS,
        include_flavor=include_flavor,
    )

    lines.extend(["", "[dataset:dataset_b]", ""])
    _render_schema_block(
        lines, DATASET_SCHEMA, DATASET_FLAVOR, DATASET_TEMPLATE_KEYS,
        include_flavor=include_flavor,
    )

    lines.extend([
        "",
        "# Repeat one [channel_group:<group_id>] section per channel comparison group.",
        "[channel_group:355_parallel]", "",
    ])
    _render_schema_block(
        lines, CHANNEL_GROUP_SCHEMA, CHANNEL_GROUP_FLAVOR, CHANNEL_GROUP_TEMPLATE_KEYS,
        include_flavor=include_flavor,
    )

    for dataset_id in ("dataset_a_reference", "dataset_b"):
        if include_flavor:
            lines.extend([
                f"# --- {dataset_id} ---",
                *(_wrap(SPECIAL_FLAVOR["dataset_atlas_channel_id"]["description"])),
                f"# Example: {SPECIAL_FLAVOR['dataset_atlas_channel_id']['example']}",
                f"{dataset_id}.atlas_channel_id =", "",
            ])
        else:
            lines.extend([f"{dataset_id}.atlas_channel_id =", ""])

    lines.extend([
        "",
        "# Repeat one [pair_group:<group_id>] section per pair comparison group.",
        "[pair_group:vldr_355]", "",
    ])
    _render_schema_block(
        lines, PAIR_GROUP_SCHEMA, PAIR_GROUP_FLAVOR, PAIR_GROUP_TEMPLATE_KEYS,
        include_flavor=include_flavor,
    )

    for dataset_id in ("dataset_a_reference", "dataset_b"):
        if include_flavor:
            lines.extend([
                f"# --- {dataset_id} ---",
                *(_wrap(SPECIAL_FLAVOR["dataset_atlas_pair_id"]["description"])),
                f"# Example: {SPECIAL_FLAVOR['dataset_atlas_pair_id']['example']}",
                f"{dataset_id}.atlas_pair_id =", "",
            ])
        else:
            lines.extend([f"{dataset_id}.atlas_pair_id =", ""])

    return "\n".join(lines).rstrip() + "\n"


def render_bare_ini() -> str:
    return render_ini(include_flavor=False)


def _dtype(meta: Mapping[str, Any]) -> str:
    name = getattr(meta["dtype"], "__name__", str(meta["dtype"]))
    return f"list[{name}]" if meta.get("is_list") else name


def _md_table(
    schema: Mapping[str, Mapping[str, Any]],
    flavor: Mapping[str, Mapping[str, Any]],
    keys: list[str],
) -> str:
    lines = [
        "| Parameter | Description | Type | Default | Allowed / limits | Example |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for key in keys:
        meta = schema[key]
        entry = flavor[key]
        values = [
            f"`{key}`",
            html.escape(str(entry["description"])).replace("|", "\\|"),
            f"`{_dtype(meta)}`",
            f"`{html.escape(_format_value(meta.get('default')))}`",
            html.escape(_allowed_display(meta)).replace("|", "\\|"),
            f"`{html.escape(str(entry.get('example', '')))}`" if entry.get("example") else "",
        ]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def render_markdown() -> str:
    validate_template_layout()
    return "\n".join([
        "# ATLAS intercomparison initialization reference",
        "",
        '!!! note "Generated reference"',
        "    This page is generated from the intercomparison parser schemas and flavor text.",
        "",
        "## Structure",
        "",
        "The file contains one `[general]` section, one optional `[plotting]` section, repeated `[dataset:<id>]` sections, and repeated `[channel_group:<id>]` and `[pair_group:<id>]` sections.",
        "",
        "Dataset sections define the exported stage, QA test, source parameters, and labels for each independently selectable comparison dataset. Different datasets may use different stages, QA tests, or source parameters, even when they originate from the same physical lidar.",
        "",
        "Exactly one dataset must set `reference = True`. Group sections map participating datasets to one `atlas_channel_id` or `atlas_pair_id`. Empty or omitted mappings exclude that dataset. A group with no provided IDs is ignored. For an active group, the reference dataset must provide its ID.",
        "",
        "## `general`", "", _md_table(GENERAL_SCHEMA, GENERAL_FLAVOR, GENERAL_TEMPLATE_KEYS), "",
        "## `plotting`", "", _md_table(PLOTTING_SCHEMA, PLOTTING_FLAVOR, PLOTTING_TEMPLATE_KEYS), "",
        "## `dataset:<dataset_id>`", "", _md_table(DATASET_SCHEMA, DATASET_FLAVOR, DATASET_TEMPLATE_KEYS), "",
        "## `channel_group:<group_id>`", "", _md_table(CHANNEL_GROUP_SCHEMA, CHANNEL_GROUP_FLAVOR, CHANNEL_GROUP_TEMPLATE_KEYS), "",
        "### Dynamic channel-group entries", "",
        f"- `<dataset_id>.atlas_channel_id = <atlas_channel_id>` selects one channel from that dataset. Empty or omitted mappings exclude the dataset. If no dataset provides an ID, the group is ignored. For an active group, the reference dataset must provide an ID. Example: `{SPECIAL_FLAVOR['dataset_atlas_channel_id']['example']}`.",
        "",
        "## `pair_group:<group_id>`", "", _md_table(PAIR_GROUP_SCHEMA, PAIR_GROUP_FLAVOR, PAIR_GROUP_TEMPLATE_KEYS), "",
        "### Dynamic pair-group entries", "",
        f"- `<dataset_id>.atlas_pair_id = <atlas_pair_id>` selects one pair from that dataset. Empty or omitted mappings exclude the dataset. If no dataset provides an ID, the group is ignored. For an active group, the reference dataset must provide an ID. Example: `{SPECIAL_FLAVOR['dataset_atlas_pair_id']['example']}`.",
        "",
        "## Deferred metadata defaults", "",
        "Dataset labels may remain empty after parsing. Background, normalization, and molecular-plot controls are resolved per group from [general]. Plotting controls are resolved per group from [plotting], except use_log_y_scale, whose group-level defaults are True for channel groups and False for pair groups.",
        "",
    ])


def _repo_root_from_this_file() -> Path:
    return Path(__file__).resolve().parents[3]


def build_outputs(repo_root: str | Path | None = None) -> dict[Path, str]:
    root = Path(repo_root).resolve() if repo_root else _repo_root_from_this_file()
    return {
        root / GENERATED_TEMPLATES_DIRECTORY / "intercomparison.ini": render_ini(),
        root / GENERATED_TEMPLATES_DIRECTORY / "intercomparison_bare.ini": render_bare_ini(),
        root / "docs" / "generated" / "intercomparison_reference.md": render_markdown(),
    }


def generate(repo_root: str | Path | None = None, *, check: bool = False) -> dict[Path, str]:
    outputs = build_outputs(repo_root)
    changed: list[Path] = []
    for path, content in outputs.items():
        if check:
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                changed.append(path)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
    if changed:
        raise SystemExit(
            "Generated intercomparison files are out of date:\n"
            + "\n".join(str(path) for path in changed)
        )
    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=None)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    outputs = generate(args.repo_root, check=args.check)
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
