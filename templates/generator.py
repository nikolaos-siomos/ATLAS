#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate ATLAS INI templates and MkDocs reference pages.

Parser schemas own technical metadata, flavor files own user-facing text, and
``template_profiles.py`` explicitly owns inclusion and ordering. This generator
can write three INI template profiles:

* full:     all schema parameters with flavor comments
* bare:     all schema parameters without per-parameter comments
* beginner: selected parameters with flavor comments, configured in
            template_profiles.py
"""

from __future__ import annotations

import argparse
import html
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

# ATLAS currently has several intra-package imports written as top-level imports
# such as ``from utils...``.  When this generator is executed as an installed
# console script, make src/atlas_actris available as a top-level import root too.
_PACKAGE_DIR = Path(__file__).resolve().parents[1]
if str(_PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(_PACKAGE_DIR))

try:
    from atlas_actris.version import __version__ as ATLAS_VERSION
except Exception:  # pragma: no cover - fallback for direct source-tree execution
    ATLAS_VERSION = "unknown"

from atlas_actris.utils.parse_init_file import SCHEMA as INIT_SCHEMA
from atlas_actris.utils.parse_config_file import SCHEMA as CONFIG_SCHEMA
from atlas_actris.utils.parse_settings_file import (
    SCHEMA as SETTINGS_SCHEMA,
    recognized_sections as SETTINGS_RECOGNIZED_SECTIONS,
)
from atlas_actris.templates.init_template_flavor import INIT_FLAVOR
from atlas_actris.templates.config_template_flavor import CONFIG_FLAVOR
from atlas_actris.templates.settings_template_flavor import SETTINGS_FLAVOR
from atlas_actris.templates.template_profiles import (
    GENERATED_TEMPLATES_DIRECTORY,
    INIT_TEMPLATE_SECTIONS,
    CONFIG_TEMPLATE_SECTIONS,
    SETTINGS_TEMPLATE_SECTIONS,
    SETTINGS_TEMPLATE_KEYS,
    BEGINNER_INIT_TEMPLATE_SECTIONS,
    BEGINNER_CONFIG_TEMPLATE_SECTIONS,
    BEGINNER_SETTINGS_TEMPLATE_KEYS,
)
from atlas_actris.templates.generate_intercomparison_template import (
    render_ini as render_intercomparison_ini,
    render_bare_ini as render_intercomparison_bare_ini,
    render_markdown as render_intercomparison_markdown,
)

AUTO_SELECTION = "empty parameter --> ignored or automatic selection"
ALLOWED_TARGETS = ("all", "ini", "docs")
ALLOWED_PROFILES = ("all", "full", "bare", "beginner")
ALLOWED_LEGACY_STATUS = {"", "new", "unchanged", "renamed", "moved"}
REQUIRED_FLAVOR_FIELDS = ("description", "example", "legacy")
REQUIRED_LEGACY_FIELDS = (
    "status",
    "introduced",
    "old_names",
    "old_location",
    "old_names_removed_in",
    "note",
)


# ---------------------------------------------------------------------------
# Generic formatting helpers
# ---------------------------------------------------------------------------


def _repo_root_from_this_file() -> Path:
    # src/atlas_actris/templates/generator.py -> repository root
    return Path(__file__).resolve().parents[3]


def _type_name(dtype: Any) -> str:
    if dtype is None:
        return "unknown"
    return getattr(dtype, "__name__", str(dtype).replace("<class '", "").replace("'>", ""))


def _type_display(meta: Mapping[str, Any]) -> str:
    dtype = _type_name(meta.get("dtype"))
    if meta.get("is_list"):
        return f"list[{dtype}]"
    return dtype


def _format_value(value: Any) -> str:
    if value is None or value == "" or value == []:
        return AUTO_SELECTION
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return AUTO_SELECTION
        return ", ".join(_format_value(v) for v in value)
    if isinstance(value, bool):
        return "True" if value else "False"
    return str(value)


def _category_display(meta: Mapping[str, Any]) -> str:
    return str(meta.get("category", "optional"))


def _allowed_range_display(meta: Mapping[str, Any]) -> str:
    parts: list[str] = []

    if "allowed" in meta:
        parts.append("Allowed values: " + ", ".join(str(v) for v in meta["allowed"]))

    has_min = "min" in meta
    has_max = "max" in meta
    if has_min and has_max:
        parts.append(f"Range: {meta['min']} to {meta['max']}")
    elif has_min:
        parts.append(f"Minimum: {meta['min']}")
    elif has_max:
        parts.append(f"Maximum: {meta['max']}")

    if "size" in meta:
        size = meta["size"]
        if isinstance(size, list):
            parts.append("Allowed list sizes: " + ", ".join(str(v) for v in size))
        else:
            parts.append(f"Required list size: {size}")

    return "; ".join(parts)


def _wrap_comment_line(text: str, width: int = 88) -> list[str]:
    import textwrap

    if not text:
        return []
    wrapped = textwrap.wrap(str(text), width=width, replace_whitespace=False)
    return [f"# {line}" for line in wrapped]


def _join_names(names: Iterable[Any]) -> str:
    return ", ".join(str(name) for name in names if str(name).strip())


def _legacy_sentence(
    legacy: Mapping[str, Any],
    *,
    for_docs: bool = False,
) -> str:
    """Return one version-history sentence for a flavor legacy dictionary.

    Empty status is treated as ``unchanged`` so templates remain readable while
    flavor files are being filled gradually.
    """

    status = str(legacy.get("status", "") or "").strip()
    if status == "":
        status = "unchanged"
    introduced = str(legacy.get("introduced", "") or "").strip()
    old_location = str(legacy.get("old_location", "") or "").strip()
    removed = str(legacy.get("old_names_removed_in", "") or "").strip()
    note = str(legacy.get("note", "") or "").strip()
    old_names = legacy.get("old_names", []) or []

    def code(value: str) -> str:
        return f"`{value}`" if for_docs else value

    def atlas_version_text(version: str) -> str:
        return f"ATLAS {code(version)}" if version else "the current ATLAS version"

    old_names_text = _join_names(code(str(name)) for name in old_names)

    if status == "new":
        sentence = f"new in {atlas_version_text(introduced)}."

    elif status == "unchanged":
        sentence = "unchanged from previous ATLAS versions."

    elif status == "renamed":
        sentence = f"renamed in {atlas_version_text(introduced)}"
        if old_names_text:
            sentence += f" from {old_names_text}"
        sentence += "."
        if removed:
            sentence += f" Old name{'s' if len(old_names) != 1 else ''} removed in {atlas_version_text(removed)}."

    elif status == "moved":
        sentence = f"moved in {atlas_version_text(introduced)}"
        if old_location:
            sentence += f" from {code(old_location)}"
        sentence += "."
        if old_names_text:
            sentence += f" Previous name{'s' if len(old_names) != 1 else ''}: {old_names_text}."
        if removed:
            sentence += f" Old name{'s' if len(old_names) != 1 else ''} removed in {atlas_version_text(removed)}."

    else:
        sentence = "unchanged from previous ATLAS versions."

    if note:
        sentence += f" {note}"

    return sentence[:1].upper() + sentence[1:]


def _version_note_line(name: str, legacy: Mapping[str, Any]) -> str:
    """Return the one-line INI version note for a parameter."""

    return f"Version note for {name}: {_legacy_sentence(legacy)}"


def _template_header(title: str, flavor_file: str | None = None) -> str:
    lines = [
        "# " + "=" * 76,
        f"# {title}",
        "# " + "=" * 76,
        "#",
        f"# Generated by ATLAS {ATLAS_VERSION}",
        f"# Template format version: {ATLAS_VERSION}",
        "#",
        "# This file was generated automatically from the ATLAS parser schema",
        "# using the explicit selection in template_profiles.py.",
    ]

    if flavor_file:
        lines.extend(
            [
                "#",
                "# Developers should not edit this generated file manually.",
                f"# To update user-facing comments and examples, edit: {flavor_file}",
                "# To change included parameters or ordering, edit: src/atlas_actris/templates/template_profiles.py",
                "# Then regenerate with: atlas-generate-templates",
            ]
        )
    else:
        lines.extend(
            [
                "#",
                "# Developers should not edit this generated file manually.",
                "# Regenerate with: atlas-generate-templates",
            ]
        )

    lines.extend(
        [
            "#",
            "# Users should copy this file and fill in the empty values.",
            "# Empty values let ATLAS use its default or automatic selection when available.",
            "#",
            "# " + "=" * 76,
            "",
        ]
    )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _validate_flavor_entry(label: str, key: str, entry: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []

    for field in REQUIRED_FLAVOR_FIELDS:
        if field not in entry:
            errors.append(f"{label}: {key}.{field} is missing")

    legacy = entry.get("legacy", {})
    if not isinstance(legacy, Mapping):
        errors.append(f"{label}: {key}.legacy must be a dictionary")
        return errors

    for field in REQUIRED_LEGACY_FIELDS:
        if field not in legacy:
            errors.append(f"{label}: {key}.legacy.{field} is missing")

    status = legacy.get("status", "")
    if status not in ALLOWED_LEGACY_STATUS:
        errors.append(
            f"{label}: {key}.legacy.status must be one of {sorted(ALLOWED_LEGACY_STATUS)}, got {status!r}"
        )

    if not isinstance(legacy.get("old_names", []), list):
        errors.append(f"{label}: {key}.legacy.old_names must be a list")

    return errors


def _section_key_list(sections: Mapping[str, Iterable[str]]) -> list[str]:
    keys: list[str] = []
    for section_keys in sections.values():
        keys.extend(list(section_keys))
    return keys


def _validate_flat_template(
    *,
    schema: Mapping[str, Mapping[str, Any]],
    sections: Mapping[str, Iterable[str]],
    flavor: Mapping[str, Mapping[str, Any]],
    label: str,
) -> None:
    schema_keys = set(schema.keys())
    flavor_keys = set(flavor.keys())
    section_keys = _section_key_list(sections)
    section_key_set = set(section_keys)
    duplicates = sorted({key for key in section_keys if section_keys.count(key) > 1})

    errors: list[str] = []
    missing_flavor = schema_keys - flavor_keys
    unknown_flavor = flavor_keys - schema_keys
    unknown_sections = section_key_set - schema_keys

    if missing_flavor:
        errors.append(f"{label}: schema keys missing flavor entries: {sorted(missing_flavor)}")
    if unknown_flavor:
        errors.append(f"{label}: flavor keys not present in schema: {sorted(unknown_flavor)}")
    if unknown_sections:
        errors.append(f"{label}: section keys not present in schema: {sorted(unknown_sections)}")
    if duplicates:
        errors.append(f"{label}: keys assigned to more than one section: {duplicates}")

    for key, entry in flavor.items():
        errors.extend(_validate_flavor_entry(label, key, entry))

    if errors:
        raise ValueError("\n".join(errors))


def _validate_subset_flat(
    *,
    schema: Mapping[str, Mapping[str, Any]],
    sections: Mapping[str, Iterable[str]],
    label: str,
) -> None:
    schema_keys = set(schema.keys())
    section_keys = _section_key_list(sections)
    section_key_set = set(section_keys)
    duplicates = sorted({key for key in section_keys if section_keys.count(key) > 1})
    unknown_sections = section_key_set - schema_keys

    errors: list[str] = []
    if unknown_sections:
        errors.append(f"{label}: section keys not present in schema: {sorted(unknown_sections)}")
    if duplicates:
        errors.append(f"{label}: keys assigned to more than one section: {duplicates}")
    if not section_key_set:
        errors.append(f"{label}: template contains no keys")

    if errors:
        raise ValueError("\n".join(errors))


def _validate_settings() -> None:
    schema_groups = set(SETTINGS_SCHEMA.keys())
    section_groups = set(SETTINGS_TEMPLATE_SECTIONS.keys())
    key_groups = set(SETTINGS_TEMPLATE_KEYS.keys())
    flavor_groups = set(SETTINGS_FLAVOR.keys())

    errors: list[str] = []
    unknown_groups = section_groups - schema_groups
    if unknown_groups:
        errors.append(f"settings: template groups not present in schema: {sorted(unknown_groups)}")
    if section_groups != key_groups:
        errors.append(
            "settings: section and parameter-selection groups do not match: "
            f"sections={sorted(section_groups)}, keys={sorted(key_groups)}"
        )
    if schema_groups != flavor_groups:
        errors.append(
            "settings: flavor groups do not match schema groups: "
            f"schema={sorted(schema_groups)}, flavor={sorted(flavor_groups)}"
        )

    for group in SETTINGS_TEMPLATE_SECTIONS:
        section_name = SETTINGS_TEMPLATE_SECTIONS[group]
        expected_section_name = SETTINGS_RECOGNIZED_SECTIONS[group]
        if section_name != expected_section_name:
            errors.append(
                f"settings: group {group} uses INI section {section_name!r}, "
                f"but parser expects {expected_section_name!r}"
            )

        schema = SETTINGS_SCHEMA[group]
        flavor = SETTINGS_FLAVOR[group]
        schema_keys = set(schema.keys())
        flavor_keys = set(flavor.keys())
        template_keys = list(SETTINGS_TEMPLATE_KEYS.get(group, ()))

        missing_flavor = schema_keys - flavor_keys
        unknown_flavor = flavor_keys - schema_keys
        if missing_flavor:
            errors.append(f"settings/{group}: schema keys missing flavor entries: {sorted(missing_flavor)}")
        if unknown_flavor:
            errors.append(f"settings/{group}: flavor keys not present in schema: {sorted(unknown_flavor)}")
        unknown_template = set(template_keys) - schema_keys
        duplicates = sorted({key for key in template_keys if template_keys.count(key) > 1})
        if unknown_template:
            errors.append(f"settings/{group}: template keys not present in schema: {sorted(unknown_template)}")
        if duplicates:
            errors.append(f"settings/{group}: duplicate template keys: {duplicates}")
        if not template_keys:
            errors.append(f"settings/{group}: template contains no keys")

        for key, entry in flavor.items():
            errors.extend(_validate_flavor_entry(f"settings/{group}", key, entry))

    if errors:
        raise ValueError("\n".join(errors))


def _validate_beginner_settings() -> None:
    errors: list[str] = []
    for group, keys in BEGINNER_SETTINGS_TEMPLATE_KEYS.items():
        if group not in SETTINGS_SCHEMA:
            errors.append(f"beginner settings: group not present in schema: {group}")
            continue
        schema_keys = set(SETTINGS_SCHEMA[group].keys())
        key_list = list(keys)
        unknown = set(key_list) - schema_keys
        duplicates = sorted({key for key in key_list if key_list.count(key) > 1})
        if unknown:
            errors.append(f"beginner settings/{group}: keys not present in schema: {sorted(unknown)}")
        if duplicates:
            errors.append(f"beginner settings/{group}: duplicate keys: {duplicates}")

    if errors:
        raise ValueError("\n".join(errors))


def validate_all() -> None:
    _validate_flat_template(
        schema=INIT_SCHEMA,
        sections=INIT_TEMPLATE_SECTIONS,
        flavor=INIT_FLAVOR,
        label="initialization",
    )
    _validate_flat_template(
        schema=CONFIG_SCHEMA,
        sections=CONFIG_TEMPLATE_SECTIONS,
        flavor=CONFIG_FLAVOR,
        label="configuration",
    )
    _validate_settings()
    _validate_subset_flat(
        schema=INIT_SCHEMA,
        sections=BEGINNER_INIT_TEMPLATE_SECTIONS,
        label="beginner initialization",
    )
    _validate_subset_flat(
        schema=CONFIG_SCHEMA,
        sections=BEGINNER_CONFIG_TEMPLATE_SECTIONS,
        label="beginner configuration",
    )
    _validate_beginner_settings()


# ---------------------------------------------------------------------------
# INI rendering
# ---------------------------------------------------------------------------


def _render_ini_entry(
    key: str,
    meta: Mapping[str, Any],
    flavor: Mapping[str, Any],
    *,
    include_flavor: bool,
) -> str:
    if not include_flavor:
        return f"{key} ="

    description = str(flavor.get("description", "") or "").strip()
    example = str(flavor.get("example", "") or "").strip()
    legacy = flavor.get("legacy", {}) or {}

    lines: list[str] = []

    if description:
        lines.extend(_wrap_comment_line(description))

        allowed_range = _allowed_range_display(meta)
        if allowed_range:
            lines.extend(_wrap_comment_line(allowed_range))

        lines.append(f"# Default: {_format_value(meta.get('default'))}")

        if example:
            lines.append(f"# Example: {example}")

        lines.extend(_wrap_comment_line(_version_note_line(key, legacy)))

    lines.append(f"{key} =")
    return "\n".join(lines)


def _render_flat_ini_template(
    *,
    title: str,
    flavor_file: str | None,
    schema: Mapping[str, Mapping[str, Any]],
    sections: Mapping[str, Iterable[str]],
    flavor: Mapping[str, Mapping[str, Any]],
    include_flavor: bool,
) -> str:
    lines = [_template_header(title, flavor_file if include_flavor else None)]

    for section, keys in sections.items():
        lines.append(f"[{section}]")
        lines.append("")
        for key in keys:
            lines.append(
                _render_ini_entry(
                    key,
                    schema[key],
                    flavor[key],
                    include_flavor=include_flavor,
                )
            )
            lines.append("")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _render_settings_ini_template(
    *,
    title: str,
    flavor_file: str | None,
    keys_by_group: Mapping[str, Iterable[str]],
    include_flavor: bool,
) -> str:
    lines = [_template_header(title, flavor_file if include_flavor else None)]

    for group, keys in keys_by_group.items():
        section_name = SETTINGS_TEMPLATE_SECTIONS[group]
        lines.append(f"[{section_name}]")
        lines.append("")
        schema = SETTINGS_SCHEMA[group]
        flavor = SETTINGS_FLAVOR[group]
        for key in keys:
            lines.append(
                _render_ini_entry(
                    key,
                    schema[key],
                    flavor[key],
                    include_flavor=include_flavor,
                )
            )
            lines.append("")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------


def _escape_md(value: Any) -> str:
    text = _format_value(value) if not isinstance(value, str) else value
    text = html.escape(str(text))
    text = text.replace("|", "\\|")
    text = text.replace("\n", "<br>")
    return text


def _markdown_header(title: str, template_name: str, flavor_file: str) -> str:
    return "\n".join(
        [
            f"# {title}",
            "",
            '!!! note "Generated reference"',
            f"    This page was generated by ATLAS `{ATLAS_VERSION}` from the parser schema, `{flavor_file}`, and the explicit selection in `src/atlas_actris/templates/template_profiles.py`.",
            f"    The corresponding generated full INI template is `{template_name}`.",
            "",
        ]
    )


def _section_title(section: str) -> str:
    """Return a readable title for an INI section name."""

    return section.replace("_", " ").title()


def _section_anchor(section: str) -> str:
    """Return the MkDocs/Markdown anchor generated from a section heading."""

    return section.replace("_", "-").lower()


def _section_intro(section: str) -> str:
    """Return short narrative text for initialization reference sections."""

    intros = {
        "configuration": (
            "Controls how ATLAS obtains or uses the system configuration, "
            "including optional SCC HOI export behaviour."
        ),
        "explicit_paths": (
            "Defines the main input, configuration, settings, radiosonde, and "
            "output paths. Relative paths are interpreted with respect to the "
            "folder containing the initialization file."
        ),
        "general_options": (
            "Selects which QA tests, profile quicklooks, background plots, and "
            "VLDR plots are produced. When process is explicitly selected, "
            "empty process_qck/process_bgd values are derived from it, while "
            "process_dedicated_dark controls automatic drk_* companions."
        ),
        "filter_channels": (
            "Restricts the channels that are processed, either by selecting "
            "explicit ATLAS channel IDs or excluding groups of channels based "
            "on parts of the channel ID."
        ),
        "trimming_options": (
            "Controls signal trimming, overflow handling, temporal averaging "
            "options, and optional time slicing/exclusion of measurements."
        ),
        "explicit_folders": (
            "Overrides the default measurement folder names inside "
            "parent_folder. Leave these empty when the standard ATLAS folder "
            "structure is used."
        ),
        "parsing_options": (
            "Controls automatic telecover file distribution and custom "
            "radiosonde parsing metadata."
        ),
    }

    return intros.get(section, "Initialization parameters for this section.")


def _schema_summary_lines(meta: Mapping[str, Any]) -> list[str]:
    """Return compact Markdown bullet lines for schema metadata."""

    lines = [
        f"- **Type:** `{_escape_md(_type_display(meta))}`",
        f"- **Category:** `{_escape_md(_category_display(meta))}`",
        f"- **Default:** `{_escape_md(_format_value(meta.get('default')))}`",
    ]

    allowed_range = _allowed_range_display(meta)
    if allowed_range:
        lines.append(f"- **Limits / allowed values:** {_escape_md(allowed_range)}")

    return lines


def _render_parameter_details(
    key: str,
    meta: Mapping[str, Any],
    flavor_entry: Mapping[str, Any],
) -> str:
    """Render one parameter as a MkDocs details block."""

    description = _escape_md(flavor_entry.get("description", ""))
    example = str(flavor_entry.get("example", "") or "").strip()
    legacy_text = _escape_md(
        _legacy_sentence(flavor_entry.get("legacy", {}), for_docs=True)
    )

    lines = [
        f'<a id="{_escape_md(key).replace("_", "-")}"></a>',
        f'??? info "`{_escape_md(key)}`"',
        "",
    ]

    if description:
        lines.append(f"    {description}")
        lines.append("")

    for item in _schema_summary_lines(meta):
        lines.append(f"    {item}")

    if example:
        lines.append(f"    - **Example:** `{_escape_md(example)}`")

    lines.append(f"    - **Version history:** {legacy_text}")
    lines.append("")
    lines.append("    ```ini")
    lines.append(f"    {key} =")
    lines.append("    ```")

    return "\n".join(lines)


def _render_section_summary_table(
    schema: Mapping[str, Mapping[str, Any]],
    flavor: Mapping[str, Mapping[str, Any]],
) -> str:
    """Render a compact section summary table for quick scanning."""

    lines = [
        "| Parameter | Type | Default | Allowed / limits |",
        "| --- | --- | --- | --- |",
    ]

    for key, meta in schema.items():
        allowed_range = _allowed_range_display(meta)
        lines.append(
            "| "
            + " | ".join(
                [
                    f"[`{_escape_md(key)}`](#{_escape_md(key).replace('_', '-')})",
                    f"`{_escape_md(_type_display(meta))}`",
                    f"`{_escape_md(_format_value(meta.get('default')))}`",
                    _escape_md(allowed_range),
                ]
            )
            + " |"
        )

    return "\n".join(lines)


def _render_initialization_markdown() -> str:
    """Render the generated MkDocs page for the initialization INI reference.

    This page is intentionally written to ``initialization_reference.md`` so
    that ``initialization.md`` can remain a hand-written overview page.
    """

    lines = [
        _markdown_header(
            "ATLAS initialization file reference",
            "call_atlas.ini",
            "src/atlas_actris/templates/init_template_flavor.py",
        ),
        "## Purpose",
        "",
        "The initialization file tells ATLAS where the input data and metadata "
        "files are located, which QA tests, quicklooks, background plots, and "
        "VLDR plots should be produced, and which optional filtering, trimming, "
        "export, and radiosonde parsing "
        "settings should be applied for one processing run.",
        "",
        "This page is generated automatically from the initialization parser "
        "schema and the user-facing template flavor text. Edit "
        "`src/atlas_actris/templates/init_template_flavor.py` to change "
        "descriptions, examples, or version notes. Edit "
        "`src/atlas_actris/templates/template_profiles.py` to change which "
        "parameters are included or their ordering. Edit "
        "`src/atlas_actris/utils/parse_init_file.py` only when the technical "
        "schema itself changes.",
        "",
        "## How values are interpreted",
        "",
        "- Empty values mean that ATLAS will use the schema default or automatic "
        "selection when available.",
        "- List values can be separated with commas or semicolons.",
        "- Boolean values must be written as `True` or `False`.",
        "- Relative explicit paths are resolved relative to the folder containing "
        "the initialization file.",
        "- Measurement folder aliases in `explicit_folders` are resolved relative "
        "to `parent_folder`.",
        "",
        "## Sections",
        "",
    ]

    for section, keys in INIT_TEMPLATE_SECTIONS.items():
        lines.append(
            f"- [`{section}`](#{_section_anchor(section)}): "
            f"{_section_intro(section)} ({len(list(keys))} parameters)."
        )

    lines.append("")
    lines.append("## Parameter reference")
    lines.append("")

    for section, keys in INIT_TEMPLATE_SECTIONS.items():
        section_schema = {key: INIT_SCHEMA[key] for key in keys}
        section_flavor = {key: INIT_FLAVOR[key] for key in keys}

        lines.append(f'<a id="{_section_anchor(section)}"></a>')
        lines.append(f"### `{section}`")
        lines.append("")
        lines.append(_section_intro(section))
        lines.append("")
        lines.append(_render_section_summary_table(section_schema, section_flavor))
        lines.append("")

        for key in keys:
            lines.append(_render_parameter_details(key, INIT_SCHEMA[key], INIT_FLAVOR[key]))
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _render_md_table(schema: Mapping[str, Mapping[str, Any]], flavor: Mapping[str, Mapping[str, Any]]) -> str:
    lines = [
        "| Parameter | Description | Type | Category | Default | Limits / allowed values | Example | Version history |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]

    for key, meta in schema.items():
        entry = flavor[key]
        legacy_text = _legacy_sentence(entry.get("legacy", {}), for_docs=True)
        example = entry.get("example", "")
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{_escape_md(key)}`",
                    _escape_md(entry.get("description", "")),
                    f"`{_escape_md(_type_display(meta))}`",
                    _escape_md(_category_display(meta)),
                    f"`{_escape_md(_format_value(meta.get('default')))}`",
                    _escape_md(_allowed_range_display(meta)),
                    f"`{_escape_md(example)}`" if example else "",
                    _escape_md(legacy_text),
                ]
            )
            + " |"
        )

    return "\n".join(lines)


def _render_flat_markdown(
    *,
    title: str,
    template_name: str,
    flavor_file: str,
    schema: Mapping[str, Mapping[str, Any]],
    sections: Mapping[str, Iterable[str]],
    flavor: Mapping[str, Mapping[str, Any]],
) -> str:
    lines = [_markdown_header(title, template_name, flavor_file)]

    for section, keys in sections.items():
        section_schema = {key: schema[key] for key in keys}
        section_flavor = {key: flavor[key] for key in keys}
        lines.append(f"## `{section}`")
        lines.append("")
        lines.append(_render_md_table(section_schema, section_flavor))
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _render_settings_markdown() -> str:
    lines = [
        _markdown_header(
            "ATLAS settings file reference",
            "settings_file.ini",
            "src/atlas_actris/templates/settings_template_flavor.py",
        )
    ]

    for group, section_name in SETTINGS_TEMPLATE_SECTIONS.items():
        lines.append(f"## `{section_name}`")
        lines.append("")
        lines.append(_render_md_table(SETTINGS_SCHEMA[group], SETTINGS_FLAVOR[group]))
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# Output planning / CLI
# ---------------------------------------------------------------------------


def _profiles_to_generate(profile: str) -> tuple[str, ...]:
    if profile not in ALLOWED_PROFILES:
        raise ValueError(f"profile must be one of {ALLOWED_PROFILES}, got {profile!r}")
    if profile == "all":
        return ("full", "bare", "beginner")
    return (profile,)


def _profile_suffix(profile: str) -> str:
    if profile == "full":
        return ""
    return f"_{profile}"


def _profile_include_flavor(profile: str) -> bool:
    return profile in ("full", "beginner")


def _profile_title(base_title: str, profile: str) -> str:
    if profile == "full":
        return base_title
    if profile == "bare":
        return base_title.replace("template", "bare template")
    if profile == "beginner":
        return base_title.replace("template", "beginner template")
    return base_title


def _profile_init_sections(profile: str) -> Mapping[str, Iterable[str]]:
    if profile == "beginner":
        return BEGINNER_INIT_TEMPLATE_SECTIONS
    return INIT_TEMPLATE_SECTIONS


def _profile_config_sections(profile: str) -> Mapping[str, Iterable[str]]:
    if profile == "beginner":
        return BEGINNER_CONFIG_TEMPLATE_SECTIONS
    return CONFIG_TEMPLATE_SECTIONS


def _profile_settings_keys(profile: str) -> Mapping[str, Iterable[str]]:
    if profile == "beginner":
        return BEGINNER_SETTINGS_TEMPLATE_KEYS
    return SETTINGS_TEMPLATE_KEYS


def build_outputs(
    repo_root: str | Path | None = None,
    target: str = "all",
    profile: str = "all",
) -> dict[Path, str]:
    if target not in ALLOWED_TARGETS:
        raise ValueError(f"target must be one of {ALLOWED_TARGETS}, got {target!r}")

    validate_all()

    root = Path(repo_root).resolve() if repo_root is not None else _repo_root_from_this_file()
    template_dir = root / GENERATED_TEMPLATES_DIRECTORY
    docs_dir = root / "docs" / "generated"

    outputs: dict[Path, str] = {}

    if target in ("all", "ini"):
        for profile_name in _profiles_to_generate(profile):
            suffix = _profile_suffix(profile_name)
            include_flavor = _profile_include_flavor(profile_name)

            outputs[template_dir / f"call_atlas{suffix}.ini"] = _render_flat_ini_template(
                title=_profile_title("ATLAS initialization file template", profile_name),
                flavor_file="src/atlas_actris/templates/init_template_flavor.py",
                schema=INIT_SCHEMA,
                sections=_profile_init_sections(profile_name),
                flavor=INIT_FLAVOR,
                include_flavor=include_flavor,
            )
            outputs[template_dir / f"config_file{suffix}.ini"] = _render_flat_ini_template(
                title=_profile_title("ATLAS configuration file template", profile_name),
                flavor_file="src/atlas_actris/templates/config_template_flavor.py",
                schema=CONFIG_SCHEMA,
                sections=_profile_config_sections(profile_name),
                flavor=CONFIG_FLAVOR,
                include_flavor=include_flavor,
            )
            outputs[template_dir / f"settings_file{suffix}.ini"] = _render_settings_ini_template(
                title=_profile_title("ATLAS settings file template", profile_name),
                flavor_file="src/atlas_actris/templates/settings_template_flavor.py",
                keys_by_group=_profile_settings_keys(profile_name),
                include_flavor=include_flavor,
            )

    if target in ("all", "ini"):
        requested_profiles = _profiles_to_generate(profile)
        if "full" in requested_profiles:
            outputs[template_dir / "intercomparison.ini"] = render_intercomparison_ini()
        if "bare" in requested_profiles:
            outputs[template_dir / "intercomparison_bare.ini"] = render_intercomparison_bare_ini()

    if target in ("all", "docs"):
        outputs[docs_dir / "initialization_reference.md"] = _render_initialization_markdown()
        outputs[docs_dir / "configuration_reference.md"] = _render_flat_markdown(
            title="ATLAS configuration file reference",
            template_name="config_file.ini",
            flavor_file="src/atlas_actris/templates/config_template_flavor.py",
            schema=CONFIG_SCHEMA,
            sections=CONFIG_TEMPLATE_SECTIONS,
            flavor=CONFIG_FLAVOR,
        )
        outputs[docs_dir / "settings_reference.md"] = _render_settings_markdown()
        outputs[docs_dir / "intercomparison_reference.md"] = render_intercomparison_markdown()

    return outputs


def generate(
    repo_root: str | Path | None = None,
    target: str = "all",
    profile: str = "all",
    check: bool = False,
) -> dict[Path, str]:
    root = Path(repo_root).resolve() if repo_root is not None else _repo_root_from_this_file()
    outputs = build_outputs(repo_root=root, target=target, profile=profile)

    if check:
        changed: list[str] = []
        for path, content in outputs.items():
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                changed.append(str(path.relative_to(root)))

        if changed:
            raise SystemExit(
                "Generated ATLAS template/docs files are out of date:\n"
                + "\n".join(f"  - {name}" for name in changed)
                + "\nRun: atlas-generate-templates"
            )
        return outputs

    for path, content in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate ATLAS INI templates and MkDocs reference pages."
    )
    parser.add_argument(
        "--repo-root",
        default=None,
        help="Repository root. Defaults to the root inferred from this module.",
    )
    parser.add_argument(
        "--target",
        choices=ALLOWED_TARGETS,
        default="all",
        help="What to generate: all, ini, or docs. Default: all.",
    )
    parser.add_argument(
        "--profile",
        choices=ALLOWED_PROFILES,
        default="all",
        help="INI template profile to generate: all, full, bare, or beginner. Default: all.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Do not write files. Fail if generated files differ from committed files.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = Path(args.repo_root).resolve() if args.repo_root else _repo_root_from_this_file()
    outputs = generate(
        repo_root=root,
        target=args.target,
        profile=args.profile,
        check=args.check,
    )

    if args.check:
        print(f"Generated ATLAS {args.target} files are up to date.")
    else:
        print(f"Generated ATLAS {args.target} files:")
        for path in outputs:
            print(f"  - {path.relative_to(root)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
