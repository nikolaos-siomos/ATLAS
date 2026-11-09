"""Render the existing application CSV in MkDocs without copying it into docs/."""
import csv
import html
from pathlib import Path


DATABASE_PATH = Path("src/atlas_actris/molecular/meteorological_database.csv")
PAGE_PATH = "meteorological_database.md"
MARKER = "<!-- ATLAS_METEOROLOGICAL_DATABASE -->"
COLUMNS = (
    "station_id", "station_name", "cloudnet_station_name",
    "rsonde_station_wmo_id", "rsonde_station_name",
    "radiosonde_distance_km", "rsonde_station_altitude",
)


def _table_cell(value):
    # Keep identifiers/precision as CSV text and escape Markdown/HTML syntax.
    text = " ".join((value or "").splitlines())
    text = html.escape(text, quote=False)
    for character in ("|", "`", "*", "_", "[", "]", "\\"):
        text = text.replace(character, f"&#{ord(character)};")
    return text


def render_database_table(path):
    """Read the current CSV on every call, preserving row order and blank cells."""
    try:
        with Path(path).open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            missing = set(COLUMNS).difference(reader.fieldnames or [])
            if missing:
                raise ValueError(f"Meteorological database {path} is missing columns: {sorted(missing)}")
            rows = list(reader)
    except OSError as exc:
        raise RuntimeError(f"Cannot read meteorological database {path}: {exc}") from exc

    lines = ["| " + " | ".join(f"`{key}`" for key in COLUMNS) + " |",
             "| " + " | ".join("---" for _ in COLUMNS) + " |"]
    lines.extend("| " + " | ".join(_table_cell(row.get(key)) for key in COLUMNS) + " |"
                 for row in rows)
    return "\n".join(lines)


def on_page_markdown(markdown, page, **kwargs):
    """MkDocs built-in hook; loaded by the hooks entry in mkdocs.yml."""
    if page.file.src_uri != PAGE_PATH:
        return markdown
    if markdown.count(MARKER) != 1:
        raise ValueError(f"{PAGE_PATH} must contain exactly one {MARKER} marker")
    project_root = Path(__file__).resolve().parents[1]
    table = render_database_table(project_root / DATABASE_PATH)
    return markdown.replace(MARKER, table)
