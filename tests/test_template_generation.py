from __future__ import annotations

from pathlib import Path

from atlas_actris.templates.generator import build_outputs, validate_all


def test_template_flavors_are_synchronized_with_schemas():
    validate_all()


def test_generated_template_files_are_up_to_date():
    repo_root = Path(__file__).resolve().parents[1]
    outputs = build_outputs(repo_root=repo_root, target="all", profile="all")

    missing_or_changed = []
    for path, expected_content in outputs.items():
        if not path.exists():
            missing_or_changed.append(str(path.relative_to(repo_root)))
            continue
        if path.read_text(encoding="utf-8") != expected_content:
            missing_or_changed.append(str(path.relative_to(repo_root)))

    assert not missing_or_changed, (
        "Generated template/docs files are missing or out of date:\n"
        + "\n".join(f"  - {name}" for name in missing_or_changed)
        + "\nRun: atlas-generate-templates"
    )
