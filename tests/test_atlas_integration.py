from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.integration
def test_atlas_smoke_test_script_runs():
    """
    Run the portable ATLAS smoke-test script.

    The script itself:
      - finds testing_pack relative to its own location,
      - removes old analysis/cache outputs,
      - runs `atlas -i ...`,
      - answers the two terminal prompts,
      - checks that plots/reports/ascii outputs were created.
    """
    repo_root = Path(__file__).resolve().parents[1]
    smoke_test = repo_root / "testing_pack" / "run_test.py"

    assert smoke_test.exists(), f"Missing smoke-test script: {smoke_test}"

    result = subprocess.run(
        [sys.executable, str(smoke_test)],
        cwd=repo_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=2400,
    )

    assert result.returncode == 0, (
        f"ATLAS smoke test failed with return code {result.returncode}\n\n"
        f"STDOUT:\n{result.stdout}\n\n"
        f"STDERR:\n{result.stderr}\n"
    )


@pytest.mark.integration
def test_atlas_smoke_test_outputs_exist_after_run():
    """
    Check the expected smoke-test outputs.

    This assumes `test_atlas_smoke_test_script_runs` has just run.
    It is intentionally light: it does not validate scientific values yet.
    """
    repo_root = Path(__file__).resolve().parents[1]
    output_root = (
        repo_root
        / "testing_pack"
        / "analysis"
        / "179_199_665_20231221"
    )

    assert output_root.exists(), f"Missing ATLAS output folder: {output_root}"

    plots = list((output_root / "plots").glob("*.png"))
    html_reports = list((output_root / "reports").glob("*.html"))
    docx_reports = list((output_root / "reports").glob("*.docx"))
    ascii_files = list((output_root / "ascii").rglob("*.txt"))

    assert len(plots) >= 10, f"Expected at least 10 plots, found {len(plots)}"
    assert len(html_reports) >= 1, "Expected at least one HTML report"
    assert len(docx_reports) >= 1, "Expected at least one DOCX report"
    assert len(ascii_files) >= 5, f"Expected at least 5 ASCII files, found {len(ascii_files)}"
