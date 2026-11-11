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
      - checks that plots/reports/ascii outputs were created,
      - removes generated outputs after checking them.

    A zero exit status includes output validation. Do not inspect the output
    folder afterward: the smoke runner deliberately cleans it in finally.
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

