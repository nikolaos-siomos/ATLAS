ATLAS atlas-test CLI and GitHub Actions update — 2026-11-12

Copy the included files into their matching paths in your 1.3.0 repository.
This patch adds the runner, its tests/documentation, its console entry point,
and updates the existing GitHub Actions test workflow. pyproject.toml also
includes the meteorological CSV packaging correction discussed previously.
Version numbers remain 1.3.0; no release or deployment is performed.

Register the new command and install pytest in your ATLAS environment:
  python -m pip install -e ".[dev]"

Run every repository test, including the smoke/integration test:
  atlas-test

List available test modules:
  atlas-test --list

Run one module:
  atlas-test --test channel_validation

Run one individual test:
  atlas-test --test tests/test_atlas_integration.py::test_atlas_smoke_test_script_runs

Run from another directory or point an installed command to the checkout:
  atlas-test --repo /path/to/atlas_dev

Pass extra pytest options after --:
  atlas-test -- --collect-only -q
  atlas-test -- -x -vv
  atlas-test -- -m "not integration"

The final example deliberately excludes the full integration smoke test;
no tests are excluded by the CLI by default. Existing dependency-based skips
remain controlled by the tests and reported by pytest.

The runner uses the current Python interpreter, starts pytest in a fresh
process with the checkout as the working directory, and passes the checkout
source paths to test children. This keeps CLI, direct Python calls and the
integration pipeline on the selected source version. Existing environment
values, including MPLBACKEND, are respected.

The repository tests and testing_pack are required; they are not added to
the installed wheel. This is a developer command for a source checkout.
The existing atlas-smoke-test and signal-viewer smoke commands remain.

GitHub Actions keeps Ubuntu 24.04 and Python 3.9–3.13. Its test job now runs
atlas-test instead of only atlas-smoke-test. The existing integration test
already invokes the smoke runner, so the smoke test runs once per matrix
entry. The documentation deployment workflow is not modified.

Validation completed:
  10 focused runner tests passed (selection, precedence, process arguments,
  exit status, dependency errors, repository discovery and environment).
  Built a wheel and verified the new console entry point, runner module and
  meteorological CSV. Invoked the runner from the unpacked wheel to list
  the full checkout's test modules. Checked workflow matrix and Python 3.9
  syntax compatibility.
  The complete numerical/integration suite could not run in this workspace:
  pytest and several ATLAS dependencies (including xarray, dask and netCDF4)
  are not installed. GitHub Actions installs .[dev] and runs the full suite.

Merge changes if you independently edited any included file after the
uploaded 1.3.0 archive. The patch does not change processing algorithms.
