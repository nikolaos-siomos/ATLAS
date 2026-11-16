# ATLAS test suite

`atlas-test` runs the repository's pytest suite. By default this includes all test modules and the full smoke/integration test, which invokes `testing_pack/run_test.py` and validates generated outputs.

## Install and run

From your ATLAS source checkout, install the development dependencies to register the command:

```bash
python -m pip install -e ".[dev]"
atlas-test
```

Use the environment containing ATLAS and its dependencies. The command invokes pytest through that environment's Python interpreter and returns its exit status, so failures also fail CI jobs.

Desktop GUI checks that create windows are skipped by default, even when a display is available. The remaining GUI/model checks still run.

The `tests/` directory and smoke-test dataset remain in the repository. A PyPI-only installation does not include them. When running elsewhere, point to the checkout:

```bash
atlas-test --repo /path/to/atlas_dev
```

## Select tests

List the available test modules:

```bash
atlas-test --list
```

Run one module using its short name:

```bash
atlas-test --test channel_validation
```

The equivalent file forms are `--test test_channel_validation.py` and `--test tests/test_channel_validation.py`. Repeat `--test` to run several modules:

```bash
atlas-test -t channel_validation -t template_generation
```

Run one individual test using its pytest node ID:

```bash
atlas-test --test tests/test_atlas_integration.py::test_atlas_smoke_test_script_runs
```

Use pytest collection to see individual test IDs without executing them:

```bash
atlas-test -- --collect-only -q
```

## Optional desktop GUI checks

A normal `atlas-test` run does not open the ATLAS editor or folder browser for
testing. To explicitly run the checks that create real desktop windows:

```bash
atlas-test --gui --test gui_compatibility
```

These checks require a working graphical display, open test windows, and close
them automatically. Their startup-folder prompt is disabled, so no interactive
file selection is required. Other GUI tests continue to run in the default suite.

`--gui` belongs to `atlas-test`, so put it before the `--` separator. The CLI
passes `ATLAS_RUN_GUI_TESTS=1` to pytest only for this explicit flag; without it,
it sets the child value to `0` even if the parent has a different value. When
invoking pytest directly, set that environment variable before collection to
opt into the desktop tests. No pytest plugin is required.

## Additional pytest options

Pass pytest options after `--`:

```bash
atlas-test -- -x -vv
atlas-test -- -k "channel"
atlas-test -- -m "not integration"
```

The last example skips tests marked `integration`, including the full smoke test. No such exclusion is applied by default. Dependency-based skips in existing tests are still reported by pytest.

## GitHub Actions

The checked-in `.github/workflows/tests.yml` invokes `atlas-test` on
Python 3.9–3.13 using Ubuntu 24.04. This runs the repository's `tests/` suite,
including the smoke integration test, PNG safety checks, HTML escaping checks,
cleanup-path checks, and the GUI compatibility checks that do not need windows.
Desktop checks remain skipped unless explicitly enabled with `--gui`.

The integration smoke test removes only its selected case's analysis output
after checking it. Other cases and caches in the measurement folders are left
untouched. Its testing dataset must be present in the checkout.
