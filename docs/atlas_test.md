# ATLAS test suite

`atlas-test` runs the repository's pytest suite. By default this includes all test modules and the full smoke/integration test, which invokes `testing_pack/run_test.py` and validates generated outputs.

## Install and run

From your ATLAS source checkout, install the development dependencies to register the command:

```bash
python -m pip install -e ".[dev]"
atlas-test
```

Use the environment containing ATLAS and its dependencies. The command invokes pytest through that environment's Python interpreter and returns its exit status, so failures also fail CI jobs.

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

## Additional pytest options

Pass pytest options after `--`:

```bash
atlas-test -- -x -vv
atlas-test -- -k "channel"
atlas-test -- -m "not integration"
```

The last example skips tests marked `integration`, including the full smoke test. No such exclusion is applied by default. Dependency-based skips in existing tests are still reported by pytest.

## GitHub Actions

The test workflow installs `.[dev]` and runs `atlas-test` on Python 3.9–3.13 using Ubuntu 24.04. The full suite includes the existing smoke/integration test, so the workflow does not launch a second separate smoke run. The documentation-deployment workflow is unchanged.

The integration smoke test removes its generated analysis output after checking it, following the existing smoke-runner behavior. Its testing dataset must be present in the checkout.
