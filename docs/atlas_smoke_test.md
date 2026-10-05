# `atlas-smoke-test`

Runs the packaged ATLAS smoke-test dataset through the main `atlas` CLI and validates that the expected plots, reports, and ASCII products are created.

## Usage

From the repository root:

```bash
atlas-smoke-test
```

The command can also be run from another directory, provided the package is installed in the active environment. The test first looks for `./testing_pack`; if that folder is not present, it falls back to the repository-level `testing_pack` associated with the editable installation.

An explicit testing-pack path can also be supplied:

```bash
atlas-smoke-test /path/to/testing_pack
```

Useful options include:

```text
--ini FILE
--case-name NAME
--timeout SECONDS
--keep-output
```

By default, the smoke test removes generated analysis files when it finishes. This keeps repeated local and CI test runs from accumulating large output folders. Cleanup also occurs after a failed run where possible.

## Keeping generated outputs

Use `--keep-output` when debugging or manually inspecting the generated files:

```bash
atlas-smoke-test --keep-output
```

With this option, existing outputs are not removed before the test and newly generated outputs are retained afterwards. This is useful when checking plots, reports, cache contents, or a failing intermediate result. Because these files can be large, omit `--keep-output` during routine testing and continuous integration.

## Longer runs

Use a longer timeout on a slower machine or CI runner:

```bash
atlas-smoke-test --timeout 3600
```

See also the [other smoke-test command](atlas_signal_viewer_smoke_test.md).
