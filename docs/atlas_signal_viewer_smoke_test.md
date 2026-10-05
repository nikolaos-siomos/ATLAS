# `atlas-signal-viewer-smoke-test`

Runs the signal-viewer workflow through the installed `atlas-signal-viewer` CLI and verifies that signal-viewer output files are created.

## Usage

```bash
atlas-signal-viewer-smoke-test
```

An explicit testing-pack path can be supplied in the same way:

```bash
atlas-signal-viewer-smoke-test /path/to/testing_pack
```

Useful options include:

```text
--ini FILE
--case-name NAME
--timeout SECONDS
--keep-output
```

By default, the signal-viewer smoke test deletes the generated `signal_viewer` and temporary `cache` folders after validation. Other ATLAS analysis products are left untouched.

To retain the generated viewer files for manual inspection, use:

```bash
atlas-signal-viewer-smoke-test --keep-output
```

This is particularly useful when checking generated HTML files, interactive plots, or viewer-specific failures. As with the main smoke test, retained outputs may consume significant disk space.

## Longer runs

Use a longer timeout on a slower machine or CI runner:

```bash
atlas-signal-viewer-smoke-test --timeout 3600
```

See also the [other smoke-test command](atlas_smoke_test.md).
