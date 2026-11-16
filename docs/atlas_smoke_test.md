# `atlas-smoke-test`

Runs the repository's ATLAS smoke-test dataset through the main `atlas` CLI and validates that the expected plots, reports, and ASCII products are created.

The `testing_pack/` dataset is distributed with the source checkout, not inside
the installed PyPI wheel. A PyPI-only installation needs a checkout or a separately
available testing pack.

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

| Option | Purpose / default |
| --- | --- |
| `--ini FILE` | Init filename within the testing pack. |
| `--case-name NAME` | A single case directory name under `analysis/`; paths are rejected. |
| `--timeout SECONDS` | Run timeout; default `1800`. |
| `--keep-output` | Skip the runner's pre/post-run removal of the selected case's analysis output. |
| `--yes-clean-cache y/n` | Answer ATLAS's cache-cleanup prompt; default `y`. |
| `--yes-export-stage y/n` | Answer its stage-export prompt; default `N`. |

The runner always answers **No** to deletion of existing stage exports during the
ATLAS run. The smoke runner's own default cleanup still removes the selected
`analysis/NAME` case folder afterwards.

By default, the smoke test removes the selected case's old output before running
and its generated output after validation, including after a failed run where
possible. Other cases under `analysis/` and caches in raw measurement folders
are preserved. The `analysis/` parent is removed only if it is empty.

Cleanup rejects traversal paths, symbolic links and directory junctions beneath
the chosen testing pack, rather than following them into other directories.

## Keeping generated outputs

Use `--keep-output` when debugging or manually inspecting the generated files:

```bash
atlas-smoke-test --keep-output
```

With this option, the smoke runner does not remove the selected case folder before
or after the test. ATLAS still receives `--yes-clean-cache y` by default, so its
cache can be deleted during processing. To retain that cache too:

```bash
atlas-smoke-test --keep-output --yes-clean-cache n
```

Plots, reports and ASCII outputs are retained when the run succeeds. This is useful when checking plots, reports, cache contents, or a failing intermediate result. Because these files can be large, omit `--keep-output` during routine testing and continuous integration.

## Longer runs

Use a longer timeout on a slower machine or CI runner:

```bash
atlas-smoke-test --timeout 3600
```

See also the [other smoke-test command](atlas_signal_viewer_smoke_test.md).
