# ATLAS ACTRIS Command-Line Interface

This page describes the current command-line interface exposed by the `atlas_actris` package during local developer installation.

The package is currently intended to be installed from the local repository with:

```bash
python -m pip install -e .
```

After installation, the main workflow is available through `atlas`. Additional installed commands provide the signal viewer and smoke tests, as documented below.

The main command expects an initialization file path:

```bash
atlas -i /path/to/call_atlas.ini
```

or:

```bash
atlas --ini_file /path/to/call_atlas.ini
```

## Current CLI entry point

The currently supported console script is defined in `pyproject.toml` under `[project.scripts]`:

```toml
[project.scripts]
atlas = "atlas_actris.cli:main"
```

The command calls the `main()` function in:

```text
src/atlas_actris/cli.py
```

The CLI wrapper then runs:

```text
src/atlas_actris/__call_atlas_interactive__.py
```

This design keeps the main ATLAS script runnable in two ways:

1. from the terminal as the installed command `atlas`, and
2. directly in Spyder by opening and running `src/atlas_actris/__call_atlas_interactive__.py`.

The second workflow is useful for development because variables created by the script remain visible in Spyder's Variable Explorer.

## Command: `atlas`

Runs the main ATLAS interactive workflow.

### Usage

```bash
atlas -i /path/to/call_atlas.ini
```

Equivalent long-option form:

```bash
atlas --ini_file /path/to/call_atlas.ini
```

### Help output

The current help output is:

```text
usage: __call_atlas_interactive__.py [-h] [-i [ini_file]]

arguments

optional arguments:
  -h, --help            show this help message and exit
  -i [ini_file], --ini_file [ini_file]
                        The path to the initialization file of the calling
                        script
```

### Arguments

- `-i`, `--ini_file`

  Path to the initialization file of the calling script.

  Example:

  ```bash
  atlas -i ./path/to/call_atlas.ini
  ```

  The current argparse help displays this option as optional because of the `[-i [ini_file]]` form. However, the ATLAS workflow expects an initialization file path for normal execution.

### Examples

Run ATLAS with an initialization file:

```bash
atlas -i ./docs/templates/call_atlas.ini
```

or:

```bash
atlas --ini_file ./docs/templates/call_atlas.ini
```

Use an absolute path if the initialization file is outside the repository:

```bash
atlas -i /home/user/atlas_runs/call_atlas.ini
```

On Windows:

```powershell
atlas -i C:\Users\YourUser\atlas_runs\call_atlas.ini
```

Check help:

```bash
atlas -h
```

or:

```bash
atlas --help
```


## Additional installed commands

The package also exposes commands for the signal viewer and for automated smoke testing:

```toml
[project.scripts]
atlas = "atlas_actris.cli:main"
atlas-signal-viewer = "atlas_actris.cli_viewer:main"
atlas-intercomparison = "atlas_actris.cli_intercomparison:main"
atlas-smoke-test = "atlas_actris.testing.smoke_test:main"
atlas-signal-viewer-smoke-test = "atlas_actris.testing.signal_viewer_smoke_test:main"
```

After changing `[project.scripts]`, reinstall the package in editable mode so the new commands are created in the active environment:

```bash
python -m pip install -e .
```


## Command: `atlas-intercomparison`

Runs the ATLAS exported-stage intercomparison workflow. The command expects one intercomparison initialization file.

### Usage

```bash
atlas-intercomparison -i /path/to/intercomparison.ini
```

Equivalent long-option form:

```bash
atlas-intercomparison --ini_file /path/to/intercomparison.ini
```

The command dispatches through:

```text
src/atlas_actris/cli_intercomparison.py
```

and runs the main interactive script:

```text
src/atlas_actris/__intercomparison_interactive__.py
```

The same main script can be opened and run directly in Spyder. Configure the Spyder run arguments with:

```text
-i /path/to/intercomparison.ini
```

The parsed initialization dictionary remains available in the Spyder namespace as `intercomparison_info`.

### Intercomparison templates

The generated full and bare templates are:

```text
src/atlas_actris/templates/intercomparison.ini
src/atlas_actris/templates/intercomparison_bare.ini
```

The full template contains parameter descriptions, defaults, allowed values, and examples, but all assignments are intentionally empty. The bare template contains the same section and parameter structure without per-parameter flavor comments.

The generated reference page is:

```text
docs/generated/intercomparison_reference.md
```

## Command: `atlas-signal-viewer`

Runs the ATLAS signal-viewer workflow using the same initialization-file interface as the main ATLAS command.

### Usage

```bash
atlas-signal-viewer -i /path/to/call_atlas.ini
```

Equivalent long-option form:

```bash
atlas-signal-viewer --ini_file /path/to/call_atlas.ini
```

The signal viewer reads and processes the configured test data, generates the requested interactive viewer outputs, and stores them under the case output folder in:

```text
analysis/<case-name>/signal_viewer/
```

The command is interactive during normal use. At the end it may ask whether temporary cache files and generated signal-viewer files should be deleted.

## Command: `atlas-smoke-test`

Runs the packaged ATLAS smoke-test dataset through the main `atlas` CLI and validates that the expected plots, reports, and ASCII products are created.

### Usage

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

### Keeping generated outputs

Use `--keep-output` when debugging or manually inspecting the generated files:

```bash
atlas-smoke-test --keep-output
```

With this option, existing outputs are not removed before the test and newly generated outputs are retained afterwards. This is useful when checking plots, reports, cache contents, or a failing intermediate result. Because these files can be large, omit `--keep-output` during routine testing and continuous integration.

## Command: `atlas-signal-viewer-smoke-test`

Runs the signal-viewer workflow through the installed `atlas-signal-viewer` CLI and verifies that signal-viewer output files are created.

### Usage

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

## Smoke-test examples

Run both smoke tests with automatic cleanup:

```bash
atlas-smoke-test
atlas-signal-viewer-smoke-test
```

Retain outputs from the main test:

```bash
atlas-smoke-test --keep-output
```

Retain signal-viewer outputs:

```bash
atlas-signal-viewer-smoke-test --keep-output
```

Use a longer timeout on a slower machine or CI runner:

```bash
atlas-smoke-test --timeout 3600
atlas-signal-viewer-smoke-test --timeout 3600
```


## Generating INI templates and reference pages

Generate all full, bare, and beginner ATLAS templates together with the generated documentation pages:

```bash
atlas-generate-templates
```

The intercomparison generator participates in the same command. With the default `--profile all`, the command writes both:

```text
src/atlas_actris/templates/intercomparison.ini
src/atlas_actris/templates/intercomparison_bare.ini
```

and:

```text
docs/generated/intercomparison_reference.md
```

Generate only full templates:

```bash
atlas-generate-templates --target ini --profile full
```

Generate only bare templates:

```bash
atlas-generate-templates --target ini --profile bare
```

Generate only the documentation pages:

```bash
atlas-generate-templates --target docs
```

Check whether committed generated files are current without rewriting them:

```bash
atlas-generate-templates --check
```

After changing parser schemas or flavor files, rerun the generator from the repository root. Editable reinstallation is only required when console-script definitions in `pyproject.toml` change.

## Developer execution from Python

The installed command is equivalent in purpose to running the main script from the package source tree.

From Spyder, open and run:

```text
src/atlas_actris/__call_atlas_interactive__.py
```

Configure the Spyder run arguments to pass the initialization file:

```text
-i /path/to/call_atlas.ini
```

From the terminal, after editable installation, use:

```bash
atlas -i /path/to/call_atlas.ini
```

or run the file directly from the repository if needed:

```bash
python src/atlas_actris/__call_atlas_interactive__.py -i /path/to/call_atlas.ini
```

The recommended terminal command is still `atlas -i ...`, because it tests the installed package entry point.

## Current `cli.py` wrapper

The current CLI wrapper can be written as:

```python
from pathlib import Path
import runpy
import sys


def main():
    script_path = Path(__file__).with_name("__call_atlas_interactive__.py")
    script_dir = str(script_path.parent.resolve())

    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    runpy.run_path(str(script_path), run_name="__main__")
```

This wrapper exists only to provide a callable target for the installed console script. The main workflow itself remains in `__call_atlas_interactive__.py`.

Because the wrapper uses `runpy.run_path(..., run_name="__main__")`, command-line arguments passed to `atlas` remain available to the underlying script through `sys.argv`. Therefore:

```bash
atlas -i /path/to/call_atlas.ini
```

is passed through to `__call_atlas_interactive__.py` as the script's command-line input.

## Legacy or planned commands

Older documentation mentioned additional commands:

- `get_config`
- `get_T_P`

These commands are not currently exposed as installed console scripts unless they are added to `[project.scripts]` in `pyproject.toml`.

If these tools should become installed commands again, add suitable callable functions in the package and expose them explicitly. For example:

```toml
[project.scripts]
atlas = "atlas_actris.cli:main"
get_config = "atlas_actris.__get_config_file_from_scc_hoi__:main"
get_T_P = "atlas_actris.__get_T_P_profiles_from_cloudnet__:main"
```

Only add entries like these if the target modules actually define a callable `main()` function and are ready to be used as stable command-line tools.

Until then, treat `get_config` and `get_T_P` as legacy or internal developer scripts, not as part of the current installed CLI.

## Verifying the CLI after installation

From the project root, install the package in editable mode:

```bash
python -m pip install -e .
```

Then verify that the package imports from the local repository:

```bash
python -c "import atlas_actris; print(atlas_actris.__file__)"
```

The printed path should point under:

```text
src/atlas_actris/
```

Then verify the command help:

```bash
atlas -h
```

Then run the command with an initialization file:

```bash
atlas -i /path/to/call_atlas.ini
```

## Troubleshooting

### `atlas: command not found`

The package is probably not installed in the active environment. Activate the correct environment and reinstall:

```bash
conda activate atlas_box
python -m pip install -e .
```

Then check:

```bash
which atlas
```

On Windows:

```bash
where atlas
```

### `ModuleNotFoundError: No module named 'atlas_actris'`

Make sure you are using the environment where the package was installed:

```bash
python -c "import sys; print(sys.executable)"
python -m pip show atlas_actris
```

If needed, reinstall from the repository root:

```bash
python -m pip install -e .
```

### The command runs but does not find the initialization file

Use an absolute path first to rule out working-directory issues:

```bash
atlas -i /absolute/path/to/call_atlas.ini
```

If that works, the problem was the relative path being interpreted from a different current working directory than expected.

### Imports work in Spyder but not in the terminal, or the reverse

Spyder and the terminal are probably using different Python environments.

Start Spyder from the activated ATLAS environment:

```bash
conda activate atlas_box
spyder
```

Then check the interpreter inside Spyder.
