# ATLAS ACTRIS Command-Line Interface

This page describes the current command-line interface exposed by the `atlas_actris` package during local developer installation.

The package is currently intended to be installed from the local repository with:

```bash
python -m pip install -e .
```

After installation, the main command is:

```bash
atlas
```

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
