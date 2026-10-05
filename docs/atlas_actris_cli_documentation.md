# ATLAS ACTRIS Command-Line Interface

This page describes the main `atlas` command during local developer installation.

The package is currently intended to be installed from the local repository with:

```bash
python -m pip install -e .
```

After installation, the main workflow is available through `atlas`. Additional commands are documented separately under [Extra tools](index.md#extra-tools).

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

### Help

```bash
atlas --help
```

The parser exposes the following options (the displayed program name and wrapping can vary by launcher and Python version):

```text
usage: atlas [-h] [-i [ini_file]] [-o [output_folder]]
             [-s slice_measurement [slice_measurement ...]]
             [-e exclude_measurement [exclude_measurement ...]]
             [-q process_qck [process_qck ...]]

arguments

options:
  -h, --help            show this help message and exit
  -i [ini_file], --ini_file [ini_file]
                        The path to the initialization file of the calling
                        script
  -o [output_folder], --output_folder [output_folder]
                        The path to the output folder. It overrides the
                        corresponding parameter provided in the initialization
                        file
  -s slice_measurement [slice_measurement ...], --slice_measurement slice_measurement [slice_measurement ...]
                        Slicing option for the QA tests. It overrides the
                        corresponding parameter provided in the initialization
                        file
  -e exclude_measurement [exclude_measurement ...], --exclude_measurement exclude_measurement [exclude_measurement ...]
                        Exclude option for the QA tests. It overrides the
                        corresponding parameter provided in the initialization
                        file
  -q process_qck [process_qck ...], --process_qck process_qck [process_qck ...]
                        Select QA test aliases for quicklook generation. It
                        overrides the corresponding parameter provided in the
                        initialization file
```

### Arguments

| Option | Value | Purpose |
| --- | --- | --- |
| `-i`, `--ini_file` | One file path | Initialization file. Required for execution; the path must exist. |
| `-o`, `--output_folder` | One folder path | Override `output_folder`. |
| `-s`, `--slice_measurement` | One or more `measurement start stop` triplets | Override the time windows selected for processing. |
| `-e`, `--exclude_measurement` | One or more `measurement start stop` triplets | Override the time windows excluded from processing. |
| `-q`, `--process_qck` | One or more quicklook aliases | Override the quicklook selection. |
| `-h`, `--help` | No value | Show help and exit. |

Although argparse displays `-i` as optional, the parser checks that an initialization file was supplied. Always supply a value with `-i` and `-o`; do not use them as bare switches.

Only the options listed here are exposed by `parse_caller_args.py`. Other initialization-schema parameters must still be set in the INI file; schema membership alone does not create a CLI option.

### Override precedence

ATLAS reads the initialization file, applies supplied CLI overrides, then performs schema conversion and validation. An omitted option leaves the INI setting in effect, subject to the usual defaults and processing rules. The initialization file itself is not modified.

List overrides **replace the entire INI value**, rather than append to it. For example, `-q drk ray` replaces the INI's `process_qck` list. The same replacement rule applies to `-s` and `-e`. Supplied overrides are listed in the terminal under `Command-line overrides:`.

Pass multiple values after a single option. Repeating an option keeps only its last occurrence. Space-separated values are easiest to read; quoted comma- or semicolon-separated lists are also accepted by the initialization parser.

### Output folder

```bash
atlas -i ./call_atlas.ini -o ./analysis_run
```

This overrides the main workflow's `output_folder` setting; it does not change any other initialization parameters.

Quote paths containing spaces, including on Windows:

```powershell
atlas -i "C:\ATLAS runs\call_atlas.ini" -o "C:\ATLAS runs\analysis"
```

### Select or exclude time windows

Each entry consists of three values: a measurement identifier, a start time, and a stop time.

Select a Rayleigh measurement window:

```bash
atlas -i ./call_atlas.ini -s ray 1200 1300
```

Select windows for multiple measurements with one `-s` option:

```bash
atlas -i ./call_atlas.ini -s ray 1200 1300 pcb 1400 1430
```

Exclude a time window:

```bash
atlas -i ./call_atlas.ini -e ray 1210 1215
```

Selection and exclusion can be used together:

```bash
atlas -i ./call_atlas.ini -s ray 1200 1300 -e ray 1210 1215
```

For an interval crossing midnight, explicit dates make the intended interval clear:

```bash
atlas -i ./call_atlas.ini -s ray 20261005_2330 20261006_0100
```

Accepted time formats are `HHMM`, `YYYYMMDD`, `YYYYMMDD_HH`, `YYYYMMDD_HHMM`, and `YYYYMMDD_HHMMSS`. Date-only values denote midnight; omitted minutes and seconds are zero. With `HHMM`, the measurement date is assigned later in processing.

Both slicing and exclusion accept individual measurement identifiers:

- `ray`, `ray_pcb`, `trg`, `dtm`, `drk`;
- `pcb_p45`, `pcb_m45`, `pcb_aux_p45`, `pcb_aux_m45`;
- `tlc_north`, `tlc_east`, `tlc_south`, `tlc_west`, `tlc_inner`, `tlc_outer`;
- `drk_ray`, `drk_pcb`, `drk_tlc`, `drk_tlc_rin`, `drk_trg`, `drk_dtm`, `drk_ray_pcb`, `drk_pcb_aux`.

The following bundle aliases apply the same interval to each listed measurement:

| Alias | Measurements |
| --- | --- |
| `pcb` | `pcb_p45`, `pcb_m45` |
| `pcb_aux` | `pcb_aux_p45`, `pcb_aux_m45` |
| `tlc` | `tlc_north`, `tlc_east`, `tlc_south`, `tlc_west` |
| `tlc_rin` | `tlc_inner`, `tlc_outer` |

For example, `-s tlc 1200 1300` applies that interval to all four telecover sectors. Bundle expansion here is specific to slicing and exclusion.

### Select quicklooks

```bash
atlas -i ./call_atlas.ini -q drk ray
```

Equivalent comma-separated form:

```bash
atlas -i ./call_atlas.ini --process_qck "drk, ray"
```

Accepted aliases are `ray`, `pcb`, `tlc`, `tlc_rin`, `ray_pcb`, `pcb_aux`, `trg`, `dtm`, `drk`, `drk_ray`, `drk_pcb`, `drk_tlc`, `drk_tlc_rin`, `drk_ray_pcb`, `drk_pcb_aux`, `drk_trg`, and `drk_dtm`.

To disable quicklooks explicitly:

```bash
atlas -i ./call_atlas.ini -q off
```

Use `off` alone. `-q` changes `process_qck`; it does not replace the `process` selection for the other QA tests. `ray` and `ray_pcb` are distinct aliases.

### Combine overrides

```bash
atlas -i ./call_atlas.ini -o ./analysis_run -s ray 1200 1300 -e ray 1210 1215 -q drk ray
```

The long-option equivalent is:

```bash
atlas --ini_file ./call_atlas.ini --output_folder ./analysis_run --slice_measurement ray 1200 1300 --exclude_measurement ray 1210 1215 --process_qck drk ray
```

See the [initialization file reference](generated/initialization_reference.md) for the corresponding schema parameters.

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
