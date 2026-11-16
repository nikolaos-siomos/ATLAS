# ATLAS ACTRIS

**[Read the ATLAS documentation](https://nikolaos-siomos.github.io/ATLAS/)** — installation, command-line tools, and INI file references.

[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/)
[![PyPI version](https://img.shields.io/pypi/v/atlas-actris.svg)](https://pypi.org/project/atlas-actris/)
[![Development status](https://img.shields.io/badge/status-development-orange)](https://github.com/nikolaos-siomos/ATLAS)
[![Tests](https://github.com/nikolaos-siomos/ATLAS/actions/workflows/tests.yml/badge.svg)](https://github.com/nikolaos-siomos/ATLAS/actions/workflows/tests.yml)
[![License](https://img.shields.io/badge/license-see%20LICENSE-lightgrey)](LICENSE)
Automated Lidar Analysis Software (ATLAS) for lidar data processing, visualization, and quality-assurance workflows.

ATLAS is intended for users involved in ACTRIS aerosol high-power lidar (AHL) quality-assurance activities. The primary users are the ACTRIS CARS group, who use ATLAS for the evaluation of quality-assurance tests of AHL systems, and National Facility PIs, who use ATLAS to check and monitor the status of their systems.


## Requirements

ATLAS currently supports Python **3.9, 3.10, 3.11, 3.12, and 3.13**.

The repository contains unit tests and a full smoke/integration dataset.
`atlas-test` runs the repository suite by default. The checked-in GitHub Actions
workflow currently runs the smoke test only; it checks the command-line workflow
and generated plots, reports, and ASCII outputs.

## Installation for users
Create and activate a Python environment. For example, using conda or miniforge:

```bash
conda create -n atlas_box python=3.13 pip
conda activate atlas_box
```

Install ATLAS through PyPI by typing the following command inside the active environment:

```bash
pip install atlas-actris
```

## Installation for developers
Create and activate a Python environment. For example, using conda or miniforge:

```bash
conda create -n atlas_box python=3.13 pip
conda activate atlas_box
```

Clone or just download and extract the repository:

```bash
git clone https://github.com/nikolaos-siomos/ATLAS.git
```

Enter the top level of the cloned/extracted project folder (the folder name can be different than in the example).
```bash
cd ATLAS
```

Install ATLAS locally by typing the following command inside the project folder:

```bash
python -m pip install -e .
```

This installs the package and exposes the `atlas` command in the active environment. Using the -e option keeps the source tree editable while working with ATLAS locally.

## Minimal usage
Activate the Python environment where ATLAS is installed. For example, using conda or miniforge:

```bash
conda activate atlas_box
```

Run ATLAS from the command line by providing the path to the initialization file:

```bash
atlas -i /path/to/call_atlas.ini
```

or equivalently:

```bash
atlas --ini_file /path/to/call_atlas.ini
```

To check the available command-line options type:

```bash
atlas -h
```

ATLAS requires an initialization file and a valid configuration. With empty init
paths, input data and `config_file.ini` are taken from beside the init; an optional
`settings_file.ini` is used there when present. Otherwise settings use built-in
defaults. Outputs default to `analysis/` beside the init. Explicit paths override
these defaults; an explicit output base gets a subfolder named after the input
parent folder. See the [quick guide](docs/quick_guide.md) for details.

To use the GUI in the same environment:

```bash
atlas-gui
```

For command-line execution, the `atlas` entry point uses a non-interactive plotting backend by default, which is suitable for terminal runs and automated tests. IDE workflows can still use the plotting backend configured in the active environment.

## Run ATLAS from an IDE
It is also possible to run atlas by opening the following directly in an IDE (e.g. Spyder):

```text
src/atlas_actris/__call_atlas_interactive__.py
```
Please make sure the IDE uses the same Python environment where ATLAS was installed.

## Smoke testing

The repository includes a smoke-test dataset in `testing_pack/`; the PyPI wheel
does not bundle this dataset. This test is intended to confirm that the complete ATLAS command-line workflow still works after code, configuration, or dependency changes.

Run the smoke test directly with:

```bash
atlas-smoke-test testing_pack
```

Alternatively, run the portable test script directly:

```bash
python testing_pack/run_test.py
```

The smoke test:

- removes previous generated analysis/cache outputs,
- runs ATLAS with the bundled initialization file,
- supplies answers for cache cleanup, stage deletion and stage export,
- checks that plots, reports, and ASCII files were created.

Generated outputs are written under:

```text
testing_pack/analysis/
```

This folder is ignored by Git and should not be committed.

## Running the test suite

Install the development dependencies first:

```bash
python -m pip install -e ".[dev]"
```

Run the whole repository suite:

```bash
atlas-test
```

Or select the integration test:

```bash
atlas-test --test atlas_integration
```

Desktop window tests are opt-in with `--gui`; the other GUI checks run normally.

The pytest integration test calls the smoke-test command and verifies that the expected output files are produced.

## Local documentation

The documentation is published at [nikolaos-siomos.github.io/ATLAS](https://nikolaos-siomos.github.io/ATLAS/). You can also preview changes locally with the same MkDocs configuration used by the deployment workflow.

The documentation source files are in:

```text
docs/
```

The automatically generated reference pages are in:

```text
docs/generated/
```

These generated pages document the initialization, configuration, settings and
intercomparison INI files. They are generated from the parser schemas and template flavor files, so they should not be edited manually.

### Simple local documentation check

From a fresh clone or extracted copy of the repository, create and activate an environment:

```bash
conda create -n atlas_box python=3.13 pip
conda activate atlas_box
```

Enter the repository root, meaning the folder that contains `mkdocs.yml`:

```bash
cd ATLAS
```

Install ATLAS with the documentation dependencies:

```bash
python -m pip install -e ".[docs]"
```

Start the local documentation server:

```bash
mkdocs serve --config-file mkdocs.yml
```

Open the local URL printed in the terminal. Usually it is:

```text
http://127.0.0.1:8000/
```

Keep the terminal open while reading the documentation. When you edit a documentation file, MkDocs usually rebuilds the page automatically and the browser can be refreshed.

Stop the documentation server with:

```text
Ctrl+C
```

