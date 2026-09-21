# ATLAS ACTRIS

[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/)
[![PyPI version](https://img.shields.io/pypi/v/atlas-actris.svg)](https://pypi.org/project/atlas-actris/)
[![Development status](https://img.shields.io/badge/status-development-orange)](https://github.com/nikolaos-siomos/atlas_dev)
[![Tests](https://github.com/nikolaos-siomos/atlas_dev/actions/workflows/tests.yml/badge.svg)](https://github.com/nikolaos-siomos/atlas_dev/actions/workflows/tests.yml)
[![License](https://img.shields.io/badge/license-see%20LICENSE-lightgrey)](LICENSE)
Automated Lidar Analysis Software (ATLAS) for lidar data processing, visualization, and quality-assurance workflows.

ATLAS is intended for users involved in ACTRIS aerosol high-power lidar (AHL) quality-assurance activities. The primary users are the ACTRIS CARS group, who use ATLAS for the evaluation of quality-assurance tests of AHL systems, and National Facility PIs, who use ATLAS to check and monitor the status of their systems.


## Requirements

ATLAS currently supports Python **3.9, 3.10, 3.11, 3.12, and 3.13**.

The package is tested with a bundled smoke-test dataset through `pytest` and GitHub Actions. The test suite checks that the command-line workflow can run to completion and produce the expected plots, reports, and ASCII outputs.

## Installation

Create and activate a Python environment. For example, using conda or miniforge:

```bash
conda create -n atlas_box python=3.13 pip
conda activate atlas_box
```

Clone or just download and extract the repository:

```bash
git clone <repository-url>
```

Enter the top level of the cloned/extracted project folder (the folder name can be different than in the example).
```bash
cd atlas_dev
```

Install ATLAS locally by typing the following command inside the project folder:

```bash
python -m pip install -e .
```

This installs the package and exposes the `atlas` command in the active environment. Using the -e option keeps the source tree editable while working with ATLAS locally.

## Developer installation

For development and testing, install ATLAS with the optional development dependencies:

```bash
python -m pip install -e ".[dev]"
```

This installs the package in editable mode and adds the tools needed to run the test suite, including `pytest`.

To include the documentation tools as well, install the optional documentation dependencies:

```bash
python -m pip install -e ".[docs]"
```

To install both development and documentation tools:

```bash
python -m pip install -e ".[dev,docs]"
```

After installation, it is useful to check that the environment has no dependency conflicts:

```bash
python -m pip check
```

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

ATLAS expects an initialization file that points to the parent folder (input lidar data), configuration file (channel/system configuration), settings file (QA test options), radiosonde or atmospheric input, and output folder where files are exported. Detailed configuration and settings documentation is provided in the MkDocs documentation.

For command-line execution, the `atlas` entry point uses a non-interactive plotting backend by default, which is suitable for terminal runs and automated tests. IDE workflows can still use the plotting backend configured in the active environment.

## Run ATLAS from an IDE
It is also possible to run atlas by opening the following directly in an IDE (e.g. Spyder):

```text
src/atlas_actris/__call_atlas_interactive__.py
```
Please make sure the IDE uses the same Python environment where ATLAS was installed.

## Smoke testing

The repository includes a bundled smoke-test dataset in `testing_pack/`. This test is intended to confirm that the complete ATLAS command-line workflow still works after code, configuration, or dependency changes.

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
- automatically answers the two terminal prompts used during the test run,
- checks that plots, reports, and ASCII files were created.

Generated outputs are written under:

```text
testing_pack/analysis/
```

This folder is ignored by Git and should not be committed.

## Running tests with pytest

Install the development dependencies first:

```bash
python -m pip install -e ".[dev]"
```

Then run the integration test:

```bash
python -m pytest -q tests/test_atlas_integration.py
```

The pytest integration test calls the smoke-test command and verifies that the expected output files are produced.

## GitHub Actions

The repository can run the smoke test automatically on GitHub Actions using the workflow in:

```text
.github/workflows/tests.yml
```

The workflow installs ATLAS, checks dependencies with `pip check`, and runs the smoke test across the supported Python versions.

## Local documentation

The documentation can be checked completely locally. This is currently the recommended way to view the documentation while the repository is private and GitHub Pages is not enabled.

The documentation source files are in:

```text
docs/
```

The automatically generated reference pages are in:

```text
docs/generated/
```

These generated pages document the initialization, configuration, and settings INI files. They are generated from the parser schemas and template flavor files, so they should not be edited manually.

### Simple local documentation check

From a fresh clone or extracted copy of the repository, create and activate an environment:

```bash
conda create -n atlas_box python=3.13 pip
conda activate atlas_box
```

Enter the repository root, meaning the folder that contains `mkdocs.yml`:

```bash
cd atlas_dev
```

Install ATLAS with the documentation dependencies:

```bash
python -m pip install -e ".[docs]"
```

Start the local documentation server:

```bash
mkdocs serve
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

### Strict documentation check

Before sharing documentation changes, run:

```bash
mkdocs build --strict
```

This checks that the documentation can be built as a static website. The `--strict` option treats warnings as errors, so broken links and missing pages are easier to catch.

A successful build ends with a message similar to:

```text
Documentation built in ... seconds
```

The generated static website is written to:

```text
site/
```

The `site/` folder is only a local build output and should not be committed to Git.

### Regenerate the automatic reference pages

If the parser schemas or template flavor files change, regenerate the INI reference pages before checking the documentation:

```bash
atlas-generate-templates --target docs
```

Then run:

```bash
mkdocs build --strict
```

The generated files that should be committed are:

```text
docs/generated/initialization_reference.md
docs/generated/configuration_reference.md
docs/generated/settings_reference.md
```

### GitHub Pages note

GitHub Pages deployment is optional for local testing. If the repository is private and GitHub Pages is not available for the current plan/settings, colleagues can still use `mkdocs serve` to view the full documentation locally.

