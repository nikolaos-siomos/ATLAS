# ATLAS Installation

This page describes the current developer installation procedure for ATLAS.

ATLAS is currently intended to be installed from a local Git checkout. The recommended installation mode is an editable pip installation:

```bash
python -m pip install -e .
```

Editable mode keeps the installed package connected to the source code in your local working tree. When you edit files under `src/atlas_actris/`, the changes are immediately used by Python without reinstalling the package.

## Current package status

The Python distribution/package name is:

```text
atlas_actris
```

The current project version is:

```text
1.0.0
```

The project requires:

```text
Python >= 3.9
```

The source code uses the `src` package layout:

```text
ATLAS/
├── pyproject.toml
├── README.md
├── docs/
├── tests/
├── tools/
└── src/
    └── atlas_actris/
        ├── __init__.py
        ├── cli.py
        ├── __call_atlas_interactive__.py
        ├── helper_functions/
        ├── processor/
        ├── readers/
        ├── templates/
        └── ...
```

The main installed command-line command is:

```bash
atlas
```

The command expects an initialization file argument:

```bash
atlas -i /path/to/call_atlas.ini
```

or equivalently:

```bash
atlas --ini_file /path/to/call_atlas.ini
```

The command is configured in `pyproject.toml` as an entry point to the package CLI wrapper. The wrapper runs the main interactive ATLAS script while still allowing developers to open and run `src/atlas_actris/__call_atlas_interactive__.py` directly in Spyder.

## Runtime dependencies

The runtime dependencies are installed automatically from `pyproject.toml` when running `python -m pip install -e .`.

Current runtime dependencies include:

- `numpy>=2.0.0`
- `pandas>=2.2.3`
- `xarray>=2023.6.0`
- `netCDF4>=1.6.3`
- `scipy>=1.13.0`
- `matplotlib>=3.8.0`
- `bokeh>=3.4.3`
- `Pillow>=11.1.0`
- `dask>=2024.8.0`
- `requests>=2.0.0`

## Recommended installation with conda or miniforge

Using conda or miniforge is recommended because it makes it easy to create a clean environment with a specific Python version.

Users with an existing conda, miniconda, anaconda, or miniforge installation can skip directly to [Create a new environment](#create-a-new-environment).

### Step 1: Download miniforge

Download the installer for your operating system from:

https://conda-forge.org/download/

### Step 2: Install miniforge

Follow the installation instructions from the miniforge webpage.

On Linux, the installer may ask whether the base environment should be activated automatically when opening a new terminal. Selecting yes is convenient for most users. If you already have another conda installation, be aware that enabling this option may make the new miniforge base environment the default one in your terminal.

On Windows, open the **Miniforge Prompt** from the Start Menu after installation. This opens a terminal with the base environment activated.

### Step 3: Activate the base environment

If your base environment is not already active, activate it manually.

Example on Linux:

```bash
source /home/your_user/miniforge3/bin/activate
```

You should see something like this at the beginning of your terminal prompt:

```text
(base)
```

To verify that the expected conda installation is being used, run:

```bash
conda info
```

or on Linux/macOS:

```bash
which conda
```

The returned path should point to your expected conda or miniforge installation.

## Create a new environment

From the base environment, create a dedicated ATLAS development environment:

```bash
conda create -n atlas_box python=3.9 pip
```

Activate it:

```bash
conda activate atlas_box
```

You can choose a different environment name if you prefer.

## Get the source code

Clone the repository, or use your existing local checkout.

Example:

```bash
git clone <repository-url>
cd ATLAS
```

If you already have the repository locally, go to the project root, where `pyproject.toml` is located:

```bash
cd path/to/ATLAS
```

The project root should contain files and folders similar to:

```text
pyproject.toml
README.md
docs/
src/
tests/
tools/
```

## Install ATLAS in editable mode

From the project root, with the `atlas_box` environment activated, run:

```bash
python -m pip install -e .
```

This installs ATLAS and its runtime dependencies from `pyproject.toml`.

Editable mode is important for development: changes made in the local source files are used immediately without reinstalling the package.

## Verify the installation

Check that the package can be imported:

```bash
python -c "import atlas_actris; print(atlas_actris.__file__)"
```

The printed path should point to your local repository, typically somewhere under:

```text
src/atlas_actris/
```

Check that the command-line entry point is available:

```bash
atlas -h
```

The help output should look similar to:

```text
usage: __call_atlas_interactive__.py [-h] [-i [ini_file]]

arguments

optional arguments:
  -h, --help            show this help message and exit
  -i [ini_file], --ini_file [ini_file]
                        The path to the initialization file of the calling
                        script
```

To run ATLAS, provide the initialization file:

```bash
atlas -i /path/to/call_atlas.ini
```

or:

```bash
atlas --ini_file /path/to/call_atlas.ini
```

Example:

```bash
atlas -i ./docs/templates/call_atlas.ini
```

Adjust the example path to point to an existing initialization file in your local checkout.

> **Note:** The current help output displays the `-i/--ini_file` option as optional because of the argparse configuration. In practice, the main workflow expects an initialization file path to run correctly.

## Install documentation dependencies

The documentation dependencies are optional and are defined in `pyproject.toml` under the `docs` extra.

To install ATLAS together with the documentation tools, run:

```bash
python -m pip install -e ".[docs]"
```

Then serve the documentation locally:

```bash
mkdocs serve
```

## Spyder development workflow

For development, it is useful to run the main script directly in Spyder so that variables remain visible in Spyder's Variable Explorer.

Install Spyder in the same environment if needed:

```bash
python -m pip install spyder
```

Start Spyder from the activated ATLAS environment:

```bash
conda activate atlas_box
spyder
```

Open and run:

```text
src/atlas_actris/__call_atlas_interactive__.py
```

When running from Spyder, configure the script arguments to provide the initialization file, for example:

```text
-i /path/to/call_atlas.ini
```

This keeps the terminal workflow and the Spyder workflow aligned:

```bash
atlas -i /path/to/call_atlas.ini
```

and:

```text
Spyder run arguments: -i /path/to/call_atlas.ini
```

should execute the same underlying script.
