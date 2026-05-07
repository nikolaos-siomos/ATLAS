# ATLAS ACTRIS

Automated Lidar Analysis Software (ATLAS) for lidar data processing, visualization, and quality-assurance workflows.

This repository is currently intended for **developers**. The recommended installation method is a local editable installation from the source tree.

## Installation

Clone the repository and enter the project folder:

```bash
git clone <repository-url>
cd ATLAS
```

Create and activate a Python environment. Python **3.9 or newer** is required.

Using conda or miniforge:

```bash
conda create -n atlas_box python=3.9 pip
conda activate atlas_box
```

Install ATLAS locally in editable mode:

```bash
python -m pip install -e .
```

This installs the package and exposes the `atlas` command in the active environment.

To include the documentation tools as well, install the optional docs dependencies:

```bash
python -m pip install -e ".[docs]"
```

## Minimal usage

Run ATLAS from the command line by providing an initialization file:

```bash
atlas -i /path/to/call_atlas.ini
```

or equivalently:

```bash
atlas --ini_file /path/to/call_atlas.ini
```

To check the available command-line options:

```bash
atlas -h
```

## Local documentation

The full documentation is built with MkDocs and is kept in the `docs/` folder.

Install the documentation dependencies if you have not already done so:

```bash
python -m pip install -e ".[docs]"
```

Start the local documentation server from the repository root:

```bash
mkdocs serve
```

Then open the local URL printed in the terminal, usually:

```text
http://127.0.0.1:8000/
```

To build the static documentation site locally:

```bash
mkdocs build
```

The generated site will be written to the `site/` folder, which should not be committed to Git.

## Development notes

The package uses a `src` layout:

```text
src/
└── atlas_actris/
```

For development in Spyder or another IDE, make sure the IDE uses the same Python environment where ATLAS was installed with:

```bash
python -m pip install -e .
```

The main interactive script can be opened directly in the IDE at:

```text
src/atlas_actris/__call_atlas_interactive__.py
```
