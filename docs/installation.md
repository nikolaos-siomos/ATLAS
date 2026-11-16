# ATLAS Installation

ATLAS is available from PyPI for users and from a source checkout for developers.
The version documented here is **{{ atlas_version }}**. The current package
supports **Python 3.9–3.13** (`>=3.9,<3.14`).

## Create a Python environment

Use an environment dedicated to ATLAS. With Conda or Miniforge, for example:

```bash
conda create -n atlas_box python=3.13 pip
conda activate atlas_box
```

Any supported Python version can be used. With an existing virtual environment,
activate it using that environment's usual command instead.

If you need Conda, the [Miniforge download page](https://conda-forge.org/download/)
provides installers. On Windows, open Miniforge Prompt after installation. Use
`conda info` to check which Conda installation is active.

## Install for normal use

```bash
python -m pip install atlas-actris
```

This installs ATLAS, its GUI launcher, and the runtime dependencies declared in
`pyproject.toml`. Tkinter itself must be available in your Python environment;
see the [GUI requirements](atlas_gui.md#tkinter-requirements).

To update an installed release:

```bash
python -m pip install --upgrade atlas-actris
```

Add `--pre` if you want to allow pre-release versions. The distribution names
`atlas-actris` and `atlas_actris` are equivalent for pip; the Python import name
is `atlas_actris`.

## Install for development

Clone the repository and enter its root, where `pyproject.toml` is located:

```bash
git clone https://github.com/nikolaos-siomos/ATLAS.git
cd ATLAS
python -m pip install -e ".[dev]"
```

If you already have a checkout, activate your environment and run the last
command from that checkout. Editable mode uses your source files on the next
run. Restart an open GUI or Python session to pick up changed imported modules.
Reinstall after changing package metadata, dependencies or console entry points.

The `dev` extra installs pytest. Run all repository tests with:

```bash
atlas-test
```

For selected tests and optional desktop checks, see [the test-suite guide](atlas_test.md).

## Pillow versions and plot safety

ATLAS selects a compatible Pillow release for each Python version:

| Python | Pillow requirement |
| --- | --- |
| 3.9 | `>=11.3,<12` |
| 3.10 | `>=12.2,<13` |
| 3.11–3.13 | `>=13,<14` |

All ATLAS plot readers explicitly accept static PNG only. Renaming a different
image format to `.png` does not make it acceptable. Animated PNGs are rejected
before Pillow opens them, protecting older Pillow versions from an APNG
allocation issue. Valid ATLAS plots, embedded PNG
metadata, report generation and colour reduction continue to work as before.

Pillow 11 and 12 remain necessary for the older Python versions, but do not
contain every security fix available in Pillow 13. Restricting decoding to PNG
avoids the unused PSD/JPEG 2000 decoder paths identified in the security review;
it does not patch those older dependencies or make every image safe. For new
environments, prefer Python 3.11 or newer so current Pillow security releases
can be installed. See Pillow's [Python support table](https://pillow.readthedocs.io/en/stable/installation/python-support.html).

After applying these dependency changes to an editable checkout, reinstall with:

```bash
python -m pip install -e ".[dev]"
```

## Verify the installation

```bash
python -c "import atlas_actris; print(atlas_actris.__file__)"
atlas --help
atlas-gui --help
```

An editable installation should report a package under your checkout's
`src/atlas_actris/`; a normal pip installation reports the environment's installed
package directory. The help commands do not start processing or open the editor.

To start the GUI:

```bash
atlas-gui
```

Or run ATLAS from the terminal using your prepared initialization file:

```bash
atlas -i /path/to/call_atlas.ini
```

The equivalent long option is `--ini_file`. See the [quick guide](quick_guide.md)
for preparing data and INIs, and the [CLI guide](atlas_actris_cli_documentation.md)
for available overrides.

The smoke-test dataset and repository tests are not bundled in a PyPI wheel.
To run the full smoke test, use a source checkout containing `testing_pack/`:

```bash
atlas-smoke-test /path/to/ATLAS/testing_pack
```

The installed command can run from any directory when that dataset is available.
See [smoke testing](atlas_smoke_test.md) for its cleanup behavior.

## Runtime dependencies

Pip installs the numerical, file-reading, download, plotting, report and viewer
dependencies specified by `pyproject.toml`. Some constraints vary by Python
version, including the NumPy, SciPy and interactive-viewer packages. Use that
file as the source of truth rather than a second, manually maintained list of
version constraints here.

## Build or preview the documentation

From the repository root:

```bash
python -m pip install -e ".[docs]"
mkdocs serve --config-file mkdocs.yml
```

Open the address printed by MkDocs. To build without starting a local server:

```bash
mkdocs build --strict --config-file mkdocs.yml
```

The root `mkdocs.yml` controls the sidebar, hooks and navigation. When applying
a documentation update, copy that file as well as `docs/`. Restart the server
if the sidebar still shows older entries. For example, **Importing stages**
should appear immediately below **Exporting stages** under **User guide**.

The documentation hooks substitute the version placeholder from
`src/atlas_actris/version.py` and render the station table from
`src/atlas_actris/molecular/meteorological_database.csv`. Generated INI reference
pages are updated with `atlas-generate-templates`, not by editing their Markdown.

## Spyder development workflow

Install and start Spyder in the same environment if you want to inspect ATLAS
variables while developing:

```bash
python -m pip install spyder
spyder
```

Open `src/atlas_actris/__call_atlas_interactive__.py` and set its run arguments
to `-i /path/to/call_atlas.ini`. The terminal entry point runs that same script.
The `atlas` CLI defaults to the noninteractive Matplotlib `Agg` backend unless
`MPLBACKEND` is already set; running the file directly in an IDE keeps the IDE's
backend setup.

## Older manuals

The PDF/ODT manuals are historical and contain earlier commands and parameter
names. Use the current pages linked above for installation and configuration;
see [the older-manuals note](legacy_manuals.md) before following an archived example.
