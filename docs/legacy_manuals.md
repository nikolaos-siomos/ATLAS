# Older PDF and ODT manuals

The current Markdown pages and generated INI references are the maintained
instructions for ATLAS {{ atlas_version }}. The files under `atlas_manual/`
and the GUI's bundled PDF were written for earlier interfaces. They remain
available as historical material, but their installation, launch, parameter
names and GUI screenshots should not be followed as current instructions.

## Earlier user and developer manuals

- [User manual PDF](atlas_manual/ATLAS_user_manual.pdf)
- [Developer manual PDF](atlas_manual/ATLAS_developer_manual.pdf)

The ODT sources and copies under `atlas_manual/source/`, `drafts/` and
`to_check/` also contain older instructions. In particular:

| Older instruction or name | Current instruction |
| --- | --- |
| `main_data_folder` in the init | Use `parent_folder`; when blank, input defaults to the init directory. |
| `rsonde_wmo_number` | Use `rsonde_station_wmo_id` in `[parsing_options]`. |
| `daq_trigger_offset` | Use `zero_bin` in `[Channels]`. |
| Launching `__master__` | Start with the installed `atlas -i ...` or `atlas-gui` command. |

Use the [quick guide](quick_guide.md), [installation guide](installation.md)
and [generated INI references](index.md#ini-file-references) for the current
configuration and commands. Scientific explanations in the older manuals may
still be useful; this note does not replace a review of that material.

## Bundled GUI PDF

The GUI-doc button opens the older PDF shipped with the GUI. It shows the
former automatic init-file chooser and separate Conda-environment selection.
Current behavior uses **Choose INI folder** when requested and always runs
ATLAS from the GUI's own installed environment. Browsing, compact saving,
draft recovery, end-of-run controls and the Linux input-method workaround
have also changed. Follow [the GUI documentation](atlas_gui.md) for those
instructions.

The historical PDF/ODT files have not been rewritten in this update.
