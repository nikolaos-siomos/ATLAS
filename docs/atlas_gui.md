# ATLAS GUI

`atlas-gui` opens the schema-driven editor for `call_atlas.ini`, `config_file.ini` and `settings_file.ini`, with a Run tab for launching ATLAS and viewing its output.

## Install and launch

The GUI is installed with ATLAS. In your active environment, install the release
or update a source checkout as described in [installation](installation.md). For
an editable checkout, run from its repository root:

```bash
python -m pip install -e .
atlas-gui
```

The following launch forms are also supported:

```bash
python -m atlas_actris.GUI
python src/atlas_actris/GUI/atlas_schema_gui.py
```

The direct-file form is for development and bypasses the `atlas-gui` launcher
options, including its automatic Linux XIM workaround. Prefer `atlas-gui` for
normal use; `python -m atlas_actris.GUI` uses that same launcher.

The GUI targets ATLAS's supported Python versions, 3.9–3.13, on Windows, Linux and macOS. It uses Tkinter and the Python standard library. Processing requires the normal ATLAS dependencies in the same environment as the GUI.

## Tkinter requirements

Tkinter must be available in the Python environment used to open the GUI. Check it with:

```bash
python -m tkinter
```

This should open a small test window. In a Conda environment, Tkinter support is normally provided by `tk`; install it in the same environment if needed:

```bash
conda install tk
```

For a distribution-managed Python on Debian/Ubuntu, install the corresponding Tkinter package, commonly `python3-tk`. A graphical desktop or configured display is required; a headless terminal alone cannot display Tk windows. Tkinter is not installed as an ATLAS pip dependency.

## INI working folder and startup

After startup, click **Choose INI folder** when ready to select the folder where
the three INI files will be stored. The folder chooser does not open automatically. The GUI
loads `call_atlas.ini` and follows its configuration/settings paths as described
below, without changing existing contents. For a new working folder, it creates
the missing standard INI files as empty drafts. A message lists
which files were loaded and created. New files are empty drafts: populate the
required configuration fields before processing. Empty settings use ATLAS defaults.
Cancelling leaves the editor open and creates nothing.

### Browser starting directory

The GUI remembers the last selected INI working directory. On the next launch,
that directory is used as the starting location when you open a browser. This
also applies to Run-tab browsing and INI/path browsing when no current file or
working folder provides a more relevant location.

This is only a browsing hint: startup leaves the editor empty. No INI file is
selected or loaded, no chooser opens automatically, and there is no reopening
checkbox. Files are loaded only after you make a selection yourself.

The directory is stored in a small per-user JSON file outside the ATLAS
installation, using Python's standard library:

| Platform | Preference file |
| --- | --- |
| Windows | `%LOCALAPPDATA%\atlas_actris\gui_preferences.json` (or the user's `AppData/Local` directory) |
| macOS | `~/Library/Application Support/atlas_actris/gui_preferences.json` |
| Linux | `$XDG_CONFIG_HOME/atlas_actris/gui_preferences.json`, or `~/.config/atlas_actris/gui_preferences.json` when unset |

Only the directory is remembered, not an INI filename, file contents, credentials
or unsaved edits. Missing/corrupt preferences or an unavailable directory use
the normal first-time browser location. Preference read/write failures are
ignored, so read-only installations and unwritable user preference folders do
not stop the GUI; the convenience feature may simply be unavailable.

### Selecting INI files and folders

Folder selection uses an ATLAS list-view browser on all platforms. Visible folders
and files are listed; `.ini` files are labeled **INI file**, so you can see
whether the expected files are present before selecting the folder. Hidden
entries are excluded by default. Enable **Show hidden files** to include dot
files/folders and native Windows hidden/system or macOS hidden entries. The
checkbox applies to that chooser and starts unchecked when a new chooser opens.
Files are
shown for reference and cannot become the selected directory. Double-click a
folder to open it, use **Up**/**Home**, or type a path and press **Go**. **New
folder** creates a folder only when requested. The output-folder chooser also
accepts a typed directory that ATLAS will create later. All file selectors
use the same ATLAS list-view browser, including **Open**, **Save as**,
the Run-tab browse button, and file-valued parameters. Hidden files start hidden
in every chooser, with **Show hidden files** available when needed. File choosers
include a file-name field and file-type filter. Double-clicking a folder opens
it; double-clicking a file selects it. Save dialogs ask before replacing an
existing file and add `.ini` when the file name has no extension. Cancelling
does not write a file.

A startup loading screen appears immediately after the Tk window is created, before fonts/styles, schemas, descriptions and editor controls are prepared. It reports each stage and the current editor filename. A loading window also appears while a folder's
INIs are read, and when opening an individual INI. It displays a message and
progress bar. Working-folder file reads run in a background worker without Tk
calls; editor pages are prepared on the main thread one at a time, with the
current filename shown. Errors close the loading window and are reported.
Loading an INI now refreshes the existing field values instead of destroying
and recreating the field widgets. The schema UI is still constructed once at
startup; this initial GUI work can take longer than parsing the small files.
File-system latency can also affect reads, especially on network storage.

The **Choose INI folder** button also lets you select another folder later. Unsaved
edits must be saved or discarded before changing folders, and changing folders
is blocked while ATLAS is running. If the selected init itself cannot be read or
parsed, an error leaves the current tabs and Run target unchanged. Problems in
associated files use the draft recovery described below.

The Run tab automatically points to the loaded initialization file. Its
**Browse** button also loads that init file and the associated configuration and
settings into all three editor tabs, after checking for unsaved changes. The
selected init file may have a different name from `call_atlas.ini`; the tab
labels identify file roles, and each tab displays the actual loaded path.

Both **Choose INI folder** (when an init file exists) and the Run-tab browser
resolve associated files using the same rules as processing:

| Init parameter | If supplied | If blank or omitted |
| --- | --- | --- |
| `atlas_configuration_file` | Use that file | Use `config_file.ini` beside the selected init |
| `atlas_settings_file` | Use that file | Use `settings_file.ini` beside the selected init |

Relative paths are resolved against the selected init file's directory, including
when `parent_folder` points elsewhere. Absolute paths and `~` are supported;
either path separator is accepted on the current platform. Foreign Windows drive
paths cannot be opened on Linux/macOS; supply the corresponding local mount path.
References are preserved without rewriting the init file.

Both selection routes keep the chosen init open if its referenced configuration
cannot be loaded. A warning explains the problem and offers two ways forward:
correct `atlas_configuration_file` in the init and select it again, or fill in
the empty configuration tab and save it. This fallback is an **unsaved draft**;
it does not overwrite the unavailable/invalid file or silently open a different
existing configuration beside the init. Missing default configurations use the
same draft behavior when an existing init is selected.

Missing settings are created as an empty file at the resolved settings path
when its parent folder exists and is writable, retaining ATLAS defaults. If
settings cannot be resolved, read, parsed or created, a warning and an unsaved
settings draft allow recovery in the same way. Existing files remain unchanged
while loading. The summary identifies full paths and unsaved drafts.

Saving a recovery draft elsewhere updates its corresponding path in the init
editor. Save the init too before processing; the GUI blocks Start while drafts
have no saved file or edited INIs remain unsaved. Saving a draft at the original
referenced path leaves the init unchanged. Processing continues to use the saved
init paths; there is no hidden runtime substitution of files.

After editing and saving the configuration/settings path entries, reselect the
init with the Run-tab browser to reload its associated files.

## Running ATLAS

Select or open `call_atlas.ini`, then click **Start**. The GUI always launches the ATLAS installation bundled with it, using its own Python interpreter. The Run tab displays that environment; there is no separate environment selector or Conda discovery step.

To use a different environment, activate it in your terminal and run `atlas-gui` there. ATLAS, its GUI, and the processing dependencies must be installed in that environment. This works with Conda and virtual environments without invoking a Conda executable from the GUI.

On Linux with an X display, `atlas-gui` sets `XMODIFIERS=@im=none` in its own
process before creating any Tk windows. This bypasses the desktop input-method
connection that can cause slow startup and keyboard disruption with large Tk
interfaces. It applies to the GUI process and its children; it does not change
the launching shell or desktop configuration. Windows and macOS are unaffected.
If you need desktop input-method composition, retain it with:

```bash
atlas-gui --system-input-method
```

This option preserves the existing `XMODIFIERS` value, including a value supplied
explicitly when launching. Basic keyboard layout input remains available with
the default workaround; input-method composition features may require the option.

The child process uses an absolute launcher path inside the GUI package and isolated Python mode. Python modules beside the INI and inherited `PYTHONPATH`/`PYTHONHOME` settings cannot select a different ATLAS. Editable installations (`pip install -e .`) still use the source checkout and see source changes on the next run. Dependencies must be installed in the environment: isolated mode excludes the per-user site-packages directory.

Save any edited INI files before starting. Output appears in the Run tab. The working directory is the selected init file's parent folder. Stop terminates the process tree using the existing Windows or POSIX implementation.

GUI runs have no interactive terminal input. Use **After processing** in the Run tab:

- **Clean cached files** is unchecked by default. Deletes reusable measurement
  caches. The next run may take longer. Recommended on your final run.
- **Default** keeps existing exported stages and does not export any stages.
- **Clean stages** deletes all exported stages in the current output folder.
  A warning appears beside the selection; stages must be exported again if needed.
- **Export stages** exports the parsed `export_stages` selection from
  `call_atlas.ini`. An omitted or blank parameter uses the existing parser default
  (`pol_cal_complete`). Existing exports with the same names are replaced;
  unrelated exported stages are kept.

These actions run after processing and report creation succeed. Export runs
before cache cleanup; an export failure leaves the measurement caches in place.
An empty parsed stage selection is rejected before processing. The controls are
disabled during a run, and choices apply to that run. They are kept for the current
GUI session and reset to their safe defaults when the GUI is reopened. Mouse-wheel
scrolling cannot change the stage selection.

GUI choices override inherited `ATLAS_*_ANSWER` environment settings. Viewer
cleanup stays disabled. Terminal runs retain their existing prompts and
environment overrides.

## Editing and documentation

The editor reads current parser schemas and flavor descriptions from the ATLAS
installation beside the GUI. Restart the GUI after editing those source files.
Validation checks declared field constraints; runtime and cross-field checks
still belong to ATLAS.

**Save** and **Save as** write compact, bare UTF-8 INIs by default:

- Only populated fields are written, regardless of whether a filter hides them.
- Empty parameters and empty sections are omitted. Explicit `False`, `0`, `off`
  and supplied values equal to defaults are retained.
- Leading/trailing whitespace is trimmed when leaving a field, checking entries
  or saving. List tokens are trimmed without removing empty positions; spaces
  inside station names and paths are preserved. Line breaks remain invalid.
- Each editor has an unchecked **Include parameter descriptions on save** box.
  Check it to add current flavor descriptions, examples and version notes to
  the parameters being saved. Uncheck it and save again for bare output.
- Saving regenerates the INI: existing comments, whitespace, BOM and line-ending
  formatting are not preserved verbatim. Populated unknown entries and legacy
  values are retained and flagged by validation. Unsupported unknown multiline
  entries cannot be compacted and are reported rather than silently flattened.
- Files changed outside the GUI are protected from accidental overwrite.

The init **configuration** tab is hidden. The parser still accepts its legacy
`scc_configuration_id` and `scc_compatible_format` values, and supplied values
are retained on Save. This does not re-enable internal SCC downloading.

The GUI has independent, hardcoded recommended-field lists in `schema_loader.py`,
initially copied from the beginner template selections:

| File | Recommended fields |
| --- | --- |
| Init | `slice_measurement`, `exclude_measurement`, `rsonde_station_name`, `rsonde_station_wmo_id`, `cloudnet_station_name` |
| Settings | Quadrant/ring telecover `normalization_region`; polarization calibration `calibration_region` and `rayleigh_region` |

These GUI lists do not change parser categories and do not automatically track
future changes to beginner templates. Configuration-file categories still come
from its parser schema. Mandatory fields always remain visible.

Path fields include a browse button for `output_folder`; its selected directory
does not need to exist yet, because ATLAS creates its output directories.

The original widget/theme font choices are retained, with the larger sizes
(12-point named fonts and 13-point headings/tabs). The GUI no longer forces one
font on all ttk widgets or overrides the native label-frame font. Button labels
use plain text, such as **Browse**, instead of special symbols or glyphs.

Mouse-wheel events over dropdowns scroll the editor page and stop the dropdown's
class binding from changing its value. This covers True/False and other choice
fields, on Windows/macOS wheel events and Linux X11 wheel buttons. Selecting a
value intentionally with the mouse or keyboard still works. No extra font or
other dependency is required.

The GUI-doc button opens the supplied `GUI-doc.pdf` through the operating system's PDF viewer. Linux needs an available `xdg-open` command and associated viewer. The PDF and example INI files are included as package data.
The bundled PDF may show older startup, save and environment controls; this
page and the GUI README describe the current behavior. See the
[older-manuals note](legacy_manuals.md) for known differences.

## Diagnosing rough or uneven text

Font family/size and antialiasing are separate issues. The GUI cannot enable a
font-rendering backend that is missing from the installed Tk library. On Linux,
Tk's Xft build supports antialiased font rendering; conda-forge distributes both
`xft` and `noxft` Tk variants. Desktop/fontconfig settings can also affect a build
that already has Xft.

From the same environment that runs the GUI:

```bash
atlas-gui --font-info
```

This reports the Tk version, windowing system, font renderer if reported by Tk,
scaling, default font and any Conda Tk build identifiers, then exits without
opening the editor or changing font settings. A working desktop display is
required. The font snapshot is taken before the GUI's point-size adjustment.

If the build identifier contains `noxft`, changing fonts in ATLAS will not add
Xft support. Select an Xft-enabled Tk build compatible with that environment.
If the renderer reports `xft`, inspect antialiasing/fontconfig/display settings
instead. Missing renderer information alone is not definitive proof of a no-Xft
build. Do not change Tcl/Tk major versions simply to alter font rendering.

References: [Tk's font-renderer query](https://core.tcl-lang.org/tips/doc/trunk/tip/205.md)
and [Conda Tk builds](https://anaconda.org/channels/conda-forge/packages/tk/files).

## Compatibility fix

Python 3.9's `tkinter.font.nametofont()` does not accept `root=`. The GUI uses `Font(root=..., name="TkDefaultFont", exists=True).copy()` instead, retaining the correct Tcl interpreter and creating a separate font before applying bold styling.

Internal imports work both as an installed package and when the script is run directly. Mouse-wheel bindings use the actual Tk windowing system: X11 button events are registered only on X11.

## Tests

The new `tests/test_gui_compatibility.py` runs through the repository test suite:

```bash
atlas-test --test gui_compatibility
```

It also invokes the supplied GUI regression tests. Normal test runs skip both
native desktop checks, even if a display is available, so they do not open the
editor or folder browser. Explicitly enable those checks with:

```bash
atlas-test --gui --test gui_compatibility
```

This opens and automatically closes test windows; it requires a desktop display.
The desktop checks disable directory memory so they do not read or write real
user preferences.
Tests construct the app with its folder prompt disabled, avoiding both unwanted
file dialogs and mock callbacks being passed to Python 3.9's Tk timer. The other
GUI checks run normally. No ATLAS processing or network downloads are launched
by the GUI compatibility checks.
