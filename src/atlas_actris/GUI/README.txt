ATLAS GUI - 16 November 2026
===========================

Install ATLAS and its GUI together in the active environment:
    python -m pip install -e .
    atlas-gui

Also supported:
    python -m atlas_actris.GUI
    python src/atlas_actris/GUI/atlas_schema_gui.py

The GUI uses Tkinter and Python's standard library. ATLAS processing requires
its usual dependencies. Supported target versions are Python 3.9–3.13 on
Windows, Linux and macOS. Check Tk with: python -m tkinter

Startup and working folder
--------------------------
Click Choose INI folder when ready; there is no automatic startup chooser.
New working folders receive the missing standard INIs as empty drafts. Existing
init files are loaded without modification and their configuration/settings
references are followed. A summary lists loaded/created files and full paths.
Fill required configuration values before processing. Cancel creates nothing.
Unsaved edits are protected; selection is blocked while processing runs.

The Run tab's Browse button loads the chosen initialization file and associated
configuration/settings into all three tabs. atlas_configuration_file and
atlas_settings_file take precedence; relative paths are based on the init file's
directory, not parent_folder. Blank/omitted entries use config_file.ini and
settings_file.ini beside the init. Existing files remain unchanged. Missing
settings are created empty at their resolved path when possible. In either
selection route, unavailable/invalid configuration files produce a warning and
an empty unsaved configuration draft while keeping the selected init open. The
warning suggests correcting the init reference or filling/saving the draft.
Unavailable/invalid settings also become an unsaved draft if they cannot be
loaded/created. No existing linked file is overwritten while loading, and no
different existing sibling is silently selected as a fallback.

Saving a recovery draft elsewhere updates the matching init entry in the GUI;
save the init too before Start. Saving at its original referenced path leaves
the init unchanged. Start is blocked for unsaved files/edits. An invalid init
itself still produces an error and preserves the previous tabs/Run target.
After editing/saving references manually, reselect the init to reload its files.

The last INI working directory is remembered only as the starting location for
browsing on later launches. No folder/file is automatically selected or loaded,
and there is no reopening checkbox. Current file/folder locations take priority
once INIs are loaded. The remembered directory is stored outside the installation
in per-user JSON (LOCALAPPDATA on Windows, Library/Application Support on macOS,
XDG_CONFIG_HOME or ~/.config on Linux, in an atlas_actris directory). Failed
preference reads/writes, corrupt JSON and unavailable paths are ignored; browsers
fall back to their normal first-time locations. No credentials or INI contents
are stored. A read-only ATLAS installation needs no preference write access.

Folder selectors use a consistent list view with Name and Type columns. Folders
and visible files are listed; INIs are labeled INI file. Hidden entries are
excluded initially. Show hidden files includes dot names and native hidden
attributes; the toggle starts unchecked each time the chooser opens. Double-click folders to open,
use Up/Home, or enter a path and Go. Files are for reference, never selected as
folders. New folder is available, and the output selector accepts a typed new
path without creating it. Browse, Open and Save as use this same list-view
browser, with file-name/type controls and overwrite confirmation for saving.

The startup screen appears before fonts/styles, schemas and controls are loaded,
with stage/file messages. A visible message and progress bar also appear during working-folder
loading and individual INI opening. Folder reads use a background worker; Tk
widgets are prepared one editor at a time on the main thread. The filename is
shown during preparation. Worker errors are reported and release the busy UI.
Existing field widgets are reused when files load, rather than recreated. The
schema controls are constructed once at startup. Parsing the INIs themselves
is lightweight; widget creation and filesystem latency can take longer.

Saving and whitespace
---------------------
Save/Save as produce bare, compact UTF-8 files by default. Only populated
parameters and their sections are written. False, 0, off and explicit defaults
are retained. Empty positions in lists are preserved, and hidden populated
fields are saved too. Existing comments/formatting are regenerated on Save.

Each editor has an unchecked Include parameter descriptions on save box.
Enable it to include current flavor descriptions, examples and version notes
only for the populated fields. Disable and save again for bare output.

Whitespace is trimmed on focus loss, validation and Save. List entries are
trimmed individually, preserving positional blanks. Internal spaces in paths
and names remain. Line breaks remain invalid. Unknown supplied entries are
retained and flagged; unsupported unknown multiline entries fail explicitly.
External file changes are checked before overwrite.

Schemas and recommendations
---------------------------
Parser schemas and flavor text are read from the installed ATLAS sources at
startup, without executing those source files. Restart after source changes.
Dynamic stage-name allowed lists remain text inputs with runtime validation.
The init configuration tab is hidden, but supplied legacy values remain in
the model and are saved; internal SCC export remains disabled.

GUI recommended lists are independent snapshots in schema_loader.py:
- Init: slice_measurement, exclude_measurement, rsonde_station_name,
  rsonde_station_wmo_id, cloudnet_station_name.
- Settings: normalization_region for quadrant/ring telecover, calibration_region
  and rayleigh_region for polarization calibration.
Parser categories and beginner template selections are not modified. Config
categories continue to follow its parser. Mandatory fields remain visible.

The output_folder field now has a directory browser. It may point to a folder
that will be created by ATLAS. Other directory/file browse behavior is retained.
The original widget/theme fonts are retained, with the larger 12-point sizes
and 13-point headings/tabs. No global ttk font override or label-frame font
substitution is applied. Browse and other button labels use plain ASCII text.
The mouse wheel scrolls the page even over a dropdown and never cycles its
True/False or other values. Deliberate mouse/keyboard choice remains available.

Running ATLAS
-------------
The Run tab points to the working folder's call_atlas.ini. The GUI always runs
the same ATLAS installation and Python interpreter as itself. To use another
environment, activate it and launch atlas-gui there. There is no Conda discovery
or separate environment selector. The isolated launcher prevents INI-folder
modules and PYTHONPATH/PYTHONHOME from replacing that installation. Editable
installs still pick up source edits. Install processing dependencies inside the
environment, not solely in user site-packages.

Save edits before Start. Output streams into Run. Relative data paths use the
INI parent directory. Stop terminates the process group/tree; Windows process
handling uses taskkill, and POSIX uses process-group signals. No interactive
stdin is provided. Use the Run tab's After processing controls for cache cleanup
and Default / Clean stages / Export stages. Defaults preserve caches and stages;
GUI choices override inherited ATLAS_*_ANSWER settings. Export precedes cache
cleanup; failed exports preserve caches.

On Linux with an X display, atlas-gui bypasses XIM in its own process before
creating Tk windows (XMODIFIERS=@im=none), avoiding input-method interaction that
can cause slow startup and keyboard disruption. Windows/macOS are unaffected.
For desktop input-method composition, launch atlas-gui --system-input-method
to retain the existing configuration. The launching shell is unchanged.

Font rendering diagnostics
--------------------------
For rough/uneven letters, inspect antialiasing separately from font type/size:
    atlas-gui --font-info
This reports renderer/build/scaling facts without changing fonts or OS settings.
It requires a working display. On Linux, a noxft Tk build cannot gain Xft support
through Python font options. Conda supplies xft and noxft builds. Confirm the
renderer/build before changing packages; if Xft is already available, desktop
fontconfig/display settings may need inspection. The diagnostic also reports
Conda tk/libtk build identifiers when metadata is available.

Documentation
-------------
See docs/atlas_gui.md for current behavior. GUI-doc opens the bundled PDF through
the operating system viewer; its older screenshots/instructions may show the
previous startup, saving or environment controls. Current Markdown/this README
supersede those parts of the historical PDF.

Verification
------------
From the project root:
    atlas-test --test gui_compatibility
    python -m unittest discover -s src/atlas_actris/GUI -v

Tests cover schemas/flavors, compact saves, whitespace, positional blanks,
provided defaults, unknown/legacy values, folder load/create/cancel behavior,
recommendations, fonts, output browsing, secure launching, streaming and Stop.
Headless checks use real temporary files plus Tk/Tcl or mocked native controls.
Desktop checks are skipped by default, even with a display, so atlas-test does
not pop up the editor/folder browser. Explicit opt-in (opens test windows):
    atlas-test --gui --test gui_compatibility
A display is required. Test windows close automatically; their startup folder
prompt is disabled. The tests no longer patch a callback with a MagicMock and
pass it to Tk.after(), which requires a function name on Python 3.9.
Native OS dialog appearance/layout still need checking on supported platforms.
