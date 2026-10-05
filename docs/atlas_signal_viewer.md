# `atlas-signal-viewer`

Runs the ATLAS signal-viewer workflow using the same initialization-file interface as the main ATLAS command.

## Usage

```bash
atlas-signal-viewer -i /path/to/call_atlas.ini
```

Equivalent long-option form:

```bash
atlas-signal-viewer --ini_file /path/to/call_atlas.ini
```

The signal viewer reads and processes the configured test data, generates the requested interactive viewer outputs, and stores them under the case output folder in:

```text
analysis/<case-name>/signal_viewer/
```

The command is interactive during normal use. At the end it may ask whether temporary cache files and generated signal-viewer files should be deleted.
