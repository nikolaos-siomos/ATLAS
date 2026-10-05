# `atlas-intercomparison`

Runs the ATLAS exported-stage intercomparison workflow. The command expects one intercomparison initialization file.

## Usage

```bash
atlas-intercomparison -i /path/to/intercomparison.ini
```

Equivalent long-option form:

```bash
atlas-intercomparison --ini_file /path/to/intercomparison.ini
```

The command dispatches through:

```text
src/atlas_actris/cli_intercomparison.py
```

and runs the main interactive script:

```text
src/atlas_actris/__intercomparison_interactive__.py
```

The same main script can be opened and run directly in Spyder. Configure the Spyder run arguments with:

```text
-i /path/to/intercomparison.ini
```

The parsed initialization dictionary remains available in the Spyder namespace as `intercomparison_info`.

## Intercomparison templates

The generated full and bare templates are:

```text
templates/intercomparison.ini
templates/intercomparison_bare.ini
```

The full template contains parameter descriptions, defaults, allowed values, and examples, but all assignments are intentionally empty. The bare template contains the same section and parameter structure without per-parameter flavor comments.

The generated reference page is:

```text
docs/generated/intercomparison_reference.md
```

Generate these templates using [atlas-generate-templates](atlas_generate_templates.md).
