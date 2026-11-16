# Importing Processing Stages

Import a saved ATLAS checkpoint into Python to inspect or reuse its results
without processing the raw measurements again. Importing returns data; it does
not restore a running `Processor` or automatically resume the processing pipeline.

First [export a processing stage](exported_processing_stages.md). Then use the
examples below to import a whole stage or only the Rayleigh profiles.

## Beginner examples: import saved results

Run these examples in Python or a notebook in the environment where ATLAS is
installed. A previous run must already have exported the stage. No active
`processor` or `caller_info` dictionary is needed.

### Set the output folder and find your stage

Use the installed package's import path:

```python
from atlas_actris.utils.export_processing_stage import (
    list_exported_stages,
    inspect_exported_stages,
    import_processing_stage,
    import_processing_entry,
)

output_folder = "/path/to/ATLAS_output"  # Replace with your actual output folder.
print(list_exported_stages(output_folder=output_folder))
```

Point to the directory **containing `exported/`**, not to `exported/` itself or
to an individual stage directory. On Windows you can use, for example,
`output_folder = r"C:\data\ATLAS_output"`.

Choose a stage name from the printed list, then inspect it:

```python
stage_name = "pol_cal_complete"
inspect_exported_stages(output_folder=output_folder, stage_name=stage_name)
```

Inspection lists the available tests and entries without loading their data.
The following examples assume this stage contains the Rayleigh test, `ray`.
Available contents depend on the measurements and processing that were performed.

### Import the whole stage

```python
stage = import_processing_stage(
    output_folder=output_folder,
    stage_name=stage_name,
)
print(stage.keys())

ray_profiles = stage["ray"]["profile"]
print(ray_profiles)
```

`stage` is a nested dictionary: **test → entry → value**. Importing the whole
stage opens every saved entry. Large numerical arrays are opened lazily; this
creates access to the data without reading all array values into memory yet.
JSON metadata entries are read when imported.

### Import only Rayleigh profiles

If you only need the profiles database for the Rayleigh test, use:

```python
ray_profiles = import_processing_entry(
    output_folder=output_folder,
    stage_name=stage_name,
    qa_test="ray",
    parameter="profile",
)
print(ray_profiles)
```

This returns the same value as `stage["ray"]["profile"]` directly. It opens
only that entry, without importing other tests or separate entries such as mean
profiles and system metadata.
You do **not** add `["ray"]["profile"]` again to this result.

The internal profiles database is called `profile_db`; its exported entry key
is the singular **`profile`**, not `profiles` or `profile_db`. Mean profiles
are a separate entry, `profile_mean`, when present. Numerical profiles are
normally xarray objects, with the coordinates and dimensions saved in the stage.

If you explicitly need the numerical values loaded into memory:

```python
ray_profiles.load()
```

This step is optional; calculations and plotting can also trigger data reads.

## Import other subsets

The examples below use `output_folder` and `stage_name` from the setup above.

### All entries for the Rayleigh test

```python
ray = import_processing_entry(
    output_folder=output_folder,
    stage_name=stage_name,
    qa_test="ray",
)
ray_profiles = ray["profile"]
```

The result is a dictionary of entries for `ray`; the outer test level is removed.

### Selected entries while keeping the nested dictionary

```python
subset = import_processing_stage(
    output_folder=output_folder,
    stage_name=stage_name,
    qa_tests=["ray"],
    parameters=["profile", "profile_mean"],
)
ray_profiles = subset["ray"]["profile"]
```

Only entries present in the export are returned. For multiple tests, pass a list
such as `qa_tests=["ray", "ray_pcb"]`, provided they exist in that stage. Always
use lists for `qa_tests` and `parameters`; use the singular shortcuts for one key:

```python
subset = import_processing_stage(
    output_folder=output_folder,
    stage_name=stage_name,
    qa_test="ray",
    parameter="profile",
)
ray_profiles = subset["ray"]["profile"]
```

Do not combine `qa_test` with `qa_tests`, or `parameter` with `parameters`.
`import_processing_stage()` filters the manifest; an unknown test can give an
empty dictionary, and unknown parameters can leave empty test dictionaries.
Use inspection to check names, or use `import_processing_entry()` for a specific
entry: it raises `KeyError` with available keys if that test/entry is missing.

### Profiles from every test that has them

```python
profiles = import_processing_entry(
    output_folder=output_folder,
    stage_name=stage_name,
    parameter="profile",
)
ray_profiles = profiles["ray"]["profile"]
```

Tests without `profile` are skipped. If no test contains it, the function raises
`KeyError`.

### Return values at a glance

| Function and selection | Result |
| --- | --- |
| `import_processing_stage(...)` | `{test: {entry: value}}`, even with filters |
| `import_processing_entry(..., qa_test="ray", parameter="profile")` | The profiles value directly |
| `import_processing_entry(..., qa_test="ray")` | `{entry: value}` for Rayleigh |
| `import_processing_entry(..., parameter="profile")` | `{test: {"profile": value}}` for matching tests |

Supplying no selectors to `import_processing_entry()` imports the whole stage.
For either importer, omitting `stage_name` opens a terminal selection prompt;
provide it explicitly in scripts and noninteractive sessions.

## Inspect available entries

```python
from atlas_actris.utils.export_processing_stage import list_exported_stage_contents

contents = list_exported_stage_contents(
    output_folder=output_folder,
    stage_name=stage_name,
)
print(contents[stage_name]["ray"])  # Entry names under ray.
```

For a detailed dictionary without printing a tree:

```python
summary = inspect_exported_stages(
    output_folder=output_folder,
    stage_name=stage_name,
    print_tree=False,
)
print(summary[stage_name]["ray"]["profile"])  # Storage kind and relative path.
```

The inspector also accepts `include_kinds=False` or `include_paths=True` to
change its printed listing. It reads only JSON manifests, not the actual Zarr
arrays or metadata entries.

## Export format and older stages

Current imports read JSON manifests and typed metadata, with numerical arrays
stored in Zarr. Older pickle-based exports are rejected; regenerate them with
current ATLAS. See [file formats and existing exports](exported_processing_stages.md#file-format-and-existing-cachesexports)
for migration details and [storage and path checks](exported_processing_stages.md#storage-and-path-checks)
for the directory requirements.

## Troubleshooting Python imports

Use the full module path:

```python
from atlas_actris.utils.export_processing_stage import import_processing_entry
```

`ModuleNotFoundError: No module named 'atlas_actris'` means the Python session
cannot find the ATLAS package; it occurs before stage import or file reading.
Check the interpreter and package discovery inside that same Python session:

```python
import sys
import importlib.util

print(sys.executable)
print(importlib.util.find_spec("atlas_actris"))
```

If discovery returns `None`, select the environment where ATLAS is installed
(for example, the correct interpreter in Spyder). If installation is needed,
activate that environment in a terminal, enter your ATLAS repository root,
and run `python -m pip install -e .`. Restart the Python kernel afterwards.
`python -m pip show atlas-actris` checks installation for the terminal's Python;
compare its interpreter with the one printed inside Spyder or the notebook.

The helpers are defined in the submodule, not re-exported from `utils`. For
example, import the batch exporter with:

```python
from atlas_actris.utils.export_processing_stage import export_processor_stages
```

Use this instead of `from atlas_actris.utils import export_processor_stages`.
