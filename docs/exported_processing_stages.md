# Exporting Processing Stages

An exported stage is a saved checkpoint of ATLAS processing. You can open its
results later in Python without processing the raw measurements again.
For whole-stage and partial-import examples, see
[Importing processing stages](importing_processing_stages.md).

## Export during a normal ATLAS run

Set the stage names in the init file, for example:

```ini
[general_options]
export_stages = pol_cal_complete
```

You can list several stages separated by commas. A blank/omitted `export_stages`
uses the current default, `pol_cal_complete`.

- **GUI:** before clicking **Start**, choose **Export stages** under **After
  processing** in the Run tab. **Default** does not export stages. Exports run
  after successful processing/reporting and replace existing exports with the
  same names. If cache cleanup is also selected, export completes first.
- **Terminal:** `atlas -i ...` asks about deleting existing exports and then
  exporting the configured stages at the end of the run. Answer **No** to
  deleting stages you want to keep, then **Yes** to export. Export replaces
  existing stages with matching names. Shell prompt-answer overrides still apply.

Exports are stored at:

```text
<output_folder>/exported/<stage_name>/
```

Here `output_folder` is the actual ATLAS output directory, which may differ from
both the raw-data folder and the folder containing the INIs. See the
[GUI documentation](atlas_gui.md) for its export and cleanup controls.

## Export from Python code

These examples are for code that already has an active ATLAS `processor`.
Most users can export through the GUI or terminal prompt instead.

```python
from atlas_actris.utils.export_processing_stage import export_processor_stage

output_folder = processor.processing_info["caller_info"]["output_folder"]

saved_path = export_processor_stage(
    processor=processor,
    output_folder=output_folder,
    stage_name="pol_cal_complete",
)
```

`output_folder` is required. By default, an existing export raises
`FileExistsError`. Add `overwrite=True` to replace it, or `ask=True` to ask
before exporting. The prompt defaults to **No**; `default_answer=True` changes
that default. A declined export returns `None`.

The exporter prints an approximate size without computing lazy arrays for that
estimate. Use `print_estimated_size=False` to hide it. Exporting normally computes
and writes entries one at a time; the estimate can differ from the final size
because of compression and metadata.

To export only Rayleigh profiles, add `qa_tests=["ray"]` and
`parameters=["profile"]`. A partial export can later import only the entries
that were saved.

To export several stages with one combined prompt:

```python
from atlas_actris.utils.export_processing_stage import export_processor_stages

saved_paths = export_processor_stages(
    processor=processor,
    output_folder=output_folder,
    stage_names=["preprocessing_complete", "pol_cal_complete"],
    ask=True,
)
```

Stages must be available in that processor. The return value maps stage names
to their saved paths, or is empty if export is declined.

If you already have a nested stage dictionary, export it directly:

```python
from atlas_actris.utils.export_processing_stage import export_processing_stage

stage_data = processor.export_test_from_stage("pol_cal_complete")
export_processing_stage(
    stage_data=stage_data,
    stage_name="pol_cal_complete",
    output_folder=output_folder,
)
```

## Delete saved stages

To delete one exported stage with confirmation:

```python
from atlas_actris.utils.export_processing_stage import (
    delete_exported_stage,
    delete_all_exported_stages,
)

delete_exported_stage(
    output_folder=output_folder,
    stage_name=stage_name,
    ask=True,
)
```

To delete everything inside the output folder's `exported/` directory:

```python
delete_all_exported_stages(output_folder=output_folder)
```

Deletion of all stages asks by default. `ask=False` disables confirmation;
`delete_exported_stage()` defaults to `ask=False`, so supply `ask=True` when
confirmation is desired. For one stage, `missing_ok=True` ignores a missing
export; deleting all already defaults to `missing_ok=True`. These operations
delete stage exports, not the original measurement files or raw caches.

## Storage and path checks

Numerical xarray entries use Zarr; mixed/object metadata use typed JSON, with
Zarr sidecars for nested numerical arrays. Imports reconstruct supported data
types, rather than arbitrary Python classes. Unknown metadata tags and unsafe
Zarr codecs are rejected. Consolidated Zarr metadata is ignored when opening
arrays so it cannot bypass codec validation.

Stage names and QA-test/entry keys must be nonempty single path components.
Empty names, `.`/`..`, separators and drive prefixes are rejected. Manifest paths
must be relative and remain inside their stage. The `exported` root and stage
folders must be real directories, not symbolic links or junctions. Links inside
a stage that escape it are rejected before reads, overwrites or deletion. The
configured output folder itself may be a symbolic link. Batch export checks all
selected targets before overwriting the first stage. These checks cover paths
present when validated; avoid modifying a stage from another process during an
operation.

## File format and existing caches/exports

New stage exports use `manifest.json` with format `atlas_processing_stage_export_v2`.
Numeric xarray entries remain `data.zarr`. Other entries use `data.json`; a
`data.json.arrays/` directory holds any referenced numeric arrays. Manifests
and metadata are committed as JSON only after their referenced arrays are written.
Numeric DataArray names are saved separately as `name.json` and restored on import.

Raw cache paths remain under `<output_folder>/cache/<measurement>/`:
signals/shots remain in `data.zarr`, metadata now use `metadata.json`, and the
cache-wide source manifest is `source_manifest.json` (`atlas_raw_cache_v2`).
The `cache` root must be a real directory; an explicitly configured output
directory may still be a symbolic link.

The first run encountering an old pickle cache reads the original measurements
and writes a fresh cache, without opening old pickle files. This run can take
longer; subsequent runs reuse the JSON cache if the measurements are unchanged.
If an original measurement path is missing, ATLAS reports that it cannot rebuild
and leaves the old cache untouched. Original measurement files are not modified.
The source-manifest comparison still detects additions, removals, size changes,
and modification times.

Old exported stages are listed and can be deleted normally, but importing or
inspecting their `manifest.pkl` raises an explicit legacy-format error. Regenerate
those exports from the original data with the updated ATLAS. There is no automatic
pickle loader or converter. Existing exported stages are not automatically deleted
or rewritten by this migration.

No new dependencies or configuration options are needed. These formats are
cross-platform, but the format change means old ATLAS releases cannot read the
new cache/export files. Numerical array and Zarr codec dependencies should still
be kept updated; data-only serialization does not guarantee arbitrary inputs
have reasonable sizes or are safe for every native decoder.
