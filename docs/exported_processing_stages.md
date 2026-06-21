# Exported Processing Stages

ATLAS can export intermediate processor stages to disk and import them again later. This is useful when you want to inspect, reuse, or pass processed data to another script without repeating the full processing chain.

The export utilities are provided by:

```python
from utils.export_processing_stage import (
    export_processing_stage,
    export_processor_stage,
    ask_export_processor_stage,
    import_processing_stage,
    import_processing_entry,
    list_exported_stages,
    inspect_exported_stages,
    list_exported_stage_contents,
    delete_exported_stage,
    delete_all_exported_stages,
)
```

The exported files are written under:

```text
caller_info["output_folder"]/exported/<stage_name>/
```

Each exported stage keeps the same nested structure used internally by ATLAS:

```text
stage_data
└── qa_test_key
    └── parameter_key
        └── value
```

For example:

```text
pol_cal_complete
├── ray
│   ├── profile
│   ├── profile_error
│   ├── profile_mean
│   ├── range
│   ├── channel_info
│   └── system_info
├── pcb_p45
│   └── ...
└── pcb_m45
    └── ...
```

Large numerical `xarray.DataArray` and `xarray.Dataset` objects are saved as Zarr. Mixed/object metadata, such as `system_info`, `channel_info`, `time_info`, and `pol_cal_info`, are saved with pickle. Zarr-backed arrays are reopened lazily when imported.

## Export a processor stage

The most convenient way to export a stage is to use `export_processor_stage()`. It extracts the stage from the active `processor` and writes it to the `exported` folder inside `caller_info["output_folder"]`.

```python
export_processor_stage(
    processor=processor,
    stage_name="pol_cal_complete",
)
```

This is equivalent to:

```python
stage_data = processor.export_test_from_stage("pol_cal_complete")
caller_info = processor.processing_info["caller_info"]

export_processing_stage(
    stage_data=stage_data,
    stage_name="pol_cal_complete",
    caller_info=caller_info,
)
```

The function returns the path of the exported stage folder.

## Export with a terminal prompt

To ask interactively before saving a stage, use `ask=True`:

```python
export_processor_stage(
    processor=processor,
    stage_name="pol_cal_complete",
    ask=True,
)
```

or use the convenience wrapper:

```python
ask_export_processor_stage(
    processor=processor,
    stage_name="pol_cal_complete",
)
```

The terminal will ask something like:

```text
Export processing stage 'pol_cal_complete' to '<output_folder>/exported/pol_cal_complete'? [y/N]
```

Pressing Enter defaults to **No**. To make Enter default to **Yes**, use:

```python
export_processor_stage(
    processor=processor,
    stage_name="pol_cal_complete",
    ask=True,
    default_answer=True,
)
```

## Estimated exported size

Before export, ATLAS prints an estimated size:

```text
-- Estimated exported size for stage 'pol_cal_complete': 512.34 MB (zarr arrays: 511.90 MB, pickle metadata: 450.12 KB, entries: 123)
```

The estimate is based on `xarray.nbytes` and does not compute lazy Dask arrays.

To disable the printout:

```python
export_processor_stage(
    processor=processor,
    stage_name="pol_cal_complete",
    print_estimated_size=False,
)
```

## Avoid overwriting existing exports

By default, exporting a stage that already exists raises an error. This prevents accidental data loss.

```python
export_processor_stage(
    processor=processor,
    stage_name="preprocessing_complete",
)
```

If the stage already exists, use `overwrite=True` to replace it:

```python
export_processor_stage(
    processor=processor,
    stage_name="preprocessing_complete",
    overwrite=True,
)
```

## Export several stages

You can save multiple stages from the same processing run. Each stage is written to a separate subfolder.

```python
export_processor_stage(
    processor=processor,
    stage_name="screening_complete",
)

export_processor_stage(
    processor=processor,
    stage_name="preprocessing_complete",
)

export_processor_stage(
    processor=processor,
    stage_name="pol_cal_complete",
)
```

The output folder will look like:

```text
<output_folder>/exported/
├── screening_complete/
├── preprocessing_complete/
└── pol_cal_complete/
```

## List exported stages

To see which stages are available:

```python
stages = list_exported_stages(caller_info)
print(stages)
```

Example output:

```python
["pol_cal_complete", "preprocessing_complete", "screening_complete"]
```

## Inspect exported stage contents before import

Before importing data, you can inspect what is stored inside the exported folder. This is useful when you do not remember the available stage names, QA-test keys, or parameter keys.

Use `inspect_exported_stages()` to print a readable tree:

```python
inspect_exported_stages(caller_info)
```

Example output:

```text
-- Exported processing stages:
pol_cal_complete/  (6 qa keys, 120 parameters)
  ray/  (18 parameters)
    - profile  [zarr_dataarray]
    - profile_mean  [zarr_dataarray]
    - profile_error  [zarr_dataarray]
    - range  [zarr_dataarray]
    - channel_info  [pickle]
    - system_info  [pickle]
    - pol_cal_info  [pickle]
  pcb_p45/  (15 parameters)
    - profile  [zarr_dataarray]
    - profile_mean  [zarr_dataarray]
    - pol_cal_ratio_mean  [zarr_dataarray]
    - pol_cal_info  [pickle]
```

The function reads only the small manifest files. It does not open the Zarr arrays and it does not load the pickle files, so it is fast even for large exported stages.

To inspect only one stage:

```python
inspect_exported_stages(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
)
```

To hide the storage type labels:

```python
inspect_exported_stages(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    include_kinds=False,
)
```

To also show the stored relative paths:

```python
inspect_exported_stages(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    include_paths=True,
)
```

You can also return the full inspection summary as a dictionary without printing anything:

```python
summary = inspect_exported_stages(
    caller_info=caller_info,
    print_tree=False,
)
```

The returned object has the structure:

```python
summary[stage_name][qa_test_key][parameter_key] = {
    "kind": "zarr_dataarray",
    "path": "ray/profile/data.zarr",
}
```

For a simpler dictionary containing only the available keys, use `list_exported_stage_contents()`:

```python
contents = list_exported_stage_contents(caller_info)
```

Example result:

```python
{
    "pol_cal_complete": {
        "ray": [
            "profile",
            "profile_mean",
            "profile_error",
            "range",
            "channel_info",
            "system_info",
        ],
        "pcb_p45": [
            "profile",
            "profile_mean",
            "pol_cal_ratio_mean",
            "pol_cal_info",
        ],
    }
}
```

To get the available keys only for one stage:

```python
contents = list_exported_stage_contents(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
)
```

The inspection functions are intended to be used before partial import. For example:

```python
inspect_exported_stages(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
)

ray_profile = import_processing_entry(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    qa_test="ray",
    parameter="profile",
)
```

## Import a full stage

To import a previously exported stage:

```python
data_pack = import_processing_stage(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
)
```

The returned object has the same nested dictionary structure:

```python
profile = data_pack["ray"]["profile"]
channel_info = data_pack["ray"]["channel_info"]
```

Zarr-backed xarray objects are opened lazily, so importing the stage does not automatically load all large arrays into memory.

## Import selected QA-test keys

To import only selected first-level QA-test keys:

```python
data_pack = import_processing_stage(
    caller_info=caller_info,
    stage_name="preprocessing_complete",
    qa_tests=["ray", "tlc_north"],
)
```

You can also use the single-key shortcut:

```python
data_pack = import_processing_stage(
    caller_info=caller_info,
    stage_name="preprocessing_complete",
    qa_test="ray",
)
```

## Import selected parameters

To import only selected parameter keys from all QA-test keys:

```python
profiles = import_processing_stage(
    caller_info=caller_info,
    stage_name="preprocessing_complete",
    parameters=["profile", "range"],
)
```

The single-parameter shortcut is:

```python
profiles = import_processing_stage(
    caller_info=caller_info,
    stage_name="preprocessing_complete",
    parameter="profile",
)
```

## Import one specific entry

For fastest access, use `import_processing_entry()` when you only need one value.

```python
profile = import_processing_entry(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    qa_test="ray",
    parameter="profile",
)
```

This opens only:

```text
<output_folder>/exported/pol_cal_complete/ray/profile/
```

and does not import the rest of the stage.

## Import all parameters for one QA-test key

To import all parameters under one QA-test key:

```python
ray_pack = import_processing_entry(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    qa_test="ray",
)
```

Then access entries as usual:

```python
profile = ray_pack["profile"]
profile_error = ray_pack["profile_error"]
```

## Import one parameter from all QA-test keys

To import one parameter wherever it exists in the stage:

```python
profiles = import_processing_entry(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    parameter="profile",
)
```

The returned object is still nested by QA-test key:

```python
ray_profile = profiles["ray"]["profile"]
pcb_p45_profile = profiles["pcb_p45"]["profile"]
```

## Delete one exported stage

To delete a specific saved stage:

```python
delete_exported_stage(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
)
```

To ask before deleting:

```python
delete_exported_stage(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    ask=True,
)
```

If you do not want an error when the stage is missing:

```python
delete_exported_stage(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    missing_ok=True,
)
```

## Delete all exported stages

To remove all exported stages from the current output folder:

```python
delete_all_exported_stages(
    caller_info=caller_info,
)
```

By default this asks for confirmation, because it deletes the complete folder:

```text
<output_folder>/exported/
```

To delete without asking:

```python
delete_all_exported_stages(
    caller_info=caller_info,
    ask=False,
)
```

## Recommended integration in an ATLAS run

A typical workflow is to export important checkpoints after recipe execution:

```python
run_linear_recipe(
    processor,
    recipe=preprocessing_recipe,
    initial_input="screening_complete",
    checkout_id="preprocessing_complete",
)

ask_export_processor_stage(
    processor=processor,
    stage_name="preprocessing_complete",
)

run_linear_recipe(
    processor,
    recipe=pol_cal_recipe,
    initial_input="preprocessing_complete",
    checkout_id="pol_cal_complete",
)

ask_export_processor_stage(
    processor=processor,
    stage_name="pol_cal_complete",
)
```

Later, in another script or notebook:

```python
from utils.export_processing_stage import import_processing_stage

pol_cal_pack = import_processing_stage(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
)
```

or only read the exact profile needed:

```python
ray_profile = import_processing_entry(
    caller_info=caller_info,
    stage_name="pol_cal_complete",
    qa_test="ray",
    parameter="profile",
)
```

## Notes and limitations

- Zarr-backed xarray objects are opened lazily with `xr.open_zarr()`.
- Object/mixed metadata are preserved with pickle.
- Pickle files should only be loaded from trusted ATLAS outputs.
- Existing stages are not overwritten unless `overwrite=True` is used.
- Partial import is faster than importing a full stage because only the requested manifest entries are opened or loaded.
- Inspection functions read only the stage manifests, not the actual Zarr arrays or pickle payloads.
- The size estimate is approximate and may differ from the final disk usage because of Zarr metadata, compression, and filesystem overhead.
