# ATLAS Quick Guide
This is a quick guide that introduces the users to ATLAS and helps them to complete their first run.

## Data preparation

In ATLAS, raw data files must be placed in specific folders whose names indicate the measurement type. The folder containing all QA-test subfolders is called the **parent folder**.

For your first run, prepare these subfolders inside the parent folder:

| Folder | Contents |
| --- | --- |
| `ray/` | Normal Rayleigh measurement files. |
| `drk/` | Common or long dark measurement files. |

After getting used to ATLAS, you can add any of the following supported subfolders. All paths below are relative to the parent folder.

| Folder | Contents |
| --- | --- |
| `tlc/north/`, `tlc/east/`, `tlc/south/`, `tlc/west/` | Quadrant telecover files, placed in the corresponding sector subfolder. |
| `tlc_rin/inner/`, `tlc_rin/outer/` | Ring telecover files, placed in the corresponding ring subfolder. |
| `pcb/p45/`, `pcb/m45/` | Δ90 polarization-calibration files, with the calibrator at +45° or −45°, respectively. |
| `ray_pcb/` | Rayleigh measurement files acquired in calibration mode. |
| `trg/` | Zero-bin test files. |
| `pcb_aux/p45/`, `pcb_aux/m45/` | Auxiliary Δ90 polarization-calibration files used to characterize the neutral-density (ND) filter used for polarization calibration. |
| `drk_xxx/` | Dedicated dark measurements. Replace `xxx` with a valid QA-test suffix, such as `ray`, `pcb`, or `tlc`, to use those dark measurements exclusively for that test. |

## ATLAS installation

ATLAS must be installed in a Python environment. This guide covers Python versions 3.9–3.13.

1. Create a new environment with a supported Python version. This guide calls the environment `atlas_box`.
2. Activate the environment.
3. Install ATLAS in `atlas_box`:

   ```bash
   pip install atlas-actris
   ```

4. Check the installation:

   ```bash
   atlas-smoke-test
   ```

The version documented here is **{{ atlas_version }}**. Add `--pre` to the
installation command if you want to include pre-release versions.

The `atlas-smoke-test` command can be executed from any directory.

## Start the graphical interface (GUI)

The GUI is included with ATLAS. It lets you edit the INI files and run ATLAS
without typing the processing command each time.

1. Activate the environment where you installed ATLAS (called `atlas_box` in
   this guide):

   ```bash
   conda activate atlas_box
   ```

2. Start the GUI from any directory:

   ```bash
   atlas-gui
   ```

3. Click **Choose INI folder** and select the folder for your INI files. The GUI
   loads available files and helps you prepare missing ones.
4. Review the three INI tabs. Fill in the required configuration values and
   save your files. Empty settings use ATLAS defaults.
5. Open the **Run** tab and click **Start**. Progress and results appear there.

On later launches, browsing starts in the last INI working directory when it is
available. You still choose the folder or file yourself; nothing loads automatically.

The GUI uses the same Python environment in which you started it. For details
about browsing, saving, recovering missing files and running ATLAS, see the
[GUI documentation](atlas_gui.md).

## Configuring ATLAS

Parameters are supplied through the following INI files:

| File | Purpose |
| --- | --- |
| `call_atlas.ini` | Initialization file. It must be provided, but its settings can be left empty. |
| `config_file.ini` | Configuration file. Mandatory. |
| `settings_file.ini` | Settings file. Optional; skip it on your first run. |

### Using the INI file templates

Generate empty templates by running the following command in `atlas_box`:

```bash
atlas-generate-templates -o <path_to_templates_folder>
```

The `atlas-generate-templates` command can be executed from any directory.

The following template types are generated:

| Filename suffix | Contents |
| --- | --- |
| `bare` | All possible parameters, without explanatory text. |
| `beginner` | The most commonly used parameters, with explanatory text. |
| No specific suffix | All possible parameters, with explanatory text. |

For your first run:

1. Copy `call_atlas_beginner.ini` into the parent folder.
2. Rename it to `call_atlas.ini`.
3. Leave its settings empty.

### Preparing the configuration file

The configuration file can be prepared manually or partially generated from the SCC Handbook of Instruments (HoI), if a corresponding SCC configuration is available.

For your first run:

1. Generate the configuration file from the SCC HoI using the instructions below.
2. Review it and fill in any missing parameters.
3. Copy it into the parent folder.
4. Rename it to `config_file.ini`.

#### Licel or PollyXT raw data files

Run the following command in `atlas_box`:

```bash
get_scc_config -c <SCC_configuration_ID> -o <path_to_output_folder>
```

The `get_scc_config` command can be executed from any directory.

You will be prompted to provide the `recorder_channel_id` for each available channel in the SCC configuration. This corresponds to the channel ID in the raw-file header:

- **Licel files:** `BT0`, `BC0`, `BT1`, and so on.
- **PollyXT files:** `1`, `2`, `3`, and so on, using ascending channel indices starting at 1.

#### SCC raw NetCDF files

Run the following command in `atlas_box`:

```bash
get_scc_config -c <SCC_configuration_ID> -o <path_to_output_folder> -s
```

The configuration file is generated using the SCC channel IDs for `recorder_channel_id`.

#### Review the generated configuration

**The automatically generated configuration must be reviewed.** As of ATLAS **{{ atlas_version }}**, automatic generation is incomplete because:

- The link between a normal configuration and a calibration configuration is not fully exposed by the SCC HoI. Fill in the `polarization_calibration` section of `config_file.ini` manually.
- Mistakes in an SCC configuration can go unnoticed.
- The `telescope_type` assignment, such as `n` for near range and `f` for far range, is based on the full-overlap distance supplied for each telescope. That distance is not always correct.
- The `channel_subtype` values `r` (reflected) and `t` (transmitted) are not directly available in the SCC. They must be inferred from the optional `Raw_signal_type` parameter used for products. A common mistake is to provide inconsistent `Raw_signal_type` values for different products, such as ELDA and HIRELPP. In that case, `get_scc_config` prints a warning; review `channel_subtype` for the affected channels.

## Running ATLAS

### First-time run

Before running ATLAS, check that the parent folder contains:

| Path relative to the parent folder | Contents |
| --- | --- |
| `ray/` | Rayleigh measurement files. |
| `drk/` | Dark measurement files. |
| `config_file.ini` | Reviewed configuration file. |
| `call_atlas.ini` | Initialization file. |

In the `atlas_box` environment, run:

```bash
atlas -i <path_to_call_atlas.ini>
```

The `atlas` command can be executed from any directory.

### Subsequent runs

After gaining some experience, you can:

- Add more QA-test folders to the parent folder, such as `tlc/` or `pcb/`.
- Use the settings file by copying `settings_file_beginner.ini` from the templates folder into the parent folder and renaming it to `settings_file.ini`.
- Consult the manual and further configure `call_atlas.ini` and `settings_file.ini` as needed.
- Supply a radiosonde file manually and configure its parsing options if necessary.

For more detail, see the [initialization reference](generated/initialization_reference.md), [settings reference](generated/settings_reference.md), and [meteorological profile selection guide](radiosonde_selection.md).

## Import saved results in Python

An exported stage is a saved processing checkpoint. Before a GUI run, select
**Export stages** in the Run tab to save the stages configured in the init.
The usual default is `pol_cal_complete`. Exports are stored inside the output
folder's `exported/` directory.

Later, open Python or a notebook in your ATLAS environment. Replace the example
path below with the **output folder containing `exported/`**:

```python
from atlas_actris.utils.export_processing_stage import (
    import_processing_stage,
    import_processing_entry,
)

output_folder = "/path/to/ATLAS_output"
stage_name = "pol_cal_complete"
```

To import the whole stage:

```python
stage = import_processing_stage(output_folder=output_folder, stage_name=stage_name)
print(stage.keys())  # Available tests, such as "ray".
```

Alternatively, import only the profiles for the Rayleigh test:

```python
ray_profiles = import_processing_entry(
    output_folder=output_folder,
    stage_name=stage_name,
    qa_test="ray",
    parameter="profile",
)
print(ray_profiles)
```

The profiles database is saved under the singular key **`profile`**. The second
example returns that value directly, without importing other tests or entries.
It requires an export containing `ray` and `profile`; stage contents depend on
the measurements and processing performed. Large numerical arrays are opened
lazily, so importing does not immediately read all their values into memory.

See [importing stages](importing_processing_stages.md) for how to list available
stages, inspect their contents, and import other subsets. For saving checkpoints,
see [exporting stages](exported_processing_stages.md).
