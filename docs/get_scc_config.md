# Generate an ATLAS configuration from SCC

The `get_scc_config` command downloads an SCC Handbook of Instruments (HOI) configuration as CSV and converts it into an ATLAS configuration INI file. Use it to prepare a configuration before running ATLAS.

## Quick start

Run the command in the Python environment where ATLAS is installed:

```bash
get_scc_config -c 665 -o ./configurations
```

Replace `665` with your SCC configuration ID. Both `-c` and `-o` are required. The output folder is created automatically, including missing parent folders. Downloading requires access to the SCC service.

By default, the exporter prompts for recorder channel IDs. Follow the terminal instructions and provide the IDs used by your raw measurement files.

For SCC-compatible input, use `-s` to use SCC channel IDs directly as recorder channel IDs and skip that prompt:

```bash
get_scc_config -c 665 -o ./configurations -s
```

To see the available options:

```bash
get_scc_config --help
```

If the command is missing after updating a local checkout, reinstall from the repository root to register the entry point:

```bash
python -m pip install -e .
```

## Arguments

| Option | Required | Description |
| --- | --- | --- |
| `-c`, `--scc_configuration_id` | Yes | SCC configuration ID to download. |
| `-o`, `--output_folder` | Yes | Destination folder for both the generated INI and downloaded CSV. Quote paths containing spaces. |
| `-s`, `--scc_compatible_format`, `--scc_format` | No | Use SCC channel IDs as recorder channel IDs, without prompting. Disabled by default. |
| `-h`, `--help` | No | Display help and exit. |

The CLI always performs the export. It has no `-e` export-mode option or separate configuration-file path option; the INI filename is generated automatically. The output-folder abbreviation is `-o`, not `-f`.

## Generated files

Both files are written directly into the folder supplied with `-o`:

```text
config_file_665_20261005_222655.ini
scc_config_665_20261005_222655.csv
```

The timestamp format is `YYYYMMDD_HHMMSS`, using the computer's local time. The INI and CSV timestamps are generated separately, so they may differ slightly. No additional `scc_hoi` subfolder is added by this CLI.

The CSV retains the downloaded SCC information. The INI contains the converted ATLAS configuration and should be reviewed before processing measurements.

## Complete the polarization calibration section

The generated INI includes these empty fields:

```ini
[polarization_calibration]
ch_r                      =
ch_t                      =
K                         =
R_to_T_transmission_ratio =
```

These values are not currently retrieved from SCC. Fill them in manually when applicable to your measurements. The `eta` field is intentionally omitted from the generated template.

Use recorder channel IDs for `ch_r` and `ch_t`. Entries at the same position describe one reflected/transmitted channel pair; provide corresponding calibration values in the same order. See the [configuration file reference](generated/configuration_reference.md) for field definitions and validation rules.

The `[System]` section uses simple `name = value` formatting. In `[Channels]` and `[polarization_calibration]`, parameter names and comma-separated values are padded with spaces to align the equals signs and commas within each section. Use a monospaced font when editing to see this alignment. Added padding does not change parsed values.

## Use the generated configuration with ATLAS

Set `atlas_configuration_file` in your ATLAS initialization file to the generated INI path. To retain your manual edits, set `export_hoi_cfg = 0` for subsequent processing runs, so ATLAS uses the prepared file without exporting it again.

The standalone CLI does not change the arguments accepted by `call_atlas`. Its existing export modes and CSV subfolder behavior remain available. Files generated through either workflow receive the empty polarization calibration fields and the updated formatting; automatic INI filename timestamps are specific to the standalone CLI.
