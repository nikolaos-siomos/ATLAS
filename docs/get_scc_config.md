# Generate an ATLAS configuration from SCC

The `get_scc_config` command downloads an SCC Handbook of Instruments (HOI) configuration as CSV and converts it into an ATLAS configuration INI file. Use it to prepare a configuration before running ATLAS.

## Quick start

Run the command in the Python environment where ATLAS is installed:

```bash
get_scc_config -c 665 -o ./configurations
```

Replace `665` with your SCC configuration ID. Both `-c` and `-o` are required. The output folder is created automatically, including missing parent folders. Downloading requires access to the SCC service.

## SCC login

Whenever a configuration download is needed, the command prompts in the terminal:

```text
SCC username:
SCC password (hidden):
```

The password is not displayed. Neither value is saved or reused across invocations. There are no username/password command-line arguments or INI settings, and no OS password-manager application or new Python dependency is required.

Use the credentials accepted by the SCC export service; these may differ from your personal SCC website login. If authentication is rejected, the command exits with an explanatory error. Run it again to retry. It does not repeatedly prompt or silently use a built-in account.

Run this CLI in an interactive terminal. If hidden password entry is unavailable, the command stops rather than falling back to visible password input. Ctrl+C cancels the export. Existing local configuration reuse through the internal helper does not trigger a prompt; the standalone CLI itself requests a fresh export on every invocation.

Downloads use certificate-verified HTTPS to the fixed SCC export endpoint, with
a 10-second connection timeout and a 60-second read timeout. ATLAS checks
redirects itself: up to three redirects are allowed only to the same
configuration export on the same HTTPS host, including addition/removal of a
trailing slash. A different host, non-default port, HTTP destination, configuration
ID, login page or URL containing user information/query/fragment is refused before
credentials are forwarded. Failed or interrupted transfers remove their temporary
file without replacing an existing completed CSV. These timeouts apply to
connecting and individual reads, not the total duration of a download.

A redirect error alone does not establish that a password is wrong. A trailing-
slash redirect is handled automatically. If SCC instead redirects to a login
page, the export may require a web/admin login or additional permission beyond
the community HTTP login; this downloader does not automate that separate login.
An HTTP downgrade requires a working HTTPS export endpoint. Error messages
identify these cases without displaying redirect URLs, cookies or tokens.

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
| `-c`, `--scc_configuration_id` | Yes | Positive integer SCC configuration ID to download. |
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

Set `atlas_configuration_file` in your ATLAS initialization file to the generated INI path, or place it beside the initialization file as `config_file.ini`. Internal SCC export during `atlas` processing is disabled: `export_hoi_cfg` is not a user-facing setting and is fixed internally to `"0"`. Entries retained in older initialization files are ignored with a warning. Run `get_scc_config` separately whenever a new export is needed.

The standalone CLI does not change the arguments accepted by `atlas`. Internal export helpers are retained for future development, but are not enabled by initialization-file values. The existing CSV output-location behavior and generated configuration formatting are retained.
