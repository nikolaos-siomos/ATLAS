# Radiosonde selection and downloading

ATLAS ACTRIS uses radiosonde profiles for processing steps that require molecular atmospheric profiles. In the interactive workflow, radiosonde handling happens after lidar data and metadata have been loaded and normalized, and before the radiosonde files are read into the meteorological dataset.

```text
parse initialization file
read lidar files
update metadata
find or select radiosonde
load radiosonde
continue processing
```

In the main interactive script this corresponds to:

```python
metadata = find_radiosonde(caller_info, metadata)
meteo, metadata = load_radiosonde(caller_info, metadata)
```

The first line decides which radiosonde file should be used and stores the decision in `metadata["radiosonde_info"]`. The second line reads the selected radiosonde file.

## Configuration parameters

Radiosonde behavior is controlled from the initialization file through `caller_info`.

| Parameter | Purpose |
| --- | --- |
| `radiosonde_file` | Optional path to a specific radiosonde file. If provided, this file is used directly. |
| `radiosonde_folder` | Folder used for automatic downloads and automatic radiosonde file selection. Ignored when `radiosonde_file` is provided. |
| `rsonde_station_wmo_id` | Optional WMO station identifier used to download a radiosonde from Wyoming. |
| `cloudnet_station_name` | Optional Cloudnet station name used to download ECMWF/Cloudnet meteorological profiles. |
| `rsonde_*` ASCII settings | Column, delimiter, header, footer, and unit settings used later when loading manually provided custom ASCII radiosondes. |

## Manual radiosonde mode

If `radiosonde_file` is provided, ATLAS uses the manually supplied file and skips all automatic discovery logic.

Manual mode does the following:

1. Validates that `radiosonde_file` points to an existing file.
2. Uses the provided file directly.
3. Skips Wyoming download.
4. Skips Cloudnet download.
5. Skips folder-based radiosonde selection.
6. Stores `radiosonde_source = "manual"` in the radiosonde metadata.
7. Detects the loadable `radiosonde_format` by trying the available radiosonde readers.

If both `radiosonde_file` and `radiosonde_folder` are provided, `radiosonde_file` takes precedence. ATLAS prints a warning that `radiosonde_folder` will be ignored.

Example initialization entry:

```ini
[System]
radiosonde_file = /data/campaign/radiosondes/20240428_1200_wyoming_bufr_16716.txt
```

Manual mode is useful when:

- the best radiosonde file is known in advance;
- the file is outside the normal `radiosonde_folder`;
- the filename does not follow the automatic station/download convention;
- automatic download should be avoided;
- reproducibility requires a fixed radiosonde file.

Manual mode does **not** reject the file based on a time-window test. If the filename contains a supported timestamp pattern, that time is stored as `radiosonde_time`. If the timestamp cannot be parsed, ATLAS stores the lidar measurement midpoint as a metadata placeholder and prints a warning.

Manual format detection is content-based. ATLAS creates a temporary `radiosonde_info` object, sets one candidate `radiosonde_format` at a time, and calls the corresponding reader. The first reader that successfully returns a valid meteorological profile determines the stored `radiosonde_format`. The temporary profile is discarded; the selected file is read again later by `load_radiosonde()`, so the rest of the processing chain keeps the usual behavior.

Typical probing order is:

```text
text-like files: custom_ascii -> wyoming -> ecmwf -> scc
NetCDF files:    ecmwf -> scc -> wyoming -> custom_ascii
```

If the filename clearly matches one of the standard ATLAS patterns, that format is tried first.

## Automatic radiosonde mode

If `radiosonde_file` is not provided, ATLAS uses automatic mode. Automatic mode works per supported QA test. Currently, radiosonde selection is applied to:

- `ray`
- `ray_pcb`

For each supported test, ATLAS computes the midpoint of the measurement time range. This midpoint is the target time used for downloads and file selection.

Automatic mode follows this sequence:

1. Compute the measurement midpoint.
2. Try to download a Wyoming radiosonde if `rsonde_station_wmo_id` is provided.
3. Try to download a Cloudnet/ECMWF profile if `cloudnet_station_name` is provided.
4. Search `radiosonde_folder` for supported radiosonde files.
5. Reject supported-looking filenames whose station identifier does not match the expected identifier for that format.
6. Prefer the closest Wyoming file within the priority time window.
7. Otherwise select the closest supported file within the fallback time window.
8. Store the selected file and metadata in `metadata["radiosonde_info"]`.

## Download sources

### Wyoming

Wyoming download is attempted when `rsonde_station_wmo_id` is set. The WMO ID is passed to the Wyoming downloader together with the midpoint date/time. If the download succeeds, the file is saved in `radiosonde_folder` and then participates in normal file selection.

Example expected downloaded filename:

```text
20240428_1200_wyoming_bufr_16716.txt
```

The selector treats this as:

```text
radiosonde_format = "wyoming"
radiosonde_source = "auto"
```

#### Wyoming reuse and failure cache

Wyoming downloads can take longer the first time they are attempted for a station and measurement time, because ATLAS may need to contact the Wyoming server for several candidate sounding times and data sources. The downloader searches nominal Wyoming sounding times within ±12 hours of the requested measurement midpoint and tries the supported Wyoming sources in priority order.

To make later runs faster, ATLAS reuses local information in two ways:

1. If a matching Wyoming file already exists in `radiosonde_folder`, it is reused directly and no new download is attempted.
2. If a previous Wyoming request failed, the failed station/time/source attempt is stored in a failure cache file named `.wyoming_download_failures.json` inside `radiosonde_folder`. On later runs, cached failures are skipped instead of contacting the Wyoming server again for the same unavailable candidate.

This means the first run can be slower, especially when no Wyoming sounding is available, while later runs are usually faster because ATLAS already knows which candidate requests failed.

#### Clearing the Wyoming failure cache

To force ATLAS to retry previously failed Wyoming requests, delete the failure cache file from the radiosonde folder:

```bash
rm /path/to/radiosonde_folder/.wyoming_download_failures.json
```

This only clears the cached information about failed Wyoming attempts. It does **not** delete downloaded radiosonde files.

To force a successful Wyoming sounding to be downloaded again, also delete the corresponding downloaded file, for example:

```bash
rm /path/to/radiosonde_folder/20240428_1200_wyoming_bufr_16716.txt
```

Use this second command only if you intentionally want to remove the local radiosonde data file.

### Cloudnet

Cloudnet download is attempted when `cloudnet_station_name` is set. The station name and midpoint date/time are passed to the Cloudnet downloader. If a file is returned, it is saved in `radiosonde_folder` and then participates in normal file selection.

Example filename:

```text
20240428_0300_ecmwf_antikythera.nc
```

The selector treats this as:

```text
radiosonde_format = "ecmwf"
radiosonde_source = "auto"
```

Cloudnet/ECMWF files are also cached locally. If the expected NetCDF file already exists in `radiosonde_folder`, ATLAS reuses it. Otherwise, the downloader tries to obtain the daily ECMWF model file from Cloudnet and saves it under the expected local filename. To force a Cloudnet/ECMWF file to be downloaded again, delete the corresponding local `.nc` file from `radiosonde_folder`.

Downloads do not automatically guarantee selection. Downloaded files are still evaluated by the filename, identifier, timestamp, and priority rules below.

## Supported filename patterns

Automatic selection only considers files whose names match a supported pattern. Other files and subfolders are ignored.

| Format | Pattern | Example |
| --- | --- | --- |
| Wyoming legacy | `YYYYMMDD_HHMM_wyoming_<wmo_id>.dat` | `20240428_1200_wyoming_16716.dat` |
| Wyoming BUFR text | `YYYYMMDD_HHMM_wyoming_bufr_<wmo_id>.txt` | `20240428_1200_wyoming_bufr_16716.txt` |
| ECMWF / Cloudnet | `YYYYMMDD_HHMM_ecmwf_<site>.nc` | `20240428_0300_ecmwf_antikythera.nc` |
| SCC | `rs_YYYYMMDD<site>HH.nc` | `rs_20240428aky12.nc` |
| Custom ASCII | `YYYYMMDD_HHMM_custom_ascii_<identifier>.<txt|dat|csv|asc>` | `20240428_1200_custom_ascii_athens.txt` |

The timestamp encoded in the filename is used to compare each radiosonde file against the measurement midpoint.

### Identifier matching in automatic mode

Automatic selection also checks the station identifier encoded in the filename. The expected identifiers are collected per format, using the same keys as `radiosonde_format`:

```python
identifiers = {
    "wyoming": rsonde_station_wmo_id,
    "ecmwf": cloudnet_station_name,
    "scc": station_id,
    "custom_ascii": rsonde_station_name,
}
```

A filename can match the regular expression and still be rejected if its parsed identifier does not match the expected identifier for that format. Identifier comparison is normalized before comparison, including case normalization, so a downloaded file such as `20251203_1159_ecmwf_barcelona.nc` can match a configuration value such as `Barcelona`. If the expected identifier for a format is missing or empty, files of that format are not selected automatically.


## Parsing selected radiosonde files

After selection, ATLAS reads the selected file with the reader that corresponds to `radiosonde_format`.

| `radiosonde_format` | Reader behavior |
| --- | --- |
| `wyoming` | Reads the ATLAS-downloaded Wyoming CSV text file using the fixed Wyoming column layout. |
| `custom_ascii` | Reads a manually provided or automatically selected custom ASCII file using the `rsonde_*` options from the initialization file. |
| `ecmwf` | Reads a Cloudnet/ECMWF NetCDF file. This follows the Cloudnet model-file structure and is normally produced by the Cloudnet downloader. |
| `scc` | Reads an SCC-style radiosonde NetCDF file using the standard SCC variable names. |

### Wyoming ASCII format

Wyoming files downloaded by ATLAS are saved as comma-separated text files, for example:

```text
20240428_1200_wyoming_bufr_16716.txt
```

The Wyoming reader expects this file to contain one header line followed by comma-separated numeric data. The reader uses fixed columns from the Wyoming CSV response:

| Quantity | Column used by ATLAS | Expected unit |
| --- | ---: | --- |
| Height above sea level | 5 | m asl |
| Pressure | 4 | hPa |
| Temperature | 6 | °C |
| Relative humidity | 9 | % |

The reader converts these internally to the units required by the molecular calculations:

- pressure: `hPa` to `Pa`
- temperature: `C` to `K`
- relative humidity: `percent` to `fraction`
- height remains `m_asl`

This fixed layout is used only for files classified as `radiosonde_format = "wyoming"`.

### Custom ASCII format

Custom ASCII files are used when `radiosonde_format = "custom_ascii"`. Their parsing is controlled by the initialization-file parameters below.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `rsonde_skip_header` | `1` | Number of lines skipped at the beginning of the file. |
| `rsonde_skip_footer` | `0` | Number of lines skipped at the end of the file. |
| `rsonde_delimiter` | `S` | `S` means whitespace-separated; `C` means comma-separated. |
| `rsonde_column_index` | `2, 1, 3, 5` | 1-based column indices for height, pressure, temperature, and optionally relative humidity. |
| `rsonde_column_units` | `m_asl, hPa, C, percent` | Units corresponding to the selected columns. |

The first three selected columns are mandatory and must represent:

```text
height, pressure, temperature
```

The fourth selected column is optional and represents relative humidity. If only three columns are provided, relative humidity is filled with missing values.

The default custom ASCII interpretation is therefore:

```text
height      -> column 2, m_asl
pressure    -> column 1, hPa
temperature -> column 3, C
humidity    -> column 5, percent
```

Allowed height units are:

```text
m_asl, m_agl, km_asl, km_agl
```

Allowed pressure units are:

```text
Pa, hPa, atm
```

Allowed temperature units are:

```text
K, C, Cx10
```

Allowed humidity units are:

```text
percent, fraction
```

If the height unit is `m_agl` or `km_agl`, `rsonde_station_altitude` must also be provided so ATLAS can convert the profile to altitude above sea level.

Example custom ASCII configuration:

```ini
[System]
rsonde_skip_header = 1
rsonde_skip_footer = 0
rsonde_delimiter = S
rsonde_column_index = 2, 1, 3, 5
rsonde_column_units = m_asl, hPa, C, percent
```

### ECMWF and SCC expected units

Cloudnet/ECMWF files are read from NetCDF files produced or reused by the Cloudnet downloader. ATLAS uses the model variables for height, pressure, temperature, and relative humidity directly from the Cloudnet file and converts the model height coordinate to `height_asl` by adding the lidar station altitude.

SCC radiosonde NetCDF files are expected to contain the standard SCC-style variables:

```text
Altitude
Pressure
Temperature
RelativeHumidity
```

For SCC files, ATLAS expects:

| Quantity | Expected unit |
| --- | --- |
| `Altitude` | m asl |
| `Pressure` | hPa |
| `Temperature` | K |
| `RelativeHumidity` | % |

### Height-coordinate cleanup

Before the meteorological profile is used for molecular calculations, ATLAS cleans the `height_asl` coordinate. This is done for all radiosonde formats after unit conversion and before number density and saturation vapour pressure are added.

The cleanup step:

1. sorts the profile by `height_asl`;
2. removes invalid height levels such as `NaN` or infinite values;
3. removes duplicate height levels, keeping the first occurrence;
4. raises an error if no valid height levels remain.

This cleanup prevents interpolation errors during molecular calculations. In particular, `xarray.interp()` requires the interpolation coordinate to be uniquely valued. Some radiosonde files contain repeated altitude levels, especially after rounding or unit conversion, so removing duplicate `height_asl` values avoids errors such as:

```text
InvalidIndexError: Reindexing only valid with uniquely valued Index objects
```

After cleanup, ATLAS adds additional atmospheric parameters used later in the processing chain:

```text
N    -> molecular number density
e_s  -> saturation vapour pressure
```

## Selection priority

Automatic mode uses two time windows:

| Rule | Default window | Description |
| --- | ---: | --- |
| Wyoming priority window | ±3 hours | If at least one Wyoming file is within this window, choose the closest Wyoming file. |
| General fallback window | ±18 hours | If no Wyoming file passes the priority rule, choose the closest supported file of any format within this window. |

This means a Wyoming file can be selected even if another non-Wyoming file is slightly closer, as long as the Wyoming file is within the priority window.

Example:

```text
Measurement midpoint: 2024-04-28 11:00 UTC
Wyoming file:         2024-04-28 12:00 UTC  -> delta = +1 h
ECMWF file:           2024-04-28 10:30 UTC  -> delta = -0.5 h
```

The Wyoming file is selected because it is within the ±3 h priority window.

If no Wyoming file is within ±3 h, the selector falls back to the closest supported file within ±18 h.

## Metadata written by selection

When a radiosonde file is selected successfully, ATLAS stores radiosonde metadata for the corresponding QA test.

Typical automatic-mode metadata:

```text
radiosonde_file     = /path/to/20240428_1200_wyoming_bufr_16716.txt
radiosonde_format   = wyoming
radiosonde_source   = auto
radiosonde_time     = 2024-04-28T12:00:00
measurement_time    = <measurement midpoint>
radiosonde_status   = 0
```

Typical manual-mode metadata:

```text
radiosonde_file     = /path/to/manual_profile.txt
radiosonde_format   = custom_ascii
radiosonde_source   = manual
radiosonde_time     = <parsed file time or measurement midpoint>
measurement_time    = <measurement midpoint>
radiosonde_status   = 0
```

The `radiosonde_format` value should remain a loadable format. The fact that the file was manually provided is stored separately in `radiosonde_source`.

## Status codes

| Status | Meaning |
| ---: | --- |
| `0` | Radiosonde file found or manually accepted. |
| `1` | Supported radiosonde files were found, but none were within the automatic selection time window. |
| `2` | No supported radiosonde file was found, or a manually provided file path was invalid. |
| `-1` | Radiosonde handling is not required for this QA test. |

If no radiosonde can be selected, processing steps that require molecular profiles are skipped or cannot be performed. The Rayleigh-fit and polarization-calibration quicklook generators also skip their molecular-dependent plots and print a warning instead of failing when molecular products are unavailable.

## Recommended usage

Use `radiosonde_file` when reproducibility is more important than automatic selection, or when the correct radiosonde file is known beforehand.

Use `radiosonde_folder` with `rsonde_station_wmo_id` and/or `cloudnet_station_name` when ATLAS should download and select the best available profile automatically. Make sure the configured station identifiers match the identifiers used in supported filenames.

Avoid mixing unrelated campaigns, stations, or dates in the same `radiosonde_folder`. The selector ignores unsupported names and supported-looking names with mismatching identifiers, but all supported filenames with matching identifiers are candidates for automatic selection.

## Behavior when molecular products are missing

If no usable radiosonde is selected, `load_radiosonde()` cannot create meteorological profiles and molecular calculations are skipped. In that case, downstream molecular products such as attenuated molecular backscatter, molecular depolarization ratio, and molecular metadata are unavailable.

The molecular-dependent quicklooks are protected against this situation:

- Rayleigh fit quicklooks are skipped for a QA test if `molecular` is missing.
- Polarization calibration quicklooks are skipped if `molecular_ratio` or `molecular_info` is missing.

ATLAS prints a warning explaining that no molecular products were found and that this usually means no usable radiosonde was selected. This prevents errors from title generation, molecular overlays, or Rayleigh-fit normalization when radiosonde-dependent data do not exist.

## Troubleshooting

### A downloaded file was not selected

A download can succeed but still fail selection. Common reasons are:

- the downloaded file timestamp is outside the selection window;
- a Wyoming file inside the priority window was preferred;
- the filename does not match a supported pattern;
- the identifier part of the filename does not match the expected identifier for that format;
- the selected QA test midpoint is different from the expected time.

### Wyoming is slow on the first run

This is expected when ATLAS has to contact the Wyoming server and test candidate sounding times/sources. Later runs may be faster because existing downloaded files are reused and previous failed attempts are skipped using `.wyoming_download_failures.json`.

Delete `.wyoming_download_failures.json` from `radiosonde_folder` when you want ATLAS to retry requests that were previously marked as unavailable or failed.

### A file is treated as custom ASCII instead of Wyoming

Wyoming BUFR text files should use the pattern:

```text
YYYYMMDD_HHMM_wyoming_bufr_<wmo_id>.txt
```

For example:

```text
20240428_1200_wyoming_bufr_16716.txt
```

This is classified as `wyoming`, not `custom_ascii`.

### Manual file timestamp warning

If a manual file does not contain a timestamp in a supported filename pattern, ATLAS cannot infer `radiosonde_time` from the name. The file can still be used if one of the readers can parse it, but the measurement midpoint is stored as the radiosonde time metadata placeholder.

### Manual file cannot be parsed by any reader

Manual mode no longer relies only on the filename extension. It tries the available radiosonde readers and accepts the file only if one reader returns a valid meteorological profile. If all readers fail, the manual file is rejected and `radiosonde_status = 2`.

Check that:

- the file path is correct;
- the file content matches one of the supported reader formats;
- for `custom_ascii`, the `rsonde_*` column, delimiter, header, footer, and unit options match the file content;
- for `custom_ascii` heights in `m_agl` or `km_agl`, `rsonde_station_altitude` is provided.

### Cloudnet file downloaded but not selected

If a Cloudnet/ECMWF file is downloaded but then rejected, check the identifier in the filename and the configured `cloudnet_station_name`. For example, `20251203_1159_ecmwf_barcelona.nc` has identifier `barcelona`. Identifier comparison is case-normalized, so `Barcelona` should match, but a genuinely different station name will not.

