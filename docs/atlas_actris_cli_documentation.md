# ATLAS ACTRIS – Command Line Tools (WIP)

This document describes the current command-line interfaces (CLI) provided by the `atlas_actris` project.


---


# Available CLI Commands

The following entrypoints are defined:

- `atlas`
- `get_config`
- `get_T_P`

---

# 1) atlas_actris

Runs the main ATLAS workflow using an initialization (`.ini`) file.

## Usage

```bash
atlas -i /path/to/call_atlas.ini
```

## Arguments

- `-i`, `--ini_file` (required)
  
  Help: The path to the initialization file of the calling script

## Behavior

- Raises an error if `--ini_file` is not provided.
- Raises an error if the file does not exist.
- Parses the INI file and runs the ATLAS processing workflow.
- Generates QA/report outputs.

## Example

```bash
atlas --ini_file ./examples/call_atlas.ini
```

# 2) get_config

Exports an ATLAS configuration file based on SCC configuration information.

## Usage

```bash
get_config -o /path/to/output/folder [OPTIONS]
```

## Arguments

- `-i`, `--scc_configuration_id` (required)
  
  Help: The SCC configuration ID
  
  Default: None
  
- `-o`, `--hoi_output_folder` (optional)
  
  Help: The path to the Handbook of Instruments csv file will be downloaded


- `-c`, `--atlas_configuration_folder` (optional)
  
  Help: The path to the folder where the ATLAS config_file.ini will be exported
  
  Default: None

- `-v`, `--verbose` / `--no-verbose` (optional)
  
  Help: Provide to print information about the SCC HOI file
  
  Default: False

## Examples

Minimal:

```bash
get_config -i 665 -o ./hoi_export
```

With explicit configuration folder:

```bash
get_config -o ./hoi_export -c ./atlas_config
```

---

# 3) get_T_P

Extracts ECMWF temperature/pressure profiles for a Cloudnet station and exports to CSV radiosonde-like files.

## Usage

```bash
get_T_P STATION DATE TIME --save-dir /path/to/cache [OPTIONS]
```

## Positional Arguments

- `station` (required)
  
  Help: e.g. "Bucharest"

- `date` (required)
  
  Help: date "dd.mm.yyyy"

- `time` (required)
  
  Help: time "hh:mm:ss" UTC

## Options

- `--save-dir` (required)
  
  Help: cache directory for .nc files

- `--src` (optional)
  
  Help: optional seed .nc to copy if cache empty
  
  Default: None

- `--nc` (optional)
  
  Help: optional seed .nc to copy if cache empty
  
  Default: None

- `--outcsv-dir` (optional)
  
  Help: CSV export folder path (defaults to the save-dir folder)
  
  Default: None (falls back to save-dir)

## Behavior

- Creates the output CSV folder if it does not exist.
- Uses cached NetCDF files if available.
- Prints model time, number of levels, and file used.

## Examples

Export CSVs into cache folder:

```bash
get_T_P Bucharest 01.01.2025 00:00:00 --save-dir ./cache
```

Export CSVs into separate folder:

```bash
get_T_P Bucharest 01.01.2025 00:00:00 --save-dir ./cache --outcsv-dir ./csv
```

---

