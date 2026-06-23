#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""User-facing flavor text for generated ATLAS INI templates.

This module is intentionally simple.  The parser schema remains the technical
source of truth for valid keys, types, defaults, mandatory/optional status, and
allowed values.  The dictionaries below provide only the comment text and
example values used by the template generator.
"""

from __future__ import annotations

try:
    from atlas_actris.version import __version__ as ATLAS_VERSION
except Exception:  # pragma: no cover - fallback for standalone editing
    ATLAS_VERSION = "unknown"

# The generated template should report this ATLAS version in its header.
TEMPLATE_ATLAS_VERSION = ATLAS_VERSION

INIT_TEMPLATE_SECTIONS = {'configuration': ['scc_compatible_format', 'export_hoi_cfg', 'scc_configuration_id'],
 'explicit_paths': ['parent_folder',
                    'atlas_configuration_file',
                    'atlas_settings_file',
                    'radiosonde_folder',
                    'radiosonde_file',
                    'output_folder'],
 'general_options': ['process',
                     'process_qck',
                     'vertical_scale',
                     'dpi',
                     'color_reduction',
                     'overwrite_output',
                     'expert_analyst',
                     'debug_signals',
                     'export_netcdf',
                     'export_all'],
 'filter_channels': ['exclude_channels',
                     'exclude_telescope_type',
                     'exclude_channel_type',
                     'exclude_acquisition_mode',
                     'exclude_channel_subtype'],
 'trimming_options': ['max_height_agl',
                      'low_shot_threshold',
                      'trim_overflows',
                      'ray_averaging_rate',
                      'ray_averaging_threshold',
                      'ray_qck_averaging_rate',
                      'ray_qck_averaging_threshold',
                      'max_adjacent_overflows',
                      'slice_measurement',
                      'exclude_measurement'],
 'explicit_folders': ['ray',
                      'nrm',
                      'pcb',
                      'tlc',
                      'tlc_rin',
                      'drk',
                      'trg',
                      'dtm',
                      'nsf',
                      'ray_pcb',
                      'nrm_pcb',
                      'pcb_aux',
                      'cam'],
 'parsing_options': ['files_per_sector',
                     'files_per_ring',
                     'files_per_set',
                     'rsonde_skip_header',
                     'rsonde_skip_footer',
                     'rsonde_delimiter',
                     'rsonde_column_index',
                     'rsonde_column_units',
                     'rsonde_station_latitude',
                     'rsonde_station_longitude',
                     'rsonde_station_altitude',
                     'cloudnet_station_name',
                     'rsonde_station_name',
                     'rsonde_station_wmo_id']}

INIT_FLAVOR = {'scc_compatible_format': {'description': 'Set to True when the input data and '
                                          'metadata should be interpreted in '
                                          'SCC-compatible mode.',
                           'example': 'True'},
 'export_hoi_cfg': {'description': 'Controls whether the SCC/HOI configuration file is '
                                   'used locally or exported/downloaded before the '
                                   'run.',
                    'example': '0'},
 'scc_configuration_id': {'description': 'SCC configuration identifier used when '
                                         'exporting or downloading a configuration '
                                         'file.',
                          'example': '665'},
 'parent_folder': {'description': 'Main folder containing the measurement subfolders '
                                  'used by ATLAS.',
                   'example': './179_199_665_20231221'},
 'atlas_configuration_file': {'description': 'Path to the ATLAS system configuration '
                                             'file.',
                              'example': './configurations/config_file_665.ini'},
 'atlas_settings_file': {'description': 'Path to the ATLAS plotting and QA-test '
                                        'settings file.',
                         'example': './settings/settings_file.ini'},
 'radiosonde_folder': {'description': 'Folder containing radiosonde or model profile '
                                      'files used for molecular calculations.',
                       'example': './radiosondes'},
 'radiosonde_file': {'description': 'Specific radiosonde or model profile file to use '
                                    'instead of automatic selection from a folder.',
                     'example': './radiosondes/20231221_1527_ecmwf_thessaloniki.nc'},
 'process': {'description': 'QA tests to process. Use off to disable processing.',
             'example': 'ray, pcb, tlc_qua, tlc_rin, drk'},
 'process_qck': {'description': 'Quicklook products to generate. Use off to disable '
                                'quicklook generation.',
                 'example': 'ray, pcb, tlc_qua, tlc_rin, drk'},
 'vertical_scale': {'description': 'Vertical coordinate used in plots and processing '
                                   'limits.',
                    'example': 'range'},
 'dpi': {'description': 'Resolution of exported figures in dots per inch.',
         'example': '300'},
 'color_reduction': {'description': 'Reduce the color palette of exported figures when '
                                    'supported.',
                     'example': 'True'},
 'output_folder': {'description': 'Base output folder where ATLAS creates analysis '
                                  'products for the current measurement.',
                   'example': './analysis'},
 'overwrite_output': {'description': 'Use stable output subfolder names instead of '
                                     'timestamped plot and ASCII subfolders.',
                      'example': 'True'},
 'expert_analyst': {'description': 'Name or identifier of the analyst responsible for '
                                   'the processing.',
                    'example': 'Name Surname'},
 'debug_signals': {'description': 'Export additional intermediate/debug signal '
                                  'information when enabled.',
                   'example': 'True'},
 'export_netcdf': {'description': 'Export NetCDF output products when enabled.',
                   'example': 'True'},
 'export_all': {'description': 'Export all available products and intermediate outputs '
                               'when enabled.',
                'example': 'True'},
 'exclude_channels': {'description': 'Recorder channel IDs to exclude from all '
                                     'selected processing.',
                      'example': '0532xcar, 1064xtax'},
 'exclude_telescope_type': {'description': 'Telescope types to exclude from all '
                                           'selected processing.',
                            'example': 'n, f, x'},
 'exclude_channel_type': {'description': 'Channel types to exclude from all selected '
                                         'processing.',
                          'example': 'p, c, t'},
 'exclude_acquisition_mode': {'description': 'Acquisition modes to exclude from all '
                                             'selected processing.',
                              'example': 'a, p, g'},
 'exclude_channel_subtype': {'description': 'Channel subtypes to exclude from all '
                                            'selected processing.',
                             'example': 'r, t, n'},
 'max_height_agl': {'description': 'Maximum height above ground level used for '
                                   'trimming or plotting.',
                    'example': '40'},
 'low_shot_threshold': {'description': 'Fractional threshold used to identify low-shot '
                                       'measurements.',
                        'example': '0.9'},
 'trim_overflows': {'description': 'Overflow trimming mode used during preprocessing.',
                    'example': '0'},
 'ray_averaging_rate': {'description': 'Averaging rate for Rayleigh measurements.',
                        'example': '30min'},
 'ray_averaging_threshold': {'description': 'Minimum valid-data fraction used for '
                                            'Rayleigh averaging.',
                             'example': '1.0'},
 'ray_qck_averaging_rate': {'description': 'Averaging rate for Rayleigh quicklooks.',
                            'example': 'raw'},
 'ray_qck_averaging_threshold': {'description': 'Minimum valid-data fraction used for '
                                                'Rayleigh quicklook averaging.',
                                 'example': '1.0'},
 'max_adjacent_overflows': {'description': 'Maximum number of adjacent overflow bins '
                                           'tolerated before trimming.',
                            'example': '1'},
 'slice_measurement': {'description': 'Optional measurement time slicing instructions '
                                      'in groups of test identifier, start time, and '
                                      'stop time.',
                       'example': 'ray, 2330, 0100'},
 'exclude_measurement': {'description': 'Optional measurement time exclusion '
                                        'instructions in groups of test identifier, '
                                        'start time, and stop time.',
                         'example': 'tlc_north, 1200, 1215'},
 'ray': {'description': 'Relative folder name for Rayleigh measurements below '
                        'parent_folder.',
         'example': 'ray'},
 'nrm': {'description': 'Legacy relative folder name for normalized/Rayleigh '
                        'measurements below parent_folder.',
         'example': 'nrm'},
 'pcb': {'description': 'Relative folder name for polarization-calibration '
                        'measurements below parent_folder.',
         'example': 'pcb'},
 'tlc': {'description': 'Relative folder name for quadrant telecover measurements '
                        'below parent_folder.',
         'example': 'tlc'},
 'tlc_rin': {'description': 'Relative folder name for ring telecover measurements '
                            'below parent_folder.',
             'example': 'tlc_rin'},
 'drk': {'description': 'Relative folder name for common dark measurements below '
                        'parent_folder.',
         'example': 'drk'},
 'trg': {'description': 'Relative folder name for trigger-delay measurements below '
                        'parent_folder.',
         'example': 'trg'},
 'dtm': {'description': 'Relative folder name for dark-current/timing measurements '
                        'below parent_folder.',
         'example': 'dtm'},
 'nsf': {'description': 'Relative folder name for near-field/far-field related '
                        'measurements below parent_folder.',
         'example': 'nsf'},
 'ray_pcb': {'description': 'Relative folder name for Rayleigh measurements associated '
                            'with polarization calibration.',
             'example': 'ray_pcb'},
 'nrm_pcb': {'description': 'Legacy relative folder name for Rayleigh/normalization '
                            'measurements associated with polarization calibration.',
             'example': 'nrm_pcb'},
 'pcb_aux': {'description': 'Relative folder name for auxiliary '
                            'polarization-calibration measurements below '
                            'parent_folder.',
             'example': 'pcb_aux'},
 'cam': {'description': 'Relative folder name for camera measurements below '
                        'parent_folder.',
         'example': 'cam'},
 'files_per_sector': {'description': 'Number of consecutive telecover files per '
                                     'quadrant sector when files must be automatically '
                                     'distributed.',
                      'example': '3'},
 'files_per_ring': {'description': 'Number of consecutive telecover files per ring '
                                   'when files must be automatically distributed.',
                    'example': '3'},
 'files_per_set': {'description': 'Number of consecutive files per measurement set '
                                  'when applicable.',
                   'example': '3'},
 'rsonde_skip_header': {'description': 'Number of header lines to skip when reading an '
                                       'ASCII radiosonde file.',
                        'example': '1'},
 'rsonde_skip_footer': {'description': 'Number of footer lines to skip when reading an '
                                       'ASCII radiosonde file.',
                        'example': '0'},
 'rsonde_delimiter': {'description': 'Delimiter type used in ASCII radiosonde files.',
                      'example': 'S'},
 'rsonde_column_index': {'description': 'Column indices for height, pressure, '
                                        'temperature, and optionally humidity in ASCII '
                                        'radiosonde files.',
                         'example': '2, 1, 3, 5'},
 'rsonde_column_units': {'description': 'Units corresponding to the radiosonde '
                                        'columns.',
                         'example': 'm_asl, hPa, C, percent'},
 'rsonde_station_latitude': {'description': 'Latitude of the radiosonde station.',
                             'example': '40.63'},
 'rsonde_station_longitude': {'description': 'Longitude of the radiosonde station.',
                              'example': '22.96'},
 'rsonde_station_altitude': {'description': 'Altitude of the radiosonde station.',
                             'example': '60'},
 'cloudnet_station_name': {'description': 'Cloudnet station name used when Cloudnet '
                                          'profiles are searched or downloaded.',
                           'example': 'thessaloniki'},
 'rsonde_station_name': {'description': 'Radiosonde station name used for profile '
                                        'selection or metadata.',
                         'example': 'Thessaloniki'},
 'rsonde_station_wmo_id': {'description': 'WMO identifier of the radiosonde station.',
                           'example': '16622'}}
