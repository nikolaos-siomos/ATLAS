#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""User-facing flavor text for generated ATLAS INI templates.

This module is intentionally simple.  The parser schema remains the technical
source of truth for valid keys, types, defaults, mandatory/optional status, and
allowed values.  The dictionaries below provide only user-facing text used by
the template/documentation generator:

- description: short explanatory comment shown above the template entry.
- example: optional example value shown in templates and MkDocs.
- legacy: structured backward-compatibility information.

The legacy dictionary uses the following fields:

- status: one of "", "new", "unchanged", "renamed", "moved".
- introduced: ATLAS version where the current parameter/location was introduced.
- old_names: previous parameter names, if any.
- old_location: previous INI file/section if the parameter was moved.
- old_names_removed_in: ATLAS version where old names stopped being accepted.
- note: optional compatibility explanation.
"""

from __future__ import annotations

from utils.parse_init_file import (
    slice_exclude_allowed_keys, 
    version_warning_nrm,
    version_warning_pcb,
    INIT_FILE_SECTIONS,
    )

try:
    from atlas_actris.version import __version__ as ATLAS_VERSION
except Exception:  # pragma: no cover - fallback for standalone editing
    ATLAS_VERSION = "unknown"

# The generated template should report this ATLAS version in its header.
TEMPLATE_ATLAS_VERSION = ATLAS_VERSION

LEGACY_STATUS_VALUES = {"", "new", "unchanged", "renamed", "moved"}

locations = {
    'init_file':'initialization file',
    'config_file':'configuration file',
    'settings_file': 'settings_file'
    }

INIT_TEMPLATE_SECTIONS = INIT_FILE_SECTIONS
# INIT_TEMPLATE_SECTIONS = {'configuration': ['scc_compatible_format', 'export_hoi_cfg', 'scc_configuration_id'],
#  'explicit_paths': ['parent_folder',
#                     'atlas_configuration_file',
#                     'atlas_settings_file',
#                     'radiosonde_folder',
#                     'radiosonde_file',
#                     'output_folder'],
#  'general_options': ['process',
#                      'process_qck',
#                      'vertical_scale',
#                      'view_mean_signal_stages',
#                      'view_signal_stages',
#                      'dpi',
#                      'color_reduction',
#                      'overwrite_output',
#                      'expert_analyst',
#                      'export_stages',
#                      'export_all'],
#  'filter_channels': ['select_channels',
#                      'exclude_telescope_type',
#                      'exclude_channel_type',
#                      'exclude_acquisition_mode',
#                      'exclude_channel_subtype'],
#  'trimming_options': ['max_height_agl',
#                       'low_shot_threshold',
#                       'trim_overflows',
#                       'low_res_averaging_period',
#                       'low_res_averaging_threshold',
#                       'high_res_averaging_period',
#                       'high_res_averaging_threshold',
#                       'max_adjacent_overflows',
#                       'slice_measurement',
#                       'exclude_measurement'],
#  'explicit_folders': ['ray',
#                       'pcb',
#                       'tlc',
#                       'tlc_rin',
#                       'drk',
#                       'trg',
#                       'dtm',
#                       'ray_pcb',
#                       'pcb_aux',
#                       'cam'],
#  'parsing_options': ['files_per_quadrant',
#                      'files_per_ring',
#                      'rsonde_skip_header',
#                      'rsonde_skip_footer',
#                      'rsonde_delimiter',
#                      'rsonde_column_index',
#                      'rsonde_column_units',
#                      'rsonde_station_altitude',
#                      'cloudnet_station_name',
#                      'rsonde_station_name',
#                      'rsonde_station_wmo_id']}

rel_path_note = 'Relative paths are provided with respect'
'to the root folder of the initialization file used to run ATLAS.'

legacy_atlas = '<=0.5.0'
last_legacy_atlas = '0.6.6'

def alias_folder_text(qa_test):
    
    text = f'Alias folder name for the {qa_test} folder. '
    'Providing will force ATLAS to read the {qa_test} input data from '
    'a folder with the given name placed in the parent folder.'
    
    return text

def not_supported_text(qa_test):

    text = f'{qa_test} measurements can be provided to inspect signals '
    'with the signal viewer but the trg test is not operation yet.'
    
    return text

INIT_FLAVOR = {
    
    'scc_compatible_format': {
        'description': 
            'Set to True if the data are in scc format and '
            'export_hoi_cfg is 1 or 3 to generate the config file by '
            'automatically connecting the scc_channel_ids and'
            'the recorder_channel_ids.',
        'example': 'True',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'export_hoi_cfg': {
        'description': 
            'Controls exporting the config_file from the SCC HOI '
            'for the provided scc_configuration_id. Select among 0, 1, or 2 | '
            '0: A local, manually created, configuration ini file will be used, '
            '1: The config_file file will be automatically created from '
            'the SCC HOI overwritting any exported config_file '
            'for the same configuration, '
            'prepared file, 1: export from HOI and overwrite existing, '
            '2: The config_file file will be automatically created from '
            'the SCC HOI if no file has already be exported for the same '
            'configuration.',
        'example': '0',
        'legacy': {
            'status': 'unchanged',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'scc_configuration_id': {
        'description': 
            'The SCC configuration ID. It is used only if  '
            'export_hoi_cfg is set to 1 or 2 to connect with the SCC HOI.',
        'example': '665',
        'legacy': {
            'status': 'unchanged',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'parent_folder': {
        'description': 
            'Absolute or relative path of the parent folder which contains '
            'the input lidar data. ' + rel_path_note
            ,
         'example': 
             './179_199_665_20231221 - ATLAS will look for a folder named '
             '179_199_665_20231221 placed in the same directory as the '
             'call_atlas.ini file which was used to run ATLAS',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
    
    'atlas_configuration_file': {
        'description': 
            'Absolute or relative path of the atlas configuration file which '
            'contains all system-, channel-, and channel-pair-related '
            'metadata. ' + rel_path_note,
         'example': 
             'my_drive/configurations/config_file_665_20231221',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'atlas_settings_file': {
        'description': 
            'Absolute or relative path of the atlas settings file which '
            'contains all system-, channel-, and channel-pair-related '
            'metadata. ' + rel_path_note,
         'example': 
             'my_drive/settings/config_file_665_20231221',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
    
    'radiosonde_folder': {
        'description': 
            'Absolute or relative path of the folder where the radiosonde '
            'files are placed. It will be ignored if radiosonde_file option '
            'is provided.' + rel_path_note,
         'example': 
             'my_drive/radiosondes/',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },

    'radiosonde_file': {
        'description': 
            'Absolute or relative path to the radiosonde file.'
            'files are placed. ' + rel_path_note,
         'example': 
             'my_drive/radiosondes/20231221_1527_ecmwf_bucharest.nc',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'process': {
        'description': 
            'The user can choose specific QA test(s) to process. '
            'Choose one or more of the allowed keys to process the '
            'corresponding QA tests. If left empty, ATLAS keeps the standard '
            'default behaviour. If one or more tests are provided, empty '
            'process_qck and process_bgd entries are derived automatically '
            'from this selection.',
        
         'example': 'ray, pcb, tlc',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'process_qck': {
         'description': 
             'Generate profile quicklooks for specific QA tests. When process '
             'is explicitly selected and process_qck is left empty, this list '
             'is derived automatically from process. An explicit process_qck '
             'value overrides that derived list. VLDR quicklooks are controlled '
             'separately by process_vldr.',
         'example': 'ray, drk',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'process_bgd': {
         'description': 
             'Generate background time-series plots for specific QA tests. '
             'When process is explicitly selected and process_bgd is left '
             'empty, this list is derived automatically from process. An '
             'explicit process_bgd value overrides that derived list.',
         'example': 'ray, drk',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.1',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'process_vldr': {
         'description':
             'Generate VLDR quicklook plots independently of process and '
             'process_qck. VLDR quicklooks are enabled by default and can be '
             'disabled explicitly with process_vldr = False.',
         'example': 'True',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.1',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },

    'process_dedicated_dark': {
         'description':
             'Control automatic dedicated-dark companions for explicitly '
             'selected QA tests. If True, process = ray automatically adds '
             'drk_ray to derived quicklook/background selections, process = '
             'tlc adds drk_tlc, and so on. If False, those automatic drk_* '
             'companions are omitted. The standalone drk test is never removed '
             'by this option when drk is selected. If process is empty, the '
             'standard independent defaults are used.',
         'example': 'True',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.1',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'vertical_scale': {
        'description': 
            'Select the vertical scale used for all QA tests. Select one of:'
            'range: Use range from the lidar in meters as the vertical scale'
            'height_agl: Use height above ground level in meters as the '
            'vertical scale'
            'height_asl: Use height above sea level in meters as the '
            'vertical scale'
            ,
        'example': 'height_agl',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': ['use_range'],
            'old_location': locations['settings_file'],
            'old_names_removed_in': '1.0.0',
            'note': ''
            }
        },

    'view_mean_signal_stages': {
        'description': 
            'Option used only providing the init file to __signal_viewer__.'
            'The user can provide one or more stages for which time-averaged '
            'signals will be ploted for each channel and each QA test.'
            'If no time-averaged signals exist for a selected stage a '
            'warning will be raised and no plots will be created.'
            ,
        'example': 'preprocessing_complete',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
        
    'view_signal_stages': {
        'description': 
            'Option used only providing the init file to __signal_viewer__.'
            'The user can provide one or more stages for which time-resolved '
            'signals will be ploted for each channel and each QA test. '
            'Each signal is plotted with a separate line.'
            ,
        'example': 'init, preprocessing_complete - this will plot all raw signals and all range-corrected signals',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
     
    'dpi': {
        'description': 
            'Resolution of exported figures in dots per inch.',
        'example': '300',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': locations['settings_file'],
            'old_names_removed_in': '1.0.0',
            'note': ''
            }
        },
        
    'color_reduction': {
        'description': 
            'If set to False, the plot colors are not reduced to decrease '
            'plot size.',
        'example': 'False',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': locations['settings_file'],
            'old_names_removed_in': '1.0.0',
            'note': 
                'Built-in color reduction performed directly in python '
                'with PIL.'
                'Imagemagick is no longer used.'
            }
        },

    'output_folder': {
        'description': 
            'Path to the folder where ATLAS output files are exported '
            '(plots, ascci, cached, and exported stage files).'
            'The files are placed under a subfolder with the same name as the '
            'parent folder.',
        'example': '/my_data/brc_analysis',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': '',
            }
        },

    'overwrite_output': {
        'description': 
            "If set to True, ATLAS writes outputs directly into the existing "
            "'plots' and 'ascii' folders, overwriting files with "
            "the same names. "
            "If set to False, ATLAS creates new timestamped subfolders "
            "for each run ('plots_<timestamp>' and 'ascii_<timestamp>') "
            "and writes the new outputs there.",
        'example': '/my_data/brc_analysis',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': '',
            }
        },

    'expert_analyst': {
        'description': 
            'Name or identifier of the analyst responsible for '
            'the processing. It is used in the filename of the generated '
            'QA reports',
        
         'example': 'ns',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },

    'export_stages': {
        'description': 
            'At the end of ATLAS excecution, the user'
            'is asked to export all export_stages. A rough '
            'calculation of the size of the data to be exported '
            'is provided. Exporting can take a while, depending '
            'on the volume of data and number of stages.',
        'example': 'preprocessing_complete',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': '',
            }
        },

    'export_all': {
        'description': 
            'If set to True, all channels are exported in the QA reports '
            '(e.g analog channels for Rayleigh fit test and photon channels '
            'for the telecover test)',
         'example': 'True',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'select_channels': {
        'description': 
            'Provide ATLAS channel IDs of the channels to be processed.'
            'All channels are processed by default. Excluding options are '
            'applied after selecting channels',
         'example': '0355xtax, 0355xtpx',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': ['channels'],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },
        
    'exclude_wavelength': {
        'description': 
            'Provide the wavelegth part of the ATLAS channel ID '
            '(first 4 characters) of the channels to be excluded '
            'before processing.'
            'All channels are processed by default.',
         'example': '0355, 1064',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'exclude_telescope_type': {
        'description': 
            'Provide the telescope type part of the ATLAS channel ID '
            '(5th character) of the channels to be excluded '
            'before processing.'
            'All channels are processed by default.',
         'example': 'n - near range telescope channels will be excluded',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },


    'exclude_channel_type': {
        'description': 
            'Provide the channel type part of the ATLAS channel ID '
            '(6th character) of the channels to be excluded '
            'before processing.'
            'All channels are processed by default.',
         'example': 'v, f - vibrational Raman and fluorescence channels will be excluded',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },
        
    'exclude_acquisition_mode': {
        'description': 
            'Provide the channel acquisition_mode part of the ATLAS channel ID '
            '(7th character) of the channels to be excluded '
            'before processing.'
            'All channels are processed by default.',
         'example': 'a - analogue channels will be excluded',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },

    'exclude_channel_subtype': {
        'description': 
            'Provide the channel subtype part of the ATLAS channel ID '
            '(8th character) of the channels to be excluded '
            'before processing.'
            'All channels are processed by default.',
         'example': 'a - analog channels will be excluded',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },

    'max_height_agl': {
        'description': 
            'The maximum height in km agl above which signal will be trimmed '
            'after the signal trimming is applied. ',
         'example': '40',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': ['vertical_trimming, vertical_limit'],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },

    'low_shot_threshold': {
        'description': 
            'The fraction of profile shots divided by the max number of shots. '
            'Profiles with shots less than low_shot_threshold times the max '
            'number of shots will be masked out by the screen low shots stage.',
         'example': '40',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
        
    'trim_overflows': {
        'description': 
            'This options determines how overflow values in the raw input '
            'files will be treated if detected. Choose among: '
            '0: the algorithm will stop and provide a diagnostic error, '
            'if overflows are found, '
            '1: the files containing at least one overflow value will be '
            'screened out, '
            '2: overflows will be interpolated from neighboring bins.'
            '3: overflows will not be masked out.',
         'example': '2',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },

    'low_res_averaging_period': {
        'description': 
            'Not applied yet for QA - '
            'Averaging period (low temporal resolution) for Rayleigh and dark '
            'measurements, specified in minutes or hours, '
            'for example 10min or 2h. '
            'If the measurement duration is shorter than the averaging period, '
            'the full dataset is averaged. '
            'Otherwise, data are averaged over equal time intervals. '
            'Averages are masked when the fraction of missing profiles within '
            'an interval exceeds low_res_averaging_threshold.',
         'example': '1h',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },  
        
    'low_res_averaging_threshold': {
        'description': 
            'Not applied for QA yet - '
            'Threshold for masking averages based on missing data within the '
            'averaging time interval. The interval is defined by '
            'low_res_averaging_period; averages are masked when the fraction of '
            'missing profiles relative to the expected number of profiles '
            'exceeds low_res_averaging_threshold.',
         'example': '1h',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },     

    'high_res_averaging_period': {
        'description': 
            'Not applied for QA yet - '
            'Averaging period (high temporal resolution) for Rayleigh and dark '
            'measurements, specified in minutes or hours, '
            'for example 10min or 2h. '
            'Decimal values are also accepted: e.g. 0.5min --> 30 seconds, '
            '0.5h --> 30 minutes'
            'If the measurement duration is shorter than the averaging period, '
            'the full dataset is averaged. '
            'Otherwise, data are averaged over equal time intervals. '
            'Averages are masked when the fraction of missing profiles within '
            'an interval exceeds high_res_averaging_threshold.',
         'example': '0.5min',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },  
        
    'high_res_averaging_threshold': {
        'description': 
            'Not applied for QA yet - '
            'Threshold for masking averages based on missing data within the '
            'averaging time interval. The interval is defined by '
            'high_res_averaging_threshold; averages are masked when the '
            'fraction of missing profiles relative to the expected number of '
            'profiles exceeds high_res_averaging_threshold.',
         'example': '1h',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },     

    'max_adjacent_overflows': {
        'description': 
            'Therehold applied when overflow interpolation is attempted: '
            'trim_overflows = 2 '
            'If the number of adjucent bins with overflow values in a '
            'measurement profiled exceeds max_adjacent_overflows, the whole '
            'profile is masked (similar to trim_overflows = 1)',
         'example': '3',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         }, 

    'slice_measurement': {
        'description': 
            'Select temporal regions to process using repeating triplets:'
            '<test_folder_alias>, <start_time>, <stop_time> '
            'The same test folder alias can appear multiple times. '
            'Accepted time formats: '
            'HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, yyyymmdd_HHMMSS '
            'For HHMM only, intervals crossing midnight are handled '
            'automatically, e.g. 2300 to 0135 means next day. '
            'For absolute formats, write the next date explicitly. '
            f'Available test folder aliases: {slice_exclude_allowed_keys}',
         'example': 'slice_measurement = drk, 2300, 0135, ray, 20251204_2300, 20251205_0600, tlc, 20251205_071510, 20251205_071545, ray_pcb, 20251205_08, 20251205_09',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': locations['settings_file'],
             'old_names_removed_in': '1.0.0',
             'note': ''
             }
         },
        
    'exclude_measurement': {
        'description': 
            'Select temporal regions to exclude using repeating triplets:'
            '<test_folder_alias>, <start_time>, <stop_time> '
            'The same test folder alias can appear multiple times. '
            'Accepted time formats: '
            'HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, yyyymmdd_HHMMSS '
            'For HHMM only, intervals crossing midnight are handled '
            'automatically, e.g. 2300 to 0135 means next day. '
            'For absolute formats, write the next date explicitly. '
            f'Available test folder aliases: {slice_exclude_allowed_keys}',
         'example': 'exclude_measurement = drk, 2300, 0135, ray, 20251204_2300, 20251205_0600, tlc, 20251205_071510, 20251205_071545, ray_pcb, 20251205_08, 20251205_09',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },

     'ray': {
         'description': alias_folder_text('ray'),
         'example': 'ray_02',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': version_warning_nrm
             }
         },
 
    'pcb': {
        'description': alias_folder_text('pcb'),
        'example': 'pcb_filter_02',
        'legacy': {
            'status': 'unchanged',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': version_warning_pcb
            }
        },

     'tlc': {
         'description': alias_folder_text('tlc'),
         'example': 'tlc',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
     
     'tlc_rin': {
         'description': alias_folder_text('tlc_rin'),
         'example': 'tlc_03',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
     
     'drk': {
         'description': alias_folder_text('drk'),
         'example': 'drk_01',
         'legacy': {
             'status': 'unchanged',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
     
     'trg': {
         'description': alias_folder_text('trg'),
         'example': 'trg_raman',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': not_supported_text('trg (zero bin measuremnt)')
             }
         },
     
     'dtm': {
         'description': alias_folder_text('dtm'),
         'example': 'dtm_02',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': not_supported_text('dtm (deat time measuremnt')
             }
         },
 
    'ray_pcb': {
        'description': alias_folder_text('dtm'),
        'example': 'ray_pcb_01',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': not_supported_text('ray_pcb (Rayleigh measuremnt in pol. cal. mode)')
            }
        },

     'pcb_aux': {
         'description': alias_folder_text('pcb_aux'),
         'example': 'pcb_aux',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': not_supported_text('pcb_aux (Pol. Cal. measuremnt for ND filter characterization)')
             }
         },
     
     'cam': {
         'description': alias_folder_text('cam'),
         'example': 'cam',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': not_supported_text('cam (camera images)')
             }
         },
     
     'files_per_quadrant': {
         'description': 
             'Number of consecutive telecover files per quadrant sector when '
             'files must be automatically distributed. The files are '
             'distributed automatically in subfolders (north, east, south, west)'
             'and are read from there for any subsequent run, automatically.',
         'example': '4',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': ['files_per_sector'],
             'old_location': 'settings_file',
             'old_names_removed_in': '1.0.0',
             'note': 'In older versions the files where not distributed into '
             'sector subfolders. This changed in version 1.0.0'
             }
         },
         
     'files_per_ring': {
         'description': 
             'Number of consecutive telecover files per quadrant sector when '
             'files must be automatically distributed. The files are '
             'distributed automatically in subfolders (outer, inner)'
             'and are read from there for any subsequent run, automatically.',
         'example': '4',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': 'settings_file',
             'old_names_removed_in': '',
             'note': 'In older versions the files where not distributed into '
             'sector subfolders. This changed in version 1.0.0'
             }
         },

     'rsonde_skip_header': {
         'description': 
             'Number of header lines to skip when reading an ASCII radiosonde '
             'file with custom format.',
         'example': '1',
         'legacy': {
             'status': 'moved',
             'introduced': legacy_atlas,
             'old_names': [],
             'old_location': 'settings_file',
             'old_names_removed_in': '',
             'note': ''
             }
         },

    'rsonde_skip_footer': {
        'description': 
            'Number of footer lines to skip when reading an ASCII radiosonde '
            'file with custom format.',
        'example': '0',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': 'settings_file',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'rsonde_delimiter': {
        'description': 
            'Delimiter type used when reading an ASCII radiosonde '
            'file with custom format.',
        'example': 'S - space delimeter',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': 'settings_file',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'rsonde_column_index': {
        'description': 
            'Column indices for height, pressure, temperature, and optionally '
            'humidity. Used when reading an ASCII radiosonde.'
            'file with custom format.',
        'example': '2, 1, 3, 5',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': 'settings_file',
            'old_names_removed_in': '',
            'note': ''
            }
        },
        
    'rsonde_column_units': {
        'description': 
            'Units corresponding to the radiosonde columns. '
            'Used when reading an ASCII radiosonde ',
        'example': 'm_asl, hPa, C, percent',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': 'settings_file',
            'old_names_removed_in': '',
            'note': ''
            }
        },

    'rsonde_station_altitude': {
        'description': 
            'Altitude of the radiosonde station. Used when reading an ASCII '
            'radiosonde and agl height units are used.',
        'example': '60.',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': ['station_altitude'],
            'old_location': 'settings_file',
            'old_names_removed_in': '',
            'note': ''
            }
        },
        
    'cloudnet_station_name': {
        'description': 
            'Cloudnet station name used when Cloudnet. If provided, ATALS will '
            'attempt downloading ECMWF meteorological files from Cloudnet to '
            'extract the temperature, pressure, and relative humidity '
            'parameters. The name of the station must be provided exactly '
            'as define in the Cloudnet API. This is not necessary the same as '
            'the display name in the web interface.',
        'example': 'garmish',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
        
    'rsonde_station_name': {
        'description': 
            'Radiosonde station name displayed in plots. If an ECMW file is '
            'used, the station name is used for displaying instead and this '
            'parameter is ignored.',
        'example': 'Thessaloniki',
        'legacy': {
            'status': 'moved',
            'introduced': legacy_atlas,
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
        
    'rsonde_station_wmo_id': {
        'description': 
            'WMO identifier of the radiosonde station. If provided, ATALS will '
            'attempt downloading Wyoming radiosonde filed to '
            'extract the temperature, pressure, and relative humidity.',
        'example': '16622',
        'legacy': {
            'status': 'moved',
            'introduced': '',
            'old_names': ['rsonde_wmo_number'],
            'old_location': '',
            'old_names_removed_in': '',
            'note': 'Before ATLAS 1.0.0 this parameter was used for '
            'display. Since ATLAS 1.0.0 it is used for automatic downloading '
            'of Wyoming radiosonde files.'
            }
        }
    }

