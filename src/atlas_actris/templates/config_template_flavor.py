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
CONFIG_FLAVOR = {
    'station_id': {
        'description': 'Station ID (e.g. SCC/EARLINET station ID)',
                'example': 'the',
        'legacy': {
            'status': '',
                           'introduced': '',
                           'old_names': [],
                           'old_location': '',
                           'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'lidar_name': {
        'description': 'Human-readable name of the lidar system.',
                'example': 'THELISYS',
        'legacy': {
            'status': '',
                           'introduced': '',
                           'old_names': [],
                           'old_location': '',
                           'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'station_name': {
        'description': 'Human-readable name of the station.',
                  'example': 'Thessaloniki',
        'legacy': {
            'status': '',
                             'introduced': '',
                             'old_names': [],
                             'old_location': '',
                             'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'lidar_id': {
        'description': 'SCC lidar ID',
              'example': '199',
        'legacy': {
            'status': '',
                         'introduced': '',
                         'old_names': [],
                         'old_location': '',
                         'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'version_name': {
        'description': 'Human-readable version name of the lidar '
                                 'configuration.',
        'example': 'default',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'version_id': {
        'description': 'SCC version ID',
                'example': '1',
        'legacy': {
            'status': '',
                           'introduced': '',
                           'old_names': [],
                           'old_location': '',
                           'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'configuration_name': {
        'description': 'Human-readable name of the configuration.',
        'example': 'Nighttime',
        'legacy': {
            'status': '',
                                   'introduced': '',
                                   'old_names': [],
                                   'old_location': '',
                                   'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'configuration_id': {
        'description': 'SCC configuration ID',
                      'example': '665',
        'legacy': {
            'status': '',
                                 'introduced': '',
                                 'old_names': [],
                                 'old_location': '',
                                 'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'station_altitude': {
        'description': 'Station altitude above sea level in meters.',
                      'example': '60',
        'legacy': {
            'status': 'renamed',
                                 'introduced': '1.0.0',
                                 'old_names': ['altitude'],
                                 'old_location': '',
                                 'old_names_removed_in': '1.0.0',
                                 'note': 'Parameter renamed during the ATLAS 1.0.0 '
                'schema cleanup.'
                }
        },
    
    'station_latitude': {
        'description': 'Station latitude in decimal degrees.',
                      'example': '40.63',
        'legacy': {
            'status': 'renamed',
                                 'introduced': '1.0.0',
                                 'old_names': ['latitude'],
                                 'old_location': '',
                                 'old_names_removed_in': '1.0.0',
                                 'note': 'Parameter renamed during the ATLAS 1.0.0 '
                'schema cleanup.'
                }
        },
    
    'station_longitude': {
        'description': 'Station longitude in decimal degrees.',
                       'example': '22.96',
        'legacy': {
            'status': 'renamed',
                                  'introduced': '1.0.0',
                                  'old_names': ['longitude'],
                                  'old_location': '',
                                  'old_names_removed_in': '1.0.0',
                                  'note': 'Parameter renamed during the ATLAS 1.0.0 '
                'schema cleanup.'
                }
        },
    
    'zenith_angle': {
        'description': 'Lidar pointing zenith angle in degrees.',
                  'example': '0',
        'legacy': {
            'status': '',
                             'introduced': '',
                             'old_names': [],
                             'old_location': '',
                             'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'azimuth_angle': {
        'description': 'Lidar pointing azimuth angle in degrees.',
                   'example': '0',
        'legacy': {
            'status': '',
                              'introduced': '',
                              'old_names': [],
                              'old_location': '',
                              'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'recorder_channel_id': {
        'description': 'Recorder channel IDs as they appear in the '
            'raw lidar files. For PollyXT raw files use ascending numbers '
            'starting from 1 as the recorder channel ID. For SCC raw files '
            'the recorder_channel_id corresponds to the scc_channel_id.',
        'example': 'BT0, BC0',
        'legacy': {
            'status': '',
                                    'introduced': '',
                                    'old_names': [],
                                    'old_location': '',
                                    'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'scc_channel_id': {
        'description': 'SCC channel IDs corresponding to each '
        'recorder_channel_id entry.',
                    'example': '1893, 1895, 1700',
        'legacy': {
            'status': '',
                               'introduced': '',
                               'old_names': [],
                               'old_location': '',
                               'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'telescope_type': {
        'description': 'ATLAS telescope type identifier for each channel. Choose among:'
            '\n#• n: first near range telescope (sort by descending expected distance of full overlap order)'
            '\n#• m: second near range telescope' 
            '\n#• l: third near range telescope'
            '\n#• f: first far range telescope (sort by ascending expected distance of full overlap order)'
            '\n#• g: second far range telescope' 
            '\n#• h: third far range telescope'
            '\n#• x: first standalone telescope '
            '\n#• y: second standalone telescope '
            '\n#• z: third standalone telescope ',
            'example': 'x, x, x, x, x, x, y, y, y, y, y, y, z, z',
        'legacy': {
            'status': '',
                               'introduced': '',
                               'old_names': [],
                               'old_location': '',
                               'old_names_removed_in': '',
            'note': ''
            }
        },

    'channel_type': {
        'description': 'ATLAS channel type identifier for each channel. Choose among:'
            '\n#• p: co-polar'
            '\n#• c: cross-polar'
            '\n#• t: total (no depolarization)'
            '\n#• v: vibrational Raman'
            '\n#• r: rotational Raman'
            '\n#• a: Cabannes (HSRL)'
            '\n#• f: fluorescence'
            '\n#• d: dual field of view',
        'example': 'c, c, p, p, v, v, c, c, p, p, r, r, t, t',
        'legacy': {
            'status': '',
                             'introduced': '',
                             'old_names': [],
                             'old_location': '',
                             'old_names_removed_in': '',
            'note': ''
            }
        },

    'acquisition_mode': {
        'description': 'ATLAS acquisition mode identifier for each channel. Use: a for analog channels, p for photon channels'
            '\n#• a: for analog channels'
            '\n#• p: for photon channels',
            'example': 'a, p, a, p, a, p, a, p, a, p, a, p, a, p, a, p, a, p, a, a',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': 'Modified identifiers with respect to pre v1.0.0 versions. In the past identifier 0 was used for analog channels and identifier 1 was used for photon channels.'
            }
        },

    'channel_subtype': {
        'description': 'ACTRIS channel subtype identifier for each channel. Choose among:'
            '\n#• r: Signal Reflected from a PBS'
            '\n#• t: Signal Transmitted through a PBS'
            '\n#• n: N2 Ramal line'
            '\n#• o: O2 Ramal line'
            '\n#• w: H2O Ramal line'
            '\n#• c: CH4 Ramal line'
            '\n#• h: High Rotational Raman'
            '\n#• l: Low Rotational Raman'
            '\n#• a: Mie (aerosol) HSRL signal'
            '\n#• m: Molecular HSRL signal'
            '\n#• b: Broadband Fluorescence'
            '\n#• s: Spectral Fluorescence'
            '\n#• x: No specific subtype',
        'example': 'r, r, t, t, n, n, t, t, r, r, x, x, x, x',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
 'zero_bin': {'description': 'Bin index corresponding to the trigger/zero range for '
                             'each channel.',
              'example': '0, 0, 0',
              'legacy': {'status': 'renamed',
                         'introduced': '1.0.0',
                         'old_names': ['daq_trigger_offset'],
                         'old_location': '',
                         'old_names_removed_in': '1.0.0',
                         'note': 'Parameter renamed during the ATLAS 1.0.0 schema '
                                 'cleanup.'}},

#    'first_signal_rangebin': {
#        'description': 'The first signal range bin index retrieved from the '
#        'zero bin test. ATLAS assumes 1 as the first bin index (0 or negative '
#        'values will not be accepted). '
#        'If the data acquisition trigger is delayed with respect to the '
#        'Q-switch pulse for a certain channel, the corresponding '
#        'trigger delay must be provided by the trigger_delay parameter in nm.'
#        'In that case, the first_signal_rangebin must be set to 1 '
#        'Decimal values are accepted but will result to slower processing '
#        'because the data needs to be interpolated along the bins dimension.',
#        'example': '2000, 2000, 2000, 2000, 2000, 2000, 10, 10, 10, 10, 10, 10, 1, 1',
#        'legacy': {
#            'status': 'new',
#            'introduced': '1.2.0',
#            'old_names': [],
#            'old_location': '',
#            'old_names_removed_in': '1.2.0',
#            'note': 'Affiliated to zero_bin parameter, introduced in v1.0.0, '
#            '(former legacy daq_trigger_offset). Since v1.2.0 '
#            'first_signal_rangebin and trigger_delay parameters have replaced '
#            'zero_bin parameter to improve compatibility with the SCC.'
#            }
#        },
    
#    'trigger_delay': {
#        'description': 'The trigger delay for channels for which the data '
#        'acquisition trigger starts later than the Q-switch pulse. '
#        'The trigger delay of channels whose data acquisition trigger '
#        'starts before the Q-switch trigger must be set to 0 '
#        'Values which result to decimal bins if divided with the '
#        'range_resolution are accepted but will result to slower processing '
#        'because the data needs to be interpolated along the bins dimension.',
#        'example': '0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 250, 275',
#        'legacy': {
#            'status': 'new',
#            'introduced': '1.2.0',
#            'old_names': [],
#            'old_location': '',
#            'old_names_removed_in': '1.2.0',
#            'note': 'Affiliated to zero_bin parameter, introduced in v1.0.0, '
#            '(former legacy daq_trigger_offset). Since v1.2.0 '
#            'first_signal_rangebin and trigger_delay parameters have replaced '
#            'zero_bin parameter to improve compatibility with the SCC.'
#            }
#        },
    

     'dead_time': {
         'description': 'Dead time of photon-counting channels in '
             'nanoseconds.',
         'example': '_, 3.7, _, 3.7, _, 3.7, _, 3.7, _, 3.7, _, 3.7, _, _',
         'legacy': {
             'status': '',
             'introduced': '',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },
             
    'background_low_bin': {
        'description': 'Lower bin index of the background region for '
            'each channel. Used for the background calculation/correction. '
            'The bin index corresponds to the original signal range bins '
            'prior to the zero bin correction. ATLAS assumes 1 as the first '
            'bin index.',
        'example': '100, 100, 100, 100, 100, 100, 100, 100,100, 100, 100, 100, 14000, 14000',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'background_high_bin': {
        'description': 'Upper bin index of the background region for '
            'each channel. Used for the background calculation/correction. '
            'The bin index corresponds to the original signal range bins '
            'prior to the zero bin correction. ATLAS assumes 1 as the first '
            'bin index.',
        'example': '900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 16000, 16000',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'emitted_wavelength': {
        'description': 'Emitted laser wavelength for each channel in nm. '
        'Defaults to 1064.14, 532.07, and 354.71 for 1064, 532, and 1064 '
        'related channels. '
        'It is recommended to either provide it with 2 decimal point accuracy '
        'or leave empty if it is not accurately known.',
        'example': '354.71, 354.71, 354.71, 354.71, 354.71, 354.71, 532.07, 532.07, 532.07, 532.07, 532.07, 532.07, 1064.14, 1064.14',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'detected_wavelength': {
        'description': 'Detected wavelength for each channel in nm. '
        'It corresponds to the central wavelength of the interference filter. '
        'It is recommended to provide it with 2 decimal point accuracy.',
        'example': '354.75, 354.75, 354.73, 354.73, 386.70, 386.70, 532.02, 532.02, 532.04, 532.04, 607.40, 607.40, 1064.10, 1064.10',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'channel_bandwidth': {
        'description': 'Interference-filter bandwidth for each channel in nm.'
        'It is recommended to provide it with 2 decimal point accuracy.',
        'example': '0.55, 0.55, 0.51, 0.51, 1.10, 1.10, 0.50, 0.50, 0.51, 0.51, 1.15, 1.15, 1.05, 1.05',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'G': {
        'description': 'Polarization cross talk parameter G for each channel. Use 1. for non-polarization channels if the cross-talk is unknown.',
        'example': '1.002, 1.002, 1.002, 1.002, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'H': {
        'description': 'Polarization cross talk H for each channel. Use 0. for non-polarization channels if the cross-talk is unknown.',
        'example': '-0.987, -0.987, 0.987, 0.987, 0., 0., -0.97, -0.97, 0.97, 0.97, 0., 0., 0., 0.',
        'legacy': {
            'status': 'new',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },

    'bins': {
        'description': 'Number of range bins for each channel.',
        'example': '16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'data_acquisition_range': {
        'description': 'Data acquisition range for each analog channel in mV. '
        'Use _ for photon channels.',
        'example': '100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 500., 500.',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'range_resolution': {
        'description': 'Range resolution for each channel in meters.',
        'example': '3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 7.5, 7.5',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },
    
    'laser_repetition_rate': {
        'description': 'Laser repetition rate for each channel in Hz.',
        'example': '20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 100, 100',
        'legacy': {
            'status': '',
            'introduced': '',
            'old_names': [],
            'old_location': '',
            'old_names_removed_in': '',
            'note': ''
            }
        },

    'ch_r': {
        'description': 'ATLAS channel ID for the reflected polarization '
            'channels. '
            'If unsure of the actual ATLAS channel ID, process once '
            'without filling up this section. The ATLAS channel ID is reported '
            ' in the header of the plots for each channel.',
        'example': '0355xcar, 0355xcpr, 0532ypar, 0532ypar',
        'legacy': {
            'status': 'moved',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': locations['settings_file'],
            'old_names_removed_in': '1.0.0',
            'note': 'Since v1.0.0 the parameters related to polarazation '
                f"channel pairs have been moved to the {locations['config_file']}."
            }
        },
    
    'ch_t': {
        'description': 'ATLAS channel ID for the transmitted polarization '
            'channels. '
            'If unsure of the actual ATLAS channel ID, process once '
            'without filling up this section. The ATLAS channel ID is reported '
            'in the header of the plots for each channel.',
        'example': '0355xpat, 0355xppt, 0532ycat, 0532ycat',
        'legacy': {
            'status': 'moved',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': locations['settings_file'],
            'old_names_removed_in': '1.0.0',
            'note': 'Since v1.0.0 the parameters related to polarazation '
                f"channel pairs have been moved to the {locations['config_file']}."
            }
        },
    
    'K': {
        'description': 'Cross talk parameter K for polarization calibration.',
        'example': '1.1520, 1.1520, 0.9804, 0.9804',
        'legacy': {
            'status': 'moved',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': locations['settings_file'],
            'old_names_removed_in': '1.0.0',
            'note': 'Since v1.0.0 the parameters related to polarazation '
                f"channel pairs have been moved to the {locations['config_file']}."
            }
        },
    
    'R_to_T_transmission_ratio': {
        'description': 'Transmission ratio between reflected '
            'and transmitted polarization paths. Provide when different ND '
            'filters are used for the polarization calibration and the '
            'Rayleigh measurements.',
        'example': '5.64, 5.64, 10.7, 10.7',
        'legacy': {
            'status': 'moved',
            'introduced': '1.0.0',
            'old_names': [],
            'old_location': locations['settings_file'],
            'old_names_removed_in': '1.0.0',
            'note': 'Since v1.0.0 the parameters related to polarazation '
                f"channel pairs have been moved to the {locations['config_file']}."
            }
        },
    
    'eta': {
        'description': 'Polarization calibration factor eta. Provide it only '
            'in case there is no dedicated polarisation calibration '
            'measurement but the calibration factor is known from a'
            'polarization calibration performed on a different day. '
            'This parameter is mainly used for testing.',
        'example': '1.0',
         'legacy': {
             'status': 'new',
             'introduced': '1.0.0',
             'old_names': [],
             'old_location': '',
             'old_names_removed_in': '',
             'note': ''
             }
         },

 'ch_n': {'description': 'Near-range channel IDs used for signal gluing.',
          'example': '0532xcnr',
          'legacy': {'status': '',
                     'introduced': '',
                     'old_names': [],
                     'old_location': '',
                     'old_names_removed_in': '',
                     'note': ''}},
 'ch_f': {'description': 'Far-range channel IDs used for signal gluing.',
          'example': '0532xcfr',
          'legacy': {'status': '',
                     'introduced': '',
                     'old_names': [],
                     'old_location': '',
                     'old_names_removed_in': '',
                     'note': ''}},
 'ch_w': {'description': 'Water-vapour Raman channel IDs.',
          'example': '0407xwpr',
          'legacy': {'status': '',
                     'introduced': '',
                     'old_names': [],
                     'old_location': '',
                     'old_names_removed_in': '',
                     'note': ''}},
 'ch_v': {'description': 'Reference channel IDs used with water-vapour channels.',
          'example': '0387xvpr',
          'legacy': {'status': '',
                     'introduced': '',
                     'old_names': [],
                     'old_location': '',
                     'old_names_removed_in': '',
                     'note': ''}},
 'wv_calibration_factor': {'description': 'Water-vapour calibration factor.',
                           'example': '1.0',
                           'legacy': {'status': '',
                                      'introduced': '',
                                      'old_names': [],
                                      'old_location': '',
                                      'old_names_removed_in': '',
                                      'note': ''}},
 'ch_h': {'description': 'High-rotational-Raman channel IDs for temperature '
                         'retrievals.',
          'example': '0353xhpr',
          'legacy': {'status': '',
                     'introduced': '',
                     'old_names': [],
                     'old_location': '',
                     'old_names_removed_in': '',
                     'note': ''}},
 'ch_l': {'description': 'Low-rotational-Raman channel IDs for temperature retrievals.',
          'example': '0353xlpr',
          'legacy': {'status': '',
                     'introduced': '',
                     'old_names': [],
                     'old_location': '',
                     'old_names_removed_in': '',
                     'note': ''}},
 'alpha_prime': {'description': 'Temperature calibration coefficient alpha prime.',
                 'example': '1.0',
                 'legacy': {'status': '',
                            'introduced': '',
                            'old_names': [],
                            'old_location': '',
                            'old_names_removed_in': '',
                            'note': ''}},
 'beta_prime': {'description': 'Temperature calibration coefficient beta prime.',
                'example': '1.0',
                'legacy': {'status': '',
                           'introduced': '',
                           'old_names': [],
                           'old_location': '',
                           'old_names_removed_in': '',
                           'note': ''}},
 'gamma_prime': {'description': 'Temperature calibration coefficient gamma prime.',
                 'example': '1.0',
                 'legacy': {'status': '',
                            'introduced': '',
                            'old_names': [],
                            'old_location': '',
                            'old_names_removed_in': '',
                            'note': ''}
                 }
 }
