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

CONFIG_TEMPLATE_SECTIONS = {'System': ['station_id',
            'lidar_name',
            'station_name',
            'lidar_id',
            'version_name',
            'version_id',
            'configuration_name',
            'configuration_id',
            'station_altitude',
            'station_latitude',
            'station_longitude',
            'zenith_angle',
            'azimuth_angle'],
 'Channels': ['recorder_channel_id',
              'scc_channel_id',
              'telescope_type',
              'channel_type',
              'channel_subtype',
              'zero_bin',
              'dead_time',
              'background_low_bin',
              'background_high_bin',
              'channel_bandwidth',
              'G',
              'H',
              'acquisition_mode',
              'emitted_wavelength',
              'detected_wavelength',
              'bins',
              'data_acquisition_range',
              'range_resolution',
              'laser_repetition_rate',
              'analog_noise_per_bin',
              'analog_noise_scaling_factor',
              'ch_n',
              'ch_f'],
 'polarization_calibration': ['ch_r', 'ch_t', 'K', 'R_to_T_transmission_ratio', 'eta'],
 'water_vapour': ['ch_w', 'ch_v', 'wv_calibration_factor'],
 'temperature': ['ch_h', 'ch_l', 'alpha_prime', 'beta_prime', 'gamma_prime']}

CONFIG_FLAVOR = {'station_id': {'description': 'ACTRIS station identifier or short station code.',
                'example': 'the'},
 'lidar_name': {'description': 'Human-readable name of the lidar system.',
                'example': 'THELISYS'},
 'station_name': {'description': 'Human-readable name of the station.',
                  'example': 'Thessaloniki'},
 'lidar_id': {'description': 'Numeric lidar identifier.', 'example': '199'},
 'version_name': {'description': 'Human-readable version name of the lidar '
                                 'configuration.',
                  'example': 'default'},
 'version_id': {'description': 'Numeric version identifier of the lidar configuration.',
                'example': '1'},
 'configuration_name': {'description': 'Human-readable name of the configuration.',
                        'example': 'standard'},
 'configuration_id': {'description': 'Numeric configuration identifier.',
                      'example': '665'},
 'station_altitude': {'description': 'Station altitude above sea level in meters.',
                      'example': '60'},
 'station_latitude': {'description': 'Station latitude in decimal degrees.',
                      'example': '40.63'},
 'station_longitude': {'description': 'Station longitude in decimal degrees.',
                       'example': '22.96'},
 'zenith_angle': {'description': 'Lidar pointing zenith angle in degrees.',
                  'example': '0'},
 'azimuth_angle': {'description': 'Lidar pointing azimuth angle in degrees.',
                   'example': '0'},
 'recorder_channel_id': {'description': 'Recorder channel IDs as they appear in the '
                                        'raw lidar files.',
                         'example': '0532xcpr, 0532xppt, 1064xtax'},
 'scc_channel_id': {'description': 'Corresponding SCC channel IDs.',
                    'example': '1893, 1895, 1700'},
 'telescope_type': {'description': 'ACTRIS telescope type code for each channel.',
                    'example': 'x, x, x'},
 'channel_type': {'description': 'ACTRIS channel type code for each channel.',
                  'example': 'c, p, t'},
 'channel_subtype': {'description': 'ACTRIS channel subtype code for each channel.',
                     'example': 'r, t, a'},
 'zero_bin': {'description': 'Bin index corresponding to the trigger/zero range for '
                             'each channel.',
              'example': '0, 0, 0'},
 'dead_time': {'description': 'Photon-counting dead time for each channel in '
                              'nanoseconds.',
               'example': '3.7, 3.7, _'},
 'background_low_bin': {'description': 'Lower bin index of the background region for '
                                       'each channel.',
                        'example': '15000, 15000, 15000'},
 'background_high_bin': {'description': 'Upper bin index of the background region for '
                                        'each channel.',
                         'example': '16000, 16000, 16000'},
 'channel_bandwidth': {'description': 'Interference-filter bandwidth for each channel.',
                       'example': '1.0, 1.0, 1.0'},
 'G': {'description': 'Polarization gain parameter G for each channel.',
       'example': '1.0, 1.0, 1.0'},
 'H': {'description': 'Polarization parameter H for each channel.',
       'example': '-1.0, 1.0, 0.0'},
 'acquisition_mode': {'description': 'Acquisition mode code for each channel.',
                      'example': 'p, p, a'},
 'emitted_wavelength': {'description': 'Emitted laser wavelength for each channel in '
                                       'nm.',
                        'example': '532.07, 532.07, 1064.14'},
 'detected_wavelength': {'description': 'Detected wavelength for each channel in nm.',
                         'example': '532.0, 532.0, 1064.0'},
 'bins': {'description': 'Number of range bins for each channel.',
          'example': '16384, 16384, 16384'},
 'data_acquisition_range': {'description': 'Data acquisition range for each channel.',
                            'example': '100000, 100000, 100000'},
 'range_resolution': {'description': 'Range resolution for each channel in meters.',
                      'example': '7.5, 7.5, 7.5'},
 'laser_repetition_rate': {'description': 'Laser repetition rate for each channel in '
                                          'Hz.',
                           'example': '20, 20, 20'},
 'analog_noise_per_bin': {'description': 'Analog noise per bin used for analog-channel '
                                         'uncertainty estimates.',
                          'example': '0.22, 0.22, 0.22'},
 'analog_noise_scaling_factor': {'description': 'Scaling factor applied to analog '
                                                'noise estimates.',
                                 'example': '0.7, 0.7, 0.7'},
 'ch_r': {'description': 'Reflected-channel IDs used in polarization calibration '
                         'channel pairs.',
          'example': '0532xcpr'},
 'ch_t': {'description': 'Transmitted-channel IDs used in polarization calibration '
                         'channel pairs.',
          'example': '0532xppt'},
 'K': {'description': 'Calibration constant K for polarization calibration.',
       'example': '1.0'},
 'R_to_T_transmission_ratio': {'description': 'Transmission ratio between reflected '
                                              'and transmitted polarization paths.',
                               'example': '1.0'},
 'eta': {'description': 'Polarization calibration eta parameter.', 'example': '1.0'},
 'ch_n': {'description': 'Near-range channel IDs used for signal gluing.',
          'example': '0532xcnr'},
 'ch_f': {'description': 'Far-range channel IDs used for signal gluing.',
          'example': '0532xcfr'},
 'ch_w': {'description': 'Water-vapour Raman channel IDs.', 'example': '0407xwpr'},
 'ch_v': {'description': 'Reference channel IDs used with water-vapour channels.',
          'example': '0387xvpr'},
 'wv_calibration_factor': {'description': 'Water-vapour calibration factor.',
                           'example': '1.0'},
 'ch_h': {'description': 'High-rotational-Raman channel IDs for temperature '
                         'retrievals.',
          'example': '0353xhpr'},
 'ch_l': {'description': 'Low-rotational-Raman channel IDs for temperature retrievals.',
          'example': '0353xlpr'},
 'alpha_prime': {'description': 'Temperature calibration coefficient alpha prime.',
                 'example': '1.0'},
 'beta_prime': {'description': 'Temperature calibration coefficient beta prime.',
                'example': '1.0'},
 'gamma_prime': {'description': 'Temperature calibration coefficient gamma prime.',
                 'example': '1.0'}}
