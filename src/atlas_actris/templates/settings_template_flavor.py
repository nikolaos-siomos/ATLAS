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

SETTINGS_TEMPLATE_SECTIONS = {'qck': 'quicklooks',
 'qck_vldr': 'quicklooks_vldr',
 'ray': 'rayleigh_fit',
 'tlc_qua': 'quadrant_telecover',
 'tlc_rin': 'ring_telecover',
 'pcb': 'polarization_calibration'}

SETTINGS_FLAVOR = {'qck': {'t_lims': {'description': 'Time limits used in the plot.',
                    'example': '2300, 0100'},
         't_tick': {'description': 'Tick spacing for the time axis.', 'example': '30'},
         'x_lims': {'description': 'Horizontal-axis limits used in the plot.',
                    'example': '0, 14'},
         'x_tick': {'description': 'Horizontal-axis tick spacing.', 'example': '1'},
         'y_lims': {'description': 'Vertical-axis limits used in the plot.',
                    'example': '0, 1'},
         'y_max_zone': {'description': 'Zone used to define or display the maximum '
                                       'signal range.',
                        'example': '0.2, 1.0'},
         'use_log_y_scale': {'description': 'Use a logarithmic vertical axis.',
                             'example': 'False'},
         'smooth': {'description': 'Apply smoothing before plotting or fitting.',
                    'example': 'True'},
         'smoothing_range': {'description': 'Height/range interval over which '
                                            'smoothing is applied.',
                             'example': '0.05, 15'},
         'smoothing_window': {'description': 'Smoothing window size.',
                              'example': '0.5'},
         'smoothing_exponential': {'description': 'Use exponential smoothing window '
                                                  'behaviour.',
                                   'example': 'False'},
         'include_channels': {'description': 'Channel IDs to include explicitly.',
                              'example': '0532xcpr, 1064xtax'},
         'exclude_wavelength': {'description': 'Wavelengths to exclude.',
                                'example': '1064'},
         'exclude_telescope_type': {'description': 'Telescope types to exclude.',
                                    'example': 'n, f'},
         'exclude_channel_type': {'description': 'Channel types to exclude.',
                                  'example': 'a, f'},
         'exclude_acquisition_mode': {'description': 'Acquisition modes to exclude.',
                                      'example': 'a'},
         'exclude_channel_subtype': {'description': 'Channel subtypes to exclude.',
                                     'example': 'w, c'}},
 'qck_vldr': {'t_lims': {'description': 'Time limits used in the plot.',
                         'example': '2300, 0100'},
              't_tick': {'description': 'Tick spacing for the time axis.',
                         'example': '30'},
              'x_lims': {'description': 'Horizontal-axis limits used in the plot.',
                         'example': '0, 14'},
              'x_tick': {'description': 'Horizontal-axis tick spacing.',
                         'example': '1'},
              'y_lims': {'description': 'Vertical-axis limits used in the plot.',
                         'example': '0, 1'},
              'y_max_zone': {'description': 'Zone used to define or display the '
                                            'maximum signal range.',
                             'example': '0.2, 1.0'},
              'use_log_y_scale': {'description': 'Use a logarithmic vertical axis.',
                                  'example': 'False'},
              'smooth': {'description': 'Apply smoothing before plotting or fitting.',
                         'example': 'True'},
              'smoothing_range': {'description': 'Height/range interval over which '
                                                 'smoothing is applied.',
                                  'example': '0.05, 15'},
              'smoothing_window': {'description': 'Smoothing window size.',
                                   'example': '0.5'},
              'smoothing_exponential': {'description': 'Use exponential smoothing '
                                                       'window behaviour.',
                                        'example': 'False'},
              'include_pairs': {'description': 'Channel/pair IDs to include '
                                               'explicitly.',
                                'example': '0532xvax, 0532xvpx'},
              'exclude_wavelength': {'description': 'Wavelengths to exclude.',
                                     'example': '1064'},
              'exclude_telescope_type': {'description': 'Telescope types to exclude.',
                                         'example': 'n, f'},
              'exclude_pair_type': {'description': 'Pair/channel types to exclude.',
                                    'example': 'a, f'},
              'exclude_acquisition_mode': {'description': 'Acquisition modes to '
                                                          'exclude.',
                                           'example': 'a'},
              'exclude_channel_subtype': {'description': 'Channel subtypes to exclude.',
                                          'example': 'w, c'}},
 'ray': {'x_lims': {'description': 'Horizontal-axis limits used in the plot.',
                    'example': '0, 14'},
         'x_tick': {'description': 'Horizontal-axis tick spacing.', 'example': '1'},
         'y_lims': {'description': 'Vertical-axis limits used in the plot.',
                    'example': '0, 1'},
         'use_lin_y_scale': {'description': 'Force a linear vertical axis.',
                             'example': 'False'},
         'normalization_region': {'description': 'Height/range interval used for '
                                                 'signal normalization.',
                                  'example': '6, 8'},
         'molecular_mask_region': {'description': 'Height/range interval searched for '
                                                  'molecular-mask candidates.',
                                   'example': '2, 34'},
         'molecular_mask_window': {'description': 'Minimum and maximum window size '
                                                  'used for molecular-mask tests.',
                                   'example': '0.5, 4'},
         'molecular_mask_window_step': {'description': 'Step size used when scanning '
                                                       'molecular-mask windows.',
                                        'example': '0.2'},
         'rsem_threshold': {'description': 'Relative standard error threshold used in '
                                           'the Rayleigh-fit mask.',
                            'example': '0.025'},
         'first_derivative_threshold': {'description': 'First-derivative threshold '
                                                       'used in the Rayleigh-fit mask.',
                                        'example': '2.0'},
         'second_derivative_threshold': {'description': 'Second-derivative threshold '
                                                        'used in the Rayleigh-fit '
                                                        'mask.',
                                         'example': '2.0'},
         'shapiro_wilk_threshold': {'description': 'Shapiro-Wilk normality threshold '
                                                   'used in the Rayleigh-fit mask.',
                                    'example': '0.05'},
         'cross_criterion_threshold': {'description': 'Cross-criterion threshold used '
                                                      'in the Rayleigh-fit mask.',
                                       'example': '1.0'},
         'durbin_watson_threshold': {'description': 'Accepted Durbin-Watson interval '
                                                    'used in the Rayleigh-fit mask.',
                                     'example': '1, 3'},
         'absolute_extinction_error_threshold': {'description': 'Maximum accepted '
                                                                'absolute extinction '
                                                                'error.',
                                                 'example': '10'},
         'isolated_point_radius': {'description': 'Neighbourhood radius used to detect '
                                                  'isolated valid mask points.',
                                   'example': '2'},
         'isolated_point_neighbour_threshold': {'description': 'Minimum number of '
                                                               'neighbours required to '
                                                               'keep a mask point.',
                                                'example': '2'},
         'smooth': {'description': 'Apply smoothing before plotting or fitting.',
                    'example': 'True'},
         'smoothing_range': {'description': 'Height/range interval over which '
                                            'smoothing is applied.',
                             'example': '0.05, 15'},
         'smoothing_window': {'description': 'Smoothing window size.',
                              'example': '0.5'},
         'include_channels': {'description': 'Channel IDs to include explicitly.',
                              'example': '0532xcpr, 1064xtax'},
         'exclude_wavelength': {'description': 'Wavelengths to exclude.',
                                'example': '1064'},
         'exclude_telescope_type': {'description': 'Telescope types to exclude.',
                                    'example': 'n, f'},
         'exclude_channel_type': {'description': 'Channel types to exclude.',
                                  'example': 'a, f'},
         'exclude_acquisition_mode': {'description': 'Acquisition modes to exclude.',
                                      'example': 'a'},
         'exclude_channel_subtype': {'description': 'Channel subtypes to exclude.',
                                     'example': 'w, c'}},
 'tlc_qua': {'plot_raw_signals': {'description': 'Plot raw telecover signals in '
                                                 'addition to processed signals.',
                                  'example': 'False'},
             'use_last_sector': {'description': 'Use the last telecover sector/ring as '
                                                'the reference sector when enabled.',
                                 'example': 'False'},
             'normalization_region': {'description': 'Height/range interval used for '
                                                     'signal normalization.',
                                      'example': '6, 8'},
             'relative_deviation_limit': {'description': 'Relative deviation limit '
                                                         'used for telecover quality '
                                                         'checks.',
                                          'example': '0.05'},
             'near_range_upper_limit': {'description': 'Upper limit of the near-range '
                                                       'interval.',
                                        'example': '2.5'},
             'smooth': {'description': 'Apply smoothing before plotting or fitting.',
                        'example': 'True'},
             'smoothing_window': {'description': 'Smoothing window size.',
                                  'example': '0.5'},
             'include_channels': {'description': 'Channel IDs to include explicitly.',
                                  'example': '0532xcpr, 1064xtax'},
             'exclude_wavelength': {'description': 'Wavelengths to exclude.',
                                    'example': '1064'},
             'exclude_telescope_type': {'description': 'Telescope types to exclude.',
                                        'example': 'n, f'},
             'exclude_channel_type': {'description': 'Channel types to exclude.',
                                      'example': 'a, f'},
             'exclude_acquisition_mode': {'description': 'Acquisition modes to '
                                                         'exclude.',
                                          'example': 'a'},
             'exclude_channel_subtype': {'description': 'Channel subtypes to exclude.',
                                         'example': 'w, c'}},
 'tlc_rin': {'plot_raw_signals': {'description': 'Plot raw telecover signals in '
                                                 'addition to processed signals.',
                                  'example': 'False'},
             'use_last_sector': {'description': 'Use the last telecover sector/ring as '
                                                'the reference sector when enabled.',
                                 'example': 'False'},
             'normalization_region': {'description': 'Height/range interval used for '
                                                     'signal normalization.',
                                      'example': '6, 8'},
             'relative_deviation_limit': {'description': 'Relative deviation limit '
                                                         'used for telecover quality '
                                                         'checks.',
                                          'example': '0.05'},
             'near_range_upper_limit': {'description': 'Upper limit of the near-range '
                                                       'interval.',
                                        'example': '2.5'},
             'smooth': {'description': 'Apply smoothing before plotting or fitting.',
                        'example': 'True'},
             'smoothing_window': {'description': 'Smoothing window size.',
                                  'example': '0.5'},
             'include_channels': {'description': 'Channel IDs to include explicitly.',
                                  'example': '0532xcpr, 1064xtax'},
             'exclude_wavelength': {'description': 'Wavelengths to exclude.',
                                    'example': '1064'},
             'exclude_telescope_type': {'description': 'Telescope types to exclude.',
                                        'example': 'n, f'},
             'exclude_channel_type': {'description': 'Channel types to exclude.',
                                      'example': 'a, f'},
             'exclude_acquisition_mode': {'description': 'Acquisition modes to '
                                                         'exclude.',
                                          'example': 'a'},
             'exclude_channel_subtype': {'description': 'Channel subtypes to exclude.',
                                         'example': 'w, c'}},
 'pcb': {'x_lims_signals': {'description': 'Horizontal-axis limits for '
                                           'polarization-calibration signal plots.',
                            'example': '0, 8'},
         'x_tick_signals': {'description': 'Horizontal-axis tick spacing for signal '
                                           'plots.',
                            'example': '1'},
         'y_lims_signals': {'description': 'Vertical-axis limits for signal plots.',
                            'example': '0, 1'},
         'x_lims_calibration': {'description': 'Horizontal-axis limits for calibration '
                                               'plots.',
                                'example': '0, 8'},
         'x_tick_calibration': {'description': 'Horizontal-axis tick spacing for '
                                               'calibration plots.',
                                'example': '1'},
         'y_lims_calibration': {'description': 'Vertical-axis limits for calibration '
                                               'plots.',
                                'example': '0, 1'},
         'x_lims_rayleigh': {'description': 'Horizontal-axis limits for Rayleigh '
                                            'comparison plots.',
                             'example': '0, 20'},
         'x_tick_rayleigh': {'description': 'Horizontal-axis tick spacing for Rayleigh '
                                            'comparison plots.',
                             'example': '2'},
         'y_lims_rayleigh': {'description': 'Vertical-axis limits for Rayleigh '
                                            'comparison plots.',
                             'example': '0, 1'},
         'calibration_region': {'description': 'Height/range interval used for '
                                               'polarization calibration.',
                                'example': '2, 4'},
         'rayleigh_region': {'description': 'Height/range interval used for Rayleigh '
                                            'reference comparison.',
                             'example': '6, 8'},
         'pldr_error_threshold': {'description': 'Allowed PLDR uncertainty threshold.',
                                  'example': '0.025'},
         'smooth': {'description': 'Apply smoothing before plotting or fitting.',
                    'example': 'True'},
         'smoothing_range': {'description': 'Height/range interval over which '
                                            'smoothing is applied.',
                             'example': '0.05, 15'},
         'smoothing_window': {'description': 'Smoothing window size.',
                              'example': '0.5'},
         'smoothing_exponential': {'description': 'Use exponential smoothing window '
                                                  'behaviour.',
                                   'example': 'False'}}}
