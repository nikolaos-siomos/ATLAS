#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Template profile configuration for generated ATLAS INI files.

The full templates are controlled by the section dictionaries in the flavor
files.  This file only controls reduced profiles, especially the beginner
profile.  Keep the lists short and practical; the generator validates all names
against the parser schemas.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Beginner initialization template
# ---------------------------------------------------------------------------
# Section names must match the strict lowercase initialization parser sections.
# Keys must exist in parse_init_file.SCHEMA.

BEGINNER_INIT_TEMPLATE_SECTIONS = {
    "configuration": [
        "scc_configuration_id",
        "export_hoi_cfg",
    ],
    "explicit_paths": [
        "parent_folder",
        "atlas_configuration_file",
        "atlas_settings_file",
        "radiosonde_folder",
        "radiosonde_file",
        "output_folder",
    ],
    "general_options": [
        "process",
        "process_qck",
        "process_bgd",
        "process_vldr",
        "process_dedicated_dark",
    ],
    "trimming_options": [
        "max_height_agl",
        "slice_measurement",
        "exclude_measurement",
    ],
    "parsing_options": [
        "rsonde_station_name",
        "rsonde_station_wmo_id",
        "cloudnet_station_name",
    ],
}


# ---------------------------------------------------------------------------
# Beginner configuration template
# ---------------------------------------------------------------------------
# Section names must match the user-facing configuration INI sections.
# Keys must exist in parse_config_file.SCHEMA.

BEGINNER_CONFIG_TEMPLATE_SECTIONS = {
    "System": [
        "station_id",
        "lidar_name",
        "station_name",
        "station_altitude",
        "station_latitude",
        "station_longitude",
        "zenith_angle",
        "azimuth_angle",
    ],
    "Channels": [
        "recorder_channel_id",
        "scc_channel_id",
        "telescope_type",
        "channel_type",
        "channel_subtype",
        "zero_bin",
        "dead_time",
        "background_low_bin",
        "background_high_bin",
        "acquisition_mode",
        "detected_wavelength",
        "emitted_wavelength",
        "channel_bandwidth",
        "bins",
        "range_resolution",
        "laser_repetition_rate",
        "G",
        "H",
    ],
    "polarization_calibration": [
        "ch_r",
        "ch_t",
        "K",
        "R_to_T_transmission_ratio",
        "eta",
    ],
    "gluing": [
        "ch_n",
        "ch_f",
    ],
    "water_vapour": [
        "ch_w",
        "ch_v",
        "wv_calibration_factor",
    ],
    "temperature": [
        "ch_h",
        "ch_l",
        "alpha_prime",
        "beta_prime",
        "gamma_prime",
    ],
}


# ---------------------------------------------------------------------------
# Beginner settings template
# ---------------------------------------------------------------------------
# Top-level keys must match parse_settings_file.SCHEMA groups.
# Values are lists of keys inside each settings schema group.

BEGINNER_SETTINGS_TEMPLATE_KEYS = {
    "qck": [
        "x_lims",
        "x_tick",
        "y_lims",
        "use_log_y_scale",
        "smooth",
        "smoothing_window",
        "select_channels",
        "exclude_wavelength",
    ],
    "qck_vldr": [
        "x_lims",
        "x_tick",
        "y_lims",
        "use_log_y_scale",
        "smooth",
        "smoothing_window",
        "include_pairs",
        "exclude_wavelength",
    ],
    "ray": [
        "x_lims",
        "x_tick",
        "y_lims",
        "normalization_region",
        "molecular_mask_region",
        "smooth",
        "smoothing_window",
        "select_channels",
        "exclude_wavelength",
    ],
    "tlc": [
        "plot_raw_signals",
        "use_last_sector",
        "normalization_region",
        "relative_deviation_limit",
        "near_range_upper_limit",
        "smooth",
        "smoothing_window",
        "select_channels",
        "exclude_wavelength",
    ],
    "tlc_rin": [
        "plot_raw_signals",
        "use_last_sector",
        "normalization_region",
        "relative_deviation_limit",
        "near_range_upper_limit",
        "smooth",
        "smoothing_window",
        "select_channels",
        "exclude_wavelength",
    ],
    "pcb": [
        "x_lims_signals",
        "x_lims_calibration",
        "x_lims_rayleigh",
        "calibration_region",
        "rayleigh_region",
        "pldr_error_threshold",
        "smooth",
        "smoothing_window",
    ],
}
