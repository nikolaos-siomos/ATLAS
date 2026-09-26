#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""User-facing flavor text for the generated intercomparison INI and docs."""

from __future__ import annotations

# Template inclusion and ordering are defined in template_profiles.py.
# This module contains user-facing descriptions/examples only.

def _legacy() -> dict[str, object]:
    return {
        "status": "new",
        "introduced": "1.0.0",
        "old_names": [],
        "old_location": "",
        "old_names_removed_in": "",
        "note": "",
    }


def _entry(description: str, example: str = "") -> dict[str, object]:
    return {"description": description, "example": example, "legacy": _legacy()}


GENERAL_FLAVOR = {
    "output_folder": _entry("Folder where intercomparison plots, tables, and cached products are written. Relative paths are resolved against the folder containing this INI file.", "./analysis"),
    "overwrite_output": _entry("If True, existing intercomparison outputs with the same names may be overwritten.", "False"),
    "default_qa_test": _entry("Universal fallback QA test. A dataset-specific qa_test takes priority.", "ray"),
    "default_signal_source": _entry("Universal fallback exported-stage parameter containing mean channel signals. A dataset-specific signal_source takes priority.", "profile_mean"),
    "default_signal_error_source": _entry("Universal fallback exported-stage parameter containing uncertainties of the mean channel signals. A dataset-specific signal_error_source takes priority.", "profile_error_mean"),
    "default_pair_source": _entry("Universal fallback exported-stage parameter containing channel-pair products. A dataset-specific pair_source takes priority.", "pol_cal_ratio_mean"),
    "default_pair_error_source": _entry("Universal fallback exported-stage parameter containing channel-pair uncertainties. A dataset-specific pair_error_source takes priority.", "pol_cal_ratio_error_mean"),
    "vertical_scale": _entry("Vertical coordinate used by vertical harmonization and later plotting. Allowed values are bins, range, height_agl, and height_asl. All four coordinates are loaded when available. Default: height_asl.", "height_asl"),
    "vertical_method": _entry("Vertical harmonization method. interpolation maps every participating entry to the coarsest native grid in that group so no entry is upscaled. vertical_binning applies conservative overlap-weighted binning onto a common regular grid. Default: interpolation.", "interpolation"),
    "vertical_bin_width": _entry("Default conservative vertical-bin width used when vertical_method=vertical_binning. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Default: 0.1 km for physical scales. A group-specific vertical_bin_width overrides this value. If the requested width is smaller than the coarsest nominal native step in a group, ATLAS warns and uses that coarsest native step instead.", "0.1"),
    "first_bin_left_edge": _entry("Left edge of the first bin in the common conservative grid. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Default: 0. Each output vertical-scale value is assigned to the center of its bin. Leading bins with no source overlap are retained and filled with NaN.", "0.0"),
    "vertical_min": _entry("Optional lower retained limit for the common binned grid. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Empty starts from first_bin_left_edge.", "0.5"),
    "vertical_max": _entry("Optional upper retained limit for the common binned grid. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Empty extends to the largest available vertical extent in the group.", "15.0"),
    "plot_native_scale": _entry("If True, vertical harmonization is skipped for plotting and every dataset is shown on its native selected vertical scale. The right-hand relative-difference panel is left without curves because the vertical samples are not aligned. Default: False.", "False"),
    "slice_measurement": _entry("Optional temporal slices as repeating start, stop pairs. Accepted formats match call_atlas.ini: HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, or yyyymmdd_HHMMSS.", "20260801_2100, 20260802_0200"),
    "exclude_measurement": _entry("Optional temporal exclusions as repeating start, stop pairs, using the same time formats as slice_measurement.", "20260801_2330, 20260801_2345"),
    "default_channel_background_correction": _entry("Default background-correction switch for channel groups. A value set directly in a [channel_group:<id>] section overrides this default.", "False"),
    "default_channel_background_region": _entry("Default channel background region in kilometres. A channel-group background_region overrides it.", "18.0, 22.0"),
    "default_channel_normalisation": _entry("Default normalization switch for channel groups. A channel-group normalisation value overrides it.", "True"),
    "default_channel_normalisation_region": _entry("Default channel normalization region in kilometres. A channel-group normalisation_region overrides it.", "4.0, 6.0"),
    "default_channel_normalise_to_molecular": _entry("Default channel normalization target. If True, channel groups are normalized to the reference entry's molecular profile; if False, they are normalized to the reference entry's measured signal.", "True"),
    "default_channel_plot_molecular": _entry("Default switch controlling whether the reference molecular channel profile is included in plots.", "True"),
    "default_pair_background_correction": _entry("Default background-correction switch for pair groups. A value set directly in a [pair_group:<id>] section overrides this default.", "False"),
    "default_pair_background_region": _entry("Default pair background region in kilometres. A pair-group background_region overrides it.", "18.0, 22.0"),
    "default_pair_normalisation": _entry("Default normalization switch for pair groups. Pair normalization is always to the reference entry's measured pair ratio, never to the molecular ratio.", "False"),
    "default_pair_normalisation_region": _entry("Default pair normalization region in kilometres. A pair-group normalisation_region overrides it.", "7.5, 9.0"),
    "default_pair_plot_molecular": _entry("Default switch controlling whether the reference molecular ratio is included in pair plots.", "True"),
    "dpi": _entry("Resolution of exported figures in dots per inch.", "150"),
    "color_reduction": _entry("If True, apply the ATLAS image color-reduction workflow to exported figures.", "False"),
}

PLOTTING_FLAVOR = {
    "channel_x_lims": _entry("Default horizontal limits for channel-group intercomparison plots. Empty means determine them automatically from the plotted vertical data. Units follow vertical_scale: km for physical scales, bin units for bins.", "0.0, 20.0"),
    "pair_x_lims": _entry("Default horizontal limits for pair-group intercomparison plots. Empty means determine them automatically from the plotted vertical data. Units follow vertical_scale: km for physical scales, bin units for bins.", "0.0, 10.0"),
    "x_tick": _entry("Default major horizontal-axis tick spacing for intercomparison plots. Units follow vertical_scale.", "1.0"),
    "channel_difference_y_lims": _entry("Default right-panel y-axis limits for channel-group relative differences. Group-specific difference_y_lims can override these values.", "-0.4, 0.4"),
    "pair_difference_y_lims": _entry("Default right-panel y-axis limits for pair-group absolute differences. Empty means determine them automatically from the SNR-filtered absolute differences. Group-specific difference_y_lims can override these values.", "-0.1, 0.1"),
    "channel_y_lims": _entry("Default left-panel y limits for channel-group plots. Empty means determine them automatically from all plotted channel signals and, when enabled, the reference molecular profile.", ""),
    "channel_smooth": _entry("Default smoothing switch for channel-group plots. Smoothing/local-STD estimation is used only for data that are not conservatively binned. Harmonized vertical_binning plots use the binned signal and propagated error directly.", "True"),
    "channel_smoothing_range": _entry("Default channel smoothing range. Units follow vertical_scale: km for physical scales, bin units for bins.", "0.05, 35.0"),
    "channel_smoothing_window": _entry("Default channel smoothing window. Units follow vertical_scale: km for physical scales, bin units for bins.", "0.5"),
    "pair_y_lims": _entry("Default left-panel y limits for pair-group plots. Empty means determine them automatically.", ""),
    "pair_smooth": _entry("Default smoothing switch for pair-group plots. Smoothing/local-STD estimation is used only for data that are not conservatively binned. Harmonized vertical_binning plots use the binned values and propagated errors directly.", "True"),
    "pair_smoothing_range": _entry("Default pair smoothing range. Units follow vertical_scale: km for physical scales, bin units for bins.", "0.05, 10.0"),
    "pair_smoothing_window": _entry("Default pair smoothing window. Units follow vertical_scale: km for physical scales, bin units for bins.", "0.5"),
    "near_range_max": _entry("Upper horizontal limit of the second-row near-range comparison. Units follow vertical_scale: km for physical scales, bin units for bins. Default: 2 km for physical scales.", "2.0"),
    "near_range_smooth": _entry("Default smoothing switch for the second-row near-range comparison. This smoothing is independent of the full-range row and is applied even when the harmonized data were conservatively binned.", "True"),
    "near_range_smoothing_window": _entry("Full smoothing window used only by the second-row near-range comparison. Units follow vertical_scale: km for physical scales, bin units for bins. The default 0.1 km corresponds to 100 m.", "0.1"),
}

DATASET_FLAVOR = {
    "stage_path": _entry("Absolute or relative path to one exported ATLAS stage directory. Different datasets may point to different stage paths, or multiple datasets may intentionally point to the same stage path with different QA/source selections.", "../dataset_a/exported/preprocessing_complete"),
    "reference": _entry("Set True for exactly one dataset. This defines the global reference dataset; every active comparison group must choose exactly one reference entry that points to it.", "True"),
    "system_label": _entry("Optional label for the physical lidar/system that produced this dataset. Multiple datasets may share the same system_label.", "Lidar A"),
    "dataset_label": _entry("Optional label describing this particular dataset or processing realization. When empty, later code may fall back to the dataset section ID.", "Rayleigh processing"),
    "qa_test": _entry("QA test used when reading this dataset. When empty, default_qa_test from [general] is used.", "ray"),
    "signal_source": _entry("Exported-stage parameter containing mean channel signals for this dataset. When empty, default_signal_source is used.", "profile_mean"),
    "signal_error_source": _entry("Exported-stage parameter containing uncertainties of the mean channel signals for this dataset. When empty, default_signal_error_source is used.", "profile_error_mean"),
    "pair_source": _entry("Exported-stage parameter containing pair products for this dataset. When empty, default_pair_source is used.", "pol_cal_ratio_mean"),
    "pair_error_source": _entry("Exported-stage parameter containing pair uncertainties for this dataset. When empty, default_pair_error_source is used.", "pol_cal_ratio_error_mean"),
}

CHANNEL_GROUP_FLAVOR = {
    "label": _entry("Optional channel-group display label. When empty, the reference entry's atlas_channel_id is used.", "355 nm parallel"),
    "background_correction": _entry("Optional channel-group override for background correction. Empty inherits default_channel_background_correction from [general].", "False"),
    "background_region": _entry("Optional channel-group override for the background interval in kilometres. Empty inherits default_channel_background_region from [general].", "18.0, 22.0"),
    "normalisation": _entry("Optional channel-group override for normalization. Empty inherits default_channel_normalisation from [general].", "True"),
    "normalisation_region": _entry("Optional channel-group override for the normalization interval in kilometres. Empty inherits default_channel_normalisation_region from [general].", "4.0, 6.0"),
    "normalise_to_molecular": _entry("Optional channel-group override. True normalizes every participating channel to the reference entry's molecular profile; False normalizes every channel to the reference entry's measured signal.", "True"),
    "plot_molecular": _entry("Optional channel-group override controlling whether the reference molecular profile is plotted. Empty inherits default_channel_plot_molecular from [general].", "True"),
    "vertical_bin_width": _entry("Optional channel-group conservative bin-width override. Empty uses [general] vertical_bin_width. If the requested width is smaller than the coarsest nominal native step in this group, ATLAS warns and uses the coarsest native step instead.", "0.03"),
    "smooth": _entry("Optional channel-group plotting smoothing override. Empty inherits channel_smooth from [plotting].", "True"),
    "smoothing_range": _entry("Optional channel-group smoothing-range override. Empty inherits channel_smoothing_range from [plotting].", "0.05, 35.0"),
    "smoothing_window": _entry("Optional channel-group smoothing-window override. Empty inherits channel_smoothing_window from [plotting].", "0.5"),
    "near_range_max": _entry("Optional group override for the second-row near-range upper limit. Empty inherits near_range_max from [plotting]. Units follow vertical_scale.", "2.0"),
    "near_range_smooth": _entry("Optional group override for second-row near-range smoothing. Empty inherits near_range_smooth from [plotting].", "True"),
    "near_range_smoothing_window": _entry("Optional group override for the full smoothing window used in the second-row near-range comparison. Empty inherits near_range_smoothing_window from [plotting].", "0.1"),
    "x_lims": _entry("Optional channel-group x-axis limits. Empty inherits channel_x_lims from [plotting]; an empty resolved value triggers automatic limits.", "0.0, 20.0"),
    "x_tick": _entry("Optional channel-group major x-axis tick spacing. Empty inherits x_tick from [plotting].", "2.0"),
    "y_lims": _entry("Optional channel-group left-panel y limits. Empty inherits channel_y_lims from [plotting]; an empty resolved value triggers automatic limits.", ""),
    "difference_y_lims": _entry("Optional channel-group right-panel y-axis limits for relative differences. Empty inherits channel_difference_y_lims from [plotting].", "-0.4, 0.4"),
    "use_log_y_scale": _entry("Channel-group logarithmic-y switch. Default: True. Set False in this group to use a linear left-panel y axis.", "True"),
}

PAIR_GROUP_FLAVOR = {
    "label": _entry("Optional pair-group display label. When empty, the reference entry's atlas_pair_id is used.", "VLDR 355 nm"),
    "background_correction": _entry("Optional pair-group override for background correction. Empty inherits default_pair_background_correction from [general].", "False"),
    "background_region": _entry("Optional pair-group override for the background interval in kilometres. Empty inherits default_pair_background_region from [general].", "18.0, 22.0"),
    "normalisation": _entry("Optional pair-group override for normalization. Pair products are normalized to the reference entry's measured pair ratio; molecular normalization is not used for pairs.", "False"),
    "normalisation_region": _entry("Optional pair-group override for the normalization interval in kilometres. Empty inherits default_pair_normalisation_region from [general].", "7.5, 9.0"),
    "plot_molecular": _entry("Optional pair-group override controlling whether the reference molecular ratio is plotted. Empty inherits default_pair_plot_molecular from [general].", "True"),
    "vertical_bin_width": _entry("Optional pair-group conservative bin-width override. Empty uses [general] vertical_bin_width. If the requested width is smaller than the coarsest nominal native step in this group, ATLAS warns and uses the coarsest native step instead.", "0.03"),
    "smooth": _entry("Optional pair-group plotting smoothing override. Empty inherits pair_smooth from [plotting].", "True"),
    "smoothing_range": _entry("Optional pair-group smoothing-range override. Empty inherits pair_smoothing_range from [plotting].", "0.05, 10.0"),
    "smoothing_window": _entry("Optional pair-group smoothing-window override. Empty inherits pair_smoothing_window from [plotting].", "0.5"),
    "near_range_max": _entry("Optional group override for the second-row near-range upper limit. Empty inherits near_range_max from [plotting]. Units follow vertical_scale.", "2.0"),
    "near_range_smooth": _entry("Optional group override for second-row near-range smoothing. Empty inherits near_range_smooth from [plotting].", "True"),
    "near_range_smoothing_window": _entry("Optional group override for the full smoothing window used in the second-row near-range comparison. Empty inherits near_range_smoothing_window from [plotting].", "0.1"),
    "x_lims": _entry("Optional pair-group x-axis limits. Empty inherits pair_x_lims from [plotting]; an empty resolved value triggers automatic limits.", "0.0, 10.0"),
    "x_tick": _entry("Optional pair-group major x-axis tick spacing. Empty inherits x_tick from [plotting].", "1.0"),
    "y_lims": _entry("Optional pair-group left-panel y limits. Empty inherits pair_y_lims from [plotting]; an empty resolved value triggers automatic limits.", ""),
    "difference_y_lims": _entry("Optional pair-group right-panel y-axis limits for absolute differences. Empty inherits pair_difference_y_lims from [plotting]; if both are empty, limits are calculated automatically.", "-0.1, 0.1"),
    "use_log_y_scale": _entry("Pair-group logarithmic-y switch. Default: False. Set True in this group only when a logarithmic left-panel y axis is explicitly desired.", "False"),
}

SPECIAL_FLAVOR = {
    "entry_dataset": _entry(
        "Dataset ID used by this named comparison entry. Multiple entries may point "
        "to the same dataset, which allows several channels or pairs from one exported "
        "dataset to appear in the same comparison group.",
        "dataset_a_reference",
    ),
    "entry_atlas_channel_id": _entry(
        "ATLAS channel ID selected by this named channel-group entry.",
        "0355xpgx",
    ),
    "entry_atlas_pair_id": _entry(
        "ATLAS pair ID selected by this named pair-group entry.",
        "0355UVAX",
    ),
    "entry_label": _entry(
        "Optional display label for this comparison entry. If empty, plotting falls "
        "back to the dataset label/system label plus the selected ATLAS product ID.",
        "Lidar A analog",
    ),
    "entry_reference": _entry(
        "Set True for exactly one named entry in each active group. That entry is used "
        "as the measured/molecular normalization target and as the reference for the "
        "difference panel. The reference entry must point to the global reference dataset.",
        "True",
    ),
}
