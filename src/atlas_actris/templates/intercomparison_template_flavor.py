#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""User-facing flavor text for the generated intercomparison INI and docs."""

from __future__ import annotations

GENERAL_TEMPLATE_KEYS = [
    "output_folder", "overwrite_output", "default_qa_test",
    "default_signal_source", "default_signal_error_source",
    "default_pair_source", "default_pair_error_source",
    "vertical_scale", "vertical_method", "vertical_binning",
    "vertical_min", "vertical_max", "slice_measurement", "exclude_measurement",
    "default_channel_background_correction", "default_channel_background_region",
    "default_channel_normalisation", "default_channel_normalisation_region",
    "default_channel_normalise_to_molecular", "default_channel_plot_molecular",
    "default_pair_background_correction", "default_pair_background_region",
    "default_pair_normalisation", "default_pair_normalisation_region",
    "default_pair_plot_molecular",
    "dpi", "color_reduction",
]

DATASET_TEMPLATE_KEYS = [
    "stage_path", "reference", "system_label", "dataset_label",
    "qa_test", "signal_source", "signal_error_source",
    "pair_source", "pair_error_source",
]

CHANNEL_GROUP_TEMPLATE_KEYS = [
    "label", "background_correction", "background_region",
    "normalisation", "normalisation_region", "normalise_to_molecular",
    "plot_molecular",
]

PAIR_GROUP_TEMPLATE_KEYS = [
    "label", "background_correction", "background_region",
    "normalisation", "normalisation_region", "plot_molecular",
]


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
    "default_signal_source": _entry("Universal fallback exported-stage parameter containing channel signals. A dataset-specific signal_source takes priority.", "profile"),
    "default_signal_error_source": _entry("Universal fallback exported-stage parameter containing channel uncertainties. A dataset-specific signal_error_source takes priority.", "profile_error"),
    "default_pair_source": _entry("Universal fallback exported-stage parameter containing channel-pair products. A dataset-specific pair_source takes priority.", "pol_cal_ratio_mean"),
    "default_pair_error_source": _entry("Universal fallback exported-stage parameter containing channel-pair uncertainties. A dataset-specific pair_error_source takes priority.", "pol_cal_ratio_error_mean"),
    "vertical_scale": _entry("Vertical coordinate to use later for harmonization and plotting. Allowed values are bins, range, height_agl, and height_asl. All four coordinates are loaded when available; this option only selects which one downstream processing will use. Default: height_asl.", "height_asl"),
    "vertical_method": _entry("Method used to place datasets on a common physical vertical grid. interpolation uses the reference-dataset grid; vertical_binning creates common altitude intervals.", "interpolation"),
    "vertical_binning": _entry("Vertical bin width in kilometres. Required only when vertical_method is vertical_binning.", "0.03"),
    "vertical_min": _entry("Optional lower comparison and plotting limit in kilometres. Empty uses the common valid overlap.", "0.5"),
    "vertical_max": _entry("Optional upper comparison and plotting limit in kilometres. Empty uses the common valid overlap.", "15.0"),
    "slice_measurement": _entry("Optional temporal slices as repeating start, stop pairs. Accepted formats match call_atlas.ini: HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, or yyyymmdd_HHMMSS.", "20260801_2100, 20260802_0200"),
    "exclude_measurement": _entry("Optional temporal exclusions as repeating start, stop pairs, using the same time formats as slice_measurement.", "20260801_2330, 20260801_2345"),
    "default_channel_background_correction": _entry("Default background-correction switch for channel groups. A value set directly in a [channel_group:<id>] section overrides this default.", "False"),
    "default_channel_background_region": _entry("Default channel background region in kilometres. A channel-group background_region overrides it.", "18.0, 22.0"),
    "default_channel_normalisation": _entry("Default normalization switch for channel groups. A channel-group normalisation value overrides it.", "True"),
    "default_channel_normalisation_region": _entry("Default channel normalization region in kilometres. A channel-group normalisation_region overrides it.", "7.5, 9.0"),
    "default_channel_normalise_to_molecular": _entry("Default channel normalization target. If True, channel groups are normalized to the reference dataset's molecular profile; if False, they are normalized to the reference measured signal.", "True"),
    "default_channel_plot_molecular": _entry("Default switch controlling whether the reference molecular channel profile is included in plots.", "True"),
    "default_pair_background_correction": _entry("Default background-correction switch for pair groups. A value set directly in a [pair_group:<id>] section overrides this default.", "False"),
    "default_pair_background_region": _entry("Default pair background region in kilometres. A pair-group background_region overrides it.", "18.0, 22.0"),
    "default_pair_normalisation": _entry("Default normalization switch for pair groups. Pair normalization is always to the reference measured pair ratio, never to the molecular ratio.", "False"),
    "default_pair_normalisation_region": _entry("Default pair normalization region in kilometres. A pair-group normalisation_region overrides it.", "7.5, 9.0"),
    "default_pair_plot_molecular": _entry("Default switch controlling whether the reference molecular ratio is included in pair plots.", "True"),
    "dpi": _entry("Resolution of exported figures in dots per inch.", "150"),
    "color_reduction": _entry("If True, apply the ATLAS image color-reduction workflow to exported figures.", "False"),
}

DATASET_FLAVOR = {
    "stage_path": _entry("Absolute or relative path to one exported ATLAS stage directory. Different datasets may point to different stage paths, or multiple datasets may intentionally point to the same stage path with different QA/source selections.", "../dataset_a/exported/preprocessing_complete"),
    "reference": _entry("Set True for exactly one dataset. Its vertical grid, molecular profile, product IDs, and metadata defaults define the comparison reference.", "True"),
    "system_label": _entry("Optional label for the physical lidar/system that produced this dataset. Multiple datasets may share the same system_label.", "Lidar A"),
    "dataset_label": _entry("Optional label describing this particular dataset or processing realization. When empty, later code may fall back to the dataset section ID.", "Rayleigh processing"),
    "qa_test": _entry("QA test used when reading this dataset. When empty, default_qa_test from [general] is used.", "ray"),
    "signal_source": _entry("Exported-stage parameter containing channel signals for this dataset. When empty, default_signal_source is used.", "profile"),
    "signal_error_source": _entry("Exported-stage parameter containing channel uncertainties for this dataset. When empty, default_signal_error_source is used.", "profile_error"),
    "pair_source": _entry("Exported-stage parameter containing pair products for this dataset. When empty, default_pair_source is used.", "pol_cal_ratio_mean"),
    "pair_error_source": _entry("Exported-stage parameter containing pair uncertainties for this dataset. When empty, default_pair_error_source is used.", "pol_cal_ratio_error_mean"),
}

CHANNEL_GROUP_FLAVOR = {
    "label": _entry("Optional channel-group display label. When empty, the reference dataset's atlas_channel_id is used.", "355 nm parallel"),
    "background_correction": _entry("Optional channel-group override for background correction. Empty inherits default_channel_background_correction from [general].", "False"),
    "background_region": _entry("Optional channel-group override for the background interval in kilometres. Empty inherits default_channel_background_region from [general].", "18.0, 22.0"),
    "normalisation": _entry("Optional channel-group override for normalization. Empty inherits default_channel_normalisation from [general].", "True"),
    "normalisation_region": _entry("Optional channel-group override for the normalization interval in kilometres. Empty inherits default_channel_normalisation_region from [general].", "7.5, 9.0"),
    "normalise_to_molecular": _entry("Optional channel-group override. True normalizes every participating channel to the reference dataset's molecular profile; False normalizes every channel to the reference measured signal.", "True"),
    "plot_molecular": _entry("Optional channel-group override controlling whether the reference molecular profile is plotted. Empty inherits default_channel_plot_molecular from [general].", "True"),
}

PAIR_GROUP_FLAVOR = {
    "label": _entry("Optional pair-group display label. When empty, the reference dataset's atlas_pair_id is used.", "VLDR 355 nm"),
    "background_correction": _entry("Optional pair-group override for background correction. Empty inherits default_pair_background_correction from [general].", "False"),
    "background_region": _entry("Optional pair-group override for the background interval in kilometres. Empty inherits default_pair_background_region from [general].", "18.0, 22.0"),
    "normalisation": _entry("Optional pair-group override for normalization. Pair products are normalized to the reference measured pair ratio; molecular normalization is not used for pairs.", "False"),
    "normalisation_region": _entry("Optional pair-group override for the normalization interval in kilometres. Empty inherits default_pair_normalisation_region from [general].", "7.5, 9.0"),
    "plot_molecular": _entry("Optional pair-group override controlling whether the reference molecular ratio is plotted. Empty inherits default_pair_plot_molecular from [general].", "True"),
}

SPECIAL_FLAVOR = {
    "dataset_atlas_channel_id": _entry("ATLAS channel ID selected from this dataset for the channel group. Leave empty or omit the parameter to exclude this dataset from the group. If all dataset channel IDs are empty, the group is ignored. For an active group, the reference dataset must provide an ID.", "0355xpgx"),
    "dataset_atlas_pair_id": _entry("ATLAS pair ID selected from this dataset for the pair group. Leave empty or omit the parameter to exclude this dataset from the group. If all dataset pair IDs are empty, the group is ignored. For an active group, the reference dataset must provide an ID.", "0355UVAX"),
}
