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
    "background_correction", "background_region", "normalisation",
    "normalisation_region", "normalise_to_molecular", "plot_molecular",
    "dpi", "color_reduction",
]

SYSTEM_TEMPLATE_KEYS = ["stage_path", "reference", "label"]
CHANNEL_TEMPLATE_KEYS = [
    "label", "force_qa_test", "background_correction", "background_region",
    "normalisation", "normalisation_region", "normalise_to_molecular",
    "plot_molecular",
]
PAIR_TEMPLATE_KEYS = ["label", "force_qa_test"]


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
    "default_qa_test": _entry("General QA-test source used by channel and pair comparisons unless force_qa_test is provided in that comparison section.", "ray"),
    "default_signal_source": _entry("Fallback exported-stage parameter containing channel signals. A source specified for the reference system in a channel bundle becomes the effective default for all systems in that bundle.", "profile"),
    "default_signal_error_source": _entry("Fallback exported-stage parameter containing channel uncertainties. A source specified for the reference system in a channel bundle becomes the effective default for all systems in that bundle.", "profile_error"),
    "default_pair_source": _entry("Fallback exported-stage parameter containing channel-pair products. A source specified for the reference system in a pair bundle becomes the effective default for all systems in that bundle.", "pol_cal_ratio_mean"),
    "default_pair_error_source": _entry("Fallback exported-stage parameter containing channel-pair uncertainties. A source specified for the reference system in a pair bundle becomes the effective default for all systems in that bundle.", "pol_cal_ratio_error_mean"),
    "vertical_scale": _entry("Physical vertical coordinate used for all comparisons. height_asl is strongly recommended because systems may be located at different station altitudes.", "height_asl"),
    "vertical_method": _entry("Method used to place systems on a common physical vertical grid. interpolation uses the reference-system grid; vertical_binning creates common altitude intervals.", "interpolation"),
    "vertical_binning": _entry("Vertical bin width in kilometres. Required only when vertical_method is vertical_binning.", "0.03"),
    "vertical_min": _entry("Optional lower comparison and plotting limit in kilometres. Empty uses the common valid overlap.", "0.5"),
    "vertical_max": _entry("Optional upper comparison and plotting limit in kilometres. Empty uses the common valid overlap.", "15.0"),
    "slice_measurement": _entry("Optional temporal slices as repeating start, stop pairs. Accepted formats match call_atlas.ini: HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, or yyyymmdd_HHMMSS.", "20260801_2100, 20260802_0200"),
    "exclude_measurement": _entry("Optional temporal exclusions as repeating start, stop pairs, using the same time formats as slice_measurement.", "20260801_2330, 20260801_2345"),
    "background_correction": _entry("General switch for subtracting the mean signal in background_region. Disabled by default and overridable per channel comparison.", "False"),
    "background_region": _entry("General fallback background interval in kilometres. A channel-specific value takes priority; when empty, later metadata loading may use the reference channel's stored region.", "18.0, 22.0"),
    "normalisation": _entry("General switch for channel normalization. It may be disabled or overridden per channel comparison. Pair products are not normalized.", "True"),
    "normalisation_region": _entry("General fallback normalization interval in kilometres. A channel-specific value takes priority; when empty, later metadata loading uses the reference channel's stored Rayleigh-fit region.", "7.5, 9.0"),
    "normalise_to_molecular": _entry("If True, normalize measured channel profiles to the reference system's molecular profile using arithmetic means over the resolved normalization region.", "True"),
    "plot_molecular": _entry("If True, include the reference system's molecular profile in channel-comparison plots.", "True"),
    "dpi": _entry("Resolution of exported figures in dots per inch.", "150"),
    "color_reduction": _entry("If True, apply the ATLAS image color-reduction workflow to exported figures.", "False"),
}

SYSTEM_FLAVOR = {
    "stage_path": _entry("Absolute or relative path to one exported ATLAS stage directory. The path may point anywhere and no common parent-folder layout is required.", "../reference_system/exported/preprocessing_complete"),
    "reference": _entry("Set True for exactly one system. Its vertical grid, molecular profile, data IDs, source parameters, and channel-specific metadata defaults define the comparison reference.", "True"),
    "label": _entry("Optional display label. When empty, metadata loading will use lidar_name and finally the system section ID as fallback.", "Reference lidar"),
}

CHANNEL_FLAVOR = {
    "label": _entry("Optional comparison label. When empty, the reference system's atlas_channel_id is used.", "355 nm parallel"),
    "force_qa_test": _entry("Optional QA-test source forced for all systems in this channel comparison. When empty, default_qa_test is used.", "ray_pcb"),
    "background_correction": _entry("Enable or disable background correction for this complete channel-comparison bundle.", "True"),
    "background_region": _entry("Channel-specific background interval in kilometres, applied to every system in this comparison. Empty falls back to the general value or later reference-channel metadata.", "18.0, 22.0"),
    "normalisation": _entry("Enable or disable normalization for this complete channel-comparison bundle.", "True"),
    "normalisation_region": _entry("Channel-specific normalization interval in kilometres, applied to every system in this comparison. Empty falls back to the general value or later reference-channel metadata.", "7.5, 9.0"),
    "normalise_to_molecular": _entry("Override whether this channel bundle is normalized to the reference molecular profile.", "True"),
    "plot_molecular": _entry("Override whether the reference molecular profile is plotted for this channel bundle.", "True"),
}

PAIR_FLAVOR = {
    "label": _entry("Optional pair-comparison label. When empty, the reference system's atlas_pair_id is used.", "VLDR 355 nm"),
    "force_qa_test": _entry("Optional QA-test source forced for all systems in this pair comparison. When empty, default_qa_test is used.", "pcb"),
}

SPECIAL_FLAVOR = {
    "system_channel_id": _entry("Map this system to an atlas_channel_id. The reference system mapping is mandatory. Other systems inherit the reference ID when omitted. Use off to exclude a non-reference system from this comparison.", "0355xpgx"),
    "system_pair_id": _entry("Map this system to an atlas_pair_id. The reference system mapping is mandatory. Other systems inherit the reference ID when omitted. Use off to exclude a non-reference system from this comparison.", "vldr355_b"),
    "system_signal_source": _entry("System-specific channel signal source. When omitted, the reference system's signal_source is inherited; if the reference also omits it, default_signal_source is used.", "derived_profile"),
    "system_signal_error_source": _entry("System-specific channel uncertainty source. When omitted, the reference system's signal_error_source is inherited; if the reference also omits it, default_signal_error_source is used.", "derived_profile_error"),
    "system_pair_source": _entry("System-specific pair-product source. When omitted, the reference system's pair_source is inherited; if the reference also omits it, default_pair_source is used.", "pol_cal_ratio_mean"),
    "system_pair_error_source": _entry("System-specific pair uncertainty source. When omitted, the reference system's pair_error_source is inherited; if the reference also omits it, default_pair_error_source is used.", "pol_cal_ratio_error_mean"),
}
