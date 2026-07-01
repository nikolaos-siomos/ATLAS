#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 23:11:27 2026

@author: nikos
"""

screening_recipe = [
    ("ranges_and_heights", "height_and_range_calculation"),
    ("sliced", "slice_and_exclude"),
    ("shots_screened", "screen_low_shots"),
    ("overflows_checked", "handling_overflows"),
    ("dark_assigned", "assign_dark"),
    # ("saturation_detected", "check_saturation"),
]

preprocessing_recipe = [
    ("averaged", "averaging_by_time"),
    ("averaged_low_res", "averaging_by_time_low_res"),
    ("averaged_high_res", "averaging_by_time_high_res"),
    ("mean_computed", "computing_mean"),
    ("photon_units_converted", "photon_units_conversion"),
    ("dark_smoothed", "smoothing_dark"),
    ("dark_corrected", "dark_correction"),
    ("background_raw_calculated", "background_calculation"),
    ("dead_time_corrected", "dead_time_correction"),
    ("background_calculated", "background_calculation"),
    ("background_corrected", "background_correction"),
    ("range_corrected", "range_correction"),
    ("vert_trimmed", "trim_vertically"),
    ("gluing_region_found", "gluing_region"),
    ("glued", "gluing"),
    ("molecular_calculated", "molecular_calculations"),
    ("mldr_generated", "mldr"),
    ("noise_calculated", "signal_noise_calculation"),
]

dark_rc_recipe = [
    ("averaged", "averaging_by_time"),
    ("averaged_low_res", "averaging_by_time_low_res"),
    ("averaged_high_res", "averaging_by_time_high_res"),
    ("mean_computed", "computing_mean"),
    ("photon_units_converted", "photon_units_conversion"),
    ("dark_smoothed", "smoothing_dark"),
    ("dark_corrected", "dark_correction"),
    ("background_raw_calculated", "background_calculation"),
    ("dead_time_corrected", "dead_time_correction"),
    ("background_calculated", "background_calculation"),
    ("background_corrected", "background_correction"),
    ("range_corrected", "range_correction"),
    ("vert_trimmed", "trim_vertically"),
    ("gluing_region_found", "gluing_region"),
    ("glued", "gluing"),
    ("molecular_calculated", "molecular_calculations"),
    ("mldr_generated", "mldr"),
    ("noise_calculated", "signal_noise_calculation"),
]


pol_cal_recipe = [
    ("gain_ratio_generated", "gain_ratio"),
    ("calibration_factor_generated", "calibration_factor"),
    ("calibrated_ratio_generated_mean", "calibrated_ratio_mean"),
    ("vldr_generated_mean", "vldr_mean"),
    ("calibrated_ratio_generated", "calibrated_ratio"),
    ("vldr_generated", "vldr"),
    ]

checkout_stages = {
    "init": "init",
    "screening": "screening_complete",
    "preprocessing": "preprocessing_complete",
    "pol_cal": "pol_cal_complete",
    }

def run_linear_recipe(
    processor,
    recipe,
    initial_input,
    checkout_id=None,
):
    input_id = initial_input

    for output_id, stage_name in recipe:
        processor.run(
            output_id=output_id,
            input_id=input_id,
            stage_name=stage_name,
        )

        input_id = output_id

    if checkout_id is not None:

        processor.checkout(
            output_id=checkout_id,
            input_id=input_id,
        )

        return checkout_id
