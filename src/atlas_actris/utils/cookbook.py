#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 23:11:27 2026

@author: nikos
"""

recipes = dict(
    # Apply screening recipe        
    screening = [
        ("ranges_and_heights", "height_and_range_calculation"),
        ("sliced", "slice_and_exclude"),
        ("shots_screened", "screen_low_shots"),
        ("overflows_checked", "handling_overflows"),
        ("dark_assigned", "assign_dark"),
        # ("saturation_detected", "check_saturation"),
    ],
    # Apply common preprocessing recipe for all signals    
    common_preprocessing = [
        ("photon_units_converted", "photon_units_conversion"),
        ("dead_time_corrected", "dead_time_correction"),
        ("averaged", "averaging_by_time"),
        ("averaged_low_res", "averaging_by_time_low_res"),
        ("averaged_high_res", "averaging_by_time_high_res"),
        ("mean_computed", "computing_mean"),
        ],
    # Apply preprocessing recipe for all QA tests except dark    
    preprocessing = [
        ("dark_smoothed", "smoothing_dark"),
        ("dark_corrected", "dark_correction"),
        ("background_raw_calculated", "background_calculation"),
        ("background_calculated", "background_calculation"),
        ("background_corrected", "background_correction"),
        ("range_corrected", "range_correction"),
        ("vert_trimmed", "trim_vertically"),
        # ("gluing_region_found", "gluing_region"),
        # ("glued", "gluing"),
        ("molecular_calculated", "molecular_calculations"),
        ("mldr_generated", "mldr"),
        ("noise_calculated", "signal_noise_calculation"),
    ],
    # Apply preprocessing recipe for dark QA test
    dark_preprocessing = [
        ("dark_background_raw_calculated", "background_calculation"),
        ("dark_background_calculated", "background_calculation"),
        ("dark_background_corrected", "background_correction"),    
        ],
    # Apply polarisation calibration recipe for pol. cal. QA test    
    pol_cal = [
        ("gain_ratio_generated", "gain_ratio"),
        ("calibration_factor_generated", "calibration_factor"),
        ("calibrated_ratio_generated_mean", "calibrated_ratio_mean"),
        ("vldr_generated_mean", "vldr_mean"),
        ("calibrated_ratio_generated", "calibrated_ratio"),
        ("vldr_generated", "vldr"),
        ],
    )

# Checkout stages for each recipe
checkout_stages = {
    "screening": "screening_complete",
    "common_preprocessing": "common_preprocessing_complete",
    "preprocessing": "preprocessing_complete",
    "dark_preprocessing": "dark_preprocessing_complete",
    "pol_cal": "pol_cal_complete",
    }

# Checkout stages for each recipe
checkin_stages = {
    "screening": "init",
    "common_preprocessing": "screening_complete",
    "preprocessing": "common_preprocessing_complete",
    "dark_preprocessing": "common_preprocessing_complete",
    "pol_cal": "preprocessing_complete",
    }

def collect_stages():
    recipe_stages = ["init"] + [
        stage
        for recipe in recipes.values()
        for stage, _ in recipe
    ]

    return recipe_stages + list(checkout_stages.values())

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
