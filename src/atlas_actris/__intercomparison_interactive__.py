#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ATLAS exported-stage intercomparison main script."""

from intercomparison.arguments import parse_intercomparison_args
from intercomparison.config import parse_intercomparison_ini
from intercomparison.stages import collect_intercomparison_stages
from intercomparison.bundles import prepare_intercomparison_bundles
from intercomparison.time_processing import filter_and_average_intercomparison_bundles
from intercomparison.background import apply_intercomparison_background_correction
from intercomparison.normalization import apply_intercomparison_normalization
from intercomparison.vertical_processing import harmonize_intercomparison_vertical
from visualizer.generate_intercomparison import generate_intercomparison


print("-----------------------------------------------")
print("1) Parsing intercomparison arguments")
print("-----------------------------------------------")
cmd_args = parse_intercomparison_args()

print("-----------------------------------------------")
print("2) Parsing intercomparison initialization file")
print("-----------------------------------------------")
intercomparison_info = parse_intercomparison_ini(filepath=cmd_args["ini_file"])

print("-----------------------------------------------")
print("3) Collecting exported processing stages")
print("-----------------------------------------------")
stage_store = collect_intercomparison_stages(intercomparison_info)

print("-----------------------------------------------")
print("4) Preparing intercomparison bundles")
print("-----------------------------------------------")
intercomparison_bundles = prepare_intercomparison_bundles(
    intercomparison_info=intercomparison_info,
    stage_store=stage_store,
)

print("-----------------------------------------------")
print("5) Time filtering and averaging")
print("-----------------------------------------------")
intercomparison_bundles = filter_and_average_intercomparison_bundles(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)

print("-----------------------------------------------")
print("6) Background correction and normalization")
print("-----------------------------------------------")
intercomparison_bundles = apply_intercomparison_background_correction(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)
intercomparison_bundles = apply_intercomparison_normalization(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)

print("-----------------------------------------------")
print("7) Vertical harmonization")
print("-----------------------------------------------")
intercomparison_bundles = harmonize_intercomparison_vertical(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)

print("-----------------------------------------------")
print("8) Generate intercomparison plots")
print("-----------------------------------------------")
intercomparison_plots = generate_intercomparison(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)

print("-----------------------------------------------")
print("9) Intercomparison complete")
print("-----------------------------------------------")
