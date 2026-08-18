#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ATLAS exported-stage intercomparison main script."""

from utils.parse_intercomparison_args import parse_intercomparison_args
from utils.parse_intercomparison_file import parse_intercomparison_ini
from utils.intercomparison_stages import collect_intercomparison_stages
from utils.intercomparison_bundles import prepare_intercomparison_bundles
from utils.intercomparison_time_processing import filter_and_average_intercomparison_bundles
from utils.intercomparison_background import apply_intercomparison_background_correction
from utils.intercomparison_normalization import apply_intercomparison_normalization


cmd_args = parse_intercomparison_args()

print("\n1) Parse intercomparison configuration")
intercomparison_info = parse_intercomparison_ini(
    filepath=cmd_args["ini_file"]
)

print("2) Initialize selective stage access")
intercomparison_stages = collect_intercomparison_stages(
    intercomparison_info=intercomparison_info
)

print("3) Load and organize requested datasets/groups")
intercomparison_bundles = prepare_intercomparison_bundles(
    intercomparison_info=intercomparison_info,
    stage_store=intercomparison_stages,
)

print("4) Apply temporal filtering and averaging")
intercomparison_bundles = filter_and_average_intercomparison_bundles(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)

print("5) Apply background correction")
intercomparison_bundles = apply_intercomparison_background_correction(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)

print("6) Apply normalization")
intercomparison_bundles = apply_intercomparison_normalization(
    intercomparison_info=intercomparison_info,
    intercomparison_bundles=intercomparison_bundles,
)

print("7) Preprocessing complete")

# Next intercomparison phases:
# 1. vertical harmonization / comparison preparation
# 2. plotting and output generation
