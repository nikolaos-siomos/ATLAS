#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ATLAS exported-stage intercomparison main script."""

from utils.parse_intercomparison_args import parse_intercomparison_args
from utils.parse_intercomparison_file import parse_intercomparison_ini


# Get the input INI file path of the intercomparison caller.
cmd_args = parse_intercomparison_args()

# Parse and validate the intercomparison initialization file.
# Optional values omitted from the INI are filled from the parser schemas.
intercomparison_info = parse_intercomparison_ini(
    filepath=cmd_args["ini_file"]
)


# Continue the intercomparison workflow directly below this point.
