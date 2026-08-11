#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Command-line argument parser for the ATLAS intercomparison script."""

import argparse
import os


def parse_intercomparison_args():
    """Collect the command-line arguments of the intercomparison script."""

    parser = argparse.ArgumentParser(
        description="ATLAS exported-stage intercomparison"
    )

    parser.add_argument(
        "-i",
        "--ini_file",
        metavar="ini_file",
        type=str,
        nargs="?",
        default=None,
        help="The path to the intercomparison initialization file",
    )

    args = vars(parser.parse_args())

    if args["ini_file"] is None:
        raise Exception("-- Error: The ini_file argument was not provided")

    if not os.path.exists(args["ini_file"]):
        raise Exception(
            "-- Error: The provided intercomparison initialization file "
            f"does not exist:\n{args['ini_file']}\n"
            "Please provide a valid filepath"
        )

    return args
