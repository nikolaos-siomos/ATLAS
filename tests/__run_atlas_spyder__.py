#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:08:02 2025

@author: nikos
"""

import sys
from atlas_actris.__call_atlas_interactive__ import main as atlas_main

if __name__ == "__main__":
    sys.argv = [
        "call_atlas_from_spyder.py",   # dummy program name
        "-i",
        "/home/nikos/Nextcloud5/CARS-Stations-QA/brc/call_atlas_brc_123_202_1073_20251203.ini",
    ]
    parser_args = atlas_main()
