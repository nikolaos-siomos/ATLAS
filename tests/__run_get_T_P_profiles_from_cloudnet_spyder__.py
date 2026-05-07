#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:08:02 2025

@author: nikos
"""

import sys
from atlas_actris.__get_T_P_profiles_from_cloudnet__ import _cli as get_T_P_profiles_from_cloudnet

if __name__ == "__main__":
    sys.argv = [
        "call_atlas_from_spyder.py",   # dummy program name
        "<station name as in Cloudnet>",   # dummy program name
        "<measurement date in yyyy.mm.dd format>",
        "<measurement time in hh:mm:ss format>",
        "--save-dir",
        "<path where the netcdf files from Cloudnet are stored>",
        "--outcsv-dir",
        "<path where the radiosond-like csv files to be used in ATLAS are saved>"
    ]
    
    get_T_P_profiles_from_cloudnet()
