#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:08:02 2025

@author: nikos
"""

import sys
from atlas_actris.__get_config_file_from_scc_hoi__ import main as get_config_file_from_scc_hoi

if __name__ == "__main__":
    sys.argv = [
        "call_atlas_from_spyder.py",   # dummy program name
        "-i",
        "<SCC_config_ID>",
        "-c",
        "<path to SCC configurations folder>",
        "-o",
        "<path to SCC HOI folder>"
    ]
    
    get_config_file_from_scc_hoi()
