#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:08:02 2025

@author: nikos
"""

import sys
from atlas_actris.__get_config_file_from_scc_hoi__ import main as get_config_file_from_scc_hoi

if __name__ == "__main__":
    # Forward any command-line args to atlas_main
    get_config_file_from_scc_hoi(sys.argv[1:])

