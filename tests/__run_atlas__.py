#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 01:08:02 2025

@author: nikos
"""

import sys
from atlas_actris.__call_atlas_interactive__ import main as atlas_main

if __name__ == "__main__":
    # Forward any command-line args to atlas_main
    atlas_main(sys.argv[1:])