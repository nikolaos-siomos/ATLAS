#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:29:11 2022

@author: nick
"""

import glob, os

# Cleaning function
def files(folder, pattern, extension):
    # Delete all existing png files within
    files = glob.glob(os.path.join(folder, f'{pattern}.{extension}'))
    for file in files:
        os.remove(file)