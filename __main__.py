#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:48:11 2022

@author: nick
"""

import warnings, os, sys
from readers.parse_config import parse_config

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Get the full path of the config_file.ini
args = parse_config()