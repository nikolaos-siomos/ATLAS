#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 18:50:15 2025

@author: nikos
"""

class FileReaderError(Exception):
    """Hard raw lidar file reading errors (wrong file format / missing header information)."""

class ConfigError(Exception):
    """Hard configuration errors (invalid/missing/consistency)."""

class DataOverflowError(Exception):
    """Errors related to overflows in raw data (too many overflows, no data left)."""

class SignalError(Exception):
    """Errors related to the signals (unrecognised values, DataArray structure and related formats)."""
