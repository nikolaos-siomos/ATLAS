#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 21:12:38 2026

@author: nikos
"""

from pathlib import Path
import runpy

VERSION_FILE = (
    Path(__file__).resolve().parents[1]
    / "src/atlas_actris/version.py"
)


def on_page_markdown(markdown, **kwargs):
    placeholder = "{{ atlas_version }}"

    if placeholder not in markdown:
        return markdown

    version = runpy.run_path(str(VERSION_FILE))["__version__"]
    return markdown.replace(placeholder, str(version))