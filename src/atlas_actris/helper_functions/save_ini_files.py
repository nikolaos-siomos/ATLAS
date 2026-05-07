#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 20:34:29 2025

@author: nikos
"""

from pathlib import Path
import shutil
from typing import Dict

def copy_ini_files(parser_args: Dict, path_call_atlas_ini: str,
                   new_extension: str = ".backup"):
    """
    Copy .ini files to dest_folder, replacing the extension.
    Example: 'file.ini' -> 'file.backup'
    """
    
    print('-----------------------------------------')
    print('Saving ini files to parent folder...')
    print('-----------------------------------------')

    dest_folder = Path(parser_args['output_folder'])
    dest_folder.mkdir(parents=True, exist_ok=True)

    src_files = [
        path_call_atlas_ini,
        parser_args['atlas_configuration_file'],
        parser_args['atlas_settings_file'],
        ]
    
    # Ensure extension starts with a dot
    if not new_extension.startswith("."):
        new_extension = "." + new_extension

    for src in src_files:
        src = Path(src)

        if not src.is_file():
            print(f"WARNING: Source file not found: {src}")
            continue

        # Replace extension
        dest = dest_folder / (src.stem + new_extension)

        shutil.copy2(src, dest)
        print(f"Copied: {src} → {dest}\n")
        