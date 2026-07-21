#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 17:58:13 2026

@author: nikos
"""

import os
import shutil

def ask_clean_cache(caller_info):
    cache_path = os.path.join(caller_info["output_folder"], "cache")

    if not os.path.exists(cache_path):
        return

    answer = os.environ.get("ATLAS_CLEAN_CACHE_ANSWER")

    if answer is None:
        answer = input(f"\nDelete temporary cache folder?\n{cache_path}\n[y/N]: ")

    answer = answer.strip().lower()

    if answer in ["y", "yes"]:
        shutil.rmtree(cache_path)
        print(f"-- Deleted cache folder: {cache_path}")
        print()
    else:
        print(f"-- Kept cache folder: {cache_path}")
        print()

def ask_clean_viewer(caller_info):
    signal_viewer_path = os.path.join(caller_info["output_folder"], "signal_viewer")

    if not os.path.exists(signal_viewer_path):
        return

    answer = os.environ.get("ATLAS_CLEAN_VIEWER_ANSWER")

    if answer is None:
        answer = input(f"\nDelete signal_viewer folder?\n{signal_viewer_path}\n[y/N]: ")

    answer = answer.strip().lower()

    if answer in ["y", "yes"]:
        shutil.rmtree(signal_viewer_path)
        print(f"-- Deleted signal_viewer folder: {signal_viewer_path}")
        print()
    else:
        print(f"-- Kept signal_viewer folder: {signal_viewer_path}")
        print()
