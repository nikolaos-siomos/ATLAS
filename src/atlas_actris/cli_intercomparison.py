#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Console-script wrapper for the ATLAS intercomparison module."""

from pathlib import Path
import os
import runpy
import sys


def main() -> None:
    # Match the normal ATLAS CLI behaviour. This does not affect direct Spyder
    # execution of __intercomparison_interactive__.py.
    os.environ.setdefault("MPLBACKEND", "Agg")

    script_path = Path(__file__).with_name("__intercomparison_interactive__.py")
    script_dir = str(script_path.parent.resolve())

    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    runpy.run_path(str(script_path), run_name="__main__")
