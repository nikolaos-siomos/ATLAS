from pathlib import Path
import os
import runpy
import sys


def main():
    # Use a non-interactive Matplotlib backend for normal terminal/batch runs.
    # This prevents TkAgg cleanup warnings in clean environments.
    #
    # Important:
    # - setdefault means users can override it, e.g.:
    #     MPLBACKEND=QtAgg atlas call_atlas.ini
    # - Spyder/debug runs that call __call_atlas_interactive__.py directly
    #   will not be affected by this CLI default.
    os.environ.setdefault("MPLBACKEND", "Agg")

    script_path = Path(__file__).with_name("__signal_viewer__.py")
    script_dir = str(script_path.parent.resolve())

    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    runpy.run_path(str(script_path), run_name="__main__")
