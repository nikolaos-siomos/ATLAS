from pathlib import Path
import runpy
import sys

def main():
    script_path = Path(__file__).with_name("__call_atlas_interactive__.py")
    script_dir = str(script_path.parent.resolve())

    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    runpy.run_path(str(script_path), run_name="__main__")
