"""Trusted subprocess entry point for the ATLAS installed alongside this GUI.

Invoked by absolute filename with Python -I, never by module lookup from a
measurement folder. Works with wheels and editable installations on Python 3.9+.
"""
from pathlib import Path
import sys


def main():
    package = Path(__file__).resolve().parents[1]
    if not (package / "__init__.py").is_file() or not (package / "cli.py").is_file():
        print("ATLAS is missing beside its GUI. Reinstall ATLAS in this environment.",
              file=sys.stderr)
        return 1
    # Derive the import location exclusively from this installed file. In an
    # editable installation this is the checkout's src directory. The data
    # directory remains cwd for relative paths, but never becomes a code path.
    sys.path.insert(0, str(package.parent))
    from atlas_actris.cli import main as run_atlas
    return run_atlas()


if __name__ == "__main__":
    raise SystemExit(main())
