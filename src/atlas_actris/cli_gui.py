"""Lightweight GUI entry point; --help works even without Tkinter installed."""
import argparse
import json
import os
from pathlib import Path
import sys


def configure_input_method(system_input_method=False):
    """Avoid the Tk/XIM slowdown in this Linux GUI process only.

    Tk uses X11 on Linux, including through Xwayland. Users who need desktop
    input-method composition can explicitly retain their existing configuration.
    This must run before creating the first Tk window.
    """
    if not system_input_method and sys.platform.startswith("linux") and os.environ.get("DISPLAY"):
        os.environ["XMODIFIERS"] = "@im=none"


def get_font_info(root):
    """Report renderer/build facts without changing font or OS settings."""
    from tkinter import TclError, font
    try:
        renderer = str(root.tk.call("tk::pkgconfig", "get", "fontsystem"))
    except TclError:
        renderer = "not reported by this Tk build"
    system = str(root.tk.call("tk", "windowingsystem"))
    result = {
        "python": sys.version.split()[0],
        "tk_version": str(root.tk.call("package", "present", "Tk")),
        "windowing_system": system,
        "font_renderer": renderer,
        "tk_scaling": float(root.tk.call("tk", "scaling")),
        "default_font": font.Font(root=root, name="TkDefaultFont", exists=True).actual(),
        "font_snapshot": "Tk defaults before ATLAS GUI size adjustments",
        "xmodifiers": os.environ.get("XMODIFIERS", ""),
    }
    builds = []
    for name in ("tk", "libtk"):
        for path in (Path(sys.prefix) / "conda-meta").glob(name + "-*.json"):
            try:
                metadata = json.loads(path.read_text(encoding="utf-8"))
                if not isinstance(metadata, dict):
                    continue
                builds.append({key: metadata.get(key) for key in ("name", "version", "build")})
            except (OSError, ValueError):
                continue
    if builds:
        result["conda_tk_builds"] = builds
    if system == "x11":
        result["note"] = (
            "Xft is available. Actual antialiasing also depends on fontconfig/display settings."
            if renderer == "xft" else
            "Xft was not confirmed. Check the Tk build; a noxft build cannot gain Xft antialiasing through Python font settings."
        )
    return result


def print_font_info():
    import tkinter as tk
    try:
        root = tk.Tk()
    except tk.TclError as error:
        print("Cannot inspect font rendering: {}. A working display is required.".format(error), file=sys.stderr)
        return 1
    try:
        root.withdraw()
        print(json.dumps(get_font_info(root), indent=2))
    finally:
        root.destroy()
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="Open the ATLAS INI editor and processing GUI.")
    parser.add_argument("--font-info", action="store_true", help="Report Tk font renderer/build details and exit.")
    parser.add_argument("--system-input-method", action="store_true",
        help="Keep the desktop input method on Linux instead of bypassing XIM (may be slower).")
    args = parser.parse_args(argv)
    configure_input_method(args.system_input_method)
    try:
        if args.font_info:
            return print_font_info()
        from .GUI.atlas_schema_gui import main as launch
    except ModuleNotFoundError as error:
        if error.name not in ("tkinter", "_tkinter"):
            raise
        print("Tkinter is unavailable in this Python environment. For Conda, install tk "
              "in that environment; for system Python, install its Tkinter package. "
              "Check with: python -m tkinter", file=sys.stderr)
        return 1
    return launch()


if __name__ == "__main__":
    raise SystemExit(main())
