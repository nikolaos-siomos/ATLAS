"""Best-effort per-user GUI preferences; never write into the installation."""
import json
import os
from pathlib import Path
import sys
import tempfile

FORMAT = "atlas_gui_preferences_v1"
MAX_BYTES = 65536


def defaults():
    return {"format": FORMAT, "last_directory": None}


def preferences_path(platform=None, environ=None, home=None):
    platform = sys.platform if platform is None else platform
    environ = os.environ if environ is None else environ
    home = Path.home() if home is None else Path(home)
    if platform.startswith("win"):
        base = environ.get("LOCALAPPDATA", "")
        base = Path(base) if base and Path(base).is_absolute() else home / "AppData" / "Local"
    elif platform == "darwin":
        base = home / "Library" / "Application Support"
    else:
        base = environ.get("XDG_CONFIG_HOME", "")
        base = Path(base) if base and Path(base).is_absolute() else home / ".config"
    return base / "atlas_actris" / "gui_preferences.json"


def _clean(values):
    result = defaults()
    if not isinstance(values, dict) or values.get("format") != FORMAT:
        return result
    for key in ("last_directory",):
        value = values.get(key)
        if (isinstance(value, str) and value and not any(ord(c) < 32 for c in value)
                and Path(value).is_absolute()):
            result[key] = value
    return result


def load_preferences(path=None):
    try:
        target = preferences_path() if path is None else Path(path)
        with target.open("rb") as handle:
            data = handle.read(MAX_BYTES + 1)
        if len(data) > MAX_BYTES:
            return defaults()
        return _clean(json.loads(data.decode("utf-8-sig")))
    except (OSError, ValueError, TypeError, RuntimeError):
        return defaults()


def save_preferences(values, path=None):
    """Return False on failure; a convenience preference must not stop the GUI."""
    temporary = None
    try:
        target = preferences_path() if path is None else Path(path)
        clean = _clean(values)
        target.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=target.parent,
                                         prefix=".gui-preferences-", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(clean, handle, ensure_ascii=False, allow_nan=False)
        os.replace(temporary, target)
        return True
    except (OSError, ValueError, TypeError, RuntimeError):
        return False
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass


def remembered_directory(preferences):
    try:
        value = preferences.get("last_directory")
        path = Path(value) if value else None
        return path if path is not None and path.is_dir() and os.access(path, os.R_OK | os.X_OK) else None
    except (OSError, ValueError, TypeError, RuntimeError):
        return None
