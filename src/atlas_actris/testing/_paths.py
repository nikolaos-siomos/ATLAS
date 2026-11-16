"""Filesystem boundaries for disposable smoke-test outputs."""
import ntpath
import os
from pathlib import Path
import shutil


def _component(value):
    if (not isinstance(value, str) or not value or value in (".", "..")
            or any(c in value for c in "/\\:")
            or any(ord(c) < 32 or ord(c) == 127 for c in value)
            or value.rstrip(" .") != value or ntpath.isabs(value)):
        raise ValueError("case-name must be a single directory name, without paths or traversal")
    return value


def _child(parent, name):
    parent = Path(parent).resolve()
    path = parent / _component(name)
    if path.is_symlink() or path.resolve() != parent / name:
        raise ValueError(f"Smoke-test output must not be a symbolic link or junction: {path}")
    return path


def case_output_root(pack, case_name):
    """Resolve a case strictly beneath the selected testing pack's analysis/."""
    case_name = _component(case_name)
    return _child(_child(Path(pack).resolve(), "analysis"), case_name)


def _validate_tree(path):
    """Fail before deletion if an existing output tree contains linked entries."""
    if not path.exists():
        return
    if not path.is_dir():
        raise ValueError(f"Smoke-test output is not a directory: {path}")
    def unreadable(error):
        raise error
    for directory, dirs, files in os.walk(path, followlinks=False, onerror=unreadable):
        for name in dirs + files:
            _child(directory, name)


def remove_test_outputs(pack, case_name, children=None, phase="generated"):
    """Delete only validated case outputs, rechecking after the subprocess.

    The caller-selected pack may itself be a link. Links below that explicit
    root are rejected. Checks cover pre-existing links, not concurrent changes
    by a process with write access to the same output tree.
    """
    root = case_output_root(pack, case_name)
    targets = [root] if children is None else [_child(root, name) for name in children]
    # Validate every target before removing the first one.
    for path in targets:
        _validate_tree(path)
    for path in targets:
        if path.exists():
            label = "analysis" if children is None else path.name
            print(f"Removing {phase} {label} folder: {path}")
            shutil.rmtree(path)
    if children is None and root.parent.exists():
        # Remove analysis/ only when empty; never recursively remove other cases.
        if not any(root.parent.iterdir()):
            root.parent.rmdir()
