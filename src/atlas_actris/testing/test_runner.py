"""Repository test-suite CLI. Run with atlas-test or python -m atlas_actris.testing.test_runner."""
from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from typing import Optional, Sequence


def _is_repository(path: Path) -> bool:
    return ((path / "pyproject.toml").is_file()
            and (path / "src" / "atlas_actris").is_dir()
            and (path / "tests").is_dir())


def resolve_repository(requested: Optional[str] = None) -> Path:
    if requested is not None:
        path = Path(requested).expanduser().resolve()
        if _is_repository(path):
            return path
        raise ValueError(f"Not an ATLAS source checkout with tests/: {path}")
    # Prefer the checkout the user is working in, then the module's checkout
    # (the latter also works with editable installs invoked from elsewhere).
    for start in (Path.cwd().resolve(), Path(__file__).resolve().parent):
        for path in (start, *start.parents):
            if _is_repository(path):
                return path
    raise ValueError(
        "Cannot find an ATLAS source checkout with tests/. "
        "Run from the repository or use --repo /path/to/atlas_dev. "
        "The repository tests and integration dataset are not bundled in the installed package."
    )


def list_test_files(repo: Path) -> list[Path]:
    return sorted(path for path in (repo / "tests").rglob("*.py")
                  if path.name.startswith("test_") or path.name.endswith("_test.py"))


def resolve_test(repo: Path, selector: str) -> str:
    """Accept short module names, file paths, and pytest file::class::method IDs."""
    name, separator, node = selector.partition("::")
    tests = list_test_files(repo)
    requested = Path(name).expanduser()
    if requested.is_absolute():
        candidates = [path for path in tests if path.resolve() == requested.resolve()]
    elif len(requested.parts) > 1:
        candidates = [path for path in tests if path.resolve() in {
            (repo / requested).resolve(), (repo / "tests" / requested).resolve(),
        }]
    else:
        filename = name if name.endswith(".py") else name + ".py"
        if not filename.startswith("test_") and not filename.endswith("_test.py"):
            filename = "test_" + filename
        candidates = [path for path in tests if path.name == filename]
    if not candidates:
        raise ValueError(f"Unknown test: {selector!r}. Use --list to see available modules.")
    if len(candidates) != 1:
        raise ValueError(f"Ambiguous test: {selector!r}. Select its path relative to the repository.")
    result = candidates[0].relative_to(repo).as_posix()
    return result + (separator + node if separator else "")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run all ATLAS repository tests, or selected tests.")
    parser.add_argument("--repo", help="ATLAS source checkout (otherwise detected automatically).")
    parser.add_argument("-t", "--test", action="append", default=[], metavar="TEST",
                        help="Module name, test-file path or pytest node ID; repeat to select several.")
    parser.add_argument("--list", action="store_true", help="List available test modules and exit.")
    parser.add_argument("pytest_args", nargs=argparse.REMAINDER,
                        help="Additional pytest arguments after --, e.g. -- -x -vv or -- -m 'not integration'.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        repo = resolve_repository(args.repo)
        files = list_test_files(repo)
        if not files:
            raise ValueError(f"No test modules found under {repo / 'tests'}")
        if args.list:
            print(f"ATLAS test modules in {repo}:")
            for path in files:
                print("  " + path.relative_to(repo).as_posix())
            return 0
        targets = [resolve_test(repo, name) for name in args.test] if args.test else ["tests"]
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if importlib.util.find_spec("pytest") is None:
        print('ERROR: pytest is required. From the repository, run: python -m pip install -e ".[dev]"',
              file=sys.stderr)
        return 2

    extra = list(args.pytest_args)
    if extra[:1] == ["--"]:
        extra = extra[1:]
    command = [sys.executable, "-m", "pytest"]
    if (repo / "pytest.ini").is_file():
        command.extend(["-c", str(repo / "pytest.ini")])
    command.extend(targets + extra)
    env = os.environ.copy()
    # Tests use both package-qualified and legacy internal imports. Also pass
    # the source paths to integration-test children so all use this checkout.
    source_paths = [str(repo / "src"), str(repo / "src" / "atlas_actris")]
    if env.get("PYTHONPATH"):
        source_paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(source_paths)
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    print(f"ATLAS test repository: {repo}", flush=True)
    print("Selected tests: " + ", ".join(targets), flush=True)
    try:
        return subprocess.run(command, cwd=repo, env=env).returncode
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
