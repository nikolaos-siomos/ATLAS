from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Union


DEFAULT_CASE_NAME = "179_199_665_20231221"
DEFAULT_INI_NAME = "call_atlas_the_179_199_665_20231221.ini"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the ATLAS signal-viewer smoke-test dataset and check that "
            "viewer outputs are created."
        )
    )

    parser.add_argument(
        "testing_pack",
        nargs="?",
        default="testing_pack",
        help=(
            "Path to the testing_pack folder. Default: ./testing_pack; "
            "if missing, use the ATLAS repository testing_pack."
        ),
    )

    parser.add_argument(
        "--ini",
        default=DEFAULT_INI_NAME,
        help=f"INI file inside testing_pack. Default: {DEFAULT_INI_NAME}",
    )

    parser.add_argument(
        "--case-name",
        default=DEFAULT_CASE_NAME,
        help=f"Expected case output folder under analysis/. Default: {DEFAULT_CASE_NAME}",
    )

    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Do not delete existing signal_viewer/cache outputs before running.",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=1800,
        help="Timeout in seconds. Default: 1800.",
    )

    return parser


def resolve_testing_pack(testing_pack: Union[str, Path]) -> Path:
    """Resolve the testing-pack directory."""
    requested_path = Path(testing_pack).expanduser()
    cwd_path = requested_path.resolve()

    if cwd_path.exists():
        return cwd_path

    if requested_path == Path("testing_pack"):
        repository_path = (
            Path(__file__).resolve().parents[3] / "testing_pack"
        ).resolve()

        if repository_path.exists():
            return repository_path

    return cwd_path


def run_smoke_test(
    testing_pack: Union[str, Path] = "testing_pack",
    ini: str = DEFAULT_INI_NAME,
    case_name: str = DEFAULT_CASE_NAME,
    keep_output: bool = False,
    timeout: int = 1800,
) -> int:
    pack_dir = resolve_testing_pack(testing_pack)
    ini_path = pack_dir / ini

    if not pack_dir.exists():
        print(f"ERROR: Testing pack does not exist: {pack_dir}", file=sys.stderr)
        return 1

    if not ini_path.exists():
        print(f"ERROR: Missing INI file: {ini_path}", file=sys.stderr)
        return 1

    output_root = pack_dir / "analysis" / case_name
    viewer_dir = output_root / "signal_viewer"
    cache_dir = output_root / "cache"

    if not keep_output:
        if viewer_dir.exists():
            print(f"Removing old signal_viewer folder: {viewer_dir}")
            shutil.rmtree(viewer_dir)

        if cache_dir.exists():
            print(f"Removing old cache folder: {cache_dir}")
            shutil.rmtree(cache_dir)

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    env["ATLAS_CLEAN_CACHE_ANSWER"] = "y"
    env["ATLAS_CLEAN_VIEWER_ANSWER"] = "N"

    print("Running ATLAS signal-viewer smoke test")
    print(f"Testing pack: {pack_dir}")
    print(f"INI file:     {ini_path}")

    try:
        result = subprocess.run(
            ["atlas-signal-viewer", "-i", str(ini_path)],
            cwd=pack_dir,
            text=True,
            env=env,
            timeout=timeout,
        )

        if result.returncode != 0:
            print(
                f"ERROR: ATLAS signal viewer failed with return code {result.returncode}",
                file=sys.stderr,
            )
            return result.returncode

        if not viewer_dir.is_dir():
            print(
                f"ERROR: Signal-viewer output folder was not created: {viewer_dir}",
                file=sys.stderr,
            )
            return 1

        output_files = [path for path in viewer_dir.rglob("*") if path.is_file()]

        print("\nOutput summary:")
        print(f"  Signal-viewer files: {len(output_files)}")

        if len(output_files) == 0:
            print("ERROR: No signal-viewer output files were created.", file=sys.stderr)
            return 1

        print("\nATLAS signal-viewer smoke test completed successfully.")
        return 0
    finally:
        if not keep_output:
            if viewer_dir.exists():
                print(f"Removing generated signal_viewer folder: {viewer_dir}")
                shutil.rmtree(viewer_dir)

            if cache_dir.exists():
                print(f"Removing generated cache folder: {cache_dir}")
                shutil.rmtree(cache_dir)


def main(argv: Union[list[str], None] = None) -> int:
    args = build_parser().parse_args(argv)

    return run_smoke_test(
        testing_pack=args.testing_pack,
        ini=args.ini,
        case_name=args.case_name,
        keep_output=args.keep_output,
        timeout=args.timeout,
    )


if __name__ == "__main__":
    raise SystemExit(main())
