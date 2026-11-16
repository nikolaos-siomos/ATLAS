from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Union
if __package__:
    from ._paths import case_output_root, remove_test_outputs
else:  # Retain direct execution of this script from its installed/source path.
    from _paths import case_output_root, remove_test_outputs


DEFAULT_CASE_NAME = "179_199_665_20231221"
DEFAULT_INI_NAME = "call_atlas_the_179_199_665_20231221.ini"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the ATLAS smoke-test dataset and check that outputs are created."
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
        help="Do not delete existing analysis/cache outputs before running.",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=1800,
        help="Timeout in seconds. Default: 1800.",
    )

    parser.add_argument(
        "--yes-clean-cache",
        default="y",
        choices=["y", "Y", "n", "N"],
        help="Answer to the first ATLAS prompt. Default: y.",
    )

    parser.add_argument(
        "--yes-export-stage",
        default="N",
        choices=["y", "Y", "n", "N"],
        help="Answer to the second ATLAS prompt. Default: N.",
    )

    return parser


def resolve_testing_pack(testing_pack: Union[str, Path]) -> Path:
    """Resolve the testing-pack directory.

    Relative paths are first resolved from the current working directory.
    When the default ``testing_pack`` path is not found there, fall back to
    the repository-level testing_pack directory located from this module.
    """
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
    yes_clean_cache: str = "y",
    yes_export_stage: str = "N",
) -> int:
    try:
        return _run_smoke_test(testing_pack, ini, case_name, keep_output, timeout,
                               yes_clean_cache, yes_export_stage)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: Smoke test stopped: {error}", file=sys.stderr)
        return 1


def _run_smoke_test(testing_pack, ini, case_name, keep_output, timeout,
                    yes_clean_cache, yes_export_stage):
    pack_dir = resolve_testing_pack(testing_pack)
    ini_path = pack_dir / ini

    if not pack_dir.exists():
        print(f"ERROR: Testing pack does not exist: {pack_dir}", file=sys.stderr)
        return 1

    if not ini_path.exists():
        print(f"ERROR: Missing INI file: {ini_path}", file=sys.stderr)
        return 1

    output_root = case_output_root(pack_dir, case_name)
    if not keep_output:
        remove_test_outputs(pack_dir, case_name, phase="old")

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

    print("Running ATLAS smoke test")
    print(f"Testing pack: {pack_dir}")
    print(f"INI file:     {ini_path}")

    env["ATLAS_CLEAN_CACHE_ANSWER"] = yes_clean_cache
    env["ATLAS_DELETE_EXPORTED_ANSWER"] = "N"
    env["ATLAS_EXPORT_STAGE_ANSWER"] = yes_export_stage

    try:
        result = subprocess.run(
            ["atlas", "-i", str(ini_path)],
            cwd=pack_dir,
            text=True,
            env=env,
            timeout=timeout,
        )

        if result.returncode != 0:
            print(f"ERROR: ATLAS failed with return code {result.returncode}", file=sys.stderr)
            return result.returncode

        output_root = case_output_root(pack_dir, case_name)

        plots = list((output_root / "plots").glob("*.png"))
        html_reports = list((output_root / "reports").glob("*.html"))
        docx_reports = list((output_root / "reports").glob("*.docx"))
        ascii_files = list((output_root / "ascii").rglob("*.txt"))

        print("\nOutput summary:")
        print(f"  Plots:        {len(plots)}")
        print(f"  HTML reports: {len(html_reports)}")
        print(f"  DOCX reports: {len(docx_reports)}")
        print(f"  ASCII files:  {len(ascii_files)}")

        if len(plots) < 10:
            print("ERROR: Too few plots were created.", file=sys.stderr)
            return 1

        if len(html_reports) < 1:
            print("ERROR: No HTML report was created.", file=sys.stderr)
            return 1

        if len(docx_reports) < 1:
            print("ERROR: No DOCX report was created.", file=sys.stderr)
            return 1

        if len(ascii_files) < 5:
            print("ERROR: Too few ASCII files were created.", file=sys.stderr)
            return 1

        print("\nATLAS smoke test completed successfully.")
        return 0
    finally:
        if not keep_output:
            remove_test_outputs(pack_dir, case_name)


def main(argv: Union[list[str], None] = None) -> int:
    args = build_parser().parse_args(argv)

    return run_smoke_test(
        testing_pack=args.testing_pack,
        ini=args.ini,
        case_name=args.case_name,
        keep_output=args.keep_output,
        timeout=args.timeout,
        yes_clean_cache=args.yes_clean_cache,
        yes_export_stage=args.yes_export_stage,
    )


if __name__ == "__main__":
    raise SystemExit(main())
