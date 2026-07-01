#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Beginner template: import one exported ATLAS processing stage.

Edit only the USER SETTINGS section below.

Typical use:
    1. Set output_folder to the same folder used by ATLAS during processing.
    2. Optionally set stage_name, qa_tests, and parameters.
    3. Run this script in the same Python environment where ATLAS is installed.

Notes:
    - If stage_name = None, the import utility can ask you to select a saved
      stage interactively, depending on your export_processing_stage.py version.
    - Zarr-backed xarray arrays are opened lazily. They are not fully loaded
      into RAM until you explicitly compute/load/use their values.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Optional

from atlas_actris.utils.export_processing_stage import (
    import_processing_stage,
    inspect_exported_stages,
    list_exported_stage_contents,
    list_exported_stages,
)


# =============================================================================
# USER SETTINGS
# =============================================================================

# Folder that contains the ATLAS exported/ folder.
# Example:
# output_folder = "/home/nikos/Big_Data/atlas_testing_dataset/analysis/179_199_665_20231221"
output_folder = "/home/nikos/Big_Data/atlas_testing_dataset/analysis/brc/123_202_1069_20251203"

# Name of the saved stage to import.
# Use None to select interactively if your module version supports it.
# Example: "preprocessing_complete", "pol_cal_complete", "rayleigh_fit_complete"
stage_name = None

# Optional QA-test filter.
# Use None to import all QA-test keys.
# Use a string for one key, for example: "ray"
# Use a list for several keys, for example: ["ray", "pcb_p45"]
qa_tests = None

# Optional parameter filter.
# Use None to import all parameters.
# Use a string for one parameter, for example: "profile"
# Use a list for several parameters, for example: ["profile", "profile_mean", "range"]
parameters = None


# =============================================================================
# SMALL HELPER FUNCTIONS
# =============================================================================

def _as_list_or_none(value: Optional[Any]) -> Optional[list[str]]:
    """Convert None/string/list-like values to the format expected by import."""

    if value is None:
        return None

    if isinstance(value, str):
        return [value]

    return [str(item) for item in value]


def _print_available_stages(output_folder: str) -> None:
    """Print all saved manifest-backed stages in output_folder/exported."""

    stages = list_exported_stages(output_folder)

    print("\nAvailable exported stages:")
    if len(stages) == 0:
        print("  No exported stages found.")
        return

    for index, stage in enumerate(stages, start=1):
        print(f"  {index}) {stage}")


def _print_stage_keys(contents: Dict[str, Dict[str, list[str]]]) -> None:
    """Print only non-empty QA-test keys and their parameter keys."""

    print("\nStored keys in selected stage(s):")

    if len(contents) == 0:
        print("  No stored keys found.")
        return

    for stage, qa_entries in contents.items():
        non_empty_qa = {
            qa_key: params
            for qa_key, params in qa_entries.items()
            if len(params) > 0
        }

        if len(non_empty_qa) == 0:
            continue

        print(f"\n{stage}/")

        for qa_key in sorted(non_empty_qa.keys()):
            params = sorted(non_empty_qa[qa_key])
            print(f"  {qa_key}/")
            for param in params:
                print(f"    - {param}")


def _print_imported_keys(data_pack: Dict[str, Dict[str, Any]]) -> None:
    """Print the keys that were actually imported."""

    print("\nImported data_pack keys:")

    if len(data_pack) == 0:
        print("  Nothing was imported. Check qa_tests/parameters filters.")
        return

    for qa_key in sorted(data_pack.keys()):
        params = sorted(data_pack[qa_key].keys())
        if len(params) == 0:
            continue

        print(f"  {qa_key}/")
        for param in params:
            print(f"    - {param}")


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main() -> Dict[str, Dict[str, Any]]:
    """Inspect available exports, import the selected stage, and return it."""

    # 1. Show all available stages before importing anything.
    _print_available_stages(output_folder)

    # 2. Show the available QA-test and parameter keys.
    #    If stage_name is None, inspect_exported_stages shows all saved stages.
    #    This does not load the actual data arrays.
    contents = list_exported_stage_contents(
        output_folder=output_folder,
        stage_name=stage_name,
    )
    _print_stage_keys(contents)

    # Optional: print the detailed manifest tree, including storage kinds.
    # Uncomment this block if you want more details.
    # inspect_exported_stages(
    #     output_folder=output_folder,
    #     stage_name=stage_name,
    #     include_kinds=True,
    #     include_paths=False,
    # )

    # 3. Import the selected stage.
    data_pack = import_processing_stage(
        output_folder=output_folder,
        stage_name=stage_name,
        qa_tests=_as_list_or_none(qa_tests),
        parameters=_as_list_or_none(parameters),
    )

    # 4. Print the keys that were actually imported after filters were applied.
    _print_imported_keys(data_pack)

    return data_pack


if __name__ == "__main__":
    data_pack = main()

    # -------------------------------------------------------------------------
    # EXAMPLES AFTER IMPORT
    # -------------------------------------------------------------------------
    # Access one imported entry like this:
    #
    # profile = data_pack["ray"]["profile"]
    # print(profile)
    #
    # If profile is a lazy xarray object and you really want to load it into RAM:
    # profile_loaded = profile.load()
