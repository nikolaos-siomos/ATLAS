#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities to export and import ATLAS processor stages.

The exported layout is:

    output_folder/exported/<stage_name>/
        manifest.pkl
        <qa_test>/<parameter>/data.zarr      # numeric xarray objects
        <qa_test>/<parameter>/data.pkl       # object/mixed metadata or non-xarray objects

The design intentionally writes each qa_test/parameter entry separately instead
of merging all entries into one xarray.Dataset. This avoids conflicts when two
entries reuse the same dimension name with different sizes, for example
``profile`` and ``profile_mean`` both having a ``time`` dimension.
"""

from __future__ import annotations

import os
import re
import pickle
import shutil
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Optional, Union

import numpy as np
import xarray as xr

_DATAARRAY_VAR = "__dataarray__"
_MANIFEST_NAME = "manifest.pkl"


def _safe_name(name: Any) -> str:
    """Return a filesystem-safe representation of a dictionary key."""

    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_")
    return safe or "unnamed"


def _export_root(output_folder: str) -> str:
    """Return output_folder/exported, creating it if needed."""

    if output_folder is None:
        raise ValueError("output_folder is required for stage export/import")

    root = os.path.join(output_folder, "exported")
    os.makedirs(root, exist_ok=True)
    return root


def _stage_dir(output_folder: str, stage_name: str) -> str:
    """Return the export folder of a specific stage."""

    return os.path.join(_export_root(output_folder), _safe_name(stage_name))


def _ask_yes_no(question: str, default: bool = False) -> bool:
    """Ask a terminal yes/no question and return the answer.

    Parameters
    ----------
    question
        Prompt text shown to the user.
    default
        Answer used when the user only presses Enter.
    """

    suffix = "[Y/n]" if default else "[y/N]"

    while True:
        answer = input(f"{question} {suffix} ").strip().lower()

        if answer == "":
            return default

        if answer in ["y", "yes"]:
            return True

        if answer in ["n", "no"]:
            return False

        print("Please answer 'y' or 'n'.")



def _print_numbered_options(title: str, options: Iterable[str]) -> list[str]:
    """Print options with ascending 1-based numbers and return them as a list."""

    option_list = [str(option) for option in options]

    if len(option_list) == 0:
        raise ValueError(f"No options available for selection: {title}")

    print(title)
    for idx, option in enumerate(option_list, start=1):
        print(f"  {idx}) {option}")

    return option_list


def _parse_number_selection(answer: str, n_options: int, allow_multiple: bool) -> list[int]:
    """Parse a numbered user selection.

    Accepts values such as ``1``, ``1,3``, ``1 3``, ``1-3`` and ``all``.
    Returned indices are zero-based and unique while preserving selection order.
    """

    text = answer.strip().lower()

    if text in ["all", "a", "*"]:
        if not allow_multiple:
            raise ValueError("Please select only one number.")
        return list(range(n_options))

    if text == "":
        raise ValueError("No selection was provided.")

    tokens = re.split(r"[ ,;]+", text)
    selected: list[int] = []

    for token in tokens:
        if token == "":
            continue

        if "-" in token:
            if not allow_multiple:
                raise ValueError("Ranges are only allowed for multiple selection.")

            parts = token.split("-", 1)
            if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isdigit():
                raise ValueError(f"Invalid range selection: {token!r}")

            start = int(parts[0])
            stop = int(parts[1])

            if start > stop:
                raise ValueError(f"Invalid descending range selection: {token!r}")

            numbers = range(start, stop + 1)
        else:
            if not token.isdigit():
                raise ValueError(f"Invalid selection: {token!r}")
            numbers = [int(token)]

        for number in numbers:
            if number < 1 or number > n_options:
                raise ValueError(
                    f"Selection {number} is outside the valid range 1-{n_options}"
                )

            idx = number - 1
            if idx not in selected:
                selected.append(idx)

    if len(selected) == 0:
        raise ValueError("No valid selection was provided.")

    if not allow_multiple and len(selected) != 1:
        raise ValueError("Please select only one number.")

    return selected


def _select_from_options(
    title: str,
    options: Iterable[str],
    allow_multiple: bool = False,
    prompt: Optional[str] = None,
) -> Union[str, list[str]]:
    """Prompt the user to select one or more values from numbered options."""

    option_list = _print_numbered_options(title, options)

    if prompt is None:
        if allow_multiple:
            prompt = "Select number(s), for example 1,3 or 1-3 or all: "
        else:
            prompt = "Select number: "

    while True:
        answer = input(prompt)

        try:
            indices = _parse_number_selection(
                answer=answer,
                n_options=len(option_list),
                allow_multiple=allow_multiple,
            )
        except ValueError as exc:
            print(f"-- {exc}")
            continue

        selected = [option_list[idx] for idx in indices]

        if allow_multiple:
            return selected

        return selected[0]


def select_exported_stage(output_folder: str) -> str:
    """Prompt the user to select one available exported stage by number."""

    return str(
        _select_from_options(
            title="-- Available exported processing stages:",
            options=list_exported_stages(output_folder),
            allow_multiple=False,
        )
    )


def select_exported_stages(output_folder: str) -> list[str]:
    """Prompt the user to select one or more available exported stages by number."""

    return list(
        _select_from_options(
            title="-- Available exported processing stages:",
            options=list_exported_stages(output_folder),
            allow_multiple=True,
        )
    )


def select_stage_names(stage_names: Iterable[str], allow_multiple: bool = True) -> Union[str, list[str]]:
    """Prompt the user to select one or more stage names from a provided list."""

    return _select_from_options(
        title="-- Available processing stages:",
        options=sorted(str(stage_name) for stage_name in stage_names),
        allow_multiple=allow_multiple,
    )


def _available_processor_stage_names(processor: Any) -> list[str]:
    """Best-effort discovery of stage names stored on an ATLAS Processor."""

    candidates: list[str] = []

    for attr_name in [
        "stages",
        "processing_stages",
        "stage_data",
        "stages_data",
        "data",
        "_stages",
    ]:
        value = getattr(processor, attr_name, None)
        if isinstance(value, dict):
            candidates.extend(str(key) for key in value.keys())

    processing_info = getattr(processor, "processing_info", None)
    if isinstance(processing_info, dict):
        for key in ["stages", "processing_stages", "stage_data", "stage_names"]:
            value = processing_info.get(key)
            if isinstance(value, dict):
                candidates.extend(str(name) for name in value.keys())
            elif isinstance(value, (list, tuple, set)):
                candidates.extend(str(name) for name in value)

    # Keep order stable but remove duplicates and obvious non-stage placeholders.
    seen: set[str] = set()
    out: list[str] = []
    for name in candidates:
        if name in seen:
            continue
        seen.add(name)
        out.append(name)

    return sorted(out)


def _resolve_processor_stage_name(processor: Any, stage_name: Optional[str]) -> str:
    """Return a processor stage name, prompting when stage_name is None."""

    if stage_name is not None:
        return str(stage_name)

    available = _available_processor_stage_names(processor)

    if len(available) == 0:
        raise ValueError(
            "stage_name was not provided and available processor stages could "
            "not be discovered automatically. Please pass stage_name explicitly."
        )

    return str(select_stage_names(available, allow_multiple=False))


def _resolve_processor_stage_names(
    processor: Any,
    stage_names: Optional[Iterable[str]],
) -> list[str]:
    """Return processor stage names, prompting when stage_names is None."""

    if stage_names is not None:
        return [str(stage_name) for stage_name in stage_names]

    available = _available_processor_stage_names(processor)

    if len(available) == 0:
        raise ValueError(
            "stage_names was not provided and available processor stages could "
            "not be discovered automatically. Please pass stage_names explicitly."
        )

    return list(select_stage_names(available, allow_multiple=True))


def _has_object_dtype_xarray(obj: Union[xr.DataArray, xr.Dataset]) -> bool:
    """Return True when an xarray object contains mixed/object data."""

    if isinstance(obj, xr.DataArray):
        return obj.dtype == object

    for var in obj.data_vars.values():
        if var.dtype == object:
            return True

    return False


def _write_pickle(obj: Any, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def _read_pickle(path: str) -> Any:
    with open(path, "rb") as f:
        return pickle.load(f)


def _human_readable_size(nbytes: Union[int, float]) -> str:
    """Return a compact human-readable byte size."""

    try:
        size = float(nbytes)
    except Exception:
        return "unknown"

    if not np.isfinite(size) or size < 0:
        return "unknown"

    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    unit = units[0]

    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            break
        size /= 1024.0

    if unit == "B":
        return f"{int(round(size))} {unit}"

    return f"{size:.2f} {unit}"


def _xarray_nbytes(obj: Union[xr.DataArray, xr.Dataset]) -> int:
    """Estimate xarray payload size without computing Dask-backed arrays."""

    total = 0

    if isinstance(obj, xr.DataArray):
        total += int(getattr(obj, "nbytes", 0))
        for coord in obj.coords.values():
            total += int(getattr(coord, "nbytes", 0))
        return total

    for var in obj.data_vars.values():
        total += int(getattr(var, "nbytes", 0))

    for coord in obj.coords.values():
        total += int(getattr(coord, "nbytes", 0))

    return total


def _pickle_nbytes(obj: Any) -> int:
    """Estimate pickle size by serializing to bytes.

    This is intended mainly for small metadata objects. If serialization fails,
    return zero rather than blocking the export.
    """

    try:
        return len(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))
    except Exception:
        return 0


def estimate_processing_stage_size(
    stage_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """Estimate the exported size of a processing stage.

    Notes
    -----
    The value is an estimate, not an exact final disk usage. Zarr metadata,
    compression, filesystem block size, and consolidated metadata can make the
    real size slightly different. Dask-backed xarray arrays are not computed.
    """

    if not isinstance(stage_data, dict):
        raise TypeError("stage_data must be a dictionary")

    zarr_bytes = 0
    pickle_bytes = 0
    entries = 0
    xarray_entries = 0
    pickle_entries = 0

    for qa_test, qa_dict in stage_data.items():
        if not isinstance(qa_dict, dict):
            continue

        for parameter, value in qa_dict.items():
            entries += 1

            if isinstance(value, (xr.DataArray, xr.Dataset)) and not _has_object_dtype_xarray(value):
                zarr_bytes += _xarray_nbytes(value)
                xarray_entries += 1
            elif isinstance(value, (xr.DataArray, xr.Dataset)):
                # Object/mixed xarray metadata are written with pickle. Use the
                # xarray payload as a cheap lower-bound estimate and avoid any
                # accidental expensive computation.
                pickle_bytes += _xarray_nbytes(value)
                pickle_entries += 1
            else:
                pickle_bytes += _pickle_nbytes(value)
                pickle_entries += 1

    total_bytes = zarr_bytes + pickle_bytes

    return {
        "total_bytes": int(total_bytes),
        "zarr_bytes": int(zarr_bytes),
        "pickle_bytes": int(pickle_bytes),
        "entries": int(entries),
        "xarray_entries": int(xarray_entries),
        "pickle_entries": int(pickle_entries),
        "total_human": _human_readable_size(total_bytes),
        "zarr_human": _human_readable_size(zarr_bytes),
        "pickle_human": _human_readable_size(pickle_bytes),
    }


def _print_stage_size_estimate(stage_name: str, estimate: Dict[str, Any]) -> None:
    """Print a compact size estimate for an exported stage."""

    print(
        f"-- Estimated exported size for stage '{stage_name}': "
        f"{estimate['total_human']} "
        f"(zarr arrays: {estimate['zarr_human']}, "
        f"pickle metadata: {estimate['pickle_human']}, "
        f"entries: {estimate['entries']})"
    )
    print()


def _dask_chunks_for_obj(obj: Union[xr.DataArray, xr.Dataset]) -> Optional[Dict[str, int]]:
    """Return safe rechunk sizes for Dask-backed xarray objects.

    Xarray/Zarr requires regular chunking where the final chunk is not larger
    than the first chunk. Some processing steps can leave arrays with chunks
    such as ``(3, 3, ..., 4)`` along time. That is valid for Dask but rejected
    by Zarr. Rechunking each affected dimension to the largest existing chunk
    preserves lazy writing while making the chunk layout Zarr-safe.
    """

    chunks: Dict[str, int] = {}

    if isinstance(obj, xr.DataArray):
        arrays = [obj]
    else:
        arrays = list(obj.data_vars.values())

    for arr in arrays:
        if arr.chunks is None:
            continue

        for dim, dim_chunks in zip(arr.dims, arr.chunks):
            finite_chunks = [int(c) for c in dim_chunks if int(c) > 0]
            if len(finite_chunks) == 0:
                continue

            chunk_size = max(finite_chunks)
            if dim not in chunks or chunk_size > chunks[dim]:
                chunks[dim] = chunk_size

    if len(chunks) == 0:
        return None

    return chunks


def _make_zarr_safe_chunks(obj: Union[xr.DataArray, xr.Dataset]) -> Union[xr.DataArray, xr.Dataset]:
    """Return obj with Dask chunks compatible with xarray.to_zarr."""

    chunks = _dask_chunks_for_obj(obj)

    if chunks is None:
        return obj

    return obj.chunk(chunks)


def _write_xarray_to_zarr(obj: Union[xr.DataArray, xr.Dataset], path: str) -> str:
    """Write one xarray object to a zarr folder and return the stored kind."""

    if os.path.isdir(path):
        shutil.rmtree(path)

    os.makedirs(os.path.dirname(path), exist_ok=True)

    obj = _make_zarr_safe_chunks(obj)

    if isinstance(obj, xr.DataArray):
        ds = obj.to_dataset(name=_DATAARRAY_VAR)
        ds.to_zarr(path, mode="w")
        return "zarr_dataarray"

    obj.to_zarr(path, mode="w")
    return "zarr_dataset"


def _read_xarray_from_zarr(path: str, kind: str) -> Union[xr.DataArray, xr.Dataset]:
    """Open one zarr entry lazily."""

    ds = xr.open_zarr(path)

    if kind == "zarr_dataarray":
        return ds[_DATAARRAY_VAR]

    if kind == "zarr_dataset":
        return ds

    raise ValueError(f"Unsupported zarr entry kind: {kind}")


def _read_export_manifest(output_folder: str, stage_name: str) -> tuple[str, Dict[str, Any]]:
    """Return exported stage directory and manifest."""

    in_dir = _stage_dir(output_folder, stage_name)
    manifest_path = os.path.join(in_dir, _MANIFEST_NAME)

    if not os.path.isfile(manifest_path):
        raise FileNotFoundError(
            f"Could not find exported stage manifest: {manifest_path}"
        )

    return in_dir, _read_pickle(manifest_path)


def _read_export_entry(in_dir: str, entry: Dict[str, Any]) -> Any:
    """Read one manifest entry only."""

    kind = entry["kind"]
    path = os.path.join(in_dir, entry["path"])

    if kind in ["zarr_dataarray", "zarr_dataset"]:
        return _read_xarray_from_zarr(path, kind)

    if kind == "pickle":
        return _read_pickle(path)

    raise ValueError(
        f"Unsupported export entry kind {kind!r} for "
        f"{entry.get('qa_test')}/{entry.get('parameter')}"
    )


def _entry_subdir(stage_dir: str, qa_test: str, parameter: str) -> str:
    return os.path.join(stage_dir, _safe_name(qa_test), _safe_name(parameter))


def export_processing_stage(
    stage_data: Dict[str, Dict[str, Any]],
    stage_name: Optional[str] = None,
    output_folder: Optional[str] = None,
    overwrite: bool = False,
    ask: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Optional[str]:
    """
    Export a full ATLAS processing stage to output_folder/exported.

    Parameters
    ----------
    stage_data
        Nested stage dictionary, normally returned by
        ``processor.export_test_from_stage(stage_name)``. Expected structure:
        ``stage_data[qa_test][parameter]``.
    stage_name
        Name of the processing stage, e.g. ``'preprocessing_complete'``.
    output_folder
        Base output folder where the ``exported`` directory is stored.
    overwrite
        If False, raise an error when the stage export already exists. If True,
        replace the existing exported stage.
    ask
        If True, show a terminal yes/no prompt before exporting. If the user
        answers no, the function returns None and does not write anything.
    default_answer
        Default answer used by the prompt when the user presses Enter.
    print_estimated_size
        If True, print an approximate exported size before writing. The estimate
        uses xarray nbytes and does not compute Dask-backed arrays.

    Returns
    -------
    str or None
        Path to the exported stage folder. Returns None when ask=True and the
        user chooses not to export.
    """

    if not isinstance(stage_data, dict):
        raise TypeError("stage_data must be a dictionary")

    if output_folder is None:
        raise ValueError("output_folder is required for stage export/import")

    if stage_name is None:
        existing_stages = list_exported_stages(output_folder)
        if len(existing_stages) == 0:
            raise ValueError(
                "stage_name was not provided and no existing exported stages "
                "are available to select. Please pass stage_name explicitly."
            )
        stage_name = select_exported_stage(output_folder)
    else:
        stage_name = str(stage_name)

    out_dir = _stage_dir(output_folder, stage_name)

    size_estimate = None
    if print_estimated_size or ask:
        size_estimate = estimate_processing_stage_size(stage_data)

    if print_estimated_size and size_estimate is not None:
        _print_stage_size_estimate(stage_name, size_estimate)

    if ask:
        do_export = _ask_yes_no(
            question=f"Export processing stage '{stage_name}' to '{out_dir}'?",
            default=default_answer,
        )

        if not do_export:
            print(f"-- Skipping export of processing stage: {stage_name}")
            return None


    if os.path.exists(out_dir):
        if not overwrite:
            raise FileExistsError(
                f"Exported stage already exists: {out_dir}. "
                "Use overwrite=True to replace it."
            )
        shutil.rmtree(out_dir)

    os.makedirs(out_dir, exist_ok=True)

    manifest = {
        "stage_name": str(stage_name),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "format": "atlas_processing_stage_export_v1",
        "entries": {},
    }

    for qa_test, qa_dict in stage_data.items():
        if not isinstance(qa_dict, dict):
            raise TypeError(
                f"stage_data[{qa_test!r}] must be a dictionary of parameters"
            )

        qa_key = str(qa_test)
        manifest["entries"].setdefault(qa_key, {})

        for parameter, value in qa_dict.items():
            param_key = str(parameter)
            entry_dir = _entry_subdir(out_dir, qa_key, param_key)

            # Numeric/lazy xarray objects are stored as zarr. Mixed object
            # arrays, such as system_info/channel_info/pol_cal_info, are stored
            # with pickle to preserve their dtype and coordinates safely.
            if isinstance(value, (xr.DataArray, xr.Dataset)) and not _has_object_dtype_xarray(value):
                rel_path = os.path.relpath(os.path.join(entry_dir, "data.zarr"), out_dir)
                kind = _write_xarray_to_zarr(value, os.path.join(out_dir, rel_path))
            else:
                rel_path = os.path.relpath(os.path.join(entry_dir, "data.pkl"), out_dir)
                _write_pickle(value, os.path.join(out_dir, rel_path))
                kind = "pickle"

            manifest["entries"][qa_key][param_key] = {
                "kind": kind,
                "path": rel_path,
                "qa_test": qa_key,
                "parameter": param_key,
            }

    _write_pickle(manifest, os.path.join(out_dir, _MANIFEST_NAME))

    return out_dir


def import_processing_stage(
    output_folder: str,
    stage_name: Optional[str] = None,
    qa_tests: Optional[Iterable[str]] = None,
    parameters: Optional[Iterable[str]] = None,
    qa_test: Optional[str] = None,
    parameter: Optional[str] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Import one exported processing stage.

    Xarray entries stored as Zarr are reopened lazily with ``xr.open_zarr``.
    Pickled entries are loaded normally.

    Parameters
    ----------
    output_folder
        Base output folder where the ``exported`` directory is stored.
    stage_name
        Name of the exported stage to read.
    qa_tests
        Optional list/set of first-level QA-test keys to import. If None, all
        QA-test keys are imported.
    parameters
        Optional list/set of parameter keys to import for each selected QA-test.
        If None, all parameters are imported.
    qa_test
        Optional single QA-test key. Convenience alternative to qa_tests=[qa_test].
    parameter
        Optional single parameter key. Convenience alternative to
        parameters=[parameter].

    Returns
    -------
    dict
        Nested dictionary with the same ``qa_test -> parameter -> value`` layout.
    """

    if stage_name is None:
        stage_name = select_exported_stage(output_folder)
    else:
        stage_name = str(stage_name)

    in_dir, manifest = _read_export_manifest(output_folder, stage_name)

    if qa_test is not None:
        if qa_tests is not None:
            raise ValueError("Use either qa_test or qa_tests, not both")
        qa_tests = [qa_test]

    if parameter is not None:
        if parameters is not None:
            raise ValueError("Use either parameter or parameters, not both")
        parameters = [parameter]

    qa_filter = None if qa_tests is None else {str(q) for q in qa_tests}
    param_filter = None if parameters is None else {str(p) for p in parameters}

    stage_data: Dict[str, Dict[str, Any]] = {}

    for qa_test, param_entries in manifest.get("entries", {}).items():
        if qa_filter is not None and qa_test not in qa_filter:
            continue

        stage_data[qa_test] = {}

        for parameter, entry in param_entries.items():
            if param_filter is not None and parameter not in param_filter:
                continue

            stage_data[qa_test][parameter] = _read_export_entry(in_dir, entry)

    return stage_data


def import_processing_entry(
    output_folder: str,
    stage_name: Optional[str] = None,
    qa_test: Optional[str] = None,
    parameter: Optional[str] = None,
) -> Any:
    """Import only one selected part of an exported processing stage.

    This is faster than importing a full stage because it opens/loads only the
    requested manifest entries. Zarr-backed xarray arrays are still opened
    lazily.

    Parameters
    ----------
    output_folder
        Base output folder where the ``exported`` directory is stored.
    stage_name
        Name of the exported stage to read.
    qa_test
        Optional first-level QA-test key, e.g. ``"ray"``.
    parameter
        Optional second-level parameter key, e.g. ``"profile"``.

    Returns
    -------
    Any
        If both qa_test and parameter are provided, returns the single value.
        If only qa_test is provided, returns ``{parameter: value}``.
        If only parameter is provided, returns ``{qa_test: {parameter: value}}``
        for all QA tests that contain that parameter.
        If neither is provided, returns the full stage via import_processing_stage.
    """

    if stage_name is None:
        stage_name = select_exported_stage(output_folder)
    else:
        stage_name = str(stage_name)

    if qa_test is None and parameter is None:
        return import_processing_stage(output_folder=output_folder, stage_name=stage_name)

    in_dir, manifest = _read_export_manifest(output_folder, stage_name)
    entries = manifest.get("entries", {})

    if qa_test is not None:
        qa_key = str(qa_test)

        if qa_key not in entries:
            raise KeyError(
                f"QA-test key {qa_key!r} was not found in exported stage "
                f"{stage_name!r}. Available keys: {sorted(entries.keys())}"
            )

        if parameter is not None:
            param_key = str(parameter)

            if param_key not in entries[qa_key]:
                raise KeyError(
                    f"Parameter key {param_key!r} was not found under "
                    f"{qa_key!r} in exported stage {stage_name!r}. "
                    f"Available parameters: {sorted(entries[qa_key].keys())}"
                )

            return _read_export_entry(in_dir, entries[qa_key][param_key])

        return {
            param_key: _read_export_entry(in_dir, entry)
            for param_key, entry in entries[qa_key].items()
        }

    param_key = str(parameter)
    out: Dict[str, Dict[str, Any]] = {}

    for qa_key, param_entries in entries.items():
        if param_key not in param_entries:
            continue

        out[qa_key] = {
            param_key: _read_export_entry(in_dir, param_entries[param_key])
        }

    if len(out) == 0:
        raise KeyError(
            f"Parameter key {param_key!r} was not found in exported stage "
            f"{stage_name!r}."
        )

    return out


def delete_exported_stage(
    output_folder: str,
    stage_name: Optional[str] = None,
    ask: bool = False,
    default_answer: bool = False,
    missing_ok: bool = False,
) -> Optional[str]:
    """Delete one exported stage folder.

    Returns the deleted path. If missing_ok=True and the stage does not exist,
    returns None.
    """

    if stage_name is None:
        stage_name = select_exported_stage(output_folder)
    else:
        stage_name = str(stage_name)

    stage_path = _stage_dir(output_folder, stage_name)

    if not os.path.exists(stage_path):
        if missing_ok:
            print(f"-- Exported stage does not exist, nothing to delete: {stage_path}")
            return None
        raise FileNotFoundError(f"Exported stage does not exist: {stage_path}")

    if ask:
        do_delete = _ask_yes_no(
            question=f"Delete exported processing stage '{stage_name}' at '{stage_path}'?",
            default=default_answer,
        )

        if not do_delete:
            print(f"-- Skipping deletion of exported stage: {stage_name}")
            return None

    shutil.rmtree(stage_path)
    print(f"-- Deleted exported processing stage: {stage_path}")
    print()

    return stage_path


def delete_all_exported_stages(
    output_folder: str,
    ask: bool = True,
    default_answer: bool = False,
    missing_ok: bool = True,
) -> Optional[str]:
    """Delete output_folder/exported and all saved stages.

    By default this asks for confirmation because it removes every saved stage.
    Returns the deleted root path. If the folder is missing and missing_ok=True,
    returns None.
    """

    root = _export_root(output_folder)

    if not os.path.exists(root):
        if missing_ok:
            print(f"-- Exported stages folder does not exist, nothing to delete: {root}")
            return None
        raise FileNotFoundError(f"Exported stages folder does not exist: {root}")

    if ask:
        stages = list_exported_stages(output_folder)
        stage_text = ", ".join(stages) if len(stages) else "no manifest-backed stages"
        do_delete = _ask_yes_no(
            question=(
                f"Delete ALL exported processing stages at '{root}' "
                f"({stage_text})?"
            ),
            default=default_answer,
        )

        if not do_delete:
            print("-- Skipping deletion of all exported processing stages")
            return None

    shutil.rmtree(root)
    print(f"-- Deleted all exported processing stages: {root}")
    print()

    return root



def estimate_processing_stages_size(
    stages_data: Dict[str, Dict[str, Dict[str, Any]]],
) -> Dict[str, Any]:
    """Estimate the combined exported size of multiple processing stages.

    Parameters
    ----------
    stages_data
        Dictionary with layout ``stages_data[stage_name][qa_test][parameter]``.

    Returns
    -------
    dict
        Combined totals plus a per-stage estimate dictionary.
    """

    if not isinstance(stages_data, dict):
        raise TypeError("stages_data must be a dictionary")

    per_stage: Dict[str, Dict[str, Any]] = {}
    total_bytes = 0
    zarr_bytes = 0
    pickle_bytes = 0
    entries = 0
    xarray_entries = 0
    pickle_entries = 0

    for stage_name, stage_data in stages_data.items():
        estimate = estimate_processing_stage_size(stage_data)
        per_stage[str(stage_name)] = estimate

        total_bytes += estimate["total_bytes"]
        zarr_bytes += estimate["zarr_bytes"]
        pickle_bytes += estimate["pickle_bytes"]
        entries += estimate["entries"]
        xarray_entries += estimate["xarray_entries"]
        pickle_entries += estimate["pickle_entries"]

    return {
        "total_bytes": int(total_bytes),
        "zarr_bytes": int(zarr_bytes),
        "pickle_bytes": int(pickle_bytes),
        "entries": int(entries),
        "xarray_entries": int(xarray_entries),
        "pickle_entries": int(pickle_entries),
        "total_human": _human_readable_size(total_bytes),
        "zarr_human": _human_readable_size(zarr_bytes),
        "pickle_human": _human_readable_size(pickle_bytes),
        "per_stage": per_stage,
    }


def _print_stages_size_estimate(estimate: Dict[str, Any]) -> None:
    """Print a compact size estimate for multiple exported stages."""

    per_stage = estimate.get("per_stage", {})
    print(
        f"-- Estimated exported size for {len(per_stage)} stage(s): "
        f"{estimate['total_human']} "
        f"(zarr arrays: {estimate['zarr_human']}, "
        f"pickle metadata: {estimate['pickle_human']}, "
        f"entries: {estimate['entries']})"
    )

    for stage_name in sorted(per_stage.keys()):
        stage_estimate = per_stage[stage_name]
        print(
            f"   - {stage_name}: {stage_estimate['total_human']} "
            f"(entries: {stage_estimate['entries']})"
        )

    print()


def _check_stage_export_targets(
    output_folder: str,
    stage_names: Iterable[str],
    overwrite: bool,
) -> None:
    """Validate output targets before a batch export writes anything."""

    safe_to_original: Dict[str, str] = {}
    existing: list[str] = []

    for stage_name in stage_names:
        stage_key = str(stage_name)
        safe = _safe_name(stage_key)

        if safe in safe_to_original and safe_to_original[safe] != stage_key:
            raise ValueError(
                "Two stage names resolve to the same export folder after "
                f"filesystem-safe conversion: {safe_to_original[safe]!r} and "
                f"{stage_key!r} -> {safe!r}"
            )

        safe_to_original[safe] = stage_key
        out_dir = _stage_dir(output_folder, stage_key)

        if os.path.exists(out_dir):
            existing.append(out_dir)

    if existing and not overwrite:
        raise FileExistsError(
            "One or more exported stages already exist. Use overwrite=True to "
            "replace them:\n" + "\n".join(f"  - {path}" for path in existing)
        )


def export_processing_stages(
    stages_data: Dict[str, Dict[str, Dict[str, Any]]],
    output_folder: str,
    stage_names: Optional[Iterable[str]] = None,
    overwrite: bool = False,
    ask: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Dict[str, str]:
    """Export multiple ATLAS processing stages with one optional prompt.

    Parameters
    ----------
    stages_data
        Dictionary with layout ``stages_data[stage_name][qa_test][parameter]``.
    output_folder
        Base output folder where the ``exported`` directory is stored.
    overwrite
        If False, raise an error when any target stage export already exists.
        If True, replace existing exported stages.
    ask
        If True, show one terminal yes/no prompt before exporting any stage.
    default_answer
        Default answer used by the prompt when the user presses Enter.
    print_estimated_size
        If True, print one combined approximate exported size before writing.

    Returns
    -------
    dict
        Mapping ``stage_name -> exported stage folder``. Returns an empty dict
        when ask=True and the user chooses not to export.
    """

    if not isinstance(stages_data, dict):
        raise TypeError("stages_data must be a dictionary")

    available_stage_names = [str(stage_name) for stage_name in stages_data.keys()]

    if len(available_stage_names) == 0:
        raise ValueError("stages_data must contain at least one stage")

    if stage_names is None:
        stage_names = list(select_stage_names(available_stage_names, allow_multiple=True))
    else:
        stage_names = [str(stage_name) for stage_name in stage_names]

    missing_stage_names = [
        stage_name for stage_name in stage_names
        if stage_name not in stages_data
    ]
    if len(missing_stage_names) > 0:
        raise KeyError(
            "Selected stage name(s) are not available in stages_data: "
            + ", ".join(missing_stage_names)
        )

    _check_stage_export_targets(
        output_folder=output_folder,
        stage_names=stage_names,
        overwrite=overwrite,
    )

    selected_stages_data = {
        stage_name: stages_data[stage_name]
        for stage_name in stage_names
    }

    size_estimate = None
    if print_estimated_size or ask:
        size_estimate = estimate_processing_stages_size(selected_stages_data)

    if print_estimated_size and size_estimate is not None:
        _print_stages_size_estimate(size_estimate)

    if ask:
        root = _export_root(output_folder)
        stage_text = ", ".join(stage_names)
        size_text = "unknown"

        if size_estimate is not None:
            size_text = size_estimate["total_human"]

        do_export = _ask_yes_no(
            question=(
                f"Export {len(stage_names)} processing stage(s) "
                f"({stage_text}) to '{root}'? "
                f"Estimated total size: {size_text}."
            ),
            default=default_answer,
        )

        if not do_export:
            print("-- Skipping export of processing stages: " + stage_text)
            return {}

    exported: Dict[str, str] = {}

    for stage_name in stage_names:
        out_dir = export_processing_stage(
            stage_data=stages_data[stage_name],
            stage_name=stage_name,
            output_folder=output_folder,
            overwrite=overwrite,
            ask=False,
            default_answer=default_answer,
            print_estimated_size=False,
        )

        if out_dir is not None:
            exported[stage_name] = out_dir

    return exported


def export_processor_stage(
    processor: Any,
    output_folder: 'str',
    stage_name: Optional[str] = None,
    overwrite: bool = False,
    ask: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Optional[str]:
    """
    Convenience wrapper that exports a stage directly from an ATLAS Processor.

    Set ask=True to show a terminal yes/no prompt before writing anything.
    """

    stage_name = _resolve_processor_stage_name(processor, stage_name)
    stage_data = processor.export_test_from_stage(stage_name)

    return export_processing_stage(
        stage_data=stage_data,
        stage_name=stage_name,
        output_folder=output_folder,
        overwrite=overwrite,
        ask=ask,
        default_answer=default_answer,
        print_estimated_size=print_estimated_size,
    )



def export_processor_stages(
    processor: Any,
    output_folder: str,
    stage_names: Optional[Iterable[str]] = None,
    overwrite: bool = False,
    ask: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Dict[str, str]:
    """Convenience wrapper that exports multiple stages from an ATLAS Processor.

    Set ask=True to show exactly one terminal yes/no prompt before writing any
    of the requested stages. The printed size estimate is the combined total for
    all stages.
    """

    stage_names = _resolve_processor_stage_names(processor, stage_names)

    if len(stage_names) == 0:
        raise ValueError("stage_names must contain at least one stage")

    stages_data = {
        stage_name: processor.export_test_from_stage(stage_name)
        for stage_name in stage_names
    }

    return export_processing_stages(
        stages_data=stages_data,
        output_folder=output_folder,
        stage_names=stage_names,
        overwrite=overwrite,
        ask=ask,
        default_answer=default_answer,
        print_estimated_size=print_estimated_size,
    )


def ask_export_processor_stages(
    processor: Any,
    stage_names: Optional[Iterable[str]] = None,
    overwrite: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Dict[str, str]:
    """Prompt once in the terminal before exporting multiple stages."""

    return export_processor_stages(
        processor=processor,
        stage_names=stage_names,
        overwrite=overwrite,
        ask=True,
        default_answer=default_answer,
        print_estimated_size=print_estimated_size,
    )

def list_exported_stages(output_folder: str) -> list[str]:
    """Return available exported stage folder names."""

    root = _export_root(output_folder)

    return sorted(
        name
        for name in os.listdir(root)
        if os.path.isdir(os.path.join(root, name))
        and os.path.isfile(os.path.join(root, name, _MANIFEST_NAME))
    )


def ask_export_processor_stage(
    processor: Any,
    stage_name: Optional[str] = None,
    overwrite: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Optional[str]:
    """Prompt in the terminal before exporting one processor stage."""

    return export_processor_stage(
        processor=processor,
        stage_name=stage_name,
        overwrite=overwrite,
        ask=True,
        default_answer=default_answer,
        print_estimated_size=print_estimated_size,
    )


def inspect_exported_stages(
    output_folder: str,
    stage_name: Optional[str] = None,
    print_tree: bool = True,
    include_kinds: bool = True,
    include_paths: bool = False,
) -> Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]:
    """Inspect saved processing stages without loading the exported data.

    The function reads only the small stage manifest files. It does not open
    Zarr arrays and it does not unpickle the stored data entries, so it is fast
    even when the exported stages contain large lazy arrays.

    Parameters
    ----------
    output_folder
        Base output folder where the ``exported`` directory is stored.
    stage_name
        Optional stage name to inspect. If None, all exported stages are
        inspected.
    print_tree
        If True, print a readable tree of stages, QA-test keys, and parameter
        keys to the terminal.
    include_kinds
        If True, include the storage kind next to each parameter in the printed
        tree, for example ``zarr_dataarray`` or ``pickle``.
    include_paths
        If True, include the relative stored path next to each parameter in the
        printed tree.

    Returns
    -------
    dict
        Nested summary dictionary with this layout::

            summary[stage_name][qa_test][parameter] = {
                "kind": "zarr_dataarray" | "zarr_dataset" | "pickle",
                "path": "relative/path/to/data",
            }
    """

    if stage_name is None:
        stage_names = list_exported_stages(output_folder)
    else:
        stage_names = [str(stage_name)]

    summary: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = {}

    for stage in stage_names:
        _, manifest = _read_export_manifest(output_folder, stage)
        stage_entries = manifest.get("entries", {})
        summary[stage] = {}

        for qa_key, param_entries in stage_entries.items():
            summary[stage][qa_key] = {}

            for param_key, entry in param_entries.items():
                summary[stage][qa_key][param_key] = {
                    "kind": entry.get("kind"),
                    "path": entry.get("path"),
                }

    if print_tree:
        print_exported_stage_tree(
            summary,
            include_kinds=include_kinds,
            include_paths=include_paths,
        )

    return summary


def print_exported_stage_tree(
    summary: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]],
    include_kinds: bool = True,
    include_paths: bool = False,
) -> None:
    """Print a readable tree returned by ``inspect_exported_stages``."""

    if len(summary) == 0:
        print("-- No exported processing stages found.")
        return

    print("-- Exported processing stages:")

    for stage in sorted(summary.keys()):
        qa_entries = summary[stage]
        n_params = sum(len(params) for params in qa_entries.values())
        print(f"{stage}/  ({len(qa_entries)} qa keys, {n_params} parameters)")

        for qa_key in sorted(qa_entries.keys()):
            params = qa_entries[qa_key]
            print(f"  {qa_key}/  ({len(params)} parameters)")

            for param_key in sorted(params.keys()):
                entry = params[param_key]
                suffix_parts = []

                if include_kinds and entry.get("kind") is not None:
                    suffix_parts.append(str(entry.get("kind")))

                if include_paths and entry.get("path") is not None:
                    suffix_parts.append(str(entry.get("path")))

                suffix = ""
                if len(suffix_parts) > 0:
                    suffix = "  [" + "; ".join(suffix_parts) + "]"

                print(f"    - {param_key}{suffix}")


def list_exported_stage_contents(
    output_folder: str,
    stage_name: Optional[str] = None,
) -> Dict[str, Dict[str, list[str]]]:
    """Return a compact stage -> qa_test -> parameter list summary.

    This is useful in scripts when you only need the available keys and not the
    storage details.
    """

    detailed = inspect_exported_stages(
        output_folder=output_folder,
        stage_name=stage_name,
        print_tree=False,
    )

    return {
        stage: {
            qa_key: sorted(params.keys())
            for qa_key, params in qa_entries.items()
        }
        for stage, qa_entries in detailed.items()
    }


# Short aliases, in case you prefer shorter imports.
export_stage = export_processing_stage
export_stages = export_processing_stages
export_processor_stage_list = export_processor_stages
import_stage = import_processing_stage
import_entry = import_processing_entry
delete_stage = delete_exported_stage
delete_all_stages = delete_all_exported_stages
inspect_stages = inspect_exported_stages
list_stage_contents = list_exported_stage_contents
select_stage = select_exported_stage
select_stages = select_exported_stages
