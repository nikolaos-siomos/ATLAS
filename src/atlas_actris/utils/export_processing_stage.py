#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities to export and import ATLAS processor stages.

The exported layout is:

    caller_info['output_folder']/exported/<stage_name>/
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


def _export_root(caller_info: Dict[str, Any]) -> str:
    """Return caller_info['output_folder']/exported, creating it if needed."""

    output_folder = caller_info.get("output_folder")

    if output_folder is None:
        raise ValueError("caller_info['output_folder'] is required for stage export/import")

    root = os.path.join(output_folder, "exported")
    os.makedirs(root, exist_ok=True)
    return root


def _stage_dir(caller_info: Dict[str, Any], stage_name: str) -> str:
    """Return the export folder of a specific stage."""

    return os.path.join(_export_root(caller_info), _safe_name(stage_name))


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


def _read_export_manifest(caller_info: Dict[str, Any], stage_name: str) -> tuple[str, Dict[str, Any]]:
    """Return exported stage directory and manifest."""

    in_dir = _stage_dir(caller_info, stage_name)
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
    stage_name: str,
    caller_info: Dict[str, Any],
    overwrite: bool = False,
    ask: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Optional[str]:
    """
    Export a full ATLAS processing stage to caller_info['output_folder']/exported.

    Parameters
    ----------
    stage_data
        Nested stage dictionary, normally returned by
        ``processor.export_test_from_stage(stage_name)``. Expected structure:
        ``stage_data[qa_test][parameter]``.
    stage_name
        Name of the processing stage, e.g. ``'preprocessing_complete'``.
    caller_info
        ATLAS caller_info dictionary. Must contain ``output_folder``.
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

    out_dir = _stage_dir(caller_info, stage_name)

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
    caller_info: Dict[str, Any],
    stage_name: str,
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
    caller_info
        ATLAS caller_info dictionary. Must contain ``output_folder``.
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

    in_dir, manifest = _read_export_manifest(caller_info, stage_name)

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
    caller_info: Dict[str, Any],
    stage_name: str,
    qa_test: Optional[str] = None,
    parameter: Optional[str] = None,
) -> Any:
    """Import only one selected part of an exported processing stage.

    This is faster than importing a full stage because it opens/loads only the
    requested manifest entries. Zarr-backed xarray arrays are still opened
    lazily.

    Parameters
    ----------
    caller_info
        ATLAS caller_info dictionary. Must contain ``output_folder``.
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

    if qa_test is None and parameter is None:
        return import_processing_stage(caller_info=caller_info, stage_name=stage_name)

    in_dir, manifest = _read_export_manifest(caller_info, stage_name)
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
    caller_info: Dict[str, Any],
    stage_name: str,
    ask: bool = False,
    default_answer: bool = False,
    missing_ok: bool = False,
) -> Optional[str]:
    """Delete one exported stage folder.

    Returns the deleted path. If missing_ok=True and the stage does not exist,
    returns None.
    """

    stage_path = _stage_dir(caller_info, stage_name)

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
    caller_info: Dict[str, Any],
    ask: bool = True,
    default_answer: bool = False,
    missing_ok: bool = True,
) -> Optional[str]:
    """Delete caller_info['output_folder']/exported and all saved stages.

    By default this asks for confirmation because it removes every saved stage.
    Returns the deleted root path. If the folder is missing and missing_ok=True,
    returns None.
    """

    root = os.path.join(caller_info.get("output_folder", ""), "exported")

    if not os.path.exists(root):
        if missing_ok:
            print(f"-- Exported stages folder does not exist, nothing to delete: {root}")
            return None
        raise FileNotFoundError(f"Exported stages folder does not exist: {root}")

    if ask:
        stages = list_exported_stages(caller_info)
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


def export_processor_stage(
    processor: Any,
    stage_name: str,
    overwrite: bool = False,
    ask: bool = False,
    default_answer: bool = False,
    print_estimated_size: bool = True,
) -> Optional[str]:
    """
    Convenience wrapper that exports a stage directly from an ATLAS Processor.

    Set ask=True to show a terminal yes/no prompt before writing anything.
    """

    stage_data = processor.export_test_from_stage(stage_name)
    caller_info = processor.processing_info["caller_info"]

    return export_processing_stage(
        stage_data=stage_data,
        stage_name=stage_name,
        caller_info=caller_info,
        overwrite=overwrite,
        ask=ask,
        default_answer=default_answer,
        print_estimated_size=print_estimated_size,
    )


def list_exported_stages(caller_info: Dict[str, Any]) -> list[str]:
    """Return available exported stage folder names."""

    root = _export_root(caller_info)

    return sorted(
        name
        for name in os.listdir(root)
        if os.path.isdir(os.path.join(root, name))
        and os.path.isfile(os.path.join(root, name, _MANIFEST_NAME))
    )


def ask_export_processor_stage(
    processor: Any,
    stage_name: str,
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
    caller_info: Dict[str, Any],
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
    caller_info
        ATLAS caller_info dictionary. Must contain ``output_folder``.
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
        stage_names = list_exported_stages(caller_info)
    else:
        stage_names = [str(stage_name)]

    summary: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = {}

    for stage in stage_names:
        _, manifest = _read_export_manifest(caller_info, stage)
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
    caller_info: Dict[str, Any],
    stage_name: Optional[str] = None,
) -> Dict[str, Dict[str, list[str]]]:
    """Return a compact stage -> qa_test -> parameter list summary.

    This is useful in scripts when you only need the available keys and not the
    storage details.
    """

    detailed = inspect_exported_stages(
        caller_info=caller_info,
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
import_stage = import_processing_stage
import_entry = import_processing_entry
delete_stage = delete_exported_stage
delete_all_stages = delete_all_exported_stages
inspect_stages = inspect_exported_stages
list_stage_contents = list_exported_stage_contents
