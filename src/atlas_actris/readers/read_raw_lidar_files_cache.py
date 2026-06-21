"""
@author: Nikos Siomos
"""
import re
import numpy as np
import xarray as xr
from datetime import datetime
import io, os, glob, contextlib, pickle, shutil

from typing import Callable
from typing import Any, Dict, List, Tuple

from utils.error_classes import FileReaderError
from utils.printouts import print_header, print_subsection, endpoint

from readers.read_scc_lazy import read_dataset as reader_scc
from readers.read_licel_lazy import read_dataset as reader_licel
from readers.read_tamarin_lazy import read_dataset as reader_tamarin
from readers.read_polly_xt_lazy import read_dataset as reader_polly_xt
from readers.read_licel_matlab_lazy import read_dataset as reader_licel_matlab
from readers.read_polly_xt_first_lazy import read_dataset as reader_polly_xt_first

# Example registry of formats -> reading functions
READERS: dict[str, Callable[[str], object]] = {
    "scc": reader_scc,
    "licel": reader_licel,
    "polly_xt": reader_polly_xt,
}

reader_menu: dict[str, Callable[[str], object]] = {
    "scc": reader_scc,
    "licel": reader_licel,
    "tamarin": reader_tamarin,
    "polly_xt": reader_polly_xt,
    "licel_matlab": reader_licel_matlab,
    "polly_xt_first": reader_polly_xt_first,
}


def measurement_type(meas_key: str) -> str:
    """
    Infer the raw-reader measurement type from a caller_info measurement key.

    This mirrors the parser-side measurement_type helper and replaces the old
    flat mtype_* entries in caller_info. The returned value is still the value
    expected by the raw file readers.
    """

    if meas_key.startswith("drk"):
        return "drk"

    if meas_key.startswith("ray"):
        return "nrm"

    if meas_key in ["trg", "dtm", "nsf"]:
        return "nrm"

    if meas_key in ["tlc_north", "tlc_east", "tlc_south", "tlc_west"]:
        return "tlc_qua"

    if meas_key in ["tlc_inner", "tlc_outer"]:
        return "tlc_rin"

    if meas_key == "pcb_p45":
        return "pcb_p45"

    if meas_key == "pcb_m45":
        return "pcb_m45"

    if meas_key == "pcb_aux_p45":
        return "pcb_p45"

    if meas_key == "pcb_aux_m45":
        return "pcb_m45"

    if meas_key.startswith("cam"):
        return "cam"

    raise FileReaderError(f"Could not infer measurement type for {meas_key}")


def _iter_physical_paths(d: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Return physical measurement folders from caller_info['paths']."""

    paths = d.get("paths", {})

    if paths is None:
        paths = {}

    if not isinstance(paths, dict):
        raise FileReaderError("caller_info['paths'] must be a dictionary")

    return [
        (meas_key, filepath)
        for meas_key, filepath in paths.items()
        if filepath is not None and not meas_key.startswith("cam")
    ]

def downcast_float_safe(da: xr.DataArray, tol: float = 1e-6) -> xr.DataArray:
    """
    Downcast float64 DataArray to float32 when possible.

    For NumPy-backed arrays, a safety check is performed before downcasting.
    For Dask-backed arrays, the downcast is returned lazily, without computing
    the full array. This preserves the RAM benefit of the lazy reader.
    """

    if not (np.issubdtype(da.dtype, np.floating) and da.dtype == np.float64):
        return da

    # Dask-backed arrays expose ``chunks``. Performing a full allclose check here
    # would force a full read/compute of the raw signal, defeating lazy loading.
    if da.chunks is not None:
        return da.astype(np.float32)

    arr = da.values
    arr32 = arr.astype(np.float32)

    if np.allclose(arr, arr32.astype(np.float64), rtol=tol, atol=tol, equal_nan=True):
        return da.astype(np.float32)

    return da


def _raw_cache_root(d: Dict[str, Any]) -> str:
    """Return the folder used for temporary raw-lidar Zarr caches.

    The cache is always placed directly under caller_info["output_folder"] / "cache".
    """

    output_folder = d.get("output_folder")

    if output_folder is None:
        raise FileReaderError(
            "caller_info['output_folder'] is required for raw-lidar caching"
        )

    return os.path.join(output_folder, "cache")

def _safe_cache_key(meas_key: str) -> str:
    """Convert a measurement key into a filesystem-safe cache key."""

    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(meas_key))


def _cache_entry_dir(cache_root: str, meas_key: str) -> str:
    return os.path.join(cache_root, _safe_cache_key(meas_key))


def _cache_signal_path(cache_root: str, meas_key: str) -> str:
    return os.path.join(_cache_entry_dir(cache_root, meas_key), "data.zarr")


def _cache_metadata_path(cache_root: str, meas_key: str) -> str:
    return os.path.join(_cache_entry_dir(cache_root, meas_key), "metadata.pkl")


def _cache_entry_exists(cache_root: str, meas_key: str) -> bool:
    """Check whether one measurement cache entry is complete enough to read."""

    return (
        os.path.isdir(_cache_signal_path(cache_root, meas_key))
        and os.path.isfile(_cache_metadata_path(cache_root, meas_key))
    )


def _cache_is_full(cache_root: str, physical_paths: List[Tuple[str, str]]) -> bool:
    """Return True only if every physical measurement has a cache entry."""

    if len(physical_paths) == 0:
        return False

    return all(_cache_entry_exists(cache_root, meas_key) for meas_key, _ in physical_paths)


def _clear_raw_cache(cache_root: str) -> None:
    """Remove the complete raw cache folder and recreate it empty."""

    if os.path.isdir(cache_root):
        shutil.rmtree(cache_root)

    os.makedirs(cache_root, exist_ok=True)


def _write_cache_entry(
    cache_root: str,
    meas_key: str,
    system_info: Any,
    channel_info: Any,
    time_info: Any,
    sig_raw: xr.DataArray,
    shots: xr.DataArray,
    meas_type: str,
) -> None:
    """
    Write one measurement to cache.

    The signal and shots are stored in Zarr. If sig_raw is Dask-backed, this
    computes chunk-by-chunk and does not assemble the full raw signal in RAM.
    Small metadata objects are stored in pickle next to the Zarr store.
    """

    entry_dir = _cache_entry_dir(cache_root, meas_key)

    if os.path.isdir(entry_dir):
        shutil.rmtree(entry_dir)

    os.makedirs(entry_dir, exist_ok=True)

    ds = xr.Dataset(
        data_vars={
            "signal": sig_raw,
            "shots": shots,
        }
    )

    ds.to_zarr(_cache_signal_path(cache_root, meas_key), mode="w")

    metadata = {
        "system_info": system_info,
        "channel_info": channel_info,
        "time_info": time_info,
        "mtype": meas_type,
    }

    with open(_cache_metadata_path(cache_root, meas_key), "wb") as f:
        pickle.dump(metadata, f, protocol=pickle.HIGHEST_PROTOCOL)


def _read_cache_entry(cache_root: str, meas_key: str) -> Tuple[Any, Any, Any, xr.DataArray, xr.DataArray, str]:
    """Read one cached measurement lazily from Zarr plus its pickled metadata."""

    if not _cache_entry_exists(cache_root, meas_key):
        raise FileReaderError(f"Raw lidar cache entry is missing or incomplete for {meas_key}")

    ds = xr.open_zarr(_cache_signal_path(cache_root, meas_key))

    with open(_cache_metadata_path(cache_root, meas_key), "rb") as f:
        metadata = pickle.load(f)

    return (
        metadata["system_info"],
        metadata["channel_info"],
        metadata["time_info"],
        ds["signal"],
        ds["shots"].compute(),
        metadata["mtype"],
    )

def infer_format(d: Dict[str, Any], station_id: str, debug: bool = False) -> str:
    """
    Infer the raw file format by trying the registered readers on the physical
    folders listed in caller_info["paths"].

    The new parser stores only folders that should be read physically under
    caller_info["paths"].
    """

    print_header("Infering raw file format")

    # Exceptional readers are detected based on the station ID.
    if station_id == "evo":
        return "polly_xt_first"

    if station_id in ["brc", "run"]:
        return "licel_matlab"

    if station_id == "tam":
        return "tamarin"

    physical_paths = _iter_physical_paths(d)

    if len(physical_paths) == 0:
        endpoint(1)

    # For the rest of the stations, infer the reader by trial and error.
    for meas_key, filepath in physical_paths:
        meas_type = measurement_type(meas_key)

        print(f"Trying readers with {meas_key} dataset:")

        for fmt, reader in READERS.items():
            print(f"--{fmt}: ", end="")

            try:
                if debug:
                    with contextlib.redirect_stdout(io.StringIO()):
                        _, _, _, sig_raw, _ = reader(filepath, meas_type=meas_type)
                else:
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        _, _, _, sig_raw, _ = reader(filepath, meas_type=meas_type)

                if isinstance(sig_raw, list) and len(sig_raw) == 0:
                    print("Wrong reader")
                    continue

                print("Correct reader")
                print(f"File format: {fmt}")
                return fmt

            except Exception as e:
                print("Wrong reader")
                if debug:
                    print(e)
                continue

    raise FileReaderError("The raw QA file format is not supported")


def flexible_reader(
    d: Dict[str, Any],
    downscale: bool = True,
    chunk: bool = True,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, str]]:
    """
    Read only physical lidar signal folders from caller_info["paths"].

    Cache behavior is now unconditional and matches the old quick_run=True mode:

    A1. If the raw cache is complete:
        Read profiles and metadata directly from the cache. Do not write again.

    A2. If the raw cache is missing/incomplete:
        Read raw files normally and write a fresh cache.

    The cache stores signals and shots as Zarr, so reopening from cache remains
    lazy and RAM-safe. Metadata are stored as small pickle files.
    """

    print_header("Reading lidar signals...")

    file_format = d["raw_file_format"]
    cache_root = _raw_cache_root(d)

    profiles: Dict[str, Any] = {}

    metadata: Dict[str, Any] = {
        "system_info": {},
        "channel_info": {},
        "time_info": {},
        "shots": {},
        "mtype": {},
    }

    physical_paths = _iter_physical_paths(d)

    if len(physical_paths) == 0:
        endpoint(1)

    read_from_cache = _cache_is_full(cache_root, physical_paths)

    if read_from_cache:
        print(f"-- Raw cache is complete. Reading from cache:\n   {cache_root}")

        for meas_key, _ in physical_paths:
            print_subsection(f"{meas_key} dataset")

            system_info, channel_info, time_info, sig_raw, shots, meas_type = _read_cache_entry(
                cache_root=cache_root,
                meas_key=meas_key,
            )

            profiles[meas_key] = sig_raw
            metadata["system_info"][meas_key] = system_info
            metadata["channel_info"][meas_key] = channel_info
            metadata["time_info"][meas_key] = time_info
            metadata["shots"][meas_key] = shots
            metadata["mtype"][meas_key] = meas_type

        if profiles == {}:
            endpoint(1)

        return profiles, metadata

    print(f"-- Raw cache is missing or incomplete. Reading raw files and replacing raw cache:\n   {cache_root}")

    # For missing/incomplete cache, rebuild the cache so that future runs are
    # consistent and complete.
    _clear_raw_cache(cache_root)

    for meas_key, filepath in physical_paths:
        meas_type = measurement_type(meas_key)

        print_subsection(f"{meas_key} dataset")

        system_info, channel_info, time_info, sig_raw, shots = reader_menu[file_format](
            filepath,
            meas_type=meas_type,
        )

        if isinstance(sig_raw, list):
            continue

        if chunk and sig_raw.chunks is None:
            time_chunks = min(5, max(10, sig_raw.sizes["time"]))
            channel_chunks = 5
            
            chunks = {
                "time": time_chunks,
                "channel": channel_chunks,
                "bins": -1,
            }

            sig_raw = sig_raw.chunk(chunks)

        if downscale:
            sig_raw = downcast_float_safe(sig_raw)

        profiles[meas_key] = sig_raw

        metadata["system_info"][meas_key] = system_info
        metadata["channel_info"][meas_key] = channel_info
        metadata["time_info"][meas_key] = time_info
        metadata["shots"][meas_key] = shots
        metadata["mtype"][meas_key] = meas_type

        _write_cache_entry(
            cache_root=cache_root,
            meas_key=meas_key,
            system_info=system_info,
            channel_info=channel_info,
            time_info=time_info,
            sig_raw=sig_raw,
            shots=shots,
            meas_type=meas_type,
        )

        # Reopen from cache so downstream computations use the fast chunked
        # temporary store rather than the original raw-file delayed tasks.
        system_info, channel_info, time_info, sig_raw_cached, shots_cached, meas_type = _read_cache_entry(
            cache_root=cache_root,
            meas_key=meas_key,
        )

        profiles[meas_key] = sig_raw_cached
        metadata["system_info"][meas_key] = system_info
        metadata["channel_info"][meas_key] = channel_info
        metadata["time_info"][meas_key] = time_info
        metadata["shots"][meas_key] = shots_cached
        metadata["mtype"][meas_key] = meas_type

    if profiles == {}:
        endpoint(1)

    return profiles, metadata
