#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lazy overflow handling for ATLAS profiles.

Design notes
------------
- get_overflow_mask() only builds a lazy boolean mask. It does not compute.
- Each trim_overflows method computes only the reductions it needs.
- trim_overflows = 3 does not check anything; it leaves the data unchanged.
- channel_info is still broadcast_like(profile) before mask construction to
  preserve the existing alignment/dimension-order behavior.
"""

import copy
from typing import Any, Dict

import dask
import xarray as xr

from utils.printouts import endpoint, print_subsection
from utils.error_classes import DataOverflowError, CustomWarning


def _has_dask_chunks(da: xr.DataArray) -> bool:
    """Return True if a DataArray is Dask-backed/chunked."""
    return hasattr(da.data, "chunks")


def _bins_are_single_chunk(da: xr.DataArray) -> bool:
    """Return True if bins is absent or already stored as one chunk."""
    if not _has_dask_chunks(da):
        return True
    if "bins" not in da.chunksizes:
        return True
    return len(da.chunksizes["bins"]) == 1


def apply_time_mask(qa_tests, mask_time, output_data):
    """
    Apply per-measurement time masks to profile/time_info/shots.

    mask_time entries are expected to be small 1D boolean DataArrays with
    dimension time. They should already be computed, avoiding accidental
    computation of full profile-sized arrays here.
    """

    del_keys = []
    for key in qa_tests:
        if not mask_time[key].any().item():
            del_keys.append(key)

    for key in qa_tests:
        if key in del_keys:
            CustomWarning(
                f"{key} measurement will not be processed because masking "
                "removes all profiles"
            )
            del output_data["profile"][key]
            if key in output_data["time_mask"]:
                del output_data["time_mask"][key]
            del output_data["time_info"][key]
            del output_data["shots"][key]

        else:
            output_data["profile"][key] = (
                output_data["profile"][key].where(mask_time[key]).copy()
            )

            output_data["time_mask"][key] = mask_time[key]

            output_data["time_info"][key] = (
                output_data["time_info"][key].where(mask_time[key]).copy()
            )

            output_data["shots"][key] = (
                output_data["shots"][key].where(mask_time[key]).copy()
            )

    return output_data


def get_overflow_mask(sig, acquisition_mode, daq_range=None, max_count=2.0 ** 15):
    """
    Build a lazy 3D photon-overflow mask.

    NaN values in sig are naturally treated as False by the comparison
    sig >= max_count, so no fillna() is needed.

    Parameters
    ----------
    sig : xr.DataArray
        Signal with dims time, channel, bins.
    acquisition_mode : xr.DataArray
        Broadcast/aligned acquisition mode. Photon channels use "p".
    daq_range : xr.DataArray, optional
        Kept in the signature for compatibility. Not used because the current
        correction policy only handles photon overflows.
    max_count : float
        Maximum allowed photon count.
    """

    mask = (acquisition_mode == "p") & (sig >= max_count)

    # Keep the same dimension order as the input signal.
    return mask.transpose(*sig.dims)


def _adjacent_overflow_time_mask(
    mask: xr.DataArray,
    max_adjacent_overflows: int,
) -> xr.DataArray:
    """
    Lazily detect time profiles with a run of adjacent overflow bins.

    Returns a 1D time mask. True means at least one channel in that time
    profile has max_adjacent_overflows consecutive overflow bins.

    This avoids xarray.rolling(...).sum(), which can be expensive on large
    Dask-backed arrays. It assumes bins are intact or at least that a run does
    not need to be detected across a Dask chunk boundary. With your preferred
    chunking (bins=-1), this is safe.
    """

    if max_adjacent_overflows is None or max_adjacent_overflows <= 1:
        return mask.any(dim=("bins", "channel"))

    run = mask
    for shift in range(1, max_adjacent_overflows):
        run = run & mask.shift(bins=-shift, fill_value=False)

    return run.any(dim=("bins", "channel"))


def _fill_overflows_nearest_bins(sig: xr.DataArray, mask: xr.DataArray) -> xr.DataArray:
    """
    Faster but cruder overflow replacement using nearest valid bins.

    Overflow bins are first set to NaN, then filled forward and backward along
    bins. Linear interpolation is usually scientifically preferable, but this
    can be useful if interpolation is too slow.
    """

    sig_nan = sig.where(~mask)
    return sig_nan.ffill("bins").bfill("bins")


def overflow_method_0(mask, filename):
    """
    trim_overflows = 0.

    If no overflows exist, do nothing.
    If overflows exist, print filename, channel, and overflow bins, then raise.

    The full 3D mask is never computed for all data. First, small reduced
    masks are computed. Then only the affected time/channel subset is computed
    for exact bin reporting.
    """

    mask_time, mask_ch_time = dask.compute(
        mask.any(dim=("bins", "channel")),
        mask.any(dim="bins"),
    )

    if not bool(mask_time.any().item()):
        return xr.ones_like(mask["time"], dtype=bool)

    print("")
    print("-- At least one bin with an overflow was detected ")
    print("-- Please revise the following bins: ")

    time = mask.time.values
    channel = mask.channel.values
    bins = mask.bins.values

    bad_times = time[mask_time.values]
    bad_channels = channel[mask_ch_time.any(dim="time").values]

    # Compute only the affected time/channel subset.
    mask_bad = mask.sel(time=bad_times, channel=bad_channels).compute()

    for t in bad_times:
        ch_ovf = bad_channels[mask_bad.sel(time=t).any(dim="bins").values]

        for ch in ch_ovf:
            bins_ovf = bins[mask_bad.sel(time=t, channel=ch).values]
            print(
                f"    file: {filename.loc[t].values} | "
                f"ch: {ch} | bins: {bins_ovf}"
            )

    raise DataOverflowError(
        "trim_overflows = 0 -> Overflows detected! In order to continue with "
        "an automated overflow removal use the trim_overflow argument with "
        "value 1 or 2 (default is 0) "
    )


def overflow_method_1(sig, shots, time_info, mask, filename):
    """
    trim_overflows = 1.

    Remove full time profiles if there is at least one overflow in any
    channel/bin. If there are no overflows, do nothing.
    """

    mask_time = mask.any(dim=("bins", "channel")).compute()

    if not bool(mask_time.any().item()):
        valid_time = xr.ones_like(sig["time"], dtype=bool)
        return sig, shots, time_info, valid_time

    time = mask_time.time.values
    time_ovf = time[mask_time.values]
    time_cor = time[~mask_time.values]

    print("")
    print(
        f"-- Warning: trim_overflows = 1 -> Removing {time_ovf.size} "
        "profiles with at least one overflow:"
    )

    for t in time_ovf:
        print(f"    {filename.loc[t].values} ")

    sig = sig.sel(time=time_cor)
    shots = shots.sel(time=time_cor)
    time_info = time_info.sel(time=time_cor)

    valid_time = xr.ones_like(sig["time"], dtype=bool)

    return sig, shots, time_info, valid_time


def overflow_method_2(
    sig,
    shots,
    time_info,
    mask,
    filename,
    max_adjacent_overflows,
    max_overflows_per_profile=100,
    fill_method="linear",
):
    """
    trim_overflows = 2.

    If there are no overflows, do nothing.
    If isolated/safe overflows exist, replace them lazily.
    If too many adjacent overflow bins exist in a profile, remove that profile.
    If too many total overflow bins exist in any time/channel profile, raise.

    fill_method can be:
    - "linear": use interpolate_na along bins.
    - "nearest": use ffill/bfill along bins.
    """

    # First do the cheapest useful check. If no overflow exists, skip all other
    # diagnostics and avoid ovfs/adjacent-run scans.
    mask_time = mask.any(dim=("bins", "channel")).compute()

    if not bool(mask_time.any().item()):
        valid_time = xr.ones_like(sig["time"], dtype=bool)
        return sig, shots, time_info, valid_time

    # Overflows exist. Compute the method-2 safety diagnostics together.
    ovfs_lazy = mask.sum(dim="bins")
    mask_adjacent_lazy = _adjacent_overflow_time_mask(
        mask=mask,
        max_adjacent_overflows=max_adjacent_overflows,
    )

    ovfs, mask_adjacent = dask.compute(ovfs_lazy, mask_adjacent_lazy)

    if bool((ovfs > max_overflows_per_profile).any().item()):
        print("")
        raise DataOverflowError(
            f"More than {max_overflows_per_profile} overflowed bins "
            "encountered in single time/channel profiles. Filling is too "
            "risky, please revise the input files or consider setting "
            "trim_overflows = 1."
        )

    time = mask_time.time.values
    time_ovf = time[mask_time.values]
    overflow_count = int(ovfs.sum().item())

    if bool(mask_adjacent.any().item()):
        time_adjacent = time[mask_adjacent.values]
        time_cor = time[~mask_adjacent.values]

        print("")
        print(
            f"-- Warning: At least {max_adjacent_overflows} adjacent "
            "overflow bins were encountered in a single profile. Filling "
            "is too risky, these files will be removed:"
        )

        for t in time_adjacent:
            print(f"    {filename.loc[t].values} ")

        sig = sig.sel(time=time_cor)
        shots = shots.sel(time=time_cor)
        time_info = time_info.sel(time=time_cor)
        mask = mask.sel(time=time_cor)
        mask_time = mask_time.sel(time=time_cor)

        # Print/correct only the remaining isolated-overflow profiles.
        time = mask_time.time.values
        time_ovf = time[mask_time.values]

    print("")
    print(
        f"-- Warning: trim_overflows = 2 -> Replacing overflows in "
        f"{time_ovf.size} profiles:"
    )

    original_chunks = sig.chunksizes if _has_dask_chunks(sig) else None

    # interpolate_na/ffill/bfill along bins need the interpolation/fill axis
    # available inside each block. With the preferred reader chunking, bins
    # should already be one chunk. Only rechunk if necessary.
    if not _bins_are_single_chunk(sig):
        sig = sig.chunk({"bins": -1})
        mask = mask.chunk({"bins": -1})

    if fill_method == "nearest":
        sig = _fill_overflows_nearest_bins(sig, mask)
        fill_text = "nearest-bin filling"
    elif fill_method == "linear":
        sig = sig.where(~mask).interpolate_na(dim="bins", method="linear")
        fill_text = "linear interpolation across the bins"
    else:
        raise ValueError(
            "fill_method must be either 'linear' or 'nearest', "
            f"got {fill_method!r}"
        )

    if original_chunks is not None:
        sig = sig.chunk(original_chunks)

    print(f"{overflow_count} overflows have been replaced by {fill_text}\n")

    valid_time = xr.ones_like(sig["time"], dtype=bool)

    return sig, shots, time_info, valid_time


def overflow_method_3(sig, shots, time_info):
    """
    trim_overflows = 3.

    Do nothing even if overflows are present. No overflow check is performed.
    """

    valid_time = xr.ones_like(sig["time"], dtype=bool)
    return sig, shots, time_info, valid_time


def compute_check_for_overflows(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    """
    Detect and handle lidar signal overflows while preserving lazy profile data.

    The full profile-sized overflow mask remains lazy. Method-specific eager
    computations are reduced to the minimum needed by the selected method.
    """

    output_data = copy.deepcopy(input_data)

    method = processing_info["caller_info"]["trim_overflows"]
    max_adjacent_overflows = processing_info["caller_info"].get(
        "max_adjacent_overflows",
        5,
    )
    max_overflows_per_profile = processing_info["caller_info"].get(
        "max_overflows_per_profile",
        100,
    )
    fill_method = processing_info["caller_info"].get(
        "overflow_fill_method",
        "linear",
    )

    profiles = output_data["profile"]
    time_info = output_data["time_info"]
    shots = output_data["shots"]
    channel_info = output_data["channel_info"]

    qa_tests = list(profiles.keys())

    # Start with full-True masks for all measurements. These are tiny 1D time
    # masks, not full profile masks.
    mask_time = {
        key: xr.ones_like(profiles[key]["time"], dtype=bool)
        for key in qa_tests
    }

    mask_ovf = {}

    for key in qa_tests:

        print_subsection(f"{key} dataset")

        # trim_overflows = 3 means: do nothing, do not even check.
        if method == 3:
            profiles[key], shots[key], time_info[key], mask_time[key] = \
                overflow_method_3(
                    sig=profiles[key],
                    shots=shots[key],
                    time_info=time_info[key],
                )
            continue

        # Keep the original broadcast_like behavior. This preserves alignment
        # with profile dimensions/order before comparisons are made.
        aq_mode = (
            channel_info[key]
            .sel({"parameters": "acquisition_mode"})
            .drop_vars("parameters")
            .broadcast_like(profiles[key])
        )

        daq_range = (
            channel_info[key]
            .sel({"parameters": "data_acquisition_range"})
            .drop_vars("parameters")
            .broadcast_like(profiles[key])
        )

        filename = time_info[key].sel({"parameters": "filename"}).compute()

        # Lazy 3D overflow mask.
        mask_ovf[key] = get_overflow_mask(
            profiles[key],
            acquisition_mode=aq_mode,
            daq_range=daq_range,
        )

        if method == 0:
            mask_time[key] = overflow_method_0(
                mask=mask_ovf[key],
                filename=filename,
            )

        elif method == 1:
            profiles[key], shots[key], time_info[key], mask_time[key] = \
                overflow_method_1(
                    sig=profiles[key],
                    shots=shots[key],
                    time_info=time_info[key],
                    mask=mask_ovf[key],
                    filename=filename,
                )

        elif method == 2:
            profiles[key], shots[key], time_info[key], mask_time[key] = \
                overflow_method_2(
                    sig=profiles[key],
                    shots=shots[key],
                    time_info=time_info[key],
                    mask=mask_ovf[key],
                    filename=filename,
                    max_adjacent_overflows=max_adjacent_overflows,
                    max_overflows_per_profile=max_overflows_per_profile,
                    fill_method=fill_method,
                )

        else:
            raise ValueError(
                "trim_overflows must be one of 0, 1, 2, or 3, "
                f"got {method!r}"
            )

    output_data["time_mask"] = mask_time
    output_data["profile_mask"] = mask_ovf

    apply_time_mask(qa_tests, mask_time, output_data)

    if output_data["profile"] == {}:
        endpoint(5)

    return output_data
