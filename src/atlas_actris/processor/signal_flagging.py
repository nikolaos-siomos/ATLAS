#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lazy saturation flagging for ATLAS profiles.

Optimized version:

1. The full time x channel x bins saturation mask is kept lazy and is built
   only if requested through store_saturation_profile_mask.
2. Channel-level warning flags are computed by first selecting analog/photon
   channels. This avoids applying analog checks to photon channels and photon
   checks to analog channels.
3. Photon saturation uses multiplication instead of division:
      sig > photon_max_counts * shots
   instead of:
      sig / shots > photon_max_counts
4. Warning output is grouped by saturation type, instead of one warning per
   channel.
5. NaN handling relies on comparison behavior: comparisons with NaN evaluate
   to False, so fillna(0.) is not needed for these masks.
"""

import copy
from typing import Any, Dict, Iterable, Sequence

import dask
import xarray as xr

from utils.error_classes import CustomWarning
from utils.printouts import print_subsection


def _reduce_dims_except_channel(da: xr.DataArray) -> Sequence[str]:
    """Return all dimensions except channel, preserving dimension order."""
    return tuple(dim for dim in da.dims if dim != "channel")


def _empty_channel_bool_like(channel_template: xr.DataArray) -> xr.DataArray:
    """Return a computed all-False boolean DataArray on the channel dimension."""
    return xr.full_like(channel_template, False, dtype=bool)


def _channel_values(mask_ch: xr.DataArray) -> list:
    """
    Convert a computed 1D channel boolean mask to a list of channel labels.

    Parameters
    ----------
    mask_ch : xr.DataArray
        Boolean DataArray with dimension channel. It is expected to be already
        computed and small.
    """
    if "channel" not in mask_ch.dims:
        return []

    return list(mask_ch["channel"].values[mask_ch.values])


def _format_channel_list(channels: Iterable[Any], indent: str = "    ") -> str:
    """Format channel labels as a compact indented list."""
    channels = list(channels)

    if len(channels) == 0:
        return ""

    return "\n".join(f"{indent}{ch}" for ch in channels)


def _print_grouped_saturation_warnings(
    analog_channels: Iterable[Any],
    photon_channels: Iterable[Any],
    analog_fraction: float,
    photon_max_counts: float,
) -> None:
    """Print grouped saturation warnings."""
    analog_channels = list(analog_channels)
    photon_channels = list(photon_channels)

    if len(analog_channels) == 0 and len(photon_channels) == 0:
        return

    print("")
    CustomWarning("Saturation was detected in one or more channels.")

    if len(analog_channels) > 0:
        print("")
        print(
            f"-- Analog channels with signal values above "
            f"{100.0 * analog_fraction:g}% of the data acquisition range:"
        )
        print(_format_channel_list(analog_channels))

    if len(photon_channels) > 0:
        print("")
        print(
            f"-- Photon channels with signal values above "
            f"{photon_max_counts:g} count per bin per shot:"
        )
        print(_format_channel_list(photon_channels))

    print("")


def _channels_by_mode(acquisition_mode: xr.DataArray) -> Dict[str, Any]:
    """
    Return channel labels split by acquisition mode.

    acquisition_mode is expected to be a small 1D DataArray with dim channel.
    """
    if "channel" not in acquisition_mode.dims:
        raise ValueError("acquisition_mode must contain a 'channel' dimension")

    modes = acquisition_mode.values
    channels = acquisition_mode["channel"].values

    return {
        "analog": channels[modes == "a"],
        "photon": channels[modes == "p"],
    }


def build_saturation_profile_mask(
    sig: xr.DataArray,
    shots: xr.DataArray,
    acquisition_mode: xr.DataArray,
    daq_range: xr.DataArray,
    analog_fraction: float = 0.80,
    photon_max_counts: float = 1.0,
) -> xr.DataArray:
    """
    Build a lazy full profile-sized saturation mask.

    This function does not compute the mask. It is only used if the caller
    asks to store output_data["profile_mask"][key].
    """
    max_daq_range = analog_fraction * daq_range

    mask_analog_saturated = (
        (acquisition_mode == "a")
        & (sig > max_daq_range)
    )

    # Avoid division over the full signal. This is equivalent to
    # sig / shots > photon_max_counts when shots are positive.
    mask_photon_saturated = (
        (acquisition_mode == "p")
        & (sig > photon_max_counts * shots)
    )

    mask_saturated = mask_analog_saturated | mask_photon_saturated

    return mask_saturated.transpose(*sig.dims)


def compute_saturation_channel_flags(
    sig: xr.DataArray,
    shots: xr.DataArray,
    acquisition_mode: xr.DataArray,
    daq_range: xr.DataArray,
    analog_fraction: float = 0.80,
    photon_max_counts: float = 1.0,
) -> Dict[str, xr.DataArray]:
    """
    Compute only channel-level saturation flags.

    This avoids computing the full 3D saturation mask. It also avoids applying
    analog checks to photon channels and photon checks to analog channels.

    Returned arrays are computed, small boolean DataArrays with dim channel
    and the same channel coordinates as acquisition_mode. Channels that do not
    belong to the corresponding acquisition type are False.
    """
    mode_channels = _channels_by_mode(acquisition_mode)
    analog_channels = mode_channels["analog"]
    photon_channels = mode_channels["photon"]

    analog_full = _empty_channel_bool_like(acquisition_mode)
    photon_full = _empty_channel_bool_like(acquisition_mode)

    lazy_results = []
    result_names = []

    if len(analog_channels) > 0:
        sig_a = sig.sel(channel=analog_channels)
        daq_a = daq_range.sel(channel=analog_channels)
        reduce_dims_a = _reduce_dims_except_channel(sig_a)

        analog_ch_lazy = (
            sig_a > analog_fraction * daq_a
        ).any(dim=reduce_dims_a)

        lazy_results.append(analog_ch_lazy)
        result_names.append("analog")

    if len(photon_channels) > 0:
        sig_p = sig.sel(channel=photon_channels)
        shots_p = shots.sel(channel=photon_channels)
        reduce_dims_p = _reduce_dims_except_channel(sig_p)

        # Avoid division over the full signal.
        photon_ch_lazy = (
            sig_p > photon_max_counts * shots_p
        ).any(dim=reduce_dims_p)

        lazy_results.append(photon_ch_lazy)
        result_names.append("photon")

    computed_values = dask.compute(*lazy_results) if lazy_results else []
    computed = dict(zip(result_names, computed_values))

    if "analog" in computed:
        analog_full.loc[dict(channel=analog_channels)] = computed["analog"]

    if "photon" in computed:
        photon_full.loc[dict(channel=photon_channels)] = computed["photon"]

    return {
        "analog": analog_full,
        "photon": photon_full,
    }


def compute_detect_saturation(
    caller_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    """
    Detect lidar profiles with saturated analog/photon signal levels.

    The full saturation mask is kept lazy in output_data["profile_mask"] only
    when store_saturation_profile_mask=True. Only small channel-level masks are
    computed eagerly for warnings.

    Parameters
    ----------
    caller_info : dict
        Runtime settings. Optional keys used here:
        - saturation_analog_fraction: default 0.80
        - saturation_photon_max_counts: default 1.0
        - store_saturation_profile_mask: default True

    input_data : dict
        Nested ATLAS data dictionary. Required top-level keys are profile,
        shots, and channel_info.
    """

    output_data = copy.deepcopy(input_data)
    output_data.setdefault("profile_mask", {})
    output_data.setdefault("saturation_channel_mask", {})

    profiles = output_data["profile"]
    shots = output_data["shots"]
    channel_info = output_data["channel_info"]

    qa_tests = list(profiles.keys())

    analog_fraction = caller_info.get("saturation_analog_fraction", 0.80)
    photon_max_counts = caller_info.get("saturation_photon_max_counts", 1.0)
    store_profile_mask = caller_info.get("store_saturation_profile_mask", True)

    for key in qa_tests:

        print_subsection(f"{key} dataset")

        aq_mode = (
            channel_info[key]
            .sel({"parameters": "acquisition_mode"})
            .drop_vars("parameters")
        )

        daq_range = (
            channel_info[key]
            .sel({"parameters": "data_acquisition_range"})
            .drop_vars("parameters")
        )

        channel_flags = compute_saturation_channel_flags(
            sig=profiles[key],
            shots=shots[key],
            acquisition_mode=aq_mode,
            daq_range=daq_range,
            analog_fraction=analog_fraction,
            photon_max_counts=photon_max_counts,
        )

        analog_channels = _channel_values(channel_flags["analog"])
        photon_channels = _channel_values(channel_flags["photon"])

        _print_grouped_saturation_warnings(
            analog_channels=analog_channels,
            photon_channels=photon_channels,
            analog_fraction=analog_fraction,
            photon_max_counts=photon_max_counts,
        )

        output_data["saturation_channel_mask"][key] = channel_flags

        if store_profile_mask:
            output_data["profile_mask"][key] = build_saturation_profile_mask(
                sig=profiles[key],
                shots=shots[key],
                acquisition_mode=aq_mode,
                daq_range=daq_range,
                analog_fraction=analog_fraction,
                photon_max_counts=photon_max_counts,
            )

    return output_data
