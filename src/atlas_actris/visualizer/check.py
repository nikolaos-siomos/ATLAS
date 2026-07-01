#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:27:01 2022

@author: nick
"""

import numpy as np
import xarray as xr
from utils.error_classes import CustomWarning

def check_channels(all_channels, settings):

    sel_channels = settings["select_channels"]
    exclude_wavelength = settings["exclude_wavelength"]
    exclude_telescope_type = settings["exclude_telescope_type"]
    exclude_channel_type = settings["exclude_channel_type"]
    exclude_acquisition_mode = settings["exclude_acquisition_mode"]
    exclude_channel_subtype = settings["exclude_channel_subtype"]

    all_channels = np.asarray(all_channels, dtype=str)

    # Include channels
    if len(sel_channels) == 0:
        channels = all_channels
    else:
        channels = np.asarray(
            [ch for ch in sel_channels if ch in all_channels],
            dtype=str,
        )

    # Missing channels
    missing_channels = [ch for ch in sel_channels if ch not in all_channels]

    if missing_channels:
        CustomWarning(
            "Channels provided in select_channels do not exist: "
            f"{missing_channels}\nPlease select one of: "
            f"{all_channels.tolist()}"
        )
        print()

    mask = np.array([
        ch[:4] not in exclude_wavelength and
        ch[4] not in exclude_telescope_type and
        ch[5] not in exclude_channel_type and
        ch[6] not in exclude_acquisition_mode and
        ch[7] not in exclude_channel_subtype
        for ch in channels
    ], dtype=bool)

    if not np.any(mask):
        CustomWarning(
            "The provided channel filtering arguments are too strict "
            "and exclude all channels. Please revise the following arguments: "
            "exclude_wavelength, exclude_telescope_type, exclude_channel_type, "
            "exclude_acquisition_mode, exclude_channel_subtype, select_channels"
        )
        print()


    channels = channels[mask]

    return channels

def check_vldr_pairs(pol_cal_ratio, pol_cal_info, settings=None):
    """
    Return VLDR pair ids and their reflected/transmitted source channels.

    The preferred selection comes from pol_cal_info ratio_type == 'vldr'.  As a
    fallback, pair ids whose sixth character is 'v' are accepted, but ch_r/ch_t
    still need to be available in pol_cal_info because the vertical scale is
    channel-based.

    Supports the same include/exclude option names used by check_pairs:
        include_pairs
        exclude_wavelength
        exclude_telescope_type
        exclude_pairl_type
        exclude_acquisition_mode
        exclude_pair_subtype

    The existing _select_vldr_pairs options are also kept:
        pairs
        vldr_pairs
        channels
    """

    settings = settings or {}

    if "pair" not in pol_cal_ratio.dims:
        raise ValueError(
            "The VLDR quicklook input must have a 'pair' dimension. "
            f"Found dimensions: {pol_cal_ratio.dims}."
        )

    pair_values = [str(pair) for pair in pol_cal_ratio.pair.values]

    selected = []

    if (
        pol_cal_info is not None
        and "parameters" in pol_cal_info.dims
        and "pair" in pol_cal_info.dims
        and "ratio_type" in pol_cal_info.parameters.values
    ):
        ratio_type = pol_cal_info.sel(parameters="ratio_type")
        for pair, value in zip(ratio_type.pair.values, ratio_type.values):
            pair_str = str(pair)
            if pair_str not in pair_values:
                continue
            if str(value).strip().lower() == "vldr":
                selected.append(pair_str)

    if len(selected) == 0:
        selected = [pair for pair in pair_values if len(pair) == 8 and pair[5] == "v"]

    # -------------------------------------------------------------------------
    # Include-pair filtering
    # -------------------------------------------------------------------------
    include_pairs = settings.get("include_pairs", [])

    if include_pairs is None:
        include_pairs = []

    include_pairs = [str(pair) for pair in include_pairs]

    if len(include_pairs) > 0:
        available_pairs = set(selected)
        missing_pairs = [pair for pair in include_pairs if pair not in available_pairs]

        if len(missing_pairs) > 0:
            CustomWarning(
                "Pairs provided in include_pairs do not exist or are "
                "not available VLDR pairs:\n"
                f"{missing_pairs}\n"
                "Please select one of:\n"
                f"{selected}"
            )

        selected = [pair for pair in include_pairs if pair in available_pairs]

    else:
        # Keep the previous _select_vldr_pairs behavior for pairs/vldr_pairs.
        requested_pairs = settings.get("pairs", settings.get("vldr_pairs", []))

        if requested_pairs is None:
            requested_pairs = []

        requested_pairs = [str(pair) for pair in requested_pairs]

        if len(requested_pairs) > 0:
            requested_pairs = set(requested_pairs)
            selected = [pair for pair in selected if pair in requested_pairs]

    # -------------------------------------------------------------------------
    # Exclude filtering copied from check_pairs logic
    # -------------------------------------------------------------------------
    exclude_wavelength = settings.get("exclude_wavelength", [])
    exclude_telescope_type = settings.get("exclude_telescope_type", [])
    exclude_pair_type = settings.get("exclude_pairl_type", [])
    exclude_acquisition_mode = settings.get("exclude_acquisition_mode", [])
    exclude_pair_subtype = settings.get("exclude_pair_subtype", [])

    if exclude_wavelength is None:
        exclude_wavelength = []
    if exclude_telescope_type is None:
        exclude_telescope_type = []
    if exclude_pair_type is None:
        exclude_pair_type = []
    if exclude_acquisition_mode is None:
        exclude_acquisition_mode = []
    if exclude_pair_subtype is None:
        exclude_pair_subtype = []

    exclude_wavelength = {str(value) for value in exclude_wavelength}
    exclude_telescope_type = {str(value) for value in exclude_telescope_type}
    exclude_pair_type = {str(value) for value in exclude_pair_type}
    exclude_acquisition_mode = {str(value) for value in exclude_acquisition_mode}
    exclude_pair_subtype = {str(value) for value in exclude_pair_subtype}

    filtered = []

    for pair in selected:
        pair = str(pair)

        # check_pairs assumes 8-character pair IDs.  Here we keep malformed
        # pairs instead of crashing, so the existing ch_r/ch_t validation below
        # can still provide the more relevant error.
        if len(pair) < 8:
            filtered.append(pair)
            continue

        keep_pair = (
            pair[:4] not in exclude_wavelength and
            pair[4] not in exclude_telescope_type and
            pair[5] not in exclude_pair_type and
            pair[6] not in exclude_acquisition_mode and
            pair[7] not in exclude_pair_subtype
        )

        if keep_pair:
            filtered.append(pair)

    if len(selected) > 0 and len(filtered) == 0:
        print(
            "Warning: The provided pair filtering arguments are too strict "
            "and exclude all VLDR pairs. Please revise the following arguments: "
            "include_pairs, exclude_wavelength, exclude_telescope_type, "
            "exclude_pairl_type, exclude_acquisition_mode, "
            "exclude_pair_subtype, pairs, vldr_pairs."
        )

    selected = filtered

    # -------------------------------------------------------------------------
    # Existing channel filtering
    # -------------------------------------------------------------------------
    requested_channels = settings.get("channels", [])

    if requested_channels is None:
        requested_channels = []

    requested_channels = {str(channel) for channel in requested_channels}

    pair_records = []

    for pair in selected:
        ch_r = _get_info_value(pol_cal_info, "ch_r", pair)
        ch_t = _get_info_value(pol_cal_info, "ch_t", pair)

        if ch_r is None:
            raise ValueError(
                f"Cannot plot VLDR pair {pair}: missing 'ch_r' in pol_cal_info. "
                "The vertical scale is channel-based, so ch_r is needed."
            )

        if len(requested_channels) > 0:
            channels_for_pair = {str(ch_r)}

            if ch_t is not None:
                channels_for_pair.add(str(ch_t))

            if channels_for_pair.isdisjoint(requested_channels):
                continue

        pair_records.append(
            {
                "pair": pair,
                "ch_r": ch_r,
                "ch_t": ch_t,
            }
        )

    return pair_records

def _get_info_value(pol_cal_info, parameter, pair, default=None):
    """Safely read one value from pol_cal_info(parameters, pair)."""

    if pol_cal_info is None:
        return default

    if "parameters" not in pol_cal_info.dims or "pair" not in pol_cal_info.dims:
        return default

    if parameter not in pol_cal_info.parameters.values:
        return default

    if pair not in pol_cal_info.pair.values:
        return default

    value = pol_cal_info.sel(parameters=parameter, pair=pair).values

    value = np.asarray(value)
    if value.size == 0:
        return default

    value = value.item()

    if not _is_valid_info_value(value):
        return default

    return str(value)

def _is_valid_info_value(value):
    """Return False for empty, None-like, and NaN-like metadata values."""

    if value is None:
        return False

    try:
        if bool(np.isnan(value)):
            return False
    except Exception:
        pass

    value_str = str(value).strip()

    return value_str.lower() not in ["", "nan", "none"]


def find_rt_channels(ch_r, ch_t, channels):
    
    if ch_r == None or ch_t == None:
        channels_r = []
        channels_t = []
        ch_r_all = np.array([ch for ch in channels if ch[7] == 'r'])
        ch_t_all = np.array([ch for ch in channels if ch[7] == 't'])
        if len(ch_r_all) == 0:
            print("-- Warning: No relfected channels were detected. Please make sure that the channel_subtype in set correctly in the configuration file")
        if len(ch_t_all) == 0:
            print("-- Warning: No transmitted channels were detected. Please make sure that the channel_subtype in set correctly in the configuration file")
        
        for ch_r_i in ch_r_all:
            for ch_t_i in ch_t_all:
                if ch_r_i[4]  == ch_t_i[4] and ch_r_i[6]  == ch_t_i[6] and \
                    ch_r_i[:4]  == ch_r_i[:4]:
                        channels_r.extend([ch_r_i])
                        channels_t.extend([ch_t_i])
    else:
        channels_r = ch_r
        channels_t = ch_t
        
    return(channels_r, channels_t)

def _valid_channel_value(value):
    if value is None:
        return False

    try:
        if bool(np.isnan(value)):
            return False
    except Exception:
        pass

    value = str(value).strip()

    return value.lower() not in ["", "nan", "none"]


def check_pairs(pol_cal_info, caller_info):
    """Filter pol_cal_info pairs using the ch_r/ch_t source channels."""

    if not isinstance(pol_cal_info, xr.DataArray):
        return pol_cal_info

    if "pair" not in pol_cal_info.dims or "parameters" not in pol_cal_info.dims:
        return pol_cal_info

    if "ch_r" not in pol_cal_info.parameters.values:
        return pol_cal_info

    channel_values = []

    for parameter in ["ch_r", "ch_t"]:
        if parameter not in pol_cal_info.parameters.values:
            continue

        values = np.asarray(
            pol_cal_info.sel(parameters=parameter).values,
            dtype=object,
        ).ravel()

        channel_values.extend(
            str(value).strip()
            for value in values
            if _valid_channel_value(value)
        )

    if len(channel_values) == 0:
        return pol_cal_info

    allowed_channels = set(check_channels(channel_values, caller_info))
    keep_pairs = []

    for pair in pol_cal_info.pair.values:
        pair_channels = []

        for parameter in ["ch_r", "ch_t"]:
            if parameter not in pol_cal_info.parameters.values:
                continue

            value = pol_cal_info.sel(parameters=parameter, pair=pair).values
            value = np.asarray(value).item()

            if _valid_channel_value(value):
                pair_channels.append(str(value).strip())

        if len(pair_channels) > 0 and all(
            channel in allowed_channels for channel in pair_channels
        ):
            keep_pairs.append(pair)

    pol_cal_info = pol_cal_info.sel({"pair": keep_pairs})

    return pol_cal_info
