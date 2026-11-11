import xarray as xr
from utils.printouts import endpoint, print_header
from visualizer.check import check_channels, check_pairs
from processor.channel_validation import drop_measurement, warn_once


def filter_channels(caller_info, profiles, metadata):
    """Filter channels per measurement while keeping signal metadata aligned."""
    print_header("Filtering out channels")
    for key in list(profiles):
        arr = profiles[key]
        channels = check_channels(arr.channel.values, caller_info)
        if len(channels) == 0:
            warn_once(caller_info, f"{key}: no selected channels; measurement skipped.")
            drop_measurement(key, profiles, metadata, caller_info)
            continue
        profiles[key] = arr.sel(channel=channels)
        for name, group in metadata.items():
            if not isinstance(group, dict) or key not in group:
                continue
            value = group[key]
            if not isinstance(value, xr.DataArray):
                continue
            if "channel" in value.dims:
                group[key] = value.sel(channel=channels)
            elif name == "pol_cal_info" and "pair" in value.dims:
                group[key] = check_pairs(value, caller_info)
    if not profiles:
        endpoint(6)
    return profiles, metadata
