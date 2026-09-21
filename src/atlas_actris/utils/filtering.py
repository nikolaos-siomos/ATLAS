import numpy as np
import xarray as xr
from utils.printouts import endpoint
from visualizer.check import check_channels, check_pairs
from utils.printouts import print_header, print_entry


def filter_channels(caller_info, profiles, metadata):
    """
    Apply channel filtering to profiles and metadata.

    profiles is expected to be a dictionary of xarray objects.
    metadata is expected to be a dictionary whose values are dictionaries.
    Non-xarray values are returned unchanged.
    """

    print_header("Filtering out channels")
    
    for qa_test, arr in profiles.items():
        if isinstance(arr, xr.DataArray):
            if "channel" in arr.dims:
                channels = arr.channel.values
                filtered_channels = check_channels(channels, caller_info)
                profiles[qa_test] = arr.sel({"channel": filtered_channels})
                
                # if profiles[qa_test].channel.size == 0.:
                #     endpoint(6)

    for key, val in metadata.items():
        if isinstance(val, dict):
            for qa_test, arr in val.items():
                if isinstance(arr, xr.DataArray):
                    if "channel" in arr.dims:
                        channels = arr.channel.values
                        filtered_channels = check_channels(channels, caller_info)

                        metadata[key][qa_test] = arr.sel({"channel": filtered_channels})

                        if metadata[key][qa_test].channel.size == 0.:
                            endpoint(6)

                    elif key == "pol_cal_info" and "pair" in arr.dims:
                        metadata[key][qa_test] = check_pairs(
                            arr,
                            caller_info,
                        )
                        print(key, qa_test)
                        # if metadata[key][qa_test].pair.size == 0.:
                        #     endpoint(6)

    return profiles, metadata
