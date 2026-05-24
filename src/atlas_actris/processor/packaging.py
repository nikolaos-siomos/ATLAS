#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:49:05 2026

@author: nikos
"""


import xarray as xr
from processor.selectors import get_source_map
from utils.dataarray_utils import shallow_copy
from utils.printouts import print_entry

packing_map = {
    "tlc_qua": [
        "tlc_north",
        "tlc_east",
        "tlc_south",
        "tlc_west",
    ],
    "tlc_rin": [
        "tlc_outer",
        "tlc_inner",
    ],
    "pcb": [
        "pcb_p45",
        "pcb_m45",
    ],
}

def combine_QA_pack(input_data):
    
    output_data = shallow_copy(input_data)

    source_map = get_source_map("dummy")

    for new_key, pack in packing_map.items():
        
        for source in source_map.keys():
    
            da = input_data[source]
    
            pack_list = []
    
            for key in pack:
                if key in da:
                    pack_list.append(da[key])
    
            if len(pack_list) > 0:
                
                if type(pack_list[0]) is type(xr.DataArray()):
                    if "time" in pack_list[0].dims:
                        output_data[source][new_key] = xr.concat(
                            pack_list,
                            dim="time",
                        ).sortby("time")
                    else:
                        output_data[source][new_key] = pack_list[0]
                  
    print_entry("Packaging complete")

    return output_data

def to_plain_value(x):
    """Convert an xarray scalar/vector to a plain Python value/list."""
    values = x.values
    if values.shape == ():
        return values.item()
    return values.tolist()

def collect_metadata(data_arrays, atlas_channel_id):
    """
    Collect metadata parameters for one channel.

    Parameters
    ----------
    data_arrays : dict[str, xarray.DataArray]
        Dictionary containing DataArrays.
    channel_id :
        Channel label used for slicing arrays with a "channel" dimension.

    Returns
    -------
    metadata : dict
        Dictionary containing parameter values for the selected channel.
    """

    metadata = {
        "atlas_channel_id": atlas_channel_id
    }
    
    for array_name, da in data_arrays.items():

        # Skip nested dictionaries, scalars, lists, etc.
        if not isinstance(da, xr.DataArray):
            continue
        
        # Ignore arrays without a parameters dimension
        if "parameters" not in da.dims:
            continue

        selected = da

        # If this array depends on channel, select the requested channel
        if "channel" in selected.dims:
            selected = selected.sel(channel=atlas_channel_id)

        # If this array depends on time, keep first and last time entry
        if "time" in selected.dims:
            time_slices = {
                "first": selected.isel(time=0),
                "last": selected.isel(time=-1),
            }
        else:
            time_slices = {
                "": selected
            }

        # Save every parameter
        for time_label, sliced in time_slices.items():
            for parameter_name in sliced["parameters"].values:

                value = sliced.sel(parameters=parameter_name)
                plain_value = to_plain_value(value)

                key = str(parameter_name)

                # Add suffix for first/last time values
                if time_label:
                    key = f"{key}_{time_label}"

                # Avoid overwriting if the same parameter name appears twice
                if key in metadata:
                    key = f"{array_name}_{key}"

                metadata[key] = plain_value

    return dict(sorted(metadata.items()))
