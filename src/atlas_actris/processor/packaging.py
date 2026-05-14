#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:49:05 2026

@author: nikos
"""


import xarray as xr
from processor.selectors import get_source_map
from utils.dataarray_utils import shallow_copy
from helper_functions.printouts import print_entry

packing_map = {
    "tlc": [
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
                if "time" in pack_list[0].dims:
                    output_data[source][new_key] = xr.concat(
                        pack_list,
                        dim="time",
                    ).sortby("time")
                else:
                    output_data[source][new_key] = pack_list[0]
                  
    print_entry("Packaging complete")

    return output_data