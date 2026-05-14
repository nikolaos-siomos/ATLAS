"""
@authors: N. Siomos, P. Paschou, Ioannis Binietoglou 
based on SULA project (https://react-gitlab.space.noa.gr/ReACT/eve/data-processing)
and also on https://gitlab.com/ioannis_binietoglou/lidar-processing/

Processing routines for signals 

=================================
Signal in 3D xarray dataset with dimensions [time, channel, bin/range]

Fucntions
 -- average_by_time: Average signals across the timeframes
 -- background_calculation: Calculates the solar background per timeframe and channel 
 -- background_correction: Performs the background correction on signals
 -- dark_correction: Removes the dark signals from the normal ananlog signals
 -- dead time correction: Performs the dead time correction onphoton channels
 -- detect_saturation: Identifies regions where signals are saturated
 -- height_calculation: Calculates the height above the lidar values per bin and channel
 -- range calculation: Calculates the range above the lidar values per bin and channel
 -- range_correction: Performs the range correction on signals
 -- smoothing: Smooths the signals (sliding average)
 -- trigger_correction: Perform the trigger correction per channel
 -- trim_vertically: Trim channels up to a maximum altitude
 -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels

"""

import copy
import numpy as np
import xarray as xr

from typing import Any, Dict
from utils.dataarray_utils import shallow_copy
from helper_functions.printouts import print_entry

def compute_gluing_region(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
        
    output_data = shallow_copy(input_data)
            
    ranges = output_data['range']

    background = output_data['background']
    
    profiles = output_data['profile']
    
    channel_info = output_data['channel_info']
        
    qa_tests = list(profiles.keys())

    for key in qa_tests:
        
        continue # Replace with your functions here
        # Signal is 3D data array with dimensions: (time, channel, bins)
        # Probably you don't need ranges. Check if you can find the gluing region in bins
        # The function should operate on profiles that have or have not the time dimension
        # --> average along the time dimension if it exists 
        # Save gluing map to channel_info and glued signals to signal databases of output data
        # Example: output_data['profile'][key] = sig_gl
        # Example: output_data['gluing_info'][key] = gluing_info
        # gluing_info can be a DataArray or a dictionary - maybe dictionary is more practical
            
    print_entry('Gluing region succesfully identified!')

    return output_data

def compute_gluing(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
        
    output_data = copy.deepcopy(input_data)
            
    ranges = output_data['range']
    
    profiles = output_data['profile']
    
    channel_info = output_data['channel_info']
        
    qa_tests = list(profiles.keys())

    if 'gluing_info' in output_data:
        for key in qa_tests:
            
            continue # Replace with your functions here
            # Signal is 3D data array with dimensions: (time, channel, bins) 
            # Calculate gluing factor within the gluing region and glue
            # Save gluing factor to gluing_info and glued signals to signal databases of output data
            # Example: output_data['profile'][key] = sig_gl
            # Example: output_data['gluing_info'][key] = gluing_info
            # gluing_info should be a DataArray with dims: (parameters, channel)
            # gluing_info can be a DataArray or a dictionary - maybe dictionary is more practical

        print_entry('Gluing region succesfully identified!')
    else:
        print_entry('Gluing could not be performed! Gluing_info data not found in input stage')


    return output_data
