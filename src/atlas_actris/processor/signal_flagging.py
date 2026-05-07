#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:02:56 2026

@author: nikos
"""
import copy

from typing import Any, Dict
from utils.error_classes import CustomWarning
from helper_functions.printouts import print_subsection

def compute_detect_saturation(
    caller_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    """
    General:
        Detects lidar profiles with photon values close to the
        maximum allowed countrate or analog values close to the data 
        acquisition range
        
    Input:
        sig: 
            A 2D or 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, ...). 

        shots: 
            A 2D xarray with the laser shots per channel and timeframe.
            It should include the following dimensions: (time, channel, ...) 
            The index should correspond to the channel dimension of sig            

        metadata: A pandas Dataframe with at least the following columns and the 
        atlas_channel_id as index:
            
            resol: 
                The raw range resolution per channel in meters
                
            dead_time: 
                The dead time per channel in nanoseconds. 
    
            DAQ_range: 
                The data acquisition range of the analog
                channels.
            
    """
    
    
    output_data = copy.deepcopy(input_data)
    # output_data.setdefault("time_mask", {})

    profiles = output_data["profile"]
    shots = output_data["shots"]
    channel_info = output_data["channel_info"]

    qa_tests = list(profiles.keys())

    for key in qa_tests:
            
        print_subsection(f"{key} dataset")
                
        channels = profiles[key].channel.values
    
        aq_mode = channel_info[key].sel({"parameters": "acquisition_mode"}).drop_vars("parameters")

        daq_range = channel_info[key].sel({"parameters": "data_acquisition_range"}).drop_vars("parameters")
        
        sig_clean = profiles[key].fillna(0.)
        daq_clean = daq_range.fillna(0.)
            
        max_counts = 1.
        
        max_daq_range = 0.80 * daq_clean   
        
        mask_analog_saturated = (sig_clean > max_daq_range) & (aq_mode == "a")
        mask_photon_saturated = (sig_clean / shots[key] > max_counts) & (aq_mode == "p")
        
        mask_saturated = (mask_analog_saturated | mask_photon_saturated).compute()
        
        mask_saturated_ch = mask_saturated.any("bins").any("time")
        
        for ch in channels:
                    
            if mask_saturated_ch.sel(channel=ch) and \
                aq_mode.sel(channel=ch) == "a":
                print("")
                CustomWarning(f"Channel {ch} - Analog signal mV values above 80% of the data acqusition range were detected! ")   
    
            if mask_saturated_ch.sel(channel=ch) and \
                aq_mode.sel(channel=ch) == "p":
                print("")
                CustomWarning(f"Channel {ch} - Photon signal values above 1 count per bin per shot were detected! ")

        print("")

        output_data["profile_mask"][key] = mask_saturated
        
    return output_data

