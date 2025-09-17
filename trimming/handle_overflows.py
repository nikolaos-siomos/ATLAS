#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 21:51:02 2025

@author: nikos
"""

import numpy as np
from helper_functions.printouts import print_header
from utils.error_classes import DataOverflowError
from helper_functions.printouts import print_subsection
from typing import Any, Dict, Tuple
import sys

def get_overflow_mask(sig, acquisition_mode, daq_range, max_count = 2.**15):
    
    sig_clean = sig.fillna(0.)
    daq_clean = daq_range.fillna(0.)
    
    mask_p = (acquisition_mode == "p") & (sig_clean >= max_count)

    # daq_range = daq_range.where(daq_range.isnull(), max_count) # avoids nan comparison warnings that occur for photon channels
    mask_a_h = (acquisition_mode == "a") & (sig_clean >= daq_clean)

    mask_a_l = (acquisition_mode == "a") & (sig_clean < 0.)

    mask = (mask_a_h) | (mask_a_l) | (mask_p)
    
    reduce_dims = set(sig.dims) - {"channel"}
    
    a_h_ch = mask_a_h.any(dim=reduce_dims)
    a_l_ch = mask_a_l.any(dim=reduce_dims)
    p_ch = mask_p.any(dim=reduce_dims)
    
    p_ch   = p_ch.compute()
    a_h_ch = a_h_ch.compute()
    a_l_ch = a_l_ch.compute()
    
    mask = mask.compute()

    ch_a_h = a_h_ch["channel"].values[a_h_ch.values]
    ch_a_l = a_l_ch["channel"].values[a_l_ch.values]
    ch_p = p_ch["channel"].values[p_ch.values]
    
    # ch_a_h = mask_a_h_ch.channel.values[mask_a_h_ch.values]
    # ch_a_l = mask_p_ch.channel.values[mask_p_ch.values]
    # ch_p   = mask_p_ch.channel.values[mask_p_ch.values]

    if len(ch_a_h) > 0:
        print("")
        print("-- Warning: Analog signal mV values above the data acqusition range detected:")
        
        for ch in ch_a_h:
            print(f"    {ch} ")

    if len(ch_a_l) > 0:
        print("")
        print("-- Warning: Negative analog signal values detected:")
        
        for ch in ch_a_l:
            print(f"    {ch} ")

    if len(ch_p) > 0:
        print("")
        print("-- Warning: Photon signal count values above the maximum allowed summed counts detected:")
        
        for ch in ch_p:
            print(f"    {ch} ")
    
    return(mask)

def overflow_method_0(mask, filename):

    print("")
    print("-- At least one bin with an overflow was detected ")
    print("-- Please revise the following bins: ")

    mask = mask.compute()

    bins = mask.bins.values
    time = mask.time.values
    channel = mask.channel.values
    
    mask_t = mask.any(dim = 'bins').any(dim = 'channel').compute()

    mask_ch = mask.any(dim = 'bins').compute()
        
    time_ovf = time[mask_t.values]
    
    for t in time_ovf:
                        
        ch_ovf = channel[mask_ch.sel({"time": t}).values]
        
        for ch in ch_ovf:
            
            bins_ovf = bins[mask.sel({"time": t, "channel": ch}).values]
            print(f"    file: {filename.loc[t].values} | ch: {ch} | bins: {bins_ovf}")
    
    DataOverflowError("trim_overflows = 0 -> Overflows detected! In order to continue with an automated overflow removal use the trim_overflow argument with value 1 or 2 (default is 0) ")
    
    return()

def overflow_method_1(sig, shots, time_info, mask, filename):

    time = mask.time.compute().values
    
    mask_t = mask.any(dim = 'bins').any(dim = 'channel')
                
    time_ovf = time[mask_t]

    time_cor = time[~mask_t]
    
    print("")
    print(f"-- Warning: trim_overflows = 1 -> Removing {time_ovf.size} profiles with at least one bin with overflows: ")
            
    for t in time_ovf:        
        print(f"    {filename.loc[t].values} ")
    
    sig = sig.sel({"time": time_cor}) 

    shots = shots.sel({"time": time_cor}) 
    
    time_info = time_info.sel({"time": time_cor}) 
    
    return(sig, shots, time_info)

def overflow_method_2(sig, shots, time_info, mask, filename, max_adjacent_overflows):

    ovfs = mask.sum(dim = 'bins')

    time = mask.time.compute().values
    
    mask_t = mask.any(dim = 'bins').any(dim = 'channel')
    
    summed = mask.rolling({"bins": max_adjacent_overflows}).sum()
    
    mask_rol = (summed >= max_adjacent_overflows).any(dim = 'bins').any(dim = 'channel')
                
    time_ovf = time[mask_t]
    
    if mask_rol.any():
                                    
        time_ovf = time[mask_rol]

        time_cor = time[~mask_rol]
        
        print("")
        print(f"-- Warning: More that {max_adjacent_overflows} adjucent bins with overflows encountered in a single profile. Interpolation is too risky, these files will be removed:")
                
        for t in time_ovf:        
            print(f"    {filename.loc[t].values} ")
        
        sig = sig.sel({"time": time_cor}) 

        shots = shots.sel({"time": time_cor}) 
        
        time_info = time_info.sel({"time": time_cor}) 
    
    if (ovfs > 100).any():
        print("")
        DataOverflowError("More that 100 overflowed bins encountered in single profiles. Interpolation is too risky, please revise the input files or consider setting trim_overflows = 1")
    
    print("")
    print(f"-- Warning: trim_overflows = 2 -> Replacing overflows in {time_ovf.size} profiles: ")
            
    sig = sig.where(~mask).interpolate_na(dim = "bins", method = "linear")

    print(f"{np.sum(mask).values} overflows have been replaced by interpolating across the bins\n")
    
    return(sig)

def check_for_overflows(caller_info: Dict[str, Any], profiles: Dict[str, Any], 
                        metadata: Dict[str, Dict[str, Any]], 
                        profile_masks: Dict[str, Dict[str, Any]]) \
    -> Tuple[Dict[str, Any], Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]:

    """
    General:
        Detects and removes lidar profiles with photon values above the
        maximum allowed countrate or analog values above the data 
        acquisition range
        
    Input:
        sig: 
            A 3D xarray with the lidar signals, it should include the 
            following dimensions: (time, channel, bins). 

        shots:
            A 2D xarray with the laser shots per timeframe and channel, 
            it should include the following dimensions: (time, channel).  
            
        channel_info:
            A pandas Dataframe with exactly one index entry per lidar channel.
            the following columns must be included:
                
            acquisition_mode: 
                The acquision_mode values per channel 
                (0 for analog, 1 for photon)
                
            data_acquisition_range: 
                The data acquisition range of the analog channels.

        time_info: 
            A pandas dataframe with exactly one index entry per 
            lidar timeframe. The index should correspond to the time 
            dimension of sig. The following column must be included         
            
            filename: 
                The raw file filename. 
        
        meas_type: 
            A 3 letter identifier that specifies the measurement type,
            it can be one of ray, tlc, pcb, drk        
        
        method:
            An integer. If set to 0 only the check for overflows will be 
            performed and an error will be raised if a single value is 
            encountered. If set to 1 all files with overflows will be
            discarded. If set to 2 overflows will be removed and 
            interpolated instead
            
            
    Returns:
        
        sig, shots, time:
            If no overflows are detected then the 3 variables reuturn intact
            
            If overflows are detected:
                method = 0 : the code exits with a diagnostic error
            
                method = 1 : timeframes with at least one overflowed bin are
                             removed
                
                method = 2 : only sig changes, overflowed values are replaced
                             by interpolated values from the surrounding bins

                method = 3 : do nothing about it, use only while debugging            
    """
        
    print_header('Handling overflow values')

    method = caller_info["trim_overflows"]
    
    max_adjacent_overflows = caller_info["max_adjacent_overflows"]

    debug_signals = caller_info["debug_signals"]

    empty_keys = []
    
    for key in profiles.keys():
            
        print_subsection(f"{key} dataset")
                    
        sig = profiles[key]
        shots = metadata[key]["shots"]
        channel_info = metadata[key]["channel_info"]
        time_info = metadata[key]["time_info"]
    
        aq_mode = channel_info.sel({"parameters": "acquisition_mode"})
        
        daq_range = channel_info.sel({"parameters": "data_acquisition_range"})
        
        filename = time_info.sel({"parameters": "filename"}).compute()
            
        # Get an overflow mask for each bin
        mask_ovf = get_overflow_mask(sig, 
                                     acquisition_mode = aq_mode, 
                                     daq_range = daq_range)
    
            
        if method == 0 and mask_ovf.any(): # Detect the problematic profiles and raise error
    
            overflow_method_0(mask = mask_ovf, filename = filename)
                
        elif method == 1 and mask_ovf.any(): # Remove the problematic profiles
        
            sig, shots, time_info = overflow_method_1(sig = sig, 
                                                      shots = shots,
                                                      time_info = time_info, 
                                                      mask = mask_ovf, 
                                                      filename = filename)
            
        elif method == 2 and mask_ovf.any(): # Replace the overflowed values with interpolated ones from the nearby bins
            
            sig = overflow_method_2(sig = sig, 
                                    shots = shots, 
                                    time_info = time_info,
                                    mask = mask_ovf, 
                                    filename = filename,
                                    max_adjacent_overflows = max_adjacent_overflows)
        
        elif method == 3 and mask_ovf.any():
    
            print("-- Warning: Overflows were detected but no action has been performed! Set trim_overflows = 3 only when debugging!\n")
    
        if any(length > 0 for length in sig.sizes.values()):
            profiles[key] = sig
            metadata[key]["shots"] = shots
            metadata[key]["time_info"]
            if debug_signals:
                profile_masks[key]["overflows"] = mask_ovf
        else:
            empty_keys.append(key)
            
    for key in empty_keys:
        del profiles[key]
        del metadata[key]["shots"]
        del profile_masks[key]
        
    if profiles == {}:
        sys.exit("Endpoint 4: Handing overflows removed all measurements. No signals to process. ATLAS terminates here")

    return(profiles, metadata, profile_masks)