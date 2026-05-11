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
from helper_functions.printouts import print_header, print_subsection, print_entry

from helper_functions.signal_utils import temporal_averaging

from utils.error_classes import CustomWarning
from utils.dataarray_utils import shallow_copy


def compute_height_and_range_calculation(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)
    
    system_info = output_data["system_info"]
    channel_info = output_data["channel_info"]
        
    qa_tests = list(channel_info.keys())
    
    for key in qa_tests:
        
        max_bins = int(channel_info[key].loc["bins", :].max().values)

        resolution = (
            channel_info[key]
            .loc["range_resolution", :]
            .astype("float32")
            .reset_coords(drop=True)
        )

        zero_bin = (
            channel_info[key]
            .loc["zero_bin", :]
            .astype("float32")
            .reset_coords(drop=True)
        )

        zenith_angle = float(system_info[key].loc["zenith_angle"].values)
        station_altitude = system_info[key].loc["station_altitude"].values
        
        bin_arr = np.arange(max_bins, dtype=np.int32)

        bins = xr.DataArray(
            bin_arr.astype(np.float32),
            dims=["bins"],
            coords={"bins": bin_arr},
            name="bins",
        )
        
        zenith_angle_rad = np.pi * zenith_angle / 180.0
        
        ranges = (
            resolution * (bins + 0.5 + zero_bin)
        ).reset_coords(drop=True)

        heights_agl = (
            ranges * np.cos(zenith_angle_rad)
        ).reset_coords(drop=True)

        output_data["range"][key] = ranges
        output_data["height_agl"][key] = heights_agl
        
        if station_altitude is not None:
            output_data["height_asl"][key] = (
                heights_agl + float(station_altitude)
            ).reset_coords(drop=True)
    
    print_entry("Ranges/heights calculated sucessfully")
    
    return output_data

def compute_unit_conv_counts_to_MHz(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    channel_info = output_data["channel_info"]
    shots = output_data["shots"]
    profiles = output_data["profile"]

    for key in channel_info:
        ci = channel_info[key]
        sig = profiles[key]

        acquisition_mode = ci.sel(parameters="acquisition_mode")
        range_resolution = ci.sel(parameters="range_resolution").astype(sig.dtype)

        photon_channels = acquisition_mode.channel.values[
            acquisition_mode.values == "p"
        ]

        if len(photon_channels) == 0:
            output_data["profile"][key] = sig
            continue

        factor = (150.0 / range_resolution.sel(channel=photon_channels)) / shots[key]

        sig_out = sig.copy()
        sig_out.loc[dict(channel=photon_channels)] = (
            sig.sel(channel=photon_channels) * factor
        )

        output_data["profile"][key] = sig_out

    print_entry("Unit conversion (counts to countrate in MHz) for photon channels complete!")
    return output_data

def compute_dead_time_correction(processing_info, input_data):
    
    output_data = shallow_copy(input_data)

    channel_info = output_data["channel_info"]
    profiles = output_data["profile"]

    for key in channel_info:
        ci = channel_info[key]
        sig = profiles[key]

        acquisition_mode = ci.sel(parameters="acquisition_mode")
        dead_time = ci.sel(parameters="dead_time").astype(sig.dtype)

        photon_channels = acquisition_mode.channel.values[
            acquisition_mode.values == "p"
        ]

        if len(photon_channels) == 0:
            output_data["profile"][key] = sig
            continue

        sig_p = sig.sel(channel=photon_channels)
        dt_p = dead_time.sel(channel=photon_channels)

        denom = 1.0 - sig_p * dt_p * 1e-3
        sig_p_corr = sig_p / denom
        sig_p_corr = sig_p_corr.where(denom != 0)

        sig_out = sig.copy()
        sig_out.loc[dict(channel=photon_channels)] = sig_p_corr

        output_data["profile"][key] = sig_out

    print_entry("Dead time correction succesfully performed!")
    
    return output_data

def compute_background_calculation(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    channel_info = output_data['channel_info']
        
    profiles = output_data['profile']

    qa_tests = list(channel_info.keys())
    
    for key in qa_tests:
        # Create a single average from all measurements for each QA test
        
        background_low_bin = channel_info[key]\
            .sel(parameters = 'background_low_bin')
        
        background_high_bin = channel_info[key]\
            .sel(parameters = 'background_high_bin')
            
        sig = profiles[key]
        
        mask_bins = (sig["bins"] >= background_low_bin) & \
            (sig["bins"] <= background_high_bin)
        
        sig_masked = sig.where(mask_bins)
        
        bg_mean_bins = sig_masked.mean("bins", skipna = True)
        bg_sdev_bins = sig_masked.std("bins", skipna = True)
                    
        output_data['background'][key] = bg_mean_bins
        output_data['background_error'][key] = bg_sdev_bins

    print_entry('Background calculated sucessfully')
    
    return output_data

                    
def compute_averaging_by_time_single(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)
            
    profiles = output_data['profile']
    
    qa_tests = list(profiles.keys())
    
    for key in qa_tests:
        
        sig = profiles[key]
            
        sig_avg = sig.mean(dim="time", keepdims = True)
        
        output_data['profile'][key] = sig_avg

    print_entry('Single mean profile per QA test produced!')

    return output_data

    
def compute_averaging_by_time_low_res(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]: 

    allowed_qa_tests = ['ray', 'drk']

    ray_averaging_rate = processing_info['caller_info']["ray_averaging_rate"]
    ray_averaging_threshold = processing_info['caller_info']["ray_averaging_threshold"]
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    profiles = output_data['profile']
    
    qa_tests = list(profiles.keys())
    
    for key in qa_tests:
        # Change the temporal resolution average from all measurements for each QA test

        if key in allowed_qa_tests:
            
            sig = profiles[key]
                
            # Averaging the measurement
            if ray_averaging_rate == None:
                sig_avg = sig.mean(dim="time", keepdims = True)
                mask_incomplete_avg = xr.full_like(
                    sig_avg, 
                    fill_value = False, 
                    dtype=bool
                    )
            else:           
                sig_avg, mask_incomplete_avg = temporal_averaging(
                    sig = sig, 
                    averaging_rate = ray_averaging_rate, 
                    averaging_threshold = ray_averaging_threshold
                    )
               
            output_data['profile'][key] = sig_avg
            output_data['profile_mask'][key] = mask_incomplete_avg

    print_entry('Low resolution averaging for the rayleigh measurement complete!')

    return output_data

def compute_averaging_by_time_high_res(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    
    allowed_qa_tests = ['ray', 'drk']

    ray_qck_averaging_rate = processing_info['caller_info']["ray_qck_averaging_rate"]
    ray_qck_averaging_threshold = processing_info['caller_info']["ray_qck_averaging_threshold"]
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    profiles = output_data['profile']
    
    qa_tests = list(profiles.keys())
    
    for key in qa_tests:
        # Change the temporal resolution average from all measurements for each QA test

        if key in allowed_qa_tests:
            
            sig = profiles[key]
                
            # Averaging the measurement
            if ray_qck_averaging_rate == None:
                sig_avg = sig.mean(dim="time", keepdims = True)
                mask_incomplete_avg = xr.full_like(
                    sig_avg, 
                    fill_value = False, 
                    dtype=bool
                    )
            else:           
                sig_avg, mask_incomplete_avg = temporal_averaging(
                    sig = sig, 
                    averaging_rate = ray_qck_averaging_rate, 
                    averaging_threshold = ray_qck_averaging_threshold
                    )
               
            output_data['profile'][key] = sig_avg
            output_data['profile_mask'][key] = mask_incomplete_avg

    print_entry('High resolution averaging for quicklooks complete!')

    return output_data

def compute_trim_vertically(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    max_height_agl = 1E3 * processing_info['caller_info']['max_height_agl']
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    heights_agl = output_data['height_agl']
    heights_asl = output_data['height_asl']
    ranges = output_data['range']

    profiles = output_data['profile']
    
    channel_info = output_data['channel_info']
    
    qa_tests = list(profiles.keys())
    
    if heights_agl:

        for key in qa_tests:
         
            sig = profiles[key]

            bins = sig.bins
            
            z_agl = heights_agl[key]
            z_asl = heights_asl[key]
            z_rng = ranges[key]
            
            zero_bin = channel_info[key]\
                .sel(parameters = 'zero_bin')
                        
            max_bin = z_agl.where(z_agl <= max_height_agl).argmax('bins')
            
            mask_bins = (bins >= zero_bin) & (bins <= max_bin)

            sig_trm = sig.where(mask_bins, drop = True)
            
            z_agl_trm = z_agl.where(mask_bins, drop = True)
            z_asl_trm = z_asl.where(mask_bins, drop = True)
            z_rng_trm = z_rng.where(mask_bins, drop = True)
            
            output_data['profile'][key] = sig_trm.reset_coords(drop=True)
            
            output_data['height_agl'][key] = z_agl_trm.reset_coords(drop=True)
            output_data['height_asl'][key] = z_asl_trm.reset_coords(drop=True)
            output_data['range'][key] = z_rng_trm.reset_coords(drop=True)
    
        print_entry('Vertical trimming succesfully performed!')
        
    else:
        print_entry('Vertical trimming could not be performed! Range data not found in input stage')

    return output_data
      
def compute_background_correction(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
        
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    profiles = output_data['profile']

    background = output_data['background']
    
    qa_tests = list(profiles.keys())

    if background:

        # Create a single average from all measurements for each QA test
        for key in qa_tests:
            sig = profiles[key]
            bc = background[key]
            
            sig_bc = sig - bc
            
            output_data['profile'][key] = sig_bc
            
        print_entry('Background correction succesfully performed!')

    else:
        print_entry('Background correction could not be performed! Background data not found in input stage')
                
    return output_data

def compute_range_correction(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    profiles = output_data['profile']

    ranges = output_data['range']
    
    qa_tests = list(profiles.keys())

    if ranges:

        for key in qa_tests:
    
            z_rng = ranges[key]#.astype("float32")
    
            sig = profiles[key]
            
            sig_rc = sig * z_rng**2
            
            mask_below_zero_range = (z_rng >= 0.)
            
            sig_rc = sig_rc.where(mask_below_zero_range, sig)
            
            output_data['profile'][key] = sig_rc
    
        print_entry('Range correction succesfully performed!')
    
    else:
        print_entry('Range correction could not be performed! Height (ASL) data not found in input stage')

    return output_data

def compute_dark_correction(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    assign_drk = {
        'ray':'drk_ray',
        'ray_pcb':'drk_ray_pcb',
        'pcb_p45':'drk_pcb',
        'pcb_m45':'drk_pcb',
        'pcb_aux_p45':'drk_pcb_aux',
        'pcb_aux_m45':'drk_pcb_aux',
        'tlc_north':'drk_tlc',
        'tlc_east':'drk_tlc',
        'tlc_south':'drk_tlc',
        'tlc_west':'drk_tlc',
        'tlc_inner':'drk_tlc',
        'tlc_outer':'drk_tlc',
        'trg':'drk_trg',
        'dtm_fi':'drk_dtm',
        'dtm_fo':'drk_dtm',
        'dtm':'drk_dtm',
        'nsf':'drk_nsf',
        'cam':'drk_cam',
        }
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    loading_map = processing_info["loading_map"]

    profiles = output_data['profile']

    channel_info = output_data['channel_info']
            
    for key in assign_drk:
        
        drk_key_alias = assign_drk[key]
        
        if drk_key_alias in loading_map and key in profiles:
        
            drk_key = loading_map[drk_key_alias]
        
            drk = profiles[drk_key].mean(dim='time', skipna = True)
            
            acquisition_mode = channel_info[key]\
                .sel(parameters = 'acquisition_mode')
                                
            mask_p = (acquisition_mode == "a")
                        
            drk = xr.where(mask_p, drk, 0.)
            
            sig_drc = profiles[key] - drk
            
            output_data['profile'][key] = sig_drc

    print_entry('Dark correction succesfully performed!')

    return output_data


def slice_along_bins(arr, indices, min_ind, max_ind, drop=False, min_bins=10):
    """
    Slice an array within an interval along the bins dimension.

    Parameters
    ----------
    arr : xarray.DataArray (time, signal, bins)
        3D data where the first dimension is time and the second
        dimension corresponds to range bins.

    indices : xarray.DataArray (bins,)
        1D array containing the coordinate associated with the bins.

    min_ind : scalar or xarray.DataArray (channel,)
        Lower limit used to select the desired interval.
        Values must have the same units as indices.
    
    max_ind : scalar or xarray.DataArray (channel,)
        Upper limit used to select the desired interval.
        Values must have the same units as indices.

    drop : bool, optional (default=False)
        Controls how values outside the selected range are handled:
        - False: keep original dimensions and replace values outside
          the region with NaN.
        - True: remove bins outside the region.

    min_bins : int, optional (default=10)
        Minimum number of bins required when `drop=True`. If fewer bins
        satisfy the condition, the original arrays are returned and a
        warning is issued.

    Returns
    -------
    arr_slice : xarray.DataArray
        Signal array restricted to the specified range region.

    indices_slice : xarray.DataArray
        Corresponding range coordinate after applying the same selection.
    """

    # create boolean mask    
    mask = (indices >= min_ind) & (indices <= max_ind)
    
    # apply mask
    arr_slice = arr.where(mask, drop=drop)
    indices_slice = indices.where(mask, drop=drop)

    return arr_slice, indices_slice

def slice_signal_by_range(sig, ranges, region, drop=False, min_bins=10):
    """
    Slice a 2D signal array by a range interval.

    Parameters
    ----------
    sig : xarray.DataArray (time, signal, bins)
        3D signal data where the first dimension is time and the second
        dimension corresponds to range bins.

    ranges : xarray.DataArray (bins,)
        1D array containing the range coordinate associated with the bins.

    region : list or tuple [min_range, max_range]
        Range limits used to select the desired interval.
        Values must have the same units as `ranges`.

    drop : bool, optional (default=False)
        Controls how values outside the selected range are handled:
        - False: keep original dimensions and replace values outside
          the region with NaN.
        - True: remove bins outside the region.

    min_bins : int, optional (default=10)
        Minimum number of bins required when `drop=True`. If fewer bins
        satisfy the condition, the original arrays are returned and a
        warning is issued.

    Returns
    -------
    sig_slice : xarray.DataArray
        Signal array restricted to the specified range region.

    ranges_slice : xarray.DataArray
        Corresponding range coordinate after applying the same selection.
    """

    # unpack region limits
    min_r, max_r = region

    # create boolean mask
    mask = (ranges >= min_r) & (ranges <= max_r)

    # check mask size if dropping bins
    if drop:
        n_selected = int(mask.sum())

        if n_selected < min_bins:
            CustomWarning(
                f"Selected region [{min_r}, {max_r}] contains only "
                f"{n_selected} bins (< {min_bins}). Returning original arrays."
            )
            return sig, ranges

    # apply mask
    sig_slice = sig.where(mask, drop=drop)
    ranges_slice = ranges.where(mask, drop=drop)

    return sig_slice, ranges_slice
    
