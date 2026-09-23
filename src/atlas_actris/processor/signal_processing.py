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

import numpy as np
import pandas as pd
import xarray as xr

from typing import Any, Dict
from utils.printouts import print_header, print_subsection, print_entry

from processor.definitions import (
    assign_drk, 
    profile_instances, 
    background_map, 
    background_error_map,
    )

from utils.signal_utils import (
    temporal_averaging, 
    temporal_averaging_error, 
    fast_rolling_mean_range,
    fast_rolling_noise,
    fast_rolling_mean,
    _rebin_and_trim_all_binned_arrays,
    _true_bin_grid,
    )

from utils.error_classes import CustomWarning
from utils.dataarray_utils import shallow_copy

def _drop_or_mean_time(da: xr.DataArray) -> xr.DataArray:
    if "time" not in da.dims:
        return da
    if da.sizes["time"] == 1:
        return da.squeeze("time", drop=True)
    return da.mean("time", skipna=True)

def _restore_dim_order(da: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    preferred = [dim for dim in template.dims if dim in da.dims]
    extra = [dim for dim in da.dims if dim not in preferred]
    return da.transpose(*preferred, *extra)

def compute_height_and_range_calculation(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)
    
    system_info = output_data["system_info"]
    channel_info = output_data["channel_info"]
    
    profiles = output_data["profile"]
        
    qa_tests = list(channel_info.keys())
    
    for key in qa_tests:
        
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
        
        bins = profiles[key].bins
        
        zenith_angle_rad = np.pi * zenith_angle / 180.0

        corrected_bins = (bins + 0.5 + zero_bin)
        
        ranges = (
            resolution * corrected_bins
        ).reset_coords(drop=True)
        
        height_agl = (
            ranges * np.cos(zenith_angle_rad)
        ).reset_coords(drop=True)

        output_data["bins"][key] = corrected_bins
        output_data["range"][key] = ranges
        output_data["height_agl"][key] = height_agl
        
        if station_altitude is not None:
            output_data["height_asl"][key] = (
                height_agl + float(station_altitude)
            ).reset_coords(drop=True)
    
    print_entry("Ranges/heights calculated sucessfully")
    
    return output_data

def compute_unit_conv_counts_to_MHz(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)
    channel_info = output_data["channel_info"]

    for prof_key in profile_instances:
        profiles = output_data[prof_key]
        if not profiles:
            continue

        for key, sig in profiles.items():

            if key not in channel_info:
                continue

            ci = channel_info[key]

            acquisition_mode = ci.sel(parameters="acquisition_mode")
            range_resolution = ci.sel(parameters="range_resolution").astype(sig.dtype)

            is_photon = acquisition_mode == "p"

            if not bool(is_photon.any()):
                output_data[prof_key][key] = sig
                continue

            if prof_key == "profile":
                shots = output_data["shots"][key]
            else:
                shots = output_data["shots"][key].median("time")

            factor = (150.0 / range_resolution) / shots

            sig_out = xr.where(is_photon, sig * factor, sig)
            sig_out = _restore_dim_order(sig_out, sig)
            
            output_data[prof_key][key] = sig_out
            
    print_entry("Unit conversion (counts to countrate in MHz) for photon channels complete!")
    return output_data

def compute_dead_time_correction(processing_info, input_data):

    output_data = shallow_copy(input_data)
    channel_info = output_data["channel_info"]

    for prof_key in profile_instances:
        profiles = output_data[prof_key]
        background = output_data['background_mean']

        if not profiles:
            continue

        for key, sig in profiles.items():

            if key not in channel_info:
                continue
            
            ci = channel_info[key]

            acquisition_mode = ci.sel(parameters="acquisition_mode")
            dead_time = ci.sel(parameters="dead_time").astype(sig.dtype)

            is_photon = acquisition_mode == "p"

            if not bool(is_photon.any()):
                output_data[prof_key][key] = sig
                continue
            
            if key in background:
                bg = background[key]
            
                mask_saturated = (bg > 60.0) & is_photon
            
                if bool(mask_saturated.any()):
                    channel_mask = mask_saturated.any(
                        dim=[dim for dim in mask_saturated.dims if dim != "channel"]
                    ).compute()
            
                    warned_channels = bg["channel"].where(
                        channel_mask,
                        drop=True,
                    ).values
            
                    dead_time = dead_time.where(~channel_mask, 0.0)
            
                    print_entry(key)
                    print(
                        "Warning: Dead-time correction was deactivated for channels "
                        f"with background > 60 MHz: {', '.join(sorted(warned_channels))}"
                    )
                    print()

            denom = 1.0 - sig * dead_time * 1e-3
            sig_corr = (sig / denom).where(denom != 0)

            sig_out = xr.where(is_photon, sig_corr, sig)
            sig_out = _restore_dim_order(sig_out, sig)
            
            output_data[prof_key][key] = sig_out

    print_entry("Dead time correction succesfully performed!")
    return output_data

def compute_background_calculation(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)
    channel_info = output_data["channel_info"]

    for prof_key in profile_instances:
        profiles = output_data[prof_key]
        if not profiles:
            continue

        back_key = background_map[prof_key]
        back_error_key = background_error_map[prof_key]

        for key, sig in profiles.items():

            if key not in channel_info:
                continue

            ci = channel_info[key]

            background_low_bin = ci.sel(parameters="background_low_bin")
            background_high_bin = ci.sel(parameters="background_high_bin")

            mask_bins = (
                (sig["bins"] >= background_low_bin)
                & (sig["bins"] <= background_high_bin)
            )

            sig_bg = sig.where(mask_bins)

            bg_mean = sig_bg.mean("bins", skipna=True)
            bg_std = sig_bg.std("bins", skipna=True)
            n_bins = sig_bg.notnull().sum("bins")

            output_data[back_key][key] = bg_mean
            output_data[back_error_key][key] = bg_std / np.sqrt(n_bins)

    print_entry("Background calculated sucessfully")
    return output_data
                    
def compute_averaging_by_time_single(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)
            
    profile = output_data['profile']
    profile_error = output_data['profile_error']
    background = output_data['background']
    background_error = output_data['background_error']
    
    qa_tests = list(profile.keys())
    
    for key in qa_tests:
        
        sig = profile[key]            
        sig_avg = sig.mean(dim="time", keepdims = True)
        
        output_data['profile_mean'][key] = sig_avg
        
        if key in background:
            bgd = background[key]
            
            bgd_avg = bgd.mean(dim="time", keepdims = True)

            output_data['background_mean'][key] = bgd_avg
        
        if key in profile_error:
            sig_err = profile_error[key]

            N = sig_err.notnull().sum(dim="time")
            sig_avg_err = sig_err.mean(dim="time", keepdims = True) / np.sqrt(N) 

            output_data['profile_error_mean'][key] = sig_avg_err

        if key in background_error:
            bgd_err = background_error[key]

            N = bgd_err.notnull().sum(dim="time")
            bgd_avg_err = bgd_err.mean(dim="time", keepdims = True) / np.sqrt(N) 

            output_data['background_error_mean'][key] = bgd_avg_err 

    print_entry('Single mean profile per QA test produced!')

    return output_data

def compute_averaging_by_time_low_res(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]: 

    allowed_qa_tests = ['ray', 'drk']

    low_res_averaging_period = processing_info['caller_info']["low_res_averaging_period"]
    low_res_averaging_threshold = processing_info['caller_info']["low_res_averaging_threshold"]
    
    output_data = shallow_copy(input_data)

    profile = output_data['profile']
    background = output_data['background']
    
    profile_error = output_data['profile_error']
    background_error = output_data['background_error']
    
    profile_mask = output_data['profile_mask']
    background_mask = output_data['background_mask']
    
    output_data["profile_low_res"] = _initialize_from(
        output_data["profile"]
        )
    
    output_data["profile_error_low_res"] = _initialize_from(
        output_data["profile_error"]
    )
    output_data["background_low_res"] = _initialize_from(
        output_data["background"]
    )
    output_data["background_error_low_res"] = _initialize_from(
        output_data["background_error"]
    )
    
    output_data["profile_mask_low_res"] = _initialize_from(
        output_data["profile_mask"]
    )
    output_data["background_mask_low_res"] = _initialize_from(
        output_data["background_mask"]
    )
    
    qa_tests = list(profile.keys())
    
    for key in qa_tests:
        # Change the temporal resolution average from all measurements for each QA test

        if key in allowed_qa_tests:
            
            sig = profile[key]
                                
            # Averaging the measurement
            if low_res_averaging_period:       
                sig_avg, sig_avg_mask = temporal_averaging(
                    sig = sig, 
                    averaging_period = low_res_averaging_period, 
                    averaging_threshold = low_res_averaging_threshold
                    )

                # The destination arrays were initialized from the native
                # resolution. If averaging was skipped, leave them unchanged.
                if sig_avg is None:
                    continue

                output_data['profile_low_res'][key] = sig_avg
                
                if key in background:
                    bgd = background[key]

                    bgd_avg, bgd_avg_mask = temporal_averaging(
                        sig = bgd, 
                        averaging_period = low_res_averaging_period, 
                        averaging_threshold = low_res_averaging_threshold
                        )
                    
                    if bgd_avg is not None:
                        output_data['background_low_res'][key] = bgd_avg
                
                if key in profile_error:
                    sig_err = profile_error[key]

                    sig_avg_err, _ = temporal_averaging_error(
                        sig_err = sig_err, 
                        averaging_period = low_res_averaging_period, 
                        averaging_threshold = low_res_averaging_threshold
                        )

                    if sig_avg_err is not None:
                        output_data['profile_error_low_res'][key] = sig_avg_err

                if key in background_error:
                    bgd_err = background_error[key]

                    bgd_avg_err, _ = temporal_averaging_error(
                        sig_err = bgd_err, 
                        averaging_period = low_res_averaging_period, 
                        averaging_threshold = low_res_averaging_threshold
                        )
           
            
                    if bgd_avg_err is not None:
                        output_data['background_error_low_res'][key] = bgd_avg_err
                
                if key in profile_mask:    
                    output_data['profile_mask_low_res'][key] = sig_avg_mask
                                        
                if key in background_mask:    
                    output_data['background_mask_low_res'][key] = bgd_avg_mask

    print_entry('Low resolution averaging for the rayleigh measurement complete!')

    return output_data

def compute_averaging_by_time_high_res(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    
    allowed_qa_tests = ['ray', 'drk']

    high_res_averaging_period = processing_info['caller_info']["high_res_averaging_period"]
    high_res_averaging_threshold = processing_info['caller_info']["high_res_averaging_threshold"]
    
    output_data = shallow_copy(input_data)

    profile = output_data['profile']
    background = output_data['background']
    
    profile_error = output_data['profile_error']
    background_error = output_data['background_error']
    
    profile_mask = output_data['profile_mask']
    background_mask = output_data['background_mask']
    
    output_data["profile_high_res"] = _initialize_from(
        output_data["profile"]
        )
    
    output_data["profile_error_high_res"] = _initialize_from(
        output_data["profile_error"]
    )
    output_data["background_high_res"] = _initialize_from(
        output_data["background"]
    )
    output_data["background_error_high_res"] = _initialize_from(
        output_data["background_error"]
    )
    
    output_data["profile_mask_high_res"] = _initialize_from(
        output_data["profile_mask"]
    )
    output_data["background_mask_high_res"] = _initialize_from(
        output_data["background_mask"]
    )
    
    qa_tests = list(profile.keys())
    
    for key in qa_tests:
        # Change the temporal resolution average from all measurements for each QA test

        if key in allowed_qa_tests:
            
            sig = profile[key]
                                
            # Averaging the measurement
            if high_res_averaging_period:   

                sig_avg, sig_avg_mask = temporal_averaging(
                    sig = sig, 
                    averaging_period = high_res_averaging_period, 
                    averaging_threshold = high_res_averaging_threshold
                    )

                # The destination arrays were initialized from the native
                # resolution. If averaging was skipped, leave them unchanged.
                if sig_avg is None:
                    continue

                output_data['profile_high_res'][key] = sig_avg
                
                if key in background:
                    bgd = background[key]

                    if bgd:
                        bgd_avg, bgd_avg_mask = temporal_averaging(
                            sig = bgd, 
                            averaging_period = high_res_averaging_period, 
                            averaging_threshold = high_res_averaging_threshold
                            )
                        
                        if bgd_avg is not None:
                            output_data['background_high_res'][key] = bgd_avg
                
                if key in profile_error:
                    sig_err = profile_error[key]

                    if sig_err:
                        sig_avg_err, _ = temporal_averaging_error(
                            sig_err = sig_err, 
                            averaging_period = high_res_averaging_period, 
                            averaging_threshold = high_res_averaging_threshold
                            )

                        if sig_avg_err is not None:
                            output_data['profile_error_high_res'][key] = sig_avg_err
                

                if key in background_error:
                    bgd_err = background_error[key]

                    if bgd_err:
                        bgd_avg_err, _ = temporal_averaging_error(
                            sig_err = bgd_err, 
                            averaging_period = high_res_averaging_period, 
                            averaging_threshold = high_res_averaging_threshold
                            )
               
                
                        if bgd_avg_err is not None:
                            output_data['background_error_high_res'][key] = bgd_avg_err
                
                if key in profile_mask:    
                    output_data['profile_mask_high_res'][key] = sig_avg_mask
                                        
                if key in background_mask:    
                    output_data['background_mask_high_res'][key] = bgd_avg_mask

    print_entry('High resolution averaging for quicklooks complete!')

    return output_data

def compute_trim_vertically_old(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    max_height_agl = 1E3 * processing_info["caller_info"]["max_height_agl"]

    output_data = shallow_copy(input_data)

    height_agl = output_data["height_agl"]
    height_asl = output_data["height_asl"]
    ranges = output_data["range"]
    channel_info = output_data["channel_info"]

    if not height_agl:
        print_entry(
            "Vertical trimming could not be performed! "
            "Range data not found in input stage"
        )
        return output_data

    qa_tests = list(ranges.keys())

    for key in qa_tests:

        if key not in height_agl:
            continue
        if key not in height_asl:
            continue
        if key not in channel_info:
            continue

        z_agl = height_agl[key]
        z_asl = height_asl[key]
        z_rng = ranges[key]

        zero_bin = channel_info[key].sel(parameters="zero_bin")

        # Use the coordinate of the vertical arrays as the reference.
        bins = z_agl.bins

        # Direct mask. Avoid argmax because it returns a position index,
        # not necessarily the actual bins coordinate value.
        mask_bins = (bins >= zero_bin) & (z_agl <= max_height_agl)

        # Trim all profile instances using the same mask.
        for prof_key in profile_instances:

            profiles = output_data[prof_key]

            if not profiles:
                continue
            if key not in profiles:
                continue

            sig = profiles[key]

            sig_trm = sig.where(mask_bins, drop=True)

            output_data[prof_key][key] = sig_trm.reset_coords(drop=True)

        # Trim vertical coordinate stores only once per QA key,
        # not once per profile instance.
        z_agl_trm = z_agl.where(mask_bins, drop=True)
        z_asl_trm = z_asl.where(mask_bins, drop=True)
        z_rng_trm = z_rng.where(mask_bins, drop=True)

        output_data["height_agl"][key] = z_agl_trm.reset_coords(drop=True)
        output_data["height_asl"][key] = z_asl_trm.reset_coords(drop=True)
        output_data["range"][key] = z_rng_trm.reset_coords(drop=True)

    print_entry("Vertical trimming succesfully performed!")

    return output_data
      
def compute_trim_vertically(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Re-bin and vertically trim all binned arrays lazily.

    Only the small channel/bin coordinate arrays are evaluated. Large profile
    arrays remain Dask-backed and are reduced with positional ``isel`` slices,
    avoiding repeated ``where(..., drop=True)`` operations and large task
    graphs.
    """

    max_height_agl = 1.0e3 * processing_info["caller_info"]["max_height_agl"]
    output_data = shallow_copy(input_data)

    ranges = output_data["range"]
    channel_info = output_data["channel_info"]
    system_info = output_data["system_info"]

    if not ranges:
        print_entry(
            "Vertical trimming could not be performed! "
            "Range data not found in input stage"
        )
        return output_data

    for key, range_ref in ranges.items():
        if key not in channel_info or key not in system_info:
            continue

        # The existing range array is the cheapest and most reliable source of
        # the original bins length; avoid scanning every store for a reference.
        if not isinstance(range_ref, xr.DataArray) or "bins" not in range_ref.dims:
            continue

        ci = channel_info[key]
        zero_bin = ci.sel(parameters="zero_bin").astype("float32").reset_coords(drop=True)
        resolution = (
            ci.sel(parameters="range_resolution")
            .astype("float32")
            .reset_coords(drop=True)
        )

        rename_dims = {}
        if "channels" in zero_bin.dims and "channel" not in zero_bin.dims:
            rename_dims["channels"] = "channel"
        if rename_dims:
            zero_bin = zero_bin.rename(rename_dims)

        rename_dims = {}
        if "channels" in resolution.dims and "channel" not in resolution.dims:
            rename_dims["channels"] = "channel"
        if rename_dims:
            resolution = resolution.rename(rename_dims)

        sys = system_info[key]
        zenith_angle_rad = np.deg2rad(float(sys.loc["zenith_angle"].values))
        cos_zenith = np.cos(zenith_angle_rad)
        station_altitude = sys.loc["station_altitude"].values

        # Build the complete true-bin grid available in the input data. Do not
        # use the requested AGL cutoff to limit this helper call: the requested
        # limit applies to height above ground, not to slant range. The exact
        # cutoff is applied below after ``z_agl`` has been calculated.
        available_height_agl = range_ref * cos_zenith
        available_max_height_agl = float(
            available_height_agl.max(skipna=True).values
        )

        target_bins = _true_bin_grid(
            zero_bin=zero_bin,
            n_old_bins=range_ref.sizes["bins"],
            resolution=resolution,
            zenith_angle_rad=zenith_angle_rad,
            max_height_agl=available_max_height_agl,
        )

        if len(target_bins) == 0:
            continue

        bins = xr.DataArray(
            target_bins,
            dims=("bins",),
            coords={"bins": target_bins},
        )
        based_bins = bins + np.float32(0.5)
        z_rng = (resolution * based_bins).reset_coords(drop=True)
        z_agl = (z_rng * cos_zenith).reset_coords(drop=True)

        # Evaluate only this small metadata mask. Determine one shared trailing
        # slice and apply it before touching any large profile array.
        mask_bins = z_agl <= max_height_agl
        non_bin_dims = [dim for dim in mask_bins.dims if dim != "bins"]
        valid_any = mask_bins.any(dim=non_bin_dims) if non_bin_dims else mask_bins
        valid_any_np = np.asarray(valid_any.values, dtype=bool)

        valid_positions = np.flatnonzero(valid_any_np)
        if valid_positions.size == 0:
            continue

        stop = int(valid_positions[-1]) + 1
        bin_slice = slice(0, stop)

        target_bins = target_bins[bin_slice]
        mask_bins = mask_bins.isel(bins=bin_slice)
        z_rng = z_rng.isel(bins=bin_slice)
        z_agl = z_agl.isel(bins=bin_slice)

        # The helper performs the channel-dependent true-bin alignment. Giving
        # it the already-trimmed target grid keeps the resulting Dask graph and
        # every downstream chunk substantially smaller.
        _rebin_and_trim_all_binned_arrays(
            output_data=output_data,
            key=key,
            zero_bin=zero_bin,
            target_bins=target_bins,
            mask_bins=mask_bins,
        )

        # Positional trimming is cheap and lazy. ``where`` is only required for
        # channel-specific tail differences; do not use ``drop=True`` here.
        all_valid = bool(np.asarray(mask_bins.all().values))
        if not all_valid:
            z_rng = z_rng.where(mask_bins)
            z_agl = z_agl.where(mask_bins)

        output_data["range"][key] = z_rng.reset_coords(drop=True)
        output_data["height_agl"][key] = z_agl.reset_coords(drop=True)

        if station_altitude is not None:
            z_asl = z_agl + float(station_altitude)
            output_data["height_asl"][key] = z_asl.reset_coords(drop=True)

    print_entry("Vertical trimming succesfully performed!")
    return output_data

def compute_background_correction(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
        
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)
    
    qa_tests = list(output_data['channel_info'].keys())

    if output_data['background']:
        for prof_key in profile_instances:
            profiles = output_data[prof_key]
            back_key = background_map[prof_key]
            background = output_data[back_key]
            
            if profiles:
                # Create a single average from all measurements for each QA test
                for key in qa_tests:
                    sig = profiles[key]
                    bc = background[key]
                    
                    sig_bc = sig - bc
                    
                    output_data[prof_key][key] = sig_bc
            
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

    ranges = output_data['range']
    
    qa_tests = list(ranges.keys())

    if ranges:
        for prof_key in profile_instances:
            profiles = output_data[prof_key]
            if profiles:
                for key in qa_tests:
            
                    z_rng = ranges[key]#.astype("float32")
            
                    sig = profiles[key]
                    
                    sig_rc = sig * z_rng**2
                    
                    mask_below_zero_range = (z_rng >= 0.)
                    
                    sig_rc = sig_rc.where(mask_below_zero_range, sig)
                    
                    output_data[prof_key][key] = sig_rc
            
        print_entry('Range correction succesfully performed!')
    
    else:
        print_entry('Range correction could not be performed! Height (ASL) data not found in input stage')

    return output_data

def compute_smoothing_dark(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)

    ranges = output_data["range"]

    dark_profiles = output_data["profile_mean"]

    for key in dark_profiles:

        if not key.startswith("drk"):
            continue

        z_rng = ranges[key]

        drk = dark_profiles[key]

        drk_sm = fast_rolling_mean_range(
            da=drk,
            ranges=z_rng,
            window=1000,
            smooth_above=2000.0,
            dim="bins",
        )
        
        output_data["profile_mean"][key] = drk_sm

    print_entry("Dark smoothing succesfully performed!")

    return output_data

def compute_dark_correction(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    
    loading_map = processing_info['caller_info']['loading_map']
    
    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    channel_info = output_data['channel_info']
    
    dark_profiles = output_data['profile_mean']
            
    for key, drk_key_alias in assign_drk.items():
        for prof_key in profile_instances:
            profiles = output_data[prof_key]
            if profiles:

                if key not in profiles:
                    continue
        
                if drk_key_alias in dark_profiles:
                    drk_key = drk_key_alias
                
                elif drk_key_alias in loading_map:
                    candidate = loading_map[drk_key_alias]
                
                    if candidate not in dark_profiles:
                        print_entry(
                            f"Dark measurement {candidate!r} was configured but "
                            "no dark profiles were successfully read. Skipping correction."
                        )
                        continue
                
                    drk_key = candidate
                
                else:
                    continue
                
                drk = dark_profiles[drk_key]
            
                drk = dark_profiles[drk_key]                   
                
                acquisition_mode = channel_info[key].sel(parameters="acquisition_mode")
                                
                mask_a = acquisition_mode == "a"

                # Apply only to analog channels.
                drk = xr.where(mask_a, drk, 0.0)
            
                sig_drc = profiles[key] - drk.squeeze("time", drop=True)
                
                output_data[prof_key][key] = sig_drc

    print_entry('Dark correction succesfully performed!')

    return output_data

def compute_signal_noise(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)

    profiles = output_data["profile"]
    profiles_m = output_data["profile_mean"]

    caller_info = processing_info["caller_info"]

    noise_smooth_window = caller_info.get("noise_smooth_window", 7)
    noise_window = caller_info.get("noise_window", 31)

    # Chunk size used only for the synthetic high-resolution error array.
    # Keeping time chunks small avoids large materialized blocks later.
    noise_time_chunk = caller_info.get("noise_time_chunk", 1)

    for key, sig in profiles.items():

        if key not in profiles_m:
            continue

        sig_m = profiles_m[key]

        sig_m_err = fast_rolling_noise(
            sig_m,
            smooth_window=noise_smooth_window,
            noise_window=noise_window,
        )

        sig_m_err = _drop_or_mean_time(sig_m_err)

        # Mean-profile noise. This stays small.
        output_data["profile_error_mean"][key] = sig_m_err.broadcast_like(sig_m)

        # Number of high-resolution profiles used in the mean.
        # Metadata only.
        N = sig.sizes["time"]

        # Convert the small mean-profile error to a Dask-backed array.
        # This prevents xarray from eagerly broadcasting it over time.
        sig_m_err_dask = (sig_m_err * np.sqrt(N)).chunk({
            dim: sig.chunksizes.get(dim, sig_m_err.sizes[dim])
            for dim in sig_m_err.dims
            if dim in sig_m_err.sizes
        })

        # Create a Dask-backed 1D time template.
        # Important: chunk this before multiplying.
        time_template = xr.ones_like(
            sig["time"],
            dtype=sig_m_err.dtype,
        ).chunk({
            "time": noise_time_chunk,
        })

        # This broadcast is now lazy because at least one operand is Dask-backed.
        sig_err = time_template * sig_m_err_dask

        sig_err = _restore_dim_order(sig_err, sig)

        output_data["profile_error"][key] = sig_err

    print_entry("Signal error calculation succesfully performed!")
    return output_data

def compute_signal_smoothing(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    output_data = shallow_copy(input_data)

    profiles_m = output_data["profile_mean"]
    profiles_m_error = output_data["profile_error_mean"]

    ranges = output_data["range"]

    caller_info = processing_info["caller_info"]

    smoothing_window_bins = caller_info.get("smoothing_window", 133)
    
    for key, sig in profiles_m.items():

        if key not in profiles_m:
            continue
        
        if key not in profiles_m_error:
            continue

        sig_m = profiles_m[key]
        sig_m_err = profiles_m_error[key]
        
        sig_m_sm = fast_rolling_mean(
            sig_m,
            window = smoothing_window_bins,
        ) 
        
        sig_m_sm_err = sig_m_err / np.sqrt(smoothing_window_bins)

        # Mean-profile noise. This stays small.
        output_data["profile_mean"][key] = sig_m_sm.broadcast_like(sig_m)

        output_data["profile_error_mean"][key] = sig_m_sm_err.broadcast_like(sig_m_err)

    print_entry("Signal error calculation succesfully performed!")
    return output_data

def compute_mean_arrays(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:

    # output_data = copy.deepcopy(input_data)
    output_data = shallow_copy(input_data)

    profiles = output_data['profile_mean']
        
    profile_error = output_data['profile_error_mean']
    
    background = output_data['background_mean']

    background_error = output_data['background_error_mean']
    
    qa_tests = list(profiles.keys())

    for key in qa_tests:

        if key in profiles:
            output_data['profile_mean'][key] = profiles[key].persist()

        if key in profile_error:
            output_data['profile_error_mean'][key] = profile_error[key].persist()
        
        if key in background:
            output_data['background_mean'][key] = background[key].persist()

        if key in background_error:
            output_data['background_error_mean'][key] = background_error[key].persist()

    print_entry('Mean arrays succesfully computed!')
    
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

    range_slice : xarray.DataArray
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
    range_slice = ranges.where(mask, drop=drop)

    return sig_slice, range_slice
    

def _initialize_from(source):
    return source.copy() if isinstance(source, dict) else source