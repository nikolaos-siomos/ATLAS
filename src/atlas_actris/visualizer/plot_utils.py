#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:26:22 2022

@author: nick
"""

import os
import numpy as np
import pandas as pd
from PIL import Image
from xarray import DataArray
from PIL import PngImagePlugin
from utils.toolbox import round_it
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from visualizer.smoothing import sliding_average_2D
from visualizer.smoothing import sliding_average_2D_fast
from visualizer.smoothing import sliding_average_1D_fast as smooth_1D


def add_plot_metadata(plot_path, plot_metadata, plot_metadata_extra=None):
    """
    Add metadata to a PNG file.

    Notes
    -----
    PNG metadata values must be text-like. Therefore:
    - None values are skipped
    - keys are converted to strings
    - values are converted to strings
    """

    im = Image.open(plot_path)
    meta = PngImagePlugin.PngInfo()

    def add_metadata_dict(metadata_dict):
        for key, value in metadata_dict.items():

            # Skip None values
            if value is None:
                continue

            if isinstance(value, list) and len(value) == 0:
                value = ""
            elif isinstance(value, tuple) and len(value) == 0:
                value = ""
            elif isinstance(value, dict) and len(value) == 0:
                value = ""
            elif isinstance(value, np.ndarray) and value.size == 0:
                value = ""

            meta.add_text(str(key), str(value))

    add_metadata_dict(plot_metadata)

    if plot_metadata_extra is not None:
        add_metadata_dict(plot_metadata_extra)

    im.save(plot_path, "png", pnginfo=meta)

    return None
    
def perform_color_reduction(color_reduction, plot_path):
    
    if color_reduction == True:
        im = Image.open(plot_path)
        im = im.convert('P', palette = Image.ADAPTIVE, colors = 255) 
        im.save(plot_path)

    return()

def export_plot(fig, args):
    
    dpi_val = args["dpi"]
    dir_out = os.path.join(args["plot_folder"])

    os.makedirs(dir_out, exist_ok=True)

    fpath = os.path.join(dir_out, f"{args['filename']}.png")
            
    fig.savefig(fpath, dpi=dpi_val)
    
    fig.clf()
    plt.close(fig)
    
    return fpath

def clean_plots(plot_dir, pattern, exclude_pattern=None):
    
    if not os.path.isdir(plot_dir):
        return

    if exclude_pattern is None:
        exclude_patterns = []
    elif isinstance(exclude_pattern, str):
        exclude_patterns = [exclude_pattern]
    else:
        exclude_patterns = list(exclude_pattern)

    for filename in os.listdir(plot_dir):
        if pattern not in filename:
            continue

        if any(excl in filename for excl in exclude_patterns):
            continue

        file_path = os.path.join(plot_dir, filename)

        if os.path.isfile(file_path):
            os.remove(file_path)
            

def prepare_folder(caller_info, pattern, exclude_pattern=None):
    
    plot_dir = os.path.join(caller_info["output_folder"], "plots")
            
    os.makedirs(plot_dir, exist_ok=True)
            
    clean_plots(
        plot_dir=plot_dir,
        pattern=pattern,
        exclude_pattern=exclude_pattern,
    )

def collect_dict(data_list, data_keys, add_dicts = []):
    
    info = {}
    
    for i in range(len(data_keys)):
        info[data_keys[i]] = data_list[i]
    
    for d in add_dicts:
        info = d | info 
    
    return dict(sorted(info.items()))

def pass_to_args(args, data_list, data_keys):
    
    for i in range(len(data_keys)):
        args[data_keys[i]] = data_list[i]
        
    return(args)
    
def smoothing(args, x_vals, y_vals, err_type = "std"):
    
    if args['smooth'] and args['smoothing_window']:

        y_vals_sm, y_errs = \
            smooth_1D(y_vals = y_vals, 
                      x_vals = x_vals,
                      x_sm_lims = args['smoothing_range'],
                      x_sm_win = 1E3 * args['smoothing_window'],
                      expo = False,
                      err_type = err_type)
    
    else:
        y_vals_sm = y_vals.copy()
        y_errs = np.nan * y_vals.copy()
        
    return(y_vals_sm, y_errs)

def smoothing_2D(args, x_vals, y_vals, err_type = "std"):
    
    if args['smooth'] and args['smoothing_window']:
        
        if isinstance(args['smoothing_window'],list):
            smooth_2D = sliding_average_2D
        else:
            smooth_2D = sliding_average_2D_fast

        y_vals_sm, y_errs = smooth_2D(
            z_vals = y_vals, 
            y_vals = x_vals,
            y_sm_lims = args['smoothing_range'],
            y_sm_win = args['smoothing_window'],
            expo = args['smooth_exponential']
            )
    else:
        y_vals_sm = y_vals
        y_errs = np.nan * np.zeros_like(y_vals)
        
    return(y_vals_sm, y_errs)

def slice_by_vertical_scale(
    da,
    vertical_scale,
    x_lims,
    time_dim="time",
    bin_dim="bins",
):
    """
    Slice a DataArray along bins using vertical_scale values.

    Works for:
    - 1D arrays with dims ("bins",)
    - 2D arrays with dims ("time", "bins")

    Bins are kept if:
    - vertical_scale is finite and within x_lims
    - data are valid:
        - 1D: bin value is not NaN
        - 2D: bin has at least one non-NaN value over time
    """

    if bin_dim not in da.dims:
        raise ValueError(f"da must contain the '{bin_dim}' dimension.")

    if bin_dim not in vertical_scale.dims:
        raise ValueError(f"vertical_scale must contain the '{bin_dim}' dimension.")

    if time_dim in da.dims:
        valid_bins = da.notnull().any(dim=time_dim)
    else:
        valid_bins = da.notnull()

    scale_mask = (
        vertical_scale.notnull()
        & (vertical_scale >= x_lims[0])
        & (vertical_scale <= x_lims[1])
    )

    mask = valid_bins & scale_mask

    # Required when da/vertical_scale are Dask-backed and drop=True is used
    mask = mask.compute()

    da_sliced = da.where(mask, drop=True)
    vertical_scale_sliced = vertical_scale.where(mask, drop=True)

    return da_sliced, vertical_scale_sliced, mask

def add_extra_plot_metadata(plot_metadata, norm_region_flag, 
                            stats_norm_region, maximum_channel_height):
    
    plot_metadata['norm_region_flag'] = f"{norm_region_flag}"
    for key in stats_norm_region.keys():
        plot_metadata[f"stats_{key}"] = f"{stats_norm_region[key]}"
    for key in stats_norm_region.keys():
        plot_metadata[f"masks_{key}"] = f"{stats_norm_region[key]}"
        
    plot_metadata['maximum_channel_height'] = f"{maximum_channel_height}"
    
    return(plot_metadata)

def convert_m_to_km(x):
   
    return 1E-3 * x

def get_vertical_axis_label(vertical_scale):
   
    # Convert meters to kilometers and select ranges or heights for the x axis depending on the use_dis value 
    if vertical_scale == 'range':
        x_label = "Range from the lidar [km]"
        
    elif vertical_scale == 'height_agl':
        x_label = "Height [km asl]"
        
    elif vertical_scale == 'height_asl':
        x_label = "Height [km agl]"
        
    else:
        raise Exception("Vertical scale '{vertical_scale}' not supported")

    return x_label 

def multiply_y_values(sig, sig_err, coef):
    
    # Multiply the   
    y_vals  = coef * sig.copy()
    
    y_errs = coef * sig_err.copy()
    
    return(y_vals, y_errs)


def get_normalization_factor(sig_1, sig_2, x_vals, region):
    
    mask = (x_vals >= region[0]) & (x_vals <= region[1])
    
    coef =  np.nanmean(sig_2[mask]) / np.nanmean(sig_1[mask])
    
    return coef

def insert_nan_time_gaps(sig, dim="time", gap_factor=1.5):
    """
    Insert NaN profiles around large time gaps while keeping original time values.

    Parameters
    ----------
    sig : xr.DataArray
        Input DataArray with a time dimension.
    dim : str
        Name of the time dimension.
    gap_factor : float
        A gap is detected when the time difference between two adjacent profiles
        is larger than gap_factor * minimum_time_difference.

    Returns
    -------
    sig_gap : xr.DataArray
        DataArray with extra NaN profiles inserted if gaps were detected.
        If no gaps were detected, the original sorted array is returned.
    time_res : pandas.Timedelta or None
        Inferred minimum temporal resolution.
    changed : bool
        True if NaN profiles were inserted, False otherwise.
    """

    sig = sig.sortby(dim)

    time = pd.DatetimeIndex(sig[dim].values)

    if len(time) < 2:
        return sig, None, False

    dt = time.to_series().diff().dropna()
    time_res = dt.min()

    gap_limit = gap_factor * time_res

    extra_times = []

    for t0, t1 in zip(time[:-1], time[1:]):
        if (t1 - t0) > gap_limit:
            extra_times.append(t0 + time_res)
            extra_times.append(t1 - time_res)

    if len(extra_times) == 0:
        return sig, time_res, False

    extra_times = pd.DatetimeIndex(extra_times)

    # Keep only valid extra times inside the measurement period
    extra_times = extra_times[
        (extra_times > time[0]) &
        (extra_times < time[-1])
    ]

    # Remove possible duplicates
    extra_times = extra_times.difference(time)

    if len(extra_times) == 0:
        return sig, time_res, False

    new_time = time.union(extra_times).sort_values()

    sig_gap = sig.reindex({dim: new_time})

    return sig_gap, time_res, True

def slice_time(da, t_lims, dim="time"):
    """
    Slice an xarray DataArray/Dataset along time and report if slicing happened.

    Parameters
    ----------
    da : xr.DataArray or xr.Dataset
        Input data with a time coordinate.
    t_lims : list
        [] or None means no slicing.
        Otherwise: ["yyyymmdd_hhmm", "yyyymmdd_hhmm"].
    dim : str
        Name of the time dimension.

    Returns
    -------
    da_out : xr.DataArray or xr.Dataset
        Sliced or original data.
    sliced : bool
        True if the returned array was actually sliced.
        False if the original array was returned.
    """

    if t_lims is None or len(t_lims) == 0:
        return da, False

    if len(t_lims) != 2:
        raise ValueError("t_lims must be [], None, or [start, end].")

    da_sorted = da.sortby(dim)

    time_min = pd.Timestamp(da_sorted[dim].values[0])
    time_max = pd.Timestamp(da_sorted[dim].values[-1])

    t_start = pd.to_datetime(t_lims[0], format="%Y%m%d_%H%M")
    t_end = pd.to_datetime(t_lims[1], format="%Y%m%d_%H%M")

    if t_start > t_end:
        raise ValueError("The start time in t_lims must be before the end time.")

    # Requested limits fully include the data, so nothing needs to be sliced
    if t_start <= time_min and t_end >= time_max:
        return da_sorted, False

    # Requested limits do not overlap the data at all, so keep original data
    if t_end < time_min or t_start > time_max:
        return da_sorted, False

    da_sliced = da_sorted.sel({dim: slice(t_start, t_end)})

    # If the number of time entries did not change, consider it untouched
    sliced = da_sliced.sizes[dim] != da_sorted.sizes[dim]

    return da_sliced, sliced

def add_fitting_suptitle(
    fig,
    title,
    y=1.,
    max_fontsize=12,
    min_fontsize=7,
    margin=0.04,
):
    """
    Add a suptitle and shrink its font size if it is too wide.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    title : str
        Title text.
    y : float
        Vertical position of the title in figure coordinates.
    max_fontsize : int
        Starting font size.
    min_fontsize : int
        Smallest allowed font size.
    margin : float
        Fractional horizontal margin on each side.

    Returns
    -------
    title_obj : matplotlib.text.Text
        The title text object.
    """

    title_obj = fig.suptitle(title, y=y, fontsize=max_fontsize)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    fig_width = fig.bbox.width
    allowed_width = fig_width * (1.0 - 2.0 * margin)

    fontsize = max_fontsize

    while fontsize > min_fontsize:
        title_width = title_obj.get_window_extent(renderer=renderer).width

        if title_width <= allowed_width:
            break

        fontsize -= 1
        title_obj.set_fontsize(fontsize)

        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

    return title_obj