#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:02:25 2022

@author: nick
"""

import warnings
import numpy as np
import xarray as xr
from version import __version__
from visualizer import plot_dark
from collections import defaultdict
from utils.printouts import print_header
from visualizer.check import check_channels
from scipy.stats import linregress, shapiro
from utils.error_classes import CustomWarning
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.molecular_signal_simulator import get_molecular_profile
from visualizer.plot_utils import (
    prepare_folder, 
    smoothing,
    smoothing_2D,
    collect_dict,
    convert_m_to_km, 
    slice_by_vertical_scale,
    perform_color_reduction, 
    add_plot_metadata,
    )

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

key_translation = {
    "molecular_mask_window": "fit_mask_window",
    "molecular_mask_window_step": "fit_mask_window_step",
    "molecular_mask_region": "fit_mask_region",
    "rsem_threshold": "rsem_threshold",
    "first_derivative_threshold": "first_derivative_threshold",
    "second_derivative_threshold": "second_derivative_threshold",
    "shapiro_wilk_threshold": "shapiro_wilk_threshold",
    "cross_criterion_threshold": "cross_criterion_threshold",
    "durbin_watson_threshold": "durbin_watson_threshold",
}

def _dark_title_key(key):
    """Return the base QA-test category used only for title generation.

    Dark-test variants keep their complete key everywhere else, including the
    returned metadata dictionary, PNG metadata, filenames, and ASCII output.
    """

    if key == "drk" or key.startswith("drk_"):
        return "drk"

    raise ValueError(f"Expected a dark-test key, received {key!r}.")

# def extract_arrays(ch_d, key, data_pack, data_pack_bc, data_pack_rc,
def extract_arrays(ch_d, key, data_pack_bc,
                   caller_info, settings):
    """Extract raw, background-corrected, smoothed, and raw-mean signals.

    The range-corrected signal is intentionally not used by the dark-test
    plots anymore. The corresponding extraction code is retained below as
    comments for reference.
    """
    # vertical_scale = data_pack[key][caller_info["vertical_scale"]].sel(ch_d)
    vertical_scale_bc = data_pack_bc[key][caller_info["vertical_scale"]].sel(ch_d)

    # RC signal is currently unused.
    # vertical_scale_rc = data_pack_rc[key][caller_info["vertical_scale"]].sel(ch_d)

    averaging_period = settings["averaging_period"]

    if averaging_period == "raw" or key != "drk":
        profiles_bc = data_pack_bc[key]["profile"].sel(ch_d)
        background = data_pack_bc[key]["background"].sel(ch_d)

    elif averaging_period in ("low_res", "high_res") and key == "drk":
        profiles_bc = data_pack_bc[key][f"profile_{averaging_period}"].sel(ch_d)
        background = data_pack_bc[key][f"background_{averaging_period}"].sel(ch_d)

        profiles_all_nan = not np.isfinite(_to_numpy(profiles_bc)).any()
        background_all_nan = not np.isfinite(_to_numpy(background)).any()

        if profiles_all_nan or background_all_nan:
            print()
            CustomWarning(
                f"{averaging_period} dark data are unavailable for "
                f"channel {ch_d['channel']}; raw data will be used instead."
            )
            print()
            profiles_bc = data_pack_bc[key]["profile"].sel(ch_d)
            background = data_pack_bc[key]["background"].sel(ch_d)

    else:
        raise ValueError(
            f"Unsupported averaging_period {averaging_period!r} for key {key!r}."
        )

    # Mean raw profile at the original raw temporal resolution.
    profiles_ra = data_pack_bc[key]["profile_mean"].sel(ch_d).squeeze(drop=True)

    y_dict = {
        "av": profiles_bc + background,
        "bc": profiles_bc,
        "sm": profiles_bc,
        "ra": profiles_ra,
        # "rc": profiles_rc,
    }

    x_dict = {
        "av": vertical_scale_bc,
        "bc": vertical_scale_bc,
        "sm": vertical_scale_bc,
        "ra": vertical_scale_bc,
        # "rc": vertical_scale_rc,
    }

    bins_dict = {
        "av": vertical_scale_bc.bins,
        "bc": vertical_scale_bc.bins,
        "sm": vertical_scale_bc.bins,
        "ra": vertical_scale_bc.bins,
        # "rc": vertical_scale_rc.bins,
    }

    return y_dict, x_dict, bins_dict

def insert_mol_sig(ch, ch_info, height, background, x_dict):
    """Create molecular reference profiles aligned with the signal grids."""

    emitted_wavelength = float(
        ch_info.sel(parameters="emitted_wavelength").item()
    )

    if ch[6] == "a":
        max_signal_ch = float(
            ch_info.sel(parameters="data_acquisition_range").item()
        )
    else:
        max_signal_ch = np.nan

    height_vals = _to_numpy(height).astype(float).squeeze()
    if height_vals.ndim != 1:
        raise ValueError(
            f"Expected one-dimensional height_agl, got shape {height_vals.shape}."
        )

    background_value = float(_to_numpy(background).squeeze())

    _, sig_mol, _ = get_molecular_profile(
        wavelength=emitted_wavelength,
        sig_nr=max_signal_ch * 0.4,
        sig_d=background_value,
        channel_type=ch[5],
        channel_subtype=ch[7],
        altitude_agl=height_vals,
    )

    sm_template = x_dict["sm"]
    if sig_mol.size != sm_template.size:
        raise ValueError(
            "Molecular-profile/smoothed-grid mismatch: "
            f"{sig_mol.size} molecular points versus {sm_template.size} signal bins."
        )

    if isinstance(sm_template, xr.DataArray):
        sig_mol = xr.DataArray(
            sig_mol,
            coords=sm_template.coords,
            dims=sm_template.dims,
            name="molecular_signal",
        )

    return {
        "av": xr.full_like(x_dict["av"], np.nan),
        "bc": sig_mol,
        "sm": sig_mol,
        "ra": xr.full_like(x_dict["ra"], np.nan),
        # "rc": xr.full_like(x_dict["rc"], np.nan),
    }

def _scale_range_km_to_m(region):
    """Convert a two-element range from kilometers to meters."""
    if region is None or len(region) == 0:
        return []
    return [None if value is None else 1e3 * value for value in region]


def smooth_arrays(x_dict, y_dict, settings):
    """Smooth the background-corrected signal while retaining xarray metadata."""

    y_dict["sm"], _ = smoothing_2D(
        args=settings,
        x_vals=x_dict["sm"],
        y_vals=y_dict["sm"],
        err_type="std",
    )

    # RC smoothing is intentionally disabled because RC is no longer used by
    # the dark-test plots.
    # y_dict["rc"], _ = smoothing_2D(
    #     args=settings,
    #     x_vals=x_dict["rc"],
    #     y_vals=y_dict["rc"],
    #     err_type="std",
    #     smooth_alias="smooth_rc",
    #     range_alias="smoothing_range_rc",
    #     window_alias="smoothing_window_rc",
    # )

    if not isinstance(y_dict["sm"], xr.DataArray):
        raise TypeError(
            f"smoothing_2D returned {type(y_dict['sm']).__name__} for "
            "y_dict['sm']. Install/use the xarray-aware versions of "
            "visualizer.smoothing and visualizer.plot_utils."
        )

    return y_dict

def convert_x_dict_to_km(x_dict):
    """Return a copy of the vertical-coordinate dictionary in kilometers."""
    return {name: convert_m_to_km(values) for name, values in x_dict.items()}


def generate_dark(
        # data_pack, data_pack_bc, data_pack_rc, caller_info, settings_info
        data_pack, caller_info, settings_info
        ):

    qa_test_info = defaultdict(dict)
    
    if 'drk' not in caller_info['process']:
        return
    
    for key in data_pack:
        
        if key.startswith('drk'):
            
            print_header(f'Initializing the Dark test ({key})')
            
            # Prepare folders
            prepare_folder(
                caller_info, pattern = key, 
                exclude_patterns = [f'_qck_{key}_']
                )
    
            # Load settings
            settings = settings_info.copy()
            
            # Load arrays
            system_info = data_pack[key]["system_info"]
            channel_info = data_pack[key]["channel_info"]
            
            shots = data_pack[key]["shots"]
            background = data_pack[key]["background_mean"]
        
            profiles = data_pack[key]['profile']
            
            system_info_d = dict(
                zip(
                    system_info.parameters.values,
                    system_info.values,
                )
            )
            
            # Check if the parsed channels exist and apply exclusion options
            channels = check_channels(
                all_channels = profiles.channel.values,
                settings = settings
                )
                    
            # iterate over the channels
            for ch in channels:
                print(f"-- channel: {ch}")
                
                ch_d = dict(channel = ch)

                qa_test_info.setdefault(key, {})
                qa_test_info[key].setdefault(ch, {})

                ch_info = channel_info.sel(ch_d)  
                ch_info_d = dict(zip(ch_info.parameters.values, ch_info.values))
                
                channel_settings = settings.copy()

                shots_ch = shots.sel(ch_d)
                all_shots_ch = float(
                    _to_scalar(shots_ch.sum(dim="time", skipna=True))
                )
                n_raw_profiles = int(shots_ch.sizes.get("time", 0))
                
                y_dict, x_dict, bins_dict = \
                    extract_arrays(
                        ch_d = ch_d, 
                        key = key, 
                        # data_pack = data_pack, 
                        data_pack_bc = data_pack, 
                        # data_pack_rc = data_pack_rc, 
                        caller_info = caller_info, 
                        settings = channel_settings
                        )
                
                # The extended molecular analysis is reserved for the long
                # dark measurement only: exact key ``drk`` and analog mode.
                # Photon-counting channels and all auxiliary ``drk_*``
                # measurements use the compact three-panel layout.
                extended_dark_analysis = key == "drk" and ch[6] == "a"

                if extended_dark_analysis:
                    m_dict = insert_mol_sig(
                        ch=ch,
                        ch_info=ch_info,
                        height=data_pack[key]["height_agl"].sel(ch_d),
                        background=background.sel(ch_d),
                        x_dict=x_dict,
                    )
                else:
                    m_dict = None

                
                # All downstream statistics, limits, exports, and plots use km.
                x_dict = convert_x_dict_to_km(x_dict)

                y_dict = smooth_arrays(
                    x_dict=x_dict,
                    y_dict=y_dict,
                    settings=channel_settings,
                    )
                
                zero_bin = get_zero_bin(ch_info)
            
                if channel_settings['far_range_region'][1] is None:
                    channel_settings['far_range_region'][1] = float(x_dict['av'][-1].values)
                
                # Slice the reconstructed averaged signal over the statistics
                # region. Its regional mean is used as the baseline offset,
                # because the added background may have been estimated over a
                # different range interval.
                av_region, _, _ = slice_by_vertical_scale(
                    da=y_dict["av"],
                    vertical_scale=x_dict["av"],
                    x_lims=channel_settings["stats_range"],
                )

                # Slice the mean raw profile over the statistics region for
                # noise, vertical-trend, and Gaussian-noise statistics.
                ra_region, x_region, _ = slice_by_vertical_scale(
                    da=y_dict["ra"],
                    vertical_scale=x_dict["ra"],
                    x_lims=channel_settings["stats_range"],
                )

                # Use the time-resolved background-corrected signal over the
                # same region for the temporal-trend calculation.
                bc_region, _, _ = slice_by_vertical_scale(
                    da=y_dict["bc"],
                    vertical_scale=x_dict["bc"],
                    x_lims=channel_settings["stats_range"],
                )

                qa_test_info[key][ch] = calculate_statistics(
                    av=av_region,
                    ra=ra_region,
                    bc=bc_region,
                    ranges=x_region,
                    all_shots_ch=all_shots_ch,
                    n_raw_profiles=n_raw_profiles,
                    stats=qa_test_info[key][ch],
                )

                # Gather the metadata that are common for all QA tests in a dictonary
                metadata = collect_metadata(data_pack[key], atlas_channel_id = ch)
                            
                plot_metadata = (
                    {
                        **system_info_d,
                        **ch_info_d,
                        **settings,
                        "atlas_channel_id": ch,
                        "ATLAS_version": __version__,
                        "QA_test_ID": key,
                        "QA_test_type": "drk",
                    }
                )
                
                plot_metadata = dict(sorted(plot_metadata.items()))

#------------------------------------------------------------------------------
# Dark
#------------------------------------------------------------------------------  
        
#------------------------------------------------------------------------------  
# Text
                # Load libraris
                lib = Libraries(
                    caller_info = caller_info,
                    metadata = metadata,
                    extra_metadata = {},
                    settings = settings,
                    qa_test_info = {
                        "qa_test": _dark_title_key(key),
                        "qa_test_id": key,
                    }
                    )
                
                # Call GenerateText class
                text_generator = GenerateText(lib = lib)
                
                # Make titles
                qa_test_info[key][ch]['title'] = \
                    text_generator.make_dark_title()
                
                # Make filenames
                qa_test_info[key][ch]['filename'] = text_generator.make_filename(
                    qa_test = key
                    )
    
                # Ascii header            
                # ascii_header = text_generator.make_header_dark()
    
            
                background_range = [ch_info_d['background_low_bin'], ch_info_d['background_high_bin']]

                # Raw signal plot x and y axis limits
                xlims_av, xlims_range_av, ylims_av = \
                    raw_lims(
                        y_dict['av'], 
                        bins = bins_dict['av'],
                        ranges = x_dict['av'],
                        region = background_range
                        )     
                    
                # Far-range zoom. The configured limits are expressed in
                # the same units as x_dict (km), and the corresponding bin
                # limits are derived from the selected x_dict interval.
                xlims_bc, xlims_range_bc, ylims_bc = \
                    pretrig_lims(
                        y_dict['bc'],
                        bins=bins_dict['bc'],
                        ranges=x_dict['bc'],
                        zero_bin=zero_bin,
                        region=channel_settings['far_range_region'],
                    )
                    
                # Zero bin zoomed
                xlims_zb, xlims_range_zb, ylims_zb = \
                    zero_bin_lims(
                        y_dict['bc'], 
                        bins = bins_dict['bc'],
                        ranges = x_dict['bc'],
                        zero_bin = zero_bin,
                        )    
                    
                if extended_dark_analysis:
                    # The smoothed and molecular-normalized panels are used
                    # only for the long analog ``drk`` measurement.
                    xlims_sm, xlims_range_sm, ylims_sm = \
                        smoothed_lims(
                            y_dict['sm'],
                            bins=bins_dict['sm'],
                            ranges=x_dict['sm'],
                            zero_bin=zero_bin,
                            region=channel_settings['far_range_region'],
                        )

                    # Time-resolved smoothed-BC deviation normalized by the
                    # corresponding molecular reference profile.
                    xlims_rc, ylims_rc, max_channel_vertical_scale = \
                        normalized_sm_deviation_lims(
                            sig=y_dict["sm"],
                            molecular=m_dict["sm"],
                            ranges=x_dict["sm"],
                            relative_limit=channel_settings[
                                "relative_molecular_deviation"
                            ],
                        )
                    qa_test_info[key][ch][
                        "max_channel_vertical_scale"
                    ] = max_channel_vertical_scale
                    qa_test_info[key][ch]["vertical_scale_alias"] = (
                        caller_info.get(
                            "vertical_scale_alias",
                            caller_info.get("vertical_scale", "range"),
                        )
                    )
                else:
                    # These panels are omitted entirely for photon-counting
                    # channels. Keep placeholders so the shared args assembly
                    # and plotting call remain backward compatible.
                    xlims_sm = None
                    xlims_range_sm = None
                    ylims_sm = None
                    xlims_rc = None
                    ylims_rc = None
                    qa_test_info[key][ch]["max_channel_vertical_scale"] = np.nan
                    qa_test_info[key][ch]["vertical_scale_alias"] = (
                        caller_info.get(
                            "vertical_scale_alias",
                            caller_info.get("vertical_scale", "range"),
                        )
                    )
    
                # Pass all generated scalar or list parameters relevant to the plots to the args dictionary
                qa_test_info[key][ch] = pass_to_args(
                    args = qa_test_info[key][ch], 
                    data_list = [
                        xlims_av,
                        xlims_range_av,
                        ylims_av,
                        xlims_bc,
                        xlims_range_bc,
                        ylims_bc,
                        xlims_zb,
                        xlims_range_zb,
                        ylims_zb,
                        xlims_sm,
                        xlims_range_sm,
                        ylims_sm,
                        xlims_rc,
                        ylims_rc,
                        ch[6],
                        extended_dark_analysis,
                        ],
                    data_keys = [
                        'xlims_av',
                        'xlims_range_av',
                        'ylims_av',
                        'xlims_bc',
                        'xlims_range_bc',
                        'ylims_bc',
                        'xlims_zb',
                        'xlims_range_zb',
                        'ylims_zb',
                        'xlims_sm',
                        'xlims_range_sm',
                        'ylims_sm',
                        'xlims_rc',
                        'ylims_rc',
                        'channel_mode',
                        'extended_dark_analysis',
                        ]
                    )
                
                # Make the plot
                qa_test_info[key][ch]['drk_plot_path'] = plot_dark.generate_plot(
                    bins_dict = bins_dict,
                    x_dict=x_dict,
                    y_dict=y_dict,
                    m_dict=m_dict,
                    args=channel_settings | qa_test_info[key][ch] | caller_info,
                    )  
            
                # Perform color reduction        
                perform_color_reduction(
                    color_reduction = True, 
                    plot_path = qa_test_info[key][ch]['drk_plot_path']
                    )
    
                # Add the metadata to the plot 
                add_plot_metadata(
                    plot_path = qa_test_info[key][ch]['drk_plot_path'], 
                    plot_metadata = plot_metadata
                    )
    
                # # Export to ascii (Volker's format)        
                # export_dark_ascii_blocks(
                #     dir_out=caller_info['ascii_folder'],
                #     fname=f"{qa_test_info[key][ch]['filename']}.txt",
                #     x_dict=x_dict,
                #     y_dict=y_dict,
                #     header=ascii_header,
                #     )

    return qa_test_info
                        

def slice_signal_by_range(sig, ranges, region):
    """
    Slice a 2D signal array by a range interval.

    Parameters
    ----------
    sig : ndarray (time, bins)
        Signal array.
    ranges : ndarray (bins,)
        Range coordinate corresponding to bins.
    region : list or tuple [min_range, max_range]
        Range limits (same units as ranges).

    Returns
    -------
    sig_slice : ndarray
        Sliced signal array.
    ranges_slice : ndarray
        Corresponding sliced ranges.
    """

    mask = (ranges >= region[0]) & (ranges <= region[1])

    sig_slice = sig.where(mask, drop=True)
    ranges_slice = ranges.where(mask, drop=True)

    return sig_slice, ranges_slice

def _to_numpy(value):
    """Return an in-memory NumPy representation of an xarray/NumPy object."""

    if hasattr(value, "compute"):
        value = value.compute()

    if hasattr(value, "values"):
        value = value.values

    return np.asarray(value)


def _to_scalar(value):
    """Convert a scalar-like xarray/NumPy value to a Python scalar."""

    array = _to_numpy(value)

    if array.size != 1:
        raise ValueError(
            f"Expected a scalar value, received shape {array.shape}."
        )

    return array.item()


def _as_1d_numpy(value, name):
    """Return a one-dimensional NumPy array."""

    array = _to_numpy(value).astype(float).squeeze()

    if array.ndim != 1:
        raise ValueError(
            f"Expected '{name}' to be one-dimensional, got shape {array.shape}."
        )

    return array


def _profile_mean_numpy(value, name):
    """Return a one-dimensional mean profile from 1-D or 2-D input."""

    array = _to_numpy(value).astype(float)

    if array.ndim == 1:
        return array

    if array.ndim == 2:
        return np.nanmean(array, axis=0)

    raise ValueError(
        f"Expected '{name}' to be one- or two-dimensional, got shape {array.shape}."
    )


def export_dark_ascii_blocks(dir_out, fname, header, x_dict, y_dict, y_err_dict):
    """Export dark-test data on their independent vertical grids.

    Raw/background/smoothed profiles and range-corrected profiles may have
    different bin counts. They are written as two sections in the same file
    instead of being forced into one rectangular array.
    """

    import os

    os.makedirs(dir_out, exist_ok=True)
    path = os.path.join(dir_out, fname)

    raw_columns = [
        _as_1d_numpy(x_dict["av"], "x_av"),
        _profile_mean_numpy(y_dict["av"], "y_av"),
        _as_1d_numpy(x_dict["bc"], "x_bc"),
        _profile_mean_numpy(y_dict["bc"], "y_bc"),
        _as_1d_numpy(x_dict["sm"], "x_sm"),
        _profile_mean_numpy(y_dict["sm"], "y_sm"),
    ]

    raw_lengths = {column.size for column in raw_columns}
    if len(raw_lengths) != 1:
        raise ValueError(
            "Raw, background-corrected, and smoothed ASCII columns must "
            f"share one grid. Received lengths: {sorted(raw_lengths)}"
        )

    # Range-corrected export is intentionally disabled because RC signals are
    # no longer collected by this dark-test implementation.
    # rc_columns = [...]
    # rc_body = np.column_stack(rc_columns)

    raw_body = np.column_stack(raw_columns)
    header_text = "" if header is None else str(header).rstrip()

    with open(path, "w", encoding="utf-8") as stream:
        if header_text:
            stream.write(header_text)
            stream.write("\n")

        stream.write("# RAW_BACKGROUND_SMOOTHED\n")
        stream.write(
            "# range_av_km mean_av range_bc_km mean_bc "
            "range_sm_km mean_sm\n"
        )
        np.savetxt(stream, raw_body, fmt="%.8e")

    return path


def _align_signal_with_vertical_scale(sig, vertical_scale):
    """Align a signal using its own coordinates and validate its vertical scale.

    The signal must be an xarray DataArray containing a ``time`` dimension and
    exactly one non-time dimension. Its own ``time`` coordinate is the source
    of truth; no separate time container is required.
    """

    if not isinstance(sig, xr.DataArray):
        raise TypeError(
            "Signal alignment requires an xarray.DataArray so dimension "
            f"names and coordinates can be used safely. Received {type(sig).__name__}."
        )

    if "time" not in sig.dims:
        raise ValueError(
            "Statistics signal does not contain a 'time' dimension. "
            f"Found dimensions: {sig.dims}."
        )

    if "time" not in sig.coords:
        raise ValueError("Statistics signal has no coordinate values for 'time'.")

    vertical_dims = [dim for dim in sig.dims if dim != "time"]

    if len(vertical_dims) != 1:
        raise ValueError(
            "Expected exactly one vertical dimension in addition to 'time'. "
            f"Found dimensions: {sig.dims}."
        )

    vertical_dim = vertical_dims[0]
    sig_aligned = sig.transpose("time", vertical_dim)
    n_vertical = int(vertical_scale.size)

    if sig_aligned.sizes[vertical_dim] != n_vertical:
        raise ValueError(
            "Signal/vertical-scale mismatch for the selected data entry: "
            f"signal has {sig_aligned.sizes[vertical_dim]} vertical bins, "
            f"but its x_dict entry has {n_vertical} values."
        )

    return sig_aligned

def calculate_statistics(
    av,
    ra,
    bc,
    ranges,
    all_shots_ch,
    n_raw_profiles,
    stats=None,
):
    """Calculate dark statistics in the configured statistics region.

    ``av`` is the reconstructed averaged signal (background-corrected signal
    plus the previously calculated background). Its mean within the statistics
    region defines the baseline offset, independently of the range used for the
    original background calculation.

    ``ra`` is the mean raw signal (``profile_mean``) restricted to the
    configured statistics region. Noise is converted to a single-shot
    equivalent by multiplying the standard deviation of this averaged profile
    by ``sqrt(all_shots_ch)``.
    The temporal trend is calculated independently from the time-resolved
    background-corrected signal ``bc`` over the same region.
    """
    if stats is None:
        stats = {}

    av_vals = _to_numpy(av).astype(float)
    ra_vals = _to_numpy(ra).astype(float).squeeze()
    ranges_vals = _to_numpy(ranges).astype(float).squeeze()

    if av_vals.ndim not in (1, 2):
        raise ValueError(
            "Expected the averaged signal in the statistics region to be "
            f"one- or two-dimensional, got {av_vals.shape}."
        )
    if not np.any(np.isfinite(av_vals)):
        raise ValueError(
            "The averaged signal contains no finite values in the statistics region."
        )

    if ra_vals.ndim != 1:
        raise ValueError(
            f"Expected the mean raw profile to be one-dimensional, got {ra_vals.shape}."
        )
    if ranges_vals.ndim != 1 or ranges_vals.size != ra_vals.size:
        raise ValueError(
            "Mean raw profile/range mismatch in statistics region: "
            f"signal shape {ra_vals.shape}, range shape {ranges_vals.shape}."
        )

    valid = np.isfinite(ra_vals) & np.isfinite(ranges_vals)
    ra_vals = ra_vals[valid]
    ranges_vals = ranges_vals[valid]

    if ra_vals.size <= 5:
        print(
            "--Warning: Insufficient number of points for the dark test "
            "statistics. Please check the provided stats_range parameter."
        )
        return stats

    if not np.isfinite(all_shots_ch) or all_shots_ch <= 0.0:
        raise ValueError(
            f"all_shots_ch must be finite and positive, got {all_shots_ch!r}."
        )

    av_mean = float(np.nanmean(av_vals))
    ra_mean = float(np.nanmean(ra_vals))
    ra_centered = ra_vals - ra_mean
    ra_std = float(np.nanstd(ra_vals))

    vert_fit = linregress(x=ranges_vals, y=ra_vals)

    stats["bins"] = int(ra_vals.size)
    stats["profiles"] = int(n_raw_profiles)
    stats["all_shots"] = float(all_shots_ch)
    stats["shots"] = float(all_shots_ch)  # Backward-compatible alias.
    stats["sample"] = float(all_shots_ch)  # Older table/export alias.

    stats["baseline_offset"] = av_mean
    stats["noise_per_bin_per_shot"] = float(
        ra_std * np.sqrt(all_shots_ch)
    )
    stats["noise_per_bin"] = stats["noise_per_bin_per_shot"]

    stats["vert_slope"] = float(vert_fit.slope)
    stats["vert_slope_sign"] = bool(vert_fit.pvalue <= 0.05)
    stats["gaussian_noise"] = bool(shapiro(ra_centered).pvalue > 0.05)

    stats["gaussian_noise_flag"] = (
        "Yes" if stats["gaussian_noise"] else "No"
    )
    stats["vert_slope_flag"] = (
        f"{vert_fit.slope:.2e} \u00b1 {vert_fit.stderr:.1e}"
    )

    # Temporal trend from the time-resolved background-corrected signal.
    bc = _align_signal_with_vertical_scale(bc, ranges)
    bc_vals = _to_numpy(bc).astype(float)
    time_values = _to_numpy(bc.coords["time"]).astype("datetime64[ns]")
    time_seconds = np.asarray(
        (time_values - time_values[0]) / np.timedelta64(1, "s"),
        dtype=float,
    )
    bc_mean_time = np.nanmean(bc_vals, axis=1)
    valid_time = np.isfinite(time_seconds) & np.isfinite(bc_mean_time)

    if np.count_nonzero(valid_time) > 2 and np.ptp(time_seconds[valid_time]) > 0.0:
        temp_fit = linregress(
            x=time_seconds[valid_time],
            y=bc_mean_time[valid_time],
        )
        stats["temp_slope"] = float(temp_fit.slope)
        stats["temp_slope_sign"] = bool(temp_fit.pvalue <= 0.05)
        # stats["temp_slope_flag"] = (
        #     "Signif." if stats["temp_slope_sign"] else "Insignif."
        # )
        stats["temp_slope_flag"] = f"{temp_fit.slope:.2e} \u00b1 {temp_fit.stderr:.1e}"
    else:
        stats["temp_slope"] = np.nan
        stats["temp_slope_sign"] = False
        stats["temp_slope_flag"] = "Too few valid profiles"

    return stats

def get_zero_bin(channel_info):
    """Return the canonical zero-bin value as a native Python integer."""

    zero_bin = channel_info.sel(parameters="zero_bin").item()

    return int(zero_bin)

def find_region_ind(arr, llim = None, ulim = None):
    
    zone = np.where((arr >= llim) & (arr <= llim))[0]
    
    if llim is None:
       llim = arr[0]

    if ulim is None:
       ulim = arr[-1]
      
    if not (isinstance(llim,int) or isinstance(llim,float)):
        raise Exception("llim must be integer or float")

    if not (isinstance(ulim,int) or isinstance(ulim,float)):
        raise Exception("ulim must be integer or float")

    if not llim < ulim:
        raise Exception("llim must be smaller than ulim")
        
    zone = np.where((arr >= llim) & (arr <= ulim))[0]
    
    if len(zone) > 3:
        lind = zone[0]
        uind = zone[-1]
    else:
        lind = 0
        uind = -1

    return(lind, uind)
    
def get_span(arr, lind = 0, uind = -1, edge_fraction = 0.03):

    first = arr[lind]
    last = arr[uind]
    
    span = (last - first)
    
    edge = np.ceil(edge_fraction * span)        

    return (first, last, span, edge)
    
def raw_lims(sig, bins, ranges, region):
    
    first_bin, last_bin, bin_span, bin_edge = get_span(bins.values)
    
    first_range, last_range, range_span, range_edge = get_span(ranges.values)
    
    xlims_bins = [first_bin - bin_edge, last_bin + bin_edge]
    
    xlims_range = [
        first_range - range_edge,
        last_range + range_edge
        ]
    
    min_y, max_y, mean_y, y_edge = \
        region_extrema(
            sig = sig,
            bins = bins,
            region = region
            )
    
    ylims = [mean_y - 3.3 * y_edge, mean_y + 3.3 * y_edge]
    
    return(xlims_bins, xlims_range, ylims)

def _coordinate_window_indices(coordinate, lower, upper, name):
    """Return inclusive indices for an explicitly requested coordinate window."""

    values = _to_numpy(coordinate).astype(float).squeeze()
    if values.ndim != 1 or values.size == 0:
        raise ValueError(
            f"Expected a non-empty one-dimensional {name} coordinate, got {values.shape}."
        )

    lower = float(values[0]) if lower is None else float(lower)
    upper = float(values[-1]) if upper is None else float(upper)
    if lower >= upper:
        raise ValueError(
            f"Invalid {name} limits [{lower}, {upper}]: lower must be below upper."
        )

    mask = np.isfinite(values) & (values >= lower) & (values <= upper)
    indices = np.flatnonzero(mask)
    if indices.size < 2:
        raise ValueError(
            f"The requested {name} window [{lower}, {upper}] contains fewer than "
            f"two points. Available extent is [{values[0]}, {values[-1]}]."
        )

    return int(indices[0]), int(indices[-1])


def pretrig_lims(sig, bins, ranges, zero_bin, region):
    """Return limits for the far-range background-corrected panel.

    ``region`` is interpreted directly on ``ranges`` (x_dict, in km). The
    matching bin-axis limits are obtained from the same selected indices.
    """

    lind, uind = _coordinate_window_indices(
        ranges, region[0], region[1], name="far-range"
    )

    bins_vals = _to_numpy(bins).astype(float).squeeze()
    ranges_vals = _to_numpy(ranges).astype(float).squeeze()

    xlims_bins = [float(bins_vals[lind]), float(bins_vals[uind])]
    xlims_range = [float(ranges_vals[lind]), float(ranges_vals[uind])]

    min_y, max_y, mean_y, y_edge = region_extrema(
        sig=sig,
        bins=ranges,
        region=xlims_range,
    )
    if not np.isfinite(y_edge) or y_edge == 0.0:
        y_edge = np.finfo(float).eps

    ylims = [-3.3 * y_edge, 3.3 * y_edge]
    return xlims_bins, xlims_range, ylims


def zero_bin_lims(sig, bins, ranges, zero_bin, left_bins=100, right_bins=300):
    """Return an asymmetric acquisition window around the laser-pulse position.

    ``zero_bin`` is a signed acquisition offset. A negative value means that
    recording started before the laser pulse, so the pulse occurs at stored
    bin ``-zero_bin``. A non-negative value means that recording started at or
    after the pulse; no pre-trigger samples exist and the window starts at the
    first available bin.
    """

    bins_vals = _to_numpy(bins).astype(float).squeeze()
    ranges_vals = _to_numpy(ranges).astype(float).squeeze()
    if (
        bins_vals.ndim != 1
        or ranges_vals.ndim != 1
        or bins_vals.size != ranges_vals.size
    ):
        raise ValueError(
            "Zero-bin plot requires one-dimensional bin and range coordinates "
            "with matching lengths."
        )

    finite_bins = bins_vals[np.isfinite(bins_vals)]
    if finite_bins.size == 0:
        raise ValueError("Zero-bin plot contains no finite bin coordinates.")

    bin_min = float(np.nanmin(finite_bins))
    bin_max = float(np.nanmax(finite_bins))

    if float(zero_bin) < 0.0:
        pulse_bin = -float(zero_bin)
        lower_bin = max(bin_min, pulse_bin - float(left_bins))
        upper_bin = min(bin_max, pulse_bin + float(right_bins))
    else:
        lower_bin = bin_min
        upper_bin = min(bin_max, bin_min + float(right_bins))

    lind, uind = _coordinate_window_indices(
        bins_vals, lower_bin, upper_bin, name="zero-bin"
    )

    # Keep the linked bin and vertical-scale axes aligned through the same
    # selected array indices.
    xlims_bins = [float(bins_vals[lind]), float(bins_vals[uind])]
    xlims_range = [float(ranges_vals[lind]), float(ranges_vals[uind])]

    min_y, max_y, mean_y, y_edge = region_extrema(
        sig=sig,
        bins=bins,
        region=[lower_bin, upper_bin],
    )
    if not np.isfinite(y_edge) or y_edge == 0.0:
        y_edge = np.finfo(float).eps

    ylims = [-1.1 * y_edge, 1.1 * y_edge]
    return xlims_bins, xlims_range, ylims

def smoothed_lims(sig, bins, ranges, zero_bin, region):
    """Return full-profile x limits and tightly fitted smoothed-signal y limits."""

    first_bin, last_bin, _, bin_edge = get_span(_to_numpy(bins).astype(float))
    first_range, last_range, _, range_edge = get_span(
        _to_numpy(ranges).astype(float)
    )

    xlims_bins = [first_bin - bin_edge, last_bin + bin_edge]
    xlims_range = [first_range - range_edge, last_range + range_edge]

    # Use the same mean-centred extrema approach as the raw panel, but evaluate
    # it in the explicitly configured far-range interval on x_dict.
    min_y, max_y, mean_y, y_edge = region_extrema(
        sig=sig,
        bins=ranges,
        region=region,
    )
    if not np.isfinite(y_edge) or y_edge == 0.0:
        y_edge = np.finfo(float).eps

    ylims = [mean_y - 3.3 * y_edge, mean_y + 3.3 * y_edge]
    return xlims_bins, xlims_range, ylims

def normalized_sm_deviation_lims(
    sig, molecular, ranges, relative_limit, padding=0.10
):
    """Return plot limits and the first threshold-exceedance range."""

    sig = _align_signal_with_vertical_scale(sig, ranges)
    vertical_dim = next(dim for dim in sig.dims if dim != "time")

    if not isinstance(molecular, xr.DataArray):
        molecular = xr.DataArray(
            _to_numpy(molecular).astype(float).squeeze(),
            dims=(vertical_dim,),
            coords={vertical_dim: sig[vertical_dim]},
        )
    else:
        molecular = molecular.squeeze(drop=True)
        if molecular.ndim != 1:
            raise ValueError(
                "Expected a one-dimensional molecular profile, "
                f"got dimensions {molecular.dims}."
            )
        molecular = molecular.rename({molecular.dims[0]: vertical_dim})
        molecular = molecular.assign_coords({vertical_dim: sig[vertical_dim]})

    if molecular.sizes[vertical_dim] != sig.sizes[vertical_dim]:
        raise ValueError(
            "Molecular-profile/signal mismatch while calculating limits: "
            f"{molecular.sizes[vertical_dim]} molecular points versus "
            f"{sig.sizes[vertical_dim]} signal bins."
        )

    deviation = sig - sig.mean(dim="time", skipna=True)
    valid_molecular = np.isfinite(molecular) & (molecular != 0.0)
    normalized = (deviation / molecular).where(valid_molecular)

    ranges_vals = _to_numpy(ranges).astype(float).squeeze()
    normalized_vals = _to_numpy(normalized).astype(float)
    molecular_vals = _to_numpy(molecular).astype(float).squeeze()

    valid_columns = (
        np.isfinite(ranges_vals)
        & np.isfinite(molecular_vals)
        & (molecular_vals != 0.0)
        & np.any(np.isfinite(normalized_vals), axis=0)
    )
    if not np.any(valid_columns):
        raise ValueError(
            "No finite normalized smoothed-BC deviations are available. "
            "Check the molecular reference profile and signal grid."
        )

    valid_ranges = ranges_vals[valid_columns]
    first_range, last_range, _, range_edge = get_span(valid_ranges)
    xlims_range = [first_range - range_edge, last_range + range_edge]

    relative_limit = float(relative_limit)
    if not np.isfinite(relative_limit) or relative_limit <= 0.0:
        raise ValueError(
            "relative_limit must be a finite positive float, "
            f"got {relative_limit!r}."
        )

    exceeds_by_column = np.any(
        np.isfinite(normalized_vals)
        & (np.abs(normalized_vals) > relative_limit),
        axis=0,
    )
    crossing_indices = np.flatnonzero(valid_columns & exceeds_by_column)

    if crossing_indices.size:
        # "First" follows the stored x-axis order used by the plot.
        max_channel_vertical_scale = float(
            ranges_vals[int(crossing_indices[0])]
        )
    else:
        max_channel_vertical_scale = np.nan

    finite_values = normalized_vals[:, valid_columns]
    finite_values = finite_values[np.isfinite(finite_values)]
    max_abs = float(np.nanmax(np.abs(finite_values)))
    if max_abs == 0.0:
        max_abs = np.finfo(float).eps

    limit = (1.0 + float(padding)) * max_abs
    return xlims_range, [-limit, limit], max_channel_vertical_scale

def region_extrema(sig, bins, region):

    bins_vals = np.asarray(
        bins.values if hasattr(bins, "values") else bins,
        dtype=float,
    )

    sig_vals = np.asarray(
        sig.values if hasattr(sig, "values") else sig,
        dtype=float,
    )

    if sig_vals.ndim == 1:
        sig_vals = sig_vals[np.newaxis, :]

    mask_bins = (
        (bins_vals >= region[0])
        & (bins_vals <= region[1])
    )

    if not np.any(mask_bins):
        raise ValueError(
            f"No bins found inside region {region}. "
            f"Available bin range is {bins_vals[0]} to {bins_vals[-1]}."
        )

    sig_region = sig_vals[:, mask_bins]

    min_y = np.nanmin(sig_region)
    max_y = np.nanmax(sig_region)
    mean_y = np.nanmean(sig_region)

    edge = max(mean_y - min_y, max_y - mean_y)

    return min_y, max_y, mean_y, edge
    
    
def pass_to_args(args, data_list, data_keys):
    
    for i in range(len(data_keys)):
        args[data_keys[i]] = data_list[i]
        
    return(args)

def add_extra_plot_metadata(plot_metadata, norm_region_flag, 
                            stats_norm_region, maximum_channel_height):
    
    plot_metadata['norm_region_flag'] = f"{norm_region_flag}"
    for key in stats_norm_region.keys():
        plot_metadata[f"stats_{key}"] = f"{stats_norm_region[key]}"
    for key in stats_norm_region.keys():
        plot_metadata[f"masks_{key}"] = f"{stats_norm_region[key]}"
        
    plot_metadata['maximum_channel_height'] = f"{maximum_channel_height}"
    
    return(plot_metadata)
