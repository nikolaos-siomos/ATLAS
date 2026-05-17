#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modern ATLAS telecover-quadrants driver.

This module follows the newer quicklook/rayleigh architecture:

    generate_telecover_quadrants(data_pack, caller_info, settings_info)
        -> qa_test_info

The telecover input is expected to be split by sector in data_pack using keys:

    tlc_north, tlc_east, tlc_south, tlc_west

Some sector keys may be missing. The function processes the available sectors
and passes sector dictionaries to external plotting/export helpers.

Low-level operations such as axis construction, y-limit selection, plotting,
text creation, and ASCII formatting are delegated to external project
functions/classes, assumed to exist.
"""

import warnings
import numpy as np
from collections import defaultdict

from version import __version__
from visualizer import export_ascii
from utils.printouts import print_header
from visualizer.check import check_channels
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.tools import sector
from visualizer.plot_utils import (
    prepare_folder,
    slice_by_vertical_scale,
    collect_dict,
    convert_m_to_km,
    perform_color_reduction,
    add_plot_metadata,
)
from visualizer.telecover_sector_processor import TelecoverSectorProcessor

# External functions/modules assumed to exist in your modern structure.
from visualizer.generate_telecover_axis import (
    get_telecover_x_axis,
    get_telecover_y_limits,
)
from visualizer import plot_telecover

warnings.filterwarnings("ignore")


SECTOR_KEYS = {
    "N": "tlc_north",
    "E": "tlc_east",
    "S": "tlc_south",
    "W": "tlc_west",
}

SECTOR_NAMES = {
    "N": "north",
    "E": "east",
    "S": "south",
    "W": "west",
}


def _get_time_dim(da, sector_id):
    """Return the time dimension name for a sector DataArray."""

    preferred = [
        "time",
        f"time_{sector_id.lower()}",
        f"time_{SECTOR_NAMES[sector_id]}",
    ]

    for dim in preferred:
        if dim in da.dims:
            return dim

    # Fallback: use the first non-channel/non-bin dimension.
    for dim in da.dims:
        if dim not in ["channel", "bins"]:
            return dim

    raise ValueError(
        f"Could not identify a time dimension for sector {sector_id}. "
        f"Available dims: {da.dims}"
    )


def generate_telecover_quadrants(data_pack, caller_info, settings_info):
    """
    Generate quadrant telecover plots/ascii and return QA metadata.

    Parameters
    ----------
    data_pack : dict
        Modern ATLAS data package with sector-specific keys:
        'tlc_north', 'tlc_east', 'tlc_south', 'tlc_west'.
        Each available sector is expected to contain at least:
            data_pack[sector_key]['profile']
            data_pack[sector_key][caller_info['vertical_scale']]
    caller_info : dict
        Runtime/caller information, including output_folder and vertical_scale.
    settings_info : dict
        Either contains a 'tlc' sub-dictionary or is already the telecover
        settings dictionary.

    Returns
    -------
    qa_test_info : collections.defaultdict(dict)
        Nested metadata dictionary: qa_test_info['tlc'][channel].
    """

    qa_key = "tlc"
    qa_test_info = defaultdict(dict)

    available_sectors = {
        sector_id: sector_key
        for sector_id, sector_key in SECTOR_KEYS.items()
        if sector_key in data_pack
    }

    if len(available_sectors) == 0:
        return qa_test_info

    print_header("Initializing the Telecover quadrants test")

    # Prepare folders once. Telecover is one QA test composed of sector inputs.
    prepare_folder(caller_info, pattern="_tlc_")

    # Load settings.
    settings = settings_info.get("tlc", settings_info).copy()

    # Load sector profiles and vertical scales.
    profiles = {
        sector_id: data_pack[sector_key]["profile"]
        for sector_id, sector_key in available_sectors.items()
    }

    vertical_scales = {
        sector_id: convert_m_to_km(
            data_pack[sector_key][caller_info["vertical_scale"]]
        )
        for sector_id, sector_key in available_sectors.items()
    }

    # Use the first available sector as reference for common coordinates and
    # metadata. The sectors should normally share channel/bin geometry.
    ref_sector = next(iter(available_sectors.keys()))
    ref_key = available_sectors[ref_sector]

    # Only process channels present in all available sectors.
    common_channels = set(profiles[ref_sector].channel.values)
    for sector_id, da in profiles.items():
        common_channels &= set(da.channel.values)

    common_channels = np.array(sorted(common_channels))

    channels = check_channels(
        all_channels=common_channels,
        settings=settings,
    )

    if len(channels) > 0:
        qa_test_info[qa_key] = {}
        
    processor = TelecoverSectorProcessor(
        settings=settings,
        smoothing_func=None,  # or your external smoothing function
    )

    for ch in channels:
        print(f"-- channel: {ch}")

        ch_d = dict(channel=ch)

        # Gather common metadata from the reference sector.
        metadata = collect_metadata(data_pack[ref_key], atlas_channel_id=ch)

        # Select channel for all available sectors.
        sig_ch = {
            sector_id: da.sel(ch_d)
            for sector_id, da in profiles.items()
        }

        vertical_scale_ch = vertical_scales[ref_sector].sel(ch_d)

        # Slice reference sector by vertical scale and keep the mask.
        # The same bin mask is then applied to all other available sectors.
        ref_time_dim = _get_time_dim(sig_ch[ref_sector], ref_sector)

        sig_ch[ref_sector], vertical_scale_ch, sl_mask = slice_by_vertical_scale(
            da=sig_ch[ref_sector],
            vertical_scale=vertical_scale_ch,
            x_lims=settings["x_lims"],
            time_dim=ref_time_dim,
        )

        for sector_id in sig_ch:
            if sector_id == ref_sector:
                continue

            sig_ch[sector_id] = sig_ch[sector_id].where(sl_mask, drop=True)

        x_vals = vertical_scale_ch.values

        # Axis construction delegated to external modern axis function.
        x_axis_info = get_telecover_x_axis(
            x_vals=x_vals,
            x_lims=settings["x_lims"],
            x_tick=settings["x_tick"],
            vertical_scale=caller_info["vertical_scale"],
            telescope_type=str(ch)[4] if len(str(ch)) > 4 else None,
        )

        x_llim = x_axis_info["x_llim"]
        x_ulim = x_axis_info["x_ulim"]
        x_lbin = x_axis_info.get("x_lbin", 0)
        x_ubin = x_axis_info.get("x_ubin", len(x_vals) - 1)

        # Avoid mutating the shared settings dict across channels.
        channel_settings = settings.copy()
        channel_settings["smoothing_range"] = [x_llim, x_ulim]
        channel_settings["available_sectors"] = list(sig_ch.keys())

        # Use the minimum number of profiles available among the included
        # sectors, preserving the behavior of the old telecover routine.
        time_dims = {
            sector_id: _get_time_dim(da, sector_id)
            for sector_id, da in sig_ch.items()
        }

        iters = min(
            sig_ch[sector_id].sizes[time_dims[sector_id]]
            for sector_id in sig_ch
        )

        extra_sec = {
            sector_id: False
            for sector_id in sig_ch
        }

        processed = {}

        for sector_id, sig_sector in sector_arrays.items():
            processed[sector_id] = processor.process(
                x=x_vals,
                y=sig_sector.sel(channel=ch).values,
                iters=iters,
                x_sm_lims=[x_llim, x_ulim],
                region=channel_settings["normalization_region"],
            )

        # Y-axis limits delegated to an external modern axis function.
        # It receives only the sectors that are actually available.
        y_axis_info = get_telecover_y_limits(
            sig=[
                processed[sector_id]["y_m_sm"][slice(x_lbin, x_ubin + 1)]
                for sector_id in processed
            ],
            sig_nr=[
                processed[sector_id]["coef"]
                * processed[sector_id]["y_m_sm"][slice(x_lbin, x_ubin + 1)]
                for sector_id in processed
            ],
            y_lims=channel_settings["y_lims"],
        )

        # Collect metadata returned by the QA test.
        qa_test_info[qa_key][ch] = collect_dict(
            data_list=[
                iters,
                list(processed.keys()),
                extra_sec,
                x_axis_info,
                y_axis_info,
                channel_settings["normalization_region"],
            ],
            data_keys=[
                "iters",
                "available_sectors",
                "extra_sec",
                "x_axis_info",
                "y_axis_info",
                "norm_region",
            ],
        )

        # Collect metadata added to plot files.
        plot_metadata = collect_dict(
            data_list=[
                iters,
                list(processed.keys()),
                extra_sec,
                channel_settings["normalization_region"],
                __version__,
                qa_key,
            ],
            data_keys=[
                "iters",
                "available_sectors",
                "extra_sec",
                "norm_region",
                "ATLAS_version",
                "QA_test_ID",
            ],
            add_dicts=[channel_settings, metadata],
        )

#------------------------------------------------------------------------------
# Text
        lib = Libraries(
            caller_info=caller_info,
            metadata=metadata,
            extra_metadata={},
            settings=channel_settings,
            qa_test_info=qa_test_info[qa_key][ch],
        )

        text_generator = GenerateText(lib=lib)

        qa_test_info[qa_key][ch]["title"] = text_generator.make_telecover_title()
        qa_test_info[qa_key][ch]["filename"] = text_generator.make_filename(
            qa_test=qa_key
        )
        ascii_header = text_generator.make_header_telecover()

#------------------------------------------------------------------------------
# Plot
        # The modern plotter is assumed to accept a dictionary of processed
        # sectors, rather than hard-coded N/E/S/W arguments.
        qa_test_info[qa_key][ch]["tlc_plot_path"], dofl_x = plot_telecover.generate_plot(
            X=x_vals,
            sectors=processed,
            args=metadata | channel_settings | qa_test_info[qa_key][ch] | caller_info,
        )

        if dofl_x == dofl_x:
            qa_test_info[qa_key][ch]["minimum_channel_height"] = str(
                int(np.round(1e3 * dofl_x, -1))
            )
            plot_metadata["minimum_channel_height"] = qa_test_info[qa_key][ch][
                "minimum_channel_height"
            ]

        perform_color_reduction(
            color_reduction=caller_info["color_reduction"],
            plot_path=qa_test_info[qa_key][ch]["tlc_plot_path"],
        )

        add_plot_metadata(
            plot_path=qa_test_info[qa_key][ch]["tlc_plot_path"],
            plot_metadata=plot_metadata,
        )

#------------------------------------------------------------------------------
# ASCII
        sectors = {
            sector_id: processed[sector_id]["y_m"]
            for sector_id in processed
        }
        sectors_e = {
            sector_id: processed[sector_id]["y_extra"]
            for sector_id in processed
        }

        export_ascii.telecover(
            dir_out=caller_info["output_folder"],
            fname=f"{qa_test_info[qa_key][ch]['filename']}.txt",
            header=ascii_header,
            iters=1,
            alt=x_vals,
            sectors=sectors,
            sectors_e=sectors_e,
        )

    print("-----------------------------------------")
    print(" ")

    return qa_test_info
