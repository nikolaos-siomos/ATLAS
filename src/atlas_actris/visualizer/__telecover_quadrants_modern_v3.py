#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modern ATLAS telecover-quadrants driver.

Expected input data_pack keys are sector-specific:

    tlc_north, tlc_east, tlc_south, tlc_west

There is no parent ``tlc`` data_pack entry. Missing sector keys are allowed.
The function treats telecover as one QA test, returned under qa_test_info['tlc'].

This module intentionally keeps the high-level workflow here, while delegating
low-level operations such as axis creation, sector processing, text creation,
plotting, ASCII export, color reduction, and PNG metadata writing to external
project helpers.
"""

import warnings
from collections import defaultdict

import numpy as np

from version import __version__
from visualizer import export_ascii, plot_telecover
from utils.printouts import print_header
from visualizer.check import check_channels
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.telecover_sector_processor import TelecoverSectorProcessor
from visualizer.plot_utils import (
    prepare_folder,
    slice_by_vertical_scale,
    collect_dict,
    convert_m_to_km,
    perform_color_reduction,
    add_plot_metadata,
)

# External modern axis helpers assumed to exist.
from visualizer.generate_telecover_axis import (
    get_telecover_x_axis,
    get_telecover_y_limits,
)

warnings.filterwarnings("ignore")


QA_KEY = "tlc"

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


def _available_sector_keys(data_pack):
    """Return mapping like {'N': 'tlc_north'} for sectors present in data_pack."""

    return {
        sector_id: data_key
        for sector_id, data_key in SECTOR_KEYS.items()
        if data_key in data_pack
    }


def _get_time_dim(da, sector_id):
    """Return the time dimension name of a sector DataArray."""

    candidates = [
        "time",
        f"time_{sector_id.lower()}",
        f"time_{SECTOR_NAMES[sector_id]}",
    ]

    for dim in candidates:
        if dim in da.dims:
            return dim

    for dim in da.dims:
        if dim not in ["channel", "bins"]:
            return dim

    raise ValueError(
        f"Could not identify a time dimension for sector {sector_id}. "
        f"Available dims: {da.dims}"
    )


def _common_channels(sector_profiles):
    """Return channels common to all available sector profiles."""

    first_sector = next(iter(sector_profiles.keys()))
    channels = set(sector_profiles[first_sector].channel.values)

    for da in sector_profiles.values():
        channels &= set(da.channel.values)

    return np.array(sorted(channels))


def _slice_sector_profiles(sector_profiles_ch, vertical_scale_ch, settings, ref_sector):
    """
    Slice all sector profiles with one bin mask from the reference sector.

    The reference sector is sliced using slice_by_vertical_scale. The returned
    bin mask is then applied to all other available sectors.
    """

    ref_time_dim = _get_time_dim(sector_profiles_ch[ref_sector], ref_sector)

    sector_profiles_ch[ref_sector], vertical_scale_ch, sl_mask = slice_by_vertical_scale(
        da=sector_profiles_ch[ref_sector],
        vertical_scale=vertical_scale_ch,
        x_lims=settings["x_lims"],
        time_dim=ref_time_dim,
    )

    for sector_id in sector_profiles_ch:
        if sector_id == ref_sector:
            continue

        sector_profiles_ch[sector_id] = sector_profiles_ch[sector_id].where(
            sl_mask,
            drop=True,
        )

    return sector_profiles_ch, vertical_scale_ch, sl_mask


def _sector_iters(sector_profiles_ch):
    """Return the common number of profiles among available sectors."""

    return min(
        da.sizes[_get_time_dim(da, sector_id)]
        for sector_id, da in sector_profiles_ch.items()
    )


def generate_telecover_quadrants(data_pack, caller_info, settings_info):
    """
    Generate telecover quadrant QA outputs.

    Parameters
    ----------
    data_pack : dict
        Data package containing some or all of:
        'tlc_north', 'tlc_east', 'tlc_south', 'tlc_west'.
        Each available sector entry is expected to contain:
        - ['profile']
        - [caller_info['vertical_scale']]
    caller_info : dict
        Runtime information, including output folder and vertical scale choice.
    settings_info : dict
        Either the telecover settings dictionary directly, or a parent dict
        containing settings_info['tlc'].

    Returns
    -------
    qa_test_info : collections.defaultdict(dict)
        Nested dictionary with one top-level telecover key:
        qa_test_info['tlc'][channel].
    """

    qa_test_info = defaultdict(dict)

    available_sectors = _available_sector_keys(data_pack)

    if len(available_sectors) == 0:
        return qa_test_info

    print_header("Initializing the Telecover quadrants test")

    # Telecover is one QA test, not a loop over multiple QA keys.
    prepare_folder(caller_info, pattern="_tlc_")

    settings = settings_info.get(QA_KEY, settings_info).copy()

    # Collect available sector DataArrays into one explicit dictionary:
    # {'N': data_pack['tlc_north']['profile'], ...}
    sector_profiles = {
        sector_id: data_pack[data_key]["profile"]
        for sector_id, data_key in available_sectors.items()
    }

    vertical_scales = {
        sector_id: convert_m_to_km(
            data_pack[data_key][caller_info["vertical_scale"]]
        )
        for sector_id, data_key in available_sectors.items()
    }

    # Use the first available sector as reference for metadata and bin geometry.
    ref_sector = next(iter(available_sectors.keys()))
    ref_key = available_sectors[ref_sector]

    channels = check_channels(
        all_channels=_common_channels(sector_profiles),
        settings=settings,
    )

    if len(channels) > 0:
        qa_test_info[QA_KEY] = {}

    for ch in channels:
        print(f"-- channel: {ch}")

        ch_d = dict(channel=ch)
        ch_key = str(ch)

        metadata = collect_metadata(data_pack[ref_key], atlas_channel_id=ch)

        sector_profiles_ch = {
            sector_id: da.sel(ch_d)
            for sector_id, da in sector_profiles.items()
        }

        vertical_scale_ch = vertical_scales[ref_sector].sel(ch_d)

        sector_profiles_ch, vertical_scale_ch, sl_mask = _slice_sector_profiles(
            sector_profiles_ch=sector_profiles_ch,
            vertical_scale_ch=vertical_scale_ch,
            settings=settings,
            ref_sector=ref_sector,
        )

        x_vals = vertical_scale_ch.values

        # Axis creation is delegated to your modern telecover-axis helper.
        x_axis_info = get_telecover_x_axis(
            x_vals=x_vals,
            x_lims=settings["x_lims"],
            x_tick=settings["x_tick"],
            vertical_scale=caller_info["vertical_scale"],
            telescope_type=ch_key[4] if len(ch_key) > 4 else None,
        )

        x_lbin = x_axis_info.get("x_lbin", 0)
        x_ubin = x_axis_info.get("x_ubin", len(x_vals) - 1)
        x_llim = x_axis_info.get("x_llim", settings["x_lims"][0])
        x_ulim = x_axis_info.get("x_ulim", settings["x_lims"][1])

        # Keep channel-specific modifications local.
        channel_settings = settings.copy()
        channel_settings["smoothing_range"] = [x_llim, x_ulim]
        channel_settings["available_sectors"] = list(sector_profiles_ch.keys())

        iters = _sector_iters(sector_profiles_ch)

        sector_processor = TelecoverSectorProcessor(
            settings=channel_settings,
        )

        processed = {}
        extra_sec = {}

        for sector_id, da in sector_profiles_ch.items():
            processed[sector_id] = sector_processor.process(
                x=x_vals,
                y=da.values.copy(),
                iters=iters,
                x_sm_lims=[x_llim, x_ulim],
                region=channel_settings["normalization_region"],
            )

            extra_sec[sector_id] = processed[sector_id]["has_extra"]

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

        # Metadata returned by the QA test.
        qa_test_info[QA_KEY][ch] = collect_dict(
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

        # Metadata added to the PNG file.
        plot_metadata = collect_dict(
            data_list=[
                iters,
                list(processed.keys()),
                extra_sec,
                channel_settings["normalization_region"],
                __version__,
                QA_KEY,
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
            qa_test_info=qa_test_info[QA_KEY][ch],
        )

        text_generator = GenerateText(lib=lib)

        qa_test_info[QA_KEY][ch]["title"] = text_generator.make_telecover_title()
        qa_test_info[QA_KEY][ch]["filename"] = text_generator.make_filename(
            qa_test=QA_KEY,
        )

        ascii_header = text_generator.make_header_telecover()

#------------------------------------------------------------------------------
# Plot
        qa_test_info[QA_KEY][ch]["tlc_plot_path"], dofl_x = plot_telecover.generate_plot(
            X=x_vals,
            sectors=processed,
            args=metadata | channel_settings | qa_test_info[QA_KEY][ch] | caller_info,
        )

        if dofl_x == dofl_x:
            qa_test_info[QA_KEY][ch]["minimum_channel_height"] = str(
                int(np.round(1e3 * dofl_x, -1))
            )
            plot_metadata["minimum_channel_height"] = qa_test_info[QA_KEY][ch][
                "minimum_channel_height"
            ]

        perform_color_reduction(
            color_reduction=caller_info["color_reduction"],
            plot_path=qa_test_info[QA_KEY][ch]["tlc_plot_path"],
        )

        add_plot_metadata(
            plot_path=qa_test_info[QA_KEY][ch]["tlc_plot_path"],
            plot_metadata=plot_metadata,
        )

#------------------------------------------------------------------------------
# ASCII
        sectors = {
            sector_id: result["y_m"]
            for sector_id, result in processed.items()
        }

        sectors_e = {
            sector_id: result["y_extra"]
            for sector_id, result in processed.items()
        }

        export_ascii.telecover(
            dir_out=caller_info["output_folder"],
            fname=f"{qa_test_info[QA_KEY][ch]['filename']}.txt",
            header=ascii_header,
            iters=1,
            alt=x_vals,
            sectors=sectors,
            sectors_e=sectors_e,
        )

    print("-----------------------------------------")
    print(" ")

    return qa_test_info
