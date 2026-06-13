#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modern ATLAS telecover-quadrants driver.

Expected input data_pack keys are sector-specific:

    tlc_outer, tlc_inner

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
from visualizer import export_ascii, plot_ring_telecover
from utils.printouts import print_header
from visualizer.check import check_channels
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.telecover_sector_processor import TelecoverSectorProcessor
from visualizer.plot_utils import (
    prepare_folder,
    collect_dict,
    convert_m_to_km,
    perform_color_reduction,
    add_plot_metadata,
)

warnings.filterwarnings("ignore")


QA_KEY = "tlc_rin"

SECTOR_KEYS = {
    "O": "tlc_outer",
    "I": "tlc_inner",
}

SECTOR_NAMES = {
    "O": "outer",
    "I": "inner",
}


def _available_sector_keys(data_pack):
    """Return mapping like {'O': 'tlc_inner'} for sectors present in data_pack."""

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

def _sector_iters(sector_profiles_ch):
    """Return the common number of profiles among available sectors."""

    return min(
        da.sizes[_get_time_dim(da, sector_id)]
        for sector_id, da in sector_profiles_ch.items()
    )


def generate_ring_telecover(data_pack, caller_info, settings_info):
    """
    Generate telecover quadrant QA outputs.

    Parameters
    ----------
    data_pack : dict
        Data package containing some or all of:
        'tlc_outer', 'tlc_inner'.
        Each available sector entry is expected to contain:
        - ['profile']
        - [caller_info['vertical_scale']]
    caller_info : dict
        Runtime information, including output folder and vertical scale choice.
    settings_info : dict
        Either the telecover settings dictionary directly, or a parent dict
        containing settings_info['tlc_rin'].

    Returns
    -------
    qa_test_info : collections.defaultdict(dict)
        Nested dictionary with one top-level telecover key:
        qa_test_info['tlc_rin'][channel].
    """

    qa_test_info = defaultdict(dict)

    available_sectors = _available_sector_keys(data_pack)

    if len(available_sectors) == 0:
        return qa_test_info

    print_header("Initializing the ring Telecover test")

    # Telecover is one QA test, not a loop over multiple QA keys.
    prepare_folder(caller_info, pattern = "_tlc_rin_", exclude_pattern = '_qck_tlc_rin_')

    settings = settings_info.copy()

    # Load common arrays
    tlc_common_key = next(iter(data_pack))
    system_info = data_pack[tlc_common_key]["system_info"]
    channel_info = data_pack[tlc_common_key]["channel_info"]
          
    # Collect available sector DataArrays into one explicit dictionary:
    # {'I': data_pack['tlc_inner']['profile'], ...}
    sector_profiles = {
        sector_id: data_pack[data_key]["profile"].persist()
        for sector_id, data_key in available_sectors.items()
    }

    vertical_scales = {
        sector_id: convert_m_to_km(
            data_pack[data_key][caller_info["vertical_scale"]]
        )
        for sector_id, data_key in available_sectors.items()
    }
    
    ranges = {
        sector_id: convert_m_to_km(
            data_pack[data_key]['range']
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

    sys_info = dict(zip(system_info.parameters.values, system_info.values))

    if len(channels) > 0:
        qa_test_info[QA_KEY] = {}

    for ch in channels:
        print(f"-- channel: {ch}")

        ch_d = dict(channel=ch)

        ch_info = channel_info.sel({'channel':ch})  
        ch_info_d = dict(zip(ch_info.parameters.values, ch_info.values))

        metadata = collect_metadata(data_pack[ref_key], atlas_channel_id=ch)

        sector_profiles_ch = {
            sector_id: da.sel(ch_d)
            for sector_id, da in sector_profiles.items()
        }

        x_vals = vertical_scales[ref_sector].sel(ch_d)
        ranges = vertical_scales[ref_sector].sel(ch_d)

        # Keep channel-specific modifications local.
        channel_settings = settings.copy()
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
                region=channel_settings["normalization_region"],
            )

            extra_sec[sector_id] = processed[sector_id]["has_extra"]

        # Metadata returned by the QA test.
        qa_test_info[QA_KEY][ch] = collect_dict(
            data_list=[
                iters,
                list(processed.keys()),
                extra_sec,
                channel_settings["normalization_region"],
            ],
            data_keys=[
                "iters",
                "available_sectors",
                "extra_sec",
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
        qa_test_info[QA_KEY][ch]["tlc_rin_plot_path"], dofl_x = \
            plot_ring_telecover.generate_ring_telecover(
                X=x_vals,
                sectors=processed,
                ranges=ranges,
                args=metadata | channel_settings | qa_test_info[QA_KEY][ch] | caller_info,
            )

        if dofl_x == dofl_x:
            qa_test_info[QA_KEY][ch]["minimum_channel_height"] = str(
                int(np.round(1e3 * dofl_x, -1))
            )
            
        minimum_channel_height = qa_test_info[QA_KEY][ch].get("minimum_channel_height")
        
        plot_metadata = (
            {
                **sys_info,
                **ch_info_d,
                **settings,
                "atlas_channel_id": ch,
                "ATLAS_version": __version__,
                "QA_test_ID": "tlc_qua",
                "minimum_channel_height": minimum_channel_height
            }
        )

        plot_metadata = dict(sorted(plot_metadata.items()))


        perform_color_reduction(
            color_reduction=caller_info["color_reduction"],
            plot_path=qa_test_info[QA_KEY][ch]["tlc_rin_plot_path"],
        )

        add_plot_metadata(
            plot_path=qa_test_info[QA_KEY][ch]["tlc_rin_plot_path"],
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
