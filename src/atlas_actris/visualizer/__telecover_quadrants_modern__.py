#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modern ATLAS telecover-quadrants driver.

This module mirrors the newer quicklook/rayleigh architecture:

    generate_telecover_quadrants(data_pack, caller_info, settings_info)
        -> qa_test_info

Notes
-----
This file intentionally keeps the high-level structure only.  Low-level
operations such as axis construction, y-limit selection, plotting, and text
creation are delegated to external project functions/classes, assumed to exist.
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

# External functions/modules assumed to exist in your modern structure.
# Keep these names aligned with your project files.
from visualizer.generate_telecover_axis import (
    get_telecover_x_axis,
    get_telecover_y_limits,
)
from visualizer import plot_telecover

warnings.filterwarnings("ignore")


def generate_telecover_quadrants(data_pack, caller_info, settings_info):
    """
    Generate quadrant telecover plots/ascii and return QA metadata.

    Parameters
    ----------
    data_pack : dict
        Modern ATLAS data package. Expected structure for each telecover key:
            data_pack[key]["sig_n"]
            data_pack[key]["sig_e"]
            data_pack[key]["sig_s"]
            data_pack[key]["sig_w"]
            data_pack[key][caller_info["vertical_scale"]]
        plus metadata arrays consumed by collect_metadata().
    caller_info : dict
        Runtime/caller information, including output_folder, process, etc.
    settings_info : dict
        Settings dictionary. Either contains a "tlc" sub-dictionary or is
        already the telecover settings dictionary.

    Returns
    -------
    qa_test_info : collections.defaultdict(dict)
        Nested metadata dictionary: qa_test_info[key][channel].
    """

    process = caller_info.get("process", [])
    process_tlc = caller_info.get("process_tlc", [])

    # Prefer an explicit telecover process list if available, otherwise infer
    # telecover-like keys from the general process/data_pack.
    if len(process_tlc) > 0:
        telecover_keys = process_tlc
    else:
        telecover_keys = [
            key for key in ["tlc", "tlc_rin", "tlc_pcb", "tlc_sec"]
            if key in process or key in data_pack
        ]

    qa_test_info = defaultdict(dict)

    for key in telecover_keys:
        if key not in data_pack:
            continue

        print_header(f"Initializing the Telecover quadrants test ({key})")

        # Prepare folders
        prepare_folder(caller_info, pattern=f"_{key}_")

        # Load arrays
        sig = {
            "N": data_pack[key]["sig_n"],
            "E": data_pack[key]["sig_e"],
            "S": data_pack[key]["sig_s"],
            "W": data_pack[key]["sig_w"],
        }

        vertical_scale = data_pack[key][caller_info["vertical_scale"]]
        vertical_scale = convert_m_to_km(vertical_scale)

        # Load settings
        settings = settings_info.get("tlc", settings_info).copy()

        # Check if the parsed channels exist
        channels = check_channels(
            all_channels=vertical_scale.channel.values,
            settings=settings,
        )

        if len(channels) > 0:
            qa_test_info[key] = {}

        for ch in channels:
            print(f"-- channel: {ch}")

            ch_d = dict(channel=ch)

            # Gather metadata common to all QA tests
            metadata = collect_metadata(data_pack[key], atlas_channel_id=ch)

            # Select channel
            sig_ch = {
                sec: sig[sec].sel(ch_d)
                for sec in ["N", "E", "S", "W"]
            }
            vertical_scale_ch = vertical_scale.sel(ch_d)

            # Slice all sectors using the same vertical-scale mask.  The first
            # sector defines the bin mask; the remaining sectors reuse it to
            # guarantee identical bins.
            sig_ch["N"], vertical_scale_ch, sl_mask = slice_by_vertical_scale(
                da=sig_ch["N"],
                vertical_scale=vertical_scale_ch,
                x_lims=settings["x_lims"],
                time_dim="time_n",
            )

            for sec, time_dim in {
                "E": "time_e",
                "S": "time_s",
                "W": "time_w",
            }.items():
                sig_ch[sec] = sig_ch[sec].where(sl_mask, drop=True)

            x_vals = vertical_scale_ch.values

            # Build/collect x-axis information externally.  This replaces the
            # old make_axis.telecover_x call and lets your modern plotter decide
            # labels/ticks/limits consistently.
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

            settings["smoothing_range"] = [x_llim, x_ulim]

            # Use the minimum number of iterations/profiles available in all
            # sectors, preserving the behavior of the old routine.
            iters = min(sig_ch[sec].sizes[sig_ch[sec].dims[0]] for sec in sig_ch)

            extra_sec = {
                "N": False,
                "E": False,
                "S": False,
                "W": False,
            }

            processed = {}

            for sec in ["N", "E", "S", "W"]:
                (
                    coef,
                    y_m,
                    y_sm,
                    y_m_sm,
                    y_l_sm,
                    y_u_sm,
                    coef_extra,
                    y_extra,
                    y_extra_sm,
                    extra_sec[sec],
                ) = sector.process(
                    x=x_vals,
                    y=sig_ch[sec].values.copy(),
                    iters=iters,
                    smooth=settings["smooth"],
                    x_sm_lims=[x_llim, x_ulim],
                    x_sm_win=settings["smoothing_window"],
                    expo=settings["smooth_exponential"],
                    region=settings["normalization_region"],
                )

                processed[sec] = {
                    "coef": coef,
                    "y_m": y_m,
                    "y_sm": y_sm,
                    "y_m_sm": y_m_sm,
                    "y_l_sm": y_l_sm,
                    "y_u_sm": y_u_sm,
                    "coef_extra": coef_extra,
                    "y_extra": y_extra,
                    "y_extra_sm": y_extra_sm,
                }

            # Y-axis limits are also delegated to a modern external function.
            y_axis_info = get_telecover_y_limits(
                sig=[
                    processed["N"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                    processed["E"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                    processed["S"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                    processed["W"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                ],
                sig_nr=[
                    processed["N"]["coef"] * processed["N"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                    processed["E"]["coef"] * processed["E"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                    processed["S"]["coef"] * processed["S"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                    processed["W"]["coef"] * processed["W"]["y_m_sm"][slice(x_lbin, x_ubin + 1)],
                ],
                y_lims=settings["y_lims"],
            )

            # Collect metadata returned by the QA test
            qa_test_info[key][ch] = collect_dict(
                data_list=[
                    iters,
                    extra_sec,
                    x_axis_info,
                    y_axis_info,
                    settings["normalization_region"],
                ],
                data_keys=[
                    "iters",
                    "extra_sec",
                    "x_axis_info",
                    "y_axis_info",
                    "norm_region",
                ],
            )

            # Collect metadata added to plot files
            plot_metadata = collect_dict(
                data_list=[
                    iters,
                    extra_sec,
                    settings["normalization_region"],
                    __version__,
                    "tlc",
                ],
                data_keys=[
                    "iters",
                    "extra_sec",
                    "norm_region",
                    "ATLAS_version",
                    "QA_test_ID",
                ],
                add_dicts=[settings, metadata],
            )

#------------------------------------------------------------------------------
# Text
            lib = Libraries(
                caller_info=caller_info,
                metadata=metadata,
                extra_metadata={},
                settings=settings,
                qa_test_info=qa_test_info[key][ch],
            )

            text_generator = GenerateText(lib=lib)

            qa_test_info[key][ch]["title"] = text_generator.make_telecover_title()
            qa_test_info[key][ch]["filename"] = text_generator.make_filename(
                qa_test="tlc"
            )
            ascii_header = text_generator.make_header_telecover()

#------------------------------------------------------------------------------
# Plot
            qa_test_info[key][ch]["tlc_plot_path"], dofl_x = plot_telecover.generate_plot(
                X=x_vals,
                YN=processed["N"],
                YE=processed["E"],
                YS=processed["S"],
                YW=processed["W"],
                args=metadata | settings | qa_test_info[key][ch] | caller_info,
            )

            if dofl_x == dofl_x:
                qa_test_info[key][ch]["minimum_channel_height"] = str(
                    int(np.round(1e3 * dofl_x, -1))
                )
                plot_metadata["minimum_channel_height"] = qa_test_info[key][ch][
                    "minimum_channel_height"
                ]

            perform_color_reduction(
                color_reduction=caller_info["color_reduction"],
                plot_path=qa_test_info[key][ch]["tlc_plot_path"],
            )

            add_plot_metadata(
                plot_path=qa_test_info[key][ch]["tlc_plot_path"],
                plot_metadata=plot_metadata,
            )

#------------------------------------------------------------------------------
# ASCII
            sectors = {
                sec: processed[sec]["y_m"]
                for sec in ["N", "E", "S", "W"]
            }
            sectors_e = {
                sec: processed[sec]["y_extra"]
                for sec in ["N", "E", "S", "W"]
            }

            export_ascii.telecover(
                dir_out=caller_info["output_folder"],
                fname=f"{qa_test_info[key][ch]['filename']}.txt",
                header=ascii_header,
                iters=1,
                alt=x_vals,
                sectors=sectors,
                sectors_e=sectors_e,
            )

    return qa_test_info
