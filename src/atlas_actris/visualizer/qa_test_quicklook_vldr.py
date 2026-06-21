#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import warnings
import numpy as np
from version import __version__
from collections import defaultdict
from utils.printouts import print_header
from processor.packaging import collect_metadata
from visualizer.check import check_vldr_pairs
from visualizer.plot_vldr import generate_plot
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder, smoothing_2D, collect_dict,
    convert_m_to_km, perform_color_reduction, 
    add_plot_metadata, insert_nan_time_gaps, slice_time
    )

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')


def _get_vertical_bin_dim(vertical_scale):
    """Return the vertical/bin dimension name from a 1D vertical scale."""

    if len(vertical_scale.dims) != 1:
        raise ValueError(
            "The selected vertical scale must be 1D after selecting one channel. "
            f"Found dimensions: {vertical_scale.dims}"
        )

    return vertical_scale.dims[0]


def _slice_vertical_scale_only(sig_ch, vertical_scale_ch, x_lims):
    """
    Slice one channel lazily using only the eager vertical scale.

    This intentionally does not inspect the profile values.  Inspecting profile
    validity with operations such as da.notnull().any(dim='time') would trigger
    a calculation over the lazy profile.  For quicklooks, filtering finite
    vertical coordinates within x_lims is enough and keeps pcolormesh happy.
    """

    bin_dim = _get_vertical_bin_dim(vertical_scale_ch)
    x_vals_all = np.asarray(vertical_scale_ch.values)

    if x_lims is None or len(x_lims) == 0:
        mask = np.isfinite(x_vals_all)
    else:
        mask = (
            np.isfinite(x_vals_all)
            & (x_vals_all >= x_lims[0])
            & (x_vals_all <= x_lims[1])
        )

    if not np.any(mask):
        selected_id = None
        for coord_name in ["channel", "pair"]:
            if coord_name in sig_ch.coords:
                try:
                    selected_id = sig_ch.coords[coord_name].values
                except Exception:
                    selected_id = sig_ch.coords[coord_name]
                break

        raise ValueError(
            "No finite vertical-scale bins were found inside x_lims "
            f"({x_lims}) for selection {selected_id}."
        )

    sig_ch = sig_ch.isel({bin_dim: mask})
    vertical_scale_ch = vertical_scale_ch.isel({bin_dim: mask})

    return sig_ch, vertical_scale_ch


def _to_numpy_selected(da):
    """
    Materialize only the already-selected quicklook slice.

    Matplotlib cannot draw a Dask-backed array directly, so a computation is
    still necessary.  The important point is that channel/time/bin selection has
    already happened before this function is called.
    """

    if hasattr(da, "compute"):
        da = da.compute()

    return np.asarray(da.values)
 
def _select_vertical_scale_for_pair(vertical_scale, pair, ch_r):
    """Select the 1D vertical scale for one VLDR pair."""

    if "pair" in vertical_scale.dims:
        return vertical_scale.sel(pair=pair)

    if "channel" in vertical_scale.dims:
        return vertical_scale.sel(channel=ch_r)

    if len(vertical_scale.dims) == 1:
        return vertical_scale

    raise ValueError(
        "Cannot select vertical scale for VLDR quicklook. Expected a 'channel' "
        "dimension, a 'pair' dimension, or an already 1D vertical scale. "
        f"Found dimensions: {vertical_scale.dims}."
    )


def _pair_or_channel_metadata(data_pack_key, pair, ch_r):
    """
    Collect metadata for plotting.

    collect_metadata is channel-based in the existing quicklook path.  For VLDR
    products we use ch_r as the channel metadata anchor and add pair metadata
    separately in plot_metadata.
    """

    try:
        return collect_metadata(data_pack_key, atlas_channel_id=ch_r)
    except Exception:
        return {}


def _safe_pair_for_filename(pair):
    """Make a pair id safe to append to a filename."""

    return str(pair).replace("/", "-").replace(" ", "_")


def generate_vldr_quicklooks(data_pack, caller_info, settings_info):
    """
    Generate quicklook plots for time-resolved VLDRs.

    This routine is analogous to generate_quicklooks(), but it reads
    data_pack[key]['pol_cal_ratio'] instead of data_pack[key]['profile'] and
    iterates over the 'pair' dimension instead of the 'channel' dimension.

    Only ratios marked as ratio_type == 'vldr' in pol_cal_info are plotted.
    The data remain lazy until a single time/pair/bin slice has been selected
    and trimmed by the eager vertical scale.
    """

    process_qck = caller_info['process_qck']

    qa_test_info = defaultdict(dict)

    for key in process_qck:
        if key not in data_pack.keys():
            continue

        if 'pol_cal_ratio' not in data_pack[key]:
            continue

        if caller_info['vertical_scale'] not in data_pack[key]:
            raise KeyError(
                f"Missing vertical scale '{caller_info['vertical_scale']}' "
                f"for VLDR quicklook key {key}."
            )

        print_header(f"Start generating VLDR quicklooks ({key})")

        # Use a distinct cleanup pattern so profile quicklooks are not removed.
        prepare_folder(caller_info, pattern="_qck_vldr_")

        ratios = data_pack[key]['pol_cal_ratio']
        pol_cal_info = data_pack[key].get('pol_cal_info')
        vertical_scale = data_pack[key][caller_info['vertical_scale']]

        settings = settings_info.copy()

        # Slice time lazily.
        ratios, time_sliced = slice_time(ratios, t_lims=settings['t_lims'])

        # Insert NaN profiles around time gaps. The new time coordinates are
        # finite; NaNs are inserted into ratio values only.
        ratios, _, has_time_gap = insert_nan_time_gaps(ratios, gap_factor=1.5)
        time = ratios.time.values

        # Convert range/height to km. vertical_scale is eager.
        vertical_scale = convert_m_to_km(vertical_scale)

        # Check if the parsed channels exist and apply exclusion options
        pair_records = check_vldr_pairs(
            pol_cal_ratio=ratios,
            pol_cal_info=pol_cal_info,
            settings=settings,
        )

        if len(pair_records) == 0:
            print(f"-- no VLDR pairs found for {key}")
            print('-----------------------------------------')
            print(' ')
            continue

        qa_test_info[key] = {}

        for record in pair_records:
            pair = record['pair']
            ch_r = record['ch_r']
            ch_t = record['ch_t']

            print(f"-- VLDR pair: {pair}")

            pair_d = dict(pair=pair)

            ratio_pair = ratios.sel(pair_d)
            vertical_scale_pair = _select_vertical_scale_for_pair(
                vertical_scale=vertical_scale,
                pair=pair,
                ch_r=ch_r,
            )

            # Trim bins before materializing the lazy ratio.  This uses only the
            # eager vertical scale and therefore avoids computing over the full
            # time-resolved pol_cal_ratio product.
            ratio_pair, vertical_scale_pair = _slice_vertical_scale_only(
                sig_ch=ratio_pair,
                vertical_scale_ch=vertical_scale_pair,
                x_lims=settings['x_lims'],
            )

            # Matplotlib and the current smoothing functions need NumPy arrays.
            # Compute only this selected 2D time/pair/bin quicklook slice.
            y_vals = _to_numpy_selected(ratio_pair)
            x_vals = np.asarray(vertical_scale_pair.values)

            y_vals_sm, _ = smoothing_2D(
                args=settings,
                x_vals=x_vals,
                y_vals=y_vals,
                err_type="std",
            )

            metadata = _pair_or_channel_metadata(
                data_pack_key=data_pack[key],
                pair=pair,
                ch_r=ch_r,
            )

            qa_test_info[key][pair] = {
                'time_sliced': time_sliced,
                'has_time_gap': has_time_gap,
                'pair': pair,
                'ch_r': ch_r,
                'ch_t': ch_t,
                'ratio_type': 'vldr',
            }


            plot_metadata = (
                {
                    # **sys_info,
                    # **ch_info_d,
                    **settings,
                    "atlas_channel_id_r": ch_r,
                    "atlas_channel_id_t": ch_t,
                    "vldr_id": pair,
                    "ATLAS_version": __version__,
                    "QA_test_ID": f"qck_{key}",
                }
            )
            
            lib = Libraries(
                caller_info=caller_info,
                metadata=metadata,
                extra_metadata={},
                settings=settings,
                qa_test_info=qa_test_info[key][pair],
            )

            text_generator = GenerateText(lib=lib)

            title = text_generator.make_vldr_title()
            # try:
            #     title = text_generator.make_quicklook_title()
            # except Exception:
            #     title = f"VLDR quicklook {key} - {pair}"

            # if ch_t is None:
            #     title = f"VLDR {pair} ({ch_r})"
            # else:
            #     title = f"VLDR {pair} ({ch_r} / {ch_t})"

            # try:
            filename = text_generator.make_filename_pair(qa_test='qck_vldr')
            # except Exception:
                # filename = 'qck_vldr'

            qa_test_info[key][pair]['title'] = title
            qa_test_info[key][pair]['filename'] = filename

            qa_test_info[key][pair]['qck_plot_path'] = generate_plot(
                T=time,
                X=x_vals,
                Y=y_vals_sm,
                args=settings | qa_test_info[key][pair] | caller_info,
            )

            perform_color_reduction(
                color_reduction=True,
                plot_path=qa_test_info[key][pair]['qck_plot_path'],
            )

            add_plot_metadata(
                plot_path=qa_test_info[key][pair]['qck_plot_path'],
                plot_metadata=plot_metadata,
            )

        print('-----------------------------------------')
        print(' ')

