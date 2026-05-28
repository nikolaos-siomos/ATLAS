#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
New-style polarization calibration visualizer.

This module consumes polarization calibration products that have already been
created by the processor, typically from a stage such as ``vldr_generated``.

Expected call pattern
---------------------
from visualizer.__polarization_calibration__ import generate_polarization_calibration

pol_cal__metadata = generate_polarization_calibration(
    data_pack=processor.export_test_from_stage("vldr_generated"),
    caller_info=processor.processing_info["caller_info"],
    settings_info=processor.settings_info,
)

Main assumptions
----------------
- data_pack["pcb"] contains combined calibration products:
    pol_cal_ratio, pol_cal_ratio_error, pol_cal_info
- data_pack["pcb_p45"] and data_pack["pcb_m45"] contain source calibration
  products and metadata/vertical scale.
- data_pack["ray_pcb"] contains Rayleigh-side calibration products and metadata.
- The left/calibration panel plots gain_ratio products, not eta products:
    gain_ratio, gain_ratio_p45, gain_ratio_m45
- eta entries are used only as the valid pair basis.
"""

from collections import defaultdict
import warnings

import numpy as np

from version import __version__
from utils.printouts import print_header
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder,
    slice_by_vertical_scale,
    smoothing,
    convert_m_to_km,
    perform_color_reduction,
    add_plot_metadata,
)
from visualizer import plot_polarization_calibration


warnings.filterwarnings("ignore")


def generate_polarization_calibration(data_pack, caller_info, settings_info):
    """
    Generate polarization calibration plots from processor products.

    Required data_pack entries
    --------------------------
    data_pack["pcb"]:
        pol_cal_ratio
        pol_cal_ratio_error
        pol_cal_info

    data_pack["pcb_p45"]:
        pol_cal_ratio
        vertical scale and metadata

    data_pack["pcb_m45"]:
        pol_cal_ratio
        vertical scale and metadata

    data_pack["ray_pcb"]:
        pol_cal_ratio
        pol_cal_ratio_error
        pol_cal_info
        vertical scale and metadata

    Optional
    --------
    data_pack["ray_pcb"]["molecular_ratio"]
    data_pack["ray_pcb"]["molecular_info"]

    Returns
    -------
    defaultdict(dict)
        QA-test metadata, indexed as qa_test_info["pcb"][eta_id].
    """

    qa_test_info = defaultdict(dict)

    required_keys = ["pcb", "pcb_m45", "pcb_p45", "ray_pcb"]
    if any(key not in data_pack for key in required_keys):
        return qa_test_info

    print_header("Initializing the Polarization Calibration test")

    prepare_folder(caller_info, pattern="_pcb_")

    settings = settings_info["pcb"].copy()
    vertical_scale_name = caller_info["vertical_scale"]

    pcb_pack = data_pack["pcb"]
    pcb_p45_pack = data_pack["pcb_p45"]
    pcb_m45_pack = data_pack["pcb_m45"]
    ray_pack = data_pack["ray_pcb"]

    pcb_ratio = pcb_pack["pol_cal_ratio"]
    pcb_info = pcb_pack["pol_cal_info"]

    pcb_p45_ratio = pcb_p45_pack["pol_cal_ratio"]
    pcb_m45_ratio = pcb_m45_pack["pol_cal_ratio"]

    ray_ratio = ray_pack["pol_cal_ratio"]
    ray_info = ray_pack["pol_cal_info"]

    molecular_ratio = ray_pack.get("molecular_ratio", None)
    molecular_info = ray_pack.get("molecular_info", None)

    # data_pack["pcb"] has only combined products, so use a source
    # calibration pack for the calibration vertical scale.
    z_cal_all = convert_m_to_km(pcb_m45_pack[vertical_scale_name])
    z_ray_all = convert_m_to_km(ray_pack[vertical_scale_name])

    eta_ids = _select_pairs_by_type(pcb_info, "eta")

    for eta_id in eta_ids:

        pair_settings = settings.copy()

        ch_r = _info_value(pcb_info, eta_id, "ch_r", default=None)
        ch_t = _info_value(pcb_info, eta_id, "ch_t", default=None)

        if ch_r is None or ch_t is None:
            continue

        print(f"-- channels: {ch_r} & {ch_t}")

        gain_ratio_id = replace_ratio_id_type(eta_id, "g")
        eta_s_f_id = replace_ratio_id_type(eta_id, "f")
        calibrated_ratio_id = replace_ratio_id_type(eta_id, "d")
        vldr_id = replace_ratio_id_type(eta_id, "v")
        mldr_id = replace_ratio_id_type(eta_id, "m")

        if not _all_pairs_exist(pcb_ratio, [gain_ratio_id, eta_s_f_id, eta_id]):
            print(f"   Skipping {eta_id}: missing pcb calibration products.")
            continue

        if not _all_pairs_exist(pcb_p45_ratio, [gain_ratio_id]):
            print(f"   Skipping {eta_id}: missing pcb_p45 gain-ratio product.")
            continue

        if not _all_pairs_exist(pcb_m45_ratio, [gain_ratio_id]):
            print(f"   Skipping {eta_id}: missing pcb_m45 gain-ratio product.")
            continue

        if not _all_pairs_exist(ray_ratio, [calibrated_ratio_id, vldr_id]):
            print(f"   Skipping {eta_id}: missing ray_pcb calibration products.")
            continue

        # ------------------------------------------------------------------
        # Calibration panel: gain-ratio profiles
        # ------------------------------------------------------------------
        z_cal = z_cal_all.sel(channel=ch_r)

        gain_ratio = pcb_ratio.sel(pair=gain_ratio_id)
        gain_ratio_p45 = pcb_p45_ratio.sel(pair=gain_ratio_id)
        gain_ratio_m45 = pcb_m45_ratio.sel(pair=gain_ratio_id)

        gain_ratio, z_cal_sliced, cal_mask = slice_by_vertical_scale(
            da=gain_ratio,
            vertical_scale=z_cal,
            x_lims=pair_settings["smoothing_range"],
        )

        gain_ratio_p45 = gain_ratio_p45.where(cal_mask, drop=True)
        gain_ratio_m45 = gain_ratio_m45.where(cal_mask, drop=True)

        x_cal = z_cal_sliced.values

        Y_cal = {}
        E_cal = {}

        for key, da in [
            ("gain_ratio", gain_ratio),
            ("gain_ratio_p45", gain_ratio_p45),
            ("gain_ratio_m45", gain_ratio_m45),
        ]:
            y_sm, y_err = smoothing(
                args=pair_settings,
                x_vals=x_cal,
                y_vals=da.values,
                err_type="std",
            )

            Y_cal[key] = y_sm
            E_cal[key] = y_err

        # ------------------------------------------------------------------
        # Rayleigh panel: time-resolved ray_pcb products -> mean over time
        # ------------------------------------------------------------------
        z_ray = z_ray_all.sel(channel=ch_r)

        calibrated_ratio = _time_mean(ray_ratio.sel(pair=calibrated_ratio_id))
        vldr = _time_mean(ray_ratio.sel(pair=vldr_id))

        calibrated_ratio, z_ray_sliced, ray_mask = slice_by_vertical_scale(
            da=calibrated_ratio,
            vertical_scale=z_ray,
            x_lims=pair_settings["smoothing_range"],
        )

        vldr = vldr.where(ray_mask, drop=True)

        x_ray = z_ray_sliced.values

        Y_ray = {}
        E_ray = {}

        for key, da in [
            ("calibrated_ratio", calibrated_ratio),
            ("vldr", vldr),
        ]:
            y_sm, y_err = smoothing(
                args=pair_settings,
                x_vals=x_ray,
                y_vals=da.values,
                err_type="std",
            )

            Y_ray[key] = y_sm
            E_ray[key] = y_err

        # Optional MLDR
        if molecular_ratio is not None and _all_pairs_exist(molecular_ratio, [mldr_id]):
            mldr = molecular_ratio.sel(pair=mldr_id)
            mldr = mldr.where(ray_mask, drop=True)

            mldr_sm, _ = smoothing(
                args=pair_settings,
                x_vals=x_ray,
                y_vals=mldr.values,
                err_type="std",
            )

            Y_ray["mldr"] = mldr_sm
            E_ray["mldr"] = np.nan * mldr_sm
        else:
            Y_ray["mldr"] = np.nan * Y_ray["vldr"]
            E_ray["mldr"] = np.nan * Y_ray["vldr"]

        # ------------------------------------------------------------------
        # Metadata / scalar diagnostics
        # ------------------------------------------------------------------
        pair_info = {
            "gain_ratio_id": gain_ratio_id,
            "eta_s_f_id": eta_s_f_id,
            "eta_id": eta_id,
            "calibrated_ratio_id": calibrated_ratio_id,
            "vldr_id": vldr_id,
            "mldr_id": mldr_id,
            "ch_r": ch_r,
            "ch_t": ch_t,
        }

        scalar_info = _collect_scalar_info(
            pcb_info=pcb_info,
            ray_info=ray_info,
            molecular_info=molecular_info,
            ids=pair_info,
            settings=pair_settings,
        )

        qa_test_info["pcb"][eta_id] = scalar_info | pair_info

        metadata_ray_r = collect_metadata(ray_pack, atlas_channel_id=ch_r)
        metadata_ray_t = collect_metadata(ray_pack, atlas_channel_id=ch_t)

        # Use source calibration packs for PCB measurement metadata.
        metadata_pcb_r = collect_metadata(pcb_m45_pack, atlas_channel_id=ch_r)
        metadata_pcb_t = collect_metadata(pcb_m45_pack, atlas_channel_id=ch_t)

        lib = Libraries(
            caller_info=caller_info,
            metadata=metadata_ray_r,
            extra_metadata=metadata_pcb_r,
            settings=pair_settings,
            qa_test_info=qa_test_info["pcb"][eta_id],
        )

        text_generator = GenerateText(lib=lib)

        # Use the polarization title function if available in your module.
        qa_test_info["pcb"][eta_id]["title"] = text_generator.make_polarization_calibration_title(
            metadata_r=metadata_ray_r,
            metadata_t=metadata_ray_t,
        )

        qa_test_info["pcb"][eta_id]["filename"] = text_generator.make_filename(
            qa_test="pcb",
            extra_metadata=metadata_ray_r,
        )

        print(scalar_info['epsilon'])
        plot_args = (
            caller_info
            | pair_settings
            | metadata_ray_r
            | qa_test_info["pcb"][eta_id]
            | {
                "vertical_scale": vertical_scale_name,
                "calibration_text": plot_polarization_calibration.make_calibration_text(
                    pair_settings | scalar_info
                ),
                "rayleigh_text": plot_polarization_calibration.make_rayleigh_text(
                    pair_settings | scalar_info
                ),
            }
        )

        plot_path = plot_polarization_calibration.generate_plot(
            X_cal=x_cal,
            Y_cal=Y_cal,
            E_cal=E_cal,
            X_ray=x_ray,
            Y_ray=Y_ray,
            E_ray=E_ray,
            args=plot_args,
        )

        qa_test_info["pcb"][eta_id]["pol_cal_plot_path"] = plot_path

        perform_color_reduction(
            color_reduction=caller_info["color_reduction"],
            plot_path=plot_path,
        )

        plot_metadata = (
            metadata_ray_r
            | {
                **pair_settings,
                **pair_info,
                **scalar_info,
                "atlas_channel_id_r": ch_r,
                "atlas_channel_id_t": ch_t,
                "ATLAS_version": __version__,
                "QA_test_ID": "pcb",
                "calibration_vertical_scale_source": "pcb_m45",
            }
        )

        add_plot_metadata(
            plot_path=plot_path,
            plot_metadata=plot_metadata,
            plot_metadata_extra=metadata_ray_t | metadata_pcb_r | metadata_pcb_t,
        )

    return qa_test_info


def replace_ratio_id_type(pair_id, new_type):
    """
    Replace the product-type character at index 5 of an 8-character pair ID.

    Examples
    --------
    0532xeax -> 0532xgax if new_type='g'
    0532xeax -> 0532xdax if new_type='d'
    """
    pair_id = str(pair_id)

    if len(pair_id) != 8:
        raise ValueError(f"Expected 8-character pair id. Got: {pair_id}")

    return f"{pair_id[:5]}{new_type}{pair_id[6:]}"


def _time_mean(da):
    """
    Average over time if a time dimension exists.
    """
    if "time" in da.dims:
        return da.mean("time")

    return da


def _select_pairs_by_type(info, ratio_type):
    """
    Select pair IDs from a pol_cal_info DataArray by ratio_type.
    """
    if info is None:
        return []

    if "parameters" not in info.dims:
        return []

    if "ratio_type" not in info.parameters.values:
        return []

    types = info.sel(parameters="ratio_type")

    return [
        str(pair)
        for pair, value in zip(types.pair.values, types.values)
        if str(value) == ratio_type
    ]


def _all_pairs_exist(da, pair_ids):
    """
    Check whether all requested pair IDs exist in a DataArray.
    """
    if da is None:
        return False

    if "pair" not in da.dims:
        return False

    available = set(str(pair) for pair in da.pair.values)

    return all(str(pair_id) in available for pair_id in pair_ids)


def _info_value(info, pair_id, parameter, default=np.nan):
    """
    Safely extract a scalar value from a pol_cal_info-like DataArray.
    """
    if info is None:
        return default

    if "parameters" not in info.dims or "pair" not in info.dims:
        return default

    if parameter not in info.parameters.values:
        return default

    if pair_id not in info.pair.values:
        return default

    value = info.sel(parameters=parameter, pair=pair_id).values

    try:
        value = np.asarray(value).item()
    except Exception:
        pass

    if value is None:
        return default

    try:
        if np.isnan(value):
            return default
    except Exception:
        pass

    return value


def _collect_scalar_info(pcb_info, ray_info, molecular_info, ids, settings):
    """
    Collect scalar diagnostics already calculated by the processor.
    """
    gain_ratio_id = ids["gain_ratio_id"]
    eta_s_f_id = ids["eta_s_f_id"]
    eta_id = ids["eta_id"]
    calibrated_ratio_id = ids["calibrated_ratio_id"]
    vldr_id = ids["vldr_id"]
    mldr_id = ids["mldr_id"]

    return {
        "gain_ratio_mean": _info_value(pcb_info, gain_ratio_id, "mean"),
        "gain_ratio_sem": _info_value(pcb_info, gain_ratio_id, "sem"),
        "eta_s_f_mean": _info_value(pcb_info, eta_s_f_id, "mean"),
        "eta_s_f_sem": _info_value(pcb_info, eta_s_f_id, "sem"),
        "eta_mean": _info_value(pcb_info, eta_id, "mean"),
        "eta_sem": _info_value(pcb_info, eta_id, "sem"),
        "epsilon": _info_value(pcb_info, eta_id, "epsilon"),
        "epsilon_sem": _info_value(pcb_info, eta_id, "epsilon_error"),
        "K": _info_value(pcb_info, eta_id, "K"),
        "calibrated_ratio_mean": _info_value(ray_info, calibrated_ratio_id, "mean"),
        "calibrated_ratio_sem": _info_value(ray_info, calibrated_ratio_id, "sem"),
        "vldr_mean": _info_value(ray_info, vldr_id, "mean"),
        "vldr_sem": _info_value(ray_info, vldr_id, "sem"),
        "vldr_residual": _info_value(ray_info, vldr_id, "vldr_residual"),
        "sr_limit": _info_value(ray_info, vldr_id, "sr_limit"),
        "G_R": _info_value(ray_info, vldr_id, "G_R"),
        "G_T": _info_value(ray_info, vldr_id, "G_T"),
        "H_R": _info_value(ray_info, vldr_id, "H_R"),
        "H_T": _info_value(ray_info, vldr_id, "H_T"),
        "mldr_mean": _info_value(molecular_info, mldr_id, "mean"),
        "pldr_error_threshold": settings.get("pldr_error_threshold", 0.025),
    }


def make_polarization_calibration_title(metadata_r, metadata_t):
    """
    Fallback title builder.

    You can replace this with a GenerateText.make_polarization_calibration_title()
    method later, following the Rayleigh-fit text architecture.
    """
    lidar_name = metadata_r.get("lidar_name", "")
    station_name = metadata_r.get("station_name", "")

    ch_r = metadata_r.get("atlas_channel_id", "")
    ch_t = metadata_t.get("atlas_channel_id", "")

    scc_r = metadata_r.get("scc_channel_id", "")
    scc_t = metadata_t.get("scc_channel_id", "")

    start = metadata_r.get("start_time_first", "")
    end = metadata_r.get("end_time_last", "")

    if hasattr(start, "strftime"):
        start_text = start.strftime("%d.%m.%Y %H:%M:%S")
    else:
        start_text = str(start)

    if hasattr(end, "strftime"):
        end_text = end.strftime("%H:%M:%S")
    else:
        end_text = str(end)

    return (
        f"{lidar_name} {station_name} {ch_r} ({scc_r}) to {ch_t} ({scc_t})\n"
        f"{start_text} to {end_text} UTC"
    ).strip()
