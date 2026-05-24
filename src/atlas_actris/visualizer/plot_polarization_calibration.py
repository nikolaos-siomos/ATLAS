#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 15:36:58 2026

@author: nikos
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from utils.toolbox import round_it
from visualizer.plot_utils import export_plot, get_vertical_axis_label, add_fitting_suptitle
from visualizer.generate_polarization_calibration_axis import (
    get_x_ticks,
    get_calibration_y_limits,
    get_rayleigh_y_limits,
)


def generate_plot(X_cal, Y_cal, E_cal, X_ray, Y_ray, E_ray, args):
    """
    Generate polarization calibration plot.

    Parameters
    ----------
    X_cal : np.ndarray
        Calibration vertical scale in km.
    Y_cal : dict
        Calibration profiles. Expected keys:
        gain_ratio, eta_s_f, eta
    E_cal : dict
        Calibration error profiles. Same keys as Y_cal, optional values may be NaN.
    X_ray : np.ndarray
        Rayleigh vertical scale in km.
    Y_ray : dict
        Rayleigh profiles. Expected keys:
        calibrated_ratio, vldr, mldr
    E_ray : dict
        Rayleigh error profiles. Same keys as Y_ray, optional values may be NaN.
    args : dict
        Plot/settings/metadata dictionary.
    """

    fig = plt.figure(figsize=(14.0, 3.5))

    add_fitting_suptitle(
        fig=fig,
        title=args["title"],
        y=0.98,
        max_fontsize=14,
        min_fontsize=10,
    )

    left_panel(
        fig=fig,
        coords=[0.045, 0.14, 0.44, 0.69],
        X=X_cal,
        Y=Y_cal,
        E=E_cal,
        args=args,
    )

    right_panel(
        fig=fig,
        coords=[0.545, 0.14, 0.44, 0.69],
        X=X_ray,
        Y=Y_ray,
        E=E_ray,
        args=args,
    )

    fpath = export_plot(fig, args)

    return fpath


def left_panel(fig, coords, X, Y, E, args):
    ax = fig.add_axes(coords)

    ax.plot(X, Y["gain_ratio"], color="tab:purple", label=r"$\eta^\star$")
    ax.plot(X, Y["gain_ratio_p45"], color="tab:red", label=r"$\eta^\star_{+45}$")
    ax.plot(X, Y["gain_ratio_m45"], color="tab:cyan", label=r"$\eta^\star_{-45}$")

    _fill_error(ax, X, Y["gain_ratio"], E.get("gain_ratio"), "tab:purple")
    # _fill_error(ax, X, Y["gain_ratio_p45"], E.get("gain_ratio_p45"), "tab:red")
    # _fill_error(ax, X, Y["gain_ratio_m45"], E.get("gain_ratio_m45"), "tab:cyan")

    x_ticks, x_labels, x_tick = get_x_ticks(
        x_lims=args["x_lims_calibration"],
        x_tick=args.get("x_tick_calibration"),
    )

    ax.set_xticks(x_ticks, labels=x_labels)
    ax.set_xlim(args["x_lims_calibration"])

    ax.set_xlabel(get_vertical_axis_label(args["vertical_scale"]))
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.0))

    y_lims = get_calibration_y_limits(
        y_values=[
            args.get("gain_ratio_m45_mean"),
            args.get("gain_ratio_p45_mean"),
            args.get("gain_ratio_mean"),
        ],
        y_lims=args.get("y_lims_calibration", []),
    )

    ax.set_ylim(y_lims)
    ax.set_ylabel(r"Gain ratio $\eta^{\star}_{f}$")

    ax.grid(which="both")

    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc="upper right")

    ax.axvspan(
        args["calibration_region"][0],
        args["calibration_region"][1],
        alpha=0.2,
        facecolor="tab:grey",
    )

    add_calibration_text(ax, args)

    return ax


def right_panel(fig, coords, X, Y, E, args):
    ax = fig.add_axes(coords)

    ax.plot(X, Y["calibrated_ratio"], color="tab:blue", label=r"$\delta^\star$")
    ax.plot(X, Y["vldr"], color="tab:orange", label=r"$\delta_v$")
    ax.plot(X, Y["mldr"], color="tab:green", label=r"$\delta_m$")

    # _fill_error(ax, X, Y["calibrated_ratio"], E.get("calibrated_ratio"), "tab:blue")
    _fill_error(ax, X, Y["vldr"], E.get("vldr"), "tab:orange")

    x_ticks, x_labels, x_tick = get_x_ticks(
        x_lims=args["x_lims_rayleigh"],
        x_tick=args.get("x_tick_rayleigh"),
    )

    ax.set_xticks(x_ticks, labels=x_labels)
    ax.set_xlim(args["x_lims_rayleigh"])

    ax.set_xlabel(get_vertical_axis_label(args["vertical_scale"]))
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 2.0))

    y_lims = get_rayleigh_y_limits(
        y_values=[
            args.get("calibrated_ratio_mean"),
            args.get("vldr_mean"),
            args.get("mldr_mean"),
        ],
        y_lims=args.get("y_lims_rayleigh", []),
    )

    ax.set_ylim(y_lims)
    ax.set_ylabel("Linear dep. ratio")

    ax.grid(which="both")

    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc="upper right")

    ax.axvspan(
        args["rayleigh_region"][0],
        args["rayleigh_region"][1],
        alpha=0.2,
        facecolor="tab:grey",
    )

    add_rayleigh_text(ax, args)

    return ax


def _fill_error(ax, X, Y, E, color):
    if E is None:
        return

    E = np.asarray(E, dtype=float)

    if E.size == 0 or np.isnan(E).all():
        return

    ax.fill_between(X, Y - E, Y + E, color=color, alpha=0.3)


def add_calibration_text(ax, args):
    text = args["calibration_text"]

    rows = [
        ("region", 0.90, "tab:green"),
        ("epsilon", 0.78, "tab:cyan"),
        ("gain_ratio", 0.66, "tab:cyan"),
        ("eta_s_f", 0.54, "tab:cyan"),
        ("eta", 0.42, text.get("eta_color", "tab:cyan")),
    ]

    for key, y, color in rows:
        if key not in text:
            continue

        ax.text(
            0.04,
            y,
            text[key],
            transform=ax.transAxes,
            bbox=dict(facecolor=color, alpha=0.22, zorder=3),
        )


def add_rayleigh_text(ax, args):
    text = args["rayleigh_text"]

    rows = [
        ("region", 0.90, "tab:green"),
        ("mldr", 0.78, "tab:cyan"),
        ("calibrated_ratio", 0.66, "tab:cyan"),
        ("vldr", 0.54, "tab:cyan"),
        ("residual", 0.42, "tab:cyan"),
        ("G", 0.30, "tab:cyan"),
        ("H", 0.18, "tab:cyan"),
        ("sr_limit", 0.06, "tab:cyan"),
    ]

    for key, y, color in rows:
        if key not in text:
            continue

        ax.text(
            0.04,
            y,
            text[key],
            transform=ax.transAxes,
            bbox=dict(facecolor=color, alpha=0.22, zorder=3),
        )


def make_calibration_text(args):
    c0, c1 = args["calibration_region"]

    eta = args.get("eta_mean", np.nan)
    eta_sem = args.get("eta_sem", np.nan)
    eta_s_f = args.get("eta_s_f_mean", np.nan)
    eta_s_f_sem = args.get("eta_s_f_sem", np.nan)
    gain_ratio = args.get("gain_ratio_mean", np.nan)
    gain_ratio_sem = args.get("gain_ratio_sem", np.nan)
    epsilon = args.get("epsilon", np.nan)

    eta_color = "tab:green"
    if np.isfinite(eta) and np.isfinite(eta_sem) and eta != 0:
        if abs(eta_sem / eta) > 0.02:
            eta_color = "tab:red"

    return {
        "region": f"cal. region: {np.round(c0, 2)} - {np.round(c1, 2)} km",
        "epsilon": r"$\epsilon$: "
        + f"{round_it(epsilon, 2)}"
        + r"$^{o}$",
        "gain_ratio": r"$\eta^\star$: "
        + f"{round_it(gain_ratio, 3)}"
        + r" $\pm$ "
        + f"{round_it(gain_ratio_sem, 2)}",
        "eta_s_f": r"$\eta^\star_f$: "
        + f"{round_it(eta_s_f, 3)}"
        + r" $\pm$ "
        + f"{round_it(eta_s_f_sem, 2)}",
        "eta": r"$\eta$: "
        + f"{round_it(eta, 3)}"
        + r" $\pm$ "
        + f"{round_it(eta_sem, 2)}",
        "eta_color": eta_color,
    }


def make_rayleigh_text(args):
    r0, r1 = args["rayleigh_region"]

    return {
        "region": f"mol. cal. region: {np.round(r0, 2)} - {np.round(r1, 2)} km",
        "mldr": r"$\delta_m$: " + f"{round_it(args.get('mldr_mean', np.nan), 4)}",
        "calibrated_ratio": r"$\delta^\star$: "
        + f"{round_it(args.get('calibrated_ratio_mean', np.nan), 4)}"
        + r" $\pm$ "
        + f"{round_it(args.get('calibrated_ratio_sem', np.nan), 2)}",
        "vldr": r"$\delta_v$: "
        + f"{round_it(args.get('vldr_mean', np.nan), 4)}"
        + r" $\pm$ "
        + f"{round_it(args.get('vldr_sem', np.nan), 2)}",
        "residual": r"$\delta_{res}$: "
        + f"{round_it(args.get('vldr_residual', np.nan), 4)}",
        "G": r"$G_R$: "
        + f"{round_it(args.get('G_R', np.nan), 4)}"
        + r", $G_T$: "
        + f"{round_it(args.get('G_T', np.nan), 4)}",
        "H": r"$H_R$: "
        + f"{round_it(args.get('H_R', np.nan), 4)}"
        + r", $H_T$: "
        + f"{round_it(args.get('H_T', np.nan), 4)}",
        "sr_limit": r"$SR$ > "
        + f"{round_it(args.get('sr_limit', np.nan), 3)}"
        + r", $\Delta\delta_p$ < "
        + f"{round_it(args.get('pldr_error_threshold', 0.025), 3)}",
    }