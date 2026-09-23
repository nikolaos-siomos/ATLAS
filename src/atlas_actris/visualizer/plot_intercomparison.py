#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-panel ATLAS intercomparison plot renderer.

The left panel shows all participating entries and, optionally, the reference
molecular profile. The right panel shows relative channel differences or
absolute pair differences to the reference when vertical grids are aligned.
When native vertical
scales are requested, the right panel is intentionally left without data.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from visualizer.plot_utils import export_plot


# Reserve tab:green for molecular products. Entry colours intentionally skip it.
ENTRY_COLORS = (
    "tab:blue",
    "tab:orange",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
    "tab:cyan",
)
MOLECULAR_COLOR = "tab:green"
ENTRY_LINESTYLES = ("-", "--", "-.", ":")


def _vertical_axis_label(vertical_scale):
    labels = {
        "bins": "Bins",
        "range": "Range from the lidar [km]",
        "height_agl": "Height [km agl]",
        "height_asl": "Height [km asl]",
    }
    if vertical_scale not in labels:
        raise ValueError(
            "Unsupported vertical_scale {!r}. Expected one of {}".format(
                vertical_scale, tuple(labels)
            )
        )
    return labels[vertical_scale]


def _finite_error(error):
    if error is None:
        return False
    arr = np.asarray(error, dtype=float)
    return np.isfinite(arr).any()


def _apply_x_axis(ax, args, *, minor_divisor=2.0):
    x_lims = args["x_lims"]
    x_tick = float(args["x_tick"])

    ax.set_xlim(x_lims)
    ax.set_xlabel(_vertical_axis_label(args["vertical_scale"]))

    if x_tick > 0:
        ax.xaxis.set_major_locator(MultipleLocator(x_tick))
        ax.xaxis.set_minor_locator(MultipleLocator(x_tick / minor_divisor))


def _shade_normalisation_region(ax, args):
    region = args.get("normalisation_region")
    if region is None or len(region) != 2:
        return
    ax.axvspan(region[0], region[1], alpha=0.15)


def left_panel(fig, ax_coords, X, Y, YE, molecular, args):
    ax = fig.add_axes(ax_coords)

    entry_colors = {}
    entry_styles = {}
    for index, (entry_id, y) in enumerate(Y.items()):
        x = np.asarray(X[entry_id], dtype=float)
        y = np.asarray(y, dtype=float)
        label = args.get("entry_labels", {}).get(entry_id, entry_id)
        color_index = index % len(ENTRY_COLORS)
        style_index = (index // len(ENTRY_COLORS)) % len(ENTRY_LINESTYLES)
        color = ENTRY_COLORS[color_index]
        linestyle = ENTRY_LINESTYLES[style_index]
        entry_colors[entry_id] = color
        entry_styles[entry_id] = linestyle

        line, = ax.plot(
            x, y, label=label, color=color, linestyle=linestyle
        )
        error = YE.get(entry_id)
        if _finite_error(error):
            error = np.asarray(error, dtype=float)
            ax.fill_between(
                x,
                y - error,
                y + error,
                alpha=0.18,
                color=line.get_color(),
            )

    if molecular is not None:
        ax.plot(
            np.asarray(molecular["x"], dtype=float),
            np.asarray(molecular["y"], dtype=float),
            linestyle="-",
            linewidth=1.6,
            color=MOLECULAR_COLOR,
            label=molecular.get("label", "molecular"),
        )

    _apply_x_axis(ax, args)
    ax.set_ylim(args["y_lims"])
    ax.set_ylabel(args.get("left_y_label", "Signal"))

    if args.get("use_log_y_scale", False):
        ax.set_yscale("log")

    _shade_normalisation_region(ax, args)
    ax.grid(which="both")

    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc="best", fontsize=8)

    args["entry_colors"] = entry_colors
    args["entry_styles"] = entry_styles
    return ax


def right_panel(fig, ax_coords, X, differences, difference_error, args):
    ax = fig.add_axes(ax_coords)
    _apply_x_axis(ax, args)
    if args.get("difference_mode") == "absolute":
        ax.set_ylabel("Absolute Diff. to Reference")
    else:
        ax.set_ylabel("Relative Diff. to Reference")
    ax.axhline(0.0, linewidth=1.0)
    _shade_normalisation_region(ax, args)

    if args.get("plot_native_scale", False):
        # Keep the Rayleigh-like two-panel layout, but deliberately do not plot
        # differences because native vertical grids are not aligned.
        ax.set_ylim(args["difference_lims"])
        ax.grid(which="both")
        return ax

    colors = args.get("entry_colors", {})
    styles = args.get("entry_styles", {})
    for index, (entry_id, diff) in enumerate(differences.items()):
        x = np.asarray(X[entry_id], dtype=float)
        diff = np.asarray(diff, dtype=float)
        label = args.get("entry_labels", {}).get(entry_id, entry_id)
        color = colors.get(entry_id, ENTRY_COLORS[index % len(ENTRY_COLORS)])
        default_style = ENTRY_LINESTYLES[
            (index // len(ENTRY_COLORS)) % len(ENTRY_LINESTYLES)
        ]
        linestyle = styles.get(entry_id, default_style)

        line, = ax.plot(
            x, diff, label=label, color=color, linestyle=linestyle
        )
        error = difference_error.get(entry_id)
        if _finite_error(error):
            error = np.asarray(error, dtype=float)
            ax.fill_between(
                x,
                diff - error,
                diff + error,
                alpha=0.18,
                color=line.get_color(),
            )

    ax.set_ylim(args["difference_lims"])
    ax.grid(which="both")

    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc="best", fontsize=8)

    return ax


def generate_plot(X, Y, YE, differences, difference_error, molecular, args):
    """Generate and export one two-panel intercomparison figure."""

    ax1_coords = [0.055, 0.17, 0.52, 0.66]
    ax2_coords = [0.625, 0.17, 0.35, 0.66]

    fig = plt.figure(figsize=(15, 3.4))
    fig.suptitle(args["title"], fontsize=11)

    left_panel(
        fig=fig,
        ax_coords=ax1_coords,
        X=X,
        Y=Y,
        YE=YE,
        molecular=molecular,
        args=args,
    )
    right_panel(
        fig=fig,
        ax_coords=ax2_coords,
        X=X,
        differences=differences,
        difference_error=difference_error,
        args=args,
    )

    return export_plot(fig, args)
