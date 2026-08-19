#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-panel ATLAS intercomparison plot renderer.

The left panel shows all participating datasets and, optionally, the reference
molecular profile.  The right panel shows relative differences to the
reference dataset when the vertical grids are aligned.  When native vertical
scales are requested, the right panel is intentionally left without data.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from visualizer.plot_utils import export_plot


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

    for dataset_id, y in Y.items():
        x = np.asarray(X[dataset_id], dtype=float)
        y = np.asarray(y, dtype=float)
        label = args.get("dataset_labels", {}).get(dataset_id, dataset_id)

        line, = ax.plot(x, y, label=label)
        error = YE.get(dataset_id)
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
            linestyle="--",
            linewidth=1.5,
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

    return ax


def right_panel(fig, ax_coords, X, relative, relative_error, args):
    ax = fig.add_axes(ax_coords)
    _apply_x_axis(ax, args)
    ax.set_ylabel("Relative Diff. to Reference")
    ax.axhline(0.0, linewidth=1.0)
    _shade_normalisation_region(ax, args)

    if args.get("plot_native_scale", False):
        # Keep the Rayleigh-like two-panel layout, but deliberately do not plot
        # differences because native vertical grids are not aligned.
        ax.set_ylim(args["relative_difference_lims"])
        ax.grid(which="both")
        return ax

    for dataset_id, rel in relative.items():
        x = np.asarray(X[dataset_id], dtype=float)
        rel = np.asarray(rel, dtype=float)
        label = args.get("dataset_labels", {}).get(dataset_id, dataset_id)

        line, = ax.plot(x, rel, label=label)
        error = relative_error.get(dataset_id)
        if _finite_error(error):
            error = np.asarray(error, dtype=float)
            ax.fill_between(
                x,
                rel - error,
                rel + error,
                alpha=0.18,
                color=line.get_color(),
            )

    ax.set_ylim(args["relative_difference_lims"])
    ax.grid(which="both")

    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc="best", fontsize=8)

    return ax


def generate_plot(X, Y, YE, relative, relative_error, molecular, args):
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
        relative=relative,
        relative_error=relative_error,
        args=args,
    )

    return export_plot(fig, args)
