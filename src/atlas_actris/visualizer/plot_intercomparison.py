#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Four-panel ATLAS intercomparison plot renderer.

The upper row shows the configured full-range comparison. The lower row shows
the same entries over the configurable near range with its own smoothing.
Molecular products always use ``tab:green``; measured entries use a deliberately
selected high-contrast categorical palette that contains no green.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from visualizer.plot_utils import export_plot


MOLECULAR_COLOR = "tab:green"
AXIS_LABEL_FONTSIZE = 10.5
TICK_LABEL_FONTSIZE = 9.5
PANEL_TITLE_FONTSIZE = 11
FIGURE_TITLE_FONTSIZE = 12
LEGEND_FONTSIZE = 9
ENTRY_LINESTYLES = ("-", "--", "-.", ":")

# High-contrast measured-entry palette.  The molecular curve keeps its dedicated
# ``tab:green`` color; measured entries may also use deliberately distinct olive
# and light-green tones.  The first colours are chosen to be very distinct on a
# white background; line styles provide an additional cue if a group contains
# more entries than colours.
ENTRY_COLORS = (
    "#000000",  # black
    "#0072B2",  # strong blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#7A7A00",  # olive
    "#E69F00",  # orange
    "#6E6E6E",  # medium grey
    "#56B4E9",  # sky blue
    "#6F42C1",  # violet
    "#009E73",  # teal
    "#8C564B",  # brown
    "#1F3A93",  # navy
)


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
    ax.set_xlabel(
        _vertical_axis_label(args["vertical_scale"]),
        labelpad=2,
        fontsize=AXIS_LABEL_FONTSIZE,
    )
    ax.tick_params(axis="x", which="major", pad=2, labelsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="y", which="major", pad=2, labelsize=TICK_LABEL_FONTSIZE)

    if x_tick > 0:
        ax.xaxis.set_major_locator(MultipleLocator(x_tick))
        ax.xaxis.set_minor_locator(MultipleLocator(x_tick / minor_divisor))


def _shade_normalisation_region(ax, args):
    region = args.get("normalisation_region")
    if region is None or len(region) != 2:
        return
    ax.axvspan(region[0], region[1], alpha=0.15)


def _entry_style(index):
    color_index = index % len(ENTRY_COLORS)
    style_index = (index // len(ENTRY_COLORS)) % len(ENTRY_LINESTYLES)
    return ENTRY_COLORS[color_index], ENTRY_LINESTYLES[style_index]


def _entry_style_from_color_index(color_index, *, entry_id):
    """Return a pinned palette color for a 1-based external color index."""
    if color_index is None:
        return None
    try:
        index = int(color_index)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Entry {entry_id!r} color_index must be an integer, got {color_index!r}"
        ) from exc
    if index < 1 or index > len(ENTRY_COLORS):
        raise ValueError(
            f"Entry {entry_id!r} color_index={index} is outside the available "
            f"palette range 1..{len(ENTRY_COLORS)}"
        )
    return ENTRY_COLORS[index - 1]


def left_panel(ax, X, Y, YE, molecular, args, *, show_legend=True, title=None):

    entry_colors = args.get("entry_colors", {})
    entry_styles = args.get("entry_styles", {})
    for index, (entry_id, y) in enumerate(Y.items()):
        x = np.asarray(X[entry_id], dtype=float)
        y = np.asarray(y, dtype=float)
        label = args.get("entry_labels", {}).get(entry_id, entry_id)
        default_color, default_style = _entry_style(index)
        pinned_color = _entry_style_from_color_index(
            args.get("entry_color_indices", {}).get(entry_id), entry_id=entry_id
        )
        color = entry_colors.get(entry_id, pinned_color or default_color)
        linestyle = entry_styles.get(entry_id, default_style)
        entry_colors[entry_id] = color
        entry_styles[entry_id] = linestyle

        line, = ax.plot(x, y, label=label, color=color, linestyle=linestyle)
        error = YE.get(entry_id)
        if _finite_error(error):
            error = np.asarray(error, dtype=float)
            ax.fill_between(
                x, y - error, y + error, alpha=0.18, color=line.get_color()
            )

    if molecular is not None:
        ax.plot(
            np.asarray(molecular["x"], dtype=float),
            np.asarray(molecular["y"], dtype=float),
            linestyle="-", linewidth=1.6, color=MOLECULAR_COLOR,
            label=molecular.get("label", "molecular"),
        )

    _apply_x_axis(ax, args)
    ax.set_ylim(args["y_lims"])
    ax.set_ylabel(
        args.get("left_y_label", "Signal"),
        labelpad=3,
        fontsize=AXIS_LABEL_FONTSIZE,
    )
    if title:
        ax.set_title(title, fontsize=PANEL_TITLE_FONTSIZE, pad=4)

    if args.get("use_log_y_scale", False):
        ax.set_yscale("log")

    _shade_normalisation_region(ax, args)
    ax.grid(which="both")

    if show_legend and ax.get_legend_handles_labels() != ([], []):
        ax.legend(
            loc="best",
            fontsize=LEGEND_FONTSIZE,
            borderaxespad=0.35,
            borderpad=0.3,
            handlelength=2.0,
            handletextpad=0.5,
            labelspacing=0.3,
        )

    args["entry_colors"] = entry_colors
    args["entry_styles"] = entry_styles
    return ax


def right_panel(
    ax, X, differences, difference_error, args, *, show_legend=False, title=None
):
    _apply_x_axis(ax, args)
    if args.get("difference_mode") == "absolute":
        ax.set_ylabel(
            "Absolute Diff. to Reference",
            labelpad=3,
            fontsize=AXIS_LABEL_FONTSIZE,
        )
    else:
        ax.set_ylabel(
            "Relative Diff. to Reference",
            labelpad=3,
            fontsize=AXIS_LABEL_FONTSIZE,
        )
    if title:
        ax.set_title(title, fontsize=PANEL_TITLE_FONTSIZE, pad=4)
    ax.axhline(0.0, linewidth=1.0)
    _shade_normalisation_region(ax, args)

    if args.get("plot_native_scale", False):
        ax.set_ylim(args["difference_lims"])
        ax.grid(which="both")
        return ax

    colors = args.get("entry_colors", {})
    styles = args.get("entry_styles", {})
    for index, (entry_id, diff) in enumerate(differences.items()):
        x = np.asarray(X[entry_id], dtype=float)
        diff = np.asarray(diff, dtype=float)
        label = args.get("entry_labels", {}).get(entry_id, entry_id)
        default_color, default_style = _entry_style(index)
        pinned_color = _entry_style_from_color_index(
            args.get("entry_color_indices", {}).get(entry_id), entry_id=entry_id
        )
        color = colors.get(entry_id, pinned_color or default_color)
        linestyle = styles.get(entry_id, default_style)

        line, = ax.plot(x, diff, label=label, color=color, linestyle=linestyle)
        error = difference_error.get(entry_id)
        if _finite_error(error):
            error = np.asarray(error, dtype=float)
            ax.fill_between(
                x, diff - error, diff + error, alpha=0.18, color=line.get_color()
            )

    ax.set_ylim(args["difference_lims"])
    ax.grid(which="both")

    if show_legend and ax.get_legend_handles_labels() != ([], []):
        ax.legend(
            loc="best",
            fontsize=LEGEND_FONTSIZE,
            borderaxespad=0.35,
            borderpad=0.3,
            handlelength=2.0,
            handletextpad=0.5,
            labelspacing=0.3,
        )

    return ax


def generate_plot(
    X, Y, YE, differences, difference_error, molecular, args,
    near_X, near_Y, near_YE, near_differences, near_difference_error,
    near_molecular, near_args,
):
    """Generate and export one compact 2x2 intercomparison figure."""

    # Use a real subplot grid instead of four manually positioned axes.  This
    # allows the labels and titles to share the available canvas efficiently
    # and avoids the large white margins produced by the previous coordinates.
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(15, 6.6),
        gridspec_kw={
            "width_ratios": (1.48, 1.0),
            "height_ratios": (1.0, 1.0),
        },
    )
    ax1, ax2 = axes[0]
    ax3, ax4 = axes[1]

    # Keep the overall figure title compact; individual axes carry their own
    # descriptive titles below it.
    fig.suptitle(args["title"], fontsize=FIGURE_TITLE_FONTSIZE, y=0.985)

    signal_name = args.get("left_y_label", "Signal")
    difference_name = (
        "Absolute difference"
        if args.get("difference_mode") == "absolute"
        else "Relative difference"
    )

    left_panel(
        ax1, X, Y, YE, molecular, args,
        show_legend=True,
        title=f"Full range - {signal_name}",
    )
    right_panel(
        ax2, X, differences, difference_error, args,
        show_legend=False,
        title=f"Full range - {difference_name}",
    )

    near_args["entry_colors"] = args.get("entry_colors", {})
    near_args["entry_styles"] = args.get("entry_styles", {})
    left_panel(
        ax3, near_X, near_Y, near_YE, near_molecular, near_args,
        show_legend=False,
        title=f"Near range - {signal_name}",
    )
    right_panel(
        ax4, near_X, near_differences, near_difference_error,
        near_args,
        show_legend=False,
        title=f"Near range - {difference_name}",
    )

    # Tight but readable spacing.  Labels stay close to their axes, while the
    # larger left column remains available for the many measured curves.
    fig.subplots_adjust(
        left=0.055,
        right=0.992,
        bottom=0.075,
        top=0.89,
        wspace=0.18,
        hspace=0.29,
    )

    return export_plot(fig, args)
