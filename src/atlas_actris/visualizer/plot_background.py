#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Background time-series plotting utilities."""

import numpy as np
import matplotlib.dates as mdates
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from visualizer.plot_utils import export_plot, add_fitting_suptitle


def generate_plot(T, Y, Y_E, args):
    """Generate a background line plot against time."""

    T, Y, Y_E = prepare_background_data(T, Y, Y_E)

    fig_coords = [0.09, 0.15, 0.88, 0.72]
    fig = plt.figure(figsize=(10.0, 5.0))

    add_fitting_suptitle(
        fig,
        args["title"],
        y=0.99,
        max_fontsize=12,
        min_fontsize=7,
    )

    background_panel(
        fig=fig,
        fig_coords=fig_coords,
        T=T,
        Y=Y,
        Y_E=Y_E,
        args=args,
    )

    return export_plot(fig, args)


def prepare_background_data(T, Y, Y_E):
    """Validate and flatten one-dimensional time-series arrays."""

    T = np.asarray(T)
    Y = np.asarray(Y).squeeze()
    Y_E = None if Y_E is None else np.asarray(Y_E).squeeze()

    if T.ndim != 1:
        raise ValueError(f"T must be one-dimensional. Found shape {T.shape}.")

    if Y.ndim != 1:
        raise ValueError(
            "Background values must be one-dimensional after selecting one "
            f"channel. Found shape {Y.shape}."
        )

    if len(T) != len(Y):
        raise ValueError(
            f"T and Y must have the same length. Found {len(T)} and {len(Y)}."
        )

    if Y_E is not None:
        if Y_E.ndim != 1:
            raise ValueError(
                "Background-error values must be one-dimensional after "
                f"selecting one channel. Found shape {Y_E.shape}."
            )
        if len(Y_E) != len(Y):
            raise ValueError(
                "Y and Y_E must have the same length. "
                f"Found {len(Y)} and {len(Y_E)}."
            )

    return T, Y, Y_E


def background_panel(fig, fig_coords, T, Y, Y_E, args):
    ax = fig.add_axes(fig_coords)

    ax.set_xlabel("Time UTC")
    
    if args['atlas_channel_id'][6] == 'a':
        y_units = 'mV'
        
        if args['input_qa_test'].startswith('drk'):
            bg_type = 'Background'
        else:
            bg_type = 'Dark corr. background'
            
    else:
        y_units = 'MHz'
        bg_type = 'Background'
        
    ax.set_ylabel(f"{bg_type} [{y_units}]")

    add_background_line(ax=ax, T=T, Y=Y, Y_E=Y_E)

    set_time_axis(
        ax=ax,
        t_tick=args.get("t_tick"),
    )

    set_y_axis(
        ax=ax,
        y_lims=args.get("y_lims", []),
    )

    add_minor_ticks(ax=ax)
    ax.grid(which="major", alpha=0.35)
    ax.grid(which="minor", alpha=0.15)

    return ax


def add_background_line(ax, T, Y, Y_E=None):
    """Plot the background signal and, when available, its uncertainty."""

    line, = ax.plot(T, Y, linewidth=1.0)

    if Y_E is not None:
        valid = np.isfinite(Y) & np.isfinite(Y_E)
        if np.any(valid):
            ax.fill_between(
                T,
                Y - Y_E,
                Y + Y_E,
                where=valid,
                color=line.get_color(),
                alpha=0.3,
                linewidth=0.0,
            )

    return ax


def set_y_axis(ax, y_lims):
    """Apply optional background-axis limits from the settings file."""

    if y_lims is not None and len(y_lims) == 2:
        ax.set_ylim(y_lims)

    try:
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    except Exception:
        pass

    return ax


def set_time_axis(ax, t_tick=None, target_ticks=8):
    """Set automatic or user-defined time tick spacing."""

    if t_tick is None:
        major_locator = mdates.AutoDateLocator(
            minticks=max(3, target_ticks - 2),
            maxticks=target_ticks + 2,
        )

        major_locator.intervald[mdates.SECONDLY] = [1, 2, 5, 10, 15, 30]
        major_locator.intervald[mdates.MINUTELY] = [1, 2, 5, 10, 15, 20, 30]
        major_locator.intervald[mdates.HOURLY] = [1, 2, 3, 4, 6, 12]
        major_locator.intervald[mdates.DAILY] = [1, 2, 3, 7, 14]
        major_locator.intervald[mdates.MONTHLY] = [1, 2, 3, 6]
        major_locator.intervald[mdates.YEARLY] = [1, 2, 5, 10]

        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(major_locator))

        try:
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        except Exception:
            pass
    else:
        major_interval_sec = max(1, int(round(60 * t_tick)))
        minor_interval_sec = max(1, int(round(15 * t_tick)))

        major_locator = mdates.SecondLocator(interval=major_interval_sec)
        major_locator.MAXTICKS = 10000
        minor_locator = mdates.SecondLocator(interval=minor_interval_sec)
        minor_locator.MAXTICKS = 10000

        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_minor_locator(minor_locator)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")
    return ax


def add_minor_ticks(ax):
    """Enable and style minor ticks on both axes."""

    try:
        ax.minorticks_on()
    except Exception:
        pass

    ax.tick_params(axis="both", which="minor", length=2.5, width=0.6)
    ax.tick_params(axis="both", which="major", length=4.0, width=0.8)
    return ax
