#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 00:22:21 2026

@author: nikos
"""

import os
import numpy as np
import matplotlib.dates as mdates
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

from visualizer.plot_utils import export_plot
from visualizer import color_lib, make_colormap
from visualizer.plot_utils import get_vertical_axis_label, add_fitting_suptitle

from visualizer.generate_quicklook_axis import (
    get_quicklook_y_limits, normalize_quicklook_y, max_quicklook_y
    )

def generate_plot(T, X, Y, args):

    # T, X, Y = prepare_quicklook_data(T, X, Y, args)
    
    fig_coords = [0.07, 0.13, 0.99, 0.74]

    fig = plt.figure(figsize=(10., 5.))
    
    add_fitting_suptitle(
        fig,
        args["title"],
        y=0.99,
        max_fontsize=12,
        min_fontsize=7,
    )

    quicklook_panel(
        fig=fig,
        fig_coords=fig_coords,
        T=T,
        X=X,
        Y=Y,
        args=args,
    )

    fpath = export_plot(fig, args)

    return fpath


def prepare_quicklook_data(T, X, Y, args):
    """
    Prepare time, height/range, and signal arrays for pcolormesh.

    Input Y is expected as:
        Y.shape == (len(T), len(X))
    """

    Y = normalize_quicklook_y(
        y_vals=Y,
        x_vals=X,
        y_max_zone=args["y_max_zone"],
    )

    return T, X, Y


def quicklook_panel(fig, fig_coords, T, X, Y, args):

    ax = fig.add_axes(fig_coords)
    
    # Get the x axis label depending on use_dis
    x_label = get_vertical_axis_label(args['vertical_scale'])

    ax.set_xlabel("Time UTC")
    ax.set_ylabel(x_label)

    set_x_axis(
        ax=ax,
        x_lims=args["x_lims"],
        x_tick=args["x_tick"],
    )

    qck = add_quicklook_mesh(
        ax=ax,
        T=T,
        X=X,
        Y=Y,
        args=args,
        )

    set_time_axis(
        ax=ax,
        t_tick=args["t_tick"],
    )
    
    if not args['has_time_gap']:
        add_time_index_axis(
            ax=ax,
            T=T,
        )

    qck = add_quicklook_mesh(
        ax=ax,
        T=T,
        X=X,
        Y=Y,
        args=args
    )

    add_colorbar(
        fig=fig,
        qck=qck,
    )

    return ax

def set_x_axis(ax, x_lims, x_tick):
    
    x_ticks = get_x_ticks(
        x_lims=x_lims,
        x_tick=x_tick,
    )

    ax.set_yticks(x_ticks, labels=x_ticks)
    ax.set_ylim(x_lims)

    return ax


def get_x_ticks(x_lims, x_tick):

    x_llim = x_lims[0]
    x_ulim = x_lims[1]
    
    x_ticks = np.arange(
        x_tick * np.ceil(x_llim / x_tick),
        x_tick * (np.floor(x_ulim / x_tick) + 1.),
        x_tick,
    )

    if np.abs(x_llim - x_ticks[0]) < x_tick * 0.25:
        x_ticks[0] = x_llim
    else:
        x_ticks = np.hstack((x_llim, x_ticks))

    if np.abs(x_ulim - x_ticks[-1]) < x_tick * 0.25:
        x_ticks[-1] = x_ulim
    else:
        x_ticks = np.hstack((x_ticks, x_ulim))

    x_ticks = np.round(x_ticks, decimals=2)

    return x_ticks


def set_time_axis(ax, t_tick=None, target_ticks=8):
    """
    Set time axis formatting.

    Parameters
    ----------
    ax : matplotlib axis
        Axis to modify.
    t_tick : float or None
        Major tick spacing in minutes.
        If None, a nice automatic tick spacing is selected.
    target_ticks : int
        Approximate number of major ticks when t_tick is None.

    Returns
    -------
    ax : matplotlib axis
        Modified axis.
    """

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

        formatter = mdates.ConciseDateFormatter(major_locator)

        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_major_formatter(formatter)

    else:
        # Manual tick spacing in minutes
        major_interval_sec = int(round(60 * t_tick))
        minor_interval_sec = int(round(15 * t_tick))

        # Avoid invalid interval=0 for very small t_tick values
        major_interval_sec = max(1, major_interval_sec)
        minor_interval_sec = max(1, minor_interval_sec)

        major_locator = mdates.SecondLocator(interval=major_interval_sec)
        major_locator.MAXTICKS = 10000

        minor_locator = mdates.SecondLocator(interval=minor_interval_sec)
        minor_locator.MAXTICKS = 10000

        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_minor_locator(minor_locator)

        ax.xaxis.set_major_formatter(
            mdates.DateFormatter("%H:%M")
        )

    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")

    return ax

def add_time_index_axis(ax, T, n_ticks=8):
    """
    Add a top twin x-axis showing nice index numbers of datetime array T.
    """

    ax_top = ax.twiny()

    # Important: copy limits after the main plot has set them
    ax_top.set_xlim(ax.get_xlim())

    locator = MaxNLocator(nbins=n_ticks, integer=True)
    idx_ticks = locator.tick_values(0, len(T) - 1)

    idx_ticks = idx_ticks.astype(int)
    idx_ticks = idx_ticks[(idx_ticks >= 0) & (idx_ticks < len(T))]
    idx_ticks = np.unique(idx_ticks)

    ax_top.set_xticks(mdates.date2num(T[idx_ticks]))
    ax_top.set_xticklabels(idx_ticks)

    return ax_top


def add_quicklook_mesh(ax, T, X, Y, args):

    # y_lims = get_quicklook_y_limits(
    #     y_vals=Y,
    #     x_vals=X,
    #     y_lims=args["y_lims"],
    #     use_log=args["use_log_y_scale"],
    # )
    
    y_max = max_quicklook_y(
        y_vals=Y,
        x_vals=X,
        y_max_zone=args["y_max_zone"],
    )

    my_cmap = get_quicklook_colormap()

    # pcolormesh expects C as (len(X), len(T))
    if Y.shape == (len(T), len(X)):
        Y_plot = Y.T
    elif Y.shape == (len(X), len(T)):
        Y_plot = Y
    else:
        raise ValueError(
            f"Unexpected Y shape {Y.shape}. Expected "
            f"({len(T)}, {len(X)}) or ({len(X)}, {len(T)})."
        )

    if args["use_log_y_scale"]:
        qck = ax.pcolormesh(
            T,
            X,
            Y_plot,
            cmap=my_cmap,
            norm=LogNorm(vmin=0., vmax=y_max),
            shading="auto",
        )
        
    else:
        qck = ax.pcolormesh(
            T,
            X,
            Y_plot,
            vmin=0.,
            vmax=y_max,
            cmap=my_cmap,
            shading="auto",
        )

    return qck

def get_quicklook_colormap():

    rgb = color_lib.volkers_rgb()

    my_cmap = make_colormap.custom_rgb(
        rgb,
        name="volkers",
    )

    return my_cmap


def add_colorbar(fig, qck):

    fig.colorbar(
        qck,
        label="Volume Linear Depolarization Ratio",
        extend="both",
        pad=0.02,
    )

    return fig
