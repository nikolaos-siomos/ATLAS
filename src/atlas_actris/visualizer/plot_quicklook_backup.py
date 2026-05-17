#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 00:22:21 2026

@author: nikos
"""

import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

from visualizer.plotting import color_lib, make_colormap


def generate_plot(
        dir_out, fname, title, dpi_val, use_log, delta_t,
        t_vals, y_vals, z_vals,
        x_lbin, x_ubin, y_lbin, y_ubin,
        y_llim, y_ulim, z_llim, z_ulim,
        y_label, t_tick, x_tick, y_tick
        ):

    args = {
        "dir_out": dir_out,
        "fname": fname,
        "title": title,
        "dpi_val": dpi_val,
        "use_log": use_log,
        "delta_t": delta_t,
        "x_lbin": x_lbin,
        "x_ubin": x_ubin,
        "y_lbin": y_lbin,
        "y_ubin": y_ubin,
        "y_llim": y_llim,
        "y_ulim": y_ulim,
        "z_llim": z_llim,
        "z_ulim": z_ulim,
        "y_label": y_label,
        "t_tick": t_tick,
        "x_tick": x_tick,
        "y_tick": y_tick,
    }

    fig_coords = [0.07, 0.13, 0.99, 0.74]

    t_vals, y_vals, z_vals, nodes, x_ubin = prepare_quicklook_data(
        t_vals=t_vals,
        y_vals=y_vals,
        z_vals=z_vals,
        delta_t=delta_t,
        x_ubin=x_ubin,
    )

    X, Y, Z = get_plot_arrays(
        t_vals=t_vals,
        y_vals=y_vals,
        z_vals=z_vals,
        x_lbin=x_lbin,
        x_ubin=x_ubin,
        y_lbin=y_lbin,
        y_ubin=y_ubin,
    )

    fig = plt.figure(figsize=(10., 5.))

    quicklook_panel(
        fig=fig,
        fig_coords=fig_coords,
        X=X,
        Y=Y,
        Z=Z,
        nodes=nodes,
        args=args,
    )

    fpath = export_quicklook(fig=fig, args=args)

    return fpath


def prepare_quicklook_data(t_vals, y_vals, z_vals, delta_t, x_ubin):
    """
    Prepare time, height/range, and signal arrays for pcolormesh.

    This function:
    - transposes z_vals to y/time orientation,
    - detects temporal gaps,
    - inserts NaNs at gaps,
    - extends the final time and y edges for pcolormesh.
    """

    z_vals = z_vals.T

    t_resol = get_time_resolution(
        t_vals=t_vals,
        delta_t=delta_t,
    )

    nodes = get_gap_nodes(
        t_vals=t_vals,
        t_resol=t_resol,
    )

    if len(nodes) > 0:
        t_vals, z_vals, x_ubin = insert_time_gaps(
            t_vals=t_vals,
            z_vals=z_vals,
            nodes=nodes,
            t_resol=t_resol,
            x_ubin=x_ubin,
        )

    t_vals = extend_time_edges(
        t_vals=t_vals,
        t_resol=t_resol,
    )

    y_vals = extend_y_edges(
        y_vals=y_vals,
    )

    return t_vals, y_vals, z_vals, nodes, x_ubin


def quicklook_panel(fig, fig_coords, X, Y, Z, nodes, args):

    ax = fig.add_axes(fig_coords)

    ax.set_title(args["title"], pad=5)

    set_axis_labels(
        ax=ax,
        y_label=args["y_label"],
    )

    set_y_axis(
        ax=ax,
        y_llim=args["y_llim"],
        y_ulim=args["y_ulim"],
        y_tick=args["y_tick"],
    )

    set_time_axis(
        ax=ax,
        t_tick=args["t_tick"],
    )

    if len(nodes) == 0:
        add_frame_axis(
            ax=ax,
            x_lbin=args["x_lbin"],
            x_ubin=args["x_ubin"],
            x_tick=args["x_tick"],
        )

    qck = add_quicklook_mesh(
        ax=ax,
        X=X,
        Y=Y,
        Z=Z,
        z_llim=args["z_llim"],
        z_ulim=args["z_ulim"],
        use_log=args["use_log"],
    )

    add_colorbar(
        fig=fig,
        qck=qck,
    )

    return ax


def get_time_resolution(t_vals, delta_t):

    dt_min = np.nanmin(t_vals[1:] - t_vals[:-1])
    t_bin = np.timedelta64(int(round(np.nanmax(delta_t))), "s")

    if dt_min > t_bin or dt_min <= 1.2 * t_bin:
        t_resol = 1.2 * dt_min
    else:
        t_resol = 1.2 * t_bin

    return t_resol


def get_gap_nodes(t_vals, t_resol):

    nodes = np.where(t_vals[1:] - t_vals[:-1] > t_resol)[0] + 1

    return nodes


def insert_time_gaps(t_vals, z_vals, nodes, t_resol, x_ubin):

    t_nodes = t_vals[nodes - 1] + t_resol

    z_vals = np.insert(
        z_vals,
        nodes,
        np.nan,
        axis=1,
    )

    t_vals = np.insert(
        t_vals,
        nodes,
        t_nodes,
    )

    x_ubin = x_ubin + len(nodes)

    return t_vals, z_vals, x_ubin


def extend_time_edges(t_vals, t_resol):

    t_vals = np.hstack(
        (
            t_vals,
            [t_vals[-1] + t_resol],
        )
    )

    return t_vals


def extend_y_edges(y_vals):

    y_vals = np.hstack(
        (
            y_vals,
            [y_vals[-1] + y_vals[1] - y_vals[0]],
        )
    )

    return y_vals


def get_plot_arrays(t_vals, y_vals, z_vals, x_lbin, x_ubin, y_lbin, y_ubin):

    X, Y = np.meshgrid(
        t_vals[slice(x_lbin, x_ubin + 2)],
        y_vals[slice(y_lbin, y_ubin + 2)],
    )

    Z = z_vals[
        slice(y_lbin, y_ubin + 1),
        slice(x_lbin, x_ubin + 1),
    ]

    return X, Y, Z


def set_axis_labels(ax, y_label):

    ax.set_xlabel("Time UTC")
    ax.set_ylabel(y_label)

    return ax


def set_y_axis(ax, y_llim, y_ulim, y_tick):

    y_ticks = get_y_ticks(
        y_llim=y_llim,
        y_ulim=y_ulim,
        y_tick=y_tick,
    )

    ax.set_yticks(y_ticks, labels=y_ticks)
    ax.set_ylim([y_llim, y_ulim])

    return ax


def get_y_ticks(y_llim, y_ulim, y_tick):

    y_ticks = np.arange(
        y_tick * np.ceil(y_llim / y_tick),
        y_tick * (np.floor(y_ulim / y_tick) + 1.),
        y_tick,
    )

    if np.abs(y_llim - y_ticks[0]) < y_tick * 0.25:
        y_ticks[0] = y_llim
    else:
        y_ticks = np.hstack((y_llim, y_ticks))

    if np.abs(y_ulim - y_ticks[-1]) < y_tick * 0.25:
        y_ticks[-1] = y_ulim
    else:
        y_ticks = np.hstack((y_ticks, y_ulim))

    y_ticks = np.round(y_ticks, decimals=2)

    return y_ticks


def set_time_axis(ax, t_tick):

    ax.xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M")
    )

    major_locator = mdates.SecondLocator(
        interval=60 * int(t_tick)
    )
    major_locator.MAXTICKS = 10000
    ax.xaxis.set_major_locator(major_locator)

    minor_locator = mdates.SecondLocator(
        interval=15 * int(t_tick)
    )
    minor_locator.MAXTICKS = 10000
    ax.xaxis.set_minor_locator(minor_locator)

    plt.xticks(rotation=35)

    return ax


def add_frame_axis(ax, x_lbin, x_ubin, x_tick):

    ax1 = ax.twiny()

    x_ticks = np.arange(
        x_lbin,
        x_ubin + x_tick,
        x_tick,
    ).astype(int)

    ax1.set_xticks(
        x_ticks,
        labels=x_ticks,
    )

    return ax1


def add_quicklook_mesh(ax, X, Y, Z, z_llim, z_ulim, use_log):

    my_cmap = get_quicklook_colormap()

    if use_log:
        qck = ax.pcolormesh(
            X,
            Y,
            Z,
            vmin=z_llim,
            vmax=z_ulim,
            cmap=my_cmap,
            norm="log",
        )
    else:
        qck = ax.pcolormesh(
            X,
            Y,
            Z,
            vmin=z_llim,
            vmax=z_ulim,
            cmap=my_cmap,
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
        label="Range-corrected Signal [A.U.]",
        extend="both",
        pad=0.02,
    )

    return fig


def export_quicklook(fig, args):

    fpath = os.path.join(
        args["dir_out"],
        args["fname"],
    )

    fig.savefig(
        fpath,
        dpi=args["dpi_val"],
    )

    fig.clf()
    plt.close()

    return fpath