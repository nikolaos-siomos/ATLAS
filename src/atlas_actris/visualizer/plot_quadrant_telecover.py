#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from visualizer.plot_utils import export_plot, get_vertical_axis_label


SECTOR_LABELS = {
    "N": "north",
    "E": "east",
    "S": "south",
    "W": "west",
}

SECTOR_COLORS = {
    "N": "tab:blue",
    "E": "tab:orange",
    "S": "tab:green",
    "W": "tab:red",
}


def generate_quadrant_telecover(X, sectors, ranges, args):
    """
    Generate telecover quadrant plot.

    Parameters
    ----------
    X : array-like
        Vertical coordinate, usually height/range in km.
    sectors : dict
        Dictionary with sector IDs as keys: "N", "E", "S", "W".
        Each value is a processed sector dictionary.
    ranges : array-like
        Range values used for non-range-corrected signal option.
    args : dict
        Plot settings, text metadata, and output settings.

    Returns
    -------
    fpath : str
        Output plot path.
    dofl_x : float
        Minimum channel height estimate. np.nan if unavailable.
    """

    data = prepare_telecover_data(
        X=X,
        sectors=sectors,
        ranges=ranges,
        args=args,
    )

    plt.rc("font", size=14)

    fig = plt.figure(figsize=(20.0, 4.0))
    fig.suptitle(args["title"])

    axes = {
        "raw_near": fig.add_axes([0.04, 0.13, 0.19, 0.7]),
        "raw_far": fig.add_axes([0.24, 0.13, 0.07, 0.7]),
        "norm_near": fig.add_axes([0.36, 0.13, 0.19, 0.7]),
        "norm_far": fig.add_axes([0.56, 0.13, 0.07, 0.7]),
        "deviation": fig.add_axes([0.69, 0.13, 0.29, 0.7]),
    }

    raw_near_panel(
        ax=axes["raw_near"],
        X=X,
        data=data,
        args=args,
    )

    raw_far_panel(
        ax=axes["raw_far"],
        X=X,
        data=data,
        args=args,
    )

    norm_near_panel(
        ax=axes["norm_near"],
        X=X,
        data=data,
        args=args,
    )

    norm_far_panel(
        ax=axes["norm_far"],
        X=X,
        data=data,
        args=args,
    )

    deviation_panel(
        ax=axes["deviation"],
        X=X,
        data=data,
        args=args,
    )

    fpath = export_plot(fig, args)

    plt.rcParams.update(plt.rcParamsDefault)

    return fpath, data["dofl_x"]


def prepare_telecover_data(X, sectors, ranges, args):
    """
    Convert processed sector dictionaries into arrays used by the panels.
    """

    X = np.asarray(X)
    ranges = np.asarray(ranges)

    norm_region = args["normalization_region"]
    use_raw = args["plot_raw_signals"]
    dofl_limit = args["relative_deviation_limit"]
    iters = args["iters"]

    raw = {}
    norm = {}
    norm_low = {}
    norm_up = {}

    for sector_id, sector in sectors.items():
        coef = sector["coef"]

        y_raw = np.asarray(sector["y_m"])
        y_sm = np.asarray(sector["y_m_sm"])
        y_low = np.asarray(sector["y_l_sm"])
        y_up = np.asarray(sector["y_u_sm"])

        if use_raw:
            raw[sector_id] = y_raw / np.power(ranges, 2)
        else:
            raw[sector_id] = y_raw

        norm[sector_id] = coef * y_sm
        norm_low[sector_id] = coef * y_low
        norm_up[sector_id] = coef * y_up

    norm_stack = np.vstack([norm[k] for k in norm])
    norm_mean = np.nanmean(norm_stack, axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        rel_dev = {
            k: (norm[k] - norm_mean) / norm_mean
            for k in norm
        }

        rmse = np.sqrt(
            np.nanmean(
                np.vstack([
                    np.power(rel_dev[k], 2)
                    for k in rel_dev
                ]),
                axis=0,
            )
        )

    dofl_ind = np.where((rmse >= dofl_limit) & (X < norm_region[0]))[0]

    if len(dofl_ind) > 0:
        dofl_x = X[dofl_ind[-1]]
    else:
        dofl_x = np.nan

    extra = get_extra_sector_data(
        sectors=sectors,
        raw=raw,
        norm=norm,
        ranges=ranges,
        use_raw=use_raw,
        iters=iters,
    )

    if use_raw:
        raw_label = "Raw BC Signals - Raw Units"
    else:
        raw_label = "RC Signals [A.U.]"

    return {
        "raw": raw,
        "norm": norm,
        "norm_low": norm_low,
        "norm_up": norm_up,
        "norm_mean": norm_mean,
        "rel_dev": rel_dev,
        "rmse": rmse,
        "extra": extra,
        "raw_label": raw_label,
        "dofl_x": dofl_x,
    }


def get_extra_sector_data(sectors, raw, norm, ranges, use_raw, iters):
    """
    Find and prepare the optional extra profile from one sector.
    """

    for sector_id in ["N", "E", "S", "W"]:
        if sector_id not in sectors:
            continue

        sector = sectors[sector_id]

        if not sector.get("has_extra", False):
            continue

        y_extra_raw = np.asarray(sector["y_extra"])
        y_extra_sm = np.asarray(sector["y_extra_sm"])
        coef_extra = sector["coef_extra"]

        if use_raw:
            y_extra_raw = y_extra_raw / np.power(ranges, 2)

        y_extra_norm = coef_extra * y_extra_sm

        if iters is None:
            label_num = ""
        else:
            label_num = str(int(iters) + 1)

        return {
            "exists": True,
            "sector_id": sector_id,
            "raw": y_extra_raw,
            "norm": y_extra_norm,
            "original_norm": norm[sector_id],
            "label": f"{SECTOR_LABELS[sector_id]}{label_num}",
            "label_short": f"{sector_id}{label_num}",
        }

    return {
        "exists": False,
        "sector_id": None,
        "raw": None,
        "norm": None,
        "original_norm": None,
        "label": "",
        "label_short": "",
    }


def raw_near_panel(ax, X, data, args):
    
    x_llim, x_ulim = get_x_lims(args)
    x_ticks, x_tick = get_x_ticks(x_llim, x_ulim, args)
    x_label = get_vertical_axis_label(args['vertical_scale'])
    y_lims = get_y_lims(X, data['raw'], args)

    for sector_id, y in data["raw"].items():
        ax.plot(
            X,
            y,
            color=SECTOR_COLORS.get(sector_id),
            label=SECTOR_LABELS.get(sector_id, sector_id),
            alpha=0.7,
        )

    if args.get("use_last_sector", False) and data["extra"]["exists"]:
        ax.plot(
            X,
            data["extra"]["raw"],
            color="tab:purple",
            label=data["extra"]["label"],
            alpha=0.7,
        )

    if ax.get_legend_handles_labels() != ([], []):
        ax.legend()

    ax.plot([X[0], X[-1]], [0.0, 0.0], color="black")

    ax.set_xticks(x_ticks, labels=x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label, loc="right")

    ax.set_ylim(y_lims)

    ax.set_ylabel(data["raw_label"])
    ax.ticklabel_format(
        axis="y",
        useMathText=True,
        style="sci",
        scilimits=(-1, 1),
    )
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.0))
    ax.grid(which="both")

    return ax


def raw_far_panel(ax, X, data, args):
    
    _, x_ulim = get_x_lims(args)

    y_lims = get_y_lims(X, data['raw'], args)

    for sector_id, y in data["raw"].items():
        ax.plot(
            X,
            y,
            color=SECTOR_COLORS.get(sector_id),
            alpha=0.7,
        )

    if args.get("use_last_sector", False) and data["extra"]["exists"]:
        ax.plot(
            X,
            data["extra"]["raw"],
            color="tab:purple",
            alpha=0.7,
        )

    ax.plot([X[0], X[-1]], [0.0, 0.0], color="black")
    ax.set_xlim([x_ulim, 20.0])
    ax.xaxis.set_minor_locator(MultipleLocator(2.0))

    ax.set_ylim(y_lims)

    ax.set(yticklabels=[])
    ax.grid(which="both")

    return ax


def norm_near_panel(ax, X, data, args):
    
    x_llim, x_ulim = get_x_lims(args)
    x_ticks, x_tick = get_x_ticks(x_llim, x_ulim, args)
    x_label = get_vertical_axis_label(args['vertical_scale'])
    
    norm_region = args["normalization_region"]
    y_lims = get_y_lims(X, data['norm'], args)

    for sector_id, y in data["norm"].items():
        ax.plot(
            X,
            y,
            color=SECTOR_COLORS.get(sector_id),
        )

    if args.get("use_last_sector", False) and data["extra"]["exists"]:
        ax.plot(
            X,
            data["extra"]["norm"],
            color="tab:purple",
        )

    for sector_id in data["norm"]:
        ax.fill_between(
            X,
            data["norm_low"][sector_id],
            data["norm_up"][sector_id],
            color=SECTOR_COLORS.get(sector_id),
            alpha=0.3,
        )

    ax.plot([X[0], X[-1]], [0.0, 0.0], color="black")
    ax.axvspan(
        norm_region[0],
        norm_region[1],
        alpha=0.2,
        facecolor="tab:grey",
    )

    ax.set_xticks(x_ticks, labels=x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label, loc="right")
    ax.set_ylim(y_lims)
    ax.set_ylabel("Normalized RC Signals [A.U.]")
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.0))
    ax.grid(which="both")

    n_llim = np.round(norm_region[0], decimals=2)
    n_ulim = np.round(norm_region[1], decimals=2)

    ax.text(
        0.30 * x_ulim,
        0.90 * y_lims[1],
        f"norm. region: {n_llim} - {n_ulim} km",
        bbox=dict(facecolor="tab:green", alpha=0.1, zorder=9),
    )

    return ax


def norm_far_panel(ax, X, data, args):
    _, x_ulim = get_x_lims(args)
    norm_region = args["normalization_region"]
    y_lims = get_y_lims(X, data['norm'], args)

    for sector_id, y in data["norm"].items():
        ax.plot(
            X,
            y,
            color=SECTOR_COLORS.get(sector_id),
        )

    if args.get("use_last_sector", False) and data["extra"]["exists"]:
        ax.plot(
            X,
            data["extra"]["norm"],
            color="tab:purple",
        )

    for sector_id in data["norm"]:
        ax.fill_between(
            X,
            data["norm_low"][sector_id],
            data["norm_up"][sector_id],
            color=SECTOR_COLORS.get(sector_id),
            alpha=0.3,
        )

    ax.plot([X[0], X[-1]], [0.0, 0.0], color="black")
    ax.axvspan(
        norm_region[0],
        norm_region[1],
        alpha=0.2,
        facecolor="tab:grey",
    )

    ax.set_xlim([x_ulim, 20.0])
    ax.xaxis.set_minor_locator(MultipleLocator(2.0))
    ax.set_ylim(y_lims)
    ax.set(yticklabels=[])
    ax.grid(which="both")

    return ax


def deviation_panel(ax, X, data, args):
    
    x_llim, x_ulim = get_x_lims(args)
    x_ticks, x_tick = get_x_ticks(x_llim, x_ulim, args)
    x_label = get_vertical_axis_label(args['vertical_scale'])
    
    norm_region = args["normalization_region"]
    dofl_limit = args["relative_deviation_limit"]

    for sector_id, y in data["rel_dev"].items():
        ax.plot(
            X,
            y,
            color=SECTOR_COLORS.get(sector_id),
            label="_nolegend_",
        )

    ax.plot(
        X,
        data["rmse"],
        color="tab:cyan",
        label="RMS Sector Diff.",
    )

    if args.get("use_last_sector", False) and data["extra"]["exists"]:
        with np.errstate(divide="ignore", invalid="ignore"):
            y_diff = (
                data["extra"]["norm"] - data["extra"]["original_norm"]
            ) / data["extra"]["original_norm"]

        ax.plot(
            X,
            y_diff,
            color="tab:purple",
            label=f"{data['extra']['label_short']} - {data['extra']['sector_id']}",
        )

    if ax.get_legend_handles_labels() != ([], []):
        ax.legend(loc="lower center")

    for sector_id in data["norm"]:
        with np.errstate(divide="ignore", invalid="ignore"):
            y_low = (
                data["norm_low"][sector_id] - data["norm_mean"]
            ) / data["norm_mean"]
            y_up = (
                data["norm_up"][sector_id] - data["norm_mean"]
            ) / data["norm_mean"]
        
        ax.fill_between(
            X,
            y_low,
            y_up,
            color=SECTOR_COLORS.get(sector_id),
            alpha=0.3,
        )

    ax.plot(X, np.zeros(X.shape), "--", color="black", zorder=10)
    ax.plot(
        X,
        -dofl_limit * np.ones(X.shape),
        "--",
        color="black",
        zorder=10,
        alpha=0.7,
    )
    ax.plot(
        X,
        dofl_limit * np.ones(X.shape),
        "--",
        color="black",
        zorder=10,
        alpha=0.7,
    )

    ax.axvspan(
        norm_region[0],
        norm_region[1],
        alpha=0.2,
        facecolor="tab:grey",
    )

    ax.set_xticks(x_ticks, labels=x_ticks)
    ax.set_xlim([x_llim, x_ulim])
    ax.set_xlabel(x_label)

    y_ticks = np.round(
        np.arange(-0.20, 0.20 + dofl_limit, dofl_limit),
        decimals=2,
    )
    ax.set_yticks(
        y_ticks,
        labels=["%.2f" % tick for tick in y_ticks],
    )
    ax.set_ylim([y_ticks[0], y_ticks[-1]])
    ax.set_ylabel("Relative Sector Deviation")
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick / 5.0))
    ax.grid(which="both")

    if np.isfinite(data["dofl_x"]):
        ax.text(
            0.30,
            0.90,
            f"min_channel_height: {data['dofl_x']:.2f} km",
            transform=ax.transAxes,
            bbox=dict(facecolor="tab:green", alpha=0.1),
            zorder=9,
        )

    return ax


def get_raw_ulim(raw, X, x_llim, x_ulim):
    
    mask = (X >= x_llim) & (X <= x_ulim)

    vals = []

    for y in raw.values():
        vals.append(y[mask])

    if len(vals) == 0:
        return np.nan

    return np.nanmax(vals)


def get_x_ticks(x_llim, x_ulim, args):
    
    # Set the x_tick depending on the telescope_type
    if args['atlas_channel_id'][4] in ['l', 'm', 'n', 'x', 'y', 'z']:
        x_tick = 0.5
    else:
        x_tick = 1.
        
    x_ticks = np.arange(
        x_tick * np.ceil(x_llim / x_tick),
        x_tick * (np.floor(x_ulim / x_tick) + 1.0),
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

    return np.round(x_ticks, decimals=2), x_tick


def get_x_lims(args):
    
    # Set the lower x limit
    x_llim = 0.
        
    # Set the upper x limit depending on the telescope_type
    if args['atlas_channel_id'][4] in ['l', 'm', 'n', 'x', 'y', 'z']:
        x_ulim = 2.5
    else:
        x_ulim = 5.
    
    return x_llim, x_ulim


def get_y_lims(X, Y, args):
    
    vals = [np.asarray(y) for y in Y.values()]

    if len(vals) == 0:
        return [0.0, 1.0]

    # Mask for min/max values
    nr_mask = (X >= 0.05) & (X < args['near_range_upper_limit'])
    
    vals = np.vstack(vals)
    y_max = np.nanmax(vals[:,nr_mask])

    if not np.isfinite(y_max) or y_max <= 0.0:
        return [0.0, 1.0]

    return [-0.1 * y_max, 1.2 * y_max]
