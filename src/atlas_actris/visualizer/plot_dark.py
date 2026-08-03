#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dark-test diagnostic plotting."""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter

from visualizer import color_lib, make_colormap
from visualizer.plot_utils import add_fitting_suptitle, export_plot


_VOLKERS_CMAP = make_colormap.custom_rgb(
    color_lib.volkers_rgb(),
    name="volkers",
)


def _values(array):
    """Return an in-memory NumPy representation of an xarray/NumPy object."""
    if hasattr(array, "compute"):
        array = array.compute()
    if hasattr(array, "values"):
        array = array.values
    return np.asarray(array)


def _elapsed_minutes(signal, n_profiles):
    """Return elapsed minutes for profile coloring and the total duration."""
    if not hasattr(signal, "coords") or "time" not in signal.coords:
        elapsed = np.arange(n_profiles, dtype=float)
        return elapsed, float(elapsed[-1]) if n_profiles else 0.0

    time = _values(signal.coords["time"]).astype("datetime64[ns]")
    if time.size != n_profiles or time.size == 0:
        elapsed = np.arange(n_profiles, dtype=float)
        return elapsed, float(elapsed[-1]) if n_profiles else 0.0

    elapsed = np.asarray(
        (time - time[0]) / np.timedelta64(1, "m"),
        dtype=float,
    )
    return elapsed, float(elapsed[-1])


def _plot_colored_profiles(ax, x, y, color_values, cmap, norm):
    """Draw all profiles as one collection without changing their appearance."""
    x = np.asarray(x, dtype=float).squeeze()
    y = np.asarray(y, dtype=float)

    if x.ndim != 1:
        raise ValueError(f"Expected one-dimensional x values, got shape {x.shape}.")
    if y.ndim != 2:
        raise ValueError(f"Expected two-dimensional profiles, got shape {y.shape}.")
    if y.shape[1] != x.size:
        raise ValueError(
            f"Profile/x mismatch: {y.shape[1]} profile points versus {x.size} x values."
        )

    segments = np.empty((y.shape[0], y.shape[1], 2), dtype=float)
    segments[:, :, 0] = x[np.newaxis, :]
    segments[:, :, 1] = y

    collection = LineCollection(
        segments,
        colors=cmap(norm(color_values)),
        linewidths=plt.rcParams["lines.linewidth"],
        antialiaseds=plt.rcParams["lines.antialiased"],
        capstyle=plt.rcParams["lines.solid_capstyle"],
        joinstyle=plt.rcParams["lines.solid_joinstyle"],
    )
    ax.add_collection(collection)
    return collection



def _scientific_y_axis(fig, ax):
    """Use Matplotlib scientific notation and return its order multiplier."""
    formatter = ScalarFormatter(useMathText=True, useOffset=False)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(formatter)

    # Let Matplotlib determine the scientific order, then move the generated
    # multiplier from the axis offset area into the y-axis label.
    fig.canvas.draw()
    offset = ax.yaxis.get_offset_text()
    offset_text = offset.get_text().replace(r"\times", r"\times\,")
    offset.set_visible(False)

    return offset_text


def generate_plot(bins_dict, x_dict, y_dict, m_dict, args):
    """Generate the dark-test summary plot.

    The full six-panel layout is reserved for the exact ``drk`` key on
    analog channels. Photon-counting channels and auxiliary ``drk_*``
    measurements use the compact raw/background/zero-bin layout.
    """
    extended_dark_analysis = bool(args.get("extended_dark_analysis", False))

    if not extended_dark_analysis:
        fig = plt.figure(figsize=(15, 3.8))
        add_fitting_suptitle(fig, args["title"], y=0.965)

        fig_x = 0.26
        fig_y = 0.57
        x_edges = (0.048, 0.37, 0.70)
        y_edge = 0.13

        ax1_coords = [x_edges[0], y_edge, fig_x, fig_y]
        ax2_coords = [x_edges[1], y_edge, fig_x, fig_y]
        ax3_coords = [x_edges[2], y_edge, fig_x, fig_y]

        # Original first-row, first-column panel.
        plot_raw_multi(
            fig, ax1_coords,
            bins_dict["av"], x_dict["av"], y_dict["av"],
            args["xlims_av"], args["xlims_range_av"], args["ylims_av"],
            "raw", args,
        )

        # Original first-row, second-column panel.
        plot_raw_multi(
            fig, ax2_coords,
            bins_dict["bc"], x_dict["bc"], y_dict["bc"],
            args["xlims_bc"], args["xlims_range_bc"], args["ylims_bc"],
            "bc", args,
        )

        # Original second-row, second-column zero-bin panel.
        plot_raw_multi(
            fig, ax3_coords,
            bins_dict["bc"], x_dict["bc"], y_dict["bc"],
            args["xlims_zb"], args["xlims_range_zb"], args["ylims_zb"],
            "bc", args, activate_colorbar=True,
        )

        return export_plot(fig, args)

    fig = plt.figure(figsize=(15, 6.4))
    add_fitting_suptitle(fig, args["title"], y=0.975)

    fig_x = 0.26
    fig_y = 0.31
    x_edges = (0.062, 0.39, 0.720)
    y_edges = (0.5, 0.075)

    ax1_coords = [x_edges[0], y_edges[0], fig_x, fig_y]
    ax2_coords = [x_edges[1], y_edges[0], fig_x, fig_y]
    ax3_coords = [x_edges[2], y_edges[0], fig_x, fig_y]
    ax4_coords = [x_edges[0], y_edges[1], fig_x, fig_y]
    ax5_coords = [x_edges[1], y_edges[1], fig_x, fig_y]
    ax6_coords = [x_edges[2], y_edges[1], fig_x, fig_y]

    plot_raw_multi(
        fig, ax1_coords,
        bins_dict["av"], x_dict["av"], y_dict["av"],
        args["xlims_av"], args["xlims_range_av"], args["ylims_av"],
        "raw", args, disable_x1_label=True,
    )

    plot_raw_multi(
        fig, ax2_coords,
        bins_dict["bc"], x_dict["bc"], y_dict["bc"],
        args["xlims_bc"], args["xlims_range_bc"], args["ylims_bc"],
        "bc", args, disable_x1_label=True,
    )

    plot_normalized_sm_deviation(
        fig=fig,
        ax_coords=ax3_coords,
        bins=bins_dict["sm"],
        ranges_km=x_dict["sm"],
        signal=y_dict["sm"],
        molecular=m_dict["sm"],
        xlims_range=args["xlims_rc"],
        ylims=args["ylims_rc"],
        args=args,
    )

    plot_raw_multi(
        fig, ax4_coords,
        bins_dict["sm"], x_dict["sm"], y_dict["sm"],
        args["xlims_sm"], args["xlims_range_sm"], args["ylims_sm"],
        "sm_bc", args, disable_x2_label=True, show_stats_region=True,
    )

    plot_raw_multi(
        fig, ax5_coords,
        bins_dict["bc"], x_dict["bc"], y_dict["bc"],
        args["xlims_zb"], args["xlims_range_zb"], args["ylims_zb"],
        "bc", args, disable_x2_label=True, activate_colorbar=True,
    )

    plot_stats_table(fig, ax6_coords, args)

    return export_plot(fig, args)

def collect_stats(args):
    y_units = "mV" if args["channel_mode"] == "a" else "MHz"
    
    z_alias = args["vertical_scale"]
    
    stats_names = [
        "Stats Region [km]",
        f"Mean offset [{y_units}]",
        f"Noise / bin / shot [{y_units}]\n(stdev of full period)",
        f"Slope over {z_alias} [{y_units} km$^{{-1}}$]",
        f"Temporal slope [{y_units} s$^{{-1}}$]",
        "Is noise Gaussian?",
        "N: bins, avg. signals, shots",
    ]

    noise_per_bin_per_shot = args.get(
        "noise_per_bin_per_shot",
        args.get("noise_per_bin", np.nan),
    )

    bins = int(args.get("bins", 0))
    profiles = int(args.get("profiles", 0))
    shots = int(round(float(
        args.get("all_shots", args.get("shots", args.get("sample", 0)))
    )))

    stats_values = [
        f"{args['stats_range'][0]:.1f} - {args['stats_range'][1]:.1f}",
        f"{args['baseline_offset']:.3e}",
        f"{noise_per_bin_per_shot:.2e}",
        str(args["vert_slope_flag"]),
        str(args["temp_slope_flag"]),
        str(args["gaussian_noise_flag"]),
        f"{bins}, {profiles}, {shots}",
    ]
    return stats_names, stats_values


def plot_normalized_sm_deviation(
    fig, ax_coords, bins, ranges_km, signal, molecular,
    xlims_range, ylims, args,
):
    """Plot ``(sm - mean_time(sm)) / molecular`` for every time profile."""
    ax = fig.add_axes(ax_coords)

    if not isinstance(signal, xr.DataArray):
        raise TypeError(
            "The smoothed background-corrected signal must remain an "
            "xarray.DataArray after smoothing so its time coordinate is retained."
        )
    if "time" not in signal.dims:
        raise ValueError(
            f"The smoothed signal has no 'time' dimension: {signal.dims}."
        )

    vertical_dims = [dim for dim in signal.dims if dim != "time"]
    if len(vertical_dims) != 1:
        raise ValueError(
            "Expected one vertical dimension in addition to 'time', "
            f"got {signal.dims}."
        )

    vertical_dim = vertical_dims[0]
    signal = signal.transpose("time", vertical_dim)
    x_bins = _values(bins).astype(float).squeeze()
    x_range = _values(ranges_km).astype(float).squeeze()

    if x_bins.ndim != 1 or x_bins.size != signal.sizes[vertical_dim]:
        raise ValueError(
            "Smoothed signal/bin mismatch: "
            f"signal shape {signal.shape}, bin shape {x_bins.shape}."
        )
    if x_range.ndim != 1 or x_range.size != signal.sizes[vertical_dim]:
        raise ValueError(
            "Smoothed signal/range mismatch: "
            f"signal shape {signal.shape}, range shape {x_range.shape}."
        )

    range_lower, range_upper = map(float, xlims_range)
    range_mask = (
        np.isfinite(x_range)
        & (x_range >= min(range_lower, range_upper))
        & (x_range <= max(range_lower, range_upper))
    )
    range_indices = np.flatnonzero(range_mask)
    if range_indices.size == 0:
        raise ValueError(
            "No vertical-scale values fall inside the normalized-panel "
            f"limits {xlims_range}."
        )

    first_index = int(range_indices[0])
    last_index = int(range_indices[-1])
    xlims_bins = [
        float(x_bins[first_index]),
        float(x_bins[last_index]),
    ]

    if isinstance(molecular, xr.DataArray):
        molecular = molecular.squeeze(drop=True)
        if molecular.ndim != 1:
            raise ValueError(
                f"Expected a one-dimensional molecular profile, got {molecular.dims}."
            )
        molecular = molecular.rename({molecular.dims[0]: vertical_dim})
        molecular = molecular.assign_coords({vertical_dim: signal[vertical_dim]})
    else:
        molecular_values = _values(molecular).astype(float).squeeze()
        molecular = xr.DataArray(
            molecular_values,
            dims=(vertical_dim,),
            coords={vertical_dim: signal[vertical_dim]},
        )

    if molecular.sizes[vertical_dim] != signal.sizes[vertical_dim]:
        raise ValueError(
            "Smoothed signal/molecular mismatch: "
            f"{signal.sizes[vertical_dim]} signal bins versus "
            f"{molecular.sizes[vertical_dim]} molecular points."
        )

    deviation = signal - signal.mean(dim="time", skipna=True)
    valid_molecular = np.isfinite(molecular) & (molecular != 0.0)
    normalized_deviation = (deviation / molecular).where(valid_molecular)
    normalized_values = _values(normalized_deviation).astype(float)

    n_time = normalized_deviation.sizes["time"]
    elapsed, duration = _elapsed_minutes(normalized_deviation, n_time)

    cmap = _VOLKERS_CMAP
    vmax = duration if duration > 0 else max(n_time - 1, 1)
    norm = Normalize(vmin=0.0, vmax=vmax)
    color_values = elapsed if duration > 0 else np.arange(n_time, dtype=float)

    _plot_colored_profiles(
        ax=ax,
        x=x_bins,
        y=normalized_values,
        color_values=color_values,
        cmap=cmap,
        norm=norm,
    )

    relative_limit = float(args["relative_molecular_deviation"])
    if not np.isfinite(relative_limit) or relative_limit <= 0.0:
        raise ValueError(
            "relative_molecular_deviation must be a finite positive float. "
            f"Received {relative_limit!r}."
        )

    ax.axhline(0.0, linewidth=1.0, color="black", alpha=0.6)
    ax.axhline(
        relative_limit,
        linewidth=1.0,
        linestyle="--",
        color="black",
        alpha=0.8,
    )
    ax.axhline(
        -relative_limit,
        linewidth=1.0,
        linestyle="--",
        color="black",
        alpha=0.8,
    )
    ax.set_xlim(xlims_bins)
    ax.set_ylim(-2.0 * relative_limit, 2.0 * relative_limit)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="x", which="minor", length=3)
    ax.tick_params(axis="y", which="minor", length=3)

    ax_top = ax.twiny()
    ax_top.set_xlim(xlims_range)

    y_offset = _scientific_y_axis(fig, ax)
    x_label_1, x_label_2, y_label, _ = make_labels(
        "normalized_sm_deviation", args["channel_mode"], y_offset
    )
    ax.set_xlabel(x_label_1)
    ax_top.set_xlabel(x_label_2, labelpad=5)
    ax.set_ylabel(y_label)
    ax.grid(which="both")

    max_vertical = float(args.get("max_channel_vertical_scale", np.nan))
    vertical_alias = args.get(
        "vertical_scale_alias",
        args.get("vertical_scale", "range"),
    )
    max_vertical_text = (
        f"Max channel {vertical_alias}: {max_vertical:.2f} km"
        if np.isfinite(max_vertical)
        else f"Max channel {vertical_alias}: Not exceeded"
    )
    text = ax.text(
        0.03,
        0.96,
        max_vertical_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        zorder=100,
        bbox=dict(
            facecolor="tab:grey",
            alpha=0.22,
            edgecolor="none",
        ),
    )
    
    text.get_bbox_patch().set_zorder(99)

    # The shared elapsed-time colorbar is displayed only on plot 2-2.
    # Plot 1-3 uses the same colormap and normalization without duplicating it.

def plot_raw_multi(
    fig, ax_coords, bins, ranges_km, signal,
    xlims, xlims_upper, ylims, signal_type, args,
    disable_x1_label=False, disable_x2_label=False,
    activate_colorbar=False, show_stats_region=False,
):
    ax = fig.add_axes(ax_coords)

    x_bins = _values(bins).astype(float)
    x_km = _values(ranges_km).astype(float)
    y = _values(signal).astype(float)

    if y.ndim == 1:
        y = y[np.newaxis, :]
    if y.ndim != 2:
        raise ValueError(f"Expected a 2-D signal for '{signal_type}', got shape {y.shape}.")
    if y.shape[-1] != x_bins.size:
        raise ValueError(
            f"Signal/bin mismatch for '{signal_type}': {y.shape[-1]} signal bins "
            f"versus {x_bins.size} coordinates."
        )

    n_time = y.shape[0]
    elapsed, duration = _elapsed_minutes(signal, n_time)

    cmap = _VOLKERS_CMAP
    vmax = duration if duration > 0 else max(n_time - 1, 1)
    norm = Normalize(vmin=0.0, vmax=vmax)
    color_values = elapsed if duration > 0 else np.arange(n_time, dtype=float)

    _plot_colored_profiles(
        ax=ax,
        x=x_bins,
        y=y,
        color_values=color_values,
        cmap=cmap,
        norm=norm,
    )

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="x", which="minor", length=3)
    ax.tick_params(axis="y", which="minor", length=3)

    ax_top = ax.twiny()
    ax_top.set_xlim(xlims_upper)

    y_offset = _scientific_y_axis(fig, ax)
    x_label_1, x_label_2, y_label, _ = make_labels(
        signal_type, args["channel_mode"], y_offset
    )

    if not disable_x1_label:
        ax.set_xlabel(x_label_1)
    if not disable_x2_label:
        ax_top.set_xlabel(x_label_2, labelpad=5)
    ax.set_ylabel(y_label)

    # The upper x-axis is already in km.
    ax_top.axvline(0.0, linewidth=2.0, color="tab:blue", alpha=0.3)

    if show_stats_region:
        ax_top.axvspan(
            args["stats_range"][0], args["stats_range"][1],
            alpha=0.2, zorder=10, facecolor="tab:grey",
        )
        ax.text(
            0.25, 0.07, "Stats region", transform=ax.transAxes,
            bbox=dict(facecolor="tab:grey", alpha=0.22, zorder=3),
        )

    ax.grid(which="both")

    if activate_colorbar:
        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.06, fraction=0.035)
        cbar.set_label("Elapsed time [min]", labelpad=5)


def plot_stats_table(fig, ax_coords, args, title=None):
    stats_names, stats_values = collect_stats(args)

    # Make the table slightly wider and extend it mainly to the left while
    # keeping its right edge aligned with the original panel.
    table_coords = list(ax_coords)
    extra_width = 0.018
    table_coords[0] -= extra_width
    table_coords[2] += extra_width

    ax = fig.add_axes(table_coords)
    ax.axis("off")

    table = ax.table(
        cellText=[[name, value] for name, value in zip(stats_names, stats_values)],
        colLabels=["Statistics on region", "Value"],
        cellLoc="left", colLoc="left", loc="center",
        # Move the separator slightly to the right by widening the first column.
        colWidths=[0.60, 0.40],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.5)

    for (row, column), cell in table.get_celld().items():
        cell.set_linewidth(0.8)

        # Keep the left-column text close to the border without touching it.
        if column == 0:
            cell.PAD = 0.04

        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#eaeaea")

    # Row 3 is the two-line noise entry (row 0 contains the headers).
    # Double the height of both cells so the second line fits comfortably.
    for column in (0, 1):
        if (3, column) in table.get_celld():
            cell = table[(3, column)]
            cell.set_height(cell.get_height() * 2.0)

    # Rows are offset by one because row 0 contains the column headers.
    slope_rows = {
        4: bool(args.get("vert_slope_sign", False)),
        5: bool(args.get("temp_slope_sign", False)),
    }
    for row, is_significant in slope_rows.items():
        if (row, 1) in table.get_celld():
            table[(row, 1)].get_text().set_color(
                "red" if is_significant else "green"
            )

    if title is not None:
        ax.set_title(title, pad=23)
    return ax


def make_labels(signal_type, channel_mode, y_offset=""):
    y_units = "mV" if channel_mode == "a" else "MHz"
    display_units = " ".join(value for value in (y_units, y_offset) if value)

    if signal_type == "raw":
        return "Bins", "Range [km]", f"Raw Signal [{display_units}]", y_units
    if signal_type == "bc":
        return "Bins", "Range [km]", f"BG Corr. Signal [{display_units}]", y_units
    if signal_type == "sm_bc":
        return "Bins", "Range [km]", f"Smoothed BG Corr. Signal\n[{display_units}]", y_units
    if signal_type == "rc":
        rc_units = " ".join(value for value in ("AU", y_offset) if value)
        return "Range [km]", "", f"RC Signal\n[{rc_units}]", "AU"
    if signal_type == "normalized_sm_deviation":
        display_units = " ".join(value for value in ("", y_offset) if value)
        return (
            "Bins",
            "Range [km]",
            f"Relative deviations from\nmolecular[{display_units}]",
            "Relative",
        )

    raise ValueError(
        "Unsupported signal_type "
        f"'{signal_type}'. Expected raw, bc, sm_bc, rc, or normalized_sm_deviation."
    )
