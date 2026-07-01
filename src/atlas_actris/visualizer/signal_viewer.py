#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:19:58 2022

@author: nick
"""

import os
import gc
import warnings
import numpy as np
import panel as pn
import holoviews as hv
from bokeh.models import (
    LinearColorMapper, ColorBar, FixedTicker,
    LinearAxis, BoxAnnotation, Span,
)
from matplotlib import colors as mpl_colors
from utils.error_classes import CustomWarning

from utils.printouts import print_header
from visualizer.check import check_channels
from visualizer import color_lib, make_colormap
from processor.packaging import collect_metadata
from visualizer.make_text import GenerateText, Libraries
from visualizer.plot_utils import (
    prepare_folder,
    convert_m_to_km, 
    )

# Ignores all warnings --> they are not printed in terminal
warnings.filterwarnings('ignore')

# Slightly enlarge Bokeh/Panel toolbar buttons in the exported HTML files.
# The exact CSS class names differ a bit between Bokeh versions, so the
# selectors below intentionally cover the common Bokeh 2.x/3.x toolbar classes.
_TOOLBAR_CSS = """
.bk-toolbar .bk-btn,
.bk-Toolbar .bk-btn,
.bk-toolbar-button .bk-btn,
.bk-tool-button .bk-btn {
    width: 42px !important;
    height: 42px !important;
    min-width: 42px !important;
    min-height: 42px !important;
}

.bk-tool-icon,
.bk-tool-button,
.bk-toolbar-button {
    width: 42px !important;
    height: 42px !important;
}

.bk-tool-icon {
    background-size: 30px 30px !important;
}

.bk-toolbar,
.bk-Toolbar {
    gap: 6px !important;
}
"""


def _ensure_panel_css():
    """Register custom CSS once for larger exported toolbar buttons."""

    if _TOOLBAR_CSS not in pn.config.raw_css:
        pn.config.raw_css.append(_TOOLBAR_CSS)


def get_quicklook_colormap():

    rgb = color_lib.volkers_rgb()

    my_cmap = make_colormap.custom_rgb(
        rgb,
        name="volkers",
    )

    return my_cmap

def _get_vertical_bin_dim(vertical_scale):
    """Return the vertical/bin dimension name from a 1D vertical scale."""

    if len(vertical_scale.dims) != 1:
        raise ValueError(
            "The selected vertical scale must be 1D after selecting one channel. "
            f"Found dimensions: {vertical_scale.dims}"
        )

    return vertical_scale.dims[0]


def _slice_vertical_scale_only(sig_ch, vertical_scale_ch, x_lims):
    """
    Slice one channel lazily using only the eager vertical scale.

    This intentionally does not inspect the profile values.  Inspecting profile
    validity with operations such as da.notnull().any(dim='time') would trigger
    a calculation over the lazy profile.  For quicklooks, filtering finite
    vertical coordinates within x_lims is enough and keeps pcolormesh happy.
    """

    bin_dim = _get_vertical_bin_dim(vertical_scale_ch)
    x_vals_all = np.asarray(vertical_scale_ch.values)

    if x_lims is None or len(x_lims) == 0:
        mask = np.isfinite(x_vals_all)
    else:
        mask = (
            np.isfinite(x_vals_all)
            & (x_vals_all >= x_lims[0])
            & (x_vals_all <= x_lims[1])
        )

    if not np.any(mask):
        selected_id = None
        for coord_name in ["channel", "pair"]:
            if coord_name in sig_ch.coords:
                try:
                    selected_id = sig_ch.coords[coord_name].values
                except Exception:
                    selected_id = sig_ch.coords[coord_name]
                break

        raise ValueError(
            "No finite vertical-scale bins were found inside x_lims "
            f"({x_lims}) for selection {selected_id}."
        )

    sig_ch = sig_ch.isel({bin_dim: mask})
    vertical_scale_ch = vertical_scale_ch.isel({bin_dim: mask})

    return sig_ch, vertical_scale_ch


def _to_numpy_selected(da):
    """
    Materialize only the already-selected quicklook slice.

    Matplotlib cannot draw a Dask-backed array directly, so a computation is
    still necessary.  The important point is that channel/time/bin selection has
    already happened before this function is called.
    """

    if hasattr(da, "compute"):
        da = da.compute()

    return np.asarray(da.values)


def _select_channel_if_present(arr, ch):
    """Select a channel from an auxiliary 1D/2D coordinate if it has one."""

    if arr is None:
        return None

    if "channel" in arr.dims or "channel" in arr.coords:
        return arr.sel(channel=ch)

    return arr


def _prepare_multiline_profile(sig_ch, vertical_scale_ch, bins_ch=None, x_name="bins"): 
    """
    Prepare one channel for an hvPlot profile plot.

    The lower x-axis uses the dedicated bins entry from the data_pack, when it
    is available.  The selected vertical scale is not used as the plotted x
    coordinate; it is only used by a Bokeh hook to add a second x-axis at the
    top of the plot.

    The output is a 2D DataArray with dimensions time x bins plus a 1D plotting
    coordinate named x_name attached along the bin dimension.
    """

    bin_dim = _get_vertical_bin_dim(vertical_scale_ch)

    if "time" not in sig_ch.dims:
        raise ValueError(
            "Cannot create a multi-line signal plot because the selected "
            f"channel has no 'time' dimension. Found dimensions: {sig_ch.dims}"
        )

    if bin_dim not in sig_ch.dims:
        raise ValueError(
            "Cannot create a multi-line signal plot because the vertical-scale "
            f"dimension {bin_dim!r} is not present in the signal dimensions "
            f"{sig_ch.dims}."
        )

    # Keep only the two dimensions that should be plotted. Any scalar channel
    # coordinate is kept as metadata, but the channel dimension has already been
    # removed by profiles.sel(channel=ch).
    sig_ch = sig_ch.transpose("time", bin_dim)

    vertical_vals = np.asarray(vertical_scale_ch.values, dtype=float)

    if vertical_vals.size != sig_ch.sizes[bin_dim]:
        raise ValueError(
            "The selected vertical scale and signal have incompatible sizes: "
            f"vertical scale has {vertical_vals.size} values, but signal dimension "
            f"{bin_dim!r} has length {sig_ch.sizes[bin_dim]}."
        )

    # Prefer the explicit data_pack['bins'] entry for the lower x-axis. This is
    # small coordinate metadata and does not compute the lazy signal values.
    if bins_ch is not None:
        if len(bins_ch.dims) != 1:
            raise ValueError(
                "The selected bins entry must be 1D after selecting one channel. "
                f"Found dimensions: {bins_ch.dims}"
            )
        bin_vals = np.asarray(bins_ch.values, dtype=float)
    elif bin_dim in sig_ch.coords:
        bin_vals = np.asarray(sig_ch[bin_dim].values, dtype=float)
    else:
        bin_vals = np.arange(sig_ch.sizes[bin_dim], dtype=float)

    if bin_vals.size != sig_ch.sizes[bin_dim]:
        raise ValueError(
            "The selected bins entry and signal have incompatible sizes: "
            f"bins has {bin_vals.size} values, but signal dimension "
            f"{bin_dim!r} has length {sig_ch.sizes[bin_dim]}."
        )

    # Use a separate coordinate name to avoid collisions with the signal's
    # existing dimension coordinate named 'bins'.
    if x_name in sig_ch.dims or x_name in sig_ch.coords:
        x_name = f"{x_name}_value"

    sig_ch = sig_ch.assign_coords({x_name: (bin_dim, bin_vals)})

    return sig_ch, x_name, bin_vals, vertical_vals




def _clean_holoviews_title(title):
    """Make ATLAS titles safe and readable in HoloViews/Bokeh.

    HoloViews applies str.format() to titles, so literal braces must be
    escaped.  Bokeh titles do not render LaTeX, so common ATLAS LaTeX-like
    snippets are converted to unicode text.
    """

    if title is None:
        return title

    title = str(title)

    # Convert common ATLAS LaTeX-like title fragments to readable unicode.
    title = title.replace(r"$\nearrow$", "↗")
    title = title.replace(r"$\searrow$", "↘")
    title = title.replace(r"$\uparrow$", "↑")
    title = title.replace(r"$\downarrow$", "↓")
    title = title.replace(r"$^{o}$", "°")
    title = title.replace(r"$^{\circ}$", "°")
    title = title.replace(r"$^\circ$", "°")

    # Remove remaining dollar signs because Bokeh titles are plain text.
    title = title.replace("$", "")

    # Escape braces because HoloViews formats titles with str.format().
    title = title.replace("{", "{{").replace("}", "}}")

    return title


def _colors_from_colormap(cmap, n_colors):
    """Return n hex colors sampled from a matplotlib-like colormap.

    The module already imports the ATLAS quicklook colormap through
    get_quicklook_colormap().  This helper samples that same colormap so the
    multi-line profiles use the imported color scale, with early profiles on
    one side of the scale and late profiles on the other.
    """

    n_colors = int(n_colors)

    if n_colors <= 0:
        return []

    # Matplotlib colormap objects are callable.
    if callable(cmap):
        values = np.linspace(0.0, 1.0, n_colors)
        return [mpl_colors.to_hex(cmap(v)) for v in values]

    # Fallback for palettes/lists of colors.
    colors = list(cmap)

    if len(colors) == n_colors:
        return [mpl_colors.to_hex(c) for c in colors]

    idx = np.linspace(0, len(colors) - 1, n_colors).astype(int)

    return [mpl_colors.to_hex(colors[i]) for i in idx]


def _format_time_tick_label(value):
    """Return a robust, short label for a time-like coordinate value.

    Numpy datetime64 values, especially datetime64[ns], must be formatted
    before calling .item().  Otherwise numpy may convert them to large integer
    nanosecond timestamps, which then appear as long numbers on the colorbar.
    """

    # Handle numpy datetime arrays/scalars first.  This preserves datetime64[ns]
    # values as datetimes instead of converting them to integer nanoseconds.
    try:
        arr = np.asarray(value)
        if np.issubdtype(arr.dtype, np.datetime64):
            return np.datetime_as_string(arr.astype("datetime64[m]"), unit="m")
    except Exception:
        pass

    # Handle Python datetime/date objects next.
    try:
        import datetime as _datetime

        if isinstance(value, _datetime.datetime):
            return value.isoformat(timespec="minutes")
        if isinstance(value, _datetime.date):
            return value.isoformat()
    except Exception:
        pass

    # Only after the datetime checks, collapse zero-dimensional/object arrays.
    try:
        value = np.asarray(value).item()
    except Exception:
        pass

    try:
        import datetime as _datetime

        if isinstance(value, _datetime.datetime):
            return value.isoformat(timespec="minutes")
        if isinstance(value, _datetime.date):
            return value.isoformat()
    except Exception:
        pass

    return str(value)


def _make_time_colorbar_hook(time_values, colors, title="time"):
    """Create a Bokeh hook that adds a time colorbar to an hvPlot figure.

    The profile colors are categorical in the line overlay, but the hook adds
    a matching visual colorbar.  The colorbar labels show actual measurement
    times at a few representative positions.

    A singleton time dimension is treated as a single averaged profile, so no
    colorbar is added.  This avoids empty-looking plots and colorbar remnants
    when the profile_mean product has shape (time: 1, channel, bins).
    """

    n_time = len(time_values)

    if n_time <= 1:
        return None

    palette = list(colors)

    # Bokeh color mappers need at least two colors.
    if len(palette) == 1:
        palette = [palette[0], palette[0]]

    tick_count = min(5, n_time)
    tick_positions = np.unique(
        np.linspace(0, n_time - 1, tick_count).round().astype(int)
    )

    tick_labels = {
        int(i): _format_time_tick_label(time_values[int(i)])
        for i in tick_positions
    }

    def add_time_colorbar(plot, element):
        mapper = LinearColorMapper(
            palette=palette,
            low=0,
            high=max(n_time - 1, 1),
        )

        colorbar = ColorBar(
            color_mapper=mapper,
            ticker=FixedTicker(ticks=[int(i) for i in tick_positions]),
            major_label_overrides=tick_labels,
            label_standoff=8,
            title=title,
            width=12,
            title_text_font_size="16pt",
            major_label_text_font_size="14pt",
        )

        plot.state.add_layout(colorbar, "right")

    return add_time_colorbar



def _make_vertical_scale_axis_hook(
        bin_values,
        vertical_values,
        axis_label="vertical scale [km]",
        label_precision=4,
        ):
    """Create a Bokeh hook that adds vertical scale as a top x-axis.

    The plotted x-coordinate remains the dedicated bins entry on the lower
    x-axis. The top x-axis uses independent tick positions, chosen from the
    actual available bin samples inside the current viewport. Its labels are
    the corresponding actual vertical-scale values from vertical_values, not a
    rounded/interpolated approximation of the lower-axis ticks.

    Gridlines remain controlled by the lower x-axis only; the upper x-axis is
    just a reference scale.
    """

    bin_values = np.asarray(bin_values, dtype=float)
    vertical_values = np.asarray(vertical_values, dtype=float)

    mask = np.isfinite(bin_values) & np.isfinite(vertical_values)

    if np.count_nonzero(mask) < 1:
        return None

    bin_values = bin_values[mask]
    vertical_values = vertical_values[mask]

    order = np.argsort(bin_values)
    bin_values = bin_values[order]
    vertical_values = vertical_values[order]

    # Remove duplicate bin positions because the JavaScript ticker/formatter
    # expects a strictly increasing x array.
    unique_mask = np.ones(bin_values.size, dtype=bool)
    unique_mask[1:] = np.diff(bin_values) != 0
    bin_values = bin_values[unique_mask]
    vertical_values = vertical_values[unique_mask]

    if bin_values.size < 1:
        return None

    max_ticks = 7
    n_ticks = min(max_ticks, bin_values.size)

    if n_ticks <= 1:
        tick_indices = np.array([0], dtype=int)
    else:
        tick_indices = np.unique(
            np.linspace(0, bin_values.size - 1, n_ticks).round().astype(int)
        )

    top_tick_positions = [float(bin_values[i]) for i in tick_indices]

    def _format_vertical_tick(value):
        value = float(value)
        abs_value = abs(value)

        if abs_value >= 100:
            out = f"{value:.1f}"
        elif abs_value >= 10:
            out = f"{value:.2f}"
        elif abs_value >= 1:
            out = f"{value:.3f}"
        else:
            out = f"{value:.{label_precision}f}"

        # Drop trailing zeros while keeping plain decimal notation.
        if "." in out:
            out = out.rstrip("0").rstrip(".")

        return out

    top_tick_labels = {
        float(bin_values[i]): _format_vertical_tick(vertical_values[i])
        for i in tick_indices
    }

    def add_vertical_scale_axis(plot, element):
        top_axis = LinearAxis(
            axis_label=axis_label,
            ticker=FixedTicker(ticks=top_tick_positions),
            major_label_overrides=top_tick_labels,
            axis_label_text_font_size="16pt",
            major_label_text_font_size="14pt",
        )
        plot.state.add_layout(top_axis, "above")

    return add_vertical_scale_axis

def _make_line_plot_style_hook(add_zero_bin_strip=True):
    """Create a Bokeh hook for grid styling and the zero-bin grey strip.

    This keeps styling in a hook because some Bokeh objects, such as minor
    gridlines and annotations, are easiest to configure after hvPlot has created the
    underlying Bokeh figure.
    """

    def apply_line_plot_style(plot, element):
        fig = plot.state

        if add_zero_bin_strip:
            zero_strip = BoxAnnotation(
                left=0.0,
                right=1.0,
                fill_color="gray",
                fill_alpha=0.22,
                line_alpha=0.0,
            )
            fig.add_layout(zero_strip)

            zero_midline = Span(
                location=0.5,
                dimension="height",
                line_color="black",
                line_width=2,
                line_alpha=1.0,
            )
            fig.add_layout(zero_midline)

        # Major gridlines are already enabled through show_grid=True.  Add minor
        # gridlines as well.  Attribute availability is Bokeh-version dependent,
        # so use guarded assignments.
        for grid in list(fig.xgrid) + list(fig.ygrid):
            for attr, value in (
                ("grid_line_alpha", 0.45),
                ("grid_line_width", 1),
                ("minor_grid_line_color", "#d9d9d9"),
                ("minor_grid_line_alpha", 0.25),
                ("minor_grid_line_width", 1),
            ):
                try:
                    setattr(grid, attr, value)
                except Exception:
                    pass

    return apply_line_plot_style



def _deduplicate_toolbar_tools_hook(plot, element):
    """Remove duplicate Bokeh toolbar buttons after HoloViews builds the plot.

    hvPlot/HoloViews may add some default tools even when a custom tool list is
    supplied.  This hook keeps only the first tool of each type/dimension pair,
    so buttons such as pan, wheel zoom, box zoom, reset, or save appear once.
    """

    try:
        toolbar = plot.state.toolbar
    except Exception:
        return

    unique_tools = []
    seen = set()

    for tool in list(getattr(toolbar, "tools", [])):
        key = (
            tool.__class__.__name__,
            getattr(tool, "dimensions", None),
            getattr(tool, "mode", None),
        )

        if key in seen:
            continue

        seen.add(key)
        unique_tools.append(tool)

    toolbar.tools = unique_tools




def _ensure_plot_name(da, fallback="signal"):
    """Return a DataArray with a non-empty name for HoloViews/hvPlot.

    HoloViews can raise a DataError for 1D unnamed xarray DataArrays,
    especially after selecting a singleton time dimension with
    isel(time=0, drop=True).  Naming the DataArray does not compute data;
    it only updates metadata.
    """

    name = getattr(da, "name", None)

    if name is None or str(name).strip() == "":
        return da.rename(fallback)

    return da


def _cleanup_plotting_state():
    """Release plotting objects that HoloViews/Panel/Bokeh may keep referenced.

    Saving many standalone HTML files in one Python process can gradually grow
    memory because renderers keep references to the last rendered plot/model.
    This helper avoids accumulating those references between channels.
    """

    try:
        renderer = hv.Store.renderers.get("bokeh")
        if renderer is not None and hasattr(renderer, "last_plot"):
            renderer.last_plot = None
    except Exception:
        pass

    try:
        pn.state._views.clear()
    except Exception:
        pass

    gc.collect()


def _save_hvplot_html(plot, output_folder, filename):
    """Save an hvPlot/HoloViews object as a left-aligned HTML file.

    Wrapping the HoloViews object in a Panel Column gives better control over
    the exported page layout than saving the HoloViews object directly.  The
    plot itself is configured as responsive, while the Panel container is
    left-aligned and stretches to the available browser width.
    """

    os.makedirs(output_folder, exist_ok=True)

    root, _ = os.path.splitext(filename)
    fpath = os.path.join(output_folder, f"{root}.html")

    _ensure_panel_css()

    page = pn.Column(
        plot,
        sizing_mode="stretch_width",
        align="start",
        margin=(0, 0, 0, 0),
    )

    try:
        page.save(
            fpath,
            resources="inline",
            title=root,
        )
    finally:
        # Drop Panel/Bokeh references immediately after each file is saved.
        page.clear()
        del page
        _cleanup_plotting_state()

    return fpath
            
def generate_line_plots(processor, stage, db):
    
    allowed_dbs = {"profile", "profile_mean"}

    if db not in allowed_dbs:
        print()
        CustomWarning(
            "No signal-viewer plots were generated because db="
            f"{db!r} is not supported."
        )
        print()

        return
    
    data_pack = processor.package_from_stage(input_id = stage)
    
    caller_info = processor.processing_info['caller_info']

    missing_db_keys = [key for key, val in data_pack.items() if db not in val]

    if missing_db_keys:
        available_by_key = {
            key: sorted(list(data_pack[key].keys()))
            for key in missing_db_keys
        }

        print()
        CustomWarning(
            "No signal-viewer plots were generated because db="
            f"{db!r} does not exist for stage {stage!r} in "
            f"the following data_pack entries: {missing_db_keys}. "
        )
        print()
        return

    for key in data_pack.keys():

        print_header(f"Start plotting signals ({key})")
        
        # Prepare folders
        prepare_folder(caller_info, pattern = f"sig_mline_{key}")

        # Load arrays
        if caller_info['vertical_scale'] not in data_pack[key]:
            print()
            CustomWarning(
                "No signal-viewer plots were generated because vertical scale "
                f"{caller_info['vertical_scale']!r} does not exist for stage "
                f"{stage!r}, data_pack entry {key!r}."
            )
            print()

            return
            
        profiles = data_pack[key][db]
        vertical_scale = data_pack[key][caller_info['vertical_scale']]
        bins = data_pack[key].get("bins")
        system_info = data_pack[key]["system_info"]
        channel_info = data_pack[key]["channel_info"]

                    
        # Convert the range/height units to km 
        vertical_scale = convert_m_to_km(vertical_scale)

        # Check if the parsed channels exist and apply exclusion options
        channels = check_channels(
            all_channels = profiles.channel.values,
            settings = caller_info
            )
        

        sys_info = dict(zip(system_info.parameters.values, system_info.values))

        # iterate over the channels
        for ch in channels:
            
            print(f"-- channel: {ch}")
    
            ch_d = dict(channel = ch)
            
            ch_info = channel_info.sel({'channel':ch})  
            
            sig_ch = profiles.sel(ch_d)
            vertical_scale_ch = vertical_scale.sel(ch_d)
            bins_ch = _select_channel_if_present(bins, ch)

            # Gather the metadata that are common for all QA tests in a dictonary
            metadata = collect_metadata(data_pack[key], atlas_channel_id = ch)

            # Load libraris
            lib = Libraries(
                caller_info = caller_info,
                metadata = metadata,
                extra_metadata = {},
                settings = {},
                qa_test_info = {'stage': stage, 'qa_test':key, 'db': db}
                )
            
            
            # Call GenerateText class
            text_generator = GenerateText(lib = lib)
            
            # Make titles
            title = _clean_holoviews_title(
                text_generator.make_sig_mlines_title()
            )
            
            # Make filenames
            filename = text_generator.make_filename_viewer(
                qa_test = f'sig_mlines_{key}'
                )
            
        
            # Build an interactive multi-line plot for this channel.
            # This is the hvPlot equivalent of xarray.plot.line(x=..., hue='time'):
            # one line per time profile, with the vertical scale on the x-axis.
            sig_ch, x_name, bin_values, vertical_values = _prepare_multiline_profile(
                sig_ch=sig_ch,
                vertical_scale_ch=vertical_scale_ch,
                bins_ch=bins_ch,
                x_name="bins",
            )

            signal_name = sig_ch.name or "signal"
            sig_ch = _ensure_plot_name(sig_ch, fallback=signal_name)

            time_cmap = get_quicklook_colormap()
            n_time = sig_ch.sizes["time"]
            time_values = sig_ch["time"].values

            hooks = []

            vertical_axis_hook = _make_vertical_scale_axis_hook(
                bin_values=bin_values,
                vertical_values=vertical_values,
                axis_label=f"{caller_info['vertical_scale']} [km]",
            )

            if vertical_axis_hook is not None:
                hooks.append(vertical_axis_hook)

            hooks.append(_make_line_plot_style_hook(add_zero_bin_strip=True))
            hooks.append(_deduplicate_toolbar_tools_hook)

            # Do not calculate min/max here.  The signal data may be Dask-backed
            # and large; calculating limits would trigger an eager computation
            # before plotting.  HoloViews padding gives extra visible space around
            # both axes while keeping the selected signal lazy until rendering.
            common_plot_kwargs = dict(
                x=x_name,
                height=720,
                responsive=True,
                legend=False,
                hover=False,
                line_width=1.0,
                title=title,
                tools=["pan", "wheel_zoom", "box_zoom", "reset", "save"],
            )

            if n_time == 1:
                # A singleton time dimension represents one averaged profile.
                # Plot it as a normal single line.  Using by="time" for a
                # single implicit integer time index can produce an empty-looking
                # overlay in HoloViews/Bokeh, and a colorbar is not meaningful.
                single_colors = _colors_from_colormap(
                    cmap=time_cmap,
                    n_colors=3,
                )
                single_color = single_colors[1] if len(single_colors) > 1 else "black"
                sig_plot = sig_ch.isel(time=0, drop=True)
                sig_plot = _ensure_plot_name(sig_plot, fallback=signal_name)

                p = sig_plot.hvplot.line(
                    color=single_color,
                    **common_plot_kwargs,
                )
            else:
                time_colors = _colors_from_colormap(
                    cmap=time_cmap,
                    n_colors=n_time,
                )
                time_colorbar_hook = _make_time_colorbar_hook(
                    time_values=time_values,
                    colors=time_colors,
                    title="time",
                )

                if time_colorbar_hook is not None:
                    hooks.append(time_colorbar_hook)

                p = sig_ch.hvplot.line(
                    by="time",
                    color=hv.Cycle(time_colors),
                    **common_plot_kwargs,
                )

            p = p.opts(
                show_grid=True,
                min_width=900,
                min_height=600,
                default_tools=[],
                tools=["pan", "wheel_zoom", "box_zoom", "reset", "save"],
                active_tools=["box_zoom", "wheel_zoom"],
                hooks=hooks,
                padding=(0.03, 0.06),
                fontsize={
                    "title": 16,
                    "labels": 16,
                    "xticks": 14,
                    "yticks": 14,
                },
                xlabel="bin",
                ylabel=signal_name,
            )
    
            dir_out = os.path.join(caller_info["output_folder"],'signal_viewer',stage)
            os.makedirs(dir_out, exist_ok=True)
            fpath = _save_hvplot_html(p, dir_out, filename)

            # Explicitly remove per-channel references before moving to the next
            # channel.  This is important when the data_pack is large and lazy.
            del p
            del sig_ch
            del vertical_scale_ch
            del bins_ch
            del bin_values
            del vertical_values
            del vertical_axis_hook
            del metadata
            del lib
            del text_generator
            del time_values
            del signal_name
            del hooks
            if "time_colors" in locals():
                del time_colors
            if "time_colorbar_hook" in locals():
                del time_colorbar_hook
            if "sig_plot" in locals():
                del sig_plot
            _cleanup_plotting_state()
                        

        print('-----------------------------------------')
        print(' ')
