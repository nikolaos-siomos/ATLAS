#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 15:36:14 2026

@author: nikos
"""

import numpy as np


def get_x_ticks(x_lims, x_tick):
    """
    Generate major x ticks and labels for a vertical-axis plot.
    """

    x_llim = x_lims[0]
    x_ulim = x_lims[1]

    if x_tick is None:
        width = x_ulim - x_llim
        if width <= 3:
            x_tick = 0.5
        elif width <= 8:
            x_tick = 1.0
        else:
            x_tick = 2.0

    if x_tick >= x_ulim - x_llim:
        raise ValueError(
            f"The x_tick ({x_tick}) must be smaller than the x-axis width "
            f"({x_llim} to {x_ulim} km)."
        )

    ticks = np.arange(
        x_tick * np.ceil(x_llim / x_tick),
        x_tick * (np.floor(x_ulim / x_tick) + 1.0),
        x_tick,
    )

    if ticks.size == 0:
        ticks = np.array([x_llim, x_ulim])

    if np.abs(x_llim - ticks[0]) < x_tick * 0.25:
        ticks[0] = x_llim
    else:
        ticks = np.hstack((x_llim, ticks))

    if np.abs(x_ulim - ticks[-1]) < x_tick * 0.25:
        ticks[-1] = x_ulim
    else:
        ticks = np.hstack((ticks, x_ulim))

    ticks = np.round(ticks, decimals=2)
    labels = [str(float(tick)).rstrip("0").rstrip(".") for tick in ticks]

    return ticks, labels, x_tick


def _as_finite_array(values):
    """
    Convert scalars/arrays/lists to one flat finite float array.
    """
    finite_values = []

    if values is None:
        return np.array([], dtype=float)

    if not isinstance(values, (list, tuple)):
        values = [values]

    for value in values:
        if value is None:
            continue

        try:
            arr = np.asarray(value, dtype=float).ravel()
        except Exception:
            continue

        arr = arr[np.isfinite(arr)]
        if arr.size > 0:
            finite_values.append(arr)

    if len(finite_values) == 0:
        return np.array([], dtype=float)

    return np.concatenate(finite_values)


def _limit_value(y_lims, index):
    """
    Return one user-provided y-limit value, allowing None/empty/partial limits.
    """
    if y_lims is None:
        return None

    if len(y_lims) <= index:
        return None

    return y_lims[index]


def get_calibration_y_limits(y_values, y_lims):
    """
    Determine y-limits for the calibration panel from scalar mean values.

    This follows the legacy polarization_calibration_cal_y logic.  The caller
    should pass mean gain-ratio diagnostics, e.g.
    [gain_ratio_m45_mean, gain_ratio_p45_mean, gain_ratio_mean].
    """

    vals = _as_finite_array(y_values)

    y_min = np.nanmin(vals) if vals.size else np.nan
    y_max = np.nanmax(vals) if vals.size else np.nan

    lower_user = _limit_value(y_lims, 0)
    upper_user = _limit_value(y_lims, 1)

    if np.isfinite(y_max) and upper_user is None:
        y_ulim = y_max * 2.0
    elif upper_user is not None:
        y_ulim = upper_user
    else:
        y_ulim = 1.0

    if np.isfinite(y_min) and lower_user is None:
        y_llim = y_min / 1.5
    elif lower_user is not None:
        y_llim = lower_user
    else:
        y_llim = 0.0

    if y_llim == y_ulim:
        margin = 0.1 * abs(y_ulim) if y_ulim != 0 else 0.1
        y_llim -= margin
        y_ulim += margin

    return [y_llim, y_ulim]


def get_rayleigh_y_limits(y_values, y_lims):
    """
    Determine y-limits for the Rayleigh/VLDR panel from scalar mean values.

    This follows the legacy polarization_calibration_ray_y logic.  The caller
    should pass scalar diagnostics, e.g.
    [calibrated_ratio_mean, vldr_mean, mldr_mean].
    """

    vals = _as_finite_array(y_values)
    y_max = np.nanmax(vals) if vals.size else np.nan

    lower_user = _limit_value(y_lims, 0)
    upper_user = _limit_value(y_lims, 1)

    if np.isfinite(y_max) and upper_user is None:
        y_ulim = y_max * 2.5
    elif upper_user is not None:
        y_ulim = upper_user
    else:
        y_ulim = 0.1

    if lower_user is None:
        y_llim = 0.0
    else:
        y_llim = lower_user

    if y_llim == y_ulim:
        margin = 0.1 * abs(y_ulim) if y_ulim != 0 else 0.1
        y_llim -= margin
        y_ulim += margin

    return [y_llim, y_ulim]
