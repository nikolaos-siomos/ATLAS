#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Telecover sector processor for the modern ATLAS visualizer architecture.

This class mirrors the legacy visualizer.sector.process(...) procedure while
returning named outputs instead of a long tuple.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

from visualizer import normalize
from visualizer.smoothing import (
    sliding_average_1D,
    sliding_average_1D_fast,
    sliding_average_2D,
    sliding_average_2D_fast,
)


@dataclass
class TelecoverSectorResult:
    """Container for one processed telecover sector."""

    coef: Any
    y_m: np.ndarray
    y_sm: np.ndarray
    y_m_sm: np.ndarray
    y_l_sm: np.ndarray
    y_u_sm: np.ndarray
    coef_extra: Any
    y_extra: Any
    y_extra_sm: Any
    has_extra: bool

    def as_dict(self) -> Dict[str, Any]:
        """Return the result as a plain dictionary."""
        return {
            "coef": self.coef,
            "y_m": self.y_m,
            "y_sm": self.y_sm,
            "y_m_sm": self.y_m_sm,
            "y_l_sm": self.y_l_sm,
            "y_u_sm": self.y_u_sm,
            "coef_extra": self.coef_extra,
            "y_extra": self.y_extra,
            "y_extra_sm": self.y_extra_sm,
            "has_extra": self.has_extra,
        }


class TelecoverSectorProcessor:
    """
    Process one telecover sector.

    This mirrors the old sector.process(...) procedure:
    - average the first ``iters`` profiles;
    - optionally smooth the averaged and unaveraged profiles;
    - use near smoothing inside ``x_sm_lims``;
    - use far smoothing from ``x_sm_lims[1]`` to 20 km;
    - keep the original signal outside the smoothed ranges;
    - use min/max of the smoothed individual profiles as variability envelope;
    - normalize with ``visualizer.normalize.to_a_point``;
    - keep one extra profile if the sector has exactly ``iters + 1`` profiles.
    """

    def __init__(self, settings: Dict[str, Any]) -> None:
        self.settings = settings

    def process(
        self,
        x: Sequence[float],
        y: np.ndarray,
        iters: Optional[int] = None,
        x_sm_lims: Optional[Sequence[float]] = None,
        region: Optional[Sequence[float]] = None,
    ) -> Dict[str, Any]:
        """
        Process one sector signal array.

        Parameters
        ----------
        x : array-like
            Vertical coordinate, already sliced and converted to km.
        y : array-like
            Sector signal. Expected shape is (time, bins). A 1D array is
            accepted and internally treated as one profile.
        iters : int or None
            Number of profiles used for the common sector mean. If None, all
            available profiles are used.
        x_sm_lims : list or tuple or None
            Smoothing range. If None, it is inferred from x.
        region : list or tuple or None
            Normalization region. If None, settings['normalization_region'] is used.

        Returns
        -------
        dict
            Dictionary equivalent to the legacy sector.process tuple outputs.
        """

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if y.ndim == 1:
            y = y[np.newaxis, :]

        if y.ndim != 2:
            raise ValueError("y must be either 1D or 2D with shape (time, bins).")

        if y.shape[1] != x.size:
            raise ValueError(
                f"The bins dimension of y ({y.shape[1]}) must match x size ({x.size})."
            )

        if iters is None:
            iters = y.shape[0]

        iters = int(min(iters, y.shape[0]))

        if iters < 1:
            raise ValueError("iters must be at least 1 after clipping to available profiles.")

        if x_sm_lims is None:
            x_sm_lims = [float(np.nanmin(x)), float(np.nanmax(x))]

        if region is None:
            region = self.settings["normalization_region"]

        smooth = bool(self.settings.get("smooth", False))
        x_sm_win = self.settings.get("smoothing_window")
        expo = bool(self.settings.get("smooth_exponential", False))

        # Legacy: average only the common number of profiles.
        y_m = np.nanmean(y[:iters, :], axis=0)

        if smooth:
            y_m_sm, y_sm = self._legacy_smoothing(
                x=x,
                y=y,
                y_m=y_m,
                x_sm_lims=x_sm_lims,
                x_sm_win=x_sm_win,
                expo=expo,
            )
        else:
            y_m_sm = y_m
            y_sm = y

        # Legacy: extra profile exists only if there is exactly one profile
        # beyond the common iteration count.
        if y_sm.shape[0] == iters + 1:
            y_extra = y[iters, :]
            y_extra_sm = y_sm[iters, :]
            has_extra = True
        else:
            y_extra = []
            y_extra_sm = []
            coef_extra = []
            has_extra = False

        # Legacy: variability envelope is min/max, not mean +/- std.
        y_l_sm = np.nanmin(y_sm[:iters, :], axis=0)
        y_u_sm = np.nanmax(y_sm[:iters, :], axis=0)

        # Legacy: normalization coefficient is calculated using the raw mean
        # signal y_m, not the smoothed mean y_m_sm.
        coef, _ = normalize.to_a_point(
            sig=y_m,
            sig_b=np.ones(x.shape),
            x_vals=x,
            region=region,
            axis=0,
        )

        if has_extra:
            coef_extra, _ = normalize.to_a_point(
                sig=y_extra,
                sig_b=np.ones(x.shape),
                x_vals=x,
                region=region,
                axis=0,
            )

        return TelecoverSectorResult(
            coef=coef,
            y_m=y_m,
            y_sm=y_sm,
            y_m_sm=y_m_sm,
            y_l_sm=y_l_sm,
            y_u_sm=y_u_sm,
            coef_extra=coef_extra,
            y_extra=y_extra,
            y_extra_sm=y_extra_sm,
            has_extra=has_extra,
        ).as_dict()

    @staticmethod
    def _legacy_smoothing(
        x: np.ndarray,
        y: np.ndarray,
        y_m: np.ndarray,
        x_sm_lims: Sequence[float],
        x_sm_win: Any,
        expo: bool,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Apply the same near/far smoothing strategy as the legacy sector.py.
        """

        # Legacy switches between variable-window and fast constant-window
        # smoothing based on whether the smoothing window is a list.
        if isinstance(x_sm_win, list):
            smooth_1d_near = sliding_average_1D
            smooth_2d_near = sliding_average_2D
        else:
            smooth_1d_near = sliding_average_1D_fast
            smooth_2d_near = sliding_average_2D_fast

        # Near smoothing: user-defined smoothing region/window.
        y_m_sm_n, _ = smooth_1d_near(
            y_vals=y_m,
            x_vals=x,
            x_sm_lims=x_sm_lims,
            x_sm_win=x_sm_win,
            expo=expo,
        )

        y_sm_n, _ = smooth_2d_near(
            z_vals=y,
            y_vals=x,
            y_sm_lims=x_sm_lims,
            y_sm_win=x_sm_win,
            expo=expo,
        )

        # Far smoothing: legacy always uses the fast functions from the upper
        # smoothing limit to 20 km with a 500 m window.
        y_m_sm_f, _ = sliding_average_1D_fast(
            y_vals=y_m,
            x_vals=x,
            x_sm_lims=[x_sm_lims[1], 20.0],
            x_sm_win=500.0,
            expo=False,
        )

        y_sm_f, _ = sliding_average_2D_fast(
            z_vals=y,
            y_vals=x,
            y_sm_lims=[x_sm_lims[1], 20.0],
            y_sm_win=500.0,
            expo=False,
        )

        y_m_sm = np.nan * np.zeros(y_m.copy().shape)

        mask_n = (x > x_sm_lims[0]) & (x < x_sm_lims[1])
        mask_f = (x >= x_sm_lims[1]) & (x < 20.0)

        y_m_sm[mask_n] = y_m_sm_n[mask_n]
        y_m_sm[mask_f] = y_m_sm_f[mask_f]

        y_sm = y.copy()
        y_sm[:, mask_n] = y_sm_n[:, mask_n]
        y_sm[:, mask_f] = y_sm_f[:, mask_f]

        return y_m_sm, y_sm
