#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Telecover sector processor for the modern ATLAS visualizer architecture.

This module is intended to replace the old functional pattern::

    sector.process(...)

with a small reusable class::

    processor = TelecoverSectorProcessor(settings)
    result = processor.process(x=x_vals, y=sig_sector_values, iters=iters)

The returned dictionary is designed to be passed directly to the modern
telecover plotting/export functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional, Sequence, Tuple

import numpy as np


SmoothFunc = Callable[..., Tuple[np.ndarray, np.ndarray]]


@dataclass
class TelecoverSectorResult:
    """Container for one processed telecover sector."""

    coef: float
    y_m: np.ndarray
    y_sm: np.ndarray
    y_m_sm: np.ndarray
    y_l_sm: np.ndarray
    y_u_sm: np.ndarray
    coef_extra: Optional[float]
    y_extra: Optional[np.ndarray]
    y_extra_sm: Optional[np.ndarray]
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

    The class mirrors the old telecover-sector processing idea while returning
    named outputs instead of a long tuple.

    Parameters
    ----------
    settings : dict
        Telecover settings dictionary. Expected keys are mostly the same as in
        the old telecover code: ``smooth``, ``smoothing_window``,
        ``smooth_exponential`` and ``normalization_region``.
    smoothing_func : callable or None
        External smoothing function. If provided, it is called as::

            smoothing_func(
                y_vals=<2D array>,
                x_vals=<1D array>,
                x_sm_lims=<list>,
                x_sm_win=<float/int>,
                expo=<bool>,
            )

        and should return ``(y_smoothed, y_spread)``. If None, a simple internal
        moving-average fallback is used.
    """

    def __init__(
        self,
        settings: Dict[str, Any],
        smoothing_func: Optional[SmoothFunc] = None,
    ) -> None:
        self.settings = settings
        self.smoothing_func = smoothing_func

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
            Sector signal. Expected shape is ``(time, bins)``. A 1D array is
            accepted and internally treated as one profile.
        iters : int or None
            Number of profiles used for the main sector mean. If None, all
            profiles are used.
        x_sm_lims : list or tuple or None
            Smoothing range. If None, it is inferred from ``x``.
        region : list or tuple or None
            Normalization region. If None, ``settings['normalization_region']``
            is used.

        Returns
        -------
        dict
            Dictionary with keys compatible with the modern telecover plotting
            structure: ``coef``, ``y_m``, ``y_sm``, ``y_m_sm``, ``y_l_sm``,
            ``y_u_sm``, ``coef_extra``, ``y_extra``, ``y_extra_sm`` and
            ``has_extra``.
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

        # Main part: only the common number of iterations/profiles.
        y_main = y[:iters, :]

        y_sm = self._smooth_profiles(
            x=x,
            y=y_main,
            x_sm_lims=x_sm_lims,
        )

        y_m = np.nanmean(y_main, axis=0)
        y_m_sm = np.nanmean(y_sm, axis=0)
        y_l_sm, y_u_sm = self._profile_spread(y_sm)
        coef = self._normalization_coef(x=x, y=y_m_sm, region=region)

        # Extra profiles: sector-specific profiles beyond the common iteration count.
        has_extra = y.shape[0] > iters

        if has_extra:
            y_extra = np.nanmean(y[iters:, :], axis=0)
            y_extra_sm = self._smooth_1d(
                x=x,
                y=y_extra,
                x_sm_lims=x_sm_lims,
            )
            coef_extra = self._normalization_coef(x=x, y=y_extra_sm, region=region)
        else:
            y_extra = None
            y_extra_sm = None
            coef_extra = None

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

    def _smooth_profiles(
        self,
        x: np.ndarray,
        y: np.ndarray,
        x_sm_lims: Sequence[float],
    ) -> np.ndarray:
        """Smooth a 2D array profile-by-profile."""

        if not self.settings.get("smooth", False):
            return y.copy()

        if not self.settings.get("smoothing_window", None):
            return y.copy()

        if self.smoothing_func is not None:
            y_sm, _ = self.smoothing_func(
                y_vals=y,
                x_vals=x,
                x_sm_lims=x_sm_lims,
                x_sm_win=self.settings["smoothing_window"],
                expo=self.settings.get("smooth_exponential", False),
            )
            return np.asarray(y_sm, dtype=float)

        return np.vstack([
            self._moving_average_1d(row, int(self.settings["smoothing_window"]))
            for row in y
        ])

    def _smooth_1d(
        self,
        x: np.ndarray,
        y: np.ndarray,
        x_sm_lims: Sequence[float],
    ) -> np.ndarray:
        """Smooth one 1D profile."""

        if not self.settings.get("smooth", False):
            return y.copy()

        if not self.settings.get("smoothing_window", None):
            return y.copy()

        if self.smoothing_func is not None:
            y_sm, _ = self.smoothing_func(
                y_vals=y[np.newaxis, :],
                x_vals=x,
                x_sm_lims=x_sm_lims,
                x_sm_win=self.settings["smoothing_window"],
                expo=self.settings.get("smooth_exponential", False),
            )
            return np.asarray(y_sm, dtype=float)[0]

        return self._moving_average_1d(y, int(self.settings["smoothing_window"]))

    @staticmethod
    def _profile_spread(y_sm: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate lower/upper variability envelope from smoothed profiles.

        This uses mean +/- standard deviation. If your legacy code used a
        different convention, replace this method only.
        """

        y_m_sm = np.nanmean(y_sm, axis=0)
        y_std = np.nanstd(y_sm, axis=0)

        return y_m_sm - y_std, y_m_sm + y_std

    @staticmethod
    def _normalization_coef(
        x: np.ndarray,
        y: np.ndarray,
        region: Sequence[float],
    ) -> float:
        """
        Calculate normalization coefficient from the selected x-region.

        The coefficient is ``1 / mean(y in region)``. This matches the common
        pattern where normalized sector curves are plotted as ``coef * y``.
        """

        mask = (x >= region[0]) & (x <= region[1]) & np.isfinite(x) & np.isfinite(y)

        if not np.any(mask):
            return 1.0

        norm_value = np.nanmean(y[mask])

        if not np.isfinite(norm_value) or norm_value == 0:
            return 1.0

        return float(1.0 / norm_value)

    @staticmethod
    def _moving_average_1d(y: np.ndarray, window: int) -> np.ndarray:
        """Simple NaN-aware moving average fallback."""

        y = np.asarray(y, dtype=float)
        window = max(1, int(window))

        if window == 1:
            return y.copy()

        kernel = np.ones(window, dtype=float)

        valid = np.isfinite(y).astype(float)
        y_safe = np.where(np.isfinite(y), y, 0.0)

        numerator = np.convolve(y_safe, kernel, mode="same")
        denominator = np.convolve(valid, kernel, mode="same")

        with np.errstate(invalid="ignore", divide="ignore"):
            y_sm = numerator / denominator

        y_sm[denominator == 0] = np.nan

        return y_sm
