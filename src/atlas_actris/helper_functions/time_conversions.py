#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 20:30:40 2025

@author: nikos
"""


from typing import Literal, Sequence, List
import pandas as pd

# ---------- 1) Anything -> ISO strings ----------
def datetimes_to_iso(
    arr: Sequence,
    tz: Literal["utc", "naive"] = "utc",
) -> List[str]:
    """
    Convert an array-like of datetime-ish values (datetime, numpy.datetime64,
    pandas.Timestamp, mixed) to ISO-8601 strings.

    tz="utc"  -> normalize to UTC and append 'Z' (recommended, consistent)
    tz="naive"-> drop timezone info (no 'Z')
    """
    # Parse everything; utc=True ensures consistent timezone handling
    s = pd.to_datetime(arr, utc=True)  # vectorized, tolerant to mixed inputs

    if tz == "naive":
        s = s.tz_localize(None)  # drop tz
        iso = (
            s.strftime("%Y-%m-%dT%H:%M:%S.%f")
             .str.rstrip("0").str.rstrip(".")
        ).tolist()
    else:  # tz == "utc"
        # Keep UTC and mark with 'Z'
        iso = (
            s.strftime("%Y-%m-%dT%H:%M:%S.%f")
             .str.rstrip("0").str.rstrip(".") + "Z"
        ).tolist()

    return iso


# ---------- 2) ISO strings -> desired type ----------
def iso_to_datetimes(
    iso_arr: Sequence[str],
    out: Literal["pandas", "numpy", "python"] = "pandas",
    tz: Literal["utc", "naive"] = "utc",
):
    """
    Convert ISO-8601 strings back to:
      - out="pandas": pandas.DatetimeIndex
      - out="numpy" : numpy.ndarray dtype datetime64[ns]
      - out="python": list[datetime.datetime]

    tz="utc"  -> keep/convert to UTC
    tz="naive"-> drop tz info
    """
    s = pd.to_datetime(iso_arr, utc=True)  # parse and normalize to UTC first

    if tz == "naive":
        s_naive = s.tz_localize(None)
        if out == "pandas":
            return pd.DatetimeIndex(s_naive)
        if out == "numpy":
            return s_naive.to_numpy(dtype="datetime64[ns]")
        if out == "python":
            return [ts.to_pydatetime() for ts in s_naive]

    else:  # tz == "utc" (tz-aware UTC)
        if out == "pandas":
            return pd.DatetimeIndex(s)
        if out == "numpy":
            # NumPy datetimes are tz-naive; convert to UTC wall time and drop tz
            return s.tz_convert("UTC").tz_localize(None).to_numpy(dtype="datetime64[ns]")
        if out == "python":
            return [ts.to_pydatetime() for ts in s.tz_convert("UTC")]

    raise ValueError("Invalid 'out' or 'tz' parameter.")