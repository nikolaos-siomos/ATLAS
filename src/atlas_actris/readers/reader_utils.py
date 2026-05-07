#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 20:55:22 2025

@author: nikos
"""

import pandas as pd
from typing import Union, Optional
import pandas as pd
import numpy as np
from typing import Union, Tuple, Optional
from pandas.api.types import is_numeric_dtype

PDObj = Union[pd.Series, pd.DataFrame]

def compare_common_with_mismatches(
    a: PDObj,
    b: PDObj,
    *,
    atol: float = 0.0,
    rtol: float = 0.0,
) -> Tuple[bool, pd.DataFrame]:
    """
    Compare two pandas objects (both Series or both DataFrames) on their
    common index/columns. Return (all_equal, mismatches_df).

    - NaNs in the same position are considered equal.
    - If both sides are numeric, differences within atol + rtol*|b| are considered equal.
    - If there's no overlap, returns (True, empty DataFrame).
    """

    if type(a) is not type(b):
        raise TypeError("Both inputs must be the same type (Series vs Series, or DataFrame vs DataFrame).")

    # --- Series case ---
    if isinstance(a, pd.Series):
        a2, b2 = a.align(b, join="inner")
        if a2.size == 0:
            return True, pd.DataFrame(columns=["left", "right"])

        # Base equality + NaN-equality
        eq = a2.eq(b2) | (a2.isna() & b2.isna())

        # Numeric tolerance (only if both are numeric dtypes)
        if is_numeric_dtype(a2) and is_numeric_dtype(b2) and (atol > 0 or rtol > 0):
            an = a2.to_numpy(dtype=float, copy=False)
            bn = b2.to_numpy(dtype=float, copy=False)
            both_num = ~(np.isnan(an) | np.isnan(bn))
            close = np.zeros_like(both_num, dtype=bool)
            close[both_num] = np.isclose(an[both_num], bn[both_num], atol=atol, rtol=rtol)
            eq = eq | pd.Series(close, index=a2.index)

        mism = ~eq
        if not mism.any():
            return True, pd.DataFrame(columns=["left", "right"])

        out = pd.DataFrame({"left": a2[mism], "right": b2[mism]})
        return False, out

    # --- DataFrame case ---
    else:
        a2, b2 = a.align(b, join="inner", axis=None)
        if a2.size == 0:  # no overlapping cells
            return True, pd.DataFrame(columns=["index", "column", "left", "right"])

        # Start with elementwise equality and NaN-equality
        eq = a2.eq(b2) | (a2.isna() & b2.isna())

        # Apply numeric tolerance per column where both sides are numeric
        if atol > 0 or rtol > 0:
            for col in a2.columns.intersection(b2.columns):
                col_a = a2[col]
                col_b = b2[col]
                if is_numeric_dtype(col_a) and is_numeric_dtype(col_b):
                    an = col_a.to_numpy(dtype=float, copy=False)
                    bn = col_b.to_numpy(dtype=float, copy=False)
                    both_num = ~(np.isnan(an) | np.isnan(bn))
                    close = np.zeros_like(both_num, dtype=bool)
                    close[both_num] = np.isclose(an[both_num], bn[both_num], atol=atol, rtol=rtol)
                    # merge tolerance result into eq for this column
                    eq[col] = eq[col] | pd.Series(close, index=a2.index)

        mism_mask = ~eq
        if not mism_mask.any().any():
            return True, pd.DataFrame(columns=["index", "column", "left", "right"])

        # Build a tidy table of mismatches
        where = np.where(mism_mask.to_numpy())
        rows = []
        idx_vals = a2.index.to_numpy()
        col_vals = a2.columns.to_numpy()
        for i, j in zip(*where):
            rows.append({
                "index": idx_vals[i],
                "column": col_vals[j],
                "left": a2.iat[i, j],
                "right": b2.iat[i, j],
            })
        out = pd.DataFrame(rows).set_index(["index", "column"])
        return False, out
    
def transfer_missing_metadata(
    dst: PDObj,
    src: Optional[PDObj] = None,
    *,
    inplace: bool = True,
    treat_falsy_as_missing: bool = False,
) -> PDObj:
    """
    Copy only metadata keys that are absent on `dst` from `src` into `dst.attrs`.

    - Works for DataFrame or Series.
    - Does NOT overwrite existing keys on `dst`.
    - If `src` is None, or has no attrs, it's a no-op.
    - If `treat_falsy_as_missing=True`, keys on `dst` with values in {None, "", []}
      are treated as missing and will be filled from `src` when available.
    """
    if src is None:
        return dst if inplace else dst.copy(deep=False)

    if not hasattr(dst, "attrs") or not hasattr(src, "attrs"):
        raise TypeError("Both dst and src must be a pandas DataFrame or Series with .attrs.")

    out = dst if inplace else dst.copy(deep=False)

    if not src.attrs:
        return out  # handles 'empty' src object case

    if treat_falsy_as_missing:
        def is_missing(v):
            return v is None or v == "" or v == []
        missing_keys = {k for k, v in out.attrs.items() if is_missing(v)}
        for k, v in src.attrs.items():
            if k not in out.attrs or k in missing_keys:
                out.attrs[k] = v
    else:
        for k, v in src.attrs.items():
            if k not in out.attrs:
                out.attrs[k] = v

    return out

def special_path_rules(caller_info):

    if caller_info["raw_file_format"] == "scc":
    
        infered = ["abs_drk_ray", "abs_drk_ray_pcb", "abs_drk_dtm", "abs_drk_trg"]
        for key in infered:
            key_m = key.replace("abs_drk_", "abs_", 1)
            if caller_info[key] == None and caller_info[key_m] != None:
                caller_info[key] = caller_info[key_m]
        
        pairs = {"abs_drk_tlc": "abs_tlc_north",
                 "abs_drk_tlc_rin": "abs_tlc_rin_inner",
                 "abs_drk_pcb": "abs_pcb_p45",
                 "abs_drk_pcb_aux": "abs_pcb_aux_p45"}
        for key, paired in pairs.items():
            if caller_info[key] == None and caller_info[paired] != None:
                caller_info[key] = caller_info[paired]
        
    if caller_info["raw_file_format"] == "polly_xt":
        
        pairs = {"abs_pcb_p45": "abs_ray",
                 "abs_pcb_m45": "abs_ray"}
        for key, paired in pairs.items():
            if caller_info[key] == None and caller_info[paired] != None:
                caller_info[key] = caller_info[paired]
    
    return(caller_info)