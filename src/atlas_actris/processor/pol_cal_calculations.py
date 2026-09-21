#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr

from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from utils.dataarray_utils import shallow_copy
from utils.printouts import print_entry
from utils.error_classes import CustomWarning


def make_ratio_id(ch_r, ch_t, type_index):
    """
    Build an 8-character ratio_id from a reflected/transmitted channel pair.

    Rule
    ----
    Digits 1-4:
        same wavelength as ch_r/ch_t

    Digit 5:
        same as ch_r/ch_t if equal, otherwise "u"

    Digit 6:
        type_index

    Digit 7:
        same as ch_r/ch_t if equal, otherwise "u"

    Digit 8:
        "x"

    The current polarization-calibration products use these type indices:
        g: gain_ratio
        e: eta
        c: calibrated_ratio
        v: vldr
        d: uncalibrated_ratio, kept for backwards compatibility
        s/f: intermediate eta_s / eta_s_f, not normally saved
    """

    allowed_indices = ["g", "e", "c", "v", "m", "d", "s", "f"]

    if type_index not in allowed_indices:
        raise ValueError(
            f"type_index must be one of: {allowed_indices}"
        )

    if len(ch_r) != 8 or len(ch_t) != 8:
        raise ValueError(
            f"Cannot build ratio_id from non-8-character channels: "
            f"{ch_r}, {ch_t}"
        )

    if ch_r[:4] != ch_t[:4]:
        raise ValueError(
            f"Cannot build ratio_id because wavelengths differ: "
            f"{ch_r[:4]} != {ch_t[:4]}"
        )

    digit_5 = ch_r[4] if ch_r[4] == ch_t[4] else "u"
    digit_7 = ch_r[6] if ch_r[6] == ch_t[6] else "u"

    return f"{ch_r[:4]}{digit_5}{type_index}{digit_7}x"


def replace_ratio_id_type(pair_id: str, old_type_index: str, new_type_index: str) -> str:
    """Return pair_id with digit 6, index 5, replaced by new_type_index.

    The helper intentionally preserves all other characters of an existing
    derived pair id so downstream products remain traceable to the product
    they were derived from.
    """

    pair_id = str(pair_id)

    if len(pair_id) != 8:
        raise ValueError(f"Cannot update non-8-character pair id: {pair_id}")

    if pair_id[5] != old_type_index:
        raise ValueError(
            f"Cannot update pair id {pair_id}: expected type index "
            f"'{old_type_index}' at position 5, found '{pair_id[5]}'."
        )

    return f"{pair_id[:5]}{new_type_index}{pair_id[6:]}"


def mean_in_region(
    da: xr.DataArray,
    z: xr.DataArray,
    averaging_range,
) -> xr.DataArray:
    """Average a DataArray over bins inside the averaging_range."""

    z_min, z_max = averaging_range
    mask = (z >= 1E3 * z_min) & (z <= 1E3 * z_max)
    return da.where(mask).mean("bins", skipna=True)


def sem_in_region(
    da: xr.DataArray,
    z: xr.DataArray,
    averaging_range,
) -> xr.DataArray:
    """
    Average an error DataArray over bins and convert it to SEM.

    This assumes da_error is already the per-bin uncertainty of the ratio.
    """

    z_min, z_max = averaging_range
    mask = (z >= 1E3 * z_min) & (z <= 1E3 * z_max)

    da_sel = da.where(mask)
    n_bins = da_sel.notnull().sum("bins")

    return da_sel.std("bins", skipna=True) / np.sqrt(n_bins)

# def sem_in_region(
#     da_error: xr.DataArray,
#     z: xr.DataArray,
#     averaging_range,
# ) -> xr.DataArray:
#     """
#     Average an error DataArray over bins and convert it to SEM.

#     This assumes da_error is already the per-bin uncertainty of the ratio.
#     """

#     z_min, z_max = averaging_range
#     mask = (z >= 1E3 * z_min) & (z <= 1E3 * z_max)

#     da_sel = da_error.where(mask)
#     n_bins = da_sel.notnull().sum("bins")

#     return da_sel.mean("bins", skipna=True) / np.sqrt(n_bins)


def channels_to_pairs(
    da: xr.DataArray,
    info: xr.DataArray,
    channel_dim: str = "channel",
    pair_dim: str = "pair",
) -> xr.DataArray:
    """
    Convert a channel-based DataArray to a pair-based DataArray by position.

    Works for arrays such as:
        da(time, channel, bins)
        da(channel, bins)
    """

    if channel_dim not in da.dims:
        raise ValueError(
            f"Input DataArray does not have dimension '{channel_dim}'. "
            f"Found dims: {da.dims}"
        )

    if pair_dim not in info.dims:
        raise ValueError(
            f"Info DataArray does not have dimension '{pair_dim}'. "
            f"Found dims: {info.dims}"
        )

    if da.sizes[channel_dim] != info.sizes[pair_dim]:
        raise ValueError(
            f"Number of channels must match number of pairs. "
            f"Got {da.sizes[channel_dim]} channels and "
            f"{info.sizes[pair_dim]} pairs."
        )

    pair_index = list(range(info.sizes[pair_dim]))

    da_pair = (
        da
        .assign_coords({channel_dim: pair_index})
        .rename({channel_dim: pair_dim})
        .assign_coords({pair_dim: info[pair_dim]})
    )

    preferred_order = [
        dim for dim in ["time", pair_dim, "bins"]
        if dim in da_pair.dims
    ]
    other_dims = [dim for dim in da_pair.dims if dim not in preferred_order]

    return da_pair.transpose(*preferred_order, *other_dims)


def simple_ratio(
    numerator: xr.DataArray,
    denominator: xr.DataArray,
    info: xr.DataArray,
) -> xr.DataArray:
    """Calculate numerator / denominator for paired channels."""

    numerator = channels_to_pairs(numerator, info)
    denominator = channels_to_pairs(denominator, info)

    ratio_da = numerator / denominator
    ratio_da = ratio_da.where(denominator != 0)
    ratio_da.name = "ratio"

    return ratio_da

def ratio_error_independent(
    numerator: xr.DataArray,
    denominator: xr.DataArray,
    numerator_error: xr.DataArray,
    denominator_error: xr.DataArray,
    ratio_values: xr.DataArray,
    info: xr.DataArray,
) -> xr.DataArray:
    """Independent error propagation for ratio = numerator / denominator."""

    numerator = channels_to_pairs(numerator, info)
    denominator = channels_to_pairs(denominator, info)
    numerator_error = channels_to_pairs(numerator_error, info)
    denominator_error = channels_to_pairs(denominator_error, info)

    rel_num = numerator_error / numerator
    rel_den = denominator_error / denominator

    ratio_error = np.abs(ratio_values) * np.sqrt(rel_num**2 + rel_den**2)
    ratio_error = ratio_error.where(denominator != 0)
    ratio_error.name = "ratio_error"

    return ratio_error


def as_ratio_dataarray(values: xr.DataArray, ratio_id: str) -> xr.DataArray:
    """Convert values(time, bins) to values(time, pair, bins)."""

    if "channel" in values.dims:
        values = values.squeeze("channel", drop=True)

    ratio_da = values.expand_dims(pair=[ratio_id])
    preferred = [dim for dim in ["time", "pair", "bins"] if dim in ratio_da.dims]
    other = [dim for dim in ratio_da.dims if dim not in preferred]
    ratio_da = ratio_da.transpose(*preferred, *other)
    ratio_da.name = "ratio"

    return ratio_da


def _as_pair_dataarray(
    values: Any,
    pairs: xr.DataArray,
    name: str,
) -> xr.DataArray:
    """Normalize scalars/sequences/DataArrays to a one-dimensional pair DataArray."""

    if isinstance(values, xr.DataArray):
        da = values
        if "pair" not in da.dims:
            if da.size == len(pairs):
                dim = da.dims[0] if da.dims else "pair"
                da = da.rename({dim: "pair"}) if dim != "pair" else da
            else:
                da = xr.DataArray(
                    np.repeat(da.item(), len(pairs)),
                    dims=["pair"],
                    coords={"pair": pairs},
                    name=name,
                )
        return da.assign_coords(pair=pairs).rename(name)

    if isinstance(values, (str, bytes)) or np.isscalar(values):
        data = [values] * len(pairs)
    else:
        data = list(values)
        if len(data) != len(pairs):
            raise ValueError(
                f"Parameter '{name}' has {len(data)} values, but {len(pairs)} pairs."
            )

    return xr.DataArray(
        data,
        dims=["pair"],
        coords={"pair": pairs},
        name=name,
    )


def add_parameter(
    pci: xr.DataArray,
    name: str,
    values: Any,
) -> xr.DataArray:
    """
    Add or replace one parameter row in pol_cal_info.

    pol_cal_info is metadata: it is not chunked and it has mixed/object dtype.
    Therefore values must be materialized before concatenation. Keeping a
    Dask-backed object row here triggers Dask's object-dtype auto-rechunking
    limitation. Full profile-like products remain lazy outside pol_cal_info.
    """

    values = _as_pair_dataarray(values, pci.pair, name)

    if "time" in values.dims:
        raise ValueError(
            f"Cannot store '{name}' in pol_cal_info because it still has a time dimension."
        )

    if "pair" not in values.dims:
        raise ValueError(
            f"Cannot store '{name}' in pol_cal_info because it has no pair dimension."
        )

    # Materialize only this small pair-level metadata row.
    # This is intentional: pol_cal_info is object metadata, not a lazy profile.
    values = values.compute()

    # Make sure the existing metadata table is also eager before object concat.
    # If it is already NumPy-backed, this is effectively a no-op.
    pci = pci.compute()

    new_row = (
        values
        .assign_coords(pair=pci.pair)
        .expand_dims(parameters=[name])
        .transpose("parameters", "pair")
        .astype(object)
        .compute()
    )

    if name in pci.parameters.values:
        pci = pci.drop_sel(parameters=name)

    return xr.concat(
        [pci.astype(object), new_row],
        dim="parameters",
        join="outer",
        combine_attrs="override",
    )

def append_or_replace_pairs(
    existing: Optional[xr.DataArray],
    new: xr.DataArray,
) -> xr.DataArray:
    """Add pair entries to an existing pair-based DataArray, replacing duplicates."""

    if existing is None:
        return new

    if "pair" not in new.dims:
        raise ValueError(f"New DataArray has no pair dimension: {new.dims}")

    if "pair" not in existing.dims:
        # Backwards compatibility for old files created with a 'ratio' dimension.
        if "ratio" in existing.dims:
            existing = existing.rename({"ratio": "pair"})
        else:
            raise ValueError(f"Existing DataArray has no pair dimension: {existing.dims}")

    drop_ids = [pid for pid in new.pair.values if pid in existing.pair.values]
    if drop_ids:
        existing = existing.drop_sel(pair=drop_ids)

    if existing.sizes.get("pair", 0) == 0:
        combined = new
    else:
        combined = xr.concat([existing, new], dim="pair", join="outer")

    preferred = [dim for dim in ["time", "pair", "bins"] if dim in combined.dims]
    other = [dim for dim in combined.dims if dim not in preferred]
    return combined.transpose(*preferred, *other)


def append_or_replace_ratio(
    existing: Optional[xr.DataArray],
    new: xr.DataArray,
    ratio_id: Optional[str] = None,
) -> xr.DataArray:
    """
    Backwards-compatible wrapper.

    The calibration module now stores products along the 'pair' dimension.  The
    old function name is kept so external callers do not break.
    """

    if "ratio" in new.dims:
        new = new.rename({"ratio": "pair"})
    return append_or_replace_pairs(existing, new)


def _filter_info_with_ratio_type(info: Optional[xr.DataArray]) -> Optional[xr.DataArray]:
    """Keep only pol_cal_info pairs with a valid, non-empty ratio_type."""

    if info is None:
        return None

    if "pair" not in info.dims:
        raise ValueError(f"pol_cal_info object has no pair dimension: {info.dims}")

    if "ratio_type" not in info.parameters.values:
        return info.isel(pair=slice(0, 0))

    ratio_type = info.sel(parameters="ratio_type")
    keep = []

    for pair, value in zip(ratio_type.pair.values, ratio_type.values):
        if value is None:
            continue

        try:
            if bool(np.isnan(value)):
                continue
        except Exception:
            pass

        value_str = str(value).strip().lower()
        if value_str in ["", "nan", "none"]:
            continue

        keep.append(pair)

    if len(keep) == 0:
        return info.isel(pair=slice(0, 0))

    return info.sel(pair=keep)


def append_or_replace_info(
    existing: Optional[xr.DataArray],
    new: xr.DataArray,
) -> xr.DataArray:
    """Add derived pair entries to pol_cal_info, replacing duplicate pair ids.

    Only pairs with a valid ratio_type are retained.  This prevents the original
    base channel-pair metadata, which can use temporary/dummy pair ids such as
    0 and 1, from being concatenated into the derived pol_cal_info products.
    """

    new = _filter_info_with_ratio_type(new)

    if new is None:
        raise ValueError("New pol_cal_info object cannot be None.")

    if existing is None:
        return new.transpose("parameters", "pair")

    existing = _filter_info_with_ratio_type(existing)

    if existing is None:
        return new.transpose("parameters", "pair")

    if "pair" not in existing.dims or "pair" not in new.dims:
        raise ValueError("Both existing and new pol_cal_info objects need a pair dimension.")

    drop_ids = [pid for pid in new.pair.values if pid in existing.pair.values]
    if drop_ids:
        existing = existing.drop_sel(pair=drop_ids)

    if existing.sizes.get("pair", 0) == 0:
        combined = new
    elif new.sizes.get("pair", 0) == 0:
        combined = existing
    else:
        combined = xr.concat(
            [existing.astype(object), new.astype(object)],
            dim="pair",
            join="outer",
            combine_attrs="override",
        )

    return combined.transpose("parameters", "pair")


def _ensure_pol_cal_output_dicts(output_data: Dict[str, Dict[str, Any]]) -> None:
    output_data.setdefault("pol_cal_ratio", {})
    output_data.setdefault("pol_cal_ratio_error", {})
    output_data.setdefault("pol_cal_info", {})


def _ensure_molecular_output_dicts(output_data: Dict[str, Dict[str, Any]]) -> None:
    output_data.setdefault("molecular_ratio", {})
    output_data.setdefault("molecular_info", {})


def _resolve_qa_alias(
    output_data: Dict[str, Dict[str, Any]],
    requested_key: str,
    loading_map: Optional[Dict[str, str]] = None,
    required_store: str = "profile",
) -> Optional[str]:
    """Return the actual key for requested_key, honoring loading_map aliases."""

    loading_map = loading_map or {}
    store = output_data.get(required_store, {})

    if requested_key in store:
        return requested_key

    mapped_key = loading_map.get(requested_key)
    if mapped_key in store:
        return mapped_key

    return None


def _with_singleton_time(
    da: xr.DataArray,
    template: Optional[xr.DataArray] = None,
) -> xr.DataArray:
    """Return da with a singleton time dimension for storage.

    Region statistics are stored in pol_cal_info and must not carry a time
    dimension.  This helper is therefore intended only for pol_cal_ratio and
    pol_cal_ratio_error products after time averaging has already happened.
    If possible, the singleton time coordinate is taken from the middle time
    coordinate of the template array.
    """

    if "time" in da.dims:
        if da.sizes.get("time", 0) == 1:
            preferred = [dim for dim in ["time", "pair", "bins"] if dim in da.dims]
            other = [dim for dim in da.dims if dim not in preferred]
            return da.transpose(*preferred, *other)

        raise ValueError(
            "Cannot add singleton time to a DataArray that already has "
            f"{da.sizes['time']} time entries."
        )

    time_value = 0
    if template is not None and "time" in template.dims and template.sizes.get("time", 0) > 0:
        mid_index = template.sizes["time"] // 2
        time_value = template["time"].values[mid_index]

    da = da.expand_dims(time=[time_value])
    preferred = [dim for dim in ["time", "pair", "bins"] if dim in da.dims]
    other = [dim for dim in da.dims if dim not in preferred]

    return da.transpose(*preferred, *other)


def _drop_singleton_time(
    da: xr.DataArray,
    name: str = "array",
) -> xr.DataArray:
    """
    Drop singleton time dimension from mean products.

    Mean products in this pipeline may still have dims:
        time, channel, bins
    or:
        time, pair, bins

    where time has size 1. The singleton time coordinate must be removed
    before combining different QA tests such as pcb_p45 and pcb_m45, because
    their singleton time coordinates are different and xarray would align
    them by time.
    """

    if "time" not in da.dims:
        return da

    if da.sizes["time"] != 1:
        raise ValueError(
            f"Cannot drop time from {name}: expected singleton time dimension, "
            f"got {da.sizes['time']} time entries."
        )

    return da.squeeze("time", drop=True)


def _drop_singleton_time_from_pair(
    values: xr.DataArray,
    errors: xr.DataArray,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Compatibility wrapper for mean products.

    The name is kept because several helpers still call this function, but
    with the current pipeline profile_mean/profile_error_mean are already
    averaged. Therefore this function no longer averages over time. It only
    drops singleton time dimensions to prevent xarray alignment issues.
    """

    values = _drop_singleton_time(values, name="values")
    errors = _drop_singleton_time(errors, name="errors")

    return values, errors

def _product_sqrt_error(
    a: xr.DataArray,
    a_error: xr.DataArray,
    b: xr.DataArray,
    b_error: xr.DataArray,
    product_sqrt: xr.DataArray,
) -> xr.DataArray:
    """Error propagation for sqrt(a * b)."""

    rel_a = a_error / a
    rel_b = b_error / b
    return np.abs(product_sqrt) * 0.5 * np.sqrt(rel_a**2 + rel_b**2)


def _division_error(
    numerator: xr.DataArray,
    numerator_error: xr.DataArray,
    denominator: xr.DataArray,
    denominator_error: xr.DataArray,
    quotient: xr.DataArray,
) -> xr.DataArray:
    """Independent error propagation for quotient = numerator / denominator."""

    rel_num = numerator_error / numerator
    rel_den = denominator_error / denominator
    return np.abs(quotient) * np.sqrt(rel_num**2 + rel_den**2)


def _numeric_parameter(
    pci: xr.DataArray,
    parameter: str,
    pair_ids: Sequence[str],
    default: Optional[float] = None,
) -> xr.DataArray:
    """Read a numeric parameter row from pol_cal_info for selected pairs."""

    if parameter in pci.parameters.values:
        out = pci.sel(parameters=parameter, pair=pair_ids).astype("float64")
    elif default is not None:
        out = xr.DataArray(
            np.full(len(pair_ids), default, dtype="float64"),
            dims=["pair"],
            coords={"pair": pair_ids},
        )
    else:
        raise KeyError(f"Missing required pol_cal_info parameter: {parameter}")

    return out


def _select_pairs_by_type(pci: xr.DataArray, ratio_type: str) -> List[str]:
    """Return pair ids whose ratio_type row matches ratio_type."""

    if "ratio_type" not in pci.parameters.values:
        return []

    types = pci.sel(parameters="ratio_type")
    return [
        str(pair)
        for pair, value in zip(types.pair.values, types.values)
        if value == ratio_type
    ]


def _pair_info_for_channels(
    pci: xr.DataArray,
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    pair_ids: Sequence[str],
    ratio_type: str,
) -> xr.DataArray:
    """Create a pol_cal_info object for derived pair products."""

    info = pci.sel(pair=pci.pair.values[:len(pair_ids)]).copy()
    info = info.assign_coords(pair=list(pair_ids))
    info = add_parameter(info, "ch_r", list(ch_r))
    info = add_parameter(info, "ch_t", list(ch_t))
    info = add_parameter(info, "ratio_type", ratio_type)
    return info



def _base_pol_cal_pairs(pci: xr.DataArray) -> xr.DataArray:
    """Return input channel-pair rows, excluding already derived ratio products."""

    if "ratio_type" not in pci.parameters.values:
        return pci

    ratio_type = pci.sel(parameters="ratio_type")
    keep = []
    for pair, value in zip(ratio_type.pair.values, ratio_type.values):
        if value is None:
            keep.append(pair)
            continue
        try:
            if bool(np.isnan(value)):
                keep.append(pair)
                continue
        except Exception:
            pass
        if str(value).lower() in ["", "nan", "none"]:
            keep.append(pair)

    if len(keep) == 0:
        return pci

    return pci.sel(pair=keep)

def _read_channels_from_info(pci: xr.DataArray) -> Tuple[List[str], List[str]]:
    ch_r = pci.sel(parameters="ch_r").values.tolist()
    ch_t = pci.sel(parameters="ch_t").values.tolist()
    return ch_r, ch_t


def _channel_parameter(
    channel_info: xr.DataArray,
    channel: str,
    parameter: str,
    default: float = 0.0,
) -> float:
    """
    Best-effort extraction of a channel parameter from channel_info.

    This supports both common layouts:
        channel_info(parameters, channel)
        channel_info(channel, parameters)
    """

    try:
        if "parameters" in channel_info.dims and "channel" in channel_info.dims:
            return float(channel_info.sel(parameters=parameter, channel=channel).values)
    except Exception:
        pass

    try:
        if "parameter" in channel_info.dims and "channel" in channel_info.dims:
            return float(channel_info.sel(parameter=parameter, channel=channel).values)
    except Exception:
        pass

    return float(default)


def _ideal_GH_for_channel(channel_id: str) -> Tuple[float, float]:
    """Return ideal/default G,H values for one channel id."""

    G = 1.0
    if channel_id[5] == "c":
        H = -1.0
    elif channel_id[5] == "p":
        H = 1.0
    else:
        H = 0.0

    return G, H

def _ideal_eta_for_channel(channel_id: str) -> Tuple[float, float]:
    """Return ideal/default eta values for one channel id."""
    """The reason that is part is needed is that the ideal calibration factor
    for a pair of total and parallel or cross channels is 2 or 0.5, NOT 1
    The reason is that an ideal analyser introduces a transmission of 1/2 
    for a parallel or cross channel but the total channel has no analyser.
    In parallel and cross channel ratios the 1/2 factors cancel out. This is 
    not the case for ratios that include the total channel."""

    if channel_id[5] == "c" or channel_id[5] == "p":
        eta = 0.5
    else:
        eta = 1.0

    return eta

def _as_1d_list(values) -> List[Any]:
    """Convert scalar, numpy scalar, DataArray values, list, tuple to a 1D list."""

    if isinstance(values, xr.DataArray):
        values = values.values

    arr = np.atleast_1d(values)

    return arr.tolist()

def _ideal_eta_for_pairs(
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    pair_ids: Sequence[str],
) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """Return ideal pair-based G_R, G_T, H_R, H_T arrays."""

    ch_r = _as_1d_list(ch_r)
    ch_t = _as_1d_list(ch_t)
    pair_ids = _as_1d_list(pair_ids)

    if not (len(ch_r) == len(ch_t) == len(pair_ids)):
        raise ValueError(
            "Cannot create ideal GH parameters because lengths differ: "
            f"len(ch_r)={len(ch_r)}, len(ch_t)={len(ch_t)}, "
            f"len(pair_ids)={len(pair_ids)}."
        )

    eta_values = []

    for r, t in zip(ch_r, ch_t):
        eta_r = _ideal_eta_for_channel(r)
        eta_t = _ideal_eta_for_channel(t)
        
        eta = eta_r / eta_t

        eta_values.append(eta)
        
    coords = {"pair": pair_ids}

    eta_arr = xr.DataArray(
        eta_values, dims=["pair"], coords=coords, name="eta"
        ).astype("float64")
    
    return eta_arr 
    

def _ideal_GH_for_pairs(
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    pair_ids: Sequence[str],
) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """Return ideal pair-based G_R, G_T, H_R, H_T arrays."""

    ch_r = _as_1d_list(ch_r)
    ch_t = _as_1d_list(ch_t)
    pair_ids = _as_1d_list(pair_ids)

    if not (len(ch_r) == len(ch_t) == len(pair_ids)):
        raise ValueError(
            "Cannot create ideal GH parameters because lengths differ: "
            f"len(ch_r)={len(ch_r)}, len(ch_t)={len(ch_t)}, "
            f"len(pair_ids)={len(pair_ids)}."
        )

    G_R_values = []
    H_R_values = []
    G_T_values = []
    H_T_values = []

    for r, t in zip(ch_r, ch_t):
        G_R, H_R = _ideal_GH_for_channel(r)
        G_T, H_T = _ideal_GH_for_channel(t)

        G_R_values.append(G_R)
        H_R_values.append(H_R)
        G_T_values.append(G_T)
        H_T_values.append(H_T)

    coords = {"pair": pair_ids}

    return (
        xr.DataArray(G_R_values, dims=["pair"], coords=coords, name="G_R").astype("float64"),
        xr.DataArray(G_T_values, dims=["pair"], coords=coords, name="G_T").astype("float64"),
        xr.DataArray(H_R_values, dims=["pair"], coords=coords, name="H_R").astype("float64"),
        xr.DataArray(H_T_values, dims=["pair"], coords=coords, name="H_T").astype("float64"),
    )

def _finite_or_default(values: xr.DataArray, defaults: xr.DataArray) -> xr.DataArray:
    """Replace non-finite channel-parameter values by defaults."""

    values = values.astype("float64")
    defaults = defaults.astype("float64")
    return values.where(np.isfinite(values), defaults)


def _channel_to_pair_GH(
    da: xr.DataArray,
    defaults: xr.DataArray,
    pair_ids: Sequence[str],
    name: str,
) -> xr.DataArray:
    """Convert a channel-selected G/H row to pair coordinates."""

    da = da.assign_coords(channel=list(pair_ids)).rename({"channel": "pair"})
    da = da.assign_coords(pair=list(pair_ids)).rename(name)

    return _finite_or_default(da, defaults.rename(name))


def _channel_info_GH_for_pairs(
    channel_info_qa: xr.DataArray,
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    pair_ids: Sequence[str],
) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """Read G/H from channel_info and align them to pair ids.

    channel_info_qa is expected to have dims (channel, parameters).  The channel
    coordinate uses atlas_channel_id values compatible with ch_r/ch_t.
    Missing or non-finite values fall back to the ideal/default GH values.
    """

    G_R_default, G_T_default, H_R_default, H_T_default = _ideal_GH_for_pairs(
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=pair_ids,
    )

    try:
        G_R = channel_info_qa.sel(channel=list(ch_r), parameters="G")
        H_R = channel_info_qa.sel(channel=list(ch_r), parameters="H")
        G_T = channel_info_qa.sel(channel=list(ch_t), parameters="G")
        H_T = channel_info_qa.sel(channel=list(ch_t), parameters="H")
    except Exception as exc:
        raise KeyError(
            "Could not read G/H from channel_info for the requested channel pairs."
        ) from exc

    return (
        _channel_to_pair_GH(G_R, G_R_default, pair_ids, "G_R"),
        _channel_to_pair_GH(G_T, G_T_default, pair_ids, "G_T"),
        _channel_to_pair_GH(H_R, H_R_default, pair_ids, "H_R"),
        _channel_to_pair_GH(H_T, H_T_default, pair_ids, "H_T"),
    )


def _GH_correct_ratio(
    ratio: xr.DataArray,
    ratio_error: xr.DataArray,
    G_R: xr.DataArray,
    G_T: xr.DataArray,
    H_R: xr.DataArray,
    H_T: xr.DataArray,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Apply the GH correction and propagate only ratio random error."""

    A = G_T + H_T
    B = G_R + H_R
    C = G_R - H_R
    D = G_T - H_T

    denominator = C - ratio * D
    corrected = (ratio * A - B) / denominator
    corrected = corrected.where(denominator != 0)

    derivative = (A * C - D * B) / denominator**2
    corrected_error = np.abs(derivative) * ratio_error
    corrected_error = corrected_error.where(denominator != 0)

    return corrected, corrected_error


def _GH_correct_values(
    ratio: xr.DataArray,
    G_R: xr.DataArray,
    G_T: xr.DataArray,
    H_R: xr.DataArray,
    H_T: xr.DataArray,
) -> xr.DataArray:
    """Apply the GH correction without error propagation."""

    A = G_T + H_T
    B = G_R + H_R
    C = G_R - H_R
    D = G_T - H_T

    denominator = C - ratio * D
    corrected = (ratio * A - B) / denominator
    corrected = corrected.where(denominator != 0)

    return corrected


def _add_GH_to_info(
    info: xr.DataArray,
    G_R: xr.DataArray,
    G_T: xr.DataArray,
    H_R: xr.DataArray,
    H_T: xr.DataArray,
) -> xr.DataArray:
    """Store GH parameters in pol_cal_info."""

    info = add_parameter(info, "G_R", G_R)
    info = add_parameter(info, "G_T", G_T)
    info = add_parameter(info, "H_R", H_R)
    info = add_parameter(info, "H_T", H_T)
    return info



def _find_eta_for_pairs(
    output_data: Dict[str, Dict[str, Any]],
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    calibrated_pair_ids: Sequence[str],
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Find scalar pcb/pcb_aux eta values and align to calibrated_pair_ids.

    The calibration-factor stage stores the calibration-region eta value and
    its SEM in output_data["pol_cal_info"] under the combined calibration keys
    ("pcb" and "pcb_aux") as ratio_type="eta" entries.  These scalar values
    are the calibration factors that ray signal ratios should be divided
    by.  The eta profiles stored in pol_cal_ratio are intentionally not used
    here, because ray must not be calibrated by a bin-resolved eta profile.
    """

    eta_pair_ids = [
        make_ratio_id(r, t, type_index="e")
        for r, t in zip(ch_r, ch_t)
    ]

    info_store = output_data.get("pol_cal_info", {})
    eta_sources = ["pcb", "pcb_aux"]
    eta_values = []
    eta_errors = []

    for eta_id, calibrated_id in zip(eta_pair_ids, calibrated_pair_ids):
        eta_value = None
        eta_error = None

        for source in eta_sources:
            if source not in info_store:
                continue

            source_info = _filter_info_with_ratio_type(info_store[source])
            if source_info is None or source_info.sizes.get("pair", 0) == 0:
                continue
            if "pair" not in source_info.dims or eta_id not in source_info.pair.values:
                continue

            source_eta_ids = _select_pairs_by_type(source_info, "eta")
            if eta_id not in source_eta_ids:
                continue

            if "mean" not in source_info.parameters.values:
                raise KeyError(
                    f"Eta entry {eta_id!r} in pol_cal_info[{source!r}] has no 'mean' row."
                )
            if "sem" not in source_info.parameters.values:
                raise KeyError(
                    f"Eta entry {eta_id!r} in pol_cal_info[{source!r}] has no 'sem' row."
                )

            eta_value = (
                source_info
                .sel(parameters="mean", pair=eta_id)
                .astype("float64")
                .reset_coords(drop=True)
            )
            eta_error = (
                source_info
                .sel(parameters="sem", pair=eta_id)
                .astype("float64")
                .reset_coords(drop=True)
            )
            break

        if eta_value is None or eta_error is None:
            continue

        eta_values.append(
            eta_value.expand_dims(pair=[calibrated_id]).rename("eta")
        )
        eta_errors.append(
            eta_error.expand_dims(pair=[calibrated_id]).rename("eta_error")
        )

    if not eta_values:
        raise KeyError(
            "No matching scalar eta entries found in pcb/pcb_aux pol_cal_info "
            "for ray channel pairs. Run compute_calibration_factor first."
        )

    eta = xr.concat(eta_values, dim="pair").astype("float64")
    eta_error = xr.concat(eta_errors, dim="pair").astype("float64")

    eta = eta.transpose("pair")
    eta_error = eta_error.transpose("pair")

    return eta, eta_error


def _collect_eta_pairs_from_pcb(
    output_data: Dict[str, Dict[str, Any]],
) -> Tuple[List[str], List[str], List[str], xr.DataArray]:
    """Collect eta channel combinations from pcb/pcb_aux pol_cal_info.

    Returns ch_r, ch_t, eta_pair_ids, eta_info.  These eta entries are the
    authoritative source of the channel combinations used by ray calibrated
    ratio and VLDR calculations.
    """

    info_store = output_data.get("pol_cal_info", {})
    eta_infos = []

    for source in ["pcb", "pcb_aux"]:
        if source not in info_store:
            continue

        source_info = _filter_info_with_ratio_type(info_store[source])
        if source_info is None or source_info.sizes.get("pair", 0) == 0:
            continue

        eta_ids = _select_pairs_by_type(source_info, "eta")
        if len(eta_ids) == 0:
            continue

        eta_infos.append(source_info.sel(pair=eta_ids))

    if len(eta_infos) == 0:
        return [], [], [], xr.DataArray()

    eta_info = xr.concat(
        [info.astype(object) for info in eta_infos],
        dim="pair",
        join="outer",
        combine_attrs="override",
    ).transpose("parameters", "pair")

    # Drop duplicate eta ids while preserving first occurrence.
    seen = set()
    keep = []
    for pair in eta_info.pair.values:
        pair_str = str(pair)
        if pair_str in seen:
            continue
        seen.add(pair_str)
        keep.append(pair)
    eta_info = eta_info.sel(pair=keep)

    ch_r, ch_t = _read_channels_from_info(eta_info)
    eta_pair_ids = [str(pair) for pair in eta_info.pair.values]

    return ch_r, ch_t, eta_pair_ids, eta_info



def pldr_error(
    delta_m: xr.DataArray,
    delta_v_err: xr.DataArray,
    delta_p_ulim: float = 0.3,
    delta_p_err_ulim: float = 0.025,
):
    """
    Calculate PLDR error lookup for each pair lazily.

    The previous implementation looped over pairs and extracted .values from
    argmax/argmin results, which forced computation. This version keeps the
    lookup table and the sr_limit selection lazy by using vectorized xarray
    operations along the R dimension.
    """

    if "pair" not in delta_m.dims:
        raise ValueError("delta_m must have dimension 'pair'.")
    if "pair" not in delta_v_err.dims:
        raise ValueError("delta_v_err must have dimension 'pair'.")

    delta_m, delta_v_err = xr.align(delta_m, delta_v_err, join="outer")

    R_values = np.arange(1.01, 3.0, 0.001)
    delta_p_values = np.arange(0.0, delta_p_ulim + 0.001, 0.001)

    R = xr.DataArray(R_values, dims=["R"], coords={"R": R_values}, name="R")
    delta_p = xr.DataArray(
        delta_p_values,
        dims=["delta_p"],
        coords={"delta_p": delta_p_values},
        name="delta_p",
    )

    sq_term_nom = (delta_v_err + delta_p) * (1.0 + delta_m) ** 2 * R**2
    sq_term_denom = (1.0 + delta_m) ** 2 * R**2

    lin_term_nom = (
        (1.0 + delta_m)
        * (
            delta_v_err * (delta_p - 2.0 * delta_m)
            - delta_p * (1.0 + delta_m)
        )
        * R
    )
    lin_term_denom = -(
        (1.0 + delta_m) * (delta_v_err + 1.0 + delta_m) * R
    )

    const_term_nom = -delta_m * (delta_p - delta_m) * delta_v_err
    const_term_denom = (delta_p - delta_m) * delta_v_err

    delta_p_err = (
        (sq_term_nom + lin_term_nom + const_term_nom)
        / (sq_term_denom + lin_term_denom + const_term_denom)
        - delta_p
    )

    delta_p_err = delta_p_err.where(np.abs(delta_p_err) <= delta_p_err_ulim)
    delta_p_err.name = "delta_p_err"
    delta_p_err = delta_p_err.transpose("delta_p", "R", "pair")

    last_delta_p_err = delta_p_err.isel(delta_p=-1)
    valid = last_delta_p_err.notnull().any("R")

    # Vectorized R-position lookup. Filling NaNs avoids argmax/argmin failures
    # on partially missing rows. Fully missing rows have no valid SR solution
    # and are masked to NaN below.
    idx_max = last_delta_p_err.fillna(-np.inf).argmax(dim="R", skipna=False)
    idx_min = last_delta_p_err.fillna(np.inf).argmin(dim="R", skipna=False)

    # Avoid R.isel(R=idx_max).  R is a 1D coordinate-backed DataArray, while
    # idx_max/idx_min can carry remaining dimensions such as ("pair",).  Some
    # xarray/pandas versions fail on that vectorized coordinate-indexing path.
    R_values = np.asarray(R.values)
    idx_max_values = np.asarray(idx_max.values).astype(int)
    idx_min_values = np.asarray(idx_min.values).astype(int)

    sr_max = xr.DataArray(
        R_values[idx_max_values],
        dims=idx_max.dims,
        coords=idx_max.coords,
        name="sr_max",
    )
    sr_min = xr.DataArray(
        R_values[idx_min_values],
        dims=idx_min.dims,
        coords=idx_min.coords,
        name="sr_min",
    )

    sr_default = xr.full_like(delta_v_err.astype("float64"), 1.01)
    sr_nan = xr.full_like(delta_v_err.astype("float64"), np.nan)

    min_bsc_ratio = xr.where(
        delta_v_err > 0.0001,
        sr_max,
        xr.where(delta_v_err < -0.0001, sr_min, sr_default),
    )

    # If the PLDR-error curve has no finite values along R, the SR limit is
    # not derivable.  Treat it like an out-of-range / too-large case and keep
    # the stored sr_limit as NaN.
    min_bsc_ratio = min_bsc_ratio.where(valid, sr_nan)
    min_bsc_ratio = min_bsc_ratio.rename("sr_limit")

    return delta_p_err, delta_p, R, min_bsc_ratio


def epsilon_angle(
    eta_p45: xr.DataArray,
    eta_m45: xr.DataArray,
    kappa: float = 1.0,
) -> xr.DataArray:
    """Calculate epsilon angle from +45 and -45 eta_s_f products."""

    psi = (eta_p45 - eta_m45) / (eta_p45 + eta_m45)
    epsilon = np.rad2deg(
        0.5 * np.arcsin(
            np.tan(0.5 * np.arcsin(psi)) / kappa
        )
    )
    return epsilon.rename("epsilon")



def add_parameters(
    info: xr.DataArray,
    parameters: Dict[str, Any],
) -> xr.DataArray:
    """Add several parameter rows to pol_cal_info."""

    for name, values in parameters.items():
        info = add_parameter(info, name=name, values=values)

    return info


def add_region_stats_to_info(
    info: xr.DataArray,
    values: xr.DataArray,
    z_pair: xr.DataArray,
    averaging_range,
    ratio_type: str,
) -> xr.DataArray:
    """Add regional mean, regional SEM, and ratio_type to pol_cal_info."""

    mean = mean_in_region(
        da=values,
        z=z_pair,
        averaging_range=averaging_range,
    )

    sem = sem_in_region(
        da=values,
        z=z_pair,
        averaging_range=averaging_range,
    )

    return add_parameters(
        info,
        {
            "mean": mean,
            "sem": sem,
            "ratio_type": ratio_type,
        },
    )


def store_pol_cal_mean_product(
    output_data: Dict[str, Dict[str, Any]],
    key: str,
    values: xr.DataArray,
    errors: xr.DataArray,
    info: xr.DataArray,
) -> None:
    """Store one mean polarization-calibration product."""

    output_data["pol_cal_ratio_mean"][key] = append_or_replace_pairs(
        output_data["pol_cal_ratio_mean"].get(key),
        values.rename("ratio"),
    )

    output_data["pol_cal_ratio_error_mean"][key] = append_or_replace_pairs(
        output_data["pol_cal_ratio_error_mean"].get(key),
        errors.rename("ratio_error"),
    )

    output_data["pol_cal_info"][key] = append_or_replace_info(
        output_data["pol_cal_info"].get(key),
        info,
    )


def pair_info_from_base(
    pci_base: xr.DataArray,
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    type_index: str,
    ratio_type: str,
) -> Tuple[List[str], xr.DataArray]:
    """Create derived pair ids and matching pol_cal_info."""

    pair_ids = [
        make_ratio_id(r, t, type_index=type_index)
        for r, t in zip(ch_r, ch_t)
    ]

    info = _pair_info_for_channels(
        pci=pci_base,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=pair_ids,
        ratio_type=ratio_type,
    )

    return pair_ids, info


def ratio_from_channel_pairs(
    sig: xr.DataArray,
    sig_err: xr.DataArray,
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    info: xr.DataArray,
    average_time: bool = True,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Create pair-based ratio and propagated ratio error."""

    sig_r = sig.sel(channel=ch_r)
    sig_t = sig.sel(channel=ch_t)
    sig_r_err = sig_err.sel(channel=ch_r)
    sig_t_err = sig_err.sel(channel=ch_t)

    if average_time:
        # Mean products are already averaged and only carry singleton time.
        # Drop it before pair operations to avoid time-coordinate alignment issues.
        sig_r, sig_r_err = _drop_singleton_time_from_pair(sig_r, sig_r_err)
        sig_t, sig_t_err = _drop_singleton_time_from_pair(sig_t, sig_t_err)

    ratio = simple_ratio(
        numerator=sig_r,
        denominator=sig_t,
        info=info,
    )

    ratio_error = ratio_error_independent(
        numerator=sig_r,
        denominator=sig_t,
        numerator_error=sig_r_err,
        denominator_error=sig_t_err,
        ratio_values=ratio,
        info=info,
    )

    return ratio.rename("ratio"), ratio_error.rename("ratio_error")


def mean_gain_ratio_from_profiles(
    qa_test: str,
    profiles: Dict[str, xr.DataArray],
    profile_errors: Dict[str, xr.DataArray],
    channel_info: Dict[str, xr.DataArray],
    pol_cal_info: Dict[str, xr.DataArray],
    vertical_scale: Dict[str, xr.DataArray],
    averaging_range,
) -> Optional[Tuple[xr.DataArray, xr.DataArray, xr.DataArray]]:
    """Build one individual +/-45 mean gain-ratio product."""

    required_stores = [profiles, profile_errors, channel_info, pol_cal_info, vertical_scale]
    if any(qa_test not in store for store in required_stores):
        return None

    pci_base = _base_pol_cal_pairs(pol_cal_info[qa_test])
    ch_r, ch_t = _read_channels_from_info(pci_base)

    if len(ch_r) == 0:
        return None

    _, info = pair_info_from_base(
        pci_base=pci_base,
        ch_r=ch_r,
        ch_t=ch_t,
        type_index="g",
        ratio_type="gain_ratio",
    )

    ratio, ratio_error = ratio_from_channel_pairs(
        sig=profiles[qa_test],
        sig_err=profile_errors[qa_test],
        ch_r=ch_r,
        ch_t=ch_t,
        info=info,
    )

    z_pair = channels_to_pairs(
        vertical_scale[qa_test].sel(channel=ch_r),
        info,
    )

    info = add_region_stats_to_info(
        info=info,
        values=ratio,
        z_pair=z_pair,
        averaging_range=averaging_range,
        ratio_type="gain_ratio",
    )

    return ratio, ratio_error, info


def mean_gain_ratio_inputs_exist(
    output_data: Dict[str, Dict[str, Any]],
    vertical_scale: Dict[str, xr.DataArray],
    p45_key: str,
    m45_key: str,
) -> bool:
    """Check whether combined gain-ratio inputs exist."""

    return (
        p45_key in output_data["pol_cal_ratio_mean"]
        and m45_key in output_data["pol_cal_ratio_mean"]
        and p45_key in output_data["pol_cal_ratio_error_mean"]
        and m45_key in output_data["pol_cal_ratio_error_mean"]
        and p45_key in output_data["pol_cal_info"]
        and m45_key in output_data["pol_cal_info"]
        and p45_key in vertical_scale
        and m45_key in vertical_scale
    )


def build_combined_gain_ratio(
    output_data: Dict[str, Dict[str, Any]],
    vertical_scale: Dict[str, xr.DataArray],
    target_key: str,
    p45_key: str,
    m45_key: str,
    averaging_range,
) -> Optional[Tuple[xr.DataArray, xr.DataArray, xr.DataArray]]:
    """Build combined pcb/pcb_aux eta_s from +45 and -45 gain ratios."""

    if not mean_gain_ratio_inputs_exist(output_data, vertical_scale, p45_key, m45_key):
        return None

    ratio_store = output_data["pol_cal_ratio_mean"]
    error_store = output_data["pol_cal_ratio_error_mean"]
    info_store = output_data["pol_cal_info"]

    pci_p45 = info_store[p45_key]
    pci_m45 = info_store[m45_key]

    p45_gain_ids = _select_pairs_by_type(pci_p45, "gain_ratio")
    m45_gain_ids = _select_pairs_by_type(pci_m45, "gain_ratio")
    gain_ids = [pair_id for pair_id in p45_gain_ids if pair_id in m45_gain_ids]

    if len(gain_ids) == 0:
        return None

    ratio_p45 = ratio_store[p45_key].sel(pair=gain_ids)
    ratio_m45 = ratio_store[m45_key].sel(pair=gain_ids)
    ratio_p45_error = error_store[p45_key].sel(pair=gain_ids)
    ratio_m45_error = error_store[m45_key].sel(pair=gain_ids)

    # Defensive compatibility with older intermediates that still carry time.
    ratio_p45, ratio_p45_error = _drop_singleton_time_from_pair(ratio_p45, ratio_p45_error)
    ratio_m45, ratio_m45_error = _drop_singleton_time_from_pair(ratio_m45, ratio_m45_error)

    info_p45 = pci_p45.sel(pair=gain_ids)
    info_m45 = pci_m45.sel(pair=gain_ids)

    ch_r_p45, ch_t_p45 = _read_channels_from_info(info_p45)
    ch_r_m45, ch_t_m45 = _read_channels_from_info(info_m45)

    if ch_r_p45 != ch_r_m45 or ch_t_p45 != ch_t_m45:
        raise ValueError(
            f"Gain-ratio channel mismatch between {p45_key} and {m45_key} "
            f"for pairs {gain_ids}."
        )

    ch_r = ch_r_p45
    ch_t = ch_t_p45

    z_pair_p45 = channels_to_pairs(
        vertical_scale[p45_key].sel(channel=ch_r),
        info_p45,
    )
    z_pair_m45 = channels_to_pairs(
        vertical_scale[m45_key].sel(channel=ch_r),
        info_m45,
    )

    if not z_pair_p45.broadcast_equals(z_pair_m45):
        CustomWarning(
            f"Skipping combined gain ratio for {target_key}: "
            f"vertical scales of {p45_key} and {m45_key} do not match."
        )
        return None

    eta_s_product = ratio_p45 * ratio_m45
    eta_s = np.sqrt(eta_s_product.where(eta_s_product >= 0))
    eta_s = eta_s.assign_coords(pair=gain_ids).rename("ratio")

    eta_s_error = _product_sqrt_error(
        ratio_p45,
        ratio_p45_error,
        ratio_m45,
        ratio_m45_error,
        eta_s,
    )
    eta_s_error = eta_s_error.assign_coords(pair=gain_ids).rename("ratio_error")

    info = _pair_info_for_channels(
        pci=info_p45,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=gain_ids,
        ratio_type="gain_ratio",
    )

    epsilon = epsilon_angle(
        eta_p45=ratio_p45.assign_coords(pair=gain_ids),
        eta_m45=ratio_m45.assign_coords(pair=gain_ids),
        kappa=1.0,
    )

    info = add_parameters(
        info,
        {
            "mean": mean_in_region(eta_s, z_pair_p45, averaging_range),
            "sem": sem_in_region(eta_s, z_pair_p45, averaging_range),
            "epsilon": mean_in_region(epsilon, z_pair_p45, averaging_range),
            "epsilon_error": sem_in_region(epsilon, z_pair_p45, averaging_range),
            "ratio_type": "gain_ratio",
        },
    )

    return eta_s, eta_s_error, info


def calibration_factor_products(
    eta_s: xr.DataArray,
    eta_s_error: xr.DataArray,
    pci_gain: xr.DataArray,
    gain_ids: Sequence[str],
) -> List[Dict[str, Any]]:
    """Build eta_s_f and eta products from combined gain-ratio eta_s."""

    ch_r, ch_t = _read_channels_from_info(pci_gain)

    eta_s_f_ids = [
        make_ratio_id(r, t, type_index="f")
        for r, t in zip(ch_r, ch_t)
    ]
    eta_ids = [
        make_ratio_id(r, t, type_index="e")
        for r, t in zip(ch_r, ch_t)
    ]

    trans_ratio = _numeric_parameter(
        pci=pci_gain,
        parameter="R_to_T_transmission_ratio",
        pair_ids=gain_ids,
        default=1.0,
    ).assign_coords(pair=gain_ids)

    K = _numeric_parameter(
        pci=pci_gain,
        parameter="K",
        pair_ids=gain_ids,
        default=1.0,
    ).assign_coords(pair=gain_ids)

    eta_s_f = eta_s / trans_ratio
    eta_s_f_error = eta_s_error / np.abs(trans_ratio)

    eta = eta_s_f / K
    eta_error = eta_s_f_error / np.abs(K)

    return [
        {
            "ratio_type": "eta_s_f",
            "pair_ids": eta_s_f_ids,
            "values": eta_s_f,
            "errors": eta_s_f_error,
            "ch_r": ch_r,
            "ch_t": ch_t,
        },
        {
            "ratio_type": "eta",
            "pair_ids": eta_ids,
            "values": eta,
            "errors": eta_error,
            "ch_r": ch_r,
            "ch_t": ch_t,
        },
    ]


def inherit_entries(
    output_data: Dict[str, Dict[str, Any]],
    target_key: str,
    source_keys: Sequence[str],
    store_keys: Sequence[str],
    overwrite: bool = False,
) -> None:
    """Copy selected store entries from the first available source key."""

    for store_key in store_keys:
        if store_key not in output_data:
            continue

        store = output_data[store_key]
        if not overwrite and target_key in store:
            continue

        for source_key in source_keys:
            if source_key in store:
                store[target_key] = store[source_key]
                break


def inherit_pcb_coordinate_entries(
    output_data: Dict[str, Dict[str, Any]],
    target_key: str,
    p45_key: str,
    m45_key: str,
) -> None:
    """Let combined pcb/pcb_aux keys inherit range/height entries."""

    inherit_entries(
        output_data=output_data,
        target_key=target_key,
        source_keys=[p45_key, m45_key],
        store_keys=["range", "height_agl", "height_asl"],
        overwrite=False,
    )


def compute_gain_ratio(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Compute mean gain ratios for pcb/pcb_aux +/-45 QA tests."""

    output_data = shallow_copy(input_data)
    vertical_scale = output_data[processing_info["caller_info"]["vertical_scale"]]
    averaging_range = processing_info["settings_info"]["pcb"]["calibration_region"]

    individual_tests = ["pcb_p45", "pcb_m45", "pcb_aux_p45", "pcb_aux_m45"]
    combined_tests = {"pcb": ("pcb_p45", "pcb_m45"), "pcb_aux": ("pcb_aux_p45", "pcb_aux_m45")}

    for qa_test in individual_tests:
        product = mean_gain_ratio_from_profiles(
            qa_test=qa_test,
            profiles=output_data["profile_mean"],
            profile_errors=output_data["profile_error_mean"],
            channel_info=output_data["channel_info"],
            pol_cal_info=output_data["pol_cal_info"],
            vertical_scale=vertical_scale,
            averaging_range=averaging_range,
        )

        if product is None:
            continue

        ratio, ratio_error, info = product
        store_pol_cal_mean_product(output_data, qa_test, ratio, ratio_error, info)

    for target_key, (p45_key, m45_key) in combined_tests.items():
        product = build_combined_gain_ratio(
            output_data=output_data,
            vertical_scale=vertical_scale,
            target_key=target_key,
            p45_key=p45_key,
            m45_key=m45_key,
            averaging_range=averaging_range,
        )

        if product is None:
            continue

        ratio, ratio_error, info = product
        store_pol_cal_mean_product(output_data, target_key, ratio, ratio_error, info)
        inherit_pcb_coordinate_entries(
            output_data=output_data,
            target_key=target_key,
            p45_key=p45_key,
            m45_key=m45_key,
        )

    print_entry("Gain ratio calculation complete!")
    return output_data


def compute_calibration_factor(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Compute eta_s_f and eta from combined pcb/pcb_aux mean gain ratios."""

    output_data = shallow_copy(input_data)
    vertical_scale = output_data[processing_info["caller_info"]["vertical_scale"]]
    averaging_range = processing_info["settings_info"]["pcb"]["calibration_region"]

    vertical_reference = {"pcb": "pcb_p45", "pcb_aux": "pcb_aux_p45"}

    for target_key in ["pcb", "pcb_aux"]:
        if target_key not in output_data["pol_cal_ratio_mean"]:
            continue
        if target_key not in output_data["pol_cal_ratio_error_mean"]:
            continue
        if target_key not in output_data["pol_cal_info"]:
            continue

        info_store = output_data["pol_cal_info"][target_key]
        gain_ids = _select_pairs_by_type(info_store, "gain_ratio")

        if len(gain_ids) == 0:
            continue

        eta_s = output_data["pol_cal_ratio_mean"][target_key].sel(pair=gain_ids).rename("ratio")
        eta_s_error = output_data["pol_cal_ratio_error_mean"][target_key].sel(pair=gain_ids).rename("ratio_error")
        pci_gain = info_store.sel(pair=gain_ids)

        eta_s, eta_s_error = _drop_singleton_time_from_pair(eta_s, eta_s_error)
        products = calibration_factor_products(eta_s, eta_s_error, pci_gain, gain_ids)

        z_pair_gain = None
        z_source = vertical_reference[target_key]

        if z_source in vertical_scale:
            ch_r, _ = _read_channels_from_info(pci_gain)
            z_pair_gain = channels_to_pairs(
                vertical_scale[z_source].sel(channel=ch_r),
                pci_gain,
            )
        else:
            CustomWarning(
                f"No vertical scale reference found for {target_key}; "
                "eta_s_f/eta mean and sem will not be added."
            )

        for product in products:
            pair_ids = product["pair_ids"]
            ratio_type = product["ratio_type"]
            values = product["values"].assign_coords(pair=pair_ids).rename("ratio")
            errors = product["errors"].assign_coords(pair=pair_ids).rename("ratio_error")

            info = _pair_info_for_channels(
                pci=pci_gain,
                ch_r=product["ch_r"],
                ch_t=product["ch_t"],
                pair_ids=pair_ids,
                ratio_type=ratio_type,
            )

            if z_pair_gain is not None:
                info = add_region_stats_to_info(
                    info=info,
                    values=values,
                    z_pair=z_pair_gain.assign_coords(pair=pair_ids),
                    averaging_range=averaging_range,
                    ratio_type=ratio_type,
                )
            else:
                info = add_parameter(info, name="ratio_type", values=ratio_type)

            store_pol_cal_mean_product(output_data, target_key, values, errors, info)

    print_entry("Calibration factor calculation complete!")
    return output_data


def store_pol_cal_product(
    output_data: Dict[str, Dict[str, Any]],
    key: str,
    values: xr.DataArray,
    errors: xr.DataArray,
    mean: bool = False,
) -> None:
    """Store one polarization-calibration ratio product."""

    ratio_key = "pol_cal_ratio_mean" if mean else "pol_cal_ratio"
    error_key = "pol_cal_ratio_error_mean" if mean else "pol_cal_ratio_error"

    output_data[ratio_key][key] = append_or_replace_pairs(
        output_data[ratio_key].get(key),
        values.rename("ratio"),
    )

    output_data[error_key][key] = append_or_replace_pairs(
        output_data[error_key].get(key),
        errors.rename("ratio_error"),
    )


def store_pol_cal_info(
    output_data: Dict[str, Dict[str, Any]],
    key: str,
    info: xr.DataArray,
) -> None:
    """Store one polarization-calibration info product."""

    output_data["pol_cal_info"][key] = append_or_replace_info(
        output_data["pol_cal_info"].get(key),
        info,
    )


def get_ray_alias(
    output_data: Dict[str, Dict[str, Any]],
    processing_info: Dict[str, Any],
    required_store: str,
) -> Optional[str]:
    """Resolve the ray key for a selected store."""

    loading_map = processing_info["caller_info"].get("loading_map", {})

    return _resolve_qa_alias(
        output_data=output_data,
        requested_key="ray",
        loading_map=loading_map,
        required_store=required_store,
    )




def get_ray_io_keys(
    output_data: Dict[str, Dict[str, Any]],
    processing_info: Dict[str, Any],
    required_store: str,
) -> Tuple[Optional[str], Optional[str]]:
    """Return input key and output key for ray-style products.

    If ray exists in the requested input store, read from and write to
    ray. Otherwise the alias/resolved key is used for both reading and
    writing, typically ray. This avoids creating synthetic ray entries
    when the original data only contains ray.
    """

    qa_alias = get_ray_alias(
        output_data=output_data,
        processing_info=processing_info,
        required_store=required_store,
    )

    if qa_alias is None:
        return None, None

    store = output_data.get(required_store, {})
    qa_store = "ray" if "ray" in store else qa_alias

    return qa_alias, qa_store


def collect_ray_eta_inputs(
    output_data: Dict[str, Dict[str, Any]],
    sig: xr.DataArray,
) -> Optional[Tuple[List[str], List[str], List[str], xr.DataArray]]:
    """Collect eta-based channel pairs that are available in a signal array."""

    ch_r, ch_t, eta_ids, eta_info = _collect_eta_pairs_from_pcb(output_data)

    if len(eta_ids) == 0:
        return None

    available_channels = set(sig.channel.values.tolist())
    keep = [
        i for i, (r, t) in enumerate(zip(ch_r, ch_t))
        if r in available_channels and t in available_channels
    ]

    if len(keep) == 0:
        return None

    ch_r = [ch_r[i] for i in keep]
    ch_t = [ch_t[i] for i in keep]
    eta_ids = [eta_ids[i] for i in keep]
    eta_info = eta_info.isel(pair=keep)

    return ch_r, ch_t, eta_ids, eta_info


def build_calibrated_ratio_product(
    output_data: Dict[str, Dict[str, Any]],
    sig: xr.DataArray,
    sig_err: xr.DataArray,
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    eta_ids: Sequence[str],
    eta_info: xr.DataArray,
    average_time: bool = False,
) -> Optional[Tuple[xr.DataArray, xr.DataArray, xr.DataArray, List[str], List[str]]]:
    """Build calibrated_ratio from signal ratios, scalar eta, and ideal GH."""

    pair_ids = [
        replace_ratio_id_type(eta_id, old_type_index="e", new_type_index="d")
        for eta_id in eta_ids
    ]

    info = _pair_info_for_channels(
        pci=eta_info,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=pair_ids,
        ratio_type="calibrated_ratio",
    )

    uncalibrated, uncalibrated_error = ratio_from_channel_pairs(
        sig=sig,
        sig_err=sig_err,
        ch_r=ch_r,
        ch_t=ch_t,
        info=info,
        average_time=average_time,
    )

    eta, eta_error = _find_eta_for_pairs(
        output_data=output_data,
        ch_r=ch_r,
        ch_t=ch_t,
        calibrated_pair_ids=pair_ids,
    )

    pair_ids_found = list(eta.pair.values)
    if len(pair_ids_found) == 0:
        return None

    uncalibrated = uncalibrated.sel(pair=pair_ids_found)
    uncalibrated_error = uncalibrated_error.sel(pair=pair_ids_found)
    info = info.sel(pair=pair_ids_found)

    ch_r_found, ch_t_found = _read_channels_from_info(info)

    calibrated = uncalibrated / eta
    calibrated_error = _division_error(
        numerator=uncalibrated,
        numerator_error=uncalibrated_error,
        denominator=eta,
        denominator_error=eta_error,
        quotient=calibrated,
    )

    G_R, G_T, H_R, H_T = _ideal_GH_for_pairs(
        ch_r=ch_r_found,
        ch_t=ch_t_found,
        pair_ids=pair_ids_found,
    )

    calibrated, calibrated_error = _GH_correct_ratio(
        ratio=calibrated,
        ratio_error=calibrated_error,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )

    info = _add_GH_to_info(
        info=info,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )

    return (
        calibrated.rename("ratio"),
        calibrated_error.rename("ratio_error"),
        info,
        ch_r_found,
        ch_t_found,
    )


def build_vldr_product(
    output_data: Dict[str, Dict[str, Any]],
    sig: xr.DataArray,
    sig_err: xr.DataArray,
    channel_info_qa: xr.DataArray,
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    eta_ids: Sequence[str],
    eta_info: xr.DataArray,
    average_time: bool = False,
) -> Optional[Tuple[xr.DataArray, xr.DataArray, xr.DataArray, List[str], List[str]]]:
    """Build VLDR from signal ratios, scalar eta, and real GH."""

    calibrated_ids = [
        replace_ratio_id_type(eta_id, old_type_index="e", new_type_index="d")
        for eta_id in eta_ids
    ]
    vldr_ids = [
        replace_ratio_id_type(calibrated_id, old_type_index="d", new_type_index="v")
        for calibrated_id in calibrated_ids
    ]

    cal_info = _pair_info_for_channels(
        pci=eta_info,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=calibrated_ids,
        ratio_type="calibrated_ratio",
    )

    uncalibrated, uncalibrated_error = ratio_from_channel_pairs(
        sig=sig,
        sig_err=sig_err,
        ch_r=ch_r,
        ch_t=ch_t,
        info=cal_info,
        average_time=average_time,
    )

    eta, eta_error = _find_eta_for_pairs(
        output_data=output_data,
        ch_r=ch_r,
        ch_t=ch_t,
        calibrated_pair_ids=calibrated_ids,
    )

    calibrated_ids_found = list(eta.pair.values)
    if len(calibrated_ids_found) == 0:
        return None

    uncalibrated = uncalibrated.sel(pair=calibrated_ids_found)
    uncalibrated_error = uncalibrated_error.sel(pair=calibrated_ids_found)
    cal_info = cal_info.sel(pair=calibrated_ids_found)

    ch_r_found, ch_t_found = _read_channels_from_info(cal_info)

    G_R, G_T, H_R, H_T = _channel_info_GH_for_pairs(
        channel_info_qa=channel_info_qa,
        ch_r=ch_r_found,
        ch_t=ch_t_found,
        pair_ids=calibrated_ids_found,
    )

    calibrated = uncalibrated / eta
    calibrated_error = _division_error(
        numerator=uncalibrated,
        numerator_error=uncalibrated_error,
        denominator=eta,
        denominator_error=eta_error,
        quotient=calibrated,
    )

    vldr, vldr_error = _GH_correct_ratio(
        ratio=calibrated,
        ratio_error=calibrated_error,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )

    vldr_ids_found = [
        replace_ratio_id_type(calibrated_id, old_type_index="d", new_type_index="v")
        for calibrated_id in calibrated_ids_found
    ]

    vldr = vldr.assign_coords(pair=vldr_ids_found).rename("ratio")
    vldr_error = vldr_error.assign_coords(pair=vldr_ids_found).rename("ratio_error")

    info = _pair_info_for_channels(
        pci=cal_info,
        ch_r=ch_r_found,
        ch_t=ch_t_found,
        pair_ids=vldr_ids_found,
        ratio_type="vldr",
    )

    info = _add_GH_to_info(
        info=info,
        G_R=G_R.assign_coords(pair=vldr_ids_found),
        G_T=G_T.assign_coords(pair=vldr_ids_found),
        H_R=H_R.assign_coords(pair=vldr_ids_found),
        H_T=H_T.assign_coords(pair=vldr_ids_found),
    )

    return vldr, vldr_error, info, ch_r_found, ch_t_found


def add_vldr_residual_to_info(
    processing_info: Dict[str, Any],
    output_data: Dict[str, Dict[str, Any]],
    info: xr.DataArray,
    vldr_mean: xr.DataArray,
    qa_test: str,
) -> xr.DataArray:
    """Add vldr_residual and sr_limit rows if matching MLDR info exists."""

    molecular_info = output_data.get("molecular_info", {})

    if qa_test not in molecular_info:
        CustomWarning(f"VLDR residual skipped: no molecular_info entry for {qa_test}.")
        return info

    mldr_info = molecular_info[qa_test]
    if "mean" not in mldr_info.parameters.values:
        CustomWarning(
            f"VLDR residual skipped: molecular_info[{qa_test!r}] has no mean row."
        )
        return info

    mldr_ids = [
        replace_ratio_id_type(vldr_id, old_type_index="v", new_type_index="m")
        for vldr_id in vldr_mean.pair.values
    ]

    mldr_mean = (
        mldr_info
        .sel(parameters="mean")
        .reindex(pair=mldr_ids)
        .assign_coords(pair=vldr_mean.pair.values)
        .astype("float64")
        .rename("mldr_mean")
    )

    residual = (vldr_mean - mldr_mean) / (1.0 - vldr_mean * mldr_mean)

    err_p = processing_info["settings_info"]["pcb"].get(
        "pldr_error_threshold",
        0.025,
    )

    _, _, _, sr_limit = pldr_error(
        delta_m=mldr_mean,
        delta_v_err=residual,
        delta_p_err_ulim=err_p,
    )
    
    residual = residual.reset_coords(drop=True)
    sr_limit = sr_limit.reset_coords(drop=True)

    return add_parameters(
        info,
        {
            "vldr_residual": residual,
            "sr_limit": sr_limit,
        },
    )


def compute_mean_calibrated_ratio(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Compute mean calibrated_ratio products and pol_cal_info metadata."""

    output_data = shallow_copy(input_data)
    vertical_scale = output_data[processing_info["caller_info"]["vertical_scale"]]
    averaging_range = processing_info["settings_info"]["pcb"]["rayleigh_region"]

    qa_alias, qa_store = get_ray_io_keys(
        output_data,
        processing_info,
        required_store="profile_mean",
    )
    if qa_alias is None:
        return output_data
    if qa_alias not in output_data["profile_error_mean"] or qa_alias not in vertical_scale:
        return output_data

    sig = output_data["profile_mean"][qa_alias]
    sig_err = output_data["profile_error_mean"][qa_alias]

    eta_inputs = collect_ray_eta_inputs(output_data, sig)
    if eta_inputs is None:
        print_entry("Mean calibrated ratio calculation skipped: no usable pcb/pcb_aux eta pairs found.")
        return output_data

    ch_r, ch_t, eta_ids, eta_info = eta_inputs
    product = build_calibrated_ratio_product(
        output_data=output_data,
        sig=sig,
        sig_err=sig_err,
        ch_r=ch_r,
        ch_t=ch_t,
        eta_ids=eta_ids,
        eta_info=eta_info,
        average_time=True,
    )
    if product is None:
        return output_data

    ratio, ratio_error, info, ch_r_found, _ = product
    z_pair = channels_to_pairs(vertical_scale[qa_alias].sel(channel=ch_r_found), info)

    info = add_region_stats_to_info(
        info=info,
        values=ratio,
        z_pair=z_pair,
        averaging_range=averaging_range,
        ratio_type="calibrated_ratio",
    )

    store_pol_cal_mean_product(output_data, qa_store, ratio, ratio_error, info)

    print_entry("Mean calibrated ratio calculation complete!")
    return output_data


def compute_calibrated_ratio(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Compute time-resolved calibrated_ratio products lazily."""

    output_data = shallow_copy(input_data)

    qa_alias, qa_store = get_ray_io_keys(
        output_data,
        processing_info,
        required_store="profile",
    )
    if qa_alias is None:
        return output_data
    if qa_alias not in output_data["profile_error"]:
        return output_data

    sig = output_data["profile"][qa_alias]
    sig_err = output_data["profile_error"][qa_alias]

    eta_inputs = collect_ray_eta_inputs(output_data, sig)
    if eta_inputs is None:
        print_entry("Calibrated ratio calculation skipped: no usable pcb/pcb_aux eta pairs found.")
        return output_data

    ch_r, ch_t, eta_ids, eta_info = eta_inputs
    product = build_calibrated_ratio_product(
        output_data=output_data,
        sig=sig,
        sig_err=sig_err,
        ch_r=ch_r,
        ch_t=ch_t,
        eta_ids=eta_ids,
        eta_info=eta_info,
        average_time=False,
    )
    if product is None:
        return output_data

    ratio, ratio_error, _, _, _ = product
    store_pol_cal_product(output_data, qa_store, ratio, ratio_error, mean=False)

    print_entry("Calibrated ratio calculation complete!")
    return output_data


def compute_mean_vldr(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Compute mean VLDR products and pol_cal_info metadata."""

    output_data = shallow_copy(input_data)
    vertical_scale = output_data[processing_info["caller_info"]["vertical_scale"]]
    averaging_range = processing_info["settings_info"]["pcb"]["rayleigh_region"]

    qa_alias, qa_store = get_ray_io_keys(
        output_data,
        processing_info,
        required_store="profile_mean",
    )
    if qa_alias is None:
        return output_data
    if qa_alias not in output_data["profile_error_mean"]:
        return output_data
    if qa_alias not in output_data["channel_info"]:
        print_entry(f"Mean VLDR calculation skipped: no {qa_alias} channel_info found.")
        return output_data
    if qa_alias not in vertical_scale:
        print_entry(f"Mean VLDR calculation skipped: no {qa_alias} vertical scale found.")
        return output_data

    sig = output_data["profile_mean"][qa_alias]
    sig_err = output_data["profile_error_mean"][qa_alias]

    eta_inputs = collect_ray_eta_inputs(output_data, sig)
    if eta_inputs is None:
        print_entry("Mean VLDR calculation skipped: no usable pcb/pcb_aux eta pairs found.")
        return output_data

    ch_r, ch_t, eta_ids, eta_info = eta_inputs
    product = build_vldr_product(
        output_data=output_data,
        sig=sig,
        sig_err=sig_err,
        channel_info_qa=output_data["channel_info"][qa_alias],
        ch_r=ch_r,
        ch_t=ch_t,
        eta_ids=eta_ids,
        eta_info=eta_info,
        average_time=True,
    )
    if product is None:
        return output_data

    ratio, ratio_error, info, ch_r_found, _ = product
    z_pair = channels_to_pairs(vertical_scale[qa_alias].sel(channel=ch_r_found), info)

    vldr_mean = mean_in_region(ratio, z_pair, averaging_range)
    vldr_sem = sem_in_region(ratio, z_pair, averaging_range)

    info = add_parameters(
        info,
        {
            "mean": vldr_mean,
            "sem": vldr_sem,
            "ratio_type": "vldr",
        },
    )
    info = add_vldr_residual_to_info(
        processing_info=processing_info,
        output_data=output_data,
        info=info,
        vldr_mean=vldr_mean,
        qa_test=qa_store,
    )

    store_pol_cal_mean_product(output_data, qa_store, ratio, ratio_error, info)

    print_entry("Mean VLDR calculation complete!")
    return output_data


def compute_vldr(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Compute time-resolved VLDR products lazily."""

    output_data = shallow_copy(input_data)

    qa_alias, qa_store = get_ray_io_keys(
        output_data,
        processing_info,
        required_store="profile",
    )
    if qa_alias is None:
        return output_data
    if qa_alias not in output_data["profile_error"]:
        return output_data
    if qa_alias not in output_data["channel_info"]:
        print_entry(f"VLDR calculation skipped: no {qa_alias} channel_info found.")
        return output_data

    sig = output_data["profile"][qa_alias]
    sig_err = output_data["profile_error"][qa_alias]

    eta_inputs = collect_ray_eta_inputs(output_data, sig)
    if eta_inputs is None:
        print_entry("VLDR calculation skipped: no usable pcb/pcb_aux eta pairs found.")
        return output_data

    ch_r, ch_t, eta_ids, eta_info = eta_inputs
    product = build_vldr_product(
        output_data=output_data,
        sig=sig,
        sig_err=sig_err,
        channel_info_qa=output_data["channel_info"][qa_alias],
        ch_r=ch_r,
        ch_t=ch_t,
        eta_ids=eta_ids,
        eta_info=eta_info,
        average_time=False,
    )
    if product is None:
        return output_data

    ratio, ratio_error, _, _, _ = product
    store_pol_cal_product(output_data, qa_store, ratio, ratio_error, mean=False)

    print_entry("VLDR calculation complete!")
    return output_data

def compute_mldr(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute MLDR for ray from molecular profiles.

    The MLDR channel combinations are created directly from the base ray
    polarization-calibration metadata, instead of depending on existing eta
    entries.  This makes MLDR independent of compute_calibration_factor while
    still producing pair ids that can be linked to VLDR by using the same
    channel-pair rules and type index "m".

    The molecular input is read from output_data["molecular"] or
    output_data["molec"] and is expected to have dims (channel, bins), after
    optional selection of opto_parameters="atten_bsc".

    The output is saved under:
        output_data["molecular_ratio"]["ray"]
        output_data["molecular_info"]["ray"]

    No random error is calculated or stored.  The regional mean is stored in
    molecular_info; no sem row is added.
    """

    output_data = shallow_copy(input_data)
    _ensure_molecular_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    qa_test_alias, qa_store = get_ray_io_keys(
        output_data=output_data,
        processing_info=processing_info,
        required_store="profile",
    )
    if qa_test_alias is None:
        return output_data

    molec_store = output_data.get("molecular", output_data.get("molec", {}))
    pol_cal_info = output_data.get("pol_cal_info", {})

    if qa_test_alias not in molec_store:
        print_entry("MLDR calculation skipped: no ray molecular profiles found.")
        return output_data
    if qa_test_alias not in vertical_scale:
        print_entry("MLDR calculation skipped: no ray vertical scale found.")
        return output_data
    qa_info_key = qa_store if qa_store in pol_cal_info else qa_test_alias
    if qa_info_key not in pol_cal_info:
        print_entry(f"MLDR calculation skipped: no {qa_store} pol_cal_info found.")
        return output_data

    molec = molec_store[qa_test_alias]
    if "opto_parameters" in molec.dims:
        molec = molec.sel(opto_parameters="atten_bsc")
    z = vertical_scale[qa_test_alias]

    if "channel" not in molec.dims:
        raise ValueError(
            f"MLDR molecular input for {qa_test_alias} must have a channel dimension. "
            f"Found dims: {molec.dims}."
        )

    pci_base = _base_pol_cal_pairs(pol_cal_info[qa_info_key])
    ch_r, ch_t = _read_channels_from_info(pci_base)

    available_channels = set(molec.channel.values.tolist())
    keep = [
        i for i, (r, t) in enumerate(zip(ch_r, ch_t))
        if r in available_channels and t in available_channels
    ]

    if len(keep) == 0:
        print_entry("MLDR calculation skipped: no base ray channel pairs found in molecular profiles.")
        return output_data

    ch_r = [ch_r[i] for i in keep]
    ch_t = [ch_t[i] for i in keep]
    pci_base = pci_base.isel(pair=keep)

    mldr_ids = [
        make_ratio_id(r, t, type_index="m")
        for r, t in zip(ch_r, ch_t)
    ]

    mldr_info = _pair_info_for_channels(
        pci=pci_base,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=mldr_ids,
        ratio_type="mldr",
    )

    molec_r = molec.sel(channel=ch_r)
    molec_t = molec.sel(channel=ch_t)

    delta_star_m = simple_ratio(
        numerator=molec_r,
        denominator=molec_t,
        info=mldr_info,
    )
    
    # The molecular ratio must be converted with the same ideal analyzer
    # response used for calibrated_ratio before comparing MLDR and VLDR.
    G_R, G_T, H_R, H_T = _ideal_GH_for_pairs(
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=mldr_ids,
    )
    
    eta = _ideal_eta_for_pairs(
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=mldr_ids,
    )

    # Correct for eta which is important for ratios that involve total channels
    delta_cor_m = delta_star_m / eta
    
    mldr = _GH_correct_values(
        ratio=delta_cor_m,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )
    
    mldr = mldr.rename("ratio")

    mldr_info = _add_GH_to_info(
        info=mldr_info,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )

    z_r = z.sel(channel=ch_r)
    z_pair = channels_to_pairs(z_r, mldr_info)

    mldr_m_bins = mean_in_region(
        da=mldr,
        z=z_pair,
        averaging_range=processing_info["settings_info"]["pcb"]["rayleigh_region"],
    )

    mldr_info = add_parameter(mldr_info, name="mean", values=mldr_m_bins)
    mldr_info = add_parameter(mldr_info, name="ratio_type", values="mldr")

    output_data["molecular_ratio"][qa_store] = append_or_replace_pairs(
        output_data["molecular_ratio"].get(qa_store),
        mldr,
    )
    output_data["molecular_info"][qa_store] = append_or_replace_info(
        output_data["molecular_info"].get(qa_store),
        mldr_info,
    )

    print_entry("MLDR calculation complete!")
    return output_data

