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
    da_error: xr.DataArray,
    z: xr.DataArray,
    averaging_range,
) -> xr.DataArray:
    """
    Average an error DataArray over bins and convert it to SEM.

    This assumes da_error is already the per-bin uncertainty of the ratio.
    """

    z_min, z_max = averaging_range
    mask = (z >= 1E3 * z_min) & (z <= 1E3 * z_max)

    da_sel = da_error.where(mask)
    n_bins = da_sel.notnull().sum("bins")

    return da_sel.mean("bins", skipna=True) / np.sqrt(n_bins)


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

    pci expected dims:
        parameters, pair

    values may be a pair DataArray, a sequence with one value per pair, or a scalar.
    String metadata such as ratio_type is preserved as object dtype.
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

    if hasattr(values.data, "compute"):
        values = values.compute()

    new_row = values.assign_coords(pair=pci.pair).expand_dims(parameters=[name])
    new_row = new_row.transpose("parameters", "pair").astype(object)

    if hasattr(pci.data, "compute"):
        pci = pci.compute()

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


def _time_mean_and_error(
    values: xr.DataArray,
    errors: xr.DataArray,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Average over time and propagate independent errors of the mean."""

    if "time" not in values.dims:
        return values, errors

    values_m_time = values.mean("time", skipna=True)
    n_time = values.count("time")
    errors_m_time = np.sqrt((errors**2).sum("time", skipna=True)) / n_time

    return values_m_time, errors_m_time


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


def _as_1d_list(values) -> List[Any]:
    """Convert scalar, numpy scalar, DataArray values, list, tuple to a 1D list."""

    if isinstance(values, xr.DataArray):
        values = values.values

    arr = np.atleast_1d(values)

    return arr.tolist()

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

    def _channel_to_pair(da: xr.DataArray, defaults: xr.DataArray, name: str) -> xr.DataArray:
        da = da.assign_coords(channel=list(pair_ids)).rename({"channel": "pair"})
        da = da.assign_coords(pair=list(pair_ids)).rename(name)
        return _finite_or_default(da, defaults.rename(name))

    return (
        _channel_to_pair(G_R, G_R_default, "G_R"),
        _channel_to_pair(G_T, G_T_default, "G_T"),
        _channel_to_pair(H_R, H_R_default, "H_R"),
        _channel_to_pair(H_T, H_T_default, "H_T"),
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
    """Find pcb/pcb_aux eta values and align them to calibrated_pair_ids.

    Only eta entries created under the combined calibration keys ("pcb" and
    "pcb_aux") are valid sources.  The returned eta arrays are time-free so
    they broadcast over time-resolved ray_pcb profiles without xarray aligning
    incompatible calibration and measurement time coordinates.
    """

    eta_pair_ids = [
        make_ratio_id(r, t, type_index="e")
        for r, t in zip(ch_r, ch_t)
    ]

    eta_sources = ["pcb", "pcb_aux"]
    eta_values = []
    eta_errors = []

    for eta_id, calibrated_id in zip(eta_pair_ids, calibrated_pair_ids):
        eta_da = None
        eta_err_da = None

        for source in eta_sources:
            ratio_store = output_data.get("pol_cal_ratio", {})
            error_store = output_data.get("pol_cal_ratio_error", {})

            if source not in ratio_store or source not in error_store:
                continue

            candidate = ratio_store[source]
            candidate_err = error_store[source]

            if "ratio" in candidate.dims:
                candidate = candidate.rename({"ratio": "pair"})
            if "ratio" in candidate_err.dims:
                candidate_err = candidate_err.rename({"ratio": "pair"})

            if "pair" not in candidate.dims or eta_id not in candidate.pair.values:
                continue

            eta_da = candidate.sel(pair=eta_id)
            eta_err_da = candidate_err.sel(pair=eta_id)

            # Calibration products must not impose their calibration time
            # coordinate on time-resolved ray_pcb products.
            if "time" in eta_da.dims:
                if eta_da.sizes.get("time", 0) == 1:
                    eta_da = eta_da.isel(time=0, drop=True)
                    eta_err_da = eta_err_da.isel(time=0, drop=True)
                else:
                    eta_da = eta_da.mean("time", skipna=True)
                    eta_err_da = np.sqrt((eta_err_da**2).sum("time", skipna=True)) / eta_err_da.count("time")

            break

        if eta_da is None or eta_err_da is None:
            continue

        eta_values.append(eta_da.expand_dims(pair=[calibrated_id]))
        eta_errors.append(eta_err_da.expand_dims(pair=[calibrated_id]))

    if not eta_values:
        raise KeyError("No matching eta entries found in pcb/pcb_aux for ray_pcb channel pairs.")

    eta = xr.concat(eta_values, dim="pair")
    eta_error = xr.concat(eta_errors, dim="pair")

    preferred = [dim for dim in ["pair", "bins"] if dim in eta.dims]
    other = [dim for dim in eta.dims if dim not in preferred]
    eta = eta.transpose(*preferred, *other)
    eta_error = eta_error.transpose(*preferred, *other)

    return eta, eta_error



def _collect_eta_pairs_from_pcb(
    output_data: Dict[str, Dict[str, Any]],
) -> Tuple[List[str], List[str], List[str], xr.DataArray]:
    """Collect eta channel combinations from pcb/pcb_aux pol_cal_info.

    Returns ch_r, ch_t, eta_pair_ids, eta_info.  These eta entries are the
    authoritative source of the channel combinations used by ray_pcb calibrated
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
    Calculate PLDR error lookup for each pair.

    delta_m and delta_v_err must be pair-based DataArrays. Missing pairs are
    retained through outer alignment and produce NaNs in the output.
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
    min_bsc_values = []

    for pair in delta_p_err.pair.values:
        row = last_delta_p_err.sel(pair=pair)

        if row.isnull().all():
            min_bsc_values.append(np.nan)
            continue

        dv = float(delta_v_err.sel(pair=pair).values)

        if dv > 0.0001:
            idx = int(row.argmax(dim="R", skipna=True).values)
            min_bsc_values.append(float(R.isel(R=idx).values))
        elif dv < -0.0001:
            idx = int(row.argmin(dim="R", skipna=True).values)
            min_bsc_values.append(float(R.isel(R=idx).values))
        else:
            min_bsc_values.append(1.01)

    min_bsc_ratio = xr.DataArray(
        min_bsc_values,
        dims=["pair"],
        coords={"pair": delta_p_err.pair.values},
        name="sr_limit",
    )

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


def compute_gain_ratio(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute gain ratios for pcb/pcb_aux +/-45 QA tests.

    The reflected and transmitted signals are averaged over time before the
    ratio is formed.  The saved pol_cal_ratio and pol_cal_ratio_error products
    are time-free arrays with dimensions (pair, bins).

    In addition to the individual +/-45 gain ratios, this function now also
    creates the combined calibration QA tests:
        pcb      from pcb_p45 / pcb_m45
        pcb_aux  from pcb_aux_p45 / pcb_aux_m45

    The combined entries are stored as ratio_type="gain_ratio" with g IDs and
    correspond to eta_s = sqrt(gain_ratio_p45 * gain_ratio_m45).  The
    calibration-factor stage can therefore start directly from pcb/pcb_aux
    gain_ratio products.
    """

    output_data = shallow_copy(input_data)
    _ensure_pol_cal_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    profiles = output_data["profile"]
    profile_errors = output_data["profile_error"]
    channel_info = output_data["channel_info"]
    pol_cal_info = output_data["pol_cal_info"]

    allowed_qa_tests = [
        "pcb_p45",
        "pcb_m45",
        "pcb_aux_p45",
        "pcb_aux_m45",
    ]

    averaging_range = processing_info["settings_info"]["pcb"]["calibration_region"]

    for qa_test in allowed_qa_tests:
        if qa_test not in profiles:
            continue
        if qa_test not in channel_info:
            continue
        if qa_test not in pol_cal_info:
            continue
        if qa_test not in vertical_scale:
            continue
        if qa_test not in profile_errors:
            continue

        z = vertical_scale[qa_test]
        sig = profiles[qa_test]
        sig_err = profile_errors[qa_test]
        pci_base = _base_pol_cal_pairs(output_data["pol_cal_info"][qa_test])

        ch_r, ch_t = _read_channels_from_info(pci_base)
        pair_ids = [
            make_ratio_id(numerator, denominator, type_index="g")
            for numerator, denominator in zip(ch_r, ch_t)
        ]

        if len(pair_ids) == 0:
            continue

        pci = _pair_info_for_channels(
            pci=pci_base,
            ch_r=ch_r,
            ch_t=ch_t,
            pair_ids=pair_ids,
            ratio_type="gain_ratio",
        )

        sig_r = sig.sel(channel=ch_r)
        sig_t = sig.sel(channel=ch_t)
        sig_r_err = sig_err.sel(channel=ch_r)
        sig_t_err = sig_err.sel(channel=ch_t)

        # Average signals first, then form the ratio. This improves the SNR of
        # photon-counting channels compared with averaging noisy ratios.
        sig_r_m_time, sig_r_err_m_time = _time_mean_and_error(sig_r, sig_r_err)
        sig_t_m_time, sig_t_err_m_time = _time_mean_and_error(sig_t, sig_t_err)

        ratio_da = simple_ratio(
            numerator=sig_r_m_time,
            denominator=sig_t_m_time,
            info=pci,
        )

        ratio_error_da = ratio_error_independent(
            numerator=sig_r_m_time,
            denominator=sig_t_m_time,
            numerator_error=sig_r_err_m_time,
            denominator_error=sig_t_err_m_time,
            ratio_values=ratio_da,
            info=pci,
        )

        z_r = z.sel(channel=ch_r)
        z_pair = channels_to_pairs(z_r, pci)

        ratio_m_bins = mean_in_region(
            da=ratio_da,
            z=z_pair,
            averaging_range=averaging_range,
        )

        ratio_error_m_bins = sem_in_region(
            da_error=ratio_error_da,
            z=z_pair,
            averaging_range=averaging_range,
        )

        pci = add_parameter(pci=pci, name="mean", values=ratio_m_bins)
        pci = add_parameter(pci=pci, name="sem", values=ratio_error_m_bins)
        pci = add_parameter(pci=pci, name="ratio_type", values="gain_ratio")

        output_data["pol_cal_ratio"][qa_test] = append_or_replace_pairs(
            output_data["pol_cal_ratio"].get(qa_test),
            ratio_da.rename("ratio"),
        )
        output_data["pol_cal_ratio_error"][qa_test] = append_or_replace_pairs(
            output_data["pol_cal_ratio_error"].get(qa_test),
            ratio_error_da.rename("ratio_error"),
        )
        output_data["pol_cal_info"][qa_test] = append_or_replace_info(
            output_data["pol_cal_info"].get(qa_test),
            pci,
        )

    # Build combined pcb / pcb_aux gain-ratio products here, so downstream
    # calibration-factor calculation only needs the combined QA test.
    assign_pcb = {
        "pcb": ("pcb_p45", "pcb_m45"),
        "pcb_aux": ("pcb_aux_p45", "pcb_aux_m45"),
    }

    for target_key, (p45, m45) in assign_pcb.items():
        if p45 not in output_data["pol_cal_ratio"] or m45 not in output_data["pol_cal_ratio"]:
            continue
        if p45 not in output_data["pol_cal_ratio_error"] or m45 not in output_data["pol_cal_ratio_error"]:
            continue
        if p45 not in output_data["pol_cal_info"] or m45 not in output_data["pol_cal_info"]:
            continue
        if p45 not in vertical_scale or m45 not in vertical_scale:
            continue

        pci_p45 = output_data["pol_cal_info"][p45]
        pci_m45 = output_data["pol_cal_info"][m45]

        p45_gain_ids = _select_pairs_by_type(pci_p45, "gain_ratio")
        m45_gain_ids = _select_pairs_by_type(pci_m45, "gain_ratio")
        gain_ids = [pid for pid in p45_gain_ids if pid in m45_gain_ids]

        if len(gain_ids) == 0:
            continue

        ratio_p45 = output_data["pol_cal_ratio"][p45].sel(pair=gain_ids)
        ratio_p45_err = output_data["pol_cal_ratio_error"][p45].sel(pair=gain_ids)
        ratio_m45 = output_data["pol_cal_ratio"][m45].sel(pair=gain_ids)
        ratio_m45_err = output_data["pol_cal_ratio_error"][m45].sel(pair=gain_ids)

        eta_s_p45, eta_s_p45_error = _time_mean_and_error(
            ratio_p45,
            ratio_p45_err,
        )
        eta_s_m45, eta_s_m45_error = _time_mean_and_error(
            ratio_m45,
            ratio_m45_err,
        )

        pci_gain = pci_p45.sel(pair=gain_ids)
        pci_gain_m45 = pci_m45.sel(pair=gain_ids)

        ch_r_p45, ch_t_p45 = _read_channels_from_info(pci_gain)
        ch_r_m45, ch_t_m45 = _read_channels_from_info(pci_gain_m45)

        if ch_r_p45 != ch_r_m45 or ch_t_p45 != ch_t_m45:
            raise ValueError(
                f"Gain-ratio channel mismatch between {p45} and {m45} "
                f"for pairs {gain_ids}."
            )

        ch_r, ch_t = ch_r_p45, ch_t_p45

        z_p45 = vertical_scale[p45].sel(channel=ch_r)
        z_m45 = vertical_scale[m45].sel(channel=ch_r)

        z_pair_p45 = channels_to_pairs(z_p45, pci_gain)
        z_pair_m45 = channels_to_pairs(z_m45, pci_gain_m45)

        if not z_pair_p45.broadcast_equals(z_pair_m45):
            CustomWarning(
                f"Skipping combined gain ratio for {target_key}: "
                f"vertical scales of {p45} and {m45} do not match."
            )
            continue

        eta_s_product = eta_s_p45 * eta_s_m45
        eta_s = np.sqrt(eta_s_product.where(eta_s_product >= 0))
        eta_s_error = _product_sqrt_error(
            eta_s_p45,
            eta_s_p45_error,
            eta_s_m45,
            eta_s_m45_error,
            eta_s,
        )

        eta_s = eta_s.assign_coords(pair=gain_ids).rename("ratio")
        eta_s_error = eta_s_error.assign_coords(pair=gain_ids).rename("ratio_error")

        combined_info = _pair_info_for_channels(
            pci=pci_gain,
            ch_r=ch_r,
            ch_t=ch_t,
            pair_ids=gain_ids,
            ratio_type="gain_ratio",
        )

        eta_s_m_bins = mean_in_region(
            da=eta_s,
            z=z_pair_p45,
            averaging_range=averaging_range,
        )
        eta_s_error_m_bins = sem_in_region(
            da_error=eta_s_error,
            z=z_pair_p45,
            averaging_range=averaging_range,
        )

        # Epsilon belongs to the combined pcb / pcb_aux gain-ratio entries.
        # It is calculated from the +45 and -45 gain ratios, treating them as
        # eta_p45 and eta_m45, and only its calibration-region mean is stored
        # in pol_cal_info. No epsilon profile is stored.
        epsilon = epsilon_angle(
            eta_p45=eta_s_p45.assign_coords(pair=gain_ids),
            eta_m45=eta_s_m45.assign_coords(pair=gain_ids),
            kappa=1.0,
        )
        epsilon_m_bins = mean_in_region(
            da=epsilon,
            z=z_pair_p45,
            averaging_range=averaging_range,
        )

        combined_info = add_parameter(combined_info, name="mean", values=eta_s_m_bins)
        combined_info = add_parameter(combined_info, name="sem", values=eta_s_error_m_bins)
        combined_info = add_parameter(combined_info, name="epsilon", values=epsilon_m_bins)
        combined_info = add_parameter(combined_info, name="ratio_type", values="gain_ratio")

        output_data["pol_cal_ratio"][target_key] = append_or_replace_pairs(
            output_data["pol_cal_ratio"].get(target_key),
            eta_s,
        )
        output_data["pol_cal_ratio_error"][target_key] = append_or_replace_pairs(
            output_data["pol_cal_ratio_error"].get(target_key),
            eta_s_error,
        )
        output_data["pol_cal_info"][target_key] = append_or_replace_info(
            output_data["pol_cal_info"].get(target_key),
            combined_info,
        )

    print_entry("Gain ratio calculation complete!")
    return output_data


def compute_calibration_factor(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute calibration products from combined pcb/pcb_aux gain ratios.

    This stage now expects compute_gain_ratio() to have already created the
    combined calibration QA tests:
        pcb
        pcb_aux

    The existing ratio_type="gain_ratio" entries in those keys are loaded as
    eta_s.  The function then continues from that point and stores:
        eta_s_f : eta_s corrected by R_to_T_transmission_ratio, using f IDs
        eta     : eta_s_f corrected by K, using e IDs

    The original +/-45 products are not modified here.
    """

    target_keys = ["pcb", "pcb_aux"]

    output_data = shallow_copy(input_data)
    _ensure_pol_cal_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    pol_cal_ratio = output_data["pol_cal_ratio"]
    pol_cal_ratio_error = output_data["pol_cal_ratio_error"]
    pol_cal_info = output_data["pol_cal_info"]

    averaging_range = processing_info["settings_info"]["pcb"]["calibration_region"]

    # Use the original +/-45 vertical scales only as a source of the bin-height
    # grid for region statistics.  The calibration quantities themselves are
    # read only from the combined pcb / pcb_aux QA-test keys.
    vertical_reference = {
        "pcb": "pcb_p45",
        "pcb_aux": "pcb_aux_p45",
    }

    def _add_region_stats(
        info: xr.DataArray,
        values: xr.DataArray,
        errors: xr.DataArray,
        z_pair: xr.DataArray,
        ratio_type: str,
    ) -> xr.DataArray:
        """Attach mean, sem, and ratio_type rows to a pol_cal_info object."""

        values_m_bins = mean_in_region(
            da=values,
            z=z_pair,
            averaging_range=averaging_range,
        )
        errors_m_bins = sem_in_region(
            da_error=errors,
            z=z_pair,
            averaging_range=averaging_range,
        )

        info = add_parameter(info, name="mean", values=values_m_bins)
        info = add_parameter(info, name="sem", values=errors_m_bins)
        info = add_parameter(info, name="ratio_type", values=ratio_type)

        return info

    def _store_product(
        key: str,
        values: xr.DataArray,
        errors: xr.DataArray,
        info: xr.DataArray,
    ) -> None:
        """Append or replace one product family in the output dictionaries."""

        output_data["pol_cal_ratio"][key] = append_or_replace_pairs(
            output_data["pol_cal_ratio"].get(key),
            values,
        )
        output_data["pol_cal_ratio_error"][key] = append_or_replace_pairs(
            output_data["pol_cal_ratio_error"].get(key),
            errors,
        )
        output_data["pol_cal_info"][key] = append_or_replace_info(
            output_data["pol_cal_info"].get(key),
            info,
        )

    for target_key in target_keys:
        if target_key not in pol_cal_ratio:
            continue
        if target_key not in pol_cal_ratio_error:
            continue
        if target_key not in pol_cal_info:
            continue

        pci_gain_store = pol_cal_info[target_key]
        gain_ids = _select_pairs_by_type(pci_gain_store, "gain_ratio")

        if len(gain_ids) == 0:
            continue

        eta_s = pol_cal_ratio[target_key].sel(pair=gain_ids).rename("ratio")
        eta_s_error = pol_cal_ratio_error[target_key].sel(pair=gain_ids).rename("ratio_error")
        pci_gain = pci_gain_store.sel(pair=gain_ids)

        # Defensive time handling for older intermediate files.  New combined
        # gain-ratio products are time-free.
        eta_s, eta_s_error = _time_mean_and_error(eta_s, eta_s_error)

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
        )
        K = _numeric_parameter(
            pci=pci_gain,
            parameter="K",
            pair_ids=gain_ids,
            default=1.0,
        )

        trans_ratio = trans_ratio.assign_coords(pair=gain_ids)
        K = K.assign_coords(pair=gain_ids)

        eta_s_f = eta_s / trans_ratio
        eta_s_f_error = eta_s_error / np.abs(trans_ratio)

        eta = eta_s_f / K
        eta_error = eta_s_f_error / np.abs(K)

        product_specs = [
            {
                "ratio_type": "eta_s_f",
                "pair_ids": eta_s_f_ids,
                "values": eta_s_f,
                "errors": eta_s_f_error,
            },
            {
                "ratio_type": "eta",
                "pair_ids": eta_ids,
                "values": eta,
                "errors": eta_error,
            },
        ]

        z_pair_gain = None
        z_source = vertical_reference.get(target_key)
        if z_source in vertical_scale:
            z_r = vertical_scale[z_source].sel(channel=ch_r)
            z_pair_gain = channels_to_pairs(z_r, pci_gain)
        else:
            CustomWarning(
                f"No vertical scale reference found for {target_key}; "
                "eta_s_f/eta mean and sem will not be added."
            )

        for spec in product_specs:
            ratio_type = spec["ratio_type"]
            pair_ids = spec["pair_ids"]

            values_for_stats = spec["values"].assign_coords(pair=pair_ids).rename("ratio")
            errors_for_stats = spec["errors"].assign_coords(pair=pair_ids).rename("ratio_error")

            info = _pair_info_for_channels(
                pci=pci_gain,
                ch_r=ch_r,
                ch_t=ch_t,
                pair_ids=pair_ids,
                ratio_type=ratio_type,
            )

            if z_pair_gain is not None:
                z_pair = z_pair_gain.assign_coords(pair=pair_ids)
                info = _add_region_stats(
                    info=info,
                    values=values_for_stats,
                    errors=errors_for_stats,
                    z_pair=z_pair,
                    ratio_type=ratio_type,
                )
            else:
                info = add_parameter(info, name="ratio_type", values=ratio_type)

            _store_product(
                key=target_key,
                values=values_for_stats,
                errors=errors_for_stats,
                info=info,
            )

    print_entry("Calibration factor calculation complete!")

    return output_data

def compute_calibrated_ratio(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute calibrated ratios for ray_pcb.

    Valid channel combinations are taken only from existing eta entries under
    the combined pcb/pcb_aux calibration keys.  If no eta products exist, this
    function leaves the input unchanged.

    Two calculation branches are used:
        1. time-resolved branch, saved in pol_cal_ratio/pol_cal_ratio_error;
        2. signal-time-averaged branch, used only for pol_cal_info mean/sem.

    The metadata branch averages reflected/transmitted signals first, then
    forms the ratio and applies the same eta and ideal-GH correction.  This
    avoids averaging noisy ratios over time.
    """

    output_data = shallow_copy(input_data)
    _ensure_pol_cal_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    qa_test = "ray_pcb"

    profiles = output_data["profile"]
    profile_errors = output_data["profile_error"]

    if qa_test not in profiles:
        return output_data
    if qa_test not in vertical_scale:
        return output_data
    if qa_test not in profile_errors:
        return output_data

    ch_r, ch_t, eta_ids, eta_info = _collect_eta_pairs_from_pcb(output_data)
    if len(eta_ids) == 0:
        print_entry("Calibrated ratio calculation skipped: no pcb/pcb_aux eta entries found.")
        return output_data

    sig = profiles[qa_test]
    sig_err = profile_errors[qa_test]
    z = vertical_scale[qa_test]

    available_channels = set(sig.channel.values.tolist())
    keep = [
        i for i, (r, t) in enumerate(zip(ch_r, ch_t))
        if r in available_channels and t in available_channels
    ]

    if len(keep) == 0:
        print_entry("Calibrated ratio calculation skipped: no eta channel pairs found in ray_pcb profiles.")
        return output_data

    ch_r = [ch_r[i] for i in keep]
    ch_t = [ch_t[i] for i in keep]
    eta_ids = [eta_ids[i] for i in keep]
    eta_info = eta_info.isel(pair=keep)

    pair_ids = [
        replace_ratio_id_type(eta_id, old_type_index="e", new_type_index="d")
        for eta_id in eta_ids
    ]

    pci = _pair_info_for_channels(
        pci=eta_info,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=pair_ids,
        ratio_type="calibrated_ratio",
    )

    sig_r = sig.sel(channel=ch_r)
    sig_t = sig.sel(channel=ch_t)
    sig_r_err = sig_err.sel(channel=ch_r)
    sig_t_err = sig_err.sel(channel=ch_t)

    # ------------------------------------------------------------------
    # 1) Time-resolved product branch: store these profiles unchanged in
    #    pol_cal_ratio/pol_cal_ratio_error.
    # ------------------------------------------------------------------
    uncalibrated = simple_ratio(
        numerator=sig_r,
        denominator=sig_t,
        info=pci,
    )

    uncalibrated_error = ratio_error_independent(
        numerator=sig_r,
        denominator=sig_t,
        numerator_error=sig_r_err,
        denominator_error=sig_t_err,
        ratio_values=uncalibrated,
        info=pci,
    )

    eta, eta_error = _find_eta_for_pairs(
        output_data=output_data,
        ch_r=ch_r,
        ch_t=ch_t,
        calibrated_pair_ids=pair_ids,
    )

    # Keep only pairs for which eta was found.
    pair_ids_found = list(eta.pair.values)
    uncalibrated = uncalibrated.sel(pair=pair_ids_found)
    uncalibrated_error = uncalibrated_error.sel(pair=pair_ids_found)
    pci = pci.sel(pair=pair_ids_found)

    ch_r_found = pci.sel(parameters="ch_r").values.tolist()
    ch_t_found = pci.sel(parameters="ch_t").values.tolist()

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

    calibrated = calibrated.rename("ratio")
    calibrated_error = calibrated_error.rename("ratio_error")

    pci = _add_GH_to_info(
        info=pci,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )

    # ------------------------------------------------------------------
    # 2) Metadata branch: average the signals first, then form/correct the
    #    ratio and finally average over bins.  This branch is pair-only and
    #    is stored only in pol_cal_info.
    # ------------------------------------------------------------------
    sig_r_m_time, sig_r_err_m_time = _time_mean_and_error(
        sig_r.sel(channel=ch_r_found),
        sig_r_err.sel(channel=ch_r_found),
    )
    sig_t_m_time, sig_t_err_m_time = _time_mean_and_error(
        sig_t.sel(channel=ch_t_found),
        sig_t_err.sel(channel=ch_t_found),
    )

    uncalibrated_m_time = simple_ratio(
        numerator=sig_r_m_time,
        denominator=sig_t_m_time,
        info=pci,
    )
    uncalibrated_error_m_time = ratio_error_independent(
        numerator=sig_r_m_time,
        denominator=sig_t_m_time,
        numerator_error=sig_r_err_m_time,
        denominator_error=sig_t_err_m_time,
        ratio_values=uncalibrated_m_time,
        info=pci,
    )

    eta_m_time = eta.sel(pair=pair_ids_found)
    eta_error_m_time = eta_error.sel(pair=pair_ids_found)

    calibrated_m_time = uncalibrated_m_time / eta_m_time
    calibrated_error_m_time = _division_error(
        numerator=uncalibrated_m_time,
        numerator_error=uncalibrated_error_m_time,
        denominator=eta_m_time,
        denominator_error=eta_error_m_time,
        quotient=calibrated_m_time,
    )

    calibrated_m_time, calibrated_error_m_time = _GH_correct_ratio(
        ratio=calibrated_m_time,
        ratio_error=calibrated_error_m_time,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )

    z_r = z.sel(channel=ch_r_found)
    z_pair = channels_to_pairs(z_r, pci)

    calibrated_m_bins = mean_in_region(
        da=calibrated_m_time,
        z=z_pair,
        averaging_range=processing_info["settings_info"]["pcb"]["rayleigh_region"],
    )
    calibrated_error_m_bins = sem_in_region(
        da_error=calibrated_error_m_time,
        z=z_pair,
        averaging_range=processing_info["settings_info"]["pcb"]["rayleigh_region"],
    )

    pci = add_parameter(pci=pci, name="mean", values=calibrated_m_bins)
    pci = add_parameter(pci=pci, name="sem", values=calibrated_error_m_bins)
    pci = add_parameter(pci=pci, name="ratio_type", values="calibrated_ratio")

    output_data["pol_cal_ratio"][qa_test] = append_or_replace_pairs(
        output_data["pol_cal_ratio"].get(qa_test),
        calibrated,
    )
    output_data["pol_cal_ratio_error"][qa_test] = append_or_replace_pairs(
        output_data["pol_cal_ratio_error"].get(qa_test),
        calibrated_error,
    )
    output_data["pol_cal_info"][qa_test] = append_or_replace_info(
        output_data["pol_cal_info"].get(qa_test),
        pci,
    )

    print_entry("Calibrated ratio calculation complete!")
    return output_data


def compute_vldr(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute VLDR for ray_pcb from time-resolved calibrated ratios.

    Two calculation branches are used:
        1. time-resolved branch from the existing calibrated_ratio entries,
           saved in pol_cal_ratio/pol_cal_ratio_error;
        2. signal-time-averaged branch, recomputed from ray_pcb signals using
           eta and both GH corrections, used only for pol_cal_info mean/sem,
           vldr_residual, and sr_limit.

    The metadata branch averages reflected/transmitted signals first, then
    forms ratios and averages over bins.  No time averaging of an already
    formed ratio/VLDR is used for pol_cal_info.
    """

    output_data = shallow_copy(input_data)
    _ensure_pol_cal_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    qa_test = "ray_pcb"

    if qa_test not in output_data.get("pol_cal_ratio", {}):
        print_entry(f"VLDR calculation skipped: no {qa_test} calibrated ratios found.")
        return output_data
    if qa_test not in output_data.get("pol_cal_ratio_error", {}):
        print_entry(f"VLDR calculation skipped: no {qa_test} calibrated ratio errors found.")
        return output_data
    if qa_test not in output_data.get("pol_cal_info", {}):
        print_entry(f"VLDR calculation skipped: no {qa_test} pol_cal_info found.")
        return output_data
    if qa_test not in vertical_scale:
        print_entry(f"VLDR calculation skipped: no {qa_test} vertical scale found.")
        return output_data
    if qa_test not in output_data.get("channel_info", {}):
        print_entry(f"VLDR calculation skipped: no {qa_test} channel_info found.")
        return output_data
    if qa_test not in output_data.get("profile", {}):
        print_entry(f"VLDR calculation skipped: no {qa_test} profiles found.")
        return output_data
    if qa_test not in output_data.get("profile_error", {}):
        print_entry(f"VLDR calculation skipped: no {qa_test} profile errors found.")
        return output_data

    ratio_store = output_data["pol_cal_ratio"][qa_test]
    error_store = output_data["pol_cal_ratio_error"][qa_test]
    pci_store = output_data["pol_cal_info"][qa_test]

    calibrated_ids = _select_pairs_by_type(pci_store, "calibrated_ratio")
    if len(calibrated_ids) == 0:
        print_entry("VLDR calculation skipped: no calibrated_ratio entries found.")
        return output_data

    calibrated = ratio_store.sel(pair=calibrated_ids)
    calibrated_error = error_store.sel(pair=calibrated_ids)
    pci_cal = pci_store.sel(pair=calibrated_ids)

    ch_r, ch_t = _read_channels_from_info(pci_cal)
    vldr_ids = [
        replace_ratio_id_type(calibrated_id, old_type_index="d", new_type_index="v")
        for calibrated_id in calibrated_ids
    ]

    vldr_info = _pair_info_for_channels(
        pci=pci_cal,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=vldr_ids,
        ratio_type="vldr",
    )

    channel_info_qa = output_data["channel_info"][qa_test]
    G_R, G_T, H_R, H_T = _channel_info_GH_for_pairs(
        channel_info_qa=channel_info_qa,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=calibrated_ids,
    )

    # ------------------------------------------------------------------
    # 1) Time-resolved product branch: apply the real GH correction to the
    #    existing time-resolved calibrated_ratio profiles.
    # ------------------------------------------------------------------
    vldr, vldr_error = _GH_correct_ratio(
        ratio=calibrated,
        ratio_error=calibrated_error,
        G_R=G_R,
        G_T=G_T,
        H_R=H_R,
        H_T=H_T,
    )

    vldr = vldr.assign_coords(pair=vldr_ids).rename("ratio")
    vldr_error = vldr_error.assign_coords(pair=vldr_ids).rename("ratio_error")

    vldr_info = _add_GH_to_info(
        info=vldr_info,
        G_R=G_R.assign_coords(pair=vldr_ids),
        G_T=G_T.assign_coords(pair=vldr_ids),
        H_R=H_R.assign_coords(pair=vldr_ids),
        H_T=H_T.assign_coords(pair=vldr_ids),
    )

    # ------------------------------------------------------------------
    # 2) Metadata branch: recompute from signal-time-averaged profiles.
    #    This mirrors compute_calibrated_ratio's metadata branch and then
    #    applies the real GH correction for VLDR.
    # ------------------------------------------------------------------
    sig = output_data["profile"][qa_test]
    sig_err = output_data["profile_error"][qa_test]
    z = vertical_scale[qa_test]

    sig_r = sig.sel(channel=ch_r)
    sig_t = sig.sel(channel=ch_t)
    sig_r_err = sig_err.sel(channel=ch_r)
    sig_t_err = sig_err.sel(channel=ch_t)

    sig_r_m_time, sig_r_err_m_time = _time_mean_and_error(sig_r, sig_r_err)
    sig_t_m_time, sig_t_err_m_time = _time_mean_and_error(sig_t, sig_t_err)

    uncalibrated_m_time = simple_ratio(
        numerator=sig_r_m_time,
        denominator=sig_t_m_time,
        info=pci_cal,
    )
    uncalibrated_error_m_time = ratio_error_independent(
        numerator=sig_r_m_time,
        denominator=sig_t_m_time,
        numerator_error=sig_r_err_m_time,
        denominator_error=sig_t_err_m_time,
        ratio_values=uncalibrated_m_time,
        info=pci_cal,
    )

    eta_m_time, eta_error_m_time = _find_eta_for_pairs(
        output_data=output_data,
        ch_r=ch_r,
        ch_t=ch_t,
        calibrated_pair_ids=calibrated_ids,
    )

    # Align the metadata branch to the calibrated ids available in eta.
    stats_pair_ids = list(eta_m_time.pair.values)
    uncalibrated_m_time = uncalibrated_m_time.sel(pair=stats_pair_ids)
    uncalibrated_error_m_time = uncalibrated_error_m_time.sel(pair=stats_pair_ids)

    pci_cal_stats = pci_cal.sel(pair=stats_pair_ids)
    ch_r_stats, ch_t_stats = _read_channels_from_info(pci_cal_stats)

    G_R_ideal, G_T_ideal, H_R_ideal, H_T_ideal = _ideal_GH_for_pairs(
        ch_r=ch_r_stats,
        ch_t=ch_t_stats,
        pair_ids=stats_pair_ids,
    )
    G_R_stats = G_R.sel(pair=stats_pair_ids)
    G_T_stats = G_T.sel(pair=stats_pair_ids)
    H_R_stats = H_R.sel(pair=stats_pair_ids)
    H_T_stats = H_T.sel(pair=stats_pair_ids)

    calibrated_m_time = uncalibrated_m_time / eta_m_time
    calibrated_error_m_time = _division_error(
        numerator=uncalibrated_m_time,
        numerator_error=uncalibrated_error_m_time,
        denominator=eta_m_time,
        denominator_error=eta_error_m_time,
        quotient=calibrated_m_time,
    )

    calibrated_m_time, calibrated_error_m_time = _GH_correct_ratio(
        ratio=calibrated_m_time,
        ratio_error=calibrated_error_m_time,
        G_R=G_R_ideal,
        G_T=G_T_ideal,
        H_R=H_R_ideal,
        H_T=H_T_ideal,
    )

    vldr_m_time, vldr_error_m_time = _GH_correct_ratio(
        ratio=calibrated_m_time,
        ratio_error=calibrated_error_m_time,
        G_R=G_R_stats,
        G_T=G_T_stats,
        H_R=H_R_stats,
        H_T=H_T_stats,
    )

    stats_vldr_ids = [
        replace_ratio_id_type(calibrated_id, old_type_index="d", new_type_index="v")
        for calibrated_id in stats_pair_ids
    ]
    vldr_m_time = vldr_m_time.assign_coords(pair=stats_vldr_ids).rename("ratio")
    vldr_error_m_time = vldr_error_m_time.assign_coords(pair=stats_vldr_ids).rename("ratio_error")

    z_r = z.sel(channel=ch_r_stats)
    z_pair_stats = channels_to_pairs(
        z_r,
        vldr_info.sel(pair=stats_vldr_ids),
    )

    vldr_m_bins_store = mean_in_region(
        da=vldr_m_time,
        z=z_pair_stats,
        averaging_range=processing_info["settings_info"]["pcb"]["rayleigh_region"],
    )
    vldr_error_m_bins_store = sem_in_region(
        da_error=vldr_error_m_time,
        z=z_pair_stats,
        averaging_range=processing_info["settings_info"]["pcb"]["rayleigh_region"],
    )

    # Reindex to the full VLDR id list so missing metadata pairs become NaN.
    vldr_m_bins_store = vldr_m_bins_store.reindex(pair=vldr_ids)
    vldr_error_m_bins_store = vldr_error_m_bins_store.reindex(pair=vldr_ids)

    vldr_info = add_parameter(vldr_info, name="mean", values=vldr_m_bins_store)
    vldr_info = add_parameter(vldr_info, name="sem", values=vldr_error_m_bins_store)

    mldr_ids = [
        replace_ratio_id_type(vldr_id, old_type_index="v", new_type_index="m")
        for vldr_id in vldr_ids
    ]

    molecular_info = output_data.get("molecular_info", {})
    if qa_test in molecular_info:
        mldr_info = molecular_info[qa_test]
        if "mean" in mldr_info.parameters.values:
            mldr_m_bins = (
                mldr_info
                .sel(parameters="mean")
                .reindex(pair=mldr_ids)
                .assign_coords(pair=vldr_ids)
                .astype("float64")
                .rename("mldr_mean")
            )

            residual = (vldr_m_bins_store - mldr_m_bins) / (
                1.0 - vldr_m_bins_store * mldr_m_bins
            )

            err_p = processing_info["settings_info"]["pcb"].get(
                "pldr_error_threshold",
                0.025,
            )
            _, _, _, sr_lim = pldr_error(
                delta_m=mldr_m_bins,
                delta_v_err=residual,
                delta_p_err_ulim=err_p,
            )

            vldr_info = add_parameter(
                vldr_info,
                name="vldr_residual",
                values=residual,
            )
            vldr_info = add_parameter(
                vldr_info,
                name="sr_limit",
                values=sr_lim,
            )
        else:
            CustomWarning(
                f"VLDR residual skipped: molecular_info[{qa_test!r}] has no mean row."
            )
    else:
        CustomWarning(
            f"VLDR residual skipped: no molecular_info entry for {qa_test}."
        )

    vldr_info = add_parameter(vldr_info, name="ratio_type", values="vldr")

    output_data["pol_cal_ratio"][qa_test] = append_or_replace_pairs(
        output_data["pol_cal_ratio"].get(qa_test),
        vldr,
    )
    output_data["pol_cal_ratio_error"][qa_test] = append_or_replace_pairs(
        output_data["pol_cal_ratio_error"].get(qa_test),
        vldr_error,
    )
    output_data["pol_cal_info"][qa_test] = append_or_replace_info(
        output_data["pol_cal_info"].get(qa_test),
        vldr_info,
    )

    print_entry("VLDR calculation complete!")
    return output_data

def compute_mldr(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute MLDR for ray_pcb from molecular profiles.

    The MLDR channel combinations are created directly from the base ray_pcb
    polarization-calibration metadata, instead of depending on existing eta
    entries.  This makes MLDR independent of compute_calibration_factor while
    still producing pair ids that can be linked to VLDR by using the same
    channel-pair rules and type index "m".

    The molecular input is read from output_data["molecular"] or
    output_data["molec"] and is expected to have dims (channel, bins), after
    optional selection of opto_parameters="atten_bsc".

    The output is saved under:
        output_data["molecular_ratio"]["ray_pcb"]
        output_data["molecular_info"]["ray_pcb"]

    No random error is calculated or stored.  The regional mean is stored in
    molecular_info; no sem row is added.
    """

    output_data = shallow_copy(input_data)
    _ensure_molecular_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    qa_test = "ray_pcb"

    molec_store = output_data.get("molecular", output_data.get("molec", {}))
    pol_cal_info = output_data.get("pol_cal_info", {})

    if qa_test not in molec_store:
        print_entry("MLDR calculation skipped: no ray_pcb molecular profiles found.")
        return output_data
    if qa_test not in vertical_scale:
        print_entry("MLDR calculation skipped: no ray_pcb vertical scale found.")
        return output_data
    if qa_test not in pol_cal_info:
        print_entry("MLDR calculation skipped: no ray_pcb pol_cal_info found.")
        return output_data

    molec = molec_store[qa_test]
    if "opto_parameters" in molec.dims:
        molec = molec.sel(opto_parameters="atten_bsc")
    z = vertical_scale[qa_test]

    if "channel" not in molec.dims:
        raise ValueError(
            f"MLDR molecular input for {qa_test} must have a channel dimension. "
            f"Found dims: {molec.dims}."
        )

    pci_base = _base_pol_cal_pairs(pol_cal_info[qa_test])
    ch_r, ch_t = _read_channels_from_info(pci_base)

    available_channels = set(molec.channel.values.tolist())
    keep = [
        i for i, (r, t) in enumerate(zip(ch_r, ch_t))
        if r in available_channels and t in available_channels
    ]

    if len(keep) == 0:
        print_entry("MLDR calculation skipped: no base ray_pcb channel pairs found in molecular profiles.")
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

    mldr = simple_ratio(
        numerator=molec_r,
        denominator=molec_t,
        info=mldr_info,
    )

    mldr = mldr.rename("ratio")

    z_r = z.sel(channel=ch_r)
    z_pair = channels_to_pairs(z_r, mldr_info)

    mldr_m_bins = mean_in_region(
        da=mldr,
        z=z_pair,
        averaging_range=processing_info["settings_info"]["pcb"]["rayleigh_region"],
    )

    mldr_info = add_parameter(mldr_info, name="mean", values=mldr_m_bins)
    mldr_info = add_parameter(mldr_info, name="ratio_type", values="mldr")

    output_data["molecular_ratio"][qa_test] = append_or_replace_pairs(
        output_data["molecular_ratio"].get(qa_test),
        mldr,
    )
    output_data["molecular_info"][qa_test] = append_or_replace_info(
        output_data["molecular_info"].get(qa_test),
        mldr_info,
    )

    print_entry("MLDR calculation complete!")
    return output_data

