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

    allowed_indices = ["g", "e", "c", "v", "d", "s", "f"]

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


def append_or_replace_info(
    existing: Optional[xr.DataArray],
    new: xr.DataArray,
) -> xr.DataArray:
    """Add pair entries to pol_cal_info, replacing duplicate pair ids."""

    if existing is None:
        return new

    if "pair" not in existing.dims or "pair" not in new.dims:
        raise ValueError("Both existing and new pol_cal_info objects need a pair dimension.")

    drop_ids = [pid for pid in new.pair.values if pid in existing.pair.values]
    if drop_ids:
        existing = existing.drop_sel(pair=drop_ids)

    if existing.sizes.get("pair", 0) == 0:
        combined = new
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


def _find_eta_for_pairs(
    output_data: Dict[str, Dict[str, Any]],
    ch_r: Sequence[str],
    ch_t: Sequence[str],
    calibrated_pair_ids: Sequence[str],
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Find eta and eta_error for the requested channel pairs and align to calibrated_pair_ids."""

    eta_pair_ids = [
        make_ratio_id(r, t, type_index="e")
        for r, t in zip(ch_r, ch_t)
    ]

    eta_sources = ["pcb", "pcb_aux", "pcb_p45", "pcb_m45", "pcb_aux_p45", "pcb_aux_m45"]
    eta_values = []
    eta_errors = []
    found_ids = []

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

            if "pair" in candidate.dims and eta_id in candidate.pair.values:
                eta_da = candidate.sel(pair=eta_id)
                eta_err_da = candidate_err.sel(pair=eta_id)
                break

        if eta_da is None or eta_err_da is None:
            continue

        eta_values.append(eta_da.expand_dims(pair=[calibrated_id]))
        eta_errors.append(eta_err_da.expand_dims(pair=[calibrated_id]))
        found_ids.append(calibrated_id)

    if not eta_values:
        raise KeyError("No matching eta entries found for ray_pcb channel pairs.")

    eta = xr.concat(eta_values, dim="pair")
    eta_error = xr.concat(eta_errors, dim="pair")

    preferred = [dim for dim in ["time", "pair", "bins"] if dim in eta.dims]
    other = [dim for dim in eta.dims if dim not in preferred]
    eta = eta.transpose(*preferred, *other)
    eta_error = eta_error.transpose(*preferred, *other)

    return eta, eta_error


def compute_gain_ratio(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute gain ratios for pcb/pcb_aux +/-45 QA tests.

    The saved products use:
        ratio_type = "gain_ratio"
        pair id type index = "g"
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

        ratio_da = simple_ratio(
            numerator=sig_r,
            denominator=sig_t,
            info=pci,
        )

        sig_r_err = sig_err.sel(channel=ch_r)
        sig_t_err = sig_err.sel(channel=ch_t)

        ratio_error_da = ratio_error_independent(
            numerator=sig_r,
            denominator=sig_t,
            numerator_error=sig_r_err,
            denominator_error=sig_t_err,
            ratio_values=ratio_da,
            info=pci,
        )

        z_r = z.sel(channel=ch_r)
        z_pair = channels_to_pairs(z_r, pci)

        ratio_m_time, ratio_error_m_time = _time_mean_and_error(
            ratio_da,
            ratio_error_da,
        )

        ratio_m_bins = mean_in_region(
            da=ratio_m_time,
            z=z_pair,
            averaging_range=averaging_range,
        )

        ratio_error_m_bins = sem_in_region(
            da_error=ratio_error_m_time,
            z=z_pair,
            averaging_range=averaging_range,
        )

        pci = add_parameter(pci=pci, name="mean", values=ratio_m_bins)
        pci = add_parameter(pci=pci, name="sem", values=ratio_error_m_bins)
        pci = add_parameter(pci=pci, name="ratio_type", values="gain_ratio")

        output_data["pol_cal_ratio"][qa_test] = append_or_replace_pairs(
            output_data["pol_cal_ratio"].get(qa_test),
            ratio_da,
        )
        output_data["pol_cal_ratio_error"][qa_test] = append_or_replace_pairs(
            output_data["pol_cal_ratio_error"].get(qa_test),
            ratio_error_da,
        )
        output_data["pol_cal_info"][qa_test] = append_or_replace_info(
            output_data["pol_cal_info"].get(qa_test),
            pci,
        )

    print_entry("Gain ratio calculation complete!")
    return output_data


def compute_calibration_factor(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute calibration products from +/-45 gain ratios.

    The p45 and m45 gain ratios are averaged over time first, because their
    time coordinates may differ.  For each pcb family, the function stores
    three ratio categories under the combined key and under the two source keys:

        gain_ratio : mean gain-ratio products, using g IDs
        eta_s_f    : gain-ratio products corrected by R_to_T, using f IDs
        eta        : final eta products corrected by K, using e IDs

    Output key groups are:
        pcb      from pcb_p45 / pcb_m45
        pcb_aux  from pcb_aux_p45 / pcb_aux_m45

    For each group, products are saved under:
        target_key, p45, m45
    for pol_cal_ratio, pol_cal_ratio_error, and pol_cal_info.
    """

    assign_pcb = {
        "pcb": ("pcb_p45", "pcb_m45"),
        "pcb_aux": ("pcb_aux_p45", "pcb_aux_m45"),
    }

    output_data = shallow_copy(input_data)
    _ensure_pol_cal_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    pol_cal_ratio = output_data["pol_cal_ratio"]
    pol_cal_ratio_error = output_data["pol_cal_ratio_error"]
    pol_cal_info = output_data["pol_cal_info"]

    averaging_range = processing_info["settings_info"]["pcb"]["calibration_region"]

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

    for target_key, (p45, m45) in assign_pcb.items():
        if p45 not in pol_cal_ratio or m45 not in pol_cal_ratio:
            continue
        if p45 not in pol_cal_ratio_error or m45 not in pol_cal_ratio_error:
            continue
        if p45 not in pol_cal_info or m45 not in pol_cal_info:
            continue
        if p45 not in vertical_scale or m45 not in vertical_scale:
            continue

        pci_p45 = pol_cal_info[p45]
        pci_m45 = pol_cal_info[m45]

        p45_gain_ids = _select_pairs_by_type(pci_p45, "gain_ratio")
        m45_gain_ids = _select_pairs_by_type(pci_m45, "gain_ratio")
        gain_ids = [pid for pid in p45_gain_ids if pid in m45_gain_ids]

        if len(gain_ids) == 0:
            continue

        ratio_p45 = pol_cal_ratio[p45].sel(pair=gain_ids)
        ratio_p45_err = pol_cal_ratio_error[p45].sel(pair=gain_ids)
        ratio_m45 = pol_cal_ratio[m45].sel(pair=gain_ids)
        ratio_m45_err = pol_cal_ratio_error[m45].sel(pair=gain_ids)

        eta_s_p45, eta_s_p45_err = _time_mean_and_error(
            ratio_p45,
            ratio_p45_err,
        )
        eta_s_m45, eta_s_m45_err = _time_mean_and_error(
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
                f"Skipping calibration factor for {target_key}: "
                f"vertical scales of {p45} and {m45} do not match."
            )
            continue

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

        eta_s_product = eta_s_p45 * eta_s_m45
        eta_s = np.sqrt(eta_s_product.where(eta_s_product >= 0))
        eta_s_err = _product_sqrt_error(
            eta_s_p45,
            eta_s_p45_err,
            eta_s_m45,
            eta_s_m45_err,
            eta_s,
        )

        # Gain-ratio family, stored with g IDs.
        gain_products = {
            target_key: (eta_s, eta_s_err),
            p45: (eta_s_p45, eta_s_p45_err),
            m45: (eta_s_m45, eta_s_m45_err),
        }

        # eta_s_f family, stored with f IDs.
        eta_s_f = eta_s / trans_ratio
        eta_s_f_error = eta_s_err / np.abs(trans_ratio)

        eta_s_f_p45 = eta_s_p45 / trans_ratio
        eta_s_f_p45_error = eta_s_p45_err / np.abs(trans_ratio)

        eta_s_f_m45 = eta_s_m45 / trans_ratio
        eta_s_f_m45_error = eta_s_m45_err / np.abs(trans_ratio)

        eta_s_f_products = {
            target_key: (eta_s_f, eta_s_f_error),
            p45: (eta_s_f_p45, eta_s_f_p45_error),
            m45: (eta_s_f_m45, eta_s_f_m45_error),
        }

        # eta family, stored with e IDs.
        eta = eta_s / K
        eta_error = eta_s_err / np.abs(K)

        eta_p45 = eta_s_p45 / K
        eta_p45_error = eta_s_p45_err / np.abs(K)

        eta_m45 = eta_s_m45 / K
        eta_m45_error = eta_s_m45_err / np.abs(K)

        eta_products = {
            target_key: (eta, eta_error),
            p45: (eta_p45, eta_p45_error),
            m45: (eta_m45, eta_m45_error),
        }

        product_specs = [
            {
                "ratio_type": "gain_ratio",
                "pair_ids": gain_ids,
                "products": gain_products,
                "z_pair": z_pair_p45,
            },
            {
                "ratio_type": "eta_s_f",
                "pair_ids": eta_s_f_ids,
                "products": eta_s_f_products,
                "z_pair": z_pair_p45.assign_coords(pair=eta_s_f_ids),
            },
            {
                "ratio_type": "eta",
                "pair_ids": eta_ids,
                "products": eta_products,
                "z_pair": z_pair_p45.assign_coords(pair=eta_ids),
            },
        ]

        for spec in product_specs:
            ratio_type = spec["ratio_type"]
            pair_ids = spec["pair_ids"]
            z_pair = spec["z_pair"]

            info = _pair_info_for_channels(
                pci=pci_gain,
                ch_r=ch_r,
                ch_t=ch_t,
                pair_ids=pair_ids,
                ratio_type=ratio_type,
            )

            for key, (values, errors) in spec["products"].items():
                values = values.assign_coords(pair=pair_ids).rename("ratio")
                errors = errors.assign_coords(pair=pair_ids).rename("ratio_error")

                info_with_stats = _add_region_stats(
                    info=info.copy(),
                    values=values,
                    errors=errors,
                    z_pair=z_pair,
                    ratio_type=ratio_type,
                )

                _store_product(
                    key=key,
                    values=values,
                    errors=errors,
                    info=info_with_stats,
                )

    print_entry("Calibration factor calculation complete!")

    return output_data


def compute_calibrated_ratio(
    processing_info: Dict[str, Any],
    input_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute calibrated ratios for ray_pcb.

    The raw reflected/transmitted ratios are divided by eta.  The saved products use:
        ratio_type = "calibrated_ratio"
        pair id type index = "c"
    """

    output_data = shallow_copy(input_data)
    _ensure_pol_cal_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    profiles = output_data["profile"]
    profile_errors = output_data["profile_error"]
    channel_info = output_data["channel_info"]
    pol_cal_info = output_data["pol_cal_info"]

    allowed_qa_tests = ["ray_pcb"]
    averaging_range = processing_info["settings_info"]["pcb"]["rayleigh_region"]

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
            make_ratio_id(numerator, denominator, type_index="c")
            for numerator, denominator in zip(ch_r, ch_t)
        ]

        if len(pair_ids) == 0:
            continue

        pci = _pair_info_for_channels(
            pci=pci_base,
            ch_r=ch_r,
            ch_t=ch_t,
            pair_ids=pair_ids,
            ratio_type="calibrated_ratio",
        )

        sig_r = sig.sel(channel=ch_r)
        sig_t = sig.sel(channel=ch_t)

        uncalibrated = simple_ratio(
            numerator=sig_r,
            denominator=sig_t,
            info=pci,
        )

        sig_r_err = sig_err.sel(channel=ch_r)
        sig_t_err = sig_err.sel(channel=ch_t)

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

        # Drop pairs for which eta was not found.
        pair_ids_found = list(eta.pair.values)
        uncalibrated = uncalibrated.sel(pair=pair_ids_found)
        uncalibrated_error = uncalibrated_error.sel(pair=pair_ids_found)
        pci = pci.sel(pair=pair_ids_found)
        ch_r_found = pci.sel(parameters="ch_r").values.tolist()

        calibrated = uncalibrated / eta
        calibrated_error = _division_error(
            numerator=uncalibrated,
            numerator_error=uncalibrated_error,
            denominator=eta,
            denominator_error=eta_error,
            quotient=calibrated,
        )

        calibrated = calibrated.rename("ratio")
        calibrated_error = calibrated_error.rename("ratio_error")

        z_r = z.sel(channel=ch_r_found)
        z_pair = channels_to_pairs(z_r, pci)

        calibrated_m_time, calibrated_error_m_time = _time_mean_and_error(
            calibrated,
            calibrated_error,
        )

        calibrated_m_bins = mean_in_region(
            da=calibrated_m_time,
            z=z_pair,
            averaging_range=averaging_range,
        )
        calibrated_error_m_bins = sem_in_region(
            da_error=calibrated_error_m_time,
            z=z_pair,
            averaging_range=averaging_range,
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
    Compute VLDR for ray_pcb from calibrated ratios.

    Cross-talk parameters g_r, g_t, h_r, h_t are read from channel_info for the
    associated channels.  The correction is intentionally left as an identity
    placeholder until the final formula is decided.

    The saved products use:
        ratio_type = "vldr"
        pair id type index = "v"
    """

    output_data = shallow_copy(input_data)
    _ensure_pol_cal_output_dicts(output_data)

    vertical_scale_name = processing_info["caller_info"]["vertical_scale"]
    vertical_scale = output_data[vertical_scale_name]

    qa_test = "ray_pcb"
    averaging_range = processing_info["settings_info"]["pcb"]["rayleigh_region"]

    if qa_test not in output_data.get("pol_cal_ratio", {}):
        print_entry("VLDR calculation skipped: no ray_pcb calibrated ratios found.")
        return output_data
    if qa_test not in output_data.get("pol_cal_ratio_error", {}):
        print_entry("VLDR calculation skipped: no ray_pcb calibrated ratio errors found.")
        return output_data
    if qa_test not in output_data.get("pol_cal_info", {}):
        print_entry("VLDR calculation skipped: no ray_pcb pol_cal_info found.")
        return output_data
    if qa_test not in vertical_scale:
        print_entry("VLDR calculation skipped: no ray_pcb vertical scale found.")
        return output_data
    if qa_test not in output_data.get("channel_info", {}):
        print_entry("VLDR calculation skipped: no ray_pcb channel_info found.")
        return output_data

    ratio_store = output_data["pol_cal_ratio"][qa_test]
    error_store = output_data["pol_cal_ratio_error"][qa_test]
    pci_store = output_data["pol_cal_info"][qa_test]
    channel_info = output_data["channel_info"][qa_test]

    calibrated_ids = _select_pairs_by_type(pci_store, "calibrated_ratio")
    if len(calibrated_ids) == 0:
        print_entry("VLDR calculation skipped: no calibrated_ratio entries found.")
        return output_data

    calibrated = ratio_store.sel(pair=calibrated_ids)
    calibrated_error = error_store.sel(pair=calibrated_ids)
    pci_cal = pci_store.sel(pair=calibrated_ids)

    ch_r, ch_t = _read_channels_from_info(pci_cal)
    vldr_ids = [
        make_ratio_id(r, t, type_index="v")
        for r, t in zip(ch_r, ch_t)
    ]

    vldr_info = _pair_info_for_channels(
        pci=pci_cal,
        ch_r=ch_r,
        ch_t=ch_t,
        pair_ids=vldr_ids,
        ratio_type="vldr",
    )

    g_r = xr.DataArray(
        [_channel_parameter(channel_info, ch, "g_r", default=0.0) for ch in ch_r],
        dims=["pair"],
        coords={"pair": calibrated_ids},
    )
    g_t = xr.DataArray(
        [_channel_parameter(channel_info, ch, "g_t", default=0.0) for ch in ch_t],
        dims=["pair"],
        coords={"pair": calibrated_ids},
    )
    h_r = xr.DataArray(
        [_channel_parameter(channel_info, ch, "h_r", default=0.0) for ch in ch_r],
        dims=["pair"],
        coords={"pair": calibrated_ids},
    )
    h_t = xr.DataArray(
        [_channel_parameter(channel_info, ch, "h_t", default=0.0) for ch in ch_t],
        dims=["pair"],
        coords={"pair": calibrated_ids},
    )

    # Placeholder for the final cross-talk correction.  Keep the parameters in
    # scope so the implementation below can be replaced directly, e.g. using
    # calibrated, g_r, g_t, h_r, h_t.
    _ = (g_r, g_t, h_r, h_t)
    vldr = calibrated.copy()
    vldr_error = calibrated_error.copy()

    vldr = vldr.assign_coords(pair=vldr_ids).rename("ratio")
    vldr_error = vldr_error.assign_coords(pair=vldr_ids).rename("ratio_error")

    z = vertical_scale[qa_test]
    z_r = z.sel(channel=ch_r)
    z_pair = channels_to_pairs(z_r, vldr_info)

    vldr_m_time, vldr_error_m_time = _time_mean_and_error(vldr, vldr_error)

    vldr_m_bins = mean_in_region(
        da=vldr_m_time,
        z=z_pair,
        averaging_range=averaging_range,
    )
    vldr_error_m_bins = sem_in_region(
        da_error=vldr_error_m_time,
        z=z_pair,
        averaging_range=averaging_range,
    )

    vldr_info = add_parameter(vldr_info, name="mean", values=vldr_m_bins)
    vldr_info = add_parameter(vldr_info, name="sem", values=vldr_error_m_bins)
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


# Backwards-compatible aliases for the previous half-made names.
compute_uncalibrated_ratio = compute_calibrated_ratio
compute_calibrated_ratio_pol_cal = compute_calibrated_ratio
compute_eta = compute_calibration_factor
