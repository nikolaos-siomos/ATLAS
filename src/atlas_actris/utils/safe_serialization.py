"""Versioned data-only metadata: JSON tags and numeric Zarr sidecars.

No module/class names from a file are imported or instantiated. Unknown types,
tags and codecs fail closed; there is deliberately no pickle compatibility path.
"""
from __future__ import annotations

import base64
import datetime as dt
import json
import math
import ntpath
import os
from pathlib import Path
import re
import tempfile
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

FORMAT = "atlas_metadata_v1"
# Fixed registry, not a list obtained from files or numcodecs.get_codec().
SAFE_CODECS = frozenset({"blosc", "zlib", "gzip", "bz2", "lz4", "zstd",
                        "delta", "astype", "shuffle", "bitround", "fixedscaleoffset",
                        "packbits", "quantize", "crc32", "adler32", "fletcher32",
                        "vlen-utf8", "vlen-bytes"})


class SerializationError(ValueError):
    """Unsupported or malformed persisted data."""


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise SerializationError("Duplicate JSON field")
        result[key] = value
    return result


def read_json(path):
    with open(path, encoding="utf-8") as stream:
        return json.load(stream, object_pairs_hook=_unique_object,
                         parse_constant=lambda value: _invalid("Nonstandard JSON number"))


def write_json(value, path):
    """Commit JSON last and atomically, after any referenced arrays are written."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=str(path.parent),
                                         prefix=".atlas_", suffix=".tmp", delete=False) as stream:
            temporary = stream.name
            json.dump(value, stream, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
        os.replace(temporary, path)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)


def _invalid(message):
    raise SerializationError(message)


def checked_reference(root, reference):
    if (not isinstance(reference, str) or not reference or ntpath.isabs(reference)
            or ntpath.splitdrive(reference)[0] or ":" in reference
            or ".." in reference.replace("\\", "/").split("/")):
        _invalid("Unsafe metadata array reference")
    parent = Path(root).resolve()
    path = parent / reference.replace("\\", "/")
    try:
        resolved = path.resolve()
        resolved.relative_to(parent)
        if resolved == parent:
            _invalid("Metadata reference addresses its containing directory")
    except (ValueError, OSError, RuntimeError) as error:
        raise SerializationError("Metadata array reference escapes its directory") from error
    return path


def validate_zarr_store(path):
    """Validate unconsolidated metadata before opening a local Zarr v2 store.

    Readers use consolidated=False so a conflicting .zmetadata cannot bypass
    these checks. Native decompression/resource limits are a separate concern.
    """
    root = Path(path).resolve()
    def fail(error):
        raise SerializationError("Cannot validate Zarr directory") from error
    for directory, dirs, files in os.walk(root, followlinks=False, onerror=fail):
        for name in dirs + files:
            checked_reference(root, os.path.relpath(os.path.join(directory, name), root))
        if ".zarray" not in files:
            continue
        metadata = read_json(Path(directory) / ".zarray")
        if metadata.get("zarr_format") != 2:
            _invalid("Only Zarr v2 arrays are supported")
        codecs = [metadata.get("compressor")] + (metadata.get("filters") or [])
        for codec in codecs:
            if codec is not None and (not isinstance(codec, dict) or codec.get("id") not in SAFE_CODECS):
                _invalid("Unsupported or unsafe Zarr codec")
        dtype = metadata.get("dtype")
        if isinstance(dtype, list):
            # Structured arrays are not emitted by ATLAS. Fail closed rather
            # than overlooking an object field nested inside a compound dtype.
            _invalid("Structured Zarr dtypes are unsupported")
        if np.dtype(dtype).hasobject:
            filters = metadata.get("filters") or []
            if not filters or filters[0].get("id") not in {"vlen-utf8", "vlen-bytes"}:
                _invalid("Object Zarr arrays require a supported string codec")


def _dtype(name):
    """Accept known pandas/numpy dtype syntax, never arbitrary extension imports."""
    if not isinstance(name, str):
        _invalid("Invalid dtype")
    if name in {"Int8", "Int16", "Int32", "Int64", "UInt8", "UInt16", "UInt32", "UInt64",
                "Float32", "Float64", "boolean", "string"}:
        return name
    if re.fullmatch(r"datetime64\[(?:s|ms|us|ns), [A-Za-z0-9_/+:. -]+\]", name):
        return name
    try:
        result = np.dtype(name)
    except (TypeError, ValueError) as error:
        raise SerializationError("Unsupported dtype") from error
    if result.fields or result.subdtype or result.kind not in "biufcMmOSU":
        _invalid("Unsupported dtype")
    return result


class _Encoder:
    def __init__(self, path=None, validate_only=False):
        self.path = None if path is None else Path(path)
        self.validate_only = validate_only
        self.counter = 0
        self.active = set()

    def encode(self, value):
        # Only exact primitive types bypass tagging (numpy subclasses do not).
        if value is None or type(value) in (bool, int, str):
            return value
        if type(value) is float:
            return value if math.isfinite(value) else {"t": "float", "v": repr(value)}
        identity = id(value)
        if identity in self.active:
            _invalid("Cyclic metadata is unsupported")
        self.active.add(identity)
        try:
            return self._encode(value)
        finally:
            self.active.remove(identity)

    def _encode(self, value):
        e = self.encode
        if value is pd.NA:
            return {"t": "pd_na"}
        if value is pd.NaT:
            return {"t": "pd_nat"}
        if isinstance(value, pd.Timestamp):
            return {"t": "timestamp", "v": value.isoformat(), "unit": value.unit,
                    "tz": getattr(value.tzinfo, "key", getattr(value.tzinfo, "zone", None))}
        if isinstance(value, pd.Timedelta):
            return {"t": "pd_timedelta", "v": str(value.value), "unit": value.unit}
        if isinstance(value, np.generic):
            if value.dtype.hasobject or value.dtype.fields or value.dtype.kind not in "biufcMmSU":
                _invalid("Unsupported numpy scalar dtype")
            return {"t": "np_scalar", "dtype": value.dtype.str,
                    "v": base64.b64encode(value.tobytes()).decode("ascii")}
        if isinstance(value, complex):
            return {"t": "complex", "v": [e(value.real), e(value.imag)]}
        if isinstance(value, (bytes, bytearray)):
            return {"t": "bytes", "v": base64.b64encode(value).decode("ascii"), "mutable": isinstance(value, bytearray)}
        if isinstance(value, dt.datetime):
            return {"t": "datetime", "v": value.isoformat(), "fold": value.fold,
                    "tz": getattr(value.tzinfo, "key", getattr(value.tzinfo, "zone", None))}
        if isinstance(value, dt.date):
            return {"t": "date", "v": value.isoformat()}
        if isinstance(value, dt.time):
            return {"t": "time", "v": value.isoformat(), "fold": value.fold}
        if isinstance(value, dt.timedelta):
            return {"t": "timedelta", "v": [value.days, value.seconds, value.microseconds]}
        if isinstance(value, Path):
            return {"t": "path", "v": str(value)}
        if isinstance(value, dict):
            return {"t": "dict", "v": [[e(k), e(v)] for k, v in value.items()]}
        if isinstance(value, (tuple, list, set, frozenset)):
            return {"t": type(value).__name__, "v": [e(v) for v in value]}
        if isinstance(value, pd.Categorical):
            return {"t": "categorical", "codes": value.codes.tolist(), "categories": e(value.categories), "ordered": value.ordered}
        if isinstance(value, pd.MultiIndex):
            return {"t": "multiindex", "levels": [e(v) for v in value.levels],
                    "codes": [v.tolist() for v in value.codes], "names": e(list(value.names))}
        if isinstance(value, pd.RangeIndex):
            return {"t": "rangeindex", "start": value.start, "stop": value.stop, "step": value.step, "name": e(value.name)}
        if isinstance(value, pd.Index):
            if isinstance(value, pd.CategoricalIndex):
                return {"t": "categoryindex", "data": e(value.array), "name": e(value.name)}
            _dtype(str(value.dtype))
            return {"t": "index", "v": e(value.tolist()), "dtype": str(value.dtype), "name": e(value.name),
                    "freq": getattr(value, "freqstr", None)}
        if isinstance(value, pd.Series):
            categorical = isinstance(value.dtype, pd.CategoricalDtype)
            if not categorical:
                _dtype(str(value.dtype))
            return {"t": "series", "index": e(value.index), "name": e(value.name),
                    "dtype": "category" if categorical else str(value.dtype),
                    "v": e(value.array if categorical else value.tolist()), "attrs": e(value.attrs)}
        if isinstance(value, pd.DataFrame):
            return {"t": "dataframe", "index": e(value.index), "columns": e(value.columns),
                    "data": [e(value.iloc[:, i]) for i in range(value.shape[1])], "attrs": e(value.attrs)}
        if isinstance(value, np.ma.MaskedArray):
            return {"t": "masked_array", "data": e(np.asarray(value.data)),
                    "mask": e(np.ma.getmaskarray(value)), "fill": e(value.fill_value)}
        if isinstance(value, np.ndarray):
            if value.dtype.fields:
                _invalid("Structured numpy arrays are unsupported")
            if value.dtype.hasobject:
                return {"t": "object_array", "shape": list(value.shape), "v": [e(v) for v in value.flat]}
            if value.dtype.kind in "SU":
                return {"t": "text_array", "shape": list(value.shape), "dtype": value.dtype.str,
                        "v": base64.b64encode(value.tobytes()).decode("ascii")}
            return self.array(value)
        # Imports are fixed in code, never supplied by the persisted document.
        try:
            import xarray as xr
        except ImportError:
            xr = None
        if xr is not None and isinstance(value, (xr.DataArray, xr.Dataset)):
            variables = value.data_vars if isinstance(value, xr.Dataset) else {"value": value.variable}
            return {"t": "dataset" if isinstance(value, xr.Dataset) else "dataarray",
                    "name": e(value.name) if isinstance(value, xr.DataArray) else None,
                    "variables": [[e(k), self.variable(v)] for k, v in variables.items()],
                    "coords": [[e(k), self.variable(v)] for k, v in value.coords.items()], "attrs": e(value.attrs)}
        try:
            import dask.array as da
        except ImportError:
            da = None
        if da is not None and isinstance(value, da.Array):
            if value.dtype.hasobject:
                return e(value.compute())
            return self.array(value, lazy=True)
        _invalid(f"Unsupported metadata type: {type(value).__module__}.{type(value).__name__}")

    def variable(self, value):
        return {"dims": list(value.dims), "data": self.encode(value.data), "attrs": self.encode(value.attrs)}

    def array(self, value, lazy=False):
        if value.dtype.kind not in "biufcMm":
            _invalid("Unsupported numeric array dtype")
        if self.validate_only:
            return {"t": "array", "path": "validation-only", "lazy": lazy}
        if self.path is None:
            _invalid("Numeric arrays require a metadata file and Zarr storage")
        reference = self.path.name + ".arrays/" + str(self.counter) + ".zarr"
        self.counter += 1
        target = checked_reference(self.path.parent, reference)
        if target.exists():
            _invalid("Metadata array target already exists; write to a fresh entry")
        target.parent.mkdir(parents=True, exist_ok=True)
        if lazy:
            import dask.array as da
            value = value.rechunk(tuple(max(parts) for parts in value.chunks))
            da.to_zarr(value, str(target))
        else:
            import zarr
            zarr.save_array(str(target), value, compressor=None)
        return {"t": "array", "path": reference, "lazy": lazy}


def _decode(value, root):
    if value is None or type(value) in (str, bool, int, float):
        return value
    if not isinstance(value, dict) or not isinstance(value.get("t"), str):
        _invalid("Invalid tagged metadata value")
    d = lambda v: _decode(v, root)
    tag = value["t"]
    v = value.get("v")
    if tag == "float":
        if v not in ("nan", "inf", "-inf"):
            _invalid("Invalid special float")
        return float(v)
    if tag in {"list", "tuple", "set", "frozenset"}:
        return {"list": list, "tuple": tuple, "set": set, "frozenset": frozenset}[tag](d(x) for x in v)
    if tag == "dict":
        result = {}
        for key, item in v:
            key = d(key)
            if key in result:
                _invalid("Duplicate metadata dictionary key")
            result[key] = d(item)
        return result
    if tag == "pd_na": return pd.NA
    if tag == "pd_nat": return pd.NaT
    if tag == "timestamp":
        result = pd.Timestamp(v)
        if value.get("tz"): result = result.tz_convert(value["tz"])
        return result.as_unit(value.get("unit", "ns"))
    if tag == "pd_timedelta": return pd.Timedelta(int(v), unit="ns").as_unit(value.get("unit", "ns"))
    if tag == "datetime":
        result = dt.datetime.fromisoformat(v)
        if value.get("tz"): result = result.astimezone(ZoneInfo(value["tz"]))
        return result.replace(fold=value.get("fold", 0))
    if tag == "date": return dt.date.fromisoformat(v)
    if tag == "time": return dt.time.fromisoformat(v).replace(fold=value.get("fold", 0))
    if tag == "timedelta": return dt.timedelta(days=v[0], seconds=v[1], microseconds=v[2])
    if tag == "path": return Path(v)
    if tag == "complex": return complex(d(v[0]), d(v[1]))
    if tag == "bytes":
        raw = base64.b64decode(v, validate=True)
        return bytearray(raw) if value.get("mutable") else raw
    if tag in {"np_scalar", "text_array"}:
        dtype = np.dtype(_dtype(value["dtype"]))
        if dtype.hasobject or (tag == "text_array" and dtype.kind not in "SU"):
            _invalid("Invalid raw-buffer dtype")
        result = np.frombuffer(base64.b64decode(v, validate=True), dtype=dtype).copy()
        return result.reshape(value["shape"]) if tag == "text_array" else result.reshape(())[()]
    if tag == "object_array":
        shape = value["shape"]
        if any(type(n) is not int or n < 0 for n in shape) or math.prod(shape) != len(v):
            _invalid("Invalid object array shape")
        result = np.empty(shape, dtype=object)
        for i, item in enumerate(v): result.flat[i] = d(item)
        return result
    if tag == "masked_array":
        return np.ma.array(d(value["data"]), mask=d(value["mask"]), fill_value=d(value["fill"]))
    if tag == "array":
        path = checked_reference(root, value["path"])
        validate_zarr_store(path)
        import zarr
        array = zarr.open_array(str(path), mode="r")
        if array.dtype.hasobject or array.dtype.kind not in "biufcMm":
            _invalid("Numeric sidecar has an invalid dtype")
        if value.get("lazy"):
            import dask.array as da
            return da.from_zarr(str(path))
        return np.asarray(array)
    if tag == "categorical":
        return pd.Categorical.from_codes(value["codes"], categories=d(value["categories"]), ordered=value["ordered"])
    if tag == "rangeindex":
        return pd.RangeIndex(value["start"], value["stop"], value["step"], name=d(value["name"]))
    if tag == "multiindex":
        return pd.MultiIndex(levels=[d(x) for x in value["levels"]], codes=value["codes"], names=d(value["names"]))
    if tag == "categoryindex": return pd.CategoricalIndex(d(value["data"]), name=d(value["name"]))
    if tag == "index":
        result = pd.Index(d(v), dtype=_dtype(value["dtype"]), name=d(value["name"]))
        if value.get("freq") is not None: result.freq = value["freq"]
        return result
    if tag == "series":
        dtype = "category" if value["dtype"] == "category" else _dtype(value["dtype"])
        result = pd.Series(d(v), index=d(value["index"]), name=d(value["name"]), dtype=dtype)
        result.attrs = d(value.get("attrs", {"t": "dict", "v": []}))
        return result
    if tag == "dataframe":
        columns = [d(x).reset_index(drop=True) for x in value["data"]]
        result = pd.concat(columns, axis=1) if columns else pd.DataFrame(index=range(len(d(value["index"]))))
        result.index, result.columns = d(value["index"]), d(value["columns"])
        result.attrs = d(value.get("attrs", {"t": "dict", "v": []}))
        return result
    if tag in {"dataset", "dataarray"}:
        import xarray as xr
        def variable(spec):
            return xr.Variable(spec["dims"], d(spec["data"]), attrs=d(spec["attrs"]))
        coords = {d(k): variable(spec) for k, spec in value["coords"]}
        variables = {d(k): variable(spec) for k, spec in value["variables"]}
        if tag == "dataset": return xr.Dataset(variables, coords=coords, attrs=d(value["attrs"]))
        return xr.DataArray(variables["value"], coords=coords, name=d(value["name"]), attrs=d(value["attrs"]))
    _invalid(f"Unsupported metadata tag: {tag!r}")


def validate_metadata(value):
    """Validate supported types before overwriting any previous export."""
    _Encoder(validate_only=True).encode(value)


def write_metadata(value, path):
    encoded = _Encoder(path).encode(value)
    write_json({"format": FORMAT, "data": encoded}, path)


def read_metadata(path):
    try:
        envelope = read_json(path)
        if not isinstance(envelope, dict) or envelope.get("format") != FORMAT:
            _invalid("Unsupported metadata format; regenerate the cache or export")
        return _decode(envelope["data"], Path(path).resolve().parent)
    except (KeyError, TypeError, OverflowError, RecursionError) as error:
        raise SerializationError("Malformed metadata document") from error
