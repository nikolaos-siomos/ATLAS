"""Data-only serialization, legacy rejection and lazy-array regressions."""
import datetime as dt
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from atlas_actris.utils import safe_serialization as codec

HAS_ARRAY_STACK = all(importlib.util.find_spec(name) is not None for name in ("xarray", "dask", "zarr"))


class MetadataTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas serialization ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.path = self.folder / "metadata.json"

    def roundtrip(self, value):
        codec.write_metadata(value, self.path)
        return codec.read_metadata(self.path)

    def test_scalar_containers_and_nonstring_dictionary_keys(self):
        value = {("tuple", 1): [True, None, b"abc", bytearray(b"def"),
                                complex(1, 2), {1, 2}, frozenset({3}), Path("relative")],
                 "looks_like_a_tag": {"t": "execute", "v": "inert text"}}
        self.assertEqual(self.roundtrip(value), value)
        self.assertEqual(json.loads(self.path.read_text())["format"], codec.FORMAT)

    def test_numpy_scalars_and_special_floats(self):
        values = [np.int16(4), np.float32(2.5), np.complex64(1+2j), np.str_("abc"),
                  np.datetime64("2025-01-02", "D"), np.timedelta64(2, "h"), np.datetime64("NaT", "ns")]
        restored = self.roundtrip(values + [float("nan"), float("inf"), float("-inf")])
        for expected, actual in zip(values, restored):
            self.assertEqual(actual.dtype, expected.dtype)
            self.assertEqual(actual.tobytes(), expected.tobytes())
        self.assertTrue(np.isnan(restored[-3]))
        self.assertEqual(restored[-2:], [float("inf"), float("-inf")])

    def test_dates_times_timezone_and_missing_values(self):
        values = [dt.datetime(2025, 10, 26, 2, 30, tzinfo=ZoneInfo("Europe/Berlin"), fold=1),
                  dt.date(2025, 1, 1), dt.time(1, 2, 3), dt.timedelta(days=-1, microseconds=2),
                  pd.Timestamp("2025-01-01T01:02:03.123456789", tz="Europe/Berlin"),
                  pd.Timedelta(123, unit="ns"), pd.NA, pd.NaT]
        result = self.roundtrip(values)
        self.assertEqual(result[:-2], values[:-2])
        self.assertEqual(result[0].fold, 1)
        self.assertIs(result[-2], pd.NA)
        self.assertIs(result[-1], pd.NaT)

    def test_mixed_object_arrays_preserve_cell_types(self):
        value = np.empty((2, 3), dtype=object)
        value.flat[:] = [None, pd.NA, pd.NaT, np.float32(2), (1, 2), "text"]
        result = self.roundtrip(value)
        self.assertEqual(result.shape, value.shape)
        self.assertEqual(result.dtype, value.dtype)
        self.assertIsNone(result[0, 0]); self.assertIs(result[0, 1], pd.NA)
        self.assertIs(result[0, 2], pd.NaT)
        self.assertIsInstance(result[1, 0], np.float32)
        self.assertEqual(result[1, 1], (1, 2))

    def test_fixed_width_text_arrays(self):
        for value in (np.array(["é", "abc"], dtype="U5"), np.array([b"a", b"b"], dtype="S4")):
            result = self.roundtrip(value)
            np.testing.assert_array_equal(result, value)
            self.assertEqual(result.dtype, value.dtype)

    def test_pandas_table_dtypes_and_missing_values(self):
        table = pd.DataFrame({"nullable": pd.Series([1, pd.NA], dtype="Int64"),
                              "float": [1.2, np.nan], "flag": pd.Series([True, pd.NA], dtype="boolean"),
                              "text": pd.Series(["a", pd.NA], dtype="string"),
                              "time": pd.date_range("2025-01-01", periods=2, tz="UTC"),
                              "category": pd.Categorical(["b", "a"], categories=["a", "b", "unused"], ordered=True)})
        table.index = pd.Index(["channel1", "channel2"], name="channel")
        table.attrs = {"source": "test"}
        actual = self.roundtrip(table)
        pd.testing.assert_frame_equal(actual, table)
        self.assertEqual(actual.attrs, table.attrs)

    def test_datetime_multiindex_and_duplicate_columns(self):
        for index in (pd.date_range("2025-01-01", periods=2, freq="h", tz="Europe/Berlin", name="time"),
                      pd.MultiIndex(levels=[["a", "b", "unused"], [1, 2]], codes=[[0, 1], [0, 1]], names=["site", "id"])):
            table = pd.DataFrame([[1, "a"], [2, "b"]], index=index, columns=["same", "same"])
            pd.testing.assert_frame_equal(self.roundtrip(table), table)

    def test_empty_tables_and_series(self):
        for table in (pd.DataFrame(), pd.DataFrame(index=pd.RangeIndex(2, 5, name="row")),
                      pd.DataFrame({"a": pd.Series([], dtype="float32")})):
            pd.testing.assert_frame_equal(self.roundtrip(table), table)
        value = pd.Series([1., 2.], index=["a", "b"], name="system_info")
        pd.testing.assert_series_equal(self.roundtrip(value), value)

    def test_unknown_objects_rejected_without_invoking_reduce(self):
        class Unsupported:
            def __reduce__(self):
                raise AssertionError("Must never invoke pickle protocol")
        with self.assertRaises(codec.SerializationError):
            self.roundtrip(Unsupported())
        self.assertFalse(self.path.exists())

    def test_cyclic_objects_rejected(self):
        value = []; value.append(value)
        with self.assertRaises(codec.SerializationError):
            self.roundtrip(value)

    def test_unknown_tags_and_object_buffers_rejected(self):
        for data in ({"t": "python_object", "module": "os", "class": "system"},
                     {"t": "np_scalar", "dtype": "O", "v": "AAAA"}):
            codec.write_json({"format": codec.FORMAT, "data": data}, self.path)
            with self.assertRaises(codec.SerializationError):
                codec.read_metadata(self.path)

    def test_duplicate_json_fields_and_nonstandard_numbers_rejected(self):
        for text in ('{"a":1,"a":2}', '{"a":NaN}'):
            self.path.write_text(text)
            with self.assertRaises(codec.SerializationError):
                codec.read_json(self.path)

    def test_escaping_array_references_rejected_without_opening_store(self):
        for reference in ("../outside.zarr", "/tmp/outside.zarr", "C:\\outside.zarr"):
            codec.write_json({"format": codec.FORMAT, "data": {"t": "array", "path": reference}}, self.path)
            with self.assertRaises(codec.SerializationError):
                codec.read_metadata(self.path)

    def test_zarr_pickle_codec_rejected_before_decoding(self):
        store = self.folder / "array.zarr"; store.mkdir()
        for field in ("compressor", "filters"):
            metadata = {"zarr_format": 2, "dtype": "|O", "compressor": None, "filters": None}
            metadata[field] = {"id": "pickle"} if field == "compressor" else [{"id": "pickle"}]
            codec.write_json(metadata, store / ".zarray")
            with self.assertRaises(codec.SerializationError):
                codec.validate_zarr_store(store)


@unittest.skipUnless(HAS_ARRAY_STACK, "Requires installed xarray, dask and zarr")
class ArrayAndStageTests(unittest.TestCase):
    setUp = MetadataTests.setUp
    roundtrip = MetadataTests.roundtrip
    def test_numpy_numeric_array_uses_zarr_and_preserves_dtype(self):
        value = np.arange(24, dtype="float32").reshape(2, 3, 4)
        result = self.roundtrip(value)
        self.assertTrue((self.folder / "metadata.json.arrays/0.zarr/.zarray").is_file())
        np.testing.assert_array_equal(result, value)
        self.assertEqual(result.dtype, value.dtype)

    def test_mixed_xarray_metadata_and_coordinates(self):
        import xarray as xr
        value = xr.DataArray(np.array([["station", 1], [None, pd.Timestamp("2025-01-01")]], dtype=object),
                             dims=("row", "field"), coords={"row": [1, 2], "field": ["a", "b"]},
                             name="metadata", attrs={"units": "mixed"})
        xr.testing.assert_identical(self.roundtrip(value), value)

    def test_numeric_stage_roundtrip_remains_lazy(self):
        import xarray as xr
        import dask.array as da
        from atlas_actris.utils import export_processing_stage as stages
        value = xr.DataArray(da.arange(12, chunks=4).reshape((3, 4)), dims=("time", "bins"),
                             coords={"time": pd.date_range("2025-01-01", periods=3)}, name="signal")
        stages.export_processing_stage({"ray": {"signal": value}}, "complete", str(self.folder), print_estimated_size=False)
        actual = stages.import_processing_stage(str(self.folder), "complete")["ray"]["signal"]
        self.assertIsInstance(actual.data, da.Array)
        xr.testing.assert_identical(actual.compute(), value.compute())

    def test_mixed_dataset_numeric_variable_stays_lazy(self):
        import xarray as xr
        import dask.array as da
        value = xr.Dataset({"metadata": ("field", np.array(["site", None], dtype=object)),
                            "numbers": ("bins", da.arange(12, chunks=4))},
                           coords={"field": ["a", "b"], "bins": np.arange(12)},
                           attrs={"source": "test"})
        actual = self.roundtrip(value)
        self.assertIsInstance(actual["numbers"].data, da.Array)
        xr.testing.assert_identical(actual.compute(), value.compute())


if __name__ == "__main__":
    unittest.main()
