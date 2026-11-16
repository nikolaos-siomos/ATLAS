"""Exercise actual cache control flow without importing format-specific readers.

Only the raw reader and Zarr I/O boundary are mocked. Function definitions are
compiled from the production source to keep these migration/security checks
independent of optional native reader dependencies.
"""
import ast
import contextlib
import io
import os
from pathlib import Path
import pickle
import re
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock

from atlas_actris.utils import safe_serialization as codec


class RawCacheMigrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas cache migration ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.output = self.folder / "output"; self.output.mkdir()
        self.raw = self.folder / "raw"; self.raw.mkdir()
        (self.raw / "measurement.dat").write_bytes(b"original measurement")
        source = Path(__file__).resolve().parents[1] / "src/atlas_actris/readers/read_raw_lidar_files_cache.py"
        tree = ast.parse(source.read_text())
        body = [ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0)]
        body += [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        scope = dict(os=os, Path=Path, re=re, shutil=shutil, FileReaderError=RuntimeError,
                     read_json=codec.read_json, write_json=codec.write_json,
                     checked_reference=codec.checked_reference, print_header=lambda *a: None,
                     print_subsection=lambda *a: None, endpoint=lambda *a: None)
        exec(compile(ast.fix_missing_locations(ast.Module(body=body, type_ignores=[])), str(source), "exec"), scope)
        self.cache = SimpleNamespace(**{key: value for key, value in scope.items() if callable(value)})
        self.scope = scope
        self.config = {"output_folder": str(self.output), "paths": {"ray": str(self.raw)}, "raw_file_format": "fixture"}
        self.fixture = ({"station": "test"}, {"channel": "a"}, {"time": "test"}, "signal", "shots")
        self.reader = Mock(return_value=self.fixture)
        scope["reader_menu"] = {"fixture": self.reader}
        def write_entry(cache_root, meas_key, **kwargs):
            entry = Path(cache_root) / meas_key
            (entry / "data.zarr").mkdir(parents=True)
            codec.write_metadata({"mtype": kwargs["meas_type"]}, entry / "metadata.json")
        scope["_write_cache_entry"] = Mock(side_effect=write_entry)
        scope["_read_cache_entry"] = Mock(return_value=self.fixture + ("nrm",))

    def run_reader(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return self.cache.flexible_reader(self.config, downscale=False, chunk=False)

    def legacy_cache(self):
        marker = self.folder / "pickle_executed"
        class Payload:
            def __reduce__(self): return (Path.touch, (marker,))
        entry = self.output / "cache/ray"; (entry / "data.zarr").mkdir(parents=True)
        payload = pickle.dumps(Payload())
        (entry / "metadata.pkl").write_bytes(payload)
        (entry.parent / "source_manifest.pkl").write_bytes(payload)
        return marker

    def test_legacy_cache_rebuilt_without_executing_payload_then_reused(self):
        marker = self.legacy_cache()
        result = self.run_reader()
        self.assertEqual(result[0], {"ray": "signal"})
        self.reader.assert_called_once()
        self.assertFalse(marker.exists())
        self.assertTrue((self.output / "cache/source_manifest.json").is_file())
        self.assertFalse((self.output / "cache/source_manifest.pkl").exists())
        self.assertEqual((self.raw / "measurement.dat").read_bytes(), b"original measurement")
        self.reader.reset_mock()
        self.run_reader()
        self.reader.assert_not_called()

    def test_missing_original_measurements_leave_legacy_cache_untouched(self):
        marker = self.legacy_cache()
        self.config["paths"]["ray"] = str(self.folder / "missing raw")
        with self.assertRaisesRegex(RuntimeError, "original measurement path is missing"):
            self.run_reader()
        self.assertFalse(marker.exists())
        self.assertTrue((self.output / "cache/source_manifest.pkl").is_file())
        self.reader.assert_not_called()

    def test_changed_measurement_rebuilds_new_cache(self):
        self.run_reader()
        self.reader.reset_mock()
        (self.raw / "measurement.dat").write_bytes(b"changed original measurement")
        self.run_reader()
        self.reader.assert_called_once()

    def test_invalid_json_manifest_forces_rebuild(self):
        self.run_reader()
        self.reader.reset_mock()
        (self.output / "cache/source_manifest.json").write_text("invalid json")
        self.run_reader()
        self.reader.assert_called_once()


if __name__ == "__main__":
    unittest.main()
