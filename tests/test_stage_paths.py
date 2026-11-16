"""Stage path regressions using temporary files and trusted local metadata."""
import contextlib
import io
from pathlib import Path
import json
import pickle
import tempfile
import unittest
from unittest.mock import patch

from atlas_actris.utils import export_processing_stage as stages


class StagePathTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas stage paths ")
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name).resolve()
        self.output = self.base / "output"
        self.output.mkdir()
        self.sentinel = self.output / "keep.txt"
        self.sentinel.write_text("keep")
        self.data = {"ray": {"metadata": {"station": "test"}}}

    def export(self, name="complete", **kwargs):
        return stages.export_processing_stage(
            self.data, stage_name=name, output_folder=str(self.output),
            print_estimated_size=False, **kwargs)

    def symlink(self, link, target, directory=False):
        try:
            link.symlink_to(target, target_is_directory=directory)
        except (OSError, NotImplementedError) as error:
            self.skipTest(f"Symlinks unavailable: {error}")

    def test_unsafe_stage_names_do_not_create_export_root_or_delete_files(self):
        for name in ("", " ", ".", "..", "../other", "a/b", "a\\b", "C:other", "C:\\other", "\\\\host\\share", "bad\x00name"):
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    self.export(name, overwrite=True)
                with self.assertRaises(ValueError):
                    stages.delete_exported_stage(str(self.output), name)
                self.assertEqual(self.sentinel.read_text(), "keep")
                self.assertFalse((self.output / "exported").exists())

    def test_normal_roundtrip_listing_overwrite_and_delete(self):
        path = Path(self.export())
        self.assertEqual(path, self.output / "exported/complete")
        self.assertEqual(stages.list_exported_stages(str(self.output)), ["complete"])
        self.assertEqual(stages.import_processing_stage(str(self.output), "complete"), self.data)
        with self.assertRaises(FileExistsError):
            self.export()
        self.data = {"ray": {"metadata": {"station": "updated"}}}
        self.export(overwrite=True)
        self.assertEqual(stages.import_processing_stage(str(self.output), "complete"), self.data)
        with contextlib.redirect_stdout(io.StringIO()):
            stages.delete_exported_stage(str(self.output), "complete")
        self.assertFalse(path.exists())
        self.assertEqual(self.sentinel.read_text(), "keep")

    def test_declining_confirmations_keeps_existing_files(self):
        path = Path(self.export())
        with patch.object(stages, "_ask_yes_no", return_value=False):
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertIsNone(self.export(overwrite=True, ask=True))
                self.assertIsNone(stages.delete_exported_stage(str(self.output), "complete", ask=True))
                self.assertIsNone(stages.delete_all_exported_stages(str(self.output)))
        self.assertTrue((path / "manifest.json").exists())

    def test_invalid_entry_keys_do_not_overwrite_existing_stage(self):
        path = Path(self.export())
        before = (path / "manifest.json").read_bytes()
        for data in ({"..": {"metadata": 1}}, {"ray": {"../metadata": 1}}, {"ray": {"..": 1}}):
            self.data = data
            with self.assertRaises(ValueError):
                self.export(overwrite=True)
            self.assertEqual((path / "manifest.json").read_bytes(), before)

    def test_batch_invalid_later_stage_does_not_overwrite_first(self):
        path = Path(self.export())
        before = (path / "manifest.json").read_bytes()
        for data in ({"complete": self.data, "..": self.data},
                     {"complete": self.data, "second": {"ray": {"..": 1}}}):
            with self.assertRaises(ValueError):
                stages.export_processing_stages(data, str(self.output),
                    stage_names=list(data), overwrite=True, print_estimated_size=False)
            self.assertEqual((path / "manifest.json").read_bytes(), before)
            self.assertFalse((self.output / "exported/second").exists())

    def test_external_export_root_link_is_rejected(self):
        outside = self.base / "outside"; outside.mkdir()
        victim = outside / "keep.txt"; victim.write_text("outside")
        self.symlink(self.output / "exported", outside, directory=True)
        for operation in (lambda: self.export(),
                          lambda: stages.delete_all_exported_stages(str(self.output), ask=False),
                          lambda: stages.list_exported_stages(str(self.output))):
            with self.assertRaises(ValueError):
                operation()
        self.assertEqual(victim.read_text(), "outside")

    def test_external_stage_link_is_rejected_for_import_overwrite_delete(self):
        root = self.output / "exported"; root.mkdir()
        outside = self.base / "outside"; outside.mkdir()
        victim = outside / "keep.txt"; victim.write_text("outside")
        self.symlink(root / "complete", outside, directory=True)
        for operation in (lambda: self.export(overwrite=True),
                          lambda: stages.import_processing_stage(str(self.output), "complete"),
                          lambda: stages.delete_exported_stage(str(self.output), "complete")):
            with self.assertRaises(ValueError):
                operation()
        self.assertEqual(victim.read_text(), "outside")

    def test_nested_escaping_link_is_rejected_before_deletion_or_overwrite(self):
        path = Path(self.export())
        outside = self.base / "outside"; outside.mkdir()
        self.symlink(path / "escape", outside, directory=True)
        for operation in (lambda: self.export(overwrite=True),
                          lambda: stages.delete_exported_stage(str(self.output), "complete"),
                          lambda: stages.delete_all_exported_stages(str(self.output), ask=False)):
            with self.assertRaises(ValueError):
                operation()
            self.assertTrue((path / "manifest.json").exists())
            self.assertTrue((path / "ray/metadata/data.json").exists())

    def test_manifest_paths_rejected_before_entry_read(self):
        path = Path(self.export())
        manifest_path = path / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        for invalid in (str(self.sentinel), "../../keep.txt", "ray/../../../keep.txt", "C:\\keep.txt", "C:keep.txt", "\\\\host\\share\\data.json", "", "."):
            manifest["entries"]["ray"]["metadata"]["path"] = invalid
            manifest_path.write_text(json.dumps(manifest))
            with self.subTest(path=invalid):
                with patch.object(stages, "_read_export_entry") as read:
                    with self.assertRaises(ValueError):
                        stages.import_processing_stage(str(self.output), "complete")
                    read.assert_not_called()
        self.assertEqual(self.sentinel.read_text(), "keep")

    def test_entry_and_manifest_file_links_cannot_escape(self):
        path = Path(self.export())
        entry = path / "ray/metadata/data.json"; entry.unlink()
        self.symlink(entry, self.sentinel)
        with self.assertRaises(ValueError):
            stages.import_processing_stage(str(self.output), "complete")
        manifest = path / "manifest.json"; manifest.unlink()
        self.symlink(manifest, self.sentinel)
        with patch.object(stages, "read_json") as read:
            with self.assertRaises(ValueError):
                stages.import_processing_stage(str(self.output), "complete")
            read.assert_not_called()

    def test_zarr_escape_is_rejected_before_opening_store(self):
        path = Path(self.export())
        store = path / "ray/metadata/data.zarr"; store.mkdir()
        self.symlink(store / "chunk", self.sentinel)
        with patch.object(stages, "_read_xarray_from_zarr") as read:
            with self.assertRaises(ValueError):
                stages._read_export_entry(str(path), {"kind": "zarr_dataset", "path": "ray/metadata/data.zarr"})
            read.assert_not_called()

    def test_explicit_output_folder_symlink_still_supported(self):
        alias = self.base / "output alias"
        self.symlink(alias, self.output, directory=True)
        stages.export_processing_stage(self.data, "complete", str(alias), print_estimated_size=False)
        self.assertEqual(stages.import_processing_stage(str(alias), "complete"), self.data)

    def test_contained_entry_link_and_windows_relative_manifest_paths(self):
        path = Path(self.export())
        original = path / "ray/metadata/data.json"
        alias = path / "alias.json"
        self.symlink(alias, original)
        self.assertEqual(stages._read_export_entry(str(path), {"kind": "json", "path": "alias.json"}), self.data["ray"]["metadata"])
        self.assertEqual(stages._read_export_entry(str(path), {"kind": "json", "path": "ray\\metadata\\data.json"}), self.data["ray"]["metadata"])

    def test_reading_missing_root_has_no_filesystem_side_effect(self):
        self.assertEqual(stages.list_exported_stages(str(self.output)), [])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertIsNone(stages.delete_all_exported_stages(str(self.output), ask=False))
        self.assertFalse((self.output / "exported").exists())

    def test_legacy_stage_refused_without_opening_pickle(self):
        marker = self.base / "executed"
        class Payload:
            def __reduce__(self): return (Path.touch, (marker,))
        stage = self.output / "exported/legacy"; stage.mkdir(parents=True)
        (stage / "manifest.pkl").write_bytes(pickle.dumps(Payload()))
        self.assertEqual(stages.list_exported_stages(str(self.output)), ["legacy"])
        with self.assertRaisesRegex(ValueError, "Legacy pickle"):
            stages.import_processing_stage(str(self.output), "legacy")
        self.assertFalse(marker.exists())

    def test_unsupported_metadata_does_not_overwrite_existing_stage(self):
        stage = Path(self.export())
        original = (stage / "manifest.json").read_bytes()
        class Unsupported:
            def __reduce__(self): raise AssertionError("Never invoke pickle")
        self.data = {"ray": {"unknown": Unsupported()}}
        with self.assertRaises(ValueError):
            self.export(overwrite=True)
        self.assertEqual((stage / "manifest.json").read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
