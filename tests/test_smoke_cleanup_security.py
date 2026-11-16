"""Smoke cleanup boundaries, using disposable fixtures and mocked processors."""
import contextlib
import io
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atlas_actris.testing import smoke_test as smoke
from atlas_actris.testing import signal_viewer_smoke_test as viewer


class SmokeCleanupSecurityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas-cleanup-test-")
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name).resolve()
        self.pack = self.base / "pack"
        self.pack.mkdir()
        (self.pack / smoke.DEFAULT_INI_NAME).write_text("[test]\n")
        self.case = "safe_case"

    def run_fixture(self, module, processor=None, **kwargs):
        with patch.object(module.subprocess, "run", side_effect=processor) as launch:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                status = module.run_smoke_test(self.pack, case_name=kwargs.pop("case_name", self.case), **kwargs)
        return status, launch

    def sentinel(self, path):
        path.mkdir(parents=True, exist_ok=True)
        target = path / "keep.txt"
        target.write_text("preserve")
        return target

    def link(self, target, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            path.symlink_to(target, target_is_directory=True)
        except (OSError, NotImplementedError) as error:
            self.skipTest(f"Directory symlinks unavailable: {error}")

    def test_unsafe_case_names_never_launch_or_delete(self):
        victim = self.sentinel(self.base / "victim" / "cache")
        for module in (smoke, viewer):
            for name in (str(victim.parent.parent), "../victim", "a/../../victim",
                         "a\\..\\victim", "C:\\victim", "C:victim", ".", "..", "", "case."):
                with self.subTest(module=module.__name__, case=name):
                    status, launch = self.run_fixture(module, case_name=name)
                    self.assertEqual(status, 1)
                    launch.assert_not_called()
                    self.assertTrue(victim.exists())

    def test_analysis_and_case_links_are_rejected(self):
        victim = self.sentinel(self.base / "victim")
        for module in (smoke, viewer):
            for link_name in ("analysis", "analysis/" + self.case):
                path = self.pack / link_name
                self.link(victim.parent, path)
                try:
                    status, launch = self.run_fixture(module)
                    self.assertEqual(status, 1)
                    launch.assert_not_called()
                    self.assertTrue(victim.exists())
                finally:
                    path.unlink()
                    analysis = self.pack / "analysis"
                    if analysis.is_dir() and not any(analysis.iterdir()):
                        analysis.rmdir()

    def test_nested_link_is_rejected_before_any_target_is_deleted(self):
        victim = self.sentinel(self.base / "victim")
        root = self.pack / "analysis" / self.case
        ordinary = self.sentinel(root / "signal_viewer")
        self.link(victim.parent, root / "cache")
        for module in (smoke, viewer):
            status, launch = self.run_fixture(module)
            self.assertEqual(status, 1)
            launch.assert_not_called()
            self.assertTrue(ordinary.exists())
            self.assertTrue(victim.exists())

    def test_main_cleanup_preserves_other_cases_and_raw_caches(self):
        unrelated = self.sentinel(self.pack / "analysis" / "other_case")
        raw_cache = self.sentinel(self.pack / "ray" / "cache")
        stale = self.sentinel(self.pack / "analysis" / self.case)
        def processor(*args, **kwargs):
            self.assertFalse(stale.exists())
            root = self.pack / "analysis" / self.case
            for name, count, suffix in (("plots", 10, ".png"), ("reports", 1, ".html"),
                                        ("reports", 1, ".docx"), ("ascii", 5, ".txt")):
                folder = root / name
                folder.mkdir(parents=True, exist_ok=True)
                for index in range(count):
                    (folder / (str(index) + suffix)).touch()
            return SimpleNamespace(returncode=0)
        status, launch = self.run_fixture(smoke, processor)
        self.assertEqual(status, 0)
        launch.assert_called_once()
        self.assertFalse(stale.parent.exists())
        self.assertTrue(unrelated.exists())
        self.assertTrue(raw_cache.exists())

    def test_viewer_cleanup_preserves_other_products(self):
        root = self.pack / "analysis" / self.case
        plot = self.sentinel(root / "plots")
        unrelated = self.sentinel(self.pack / "analysis" / "other_case")
        def processor(*args, **kwargs):
            self.sentinel(root / "signal_viewer")
            self.sentinel(root / "cache")
            return SimpleNamespace(returncode=0)
        status, _ = self.run_fixture(viewer, processor)
        self.assertEqual(status, 0)
        self.assertTrue(plot.exists())
        self.assertTrue(unrelated.exists())
        self.assertFalse((root / "signal_viewer").exists())
        self.assertFalse((root / "cache").exists())

    def test_keep_output_retains_existing_outputs(self):
        for module in (smoke, viewer):
            marker = self.sentinel(self.pack / "analysis" / self.case / "cache")
            status, _ = self.run_fixture(module, lambda *a, **k: SimpleNamespace(returncode=9), keep_output=True)
            self.assertEqual(status, 9)
            self.assertTrue(marker.exists())

    def test_cleanup_rechecks_links_created_by_the_subprocess(self):
        victim = self.sentinel(self.base / "victim")
        root = self.pack / "analysis" / self.case
        for module in (smoke, viewer):
            root.mkdir(parents=True, exist_ok=True)
            def processor(*args, **kwargs):
                self.link(victim.parent, root / "cache")
                return SimpleNamespace(returncode=9)
            status, _ = self.run_fixture(module, processor)
            self.assertEqual(status, 1)
            self.assertTrue(victim.exists())
            (root / "cache").unlink()


if __name__ == "__main__":
    unittest.main()
