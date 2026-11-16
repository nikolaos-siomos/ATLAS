"""Check GUI action isolation and preservation of caches when exports fail."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atlas_actris.GUI.atlas_runner import post_run_environment, ProcessRunner
from atlas_actris.utils.gui_post_run import run_gui_post_run, validate_gui_post_run
from atlas_actris.utils.cleaners import ask_clean_cache


class GuiPostRunTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cache = self.root / "cache"
        self.cache.mkdir()
        (self.cache / "measurement").write_text("cached measurement")
        self.info = {"output_folder": str(self.root), "export_stages": ["pol_cal_complete"]}

    def run_actions(self, clean=False, action="Default", export=None, delete=None):
        export = export if export is not None else Mock()
        delete = delete if delete is not None else Mock()
        with patch.dict(os.environ, post_run_environment(clean, action)):
            run_gui_post_run(object(), self.info, export, delete, ask_clean_cache)
        return export, delete

    def test_defaults_preserve_cache_and_stages_despite_inherited_answers(self):
        with patch.dict(os.environ, {"ATLAS_CLEAN_CACHE_ANSWER": "Y",
                                    "ATLAS_DELETE_EXPORTED_ANSWER": "Y",
                                    "ATLAS_EXPORT_STAGE_ANSWER": "Y"}):
            export, delete = self.run_actions()
        export.assert_not_called()
        delete.assert_not_called()
        self.assertTrue(self.cache.is_dir())

    def test_export_finishes_before_cache_cleanup(self):
        def export(**kwargs):
            self.assertTrue(self.cache.is_dir())
            self.assertEqual(kwargs["stage_names"], ["pol_cal_complete"])
            self.assertFalse(kwargs["ask"])
            self.assertTrue(kwargs["overwrite"])
        self.run_actions(True, "Export stages", export=export)
        self.assertFalse(self.cache.exists())

    def test_failed_export_preserves_cache_and_does_not_delete_stages(self):
        delete = Mock()
        with self.assertRaisesRegex(OSError, "disk full"):
            self.run_actions(True, "Export stages", export=Mock(side_effect=OSError("disk full")), delete=delete)
        self.assertTrue((self.cache / "measurement").is_file())
        delete.assert_not_called()

    def test_clean_stages_uses_current_output_without_export_or_cache_cleanup(self):
        export, delete = self.run_actions(action="Clean stages")
        delete.assert_called_once_with(output_folder=str(self.root), ask=False)
        export.assert_not_called()
        self.assertTrue(self.cache.is_dir())

    def test_empty_export_selection_rejected_before_actions(self):
        self.info["export_stages"] = []
        export, delete = Mock(), Mock()
        with self.assertRaisesRegex(ValueError, "nonempty export_stages"):
            self.run_actions(True, "Export stages", export=export, delete=delete)
        export.assert_not_called()
        delete.assert_not_called()
        self.assertTrue(self.cache.is_dir())

    def test_terminal_validation_is_unchanged(self):
        with patch.dict(os.environ, {}, clear=True):
            validate_gui_post_run({})

    def test_invalid_action_never_starts_a_process(self):
        with patch("atlas_actris.GUI.atlas_runner.subprocess.Popen") as launch:
            with self.assertRaises(ValueError):
                ProcessRunner().start([sys.executable], self.root, stage_action="unexpected")
        launch.assert_not_called()


if __name__ == "__main__":
    unittest.main()
