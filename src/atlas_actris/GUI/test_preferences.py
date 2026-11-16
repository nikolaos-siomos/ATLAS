"""Directory memory must not require writable application files or load INIs."""
import json
from pathlib import Path
import tempfile
from types import MethodType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from atlas_schema_gui import App, SchemaTab
from gui_common import RunTab
from preferences import (defaults, preferences_path, load_preferences,
                         save_preferences, remembered_directory, MAX_BYTES)
from schema_loader import load_schemas


class PreferenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas GUI preferences ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.path = self.folder / "preferences" / "gui_preferences.json"

    def test_platform_locations_use_user_data_and_ignore_relative_overrides(self):
        home, config = self.folder / "home", self.folder / "config"
        suffix = Path("atlas_actris/gui_preferences.json")
        self.assertEqual(preferences_path("win32", {"LOCALAPPDATA": str(config)}, home), config / suffix)
        self.assertEqual(preferences_path("darwin", {}, home), home / "Library/Application Support" / suffix)
        self.assertEqual(preferences_path("linux", {"XDG_CONFIG_HOME": str(config)}, home), config / suffix)
        self.assertEqual(preferences_path("linux", {"XDG_CONFIG_HOME": "relative"}, home), home / ".config" / suffix)
        self.assertEqual(preferences_path("win32", {"LOCALAPPDATA": "relative"}, home), home / "AppData/Local" / suffix)

    def test_only_directory_is_saved_not_ini_selection_or_other_values(self):
        values = dict(defaults(), last_directory=str(self.folder), last_init_file="custom.ini",
                      reopen_last=True, password="not-a-stored-preference")
        self.assertTrue(save_preferences(values, self.path))
        self.assertEqual(load_preferences(self.path), dict(defaults(), last_directory=str(self.folder)))
        self.assertEqual(set(json.loads(self.path.read_text())), {"format", "last_directory"})
        self.assertEqual(remembered_directory(load_preferences(self.path)), self.folder)

    def test_missing_unreadable_corrupt_or_oversize_preferences_fall_back(self):
        self.assertEqual(load_preferences(self.path), defaults())
        self.assertFalse(self.path.parent.exists())
        with patch.object(Path, "open", side_effect=PermissionError("no access")):
            self.assertEqual(load_preferences(self.path), defaults())
        self.path.parent.mkdir()
        for data in (b"broken JSON", b"[]", b'{"format":"unknown"}', b"x" * (MAX_BYTES + 1), b"\xff"):
            self.path.write_bytes(data)
            self.assertEqual(load_preferences(self.path), defaults())

    def test_invalid_missing_or_inaccessible_directories_are_not_used(self):
        self.path.parent.mkdir()
        for value in (123, "relative", "bad\x00path"):
            self.path.write_text(json.dumps(dict(defaults(), last_directory=value)))
            self.assertEqual(load_preferences(self.path), defaults())
        self.assertIsNone(remembered_directory({"last_directory": str(self.folder / "missing")}))
        with patch.object(Path, "is_dir", side_effect=PermissionError("no access")):
            self.assertIsNone(remembered_directory({"last_directory": str(self.folder)}))
        with patch("preferences.os.access", return_value=False):
            self.assertIsNone(remembered_directory({"last_directory": str(self.folder)}))

    def test_unwritable_preference_directory_does_not_raise(self):
        with patch.object(Path, "mkdir", side_effect=PermissionError("read only")):
            self.assertFalse(save_preferences(defaults(), self.path))
        self.assertFalse(self.path.exists())

    def test_failed_atomic_replacement_preserves_previous_preferences(self):
        self.assertTrue(save_preferences(defaults(), self.path))
        before = self.path.read_bytes()
        with patch("preferences.os.replace", side_effect=PermissionError("locked")):
            self.assertFalse(save_preferences(dict(defaults(), last_directory=str(self.folder)), self.path))
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(list(self.path.parent.iterdir()), [self.path])

    def test_remembering_directory_tolerates_failed_write_and_syncs_run_path(self):
        app = SimpleNamespace(preferences=defaults(), directory_memory_enabled=True,
                              _loading=False, runner=Mock())
        app.remember_ini_directory = MethodType(App.remember_ini_directory, app)
        target = self.folder / "custom run.ini"
        with patch("atlas_schema_gui.save_preferences", return_value=False):
            App.init_path_changed(app, target)
        self.assertEqual(app.preferences["last_directory"], str(self.folder))
        app.runner.set_call_atlas_path.assert_called_once_with(target)

    def test_working_folder_browser_starts_at_remembered_path_without_selecting(self):
        app = SimpleNamespace(_loading=False, runner=Mock(), working_folder=None,
                              preferences=dict(defaults(), last_directory=str(self.folder)),
                              editors=[], begin_folder_load=Mock())
        app.runner.running.return_value = False
        with patch("atlas_schema_gui.ask_directory", return_value="") as browse:
            App.select_working_folder(app)
        self.assertEqual(browse.call_args.kwargs["initialdir"], self.folder)
        app.begin_folder_load.assert_not_called()

    def test_run_browser_uses_remembered_directory_without_loading_on_cancel(self):
        tab = SimpleNamespace(engine=Mock(busy=False), get_call_atlas_path=Mock(side_effect=ValueError()),
                              browse_directory=lambda: self.folder, call_atlas_selected=Mock(),
                              set_call_atlas_path=Mock())
        with patch("gui_common.ask_open_file", return_value="") as browse:
            RunTab.select_call_atlas(tab)
        self.assertEqual(browse.call_args.kwargs["initialdir"], str(self.folder))
        tab.call_atlas_selected.assert_not_called()
        tab.set_call_atlas_path.assert_not_called()

    def test_parameter_browser_uses_directory_hint_without_changing_value(self):
        schemas = load_schemas(Path(__file__).resolve().parents[1])
        field = next(f for f in schemas["call_atlas.ini"].fields if f.key == "output_folder")
        variable = Mock()
        tab = SimpleNamespace(model=SimpleNamespace(path=None), browse_directory=lambda: self.folder,
                              variables={field.identity: variable})
        with patch("atlas_schema_gui.ask_directory", return_value="") as browse:
            SchemaTab.browse_value(tab, field)
        self.assertEqual(browse.call_args.kwargs["initialdir"], self.folder)
        variable.set.assert_not_called()


if __name__ == "__main__":
    unittest.main()
