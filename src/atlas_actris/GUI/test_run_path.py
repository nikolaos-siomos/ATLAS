"""Exercise editor-to-Run synchronization without a graphical display."""
from pathlib import Path
import tempfile
import tkinter as tk
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from atlas_schema_gui import SchemaTab
from editor_model import EditorModel
from gui_common import RunTab
from schema_loader import load_schemas

BASE = Path(__file__).resolve().parent


class RunPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = load_schemas(BASE.parent)["call_atlas.ini"]

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas run ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.first = self.folder / "first call.ini"
        self.second = self.folder / "second call.ini"
        for path in (self.first, self.second):
            path.write_bytes((BASE / "examples/call_atlas.ini").read_bytes())
        self.tcl = tk.Tcl()
        self.runner = SimpleNamespace(call_atlas_path=tk.StringVar(master=self.tcl),
                                      engine=Mock(busy=False), call_atlas_selected=None, browse_directory=None)
        self.runner.set_call_atlas_path = lambda path: RunTab.set_call_atlas_path(self.runner, path)
        self.runner.get_call_atlas_path = lambda: RunTab.get_call_atlas_path(self.runner)
        self.editor = SimpleNamespace(
            filename="call_atlas.ini", model=EditorModel(self.schema),
            path_changed=self.runner.set_call_atlas_path,
            browse_directory=None, variables={},
            include_descriptions=Mock(get=lambda: False),
            status=tk.StringVar(master=self.tcl), build=Mock(), refresh_values=Mock(), update_status=Mock(),
            _refreshing=False,
            may_discard=Mock(return_value=True),
        )
        self.editor.normalize_fields = lambda: SchemaTab.normalize_fields(self.editor)
        self.editor.notify_path_changed = lambda: SchemaTab.notify_path_changed(self.editor)

    def open_editor(self, path):
        with patch("atlas_schema_gui.ask_open_file", return_value=str(path)), \
                patch("atlas_schema_gui.LoadingDialog"):
            SchemaTab.open(self.editor)

    def test_open_and_explicit_reselection_update_run_path(self):
        self.open_editor(self.first)
        self.assertEqual(self.runner.get_call_atlas_path(), self.first)
        self.runner.call_atlas_path.set(str(self.second))
        self.open_editor(self.first)
        self.assertEqual(self.runner.get_call_atlas_path(), self.first)

    def test_run_override_survives_value_edits_and_normal_save(self):
        self.open_editor(self.first)
        self.runner.call_atlas_path.set(str(self.second))
        SchemaTab.changed(self.editor, ("configuration", "scc_configuration_id"), "713")
        self.assertTrue(SchemaTab.save(self.editor))
        self.assertEqual(self.runner.get_call_atlas_path(), self.second)
        target = self.folder / "saved as.ini"
        with patch("atlas_schema_gui.ask_save_file", return_value=str(target)):
            self.assertTrue(SchemaTab.save(self.editor, save_as=True))
        self.assertEqual(self.runner.get_call_atlas_path(), target)

    def test_new_clears_path_and_first_save_populates_it(self):
        self.open_editor(self.first)
        SchemaTab.new(self.editor)
        self.assertEqual(self.runner.call_atlas_path.get(), "")
        target = self.folder / "new call.ini"
        with patch("atlas_schema_gui.ask_save_file", return_value=str(target)):
            self.assertTrue(SchemaTab.save(self.editor))
        self.assertEqual(self.runner.get_call_atlas_path(), target)

    def test_failed_or_cancelled_selection_preserves_run_override(self):
        self.open_editor(self.first)
        self.runner.call_atlas_path.set(str(self.second))
        with patch("atlas_schema_gui.ask_open_file", return_value=""):
            SchemaTab.open(self.editor)
        with patch("atlas_schema_gui.messagebox.showerror") as error:
            self.open_editor(self.folder / "missing.ini")
            error.assert_called_once()
        with patch("atlas_schema_gui.ask_save_file", return_value=""):
            self.assertFalse(SchemaTab.save(self.editor, save_as=True))
        self.assertEqual(self.editor.model.path, self.first)
        self.assertEqual(self.runner.get_call_atlas_path(), self.second)

    def test_path_getter_uses_latest_text_and_resolves_relative_paths(self):
        self.open_editor(self.first)
        self.runner.call_atlas_path.set('"{}"'.format(self.second))
        self.assertEqual(self.runner.get_call_atlas_path(), self.second)
        with patch("gui_common.__file__", str(self.folder / "gui_common.py")):
            self.runner.call_atlas_path.set(self.first.name)
            self.assertEqual(self.runner.get_call_atlas_path(), self.first)
        for value in ("", str(self.folder), str(self.folder / "missing.ini")):
            with self.subTest(value=value):
                self.runner.call_atlas_path.set(value)
                with self.assertRaises(ValueError):
                    self.runner.get_call_atlas_path()

    def test_run_browse_changes_run_target_without_reloading_editor(self):
        self.open_editor(self.first)
        with patch("gui_common.ask_open_file", return_value=str(self.second)):
            RunTab.select_call_atlas(self.runner)
        self.assertEqual(self.runner.get_call_atlas_path(), self.second)
        self.assertEqual(self.editor.model.path, self.first)
        with patch("gui_common.ask_open_file", return_value=""):
            RunTab.select_call_atlas(self.runner)
        self.assertEqual(self.runner.get_call_atlas_path(), self.second)

    def test_run_browse_dispatches_linked_loading_without_changing_target_early(self):
        self.runner.set_call_atlas_path(self.first)
        self.runner.call_atlas_selected = Mock()
        with patch("gui_common.ask_open_file", return_value=str(self.second)):
            RunTab.select_call_atlas(self.runner)
        self.runner.call_atlas_selected.assert_called_once_with(self.second)
        self.assertEqual(self.runner.get_call_atlas_path(), self.first)

    def test_run_browse_is_blocked_during_processing(self):
        self.runner.engine.busy = True
        with patch("gui_common.ask_open_file") as browse, patch("gui_common.messagebox.showinfo"):
            RunTab.select_call_atlas(self.runner)
        browse.assert_not_called()


if __name__ == "__main__":
    unittest.main()
