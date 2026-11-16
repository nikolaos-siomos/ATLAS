"""Folder listing, asynchronous loading and mouse-wheel regression checks."""
import ast
from pathlib import Path
import queue
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from directory_browser import (directory_entries, selected_folder, selected_file,
                               matching_entries, DirectoryDialog, background_call, is_hidden)
from gui_common import ScrollFrame
from atlas_schema_gui import App, SchemaTab, configure_fonts


class BrowserTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas browser ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()

    def test_list_contains_files_including_inis_and_sorts_folders_first(self):
        (self.folder / "z folder").mkdir()
        (self.folder / "call_atlas.ini").write_text("")
        (self.folder / "notes.txt").write_text("")
        (self.folder / ".hidden").write_text("")
        folder, rows = directory_entries(self.folder, show_hidden=True)
        self.assertEqual(folder, self.folder)
        self.assertEqual(rows[0][0], "z folder")
        types = {row[0]: row[3] for row in rows}
        self.assertEqual(types["call_atlas.ini"], "INI file")
        self.assertEqual(types["notes.txt"], "File")
        self.assertIn(".hidden", types)

    def test_selecting_file_returns_current_folder_instead_of_file(self):
        file = self.folder / "call_atlas.ini"; file.touch()
        self.assertEqual(selected_folder(self.folder, str(self.folder), file), self.folder)
        with self.assertRaises(ValueError):
            selected_folder(self.folder, str(file))

    def test_hidden_files_and_directories_are_hidden_by_default(self):
        (self.folder / ".hidden.ini").touch()
        (self.folder / ".hidden folder").mkdir()
        (self.folder / "config_file.ini").touch()
        _, entries = directory_entries(self.folder)
        self.assertEqual([row[0] for row in entries], ["config_file.ini"])
        _, entries = directory_entries(self.folder, show_hidden=True)
        self.assertEqual({row[0] for row in entries}, {".hidden.ini", ".hidden folder", "config_file.ini"})

    def test_native_hidden_attributes_are_recognized(self):
        entry = SimpleNamespace(name="hidden.ini", stat=Mock(return_value=SimpleNamespace(st_file_attributes=2)))
        with patch("directory_browser.os.name", "nt"):
            self.assertTrue(is_hidden(entry))
            entry.stat.return_value.st_file_attributes = 0
            self.assertFalse(is_hidden(entry))
        entry.stat.return_value = SimpleNamespace(st_flags=0x8000)
        with patch("directory_browser.sys.platform", "darwin"), patch("directory_browser.stat.UF_HIDDEN", 0x8000, create=True):
            self.assertTrue(is_hidden(entry))

    def test_value_refresh_reuses_widgets_and_suppresses_per_field_status_updates(self):
        variables = {("section", "a"): Mock(get=lambda: "old"),
                     ("section", "b"): Mock(get=lambda: "same")}
        model = SimpleNamespace(values={("section", "a"): "new", ("section", "b"): "same"})
        tab = SimpleNamespace(variables=variables, model=model, _refreshing=False,
                              build=Mock(), update_status=Mock())
        variables[("section", "a")].set.side_effect = lambda value: SchemaTab.changed(tab, ("section", "a"), value)
        SchemaTab.refresh_values(tab)
        self.assertIs(tab.variables, variables)
        self.assertFalse(tab._refreshing)
        variables[("section", "a")].set.assert_called_once_with("new")
        variables[("section", "b")].set.assert_not_called()
        tab.update_status.assert_called_once()
        tab.build.assert_not_called()

    def test_selected_subdirectory_and_typed_paths(self):
        child = self.folder / "INI files"; child.mkdir()
        self.assertEqual(selected_folder(self.folder, str(self.folder), child), child)
        self.assertEqual(selected_folder(self.folder, "INI files"), child)
        self.assertEqual(selected_folder(self.folder, str(self.folder)), self.folder)
        missing = self.folder / "new output"
        self.assertEqual(selected_folder(self.folder, str(missing), mustexist=False), missing)
        self.assertFalse(missing.exists())
        with self.assertRaises(ValueError): selected_folder(self.folder, str(missing))

    def test_open_file_requires_existing_regular_file_and_trims_edges(self):
        path = self.folder / "config file.ini"
        path.touch()
        self.assertEqual(selected_file(self.folder, "  config file.ini  "), path)
        for name in ("", "missing.ini", str(self.folder)):
            with self.assertRaises(ValueError):
                selected_file(self.folder, name)

    def test_save_file_adds_extension_without_creating_anything(self):
        path = self.folder / "new config.ini"
        self.assertEqual(selected_file(self.folder, "new config", save=True, defaultextension=".ini"), path)
        self.assertFalse(path.exists())
        self.assertEqual(selected_file(self.folder, "config.INI", save=True, defaultextension=".ini"),
                         self.folder / "config.INI")
        for name in (str(self.folder), "missing/config.ini"):
            with self.assertRaises(ValueError):
                selected_file(self.folder, name, save=True, defaultextension=".ini")

    def test_file_filters_keep_navigation_and_respect_hidden_setting(self):
        (self.folder / "measurements").mkdir()
        (self.folder / "CONFIG.INI").touch()
        (self.folder / "sounding.nc").touch()
        (self.folder / ".private.ini").touch()
        _, entries = directory_entries(self.folder)
        self.assertEqual({row[0] for row in matching_entries(entries, "*.ini")},
                         {"measurements", "CONFIG.INI"})
        self.assertEqual({row[0] for row in matching_entries(entries, ("*.ini", "*.nc"))},
                         {"measurements", "CONFIG.INI", "sounding.nc"})
        _, entries = directory_entries(self.folder, show_hidden=True)
        self.assertIn(".private.ini", {row[0] for row in matching_entries(entries, "*.ini")})

    def test_save_replacement_requires_confirmation_and_never_writes_in_chooser(self):
        path = self.folder / "call_atlas.ini"
        path.write_text("existing contents")
        dialog = SimpleNamespace(busy=False, mode="save", current=self.folder,
                                 filename=Mock(get=lambda: path.name), defaultextension=".ini",
                                 result="", destroy=Mock())
        with patch("directory_browser.messagebox.askyesno", return_value=False) as confirm:
            DirectoryDialog.accept(dialog)
        confirm.assert_called_once()
        dialog.destroy.assert_not_called()
        self.assertEqual(dialog.result, "")
        with patch("directory_browser.messagebox.askyesno", return_value=True):
            DirectoryDialog.accept(dialog)
        self.assertEqual(dialog.result, str(path))
        dialog.destroy.assert_called_once()
        self.assertEqual(path.read_text(), "existing contents")

    def test_worker_does_not_block_and_reports_errors_as_data(self):
        release = threading.Event()
        self.addCleanup(release.set)
        result = background_call(lambda: release.wait(3))
        with self.assertRaises(queue.Empty): result.get_nowait()
        release.set()
        self.assertEqual(result.get(timeout=5), (True, True))
        def fail(): raise OSError("Cannot read folder")
        success, error = background_call(fail).get(timeout=5)
        self.assertFalse(success)
        self.assertIsInstance(error, OSError)

    def test_pending_ini_load_keeps_polling_without_hiding_message(self):
        dialog = Mock()
        app = SimpleNamespace(pending_load=queue.Queue(), after=Mock(), loading_dialog=dialog)
        app.poll_folder_load = Mock()
        App.poll_folder_load(app)
        app.after.assert_called_once_with(30, app.poll_folder_load)
        dialog.close.assert_not_called()

    def test_failed_ini_load_closes_progress_and_preserves_current_models(self):
        result = queue.Queue(); result.put((False, OSError("No permission")))
        dialog = Mock()
        app = SimpleNamespace(pending_load=result, loading_dialog=dialog,
                              working_folder=self.folder, _loading=True)
        app.end_folder_load = lambda: App.end_folder_load(app)
        with patch("atlas_schema_gui.messagebox.showerror") as show:
            App.poll_folder_load(app)
        self.assertFalse(app._loading)
        dialog.close.assert_called_once()
        self.assertEqual(app.working_folder, self.folder)
        show.assert_called_once()

    def test_mouse_wheel_scrolls_page_and_stops_combobox_class_binding(self):
        page = SimpleNamespace(canvas=Mock())
        for event, direction in ((SimpleNamespace(delta=-120), 1),
                                 (SimpleNamespace(delta=120), -1),
                                 (SimpleNamespace(num=4), -1),
                                 (SimpleNamespace(num=5), 1)):
            page.canvas.reset_mock()
            self.assertEqual(ScrollFrame.scroll_wheel(page, event), "break")
            page.canvas.yview_scroll.assert_called_once_with(direction, "units")

    def test_font_resize_does_not_override_theme_fonts_or_families(self):
        with patch("atlas_schema_gui.tkfont.Font") as font, patch("atlas_schema_gui.ttk.Style") as style:
            configure_fonts(Mock())
        for call in style.return_value.configure.call_args_list:
            self.assertNotIn("font", call.kwargs)
        for call in font.return_value.configure.call_args_list:
            self.assertNotIn("family", call.kwargs)
            self.assertNotIn("weight", call.kwargs)
        self.assertIn(unittest.mock.call(size=12), font.return_value.configure.call_args_list)

    def test_button_labels_do_not_depend_on_unicode_symbols(self):
        for name in ("atlas_schema_gui.py", "gui_common.py", "directory_browser.py"):
            tree = ast.parse(Path(__file__).with_name(name).read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "Button":
                    for argument in node.keywords:
                        if argument.arg == "text" and isinstance(argument.value, ast.Constant):
                            self.assertTrue(argument.value.value.isascii(), argument.value.value)


if __name__ == "__main__":
    unittest.main()
