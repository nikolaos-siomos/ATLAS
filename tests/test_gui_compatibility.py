"""Portable GUI entry-point checks; native window check skips without a display."""
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from atlas_actris import cli_gui

try:
    from atlas_actris.GUI import atlas_schema_gui as gui
    from atlas_actris.GUI import atlas_runner
except ModuleNotFoundError as error:
    if error.name not in ("tkinter", "_tkinter"):
        raise
    gui = None


class GuiCliTests(unittest.TestCase):
    def test_linux_display_bypasses_xim(self):
        with patch.object(cli_gui.sys, "platform", "linux"), \
                patch.dict(os.environ, {"DISPLAY": ":1", "XMODIFIERS": "@im=ibus"}, clear=True):
            cli_gui.configure_input_method()
            self.assertEqual(os.environ["XMODIFIERS"], "@im=none")

    def test_other_platforms_and_headless_runs_preserve_input_method(self):
        for platform, display in (("win32", ":1"), ("darwin", ":1"), ("linux", "")):
            with self.subTest(platform=platform, display=display), \
                    patch.object(cli_gui.sys, "platform", platform), \
                    patch.dict(os.environ, {"DISPLAY": display, "XMODIFIERS": "@im=ibus"}, clear=True):
                cli_gui.configure_input_method()
                self.assertEqual(os.environ["XMODIFIERS"], "@im=ibus")

    def test_help_needs_no_window(self):
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                cli_gui.main(["--help"])
        self.assertEqual(caught.exception.code, 0)


@unittest.skipIf(gui is None, "Tkinter is not installed in this Python environment")
class GuiCompatibilityTests(unittest.TestCase):
    def test_gui_launch_applies_input_method_choice_before_starting(self):
        for args, expected in (([], "@im=none"), (["--system-input-method"], "@im=ibus")):
            with self.subTest(args=args), patch.object(cli_gui.sys, "platform", "linux"), \
                    patch.dict(os.environ, {"DISPLAY": ":1", "XMODIFIERS": "@im=ibus"}, clear=True), \
                    patch.object(gui, "main", side_effect=lambda: os.environ["XMODIFIERS"]):
                self.assertEqual(cli_gui.main(args), expected)

    def test_default_desktop_checks_skip_before_creating_any_window(self):
        if os.environ.get("ATLAS_RUN_GUI_TESTS") == "1":
            self.skipTest("Default-mode safeguard; desktop checks were explicitly enabled")
        for name in ("test_native_window_initializes_when_display_available",
                     "test_native_folder_list_and_wheel_preserve_values_when_display_available"):
            result = unittest.TestResult()
            with patch.object(gui.tk, "Tk", side_effect=AssertionError("Unexpected window")) as create:
                GuiCompatibilityTests(name).run(result)
            create.assert_not_called()
            self.assertEqual(len(result.skipped), 1)

    def test_startup_messages_cover_schema_loading_and_editor_creation(self):
        events = []
        app = SimpleNamespace(loading_dialog=MagicMock(), select_working_folder=MagicMock(),
                              select_call_atlas=MagicMock(),
                              call_atlas_directory=MagicMock(return_value=None),
                              ready_to_run=MagicMock(), install_schemas=lambda schemas: events.append("editors"))
        app.loading_dialog.show_message.side_effect = events.append
        def schemas(source):
            events.append("schemas")
            return {}
        with patch.object(gui, "configure_fonts"), patch.object(gui, "configure_section_tabs"), \
                patch.object(gui, "load_schemas", side_effect=schemas), \
                patch.object(gui.tk, "BooleanVar"), patch.object(gui.tk, "StringVar"), \
                patch.object(gui.ttk, "Frame"), patch.object(gui.ttk, "Button"), \
                patch.object(gui.ttk, "Label"), patch.object(gui.ttk, "Notebook"), \
                patch.object(gui, "RunTab"):
            gui.App.initialize_interface(app)
        self.assertLess(events.index("Reading INI schemas and descriptions..."), events.index("schemas"))
        self.assertLess(events.index("Preparing interface fonts and styles..."), events.index("editors"))
    def test_font_info_mode_does_not_launch_editor(self):
        with patch.dict(os.environ), patch.object(cli_gui, "print_font_info", return_value=0) as inspect_fonts, \
                patch.object(gui, "main") as launch:
            self.assertEqual(cli_gui.main(["--font-info"]), 0)
        inspect_fonts.assert_called_once_with()
        launch.assert_not_called()

    def test_font_info_reports_renderer_and_conda_build_without_changing_fonts(self):
        with tempfile.TemporaryDirectory() as temporary:
            metadata = Path(temporary) / "conda-meta"
            metadata.mkdir()
            (metadata / "tk-8.6.13-noxft_test.json").write_text('{"name":"tk","version":"8.6.13","build":"noxft_test"}')
            root = MagicMock()
            values = {("tk::pkgconfig", "get", "fontsystem"): "xft", ("tk", "windowingsystem"): "x11",
                      ("package", "present", "Tk"): "8.6.13", ("tk", "scaling"): 1.333}
            root.tk.call.side_effect = lambda *args: values[args]
            with patch("tkinter.font.Font") as font, patch.object(cli_gui.sys, "prefix", temporary):
                font.return_value.actual.return_value = {"family": "Sans", "size": 12}
                result = cli_gui.get_font_info(root)
            self.assertEqual(result["font_renderer"], "xft")
            self.assertEqual(result["conda_tk_builds"][0]["build"], "noxft_test")
            font.return_value.configure.assert_not_called()

    def test_font_info_without_display_gives_actionable_error(self):
        with patch.dict(os.environ), patch("tkinter.Tk", side_effect=gui.tk.TclError("no display")), \
                contextlib.redirect_stderr(io.StringIO()) as output:
            self.assertEqual(cli_gui.main(["--font-info"]), 1)
        self.assertIn("display is required", output.getvalue())

    def test_launcher_delegates_without_changing_interpreter(self):
        with patch.object(gui, "main", return_value=0) as launch:
            self.assertEqual(cli_gui.main([]), 0)
        launch.assert_called_once_with()

    def test_python39_font_api(self):
        # Simulate an environment where nametofont cannot accept root=.
        style = MagicMock()
        style.element_names.return_value = ()
        parent = MagicMock()
        with patch.object(gui.ttk, "Style", return_value=style):
            with patch.object(gui.tkfont, "nametofont", side_effect=TypeError("unexpected keyword root")) as old:
                with patch.object(gui.tkfont, "Font") as font:
                    gui.configure_section_tabs(parent)
        old.assert_not_called()
        font.assert_called_once_with(root=parent, name="TkDefaultFont", exists=True)
        font.return_value.copy.return_value.configure.assert_called_once_with(weight="bold", size=13)

    def test_nonexistent_ini_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(ValueError):
                atlas_runner.build_current_python_command(Path(folder) / "missing.ini")

    def test_no_display_error_is_actionable(self):
        with patch.object(gui, "App", side_effect=gui.tk.TclError("no display")):
            with contextlib.redirect_stderr(io.StringIO()) as error:
                self.assertEqual(gui.main(), 1)
        self.assertIn("python -m tkinter", error.getvalue())

    def test_supplied_gui_regressions(self):
        directory = ROOT / "src/atlas_actris/GUI"
        result = subprocess.run(
            [sys.executable, "-m", "unittest", "discover", "-s", str(directory), "-v"],
            cwd=ROOT, capture_output=True, text=True, timeout=180,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    @unittest.skipUnless(os.environ.get("ATLAS_RUN_GUI_TESTS") == "1",
                         "Desktop tests are opt-in: atlas-test --gui --test gui_compatibility")
    def test_native_window_initializes_when_display_available(self):
        try:
            probe = gui.tk.Tk()
        except gui.tk.TclError as error:
            self.skipTest(f"No graphical display/Tk available: {error}")
        probe.destroy()
        app = gui.App(startup_prompt=False, remember_directory=False)
        try:
            app.withdraw()
            app.update_idletasks()
            self.assertEqual(len(app.editors), 3)
            self.assertFalse(hasattr(app.runner, "environment_combo"))
            self.assertTrue(all(not editor.include_descriptions.get() for editor in app.editors))
            init_editor = next(editor for editor in app.editors if editor.filename == "call_atlas.ini")
            titles = [init_editor.notebook.tab(tab, "text") for tab in init_editor.notebook.tabs()]
            self.assertNotIn("configuration", titles)
            self.assertEqual(gui.tkfont.Font(root=app, name="TkDefaultFont", exists=True).actual("size"), 12)
            variables = dict(init_editor.variables)
            with tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / "call_atlas.ini"
                path.write_text("[general_options]\nprocess_dedicated_dark = True\n")
                model = gui.EditorModel(init_editor.model.schema)
                model.load(path)
                init_editor.model = model
                init_editor.refresh_values()
                identity = ("general_options", "process_dedicated_dark")
                self.assertIs(init_editor.variables[identity], variables[identity])
                self.assertEqual(init_editor.variables[identity].get(), "True")
                self.assertFalse(init_editor.dirty())
        finally:
            app.destroy()

    @unittest.skipUnless(os.environ.get("ATLAS_RUN_GUI_TESTS") == "1",
                         "Desktop tests are opt-in: atlas-test --gui --test gui_compatibility")
    def test_native_folder_list_and_wheel_preserve_values_when_display_available(self):
        try:
            probe = gui.tk.Tk()
        except gui.tk.TclError as error:
            self.skipTest(f"No graphical display/Tk available: {error}")
        probe.destroy()
        app = gui.App(startup_prompt=False, remember_directory=False)
        dialog = None
        try:
            editor = next(e for e in app.editors if e.filename == "call_atlas.ini")
            app.notebook.select(editor)
            tab = next(t for t in editor.notebook.tabs() if editor.notebook.tab(t, "text") == "general_options")
            editor.notebook.select(tab)
            app.update()
            identity = ("general_options", "process_dedicated_dark")
            editor.variables[identity].set("True")
            row = next(row for field, row in editor.rows if field.identity == identity)
            combo = next(widget for cell in row.winfo_children() for widget in cell.winfo_children()
                         if isinstance(widget, gui.ttk.Combobox))
            page = app.nametowidget(tab)
            page.canvas.yview_moveto(0.2)
            app.update_idletasks()
            before = page.canvas.yview()[0]
            combo.event_generate("<MouseWheel>", delta=-120)
            app.update_idletasks()
            self.assertEqual(editor.variables[identity].get(), "True")
            self.assertGreater(page.canvas.yview()[0], before)
            from atlas_actris.GUI.directory_browser import DirectoryDialog
            with tempfile.TemporaryDirectory() as temporary:
                folder = Path(temporary).resolve()
                (folder / "call_atlas.ini").touch()
                (folder / ".hidden.ini").touch()
                (folder / "child").mkdir()
                dialog = DirectoryDialog(app, "Test folder browser", folder)
                deadline = time.monotonic() + 5
                while dialog.busy and time.monotonic() < deadline:
                    app.update()
                    time.sleep(0.01)
                self.assertFalse(dialog.busy)
                self.assertNotIn(".hidden.ini", {path.name for path, _ in dialog.rows.values()})
                dialog.hidden_toggle.invoke()
                deadline = time.monotonic() + 5
                while dialog.busy and time.monotonic() < deadline:
                    app.update()
                    time.sleep(0.01)
                self.assertFalse(dialog.busy)
                self.assertIn(".hidden.ini", {path.name for path, _ in dialog.rows.values()})
                self.assertIn("headings", str(dialog.table.cget("show")))
                file_row = next(key for key, (path, _) in dialog.rows.items() if path.name == "call_atlas.ini")
                self.assertEqual(dialog.table.item(file_row, "values")[1], "INI file")
                dialog.table.selection_set(file_row)
                dialog.accept()
                self.assertEqual(dialog.result, str(folder))
        finally:
            if dialog is not None and dialog.winfo_exists():
                dialog.destroy()
            app.destroy()

if __name__ == "__main__":
    unittest.main()
