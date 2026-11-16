"""Loading a selected init must show the same associated files as processing."""
import os
from pathlib import Path
import queue
import tempfile
from types import MethodType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from atlas_schema_gui import App, SchemaTab
from editor_model import EditorModel, load_call_atlas_file, load_working_folder
from schema_loader import load_schemas


class LinkedIniTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schemas = load_schemas(Path(__file__).resolve().parents[1])

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas linked inis ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.call = self.folder / "custom run.ini"
        self.config = self.folder / "config_file.ini"
        self.settings = self.folder / "settings_file.ini"
        self.config.write_text("[System]\nstation_id = abc\n", encoding="utf-8")

    def init(self, text=""):
        self.call.write_text(text, encoding="utf-8")

    def test_default_siblings_and_empty_settings_creation_preserve_existing_init(self):
        self.init("# existing init\n[explicit_paths]\nparent_folder = /unrelated/data\n")
        before = self.call.read_bytes()
        models, actions = load_call_atlas_file(self.call, self.schemas)
        self.assertEqual(models["call_atlas.ini"].path, self.call)
        self.assertEqual(models["config_file.ini"].path, self.config)
        self.assertEqual(models["settings_file.ini"].path, self.settings)
        self.assertEqual(self.settings.read_bytes(), b"")
        self.assertIn("Created settings_file.ini", actions)
        self.assertEqual(self.call.read_bytes(), before)

    def test_relative_and_absolute_overrides_win_over_standard_siblings(self):
        external = self.folder / "other files"
        external.mkdir()
        config = external / "station.ini"
        settings = external / "processing.ini"
        config.write_text("[System]\nstation_id = xyz\n", encoding="utf-8")
        settings.write_text("# keep settings unchanged\n", encoding="utf-8")
        self.init("[explicit_paths]\nparent_folder = elsewhere\n"
                  "atlas_configuration_file = other files/station.ini\n"
                  "atlas_settings_file = {}\n".format(settings))
        models, _ = load_call_atlas_file(self.call, self.schemas)
        self.assertEqual(models["config_file.ini"].values[("System", "station_id")], "xyz")
        self.assertEqual(models["config_file.ini"].path, config)
        self.assertEqual(models["settings_file.ini"].path, settings)
        self.assertEqual(settings.read_bytes(), b"# keep settings unchanged\n")
        self.assertFalse(self.settings.exists())

    def test_explicit_missing_settings_created_at_referenced_path(self):
        folder = self.folder / "settings"
        folder.mkdir()
        self.init("[explicit_paths]\natlas_settings_file = settings\\custom.ini\n")
        models, _ = load_call_atlas_file(self.call, self.schemas)
        target = folder / "custom.ini"
        self.assertEqual(models["settings_file.ini"].path, target)
        self.assertEqual(target.read_bytes(), b"")
        self.assertFalse(self.settings.exists())

    def test_missing_explicit_config_opens_empty_draft_without_using_existing_sibling(self):
        self.init("[explicit_paths]\natlas_configuration_file = missing.ini\n")
        before = self.config.read_bytes()
        models, _ = load_call_atlas_file(self.call, self.schemas)
        draft = models["config_file.ini"]
        self.assertIsNone(draft.path)
        self.assertEqual(draft.render(), b"")
        self.assertIn("missing.ini", draft.load_warning)
        self.assertIn("Correct atlas_configuration_file", draft.load_warning)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertTrue(self.settings.exists())
        self.assertFalse((self.folder / "missing.ini").exists())

    def test_run_and_working_folder_use_same_empty_configuration_fallback(self):
        self.config.unlink()
        self.init()
        models, _ = load_call_atlas_file(self.call, self.schemas)
        self.assertIsNone(models["config_file.ini"].path)
        self.assertFalse(self.config.exists())
        self.call.rename(self.folder / "call_atlas.ini")
        models, _ = load_working_folder(self.folder, self.schemas)
        self.assertIsNone(models["config_file.ini"].path)
        self.assertEqual(models["config_file.ini"].render(), b"")
        self.assertFalse(self.config.exists())

    def test_invalid_existing_settings_are_preserved_with_warning_and_empty_draft(self):
        self.config.unlink()
        self.init()
        before = b"[quicklooks]\na = 1\na = 2\n"
        self.settings.write_bytes(before)
        models, _ = load_call_atlas_file(self.call, self.schemas)
        self.assertIsNone(models["settings_file.ini"].path)
        self.assertIn("Duplicate key", models["settings_file.ini"].load_warning)
        self.assertEqual(self.settings.read_bytes(), before)
        self.assertFalse(self.config.exists())

    def test_missing_settings_parent_opens_draft_without_creating_folders(self):
        self.init("[explicit_paths]\natlas_settings_file = missing/settings.ini\n")
        models, _ = load_call_atlas_file(self.call, self.schemas)
        self.assertIn("destination folder does not exist", models["settings_file.ini"].load_warning)
        self.assertIsNone(models["settings_file.ini"].path)
        self.assertFalse((self.folder / "missing").exists())

    def test_aliased_roles_open_separate_empty_draft(self):
        self.init("[explicit_paths]\natlas_configuration_file = custom run.ini\n")
        before = self.call.read_bytes()
        models, _ = load_call_atlas_file(self.call, self.schemas)
        self.assertIn("different file", models["config_file.ini"].load_warning)
        self.assertIsNone(models["config_file.ini"].path)
        self.assertEqual(self.call.read_bytes(), before)

    def test_foreign_drive_path_opens_draft_with_local_save_suggestion(self):
        if os.name == "nt":
            self.skipTest("Windows accepts native drive paths")
        self.init("[explicit_paths]\natlas_configuration_file = C:\\station\\config.ini\n")
        models, _ = load_call_atlas_file(self.call, self.schemas)
        self.assertIn("Windows drive", models["config_file.ini"].load_warning)
        self.assertIsNone(models["config_file.ini"].path)
        self.assertEqual(models["config_file.ini"].suggested_save_path, self.config)

    def make_app(self):
        callbacks, selected = [], [None]
        runner = Mock()
        runner.running.return_value = False
        runner.set_call_atlas_path.side_effect = lambda path: selected.__setitem__(0, path)
        editors = []
        for filename, schema in self.schemas.items():
            editor = SimpleNamespace(filename=filename, model=EditorModel(schema),
                may_discard=Mock(return_value=True), refresh_values=Mock())
            editor.notify_path_changed = lambda e=editor: runner.set_call_atlas_path(e.model.path) if e.filename == "call_atlas.ini" else None
            editors.append(editor)
        app = SimpleNamespace(_loading=False, runner=runner, schemas=self.schemas,
            editors=editors, working_folder=None, folder_display=Mock(), notebook=Mock(),
            after=lambda delay, callback: callbacks.append(callback))
        for name in ("select_call_atlas", "begin_folder_load", "poll_folder_load",
                     "prepare_next_editor", "end_folder_load"):
            setattr(app, name, MethodType(getattr(App, name), app))
        return app, callbacks, selected

    @staticmethod
    def background(function):
        result = queue.Queue()
        try:
            result.put((True, function()))
        except Exception as error:
            result.put((False, error))
        return result

    def test_app_selection_updates_all_tabs_and_run_target_after_loading(self):
        self.init()
        app, callbacks, selected = self.make_app()
        with patch("atlas_schema_gui.LoadingDialog"), \
                patch("atlas_schema_gui.background_call", side_effect=self.background), \
                patch("atlas_schema_gui.messagebox.showinfo") as summary:
            app.select_call_atlas(self.call)
            self.assertTrue(app._loading)
            self.assertIsNone(selected[0])
            while callbacks:
                callbacks.pop(0)()
        self.assertFalse(app._loading)
        self.assertEqual(selected[0], self.call)
        self.assertEqual(app.working_folder, self.folder)
        for editor in app.editors:
            editor.refresh_values.assert_called_once()
        self.assertIn(str(self.call), summary.call_args.args[1])
        self.assertIn(str(self.settings), summary.call_args.args[1])

    def test_failed_load_preserves_tabs_and_run_target(self):
        self.init("[explicit_paths]\natlas_configuration_file = a\natlas_configuration_file = b\n")
        app, callbacks, selected = self.make_app()
        originals = [editor.model for editor in app.editors]
        selected[0] = self.folder / "previous.ini"
        with patch("atlas_schema_gui.LoadingDialog"), \
                patch("atlas_schema_gui.background_call", side_effect=self.background), \
                patch("atlas_schema_gui.messagebox.showerror") as error:
            app.select_call_atlas(self.call)
            while callbacks:
                callbacks.pop(0)()
        self.assertFalse(app._loading)
        self.assertEqual([editor.model for editor in app.editors], originals)
        self.assertEqual(selected[0], self.folder / "previous.ini")
        self.assertFalse(self.settings.exists())
        error.assert_called_once()

    def test_cancelled_unsaved_changes_do_not_load_or_create_files(self):
        self.init()
        app, callbacks, _ = self.make_app()
        app.editors[0].may_discard.return_value = False
        with patch("atlas_schema_gui.background_call") as load:
            app.select_call_atlas(self.call)
        load.assert_not_called()
        self.assertEqual(callbacks, [])
        self.assertFalse(self.settings.exists())

    def load_models_into_app(self):
        models, _ = load_call_atlas_file(self.call, self.schemas)
        app, callbacks, selected = self.make_app()
        for editor in app.editors:
            editor.model = models[editor.filename]
            editor.dirty = lambda e=editor: e.model.dirty()
            editor.variables = {identity: Mock() for identity in editor.model.values}
            editor.update_status = Mock()
        app.linked_draft_saved = MethodType(App.linked_draft_saved, app)
        app.ready_to_run = MethodType(App.ready_to_run, app)
        return app

    def test_unedited_empty_configuration_draft_blocks_run(self):
        self.init("[explicit_paths]\natlas_configuration_file = missing.ini\n")
        app = self.load_models_into_app()
        self.assertFalse(any(editor.dirty() for editor in app.editors))
        with patch("atlas_schema_gui.messagebox.showinfo") as message:
            self.assertFalse(app.ready_to_run())
        self.assertIn("config_file.ini", message.call_args.args[1])

    def test_saving_filled_fallback_elsewhere_links_init_and_requires_init_save(self):
        self.init("[explicit_paths]\natlas_configuration_file = missing.ini\n")
        before = self.call.read_bytes()
        app = self.load_models_into_app()
        init = next(e for e in app.editors if e.filename == "call_atlas.ini")
        config = next(e for e in app.editors if e.filename == "config_file.ini")
        config.model.values[("System", "station_id")] = "xyz"
        config.normalize_fields = config.model.normalize
        config.include_descriptions = Mock(get=lambda: False)
        config.status = Mock(get=lambda: "Saved")
        config.notify_path_changed = lambda: app.linked_draft_saved(config.filename, config.model.path)
        destination = self.folder / "new station.ini"
        with patch("atlas_schema_gui.ask_save_file", return_value=str(destination)), \
                patch("atlas_schema_gui.messagebox.showinfo"):
            self.assertTrue(SchemaTab.save(config))
            self.assertFalse(app.ready_to_run())
        identity = ("explicit_paths", "atlas_configuration_file")
        self.assertEqual(init.model.values[identity], str(destination))
        init.variables[identity].set.assert_called_once_with(str(destination))
        self.assertTrue(init.model.dirty())
        self.assertEqual(self.call.read_bytes(), before)
        init.model.save(self.call)
        self.assertTrue(app.ready_to_run())
        self.assertIn(str(destination), self.call.read_text())
        self.assertIn("station_id = xyz", destination.read_text())

    def test_saving_at_original_missing_path_keeps_init_unchanged(self):
        self.init("[explicit_paths]\natlas_configuration_file = missing.ini\n")
        before = self.call.read_bytes()
        app = self.load_models_into_app()
        config = next(e for e in app.editors if e.filename == "config_file.ini")
        destination = self.folder / "missing.ini"
        self.assertEqual(config.model.suggested_save_path, destination)
        config.model.values[("System", "station_id")] = "xyz"
        config.model.save(destination)
        with patch("atlas_schema_gui.messagebox.showinfo") as message:
            app.linked_draft_saved(config.filename, destination)
        message.assert_not_called()
        self.assertEqual(self.call.read_bytes(), before)
        self.assertTrue(app.ready_to_run())

    def test_config_load_failure_keeps_new_init_open_and_shows_recovery_warning(self):
        self.init("[explicit_paths]\natlas_configuration_file = missing.ini\n")
        app, callbacks, selected = self.make_app()
        with patch("atlas_schema_gui.LoadingDialog"), \
                patch("atlas_schema_gui.background_call", side_effect=self.background), \
                patch("atlas_schema_gui.messagebox.showwarning") as warning:
            app.select_call_atlas(self.call)
            while callbacks:
                callbacks.pop(0)()
        self.assertEqual(selected[0], self.call)
        config = next(e for e in app.editors if e.filename == "config_file.ini")
        self.assertIsNone(config.model.path)
        self.assertEqual(config.model.render(), b"")
        warning.assert_called_once()
        self.assertIn("Correct atlas_configuration_file", warning.call_args.args[1])
        self.assertIn("fill in the empty", warning.call_args.args[1])

    def test_invalid_configuration_contents_are_not_overwritten_or_hidden(self):
        self.init()
        before = b"[System]\nstation_id = a\nstation_id = b\n"
        self.config.write_bytes(before)
        models, _ = load_call_atlas_file(self.call, self.schemas)
        config = models["config_file.ini"]
        self.assertIsNone(config.path)
        self.assertEqual(config.render(), b"")
        self.assertIn("Duplicate key", config.load_warning)
        self.assertEqual(self.config.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
