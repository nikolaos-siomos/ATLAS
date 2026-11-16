"""Compact saves, whitespace, independent recommendations and folder startup."""
from pathlib import Path
import tempfile
from types import SimpleNamespace, MethodType
import queue
import unittest
from unittest.mock import Mock, patch

from editor_model import EditorModel, load_working_folder, normalize_value, validate_value
from schema_loader import load_schemas, RECOMMENDED_FIELDS
from atlas_schema_gui import App, SchemaTab, path_kind, configure_fonts

BASE = Path(__file__).resolve().parent


class UsabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schemas = load_schemas(BASE.parent)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas ini folder ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()

    def model(self, filename="call_atlas.ini"):
        return EditorModel(self.schemas[filename])

    def field(self, filename, section, key):
        return next(f for f in self.schemas[filename].fields if f.identity == (section, key))

    def test_bare_save_omits_blank_entries_and_empty_sections_keeps_explicit_defaults(self):
        model = self.model()
        model.values[("general_options", "process_dedicated_dark")] = " False "
        model.values[("general_options", "vertical_scale")] = "height_agl"
        model.values[("explicit_paths", "output_folder")] = "  "
        path = self.folder / "call_atlas.ini"
        model.save(path)
        self.assertEqual(path.read_text(), "[general_options]\nprocess_dedicated_dark = False\nvertical_scale = height_agl\n")
        self.assertFalse(model.dirty())

    def test_descriptions_can_be_enabled_and_removed_without_losing_values(self):
        model = self.model()
        identity = ("parsing_options", "cloudnet_station_name")
        model.values[identity] = " mindelo "
        path = self.folder / "call_atlas.ini"
        field = self.field("call_atlas.ini", *identity)
        model.save(path, include_descriptions=True)
        self.assertIn("# " + field.flavor["description"].splitlines()[0], path.read_text())
        model.save(path)
        self.assertNotIn("#", path.read_text())
        self.assertEqual(model.values[identity], "mindelo")

    def test_gui_save_as_uses_description_switch_and_normalizes_form(self):
        model = self.model()
        identity = ("parsing_options", "cloudnet_station_name")
        text = [" mindelo  "]
        model.values[identity] = text[0]
        def set_value(value):
            text[0] = value
            model.values[identity] = value
        tab = SimpleNamespace(model=model, filename="call_atlas.ini",
                              variables={identity: SimpleNamespace(get=lambda: text[0], set=set_value)},
                              browse_directory=lambda: self.folder,
                              include_descriptions=Mock(get=lambda: True),
                              update_status=Mock(), notify_path_changed=Mock(), status=Mock())
        tab.normalize_field = lambda field: SchemaTab.normalize_field(tab, field)
        tab.normalize_fields = lambda: SchemaTab.normalize_fields(tab)
        path = self.folder / "saved.ini"
        with patch("atlas_schema_gui.ask_save_file", return_value=str(path)) as ask:
            self.assertTrue(SchemaTab.save(tab, save_as=True))
        self.assertEqual(ask.call_args.kwargs["initialdir"], str(self.folder))
        self.assertEqual(text[0], "mindelo")
        self.assertIn("# ", path.read_text())
        self.assertIn("cloudnet_station_name = mindelo\n", path.read_text())
        tab.notify_path_changed.assert_called_once()

    def test_new_empty_drafts_have_no_empty_sections(self):
        for filename in self.schemas:
            model = self.model(filename)
            self.assertEqual(model.render(), b"")
            self.assertEqual(model.render(include_descriptions=True), b"")

    def test_trim_preserves_internal_spaces_and_list_positions(self):
        path = self.field("call_atlas.ini", "explicit_paths", "output_folder")
        self.assertEqual(normalize_value(path, "  /data/My station/output  "), "/data/My station/output")
        dead_time = self.field("config_file.ini", "Channels", "dead_time")
        self.assertEqual(normalize_value(dead_time, " 3.7 ,  , 4.2 , "), "3.7,,4.2,")
        self.assertEqual(validate_value(dead_time, " 3.7 ,  , 4.2 , ", "config_file.ini"), [])
        model = self.model("config_file.ini")
        model.values[dead_time.identity] = " 0 , , 4.2 , "
        self.assertIn(b"dead_time = 0,,4.2,\n", model.render())

    def test_line_breaks_still_rejected_before_overwriting(self):
        model = self.model()
        path = self.folder / "call_atlas.ini"; model.save(path)
        model.values[("explicit_paths", "output_folder")] = "bad\npath"
        with self.assertRaises(ValueError): model.save(path)
        self.assertEqual(path.read_bytes(), b"")

    def test_hidden_legacy_configuration_values_are_retained(self):
        path = self.folder / "call_atlas.ini"
        path.write_text("[configuration]\nscc_configuration_id = 713\nscc_compatible_format = False\n")
        model = self.model(); model.load(path); model.save(path)
        self.assertIn("scc_configuration_id = 713", path.read_text())
        self.assertIn("scc_compatible_format = False", path.read_text())

    def test_whitespace_from_existing_files_is_normalized_on_save(self):
        path = self.folder / "config_file.ini"
        path.write_text("# existing comments\n[System]\nstation_id = abc   \nstation_name = My station   \n[Channels]\ndead_time = 0 ,  , 4.2   \n")
        model = self.model("config_file.ini"); model.load(path); model.save(path)
        self.assertIn("station_name = My station\n", path.read_text())
        self.assertIn("dead_time = 0,,4.2\n", path.read_text())
        self.assertTrue(all(line == line.rstrip() for line in path.read_text().splitlines()))

    def test_independent_recommended_sets_do_not_change_parser_categories(self):
        for filename, expected in RECOMMENDED_FIELDS.items():
            fields = self.schemas[filename].fields
            self.assertEqual({f.identity for f in fields if f.category == "recommended"}, expected)
            self.assertTrue(all(f.meta.get("category", "optional") == "optional" for f in fields))
        self.assertEqual(len(self.schemas["config_file.ini"].fields), 47)

    def test_output_browser_does_not_require_existing_output_directory(self):
        field = self.field("call_atlas.ini", "explicit_paths", "output_folder")
        self.assertEqual(path_kind(field), "dir")
        self.assertNotIn("check_path", field.meta)
        variable = Mock()
        tab = SimpleNamespace(variables={field.identity: variable}, model=SimpleNamespace(path=None))
        target = str(self.folder / "new output")
        with patch("atlas_schema_gui.ask_directory", return_value=target) as ask:
            SchemaTab.browse_value(tab, field)
        self.assertFalse(ask.call_args.kwargs["mustexist"])
        variable.set.assert_called_once_with(target)

    def test_folder_loads_existing_files_byte_for_byte_and_creates_only_missing(self):
        path = self.folder / "call_atlas.ini"
        original = b"# untouched at startup\n[explicit_paths]\nparent_folder = /data/source\n"
        path.write_bytes(original)
        models, actions = load_working_folder(self.folder, self.schemas)
        self.assertEqual(path.read_bytes(), original)
        self.assertIn("Loaded call_atlas.ini", actions)
        self.assertIn("Prepared empty config_file.ini draft", actions)
        self.assertIn("Created settings_file.ini", actions)
        self.assertEqual(models["call_atlas.ini"].values[("explicit_paths", "parent_folder")], "/data/source")
        self.assertFalse((self.folder / "config_file.ini").exists())
        self.assertIsNone(models["config_file.ini"].path)
        self.assertEqual((self.folder / "settings_file.ini").read_bytes(), b"")
        _, actions = load_working_folder(self.folder, self.schemas)
        self.assertIn("Prepared empty config_file.ini draft", actions)

    def test_bad_existing_file_is_reported_before_creating_missing_files(self):
        (self.folder / "config_file.ini").write_text("[System]\nstation_id = a\nstation_id = b\n")
        with self.assertRaisesRegex(ValueError, "Duplicate key"):
            load_working_folder(self.folder, self.schemas)
        self.assertFalse((self.folder / "call_atlas.ini").exists())
        self.assertFalse((self.folder / "settings_file.ini").exists())

    def test_cancelled_folder_selection_has_no_effect(self):
        app = SimpleNamespace(runner=Mock(), working_folder=None, editors=[], _loading=False)
        app.runner.running.return_value = False
        with patch("atlas_schema_gui.ask_directory", return_value=""), \
                patch("atlas_schema_gui.load_working_folder") as load:
            App.select_working_folder(app)
        load.assert_not_called()
        self.assertIsNone(app.working_folder)

    def test_folder_selection_connects_models_and_run_target(self):
        editors = []
        runner = Mock()
        runner.running.return_value = False
        for filename in self.schemas:
            editor = SimpleNamespace(filename=filename, may_discard=Mock(return_value=True), build=Mock(), refresh_values=Mock())
            editor.notify_path_changed = lambda e=editor: runner.set_call_atlas_path(e.model.path) if e.filename == "call_atlas.ini" else None
            editors.append(editor)
        callbacks = []
        app = SimpleNamespace(runner=runner, working_folder=None, schemas=self.schemas,
                              editors=editors, folder_display=Mock(), notebook=Mock(), _loading=False,
                              after=lambda delay, callback: callbacks.append(callback))
        for name in ("begin_folder_load", "poll_folder_load", "prepare_next_editor", "end_folder_load"):
            setattr(app, name, MethodType(getattr(App, name), app))
        def background(function):
            result = queue.Queue(); result.put((True, function())); return result
        with patch("atlas_schema_gui.ask_directory", return_value=str(self.folder)), \
                patch("atlas_schema_gui.messagebox.showinfo") as message, \
                patch("atlas_schema_gui.LoadingDialog") as loading, \
                patch("atlas_schema_gui.background_call", side_effect=background):
            App.select_working_folder(app)
            self.assertTrue(app._loading)
            self.assertIsNone(app.working_folder)
            while callbacks:
                callbacks.pop(0)()
            loading.return_value.close.assert_called_once()
            self.assertFalse(app._loading)
        self.assertEqual(app.working_folder, self.folder)
        for editor in editors:
            editor.build.assert_not_called()
            editor.refresh_values.assert_called_once()
        runner.set_call_atlas_path.assert_called_once_with(self.folder / "call_atlas.ini")
        self.assertIn("Created call_atlas.ini", message.call_args.args[1])

    def test_font_configuration_uses_named_platform_fonts_and_point_sizes(self):
        parent = Mock()
        with patch("atlas_schema_gui.tkfont.Font") as font, patch("atlas_schema_gui.ttk.Style"):
            configure_fonts(parent)
        self.assertTrue(all(call.kwargs.get("exists") for call in font.call_args_list))
        self.assertTrue(all(call.kwargs["name"].startswith("Tk") for call in font.call_args_list))
        self.assertIn(unittest.mock.call(size=12), font.return_value.configure.call_args_list)


if __name__ == "__main__":
    unittest.main()
