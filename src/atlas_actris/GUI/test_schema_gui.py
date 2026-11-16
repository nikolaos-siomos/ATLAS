import ast
import codecs
from collections import Counter
from dataclasses import replace
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from schema_loader import load_schemas, SourceError, SourceReader, Unresolved, value_syntax
from editor_model import EditorModel, validate_value
from atlas_schema_gui import SchemaTab

BASE = Path(__file__).resolve().parent
SOURCE = BASE.parents[2]  # ATLAS project root, above src/atlas_actris/GUI.


def copy_source_fixture(project):
    # Copy only the schemas/flavor files under test, avoiding the repository's
    # datasets and recursion when TMPDIR happens to be inside the checkout.
    for subdir, names in (
        ('utils', ('parse_init_file.py', 'parse_config_file.py', 'parse_settings_file.py')),
        ('templates', ('init_template_flavor.py', 'config_template_flavor.py', 'settings_template_flavor.py')),
    ):
        destination = project / 'src/atlas_actris' / subdir
        destination.mkdir(parents=True, exist_ok=True)
        for name in names:
            shutil.copyfile(SOURCE / 'src/atlas_actris' / subdir / name, destination / name)


class SchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schemas = load_schemas(SOURCE)

    def test_source_fields_and_categories(self):
        self.assertEqual([len(s.fields) for s in self.schemas.values()], [59, 47, 127])
        fields = self.schemas['config_file.ini'].fields
        self.assertEqual(Counter(f.category for f in fields), {'mandatory': 5, 'recommended': 22, 'optional': 20})
        self.assertEqual({f.key for f in fields if f.category == 'mandatory'},
                         {'station_id', 'recorder_channel_id', 'telescope_type', 'channel_type', 'channel_subtype'})
        for schema in self.schemas.values():
            for f in schema.fields:
                self.assertTrue(f.flavor.get('description'), f.identity)
                self.assertIn(f.flavor['description'], f.help_text())

    def test_filter_combinations(self):
        fields = self.schemas['config_file.ini'].fields
        for rec, opt, expected in ((True, True, 47), (False, True, 25), (True, False, 27), (False, False, 5)):
            self.assertEqual(sum(f.visible(rec, opt) for f in fields), expected)
        self.assertEqual(sum(f.visible(True, False) for f in self.schemas['call_atlas.ini'].fields), 5)
        self.assertFalse(any(f.visible(False, False) for f in self.schemas['call_atlas.ini'].fields))
        self.assertFalse(any(f.visible(False, False) for f in self.schemas['settings_file.ini'].fields))

    def test_dynamic_choices_are_explicit(self):
        unresolved = [f for f in self.schemas['call_atlas.ini'].fields if isinstance(f.meta.get('allowed'), Unresolved)]
        self.assertEqual({f.key for f in unresolved}, {'export_stages', 'view_signal_stages', 'view_mean_signal_stages'})
        for f in unresolved:
            self.assertIn('collect_stages()', f.help_text())
            self.assertEqual(validate_value(f, 'custom_stage', 'call_atlas.ini'), [])

    def test_static_loading_does_not_execute_source(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / 'source.py'
            marker = Path(temp) / 'executed'
            p.write_text('from pathlib import Path\nPath(' + repr(str(marker)) + ').touch()\nraise RuntimeError("Do not run")\nVALUE = {"a": 1}\n')
            self.assertEqual(SourceReader(p).get('VALUE'), {'a': 1})
            self.assertFalse(marker.exists())

    def test_reload_tracks_source_changes_and_new_fields(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / 'atlas'
            copy_source_fixture(project)
            p = project / 'src/atlas_actris/utils/parse_settings_file.py'
            text = p.read_text().replace('SCHEMA_QCK: Dict[str, Dict[str, Any]] = {',
                'SCHEMA_QCK: Dict[str, Dict[str, Any]] = {\n "new_switch": {"dtype": bool, "default": False, "is_list": False, "category": "recommended"},', 1)
            p.write_text(text)
            schema = load_schemas(project)['settings_file.ini']
            self.assertEqual(len(schema.fields), 128)
            new = schema.sections['quicklooks'][0]
            self.assertEqual(new.key, 'new_switch')
            self.assertTrue(new.visible(True, False))
            self.assertFalse(new.visible(False, True))

    def test_configuration_group_addition_without_flavor_entry(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / 'atlas'
            copy_source_fixture(project)
            p = project / 'src/atlas_actris/utils/parse_config_file.py'
            text = p.read_text().replace('SCHEMA: Dict[str, Dict[str, Any]] = {',
                'SCHEMA: Dict[str, Dict[str, Any]] = {\n "new_key": {"dtype": str, "default": None, "is_list": False, "category": "optional"},', 1)
            text = text.replace('SYSTEM_KEYS = {', 'SYSTEM_KEYS = {\n "new_key",', 1)
            p.write_text(text)
            schema = load_schemas(project)['config_file.ini']
            self.assertEqual(len(schema.fields), 48)
            self.assertEqual(schema.sections['System'][0].key, 'new_key')

    def test_scalar_and_list_constraints(self):
        def check(filename, section, key, value, good):
            f = next(f for f in self.schemas[filename].fields if f.identity == (section, key))
            self.assertEqual(not validate_value(f, value, filename), good, (key, value))
        check('call_atlas.ini', 'general_options', 'process_vldr', 'True', True)
        check('call_atlas.ini', 'general_options', 'process_vldr', 'true', False)
        check('call_atlas.ini', 'general_options', 'vertical_scale', 'height_agl', True)
        check('call_atlas.ini', 'general_options', 'vertical_scale', 'bad', False)
        check('config_file.ini', 'System', 'station_latitude', '91', False)
        check('config_file.ini', 'Channels', 'dead_time', '3.7, _', True)
        check('config_file.ini', 'Channels', 'bins', '16384.0', True)
        check('config_file.ini', 'Channels', 'bins', '16384.5', False)
        check('settings_file.ini', 'quicklooks', 'x_lims', '0, 14', True)
        check('settings_file.ini', 'quicklooks', 'x_lims', '0,14', True)
        check('settings_file.ini', 'quicklooks', 'x_lims', '0;14', True)
        check('call_atlas.ini', 'parsing_options', 'rsonde_column_index', '1, 2, 3', True)
        check('call_atlas.ini', 'parsing_options', 'rsonde_column_index', '1, 2', False)

    def test_numeric_limits_follow_field_metadata(self):
        field = next(f for f in self.schemas['config_file.ini'].fields if f.key == 'laser_repetition_rate')
        # A test-owned bound exercises validation without pinning the live schema.
        field = replace(field, meta=dict(field.meta, max=10000.0))
        self.assertEqual(validate_value(field, '10000', 'config_file.ini'), [])
        self.assertTrue(validate_value(field, '10001', 'config_file.ini'))

    def test_config_optional_positions_are_preserved(self):
        f = next(f for f in self.schemas['config_file.ini'].fields if f.key == 'dead_time')
        # A size constraint makes dropped positions observable to the validator.
        f = replace(f, meta=dict(f.meta, size=3))
        for raw in ('3.7,,4', '3.7;None;4', '3.7,_,4', ',3.7,'):
            self.assertEqual(validate_value(f, raw, 'config_file.ini', required=True), [], raw)
        self.assertTrue(validate_value(f, '3.7,4', 'config_file.ini'))
        model = EditorModel(self.schemas['config_file.ini'])
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'config.ini'
            raw = b'[Channels]\ndead_time = 3.7,,None,_\n'
            path.write_bytes(raw)
            model.load(path)
            self.assertEqual(model.render(), raw)
            model.save(path)
            self.assertEqual(model.values[('Channels', 'dead_time')], '3.7,,None,_')

    def test_required_config_reports_missing_positions(self):
        fields = {f.key: f for f in self.schemas['config_file.ini'].fields}
        for raw, positions in (('a,,b', '2'), ('_,a,None', '1, 3'), (',a,', '1, 3')):
            errors = validate_value(fields['recorder_channel_id'], raw, 'config_file.ini', required=True)
            self.assertIn('positions ' + positions + ' (1-based)', errors[0])
            self.assertEqual(validate_value(fields['recorder_channel_id'], raw, 'config_file.ini'), [])
        for raw in ('_', 'None', ''):
            self.assertTrue(validate_value(fields['station_id'], raw, 'config_file.ini', required=True))
            self.assertEqual(validate_value(fields['station_id'], raw, 'config_file.ini'), [])
        settings = next(f for f in self.schemas['settings_file.ini'].fields if f.identity == ('quicklooks', 'x_lims'))
        self.assertTrue(validate_value(settings, '0,None', 'settings_file.ini'))

    def test_legacy_list_syntax_is_read_from_source(self):
        settings = next(f for f in self.schemas['settings_file.ini'].fields if f.identity == ('quicklooks', 'x_lims'))
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / 'parser.py'
            source.write_text('def _split_list(raw):\n    return [x.strip() for x in s.replace(";", ", ").split(", ") if x.strip() != ""]\n')
            legacy = replace(settings, value_syntax=value_syntax(SourceReader(source)))
            self.assertTrue(validate_value(legacy, '0,14', 'settings_file.ini'))
            self.assertEqual(validate_value(legacy, '0, 14', 'settings_file.ini'), [])
            self.assertEqual(validate_value(legacy, '0;14', 'settings_file.ini'), [])
            source.write_text('blank_tokens = ["_"]\ndef _split_list(raw):\n    return [x.strip() for x in s.replace(";", ",").split(",") if x.strip() != ""]\n')
            legacy = replace(settings, value_syntax=value_syntax(SourceReader(source)))
            self.assertEqual(validate_value(legacy, '0,,14', 'config_file.ini'), [])
            self.assertEqual(validate_value(legacy, '0,_', 'config_file.ini'), [])
            self.assertTrue(validate_value(legacy, '0,None', 'config_file.ini'))

    def test_legacy_config_layout_order_is_retained(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / 'atlas'
            copy_source_fixture(project)
            p = project / 'src/atlas_actris/templates/config_template_flavor.py'
            p.write_text(p.read_text() + '\nCONFIG_TEMPLATE_SECTIONS = {"System": ["station_name", "station_id"]}\n')
            schema = load_schemas(project)['config_file.ini']
            self.assertEqual([f.key for f in schema.sections['System'][:2]], ['station_name', 'station_id'])
            self.assertEqual(len(schema.fields), 47)

    def test_unsupported_layout_and_list_syntax_fail_explicitly(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / 'atlas'
            copy_source_fixture(project)
            p = project / 'src/atlas_actris/templates/config_template_flavor.py'
            original = p.read_text()
            p.write_text(original + '\nCONFIG_TEMPLATE_SECTIONS = unknown_layout()\n')
            with self.assertRaisesRegex(SourceError, 'template layout'):
                load_schemas(project)
            p.write_text(original)
            p = project / 'src/atlas_actris/utils/parse_settings_file.py'
            p.write_text(p.read_text().replace('.split(",")', '.split("|")'))
            with self.assertRaisesRegex(SourceError, 'list syntax'):
                load_schemas(project)

    def test_required_validation_allows_saving_drafts(self):
        model = EditorModel(self.schemas['config_file.ini'])
        self.assertEqual(len(model.validate()), 5)
        self.assertEqual(model.validate(required=False), [])
        with tempfile.TemporaryDirectory() as temp:
            model.save(Path(temp) / 'config.ini')
            self.assertEqual(len(model.validate()), 5)
            self.assertFalse(model.dirty())

    def test_real_ini_compact_roundtrip_and_edits(self):
        from gui_common import IniDocument, Setting
        def entries(raw):
            document = IniDocument(raw, inline_comments=False)
            return {(section.name, item.key): item.value for section in document.sections
                    for item in section.items if isinstance(item, Setting)}
        for filename, schema in self.schemas.items():
            model = EditorModel(schema)
            model.load(BASE / 'examples' / filename)
            before = entries(model.render())
            self.assertFalse(any(line.startswith(b'#') for line in model.render().splitlines()))
            for f in schema.fields:
                original = model.values[f.identity]
                model.values[f.identity] = 'edited'
                after = entries(model.render())
                expected = dict(before)
                expected[f.identity] = 'edited'
                self.assertEqual(after, expected)
                model.values[f.identity] = original

    def test_gui_filters_keep_values_and_do_not_recreate_fields(self):
        model = EditorModel(self.schemas['config_file.ini'])
        class Widget:
            def __init__(self): self.hidden = False
            def grid(self): self.hidden = False
            def grid_remove(self): self.hidden = True
        rows = [(f, Widget()) for f in model.schema.fields]
        model.values[('Channels', 'dead_time')] = '3.7'
        model.values[('System', 'lidar_name')] = 'Example'
        fake = SimpleNamespace(rows=rows, page_notes=[], recommended=SimpleNamespace(get=lambda:False), optional=SimpleNamespace(get=lambda:False))
        SchemaTab.apply_filters(fake)
        self.assertEqual(sum(not w.hidden for f,w in rows), 5)
        self.assertEqual(model.values[('Channels', 'dead_time')], '3.7')
        self.assertIn(b'dead_time = 3.7', model.render())
        self.assertIn(b'lidar_name = Example', model.render())
        self.assertTrue(model.dirty())

    def test_missing_keys_new_sections_and_unknown_preservation(self):
        model = EditorModel(self.schemas['config_file.ini'])
        raw = codecs.BOM_UTF8 + b'# comment\r\n[System]\r\nstation_id = abc\r\nunknown = keep\r\n'
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'config.ini'
            path.write_bytes(raw)
            model.load(path)
            model.values[('System', 'station_name')] = 'Test station'
            model.values[('Channels', 'recorder_channel_id')] = '0532xcpr'
            model.save(path)
            self.assertFalse(path.read_bytes().startswith(codecs.BOM_UTF8))
            self.assertNotIn(b'# comment', path.read_bytes())
            self.assertIn(b'unknown = keep', path.read_bytes())
            self.assertEqual(model.values[('System', 'station_name')], 'Test station')
            self.assertEqual(model.values[('Channels', 'recorder_channel_id')], '0532xcpr')
            self.assertEqual(model.unknown, [('System','unknown')])
            self.assertFalse(model.dirty())

    def test_file_changed_externally_is_not_overwritten(self):
        model = EditorModel(self.schemas['call_atlas.ini'])
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'call.ini'
            path.write_bytes(b'[configuration]\nscc_configuration_id = 713\n')
            model.load(path)
            path.write_bytes(b'# external edit\n')
            with self.assertRaises(ValueError): model.save(path)
            self.assertEqual(path.read_bytes(), b'# external edit\n')

    def test_duplicate_rejected_without_replacing_loaded_document(self):
        model = EditorModel(self.schemas['call_atlas.ini'])
        model.load(BASE / 'examples/call_atlas.ini')
        previous = model.document
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'bad.ini'
            path.write_text('[configuration]\nscc_configuration_id=713\nscc_configuration_id=714\n')
            with self.assertRaises(ValueError): model.load(path)
            self.assertIs(model.document, previous)

    def test_python_39_through_313_syntax(self):
        for filename in ('schema_loader.py', 'editor_model.py', 'gui_common.py', 'atlas_schema_gui.py', 'atlas_runner.py'):
            for minor in range(9, 14):
                ast.parse((BASE / filename).read_text(), feature_version=(3, minor))


if __name__ == '__main__':
    unittest.main()
