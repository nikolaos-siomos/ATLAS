"""The atlas init/CLI inputs cannot enable internal SCC configuration export."""
import ast
import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "src/atlas_actris"))
from atlas_actris.utils import parse_init_file as parser


class InternalSccDisabledTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ini = self.root / "call_atlas.ini"

    def read(self, value=None, cmd=None):
        text = "[configuration]\nscc_configuration_id = 713\n"
        if value is not None:
            text += f"export_hoi_cfg = {value}\n"
        self.ini.write_text(text)
        with contextlib.redirect_stdout(io.StringIO()) as output:
            result = parser.read_ini_file(str(self.ini), cmd_args=cmd)
        return result, output.getvalue()

    def test_mode_is_internal_and_absent_from_schema(self):
        self.assertNotIn("export_hoi_cfg", parser.SCHEMA)
        self.assertNotIn("export_hoi_cfg", parser.INIT_KEY_TO_SECTION)
        result, _ = self.read()
        self.assertEqual(result["export_hoi_cfg"], "0")

    def test_legacy_values_are_ignored_without_conversion(self):
        for value in ("", "0", "1", "2", "unexpected-value"):
            with self.subTest(value=value):
                result, output = self.read(value)
                self.assertEqual(result["export_hoi_cfg"], "0")
                self.assertIn("export_hoi_cfg is disabled and ignored", output)

    def test_programmatic_cli_override_cannot_enable_export(self):
        result, _ = self.read("2", cmd={"export_hoi_cfg": "1"})
        self.assertEqual(result["export_hoi_cfg"], "0")

    def test_unrelated_unknown_fields_still_fail(self):
        self.ini.write_text("[configuration]\nexport_hoi_cfg = 1\nunknown_option = 1\n")
        with self.assertRaises(parser.ConfigError):
            parser.read_ini_file(str(self.ini))

    def test_disabled_field_in_wrong_section_still_fails(self):
        self.ini.write_text("[general_options]\nexport_hoi_cfg = 1\n")
        with self.assertRaises(parser.ConfigError):
            parser.read_ini_file(str(self.ini))

    def test_local_configuration_path_used_even_for_legacy_export_mode(self):
        result, _ = self.read("1")
        config = self.root / "config_file.ini"
        config.touch()
        parser._export_hoi_cfg_check(result)
        parser._scc_configuration_id_check(result)
        parser._resolve_default_configuration_file(result, str(self.ini))
        self.assertEqual(Path(result["atlas_configuration_file"]), config)

    def test_main_does_not_import_or_call_downloader(self):
        caller = ROOT / "src/atlas_actris/__call_atlas_interactive__.py"
        tree = ast.parse(caller.read_text())
        guards = [node for node in tree.body if isinstance(node, ast.If)
                  and any(isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                          and n.func.id == "export_scc_config" for n in ast.walk(node))]
        self.assertEqual(len(guards), 1)
        self.assertFalse(any(isinstance(node, ast.ImportFrom)
                             and node.module == "utils.get_scc_config" for node in tree.body))
        for legacy in (None, "1", "2"):
            result, _ = self.read(legacy)
            code = compile(ast.Module(body=guards, type_ignores=[]), str(caller), "exec")
            with patch("builtins.__import__", side_effect=AssertionError("Downloader must not be imported")):
                exec(code, {"caller_info": result})


if __name__ == "__main__":
    unittest.main()
