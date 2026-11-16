"""Real subprocess checks for trusted imports and INI-relative data paths."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import atlas_runner


class InstalledLauncherTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas launcher ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.package = self.folder / "installed location" / "atlas_actris"
        self.gui = self.package / "GUI"
        self.gui.mkdir(parents=True)
        (self.package / "__init__.py").touch()
        shutil.copyfile(Path(atlas_runner.__file__).with_name("_run_atlas.py"),
                        self.gui / "_run_atlas.py")
        (self.package / "cli.py").write_text(
            "import json, os, sys\nfrom pathlib import Path\n"
            "def main():\n"
            " print(json.dumps(dict(package=str(Path(__file__).resolve()), "
            "python=sys.executable, cwd=os.getcwd(), args=sys.argv[1:], "
            "relative=Path('relative.txt').read_text(), isolated=sys.flags.isolated)))\n"
            " return 7\n", encoding="utf-8")
        self.data = self.folder / "received data"
        self.data.mkdir()
        self.ini = self.data / "call atlas.ini"
        self.ini.touch()
        (self.data / "relative.txt").write_text("relative data", encoding="utf-8")

    def command(self):
        with patch.object(atlas_runner, "__file__", str(self.gui / "atlas_runner.py")):
            return atlas_runner.build_current_python_command(self.ini)

    def test_ignores_local_package_modules_and_python_environment_overrides(self):
        marker = self.data / "untrusted_code_executed"
        payload = "from pathlib import Path\nPath({!r}).touch()\nraise RuntimeError('untrusted import')\n".format(str(marker))
        fake = self.data / "atlas_actris"
        fake.mkdir()
        (fake / "__init__.py").write_text(payload, encoding="utf-8")
        for name in ("sitecustomize.py", "usercustomize.py", "json.py", "pathlib.py"):
            (self.data / name).write_text(payload, encoding="utf-8")
        environment = os.environ.copy()
        environment.update(PYTHONPATH=str(self.data), PYTHONHOME=str(self.data / "not Python"))
        result = subprocess.run(self.command(), cwd=self.data, env=environment,
                                capture_output=True, text=True, encoding="utf-8", timeout=15)
        self.assertEqual(result.returncode, 7, result.stderr)
        self.assertFalse(marker.exists())
        details = json.loads(result.stdout)
        self.assertEqual(details["package"], str(self.package / "cli.py"))
        self.assertEqual(details["python"], sys.executable)
        self.assertEqual(details["cwd"], str(self.data))
        self.assertEqual(details["args"], ["-i", str(self.ini)])
        self.assertEqual(details["relative"], "relative data")
        self.assertEqual(details["isolated"], 1)

    def test_source_edits_are_seen_on_the_next_run(self):
        # The trusted package can be in an editable checkout, outside sys.prefix.
        (self.package / "cli.py").write_text("def main():\n print('updated source')\n return 0\n")
        result = subprocess.run(self.command(), cwd=self.data, capture_output=True,
                                text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "updated source")

    def test_missing_bootstrap_is_reported_before_starting(self):
        (self.gui / "_run_atlas.py").unlink()
        with self.assertRaisesRegex(FileNotFoundError, "Reinstall ATLAS"):
            self.command()

    def test_missing_atlas_never_falls_back_to_another_package(self):
        (self.package / "cli.py").unlink()
        result = subprocess.run(self.command(), cwd=self.data, capture_output=True,
                                text=True, timeout=15)
        self.assertEqual(result.returncode, 1)
        self.assertIn("Reinstall ATLAS", result.stderr)


if __name__ == "__main__":
    unittest.main()
