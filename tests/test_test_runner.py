"""CLI selection and process isolation checks; no lidar I/O or network calls."""
import contextlib
import io
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atlas_actris.testing import test_runner


class TestRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = Path(self.temp.name)
        (self.repo / "src/atlas_actris").mkdir(parents=True)
        (self.repo / "tests").mkdir()
        (self.repo / "pyproject.toml").write_text('[project]\nname="atlas_actris"\n')
        (self.repo / "pytest.ini").write_text('[pytest]\n')
        (self.repo / "tests/test_alpha.py").write_text('def test_one(): pass\n')
        (self.repo / "tests/test_beta.py").write_text('def test_two(): pass\n')
        (self.repo / "tests/__run_manual__.py").write_text('')

    def invoke(self, arguments=(), returncode=0):
        with patch.object(test_runner.importlib.util, "find_spec", return_value=object()):
            with patch.object(test_runner.subprocess, "run", return_value=SimpleNamespace(returncode=returncode)) as run:
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    status = test_runner.main(["--repo", str(self.repo), *arguments])
        return status, run

    def test_default_runs_whole_suite_with_current_interpreter(self):
        status, run = self.invoke()
        self.assertEqual(status, 0)
        args, kwargs = run.call_args
        self.assertEqual(args[0], [sys.executable, "-m", "pytest", "-c", str(self.repo / "pytest.ini"), "tests"])
        self.assertEqual(kwargs["cwd"], self.repo)
        self.assertEqual(kwargs["env"]["PYTHONPATH"].split(os.pathsep)[:2],
                         [str(self.repo / "src"), str(self.repo / "src/atlas_actris")])

    def test_module_and_individual_selection(self):
        status, run = self.invoke(["-t", "alpha", "--test", "tests/test_beta.py::test_two", "--", "-x", "-vv"])
        self.assertEqual(status, 0)
        self.assertEqual(run.call_args.args[0][-4:], ["tests/test_alpha.py", "tests/test_beta.py::test_two", "-x", "-vv"])

    def test_list_without_pytest_or_process_launch(self):
        output = io.StringIO()
        with patch.object(test_runner.importlib.util, "find_spec", return_value=None):
            with patch.object(test_runner.subprocess, "run") as run:
                with contextlib.redirect_stdout(output):
                    status = test_runner.main(["--repo", str(self.repo), "--list"])
        self.assertEqual(status, 0)
        self.assertIn("tests/test_alpha.py", output.getvalue())
        self.assertNotIn("__run_manual__", output.getvalue())
        run.assert_not_called()

    def test_failure_exit_status_is_preserved(self):
        status, _ = self.invoke(returncode=1)
        self.assertEqual(status, 1)

    def test_typo_cannot_silently_run_whole_suite(self):
        status, run = self.invoke(["-t", "missing"])
        self.assertEqual(status, 2)
        run.assert_not_called()

    def test_missing_pytest_is_actionable(self):
        errors = io.StringIO()
        with patch.object(test_runner.importlib.util, "find_spec", return_value=None):
            with contextlib.redirect_stderr(errors):
                status = test_runner.main(["--repo", str(self.repo)])
        self.assertEqual(status, 2)
        self.assertIn('python -m pip install -e ".[dev]"', errors.getvalue())

    def test_missing_checkout(self):
        with self.assertRaises(ValueError):
            test_runner.resolve_repository(str(self.repo / "missing"))

    def test_detect_checkout_from_nested_cwd(self):
        with patch.object(test_runner.Path, "cwd", return_value=self.repo / "tests"):
            self.assertEqual(test_runner.resolve_repository(), self.repo)

    def test_ambiguous_short_name_requires_path(self):
        sub = self.repo / "tests/nested"
        sub.mkdir()
        (sub / "test_alpha.py").touch()
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            test_runner.resolve_test(self.repo, "alpha")
        self.assertEqual(test_runner.resolve_test(self.repo, "tests/nested/test_alpha.py"),
                         "tests/nested/test_alpha.py")

    def test_source_and_environment_are_inherited_without_parent_mutation(self):
        with patch.dict(os.environ, {"PYTHONPATH": "custom-path", "MPLBACKEND": "custom-backend"}):
            _, run = self.invoke()
            env = run.call_args.kwargs["env"]
            self.assertEqual(env["MPLBACKEND"], "custom-backend")
            self.assertTrue(env["PYTHONPATH"].endswith(os.pathsep + "custom-path"))
            self.assertEqual(os.environ["PYTHONPATH"], "custom-path")


if __name__ == "__main__":
    unittest.main()
