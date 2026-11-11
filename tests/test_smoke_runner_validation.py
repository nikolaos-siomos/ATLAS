"""Verify output checks run before smoke-test cleanup, without running lidar I/O."""
import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from atlas_actris.testing.smoke_test import run_smoke_test, DEFAULT_CASE_NAME, DEFAULT_INI_NAME


class SmokeRunnerValidationTests(unittest.TestCase):
    def exercise(self, plot_count):
        with tempfile.TemporaryDirectory() as directory:
            pack=Path(directory)
            (pack/DEFAULT_INI_NAME).write_text('[test]\n')
            def simulated_atlas(*args,**kwargs):
                root=pack/'analysis'/DEFAULT_CASE_NAME
                for folder,count,suffix in (('plots',plot_count,'.png'),('reports',1,'.html'),
                                             ('reports',1,'.docx'),('ascii',5,'.txt')):
                    destination=root/folder
                    destination.mkdir(parents=True,exist_ok=True)
                    for i in range(count):
                        (destination/f'output_{i}{suffix}').touch()
                return SimpleNamespace(returncode=0)
            with patch('atlas_actris.testing.smoke_test.subprocess.run',side_effect=simulated_atlas):
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    with contextlib.redirect_stderr(io.StringIO()) as errors:
                        status=run_smoke_test(testing_pack=pack)
            self.assertFalse((pack/'analysis').exists())
            return status,output.getvalue(),errors.getvalue()

    def test_success_checks_outputs_then_removes_them(self):
        status,output,errors=self.exercise(10)
        self.assertEqual(status,0)
        self.assertIn('completed successfully',output)
        self.assertIn('Removing generated analysis folder',output)
        self.assertEqual(errors,'')

    def test_missing_outputs_fail_even_if_atlas_returns_zero(self):
        status,output,errors=self.exercise(0)
        self.assertEqual(status,1)
        self.assertIn('Too few plots',errors)


if __name__=='__main__':
    unittest.main()
