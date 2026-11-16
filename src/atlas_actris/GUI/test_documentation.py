"""Check PDF launching without opening external applications."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import gui_common


class DocumentationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas gui ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.pdf = self.folder / "GUI-doc.pdf"
        self.pdf.write_bytes(b"%PDF-1.4\n")
        source = patch.object(gui_common, "__file__", str(self.folder / "gui_common.py"))
        source.start()
        self.addCleanup(source.stop)
        self.parent = Mock()

    def test_missing_document_does_not_launch_viewer(self):
        self.pdf.unlink()
        with patch.object(gui_common.messagebox, "showerror") as error, \
                patch.object(gui_common.subprocess, "Popen") as launch, \
                patch.object(gui_common.os, "startfile", create=True) as startfile:
            gui_common.open_documentation(self.parent)
            launch.assert_not_called()
            startfile.assert_not_called()
            self.assertIn(str(self.folder), error.call_args.args[1])

    def test_windows_opens_document_beside_gui(self):
        with patch.object(gui_common.sys, "platform", "win32"), \
                patch.object(gui_common.os, "startfile", create=True) as startfile:
            gui_common.open_documentation(self.parent)
            startfile.assert_called_once_with(str(self.pdf))

    def test_linux_and_macos_launch_absolute_path_as_one_argument(self):
        for platform, command in (("linux", "xdg-open"), ("darwin", "open")):
            with self.subTest(platform=platform), \
                    patch.object(gui_common.sys, "platform", platform), \
                    patch.object(gui_common.subprocess, "Popen") as launch, \
                    patch.object(gui_common.messagebox, "showerror") as error:
                launch.return_value.poll.return_value = 0
                gui_common.open_documentation(self.parent)
                self.assertEqual(launch.call_args.args[0], [command, str(self.pdf)])
                self.assertFalse(launch.call_args.kwargs.get("shell", False))
                self.parent.after.call_args.args[1]()
                error.assert_not_called()

    def test_viewer_launch_errors_are_reported(self):
        with patch.object(gui_common.sys, "platform", "linux"), \
                patch.object(gui_common.subprocess, "Popen", side_effect=OSError("No viewer")), \
                patch.object(gui_common.messagebox, "showerror") as error:
            gui_common.open_documentation(self.parent)
            self.assertIn("No viewer", error.call_args.args[1])
        with patch.object(gui_common.sys, "platform", "linux"), \
                patch.object(gui_common.subprocess, "Popen") as launch, \
                patch.object(gui_common.messagebox, "showerror") as error:
            launch.return_value.poll.side_effect = [None, 3]
            gui_common.open_documentation(self.parent)
            self.parent.after.call_args.args[1]()
            error.assert_not_called()
            self.parent.after.call_args.args[1]()
            self.assertIn(str(self.pdf), error.call_args.args[1])


if __name__ == "__main__":
    unittest.main()
