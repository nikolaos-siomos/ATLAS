"""Exercise output pipes/process groups without processing measurements."""
import json
import io
import os
from pathlib import Path
import queue
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from atlas_runner import ProcessRunner, build_current_python_command
from gui_common import RunTab

FAKE_PROCESSOR = '''import json, os, pathlib, subprocess, sys, time
path = pathlib.Path(sys.argv[-1])
mode = path.read_text()
print("ARGS=" + json.dumps(sys.argv[1:]), flush=True)
print("CWD=" + os.getcwd(), flush=True)
print("UNBUFFERED=" + os.environ.get("PYTHONUNBUFFERED", ""), flush=True)
if mode == "stream":
    os.write(1, b"partial \\xc3")
    os.write(1, b"\\xa4 text")
    while not path.with_suffix(".release").exists():
        time.sleep(0.02)
    print(" stdout end", flush=True)
    print("stderr message", file=sys.stderr, flush=True)
elif mode == "fail":
    print("ATLAS failed", file=sys.stderr, flush=True)
    sys.exit(7)
elif mode == "child":
    child = subprocess.Popen([sys.executable, "-u", "-c",
        "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); print('CHILD_READY', flush=True); time.sleep(60)"])
    child.wait()
'''


class RunnerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas process ")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.fake = self.folder / "fake processor.py"
        self.fake.write_text(FAKE_PROCESSOR)
        self.ini = self.folder / "call atlas.ini"
        self.ini.write_text("stream")
        self.command = [sys.executable, "-u", str(self.fake), "-i", str(self.ini)]
        self.engine = ProcessRunner()
        self.addCleanup(self.cleanup_process)

    def cleanup_process(self):
        if self.engine.busy:
            self.engine.stop()
            deadline = time.monotonic() + 8
            while time.monotonic() < deadline:
                try:
                    kind, value = self.engine.events.get(timeout=0.1)
                except queue.Empty:
                    continue
                if kind == "done":
                    self.engine.finish()
                    break

    def collect(self, until):
        text, events = "", []
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            try:
                kind, value = self.engine.events.get(timeout=0.1)
            except queue.Empty:
                continue
            events.append((kind, value))
            if kind == "text":
                text += value
            if until(text, events):
                return text, events
        self.fail("Timed out waiting for process output: " + text)

    def test_streams_before_exit_and_merges_stdout_stderr(self):
        command = self.command
        self.engine.start(command, self.folder)
        first, _ = self.collect(lambda text, events: "partial ä text" in text)
        self.assertIsNone(self.engine.process.poll())
        args = json.loads(first.split("ARGS=", 1)[1].splitlines()[0])
        self.assertEqual(args, ["-i", str(self.ini)])
        self.assertIn("CWD=" + str(self.folder), first)
        self.assertIn("UNBUFFERED=1", first)
        self.ini.with_suffix(".release").touch()
        rest, events = self.collect(lambda text, events: any(kind == "done" for kind, _ in events))
        self.assertIn("stdout end", rest)
        self.assertIn("stderr message", rest)
        self.assertIn(("done", 0), events)
        self.assertTrue(self.engine.finish())

    def test_nonzero_exit_and_restart(self):
        self.ini.write_text("fail")
        command = self.command
        for _ in range(2):
            self.engine.start(command, self.folder)
            with self.assertRaises(RuntimeError):
                self.engine.start(command, self.folder)
            text, events = self.collect(lambda text, events: any(kind == "done" for kind, _ in events))
            self.assertIn("ATLAS failed", text)
            self.assertIn(("done", 7), events)
            self.assertTrue(self.engine.finish())

    @unittest.skipIf(sys.platform == "win32", "POSIX process-group integration test")
    def test_stop_terminates_child_even_if_wrapper_exits_first(self):
        self.ini.write_text("child")
        self.engine.start(self.command, self.folder)
        self.collect(lambda text, events: "CHILD_READY" in text)
        self.engine.stop()
        _, events = self.collect(lambda text, events: any(kind == "done" for kind, _ in events))
        self.assertFalse(self.engine.stop_failed)
        self.assertTrue(self.engine.finish())
        self.assertFalse(self.engine.busy)
        self.assertNotEqual(next(value for kind, value in events if kind == "done"), 0)

    def test_missing_executable_leaves_runner_idle(self):
        with self.assertRaises(OSError):
            self.engine.start([str(self.folder / "missing")], self.folder)
        self.assertFalse(self.engine.busy)

    def test_gui_start_reads_current_path_and_honors_unsaved_edits(self):
        first = self.ini
        second = self.folder / "different.ini"
        second.write_text("fail")
        current = [first]
        engine = Mock(busy=False)
        tab = SimpleNamespace(engine=engine,
            get_call_atlas_path=lambda: current[0],
            before_start=Mock(return_value=True), append=Mock(),
            status=Mock(), start_button=Mock(), stop_button=Mock(),
            clean_cache=Mock(get=Mock(return_value=False)),
            stage_action=Mock(get=Mock(return_value="Default")),
            cache_checkbox=Mock(), stage_selector=Mock())
        RunTab.start(tab)
        self.assertEqual(engine.start.call_args.args[0][-1], str(first))
        self.assertEqual(engine.start.call_args.kwargs, {"clean_cache": False, "stage_action": "Default"})
        tab.cache_checkbox.configure.assert_called_with(state="disabled")
        tab.stage_selector.configure.assert_called_with(state="disabled")
        current[0] = second
        tab.clean_cache.get.return_value = True
        tab.stage_action.get.return_value = "Export stages"
        RunTab.start(tab)
        self.assertEqual(engine.start.call_args.args[0][-1], str(second))
        self.assertEqual(engine.start.call_args.kwargs, {"clean_cache": True, "stage_action": "Export stages"})
        tab.before_start.return_value = False
        engine.start.reset_mock()
        RunTab.start(tab)
        engine.start.assert_not_called()

    def test_child_receives_gui_choices_without_changing_parent_environment(self):
        keys = ("ATLAS_CLEAN_CACHE_ANSWER", "ATLAS_DELETE_EXPORTED_ANSWER", "ATLAS_EXPORT_STAGE_ANSWER")
        code = "import os,json; print(json.dumps([os.environ.get(k) for k in {}]), flush=True)".format(repr(keys))
        with patch.dict(os.environ, {key: "Y" for key in keys}):
            for clean, action, expected in ((False, "Default", ["N", "N", "N"]),
                                             (True, "Export stages", ["Y", "N", "Y"]),
                                             (False, "Clean stages", ["N", "Y", "N"])):
                self.engine.start([sys.executable, "-u", "-c", code], self.folder,
                                  clean_cache=clean, stage_action=action)
                text, _ = self.collect(lambda text, events: any(kind == "done" for kind, _ in events))
                self.assertEqual(json.loads(text), expected)
                self.assertTrue(self.engine.finish())
                self.assertEqual([os.environ[key] for key in keys], ["Y", "Y", "Y"])

    def test_old_runner_reports_incomplete_update_without_callback_traceback(self):
        def old_start(command, cwd):
            self.fail("An outdated runner must not start the processing job")
        tab = SimpleNamespace(engine=SimpleNamespace(busy=False, start=old_start),
            get_call_atlas_path=lambda: self.ini, before_start=None,
            clean_cache=Mock(get=Mock(return_value=False)),
            stage_action=Mock(get=Mock(return_value="Default")))
        with patch("gui_common.messagebox.showerror") as show:
            RunTab.start(tab)
        show.assert_called_once()
        self.assertIn("different ATLAS updates", show.call_args.args[1])

    def test_windows_uses_hidden_process_and_taskkill_tree(self):
        process = Mock(pid=12345, stdout=io.BytesIO(b"Windows output\n"))
        process.wait.return_value = 0
        process.poll.return_value = 0
        with patch("atlas_runner.sys.platform", "win32"), \
                patch("atlas_runner.subprocess.CREATE_NO_WINDOW", 0x08000000, create=True), \
                patch("atlas_runner.subprocess.CREATE_NEW_PROCESS_GROUP", 0x00000200, create=True), \
                patch("atlas_runner.subprocess.Popen", return_value=process) as launch:
            self.engine.start(build_current_python_command(self.ini), self.folder)
            options = launch.call_args.kwargs
            self.assertEqual(options["creationflags"], 0x08000200)
            self.assertNotIn("shell", options)
            text, _ = self.collect(lambda text, events: any(kind == "done" for kind, _ in events))
            self.assertIn("Windows output", text)
            self.engine.finish()
            with patch("atlas_runner.subprocess.run", return_value=Mock(returncode=0)) as kill:
                self.engine._stop_tree(process)
                self.assertEqual(kill.call_args.args[0][1:], ["/PID", "12345", "/T", "/F"])



if __name__ == "__main__":
    unittest.main()
