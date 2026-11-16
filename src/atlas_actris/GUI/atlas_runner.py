"""Launch this GUI's ATLAS installation and stream its output."""
import codecs
import os
from pathlib import Path
import queue
import signal
import subprocess
import sys
import threading


def build_current_python_command(ini_path):
    """Use this GUI's interpreter and trusted bootstrap, never the INI folder."""
    path = Path(ini_path).expanduser().resolve()
    if not path.is_file():
        raise ValueError("The selected INI file does not exist:\n{}".format(path))
    bootstrap = Path(__file__).resolve().with_name("_run_atlas.py")
    if not bootstrap.is_file():
        raise FileNotFoundError("The ATLAS GUI launcher is missing. Reinstall ATLAS in this environment.")
    # -I excludes the working directory, user site, and PYTHON* overrides.
    # Explicit flags retain streaming UTF-8 output even in isolated mode.
    return [sys.executable, "-I", "-u", "-X", "utf8", str(bootstrap), "-i", str(path)]


def window_options():
    return {"creationflags": subprocess.CREATE_NO_WINDOW} if sys.platform == "win32" else {}


def post_run_environment(clean_cache=False, stage_action="Default"):
    """Explicit per-run choices override inherited terminal prompt answers."""
    if stage_action not in ("Default", "Clean stages", "Export stages"):
        raise ValueError("Unknown stage action: {}".format(stage_action))
    return {
        "ATLAS_GUI_POST_RUN": "1",
        "ATLAS_CLEAN_CACHE_ANSWER": "Y" if clean_cache else "N",
        "ATLAS_CLEAN_VIEWER_ANSWER": "N",
        "ATLAS_DELETE_EXPORTED_ANSWER": "Y" if stage_action == "Clean stages" else "N",
        "ATLAS_EXPORT_STAGE_ANSWER": "Y" if stage_action == "Export stages" else "N",
    }


class ProcessRunner:
    """Only events cross into Tk; workers never access Tkinter widgets."""
    def __init__(self):
        self.events = queue.Queue(maxsize=2048)
        self.process = None
        self.busy = False
        self.stopping = False
        self.stop_failed = False
        self._stop_complete = threading.Event()
        self._stop_complete.set()

    def start(self, command, cwd, clean_cache=False, stage_action="Default"):
        if self.busy:
            raise RuntimeError("ATLAS is already running.")
        environment = os.environ.copy()
        environment.update(PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8", PYTHONUTF8="1", NO_COLOR="1")
        environment.update(post_run_environment(clean_cache, stage_action))
        options = window_options()
        if sys.platform == "win32":
            options["creationflags"] |= subprocess.CREATE_NEW_PROCESS_GROUP
        else:
            options["start_new_session"] = True
        process = subprocess.Popen(command, cwd=str(cwd), env=environment, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=0, **options)
        self.process = process
        self.busy, self.stopping, self.stop_failed = True, False, False
        self._stop_complete.set()
        threading.Thread(target=self._read_output, args=(process,), daemon=True).start()

    def _read_output(self, process):
        decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
        try:
            with process.stdout:
                while True:
                    chunk = process.stdout.read(4096)
                    if not chunk:
                        break
                    text = decoder.decode(chunk)
                    if text:
                        self.events.put(("text", text))
                tail = decoder.decode(b"", final=True)
                if tail:
                    self.events.put(("text", tail))
        except OSError as error:
            self.events.put(("text", "\nOutput error: {}\n".format(error)))
        code = process.wait()
        self._stop_complete.wait()
        self.events.put(("done", code))

    def finish(self):
        if not self._stop_complete.is_set():
            return False
        self.busy = False
        self.process = None
        return True

    def stop(self):
        if not self.busy or self.stopping:
            return
        self.stopping, self.stop_failed = True, False
        self._stop_complete.clear()
        threading.Thread(target=self._stop_tree, args=(self.process,), daemon=True).start()

    def _stop_tree(self, process):
        try:
            if sys.platform == "win32":
                taskkill = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32/taskkill.exe"
                result = subprocess.run([str(taskkill), "/PID", str(process.pid), "/T", "/F"],
                    stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    timeout=15, **window_options())
                if result.returncode and process.poll() is None:
                    raise RuntimeError(result.stdout.decode("utf-8", errors="replace"))
            else:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    return
                # The wrapper may exit before its children. Kill the group after a grace period.
                threading.Event().wait(2)
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
            self.stop_failed = True
            self.stopping = False
            self.events.put(("stop_error", str(error)))
        finally:
            self._stop_complete.set()
