"""Shared INI I/O, Tk widgets and documentation launcher."""
import codecs
import inspect


import os


from pathlib import Path


import re
import queue
import shlex


import stat


import subprocess


import sys


import tempfile


import tkinter as tk


from tkinter import messagebox, ttk, font as tkfont


from dataclasses import dataclass, field
if __package__:
    from .atlas_runner import ProcessRunner, build_current_python_command
    from .directory_browser import ask_open_file
else:
    from atlas_runner import ProcessRunner, build_current_python_command
    from directory_browser import ask_open_file


FILES = ("call_atlas.ini", "config_file.ini", "settings_file.ini")


SECTION = re.compile(r"^\s*\[([^\]\r\n]+)\]\s*([;#].*)?$")


OPTION = re.compile(r"^(\s*[^=:#;\s][^=:]*?\s*[=:][ \t]*)(.*)$")


@dataclass
class Setting:
    line: int
    key: str
    prefix: str
    value: str
    suffix: str
    newline: str
    indent: int
    multiline: bool = False


@dataclass
class Section:
    name: str
    items: list = field(default_factory=list)
    implicit: bool = False
    heading_comment: str = ""


def comment_text(line):
    """Remove only the comment marker and one optional space."""
    text = line.lstrip()[1:]
    return text[1:] if text.startswith(" ") else text


def value_parts(text, inline_comments):
    """Recognize inline comments outside quotes, following whitespace."""
    quote = None
    end = len(text)
    escaped = False
    if inline_comments:
        for pos, char in enumerate(text):
            if char in ";#" and quote is None and (pos == 0 or text[pos - 1].isspace()):
                end = pos
                break
            if char in "\"'" and not escaped:
                if quote == char:
                    quote = None
                elif quote is None:
                    quote = char
            escaped = (char == "\\" and not escaped)
    value = text[:end].rstrip(" \t")
    return value, text[len(value):]


class IniDocument:
    """Keep the original lines; replace only values when writing."""
    def __init__(self, raw, encoding="utf-8", inline_comments=True):
        self.raw = raw
        self.encoding = encoding
        self.inline_comments = inline_comments
        self.lines = raw.decode(encoding).splitlines(keepends=True)
        self.sections = []
        section = None
        previous = None
        for index, line in enumerate(self.lines):
            body = line.rstrip("\r\n")
            newline = line[len(body):]
            stripped = body.strip()
            indent = len(body) - len(body.lstrip())
            if not stripped or stripped.startswith((";", "#")):
                if section is None:
                    section = Section("No section", implicit=True)
                    self.sections.append(section)
                section.items.append(body)
                continue
            # More deeply indented content belongs to the previous value.
            if previous is not None and indent > previous.indent:
                previous.multiline = True
                section.items.append(body)
                continue
            match = SECTION.match(body)
            if match:
                section = Section(match.group(1), heading_comment=comment_text(match.group(2)) if match.group(2) else "")
                self.sections.append(section)
                previous = None
                continue
            match = OPTION.match(body)
            if not match:
                raise ValueError("Unsupported INI line {}: {}".format(index + 1, body))
            if section is None:
                section = Section("No section", implicit=True)
                self.sections.append(section)
            prefix, tail = match.groups()
            value, suffix = value_parts(tail, inline_comments)
            key = re.split(r"[=:]", prefix, maxsplit=1)[0].strip()
            previous = Setting(index, key, prefix, value, suffix, newline, indent)
            section.items.append(previous)
        if not self.sections:
            self.sections.append(Section("No section", implicit=True))

    @classmethod
    def load(cls, path, inline_comments=True):
        raw = Path(path).read_bytes()
        # Honor explicit BOMs; otherwise use UTF-8 without guessing other encodings.
        if raw.startswith((codecs.BOM_UTF32_LE, codecs.BOM_UTF32_BE)):
            raise ValueError("UTF-32 is not supported. Please save the file as UTF-8.")
        if raw.startswith(codecs.BOM_UTF8):
            encoding = "utf-8-sig"
        elif raw.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
            # Preserve the explicit byte order, including the original BOM.
            encoding = "utf-16-le" if raw.startswith(codecs.BOM_UTF16_LE) else "utf-16-be"
            obj = cls(raw[2:], encoding, inline_comments)
            obj.raw = raw
            obj.bom = raw[:2]
            return obj
        else:
            encoding = "utf-8"
        return cls(raw, encoding, inline_comments)

    def render(self, values):
        lines = self.lines.copy()
        for section in self.sections:
            for item in section.items:
                if not isinstance(item, Setting):
                    continue
                value = values.get(item.line, item.value)
                if value == item.value:
                    continue
                if item.multiline:
                    raise ValueError("Multiline values are read-only in this version.")
                if "\n" in value or "\r" in value:
                    raise ValueError("A value must not contain a line break.")
                # Prevent new values from being interpreted as comments or sections.
                if value != value.strip(" \t"):
                    raise ValueError("The value for '{}' has leading or trailing spaces. Use quotation marks if needed.".format(item.key))
                if value_parts(value, self.inline_comments)[0] != value:
                    raise ValueError("The value for '{}' contains a comment marker. Use quotation marks or reload with inline comment detection disabled.".format(item.key))
                suffix = item.suffix
                if value and suffix.startswith((";", "#")):
                    suffix = " " + suffix
                lines[item.line] = item.prefix + value + suffix + item.newline
        return getattr(self, "bom", b"") + "".join(lines).encode(self.encoding)


def atomic_write(path, data):
    """Write completely before replacing the file; preserve its permissions."""
    path = Path(path).resolve()
    descriptor, name = tempfile.mkstemp(prefix="." + path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        if path.exists():
            os.chmod(name, stat.S_IMODE(path.stat().st_mode))
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


class ScrollFrame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        canvas = tk.Canvas(self, highlightthickness=0)
        self.canvas = canvas
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.body = ttk.Frame(canvas, padding=12)
        window = canvas.create_window((0, 0), window=self.body, anchor="nw")
        self.body.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        # Scroll only within this page; do not use global bindings.
        def bind_scroll(widget):
            sequences = ["<MouseWheel>"]
            if self.tk.call("tk", "windowingsystem") == "x11":
                sequences += ["<Button-4>", "<Button-5>"]
            for sequence in sequences:
                widget.bind(sequence, self.scroll_wheel, add="+")
            for child in widget.winfo_children():
                bind_scroll(child)
        self.bind_scroll = lambda: bind_scroll(self)


    def scroll_wheel(self, event):
        delta, number = getattr(event, "delta", 0), getattr(event, "num", None)
        if delta or number in (4, 5):
            self.canvas.yview_scroll(-1 if delta > 0 or number == 4 else 1, "units")
        # Widget bindings run before ttk class bindings. Stop the latter from
        # changing a Combobox selection after the page has already scrolled.
        return "break"


class LoadingDialog(tk.Toplevel):
    """Visible busy state while models load and editor widgets are prepared."""
    def __init__(self, parent, message):
        super().__init__(parent)
        self.title("Loading INI files")
        self.transient(parent.winfo_toplevel())
        self.resizable(False, False)
        self.protocol("WM_DELETE_WINDOW", lambda: None)
        self.message = tk.StringVar(value=message)
        ttk.Label(self, textvariable=self.message, wraplength=550, padding=16).pack(fill="x")
        self.progress = ttk.Progressbar(self, mode="indeterminate", length=420)
        self.progress.pack(fill="x", padx=16, pady=(0, 16))
        self.configure(cursor="watch")
        self.progress.start(30)
        self.update_idletasks()
        if not self.winfo_viewable():
            self.wait_visibility()
        self.grab_set()

    def show_message(self, text):
        self.message.set(text)
        self.update_idletasks()

    def close(self):
        self.progress.stop()
        self.grab_release()
        self.destroy()


class InfoWindow(tk.Toplevel):
    """A modeless, resizable window with selectable comment text."""
    def __init__(self, parent):
        super().__init__(parent)
        self.transient(parent.winfo_toplevel())
        self.geometry("760x430")
        self.minsize(420, 240)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.heading = ttk.Label(self, font=("TkDefaultFont", 11, "bold"), wraplength=700)
        self.heading.grid(row=0, column=0, columnspan=2, sticky="w", padx=16, pady=12)
        self.text = tk.Text(self, wrap="word", padx=12, pady=10, height=12,
                            relief="flat", state="disabled", font="TkDefaultFont")
        self.text.grid(row=1, column=0, sticky="nsew", padx=(16, 0))
        scrollbar = ttk.Scrollbar(self, command=self.text.yview)
        scrollbar.grid(row=1, column=1, sticky="ns", padx=(0, 16))
        self.text.configure(yscrollcommand=scrollbar.set)
        ttk.Button(self, text="Close", command=self.destroy).grid(
            row=2, column=0, columnspan=2, sticky="e", padx=16, pady=12)
        self.bind("<Escape>", lambda event: self.destroy())

    def show(self, title, text):
        self.title("Info - " + title)
        self.heading.configure(text=title)
        self.text.configure(state="normal")
        self.text.delete("1.0", "end")
        self.text.insert("1.0", text)
        self.text.configure(state="disabled")
        self.text.yview_moveto(0)
        self.deiconify()
        self.lift()
        self.text.focus_set()


def open_documentation(parent):
    """Open GUI-doc.pdf beside the GUI modules with the system PDF viewer."""
    path = Path(__file__).resolve().with_name("GUI-doc.pdf")
    if not path.is_file():
        messagebox.showerror("Documentation not found",
                             "Place GUI-doc.pdf in the GUI folder:\n{}".format(path.parent),
                             parent=parent)
        return
    try:
        if sys.platform == "win32":
            os.startfile(str(path))
            return
        command = "open" if sys.platform == "darwin" else "xdg-open"
        process = subprocess.Popen([command, str(path)], stdin=subprocess.DEVNULL,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as error:
        messagebox.showerror("Cannot open documentation", str(error), parent=parent)
        return

    def check_opened():
        result = process.poll()
        if result is None:
            parent.after(250, check_opened)
        elif result != 0:
            messagebox.showerror("Cannot open documentation",
                                 "The system PDF viewer could not open:\n{}".format(path),
                                 parent=parent)

    parent.after(250, check_opened)


class RunTab(ttk.Frame):
    def __init__(self, parent, before_start=None, call_atlas_selected=None, browse_directory=None):
        super().__init__(parent, padding=10)
        self.before_start = before_start
        self.call_atlas_selected = call_atlas_selected
        self.browse_directory = browse_directory
        self.engine = ProcessRunner()
        self._pending_exit = None
        self._ansi_pending = ""
        self._poll_id = None
        self.call_atlas_path = tk.StringVar()
        self.status = tk.StringVar(value="Ready")
        self.clean_cache = tk.BooleanVar(value=False)
        self.stage_action = tk.StringVar(value="Default")
        self.stage_warning = tk.StringVar(value="Keeps existing stages; no stage export.")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(5, weight=1)
        bar = ttk.Frame(self)
        bar.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        self.start_button = ttk.Button(bar, text="Start", command=self.start)
        self.start_button.pack(side="left")
        self.stop_button = ttk.Button(bar, text="Stop", command=self.stop, state="disabled")
        self.stop_button.pack(side="left", padx=8)
        ttk.Button(bar, text="GUI-doc", command=lambda: open_documentation(self)).pack(side="right")
        path_box = ttk.LabelFrame(self, text="call_atlas.ini", padding=8)
        path_box.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        self.path_font = tkfont.Font(root=self, name="TkDefaultFont", exists=True).copy()
        self.path_font.configure(weight="bold")
        ttk.Button(path_box, text="Browse ...", command=self.select_call_atlas).pack(side="right", padx=(8, 0))
        tk.Entry(path_box, textvariable=self.call_atlas_path, font=self.path_font,
                 background="#e8f2ff", foreground="#123e70", relief="flat").pack(fill="x", ipady=5)
        environment_box = ttk.LabelFrame(self, text="Python environment", padding=8)
        environment_box.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(environment_box, text=sys.prefix).pack(anchor="w")
        actions = ttk.LabelFrame(self, text="After processing", padding=8)
        actions.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        actions.columnconfigure(1, weight=1)
        self.cache_checkbox = ttk.Checkbutton(actions, text="Clean cached files", variable=self.clean_cache)
        self.cache_checkbox.grid(row=0, column=0, sticky="w", padx=(0, 12))
        ttk.Label(actions, text="Deletes reusable measurement caches. The next run may take longer. "
                  "Recommended on your final run.", wraplength=700).grid(row=0, column=1, sticky="w")
        ttk.Label(actions, text="Stage action").grid(row=1, column=0, sticky="w", pady=(8, 0))
        self.stage_selector = ttk.Combobox(actions, textvariable=self.stage_action, state="readonly",
            values=("Default", "Clean stages", "Export stages"), width=20)
        self.stage_selector.grid(row=1, column=1, sticky="w", pady=(8, 0))
        self.stage_selector.bind("<<ComboboxSelected>>", self.update_stage_warning)
        for event in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            self.stage_selector.bind(event, lambda event: "break")
        ttk.Label(actions, textvariable=self.stage_warning, wraplength=950).grid(
            row=2, column=0, columnspan=2, sticky="w", pady=(4, 0))
        ttk.Label(self, textvariable=self.status).grid(row=4, column=0, sticky="w", pady=(0, 8))
        frame = ttk.Frame(self)
        frame.grid(row=5, column=0, sticky="nsew")
        self.output = tk.Text(frame, wrap="word", state="disabled", height=16)
        scroll = ttk.Scrollbar(frame, command=self.output.yview)
        self.output.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        self.output.pack(fill="both", expand=True)
        self._poll_id = self.after(50, self.poll)

    def update_stage_warning(self, event=None):
        messages = {
            "Default": "Keeps existing stages; no stage export.",
            "Clean stages": "Warning: deletes all exported stages in this run's output folder. "
                            "They must be exported again if needed.",
            "Export stages": "Exports the export_stages selection from call_atlas.ini "
                             "(including its defaults). Existing stages with the same names are replaced.",
        }
        self.stage_warning.set(messages[self.stage_action.get()])

    def set_call_atlas_path(self, path):
        """Receive a file selection from the call_atlas.ini editor."""
        self.call_atlas_path.set(str(path) if path is not None else "")

    def get_call_atlas_path(self):
        """Read the current Run field when launching, without caching the editor path."""
        raw = self.call_atlas_path.get().strip()
        if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in ("'", '"'):
            raw = raw[1:-1]
        if not raw:
            raise ValueError("Select a call_atlas.ini file in the Run tab.")
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = Path(__file__).resolve().parent / path
        path = path.resolve()
        if not path.is_file():
            raise ValueError("The selected INI file does not exist:\n{}".format(path))
        return path

    def select_call_atlas(self):
        if self.engine.busy:
            messagebox.showinfo("ATLAS running", "Stop ATLAS before selecting another INI file.", parent=self)
            return
        try:
            current = self.get_call_atlas_path()
        except (OSError, ValueError, RuntimeError):
            directory = self.browse_directory() if self.browse_directory is not None else None
            current = (Path(directory) / "call_atlas.ini" if directory is not None
                       else Path(__file__).resolve().with_name("call_atlas.ini"))
        path = ask_open_file(parent=self, title="Select call_atlas.ini",
            initialdir=str(current.parent), initialfile=current.name,
            filetypes=[("INI files", "*.ini"), ("All files", "*")])
        if path:
            if self.call_atlas_selected is not None:
                self.call_atlas_selected(Path(path).resolve())
            else:
                self.set_call_atlas_path(Path(path).resolve())

    def running(self):
        return self.engine.busy

    def start(self):
        if self.engine.busy:
            return
        try:
            path = self.get_call_atlas_path()
            command = build_current_python_command(path)
            if self.before_start is not None and not self.before_start():
                return
            choices = {"clean_cache": self.clean_cache.get(), "stage_action": self.stage_action.get()}
            try:
                inspect.signature(self.engine.start).bind(command, path.parent, **choices)
            except TypeError as error:
                raise RuntimeError("The GUI and process runner are from different ATLAS updates. "
                    "Update gui_common.py and atlas_runner.py together, then restart atlas-gui.") from error
            # Resolve the INI at Start; always run this GUI's ATLAS installation.
            self.engine.start(command, path.parent, **choices)
        except (OSError, ValueError, RuntimeError) as error:
            messagebox.showerror("Cannot start ATLAS", str(error), parent=self)
            return
        self._pending_exit = None
        self._ansi_pending = ""
        display_command = subprocess.list2cmdline(command) if sys.platform == "win32" else shlex.join(command)
        self.append("\n--- Start ---\n{}\n\n".format(display_command))
        self.status.set("Running ATLAS ...")
        self.start_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.cache_checkbox.configure(state="disabled")
        self.stage_selector.configure(state="disabled")

    def stop(self):
        if self.engine.busy:
            self.engine.stop()
            self.status.set("Stopping ATLAS ...")
            self.stop_button.configure(state="disabled")

    def append(self, text):
        text = self._ansi_pending + text
        self._ansi_pending = ""
        text = re.sub(r"\x1b(?:\[[0-?]*[ -/]*[@-~]|\][^\x07]*(?:\x07|\x1b\\))", "", text)
        escape = text.rfind("\x1b")
        if escape >= 0 and len(text) - escape < 1024:
            text, self._ansi_pending = text[:escape], text[escape:]
        text = text.replace("\r\n", "\n").replace("\r", "\n")
        if not text:
            return
        follow = self.output.yview()[1] >= 0.999
        self.output.configure(state="normal")
        self.output.insert("end", text)
        length = self.output.count("1.0", "end-1c", "chars")[0]
        if length > 1000000:
            self.output.delete("1.0", "1.0+{}c".format(length - 1000000))
        if follow:
            self.output.see("end")
        self.output.configure(state="disabled")

    def poll(self):
        output = []
        for _ in range(128):
            try:
                kind, value = self.engine.events.get_nowait()
            except queue.Empty:
                break
            if kind == "text":
                output.append(value)
            elif kind == "done":
                self._pending_exit = value
            elif kind == "stop_error":
                output.append("\nCould not stop ATLAS: {}\n".format(value))
                self.status.set("Stop failed; see output")
                self.stop_button.configure(state="normal")
        if output:
            self.append("".join(output))
        if self._pending_exit is not None and self.engine.finish():
            code = self._pending_exit
            self._pending_exit = None
            result = "Stopped" if self.engine.stopping else ("Finished" if code == 0 else "Failed")
            self.status.set("{}, exit code: {}".format(result, code))
            self.append("\n--- {}, exit code: {} ---\n".format(result, code))
            self.start_button.configure(state="normal")
            self.stop_button.configure(state="disabled")
            self.cache_checkbox.configure(state="normal")
            self.stage_selector.configure(state="readonly")
        self._poll_id = self.after(50, self.poll)

    def destroy(self):
        if self._poll_id is not None:
            self.after_cancel(self._poll_id)
        super().destroy()
