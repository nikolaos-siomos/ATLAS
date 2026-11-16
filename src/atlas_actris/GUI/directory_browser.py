"""Portable list-view folder, file-open and file-save choosers."""
import fnmatch
import os
from pathlib import Path
import queue
import threading
import stat
import sys
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, font as tkfont


def background_call(function):
    """Run filesystem work only; workers never access Tk widgets."""
    result = queue.Queue(maxsize=1)
    def run():
        try:
            result.put((True, function()))
        except Exception as error:
            result.put((False, error))
    threading.Thread(target=run, daemon=True).start()
    return result


def is_hidden(entry):
    """Recognize Unix dot names and Windows hidden/system attributes."""
    if entry.name.startswith("."):
        return True
    if os.name == "nt":
        try:
            attributes = getattr(entry.stat(follow_symlinks=False), "st_file_attributes", 0)
            return bool(attributes & (getattr(stat, "FILE_ATTRIBUTE_HIDDEN", 2)
                                      | getattr(stat, "FILE_ATTRIBUTE_SYSTEM", 4)))
        except OSError:
            pass
    elif sys.platform == "darwin":
        try:
            return bool(getattr(entry.stat(follow_symlinks=False), "st_flags", 0)
                        & getattr(stat, "UF_HIDDEN", 0))
        except OSError:
            pass
    return False


def directory_entries(folder, show_hidden=False):
    """List one level, folders first; hidden entries require explicit opt-in."""
    folder = Path(folder).expanduser().resolve()
    entries = []
    with os.scandir(folder) as items:
        for item in items:
            if not show_hidden and is_hidden(item):
                continue
            try:
                is_dir = item.is_dir()
                kind = "Folder" if is_dir else "INI file" if item.name.lower().endswith(".ini") else "File"
            except OSError:
                is_dir, kind = False, "Unavailable"
            entries.append((item.name, Path(item.path), is_dir, kind))
    entries.sort(key=lambda row: (not row[2], row[0].casefold(), row[0]))
    return folder, entries


def selected_folder(current, text, selected=None, mustexist=True):
    """Honor an edited location first; clicking a file keeps its parent folder."""
    current = Path(current)
    requested = Path(text.strip()).expanduser() if text.strip() else current
    if not requested.is_absolute():
        requested = current / requested
    requested = requested.resolve()
    if requested == current.resolve() and selected is not None and Path(selected).is_dir():
        requested = Path(selected).resolve()
    if requested.exists() and not requested.is_dir():
        raise ValueError("Select a directory, not a file.")
    if mustexist and not requested.is_dir():
        raise ValueError("The selected directory does not exist.")
    return requested


def selected_file(current, text, save=False, defaultextension=""):
    """Validate a file selection without opening, creating or replacing it."""
    text = text.strip()
    if not text:
        raise ValueError("Choose a file or enter a file name.")
    requested = Path(text).expanduser()
    if not requested.is_absolute():
        requested = Path(current) / requested
    if requested.is_dir():
        raise ValueError("Choose a file, not a folder.")
    if save and defaultextension and not requested.suffix:
        requested = requested.with_name(requested.name + defaultextension)
    requested = requested.resolve()
    if requested.is_dir():
        raise ValueError("Choose a file, not a folder.")
    if save:
        if not requested.parent.is_dir():
            raise ValueError("The destination folder does not exist.")
        if requested.exists() and not requested.is_file():
            raise ValueError("The destination is not a regular file.")
    elif not requested.is_file():
        raise ValueError("The selected file does not exist or is not a regular file.")
    return requested


def matching_entries(entries, patterns):
    """Retain folders for navigation and files matching the selected filter."""
    if isinstance(patterns, str):
        patterns = patterns.split()
    return [row for row in entries if row[2] or any(
        fnmatch.fnmatchcase(row[0].casefold(), pattern.casefold()) for pattern in patterns)]


class DirectoryDialog(tk.Toplevel):
    def __init__(self, parent, title, initialdir, mustexist=True, mode="directory",
                 initialfile="", filetypes=None, defaultextension=""):
        super().__init__(parent)
        self.withdraw()
        self.title(title)
        self.transient(parent.winfo_toplevel())
        self.geometry("850x550")
        self.minsize(550, 350)
        self.result = ""
        self.mustexist = mustexist
        self.mode = mode
        self.defaultextension = defaultextension
        self.filetypes = filetypes or [("All files", "*")]
        self.filename = tk.StringVar(value=initialfile)
        self.filetype = tk.StringVar(value=self.filetypes[0][0])
        self.current = Path(initialdir or Path.home()).expanduser()
        if not self.current.is_dir():
            self.current = Path.home()
        self.current = self.current.resolve()
        self.rows = {}
        self.busy = False
        self.poll_id = None
        self.location = tk.StringVar(value=str(self.current))
        self.status = tk.StringVar(value="Files are visible for reference. Choose a folder.")
        self.show_hidden = tk.BooleanVar(value=False)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)
        location = ttk.Frame(self, padding=10)
        location.grid(row=0, column=0, sticky="ew")
        location.columnconfigure(0, weight=1)
        self.entry = ttk.Entry(location, textvariable=self.location)
        self.entry.grid(row=0, column=0, sticky="ew")
        self.entry.bind("<Return>", lambda event: self.navigate(self.location.get()))
        self.go = ttk.Button(location, text="Go", command=lambda: self.navigate(self.location.get()))
        self.go.grid(row=0, column=1, padx=(8, 0))
        toolbar = ttk.Frame(self, padding=(10, 0, 10, 8))
        toolbar.grid(row=1, column=0, sticky="ew")
        self.up = ttk.Button(toolbar, text="Up", command=lambda: self.navigate(self.current.parent))
        self.up.pack(side="left")
        self.home = ttk.Button(toolbar, text="Home", command=lambda: self.navigate(Path.home()))
        self.home.pack(side="left", padx=8)
        self.refresh = ttk.Button(toolbar, text="Refresh", command=lambda: self.navigate(self.current))
        self.refresh.pack(side="left")
        self.new = ttk.Button(toolbar, text="New folder", command=self.new_folder)
        self.new.pack(side="right")
        self.hidden_toggle = ttk.Checkbutton(toolbar, text="Show hidden files", variable=self.show_hidden,
                                            command=lambda: self.navigate(self.current))
        self.hidden_toggle.pack(side="right", padx=12)
        content = ttk.Frame(self, padding=(10, 0))
        content.grid(row=2, column=0, sticky="nsew")
        content.rowconfigure(0, weight=1)
        content.columnconfigure(0, weight=1)
        style = ttk.Style(self)
        height = tkfont.Font(root=self, name="TkDefaultFont", exists=True).metrics("linespace") + 8
        style.configure("AtlasFolders.Treeview", rowheight=height)
        self.table = ttk.Treeview(content, columns=("name", "kind"), show="headings",
                                  selectmode="browse", style="AtlasFolders.Treeview")
        self.table.heading("name", text="Name")
        self.table.heading("kind", text="Type")
        self.table.column("name", width=570, minwidth=160)
        self.table.column("kind", width=140, minwidth=90, stretch=False)
        self.table.grid(row=0, column=0, sticky="nsew")
        vertical = ttk.Scrollbar(content, orient="vertical", command=self.table.yview)
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal = ttk.Scrollbar(content, orient="horizontal", command=self.table.xview)
        horizontal.grid(row=1, column=0, sticky="ew")
        self.table.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        self.table.bind("<Double-1>", self.open_selected)
        self.table.bind("<Return>", self.open_selected)
        if self.mode != "directory":
            self.table.bind("<<TreeviewSelect>>", self.select_file_row)
        ttk.Label(self, textvariable=self.status, wraplength=800, padding=10).grid(row=3, column=0, sticky="ew")
        footer = ttk.Frame(self, padding=(10, 0, 10, 10))
        footer.grid(row=4, column=0, sticky="ew")
        if self.mode != "directory":
            fields = ttk.Frame(footer)
            fields.pack(fill="x", pady=(0, 10))
            fields.columnconfigure(1, weight=1)
            ttk.Label(fields, text="File name").grid(row=0, column=0, sticky="w", padx=(0, 8))
            self.filename_entry = ttk.Entry(fields, textvariable=self.filename)
            self.filename_entry.grid(row=0, column=1, sticky="ew")
            self.filename_entry.bind("<Return>", lambda event: self.accept())
            self.filter_selector = ttk.Combobox(fields, textvariable=self.filetype, state="readonly",
                values=[item[0] for item in self.filetypes], width=18)
            self.filter_selector.grid(row=0, column=2, padx=(8, 0))
            self.filter_selector.bind("<<ComboboxSelected>>", lambda event: self.navigate(self.current))
            for event in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                self.filter_selector.bind(event, lambda event: "break")
        ttk.Button(footer, text="Cancel", command=self.destroy).pack(side="right")
        label = "Select folder" if mode == "directory" else "Save" if mode == "save" else "Open"
        self.select = ttk.Button(footer, text=label, command=self.accept)
        self.select.pack(side="right", padx=8)
        self.bind("<Escape>", lambda event: self.destroy())
        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.deiconify()
        self.wait_visibility()
        self.grab_set()
        self.navigate(self.current)

    def controls(self, state):
        for widget in (self.entry, self.go, self.up, self.home, self.refresh, self.new, self.select, self.hidden_toggle):
            widget.configure(state=state)
        if self.mode != "directory":
            self.filename_entry.configure(state=state)
            self.filter_selector.configure(state="readonly" if state == "normal" else state)

    def navigate(self, folder):
        if self.busy:
            return "break"
        try:
            path = Path(folder).expanduser()
        except (OSError, ValueError, RuntimeError) as error:
            messagebox.showerror("Cannot open folder", str(error), parent=self)
            return "break"
        if not path.is_absolute():
            path = self.current / path
        self.busy = True
        self.controls("disabled")
        self.status.set("Loading folder: " + str(path))
        show_hidden = self.show_hidden.get()  # Read Tk variables only on the UI thread.
        self.pending = background_call(lambda: directory_entries(path, show_hidden=show_hidden))
        self.poll_id = self.after(30, self.poll_directory)
        return "break"

    def poll_directory(self):
        self.poll_id = None
        try:
            success, payload = self.pending.get_nowait()
        except queue.Empty:
            self.poll_id = self.after(30, self.poll_directory)
            return
        self.busy = False
        self.controls("normal")
        if not success:
            self.status.set("Could not open folder: " + str(payload))
            messagebox.showerror("Cannot open folder", str(payload), parent=self)
            return
        self.current, entries = payload
        if self.mode != "directory":
            patterns = next(item[1] for item in self.filetypes if item[0] == self.filetype.get())
            entries = matching_entries(entries, patterns)
        self.location.set(str(self.current))
        self.table.delete(*self.table.get_children())
        self.rows = {}
        for index, (name, path, is_dir, kind) in enumerate(entries):
            key = str(index)
            self.table.insert("", "end", iid=key, values=(name, kind))
            self.rows[key] = (path, is_dir)
        hint = ("Files are shown for reference; double-click a folder to open it." if self.mode == "directory"
                else "Double-click a folder to open it, or select a file and press {}.".format(
                    "Save" if self.mode == "save" else "Open"))
        self.status.set("{} folders, {} files. {}".format(
            sum(row[2] for row in entries), sum(not row[2] for row in entries), hint))

    def select_file_row(self, event=None):
        selection = self.table.selection()
        if not self.busy and selection and selection[0] in self.rows:
            path, is_dir = self.rows[selection[0]]
            if not is_dir:
                self.filename.set(path.name)

    def open_selected(self, event=None):
        if not self.busy:
            selection = self.table.selection()
            if selection and selection[0] in self.rows:
                path, is_dir = self.rows[selection[0]]
                if is_dir:
                    self.navigate(path)
                elif self.mode != "directory":
                    self.filename.set(path.name)
                    self.accept()
        return "break"

    def accept(self):
        if self.busy:
            return
        if self.mode != "directory":
            try:
                result = selected_file(self.current, self.filename.get(), self.mode == "save",
                                       self.defaultextension)
            except (OSError, ValueError, RuntimeError) as error:
                messagebox.showerror("Cannot select file", str(error), parent=self)
                return
            if self.mode == "save" and result.exists() and not messagebox.askyesno(
                    "Replace existing file?", "Replace this file?\n{}".format(result), parent=self):
                return
            self.result = str(result)
            self.destroy()
            return
        selection = self.table.selection()
        path = self.rows[selection[0]][0] if selection else None
        try:
            result = selected_folder(self.current, self.location.get(), path, self.mustexist)
        except (OSError, ValueError, RuntimeError) as error:
            messagebox.showerror("Cannot select folder", str(error), parent=self)
            return
        self.result = str(result)
        self.destroy()

    def new_folder(self):
        name = simpledialog.askstring("New folder", "Folder name:", parent=self)
        if name is None:
            return
        name = name.strip()
        if not name or name in (".", "..") or any(char in "/\\:" or ord(char) < 32 for char in name):
            messagebox.showerror("Cannot create folder", "Enter a folder name without path separators.", parent=self)
            return
        try:
            path = self.current / name
            path.mkdir()
        except OSError as error:
            messagebox.showerror("Cannot create folder", str(error), parent=self)
            return
        self.navigate(path)

    def destroy(self):
        if self.poll_id is not None:
            self.after_cancel(self.poll_id)
            self.poll_id = None
        super().destroy()


def ask_directory(parent, title="Select folder", initialdir=None, mustexist=True):
    dialog = DirectoryDialog(parent, title, initialdir, mustexist)
    parent.wait_window(dialog)
    return dialog.result


def ask_open_file(parent, title="Select file", initialdir=None, initialfile="", filetypes=None):
    dialog = DirectoryDialog(parent, title, initialdir, mode="open",
                             initialfile=initialfile, filetypes=filetypes)
    parent.wait_window(dialog)
    return dialog.result


def ask_save_file(parent, title="Save file", initialdir=None, initialfile="",
                  filetypes=None, defaultextension=""):
    dialog = DirectoryDialog(parent, title, initialdir, mode="save", initialfile=initialfile,
                             filetypes=filetypes, defaultextension=defaultextension)
    parent.wait_window(dialog)
    return dialog.result
