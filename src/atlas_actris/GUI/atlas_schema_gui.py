#!/usr/bin/env python3
"""English, schema-driven ATLAS INI GUI. Python 3.9-3.13; Tkinter only."""
from pathlib import Path
import queue
import sys

# Resolve local helpers and sibling ATLAS sources from this file, not the working directory.
BASE = Path(__file__).resolve().parent
if not __package__ and str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

import tkinter as tk
from tkinter import messagebox, ttk, font as tkfont
if __package__:
    from .gui_common import InfoWindow, RunTab, ScrollFrame, LoadingDialog, open_documentation
    from .directory_browser import ask_directory, ask_open_file, ask_save_file, background_call
    from .schema_loader import load_schemas, display
    from .editor_model import EditorModel, normalize_value, load_working_folder, load_call_atlas_file, linked_ini_path
    from .preferences import defaults as preference_defaults, load_preferences, save_preferences, remembered_directory
else:
    from gui_common import InfoWindow, RunTab, ScrollFrame, LoadingDialog, open_documentation
    from directory_browser import ask_directory, ask_open_file, ask_save_file, background_call
    from schema_loader import load_schemas, display
    from editor_model import EditorModel, normalize_value, load_working_folder, load_call_atlas_file, linked_ini_path
    from preferences import defaults as preference_defaults, load_preferences, save_preferences, remembered_directory

DEFAULT_SOURCE = BASE.parent  # src/atlas_actris; utils and templates are siblings of GUI.
HIDDEN_CONFIG_SECTIONS = {"water_vapour", "gluing", "temperature"}


def configure_fonts(parent):
    """Use platform-native interface fonts, at accessible point sizes."""
    parent.interface_fonts = []
    for name in ("TkDefaultFont", "TkTextFont", "TkMenuFont", "TkHeadingFont",
                 "TkCaptionFont", "TkSmallCaptionFont", "TkTooltipFont", "TkFixedFont"):
        font = tkfont.Font(root=parent, name=name, exists=True)
        font.configure(size=13 if name == "TkHeadingFont" else 12)
        parent.interface_fonts.append(font)
    style = ttk.Style(parent)
    style.configure("TButton", padding=(10, 5))
    style.configure("TEntry", padding=(6, 4))
    parent.option_add("*TCombobox*Listbox.font", "TkTextFont")


def path_kind(field):
    # output_folder is created by ATLAS; do not change parser validation to
    # require an existing directory just to display a GUI browse button.
    return "dir" if field.key == "output_folder" else field.meta.get("check_path")


def configure_section_tabs(parent):
    style = ttk.Style(parent)
    parent.tab_label_font = tkfont.Font(root=parent, name="TkDefaultFont", exists=True).copy()
    parent.tab_label_font.configure(weight="bold", size=13)
    style.configure("AtlasMain.TNotebook.Tab", font=parent.tab_label_font, padding=(12, 6))
    style.map("AtlasMain.TNotebook.Tab", foreground=[("selected", "#233c56")])
    # Use a colorable tab element even with native Windows/macOS themes.
    element = "AtlasSections.Notebook.tab"
    if element not in style.element_names():
        style.element_create(element, "from", "clam", "Notebook.tab")
    style.layout("AtlasSections.TNotebook.Tab", [
        (element, {"sticky": "nswe", "children": [
            ("Notebook.padding", {"side": "top", "sticky": "nswe", "children": [
                ("Notebook.focus", {"side": "top", "sticky": "nswe", "children": [
                    ("Notebook.label", {"side": "top", "sticky": ""})
                ]})
            ]})
        ]})
    ])
    style.configure("AtlasSections.TNotebook.Tab", padding=(12, 6),
                    background="#f5f6f7", foreground="#485361",
                    bordercolor="#d8dee5", lightcolor="#e9edf1", darkcolor="#d8dee5")
    style.map("AtlasSections.TNotebook.Tab",
              background=[("selected", "#e8edf3"), ("active", "#eff2f5"), ("!selected", "#f5f6f7")],
              foreground=[("selected", "#233c56"), ("!selected", "#485361")],
              font=[("selected", parent.tab_label_font)])


class SchemaTab(ttk.Frame):
    def __init__(self, parent, schema, recommended, optional, filter_changed,
                 path_changed=None, browse_directory=None):
        super().__init__(parent, padding=10)
        self.model = EditorModel(schema)
        self.filename = schema.filename
        self.recommended, self.optional = recommended, optional
        self.path_changed = path_changed
        self.browse_directory = browse_directory
        self.info_window = None
        self._refreshing = False
        self.variables = {}
        self.rows = []
        self.page_notes = []
        self.status = tk.StringVar()
        self.include_descriptions = tk.BooleanVar(value=False)
        self.path_display = tk.StringVar(value="No INI file loaded")
        self.path_font = tkfont.Font(root=self, name="TkDefaultFont", exists=True).copy()
        self.path_font.configure(weight="bold")
        bar = ttk.Frame(self)
        bar.pack(fill="x")
        for title, action in (("New", self.new), ("Open ...", self.open), ("Save", self.save),
                              ("Save as ...", lambda: self.save(True)), ("Check ini-file entries", self.validate)):
            ttk.Button(bar, text=title, command=action).pack(side="left", padx=(0, 6))
        ttk.Button(bar, text="GUI-doc", command=lambda: open_documentation(self)).pack(side="right")
        path_box = ttk.LabelFrame(self, text="Loaded INI file", padding=8)
        path_box.pack(fill="x", pady=(10, 0))
        tk.Entry(path_box, textvariable=self.path_display, state="readonly",
                 font=self.path_font, readonlybackground="#e8f2ff",
                 foreground="#123e70", relief="flat").pack(fill="x", ipady=5)
        filters = ttk.Frame(self)
        filters.pack(fill="x", pady=(10, 0))
        ttk.Checkbutton(filters, text="Show recommended inputs", variable=self.recommended,
                        command=filter_changed).pack(side="left")
        ttk.Checkbutton(filters, text="Show optional inputs", variable=self.optional,
                        command=filter_changed).pack(side="left", padx=16)
        ttk.Checkbutton(filters, text="Include parameter descriptions on save",
                        variable=self.include_descriptions).pack(side="left", padx=8)
        ttk.Label(self, textvariable=self.status, wraplength=1000).pack(side="bottom", fill="x", pady=(8, 0))
        self.notebook = ttk.Notebook(self, style="AtlasSections.TNotebook")
        self.notebook.pack(fill="both", expand=True, pady=(10, 0))
        self.build()

    def build(self):
        self.variables = {}
        self.rows = []
        self.page_notes = []
        for child in self.notebook.winfo_children():
            child.destroy()
        for section, fields in self.model.schema.sections.items():
            if self.filename == "call_atlas.ini" and section == "configuration":
                # Keep legacy values in the model and save them if provided.
                continue
            page = ScrollFrame(self.notebook)
            self.notebook.add(page, text=section)
            if self.filename == "config_file.ini" and section in HIDDEN_CONFIG_SECTIONS:
                self.notebook.hide(page)
            page.body.columnconfigure(0, weight=1)
            for row_number, field in enumerate(fields):
                row = ttk.Frame(page.body)
                row.grid(row=row_number, column=0, sticky="ew", pady=4)
                row.columnconfigure(1, weight=1)
                label = ttk.Frame(row, width=270)
                label.grid(row=0, column=0, sticky="nw", padx=(0, 12))
                ttk.Label(label, text=field.key, width=34, wraplength=260).pack(anchor="w")
                category = "Required" if field.category == "mandatory" else field.category.capitalize()
                ttk.Label(label, text=category, foreground="#666666").pack(anchor="w")
                cell = ttk.Frame(row)
                cell.grid(row=0, column=1, sticky="ew")
                variable = tk.StringVar(value=self.model.values[field.identity])
                self.variables[field.identity] = variable
                allowed = field.meta.get("allowed")
                if not field.meta.get("is_list") and field.meta["dtype"] is bool:
                    widget = ttk.Combobox(cell, textvariable=variable, values=("", "True", "False"))
                elif not field.meta.get("is_list") and isinstance(allowed, (list, tuple)):
                    widget = ttk.Combobox(cell, textvariable=variable, values=[""] + [str(v) for v in allowed])
                else:
                    widget = ttk.Entry(cell, textvariable=variable)
                widget.pack(fill="x")
                widget.bind("<FocusOut>", lambda event, f=field: self.normalize_field(f))
                default = display(field.meta.get("default"))
                if len(default) > 85:
                    default = default[:82] + "..."
                ttk.Label(cell, text="Default: " + default, foreground="#666666", wraplength=560).pack(anchor="w")
                controls = ttk.Frame(row)
                controls.grid(row=0, column=2, sticky="ne", padx=(8, 0))
                ttk.Button(controls, text="Info", width=5, command=lambda f=field: self.show_info(
                    "[{}] {}".format(f.section, f.key), f.help_text())).pack(side="left")
                if path_kind(field) in ("file", "dir"):
                    ttk.Button(controls, text="Browse", width=7, command=lambda f=field: self.browse_value(f)).pack(side="left", padx=(4, 0))
                variable.trace_add("write", lambda *_, f=field, v=variable: self.changed(f.identity, v.get()))
                self.rows.append((field, row))
            note = ttk.Label(page.body, text="No settings match the current filters.", wraplength=750)
            note.grid(row=len(fields), column=0, sticky="w", pady=12)
            self.page_notes.append((fields, note))
            page.bind_scroll()
        self.apply_filters()
        self.update_status()

    def changed(self, identity, value):
        if self._refreshing:
            return
        self.model.values[identity] = value
        self.update_status()

    def refresh_values(self):
        """Reuse schema widgets when only the loaded values have changed."""
        self._refreshing = True
        try:
            for identity, variable in self.variables.items():
                value = self.model.values[identity]
                if variable.get() != value:
                    variable.set(value)
        finally:
            self._refreshing = False
        self.update_status()

    def normalize_field(self, field):
        variable = self.variables[field.identity]
        try:
            value = normalize_value(field, variable.get())
        except ValueError:
            return  # Check/Save will report line breaks, without a Tk traceback.
        if variable.get() != value:
            variable.set(value)

    def normalize_fields(self):
        for field in self.model.schema.fields:
            if field.identity in self.variables:
                self.normalize_field(field)
        self.model.normalize()

    def dirty(self):
        return self.model.dirty()

    def update_status(self):
        self.path_display.set(str(self.model.path) if self.model.path else "No INI file loaded")
        self.status.set(("Loaded" if self.model.path else "New file from Python schema") +
                        ("  |  Unsaved changes" if self.dirty() else "") +
                        ("  |  {} unrecognized entries preserved; click Check ini-file entries".format(len(self.model.unknown)) if self.model.unknown else ""))

    def apply_filters(self):
        for field, row in self.rows:
            if field.visible(self.recommended.get(), self.optional.get()):
                row.grid()
            else:
                row.grid_remove()
        for fields, note in self.page_notes:
            if any(f.visible(self.recommended.get(), self.optional.get()) for f in fields):
                note.grid_remove()
            else:
                note.grid()

    def show_info(self, title, text):
        if self.info_window is None or not self.info_window.winfo_exists():
            self.info_window = InfoWindow(self)
        self.info_window.show(title, text)

    def may_discard(self):
        return not self.dirty() or messagebox.askyesno("Unsaved changes", "Discard changes to {}?".format(self.filename), parent=self)

    def notify_path_changed(self):
        if self.path_changed is not None:
            self.path_changed(self.model.path)

    def new(self):
        if self.may_discard():
            self.model = EditorModel(self.model.schema)
            self.refresh_values()
            self.notify_path_changed()

    def open(self):
        if not self.may_discard():
            return
        directory = self.browse_directory() if self.browse_directory is not None else None
        if directory is None:
            directory = self.model.path.parent if self.model.path else BASE / "examples"
        path = ask_open_file(parent=self, title="Open " + self.filename,
            initialfile=self.filename, initialdir=str(directory),
            filetypes=[("INI files", "*.ini"), ("All files", "*")])
        if path:
            loading = LoadingDialog(self, "Loading " + str(path))
            problem = None
            try:
                self.model.load(path)
                loading.show_message("Preparing " + self.filename + "...")
                self.refresh_values()
                self.notify_path_changed()
            except (OSError, ValueError) as error:
                problem = error
            finally:
                loading.close()
            if problem is not None:
                messagebox.showerror("Load failed", str(problem), parent=self)

    def save(self, save_as=False):
        previous_path = self.model.path
        target = self.model.path
        if target is None or save_as:
            suggested = target or self.model.suggested_save_path
            directory = suggested.parent if suggested is not None else None
            if directory is None and self.browse_directory is not None:
                directory = self.browse_directory()
            name = ask_save_file(parent=self, title="Save " + self.filename,
                initialdir=str(directory or Path.home()),
                initialfile=suggested.name if suggested is not None else self.filename,
                defaultextension=".ini", filetypes=[("INI files", "*.ini")])
            if not name:
                return False
            target = Path(name)
        try:
            self.normalize_fields()
            self.model.save(target, include_descriptions=self.include_descriptions.get())
            self.update_status()
            if save_as or self.model.path != previous_path:
                self.notify_path_changed()
            errors = self.model.validate()
            if errors:
                self.status.set(self.status.get() + "  |  Draft saved: {} validation issues; click Check ini-file entries".format(len(errors)))
            return True
        except (OSError, ValueError) as error:
            messagebox.showerror("Save failed", str(error), parent=self)
            return False

    def validate(self):
        try:
            self.normalize_fields()
        except ValueError as error:
            self.show_info("Validation - " + self.filename, str(error))
            return
        errors = self.model.validate()
        self.show_info("Validation - " + self.filename, "\n\n".join(errors) if errors else
                       "Declared field constraints passed. ATLAS runtime and cross-field checks are not run here.")

    def browse_value(self, field):
        directory = self.model.path.parent if self.model.path else None
        if directory is None and getattr(self, "browse_directory", None) is not None:
            directory = self.browse_directory()
        directory = directory or Path.home()
        if path_kind(field) == "dir":
            value = ask_directory(parent=self, title="Select " + field.key,
                                           initialdir=directory,
                                           mustexist=field.key != "output_folder")
        else:
            value = ask_open_file(parent=self, title="Select " + field.key,
                                  initialdir=directory)
        if value:
            self.variables[field.identity].set(value)


class App(tk.Tk):
    def __init__(self, startup_prompt=False, remember_directory=True):
        super().__init__()
        self.title("ATLAS - Schema INI Editor")
        self.geometry("1200x800")
        self.minsize(850, 550)
        self._loading = True
        # Show the busy screen before reading schemas, configuring styles, or
        # creating the thousands of controls used by the schema editors.
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.loading_dialog = LoadingDialog(self, "Starting ATLAS GUI...")
        self.directory_memory_enabled = remember_directory
        self.preferences = load_preferences() if remember_directory else preference_defaults()
        try:
            self.initialize_interface()
        finally:
            self.end_folder_load()
        if startup_prompt:
            self.after(100, self.select_working_folder)

    def initialize_interface(self):
        self.loading_dialog.show_message("Preparing interface fonts and styles...")
        configure_fonts(self)
        configure_section_tabs(self)
        self.loading_dialog.show_message("Reading INI schemas and descriptions...")
        self.schemas = load_schemas(DEFAULT_SOURCE)
        self.recommended = tk.BooleanVar(value=True)
        self.optional = tk.BooleanVar(value=True)
        self.working_folder = None
        self.folder_display = tk.StringVar(value="No INI working folder selected")
        folder_bar = ttk.Frame(self, padding=(10, 10, 10, 0))
        folder_bar.pack(fill="x")
        ttk.Button(folder_bar, text="Choose INI folder ...", command=self.select_working_folder).pack(side="left")
        ttk.Label(folder_bar, textvariable=self.folder_display, wraplength=850).pack(side="left", padx=12)
        self.notebook = ttk.Notebook(self, style="AtlasMain.TNotebook")
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)
        self.editors = []
        self.runner = RunTab(self.notebook, before_start=self.ready_to_run,
                             call_atlas_selected=self.select_call_atlas,
                             browse_directory=self.call_atlas_directory)
        self.notebook.add(self.runner, text="Run")
        self.install_schemas(self.schemas)

    def select_working_folder(self):
        if self._loading:
            return
        if self.runner.running():
            messagebox.showinfo("ATLAS running", "Stop ATLAS before changing the INI working folder.", parent=self)
            return
        folder = ask_directory(parent=self, title="Select the INI working folder",
                               initialdir=self.working_folder or remembered_directory(
                                   getattr(self, "preferences", {})) or Path.home(), mustexist=True)
        if not folder or any(not editor.may_discard() for editor in self.editors):
            return
        self.begin_folder_load(folder)

    def select_call_atlas(self, path):
        if self._loading:
            return
        if self.runner.running():
            messagebox.showinfo("ATLAS running", "Stop ATLAS before selecting another INI file.", parent=self)
            return
        if any(not editor.may_discard() for editor in self.editors):
            return
        self.begin_folder_load(Path(path).parent, call_atlas=path)

    def begin_folder_load(self, folder, call_atlas=None):
        self._loading = True
        self.loading_dialog = LoadingDialog(self, "Reading INI files from " + str(folder))
        self.pending_load = background_call(
            lambda: load_call_atlas_file(call_atlas, self.schemas) if call_atlas is not None
            else load_working_folder(folder, self.schemas))
        self._load_folder = Path(folder)
        self.after(30, self.poll_folder_load)

    def poll_folder_load(self):
        try:
            success, payload = self.pending_load.get_nowait()
        except queue.Empty:
            self.after(30, self.poll_folder_load)
            return
        if not success:
            self.end_folder_load()
            messagebox.showerror("Cannot open INI folder", str(payload), parent=self)
            return
        self._loaded_models, self._load_actions = payload
        self._editor_index = 0
        self.after(10, self.prepare_next_editor)

    def prepare_next_editor(self):
        if self._editor_index < len(self.editors):
            editor = self.editors[self._editor_index]
            self.loading_dialog.show_message("Preparing " + editor.filename + "...")
            try:
                editor.model = self._loaded_models[editor.filename]
                editor.refresh_values()
                editor.notify_path_changed()
            except Exception as error:
                self.end_folder_load()
                messagebox.showerror("Cannot prepare INI editor", str(error), parent=self)
                return
            self._editor_index += 1
            self.after(10, self.prepare_next_editor)
            return
        self.working_folder = self._load_folder
        self.folder_display.set(str(self.working_folder))
        self.notebook.select(self.editors[0])
        if getattr(self, "directory_memory_enabled", False):
            self.remember_ini_directory(self._loaded_models["call_atlas.ini"].path)
        detail = ("Created files are empty drafts; fill in any required fields before processing."
                  if any(action.startswith("Created ") for action in self._load_actions)
                  else "Existing files were loaded without modifying them.")
        self.end_folder_load()
        locations = "\n".join("{}: {}".format(name, model.path or "Unsaved draft")
                              for name, model in self._loaded_models.items())
        warnings = [model.load_warning for model in self._loaded_models.values() if model.load_warning]
        message = "\n".join(self._load_actions) + "\n\n" + locations + "\n\n" + detail
        if warnings:
            messagebox.showwarning("INI files need attention", message + "\n\n" + "\n\n".join(warnings), parent=self)
        else:
            messagebox.showinfo("INI working folder", message, parent=self)

    def end_folder_load(self):
        self.loading_dialog.close()
        self._loading = False

    def install_schemas(self, schemas):
        paths = {e.filename: e.model.path for e in self.editors}
        for editor in self.editors:
            editor.destroy()
        self.editors = []
        for index, schema in enumerate(schemas.values()):
            if self._loading:
                self.loading_dialog.show_message("Preparing " + schema.filename + "...")
            path_changed = self.init_path_changed if schema.filename == "call_atlas.ini" else None
            if schema.filename in ("config_file.ini", "settings_file.ini"):
                path_changed = lambda path, filename=schema.filename: self.linked_draft_saved(filename, path)
            browse_directory = self.call_atlas_directory
            editor = SchemaTab(self.notebook, schema, self.recommended, self.optional,
                               self.filter_changed, path_changed, browse_directory)
            path = paths.get(schema.filename)
            if path is not None:
                try:
                    editor.model.load(path)
                    editor.build()
                except (OSError, ValueError) as error:
                    messagebox.showwarning("Could not reload " + schema.filename, str(error), parent=self)
            editor.notify_path_changed()
            self.editors.append(editor)
            self.notebook.insert(index, editor, text=schema.filename)
        self.notebook.select(self.editors[0])

    def call_atlas_directory(self):
        if self.working_folder is not None:
            return self.working_folder
        for editor in self.editors:
            if editor.filename == "call_atlas.ini" and editor.model.path is not None:
                directory = editor.model.path.parent
                if directory.is_dir():
                    return directory
        # A file may have been selected only in Run, without loading the editor.
        try:
            return self.runner.get_call_atlas_path().parent
        except (OSError, ValueError, RuntimeError):
            return remembered_directory(getattr(self, "preferences", {}))

    def remember_ini_directory(self, path):
        if not getattr(self, "directory_memory_enabled", False) or path is None:
            return
        self.preferences["last_directory"] = str(Path(path).parent)
        save_preferences(self.preferences)

    def init_path_changed(self, path):
        self.runner.set_call_atlas_path(path)
        if not self._loading:
            self.remember_ini_directory(path)

    def filter_changed(self):
        for editor in self.editors:
            editor.apply_filters()

    def linked_draft_saved(self, filename, path):
        if self._loading or path is None:
            return
        editor = next(e for e in self.editors if e.filename == filename)
        if not editor.model.relink_on_save:
            return
        init = next(e for e in self.editors if e.filename == "call_atlas.ini")
        if init.model.path is None:
            return
        key = "atlas_configuration_file" if filename == "config_file.ini" else "atlas_settings_file"
        identity = ("explicit_paths", key)
        try:
            expected = linked_ini_path(init.model.path, init.model.values[identity], filename).resolve()
        except (OSError, ValueError, RuntimeError):
            expected = None
        if expected != Path(path).resolve():
            init.variables[identity].set(str(path))
            init.model.values[identity] = str(path)
            init.update_status()
            messagebox.showinfo("Save the init file", "Updated {} to the saved draft's path. "
                "Save {} before running ATLAS.".format(key, init.model.path.name), parent=self)
        editor.model.relink_on_save = False
        editor.model.load_warning = ""

    def ready_to_run(self):
        if self._loading:
            return False
        unsaved = [editor.filename for editor in self.editors if editor.model.path is None]
        if unsaved:
            messagebox.showinfo("INI draft needs saving", "These tabs have no saved file: {}. "
                "Fill in and save the drafts, or correct the file paths in the init and select it again "
                "before running ATLAS.".format(", ".join(unsaved)), parent=self)
            return False
        if any(editor.dirty() for editor in self.editors):
            messagebox.showinfo("Unsaved changes", "Save the modified INI files before starting ATLAS.", parent=self)
            return False
        return True

    def close(self):
        if self._loading:
            return
        if self.runner.running() and not messagebox.askyesno("ATLAS running",
                "Stop ATLAS and close the GUI?", parent=self):
            return
        if any(not editor.may_discard() for editor in self.editors):
            return
        if self.runner.running():
            self.runner.stop()
            self.after(100, self.finish_close)
            return
        self.destroy()

    def finish_close(self):
        if self.runner.engine.stop_failed:
            return
        if self.runner.running():
            self.after(100, self.finish_close)
        else:
            self.destroy()


def main():
    try:
        App().mainloop()
    except tk.TclError as error:
        print("Cannot start the ATLAS GUI: {}\n"
              "A graphical desktop/display and working Tcl/Tk are required. "
              "Check this Python environment with: python -m tkinter".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
