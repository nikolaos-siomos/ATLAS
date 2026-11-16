"""Schema-based values, validation and compact, normalized INI saves."""
import math
import ntpath
import os
import re
from pathlib import Path
if __package__:
    from .gui_common import IniDocument, Setting, atomic_write
    from .schema_loader import Unresolved, display
else:
    from gui_common import IniDocument, Setting, atomic_write
    from schema_loader import Unresolved, display


def normalize_value(field, raw):
    """Trim tokens without dropping empty list positions or internal spaces."""
    if "\n" in raw or "\r" in raw:
        raise ValueError("Values must not contain line breaks.")
    text = raw.strip()
    if field.meta.get("is_list") and text:
        separator = field.value_syntax.separator
        text = separator.join(token.strip() for token in text.replace(";", separator).split(separator))
    return text.rstrip()


def validate_value(field, raw, filename, required=False):
    """Check declared field constraints, not ATLAS cross-field/runtime logic."""
    if "\n" in raw or "\r" in raw:
        return ["Use a single line for this value."]
    text = normalize_value(field, raw)
    if not text:
        return ["A value is required."] if required and field.category == "mandatory" else []
    meta = field.meta
    syntax = field.value_syntax
    if meta.get("is_list"):
        tokens = [t.strip() for t in text.replace(";", syntax.separator).split(syntax.separator)]
        if syntax.drop_empty:
            tokens = [t for t in tokens if t]
    else:
        tokens = [text]
    values = []
    for token in tokens:
        if token in syntax.blank_tokens or (not token and not syntax.drop_empty):
            values.append(None)
            continue
        try:
            dtype = meta["dtype"]
            if dtype is bool:
                if token not in ("True", "False"):
                    raise ValueError("Use True or False (case-sensitive).")
                value = token == "True"
            elif dtype in (int, float):
                value = float(token)
                if not math.isfinite(value):
                    raise ValueError("Use a finite number.")
                if dtype is int:
                    if not value.is_integer():
                        raise ValueError("An integer is required.")
                    value = int(value)
            else:
                value = token
            if "min" in meta and value < meta["min"]:
                raise ValueError("Minimum: " + str(meta["min"]))
            if "max" in meta and value > meta["max"]:
                raise ValueError("Maximum: " + str(meta["max"]))
            allowed = meta.get("allowed")
            if isinstance(allowed, (list, tuple, set)) and value not in allowed:
                raise ValueError("Allowed values: " + display(allowed))
            values.append(value)
        except (ValueError, TypeError) as error:
            return ["{}: {}".format(token, error)]
    if required and field.category == "mandatory" and (not values or all(v is None for v in values)):
        return ["A value is required."]
    if required and field.category == "mandatory" and meta.get("is_list"):
        missing = [str(index + 1) for index, value in enumerate(values) if value is None]
        if missing:
            return ["Required values are missing at list positions {} (1-based).".format(", ".join(missing))]
    size = meta.get("size")
    sizes = [size] if isinstance(size, int) else size
    if meta.get("is_list") and sizes and len(values) not in sizes:
        return ["Expected list size: " + display(sizes)]
    return []


def comment_block(field, newline="\n"):
    text = field.flavor.get("description", "")
    lines = text.splitlines() + ["Category: " + field.category, "Default: " + display(field.meta.get("default"))]
    if field.flavor.get("example"):
        lines.append("Example: " + display(field.flavor["example"]))
    legacy = field.flavor.get("legacy", {})
    if isinstance(legacy, dict):
        lines.extend(name.replace("_", " ").capitalize() + ": " + display(value)
                     for name, value in legacy.items() if value)
    return "".join("# " + line + newline for line in lines)


class EditorModel:
    def __init__(self, schema):
        self.schema = schema
        self.suggested_save_path = None
        self.relink_on_save = False
        self.load_warning = ""
        self.path = None
        self.document = None
        self.values = {f.identity: "" for f in schema.fields}
        self.baseline = self.values.copy()
        self.entries = {}
        self.unknown = []

    def dirty(self):
        return self.values != self.baseline

    def load(self, path):
        document = IniDocument.load(path, inline_comments=False)
        if document.encoding not in ("utf-8", "utf-8-sig"):
            raise ValueError("ATLAS expects UTF-8. Convert this file to UTF-8 before loading it.")
        entries, unknown, section_names = {}, [], set()
        values = dict.fromkeys(self.values, "")
        for section in document.sections:
            if not section.implicit:
                if section.name in section_names:
                    raise ValueError("Duplicate section: " + section.name)
                section_names.add(section.name)
            for item in section.items:
                if not isinstance(item, Setting):
                    continue
                identity = (section.name, item.key)
                if identity in entries:
                    raise ValueError("Duplicate key: [{}] {}".format(*identity))
                entries[identity] = item
                if identity in values:
                    if item.multiline:
                        raise ValueError("Multiline value cannot be edited: [{}] {}".format(*identity))
                    field = next(f for f in self.schema.fields if f.identity == identity)
                    values[identity] = normalize_value(field, item.value)
                else:
                    unknown.append(identity)
        self.path, self.document = Path(path).resolve(), document
        self.entries, self.unknown = entries, unknown
        self.values = values
        self.baseline = values.copy()

    def validate(self, required=True):
        self.normalize()
        errors = []
        for item in self.schema.fields:
            for message in validate_value(item, self.values[item.identity], self.schema.filename, required):
                errors.append("[{}] {}: {}".format(item.section, item.key, message))
        errors += ["Unrecognized entry preserved in file: [{}] {}".format(*identity) for identity in self.unknown]
        return errors

    def normalize(self):
        for field in self.schema.fields:
            self.values[field.identity] = normalize_value(field, self.values[field.identity])

    def render(self, include_descriptions=False):
        self.normalize()
        sections = {}
        for section, fields in self.schema.sections.items():
            for field in fields:
                value = self.values[field.identity]
                if value != "":
                    comment = comment_block(field) if include_descriptions else ""
                    sections.setdefault(section, []).append(comment + field.key + " = " + value + "\n")
        # Keep supplied legacy/unknown values, but never silently flatten an
        # unsupported multiline entry or turn an implicit entry into a section.
        for identity in self.unknown:
            item = self.entries[identity]
            value = item.value.strip()
            if not value:
                continue
            if item.multiline or identity[0] == "No section":
                raise ValueError("Cannot compact an unrecognized multiline or unsectioned entry: " + item.key)
            if "\n" in value or "\r" in value:
                raise ValueError("Values must not contain line breaks.")
            sections.setdefault(identity[0], []).append(item.key + " = " + value + "\n")
        parts = []
        for section, entries in sections.items():
            parts.append("[{}]\n{}".format(section, "\n".join(entries) if include_descriptions else "".join(entries)))
        return "\n".join(parts).encode("utf-8")

    def save(self, path, include_descriptions=False):
        target = Path(path).resolve()
        if self.document is not None and target == self.path and target.read_bytes() != self.document.raw:
            raise ValueError("The file changed outside the GUI. Reload it or save under a different name.")
        data = self.render(include_descriptions=include_descriptions)
        atomic_write(target, data)
        self.load(target)


def load_working_folder(folder, schemas):
    """Read existing INIs first; create only missing files, without overwrites."""
    folder = Path(folder).expanduser().resolve()
    if not folder.is_dir():
        raise ValueError("Select an existing INI working folder.")
    call_path = folder / "call_atlas.ini"
    if call_path.exists() or call_path.is_symlink():
        return load_call_atlas_file(call_path, schemas)
    models, actions, missing = {}, [], []
    for filename, schema in schemas.items():
        model = EditorModel(schema)
        path = folder / filename
        if path.exists() or path.is_symlink():
            model.load(path)
            actions.append("Loaded " + filename)
        else:
            missing.append((filename, model, path))
        models[filename] = model
    for filename, model, path in missing:
        try:
            # Exclusive creation: never replace a file created concurrently.
            with path.open("xb"):
                pass
            actions.append("Created " + filename)
        except FileExistsError:
            actions.append("Loaded " + filename)
        model.load(path)
    return models, actions


def linked_ini_path(call_path, value, default_name):
    """Match runtime path rules: ~ and either separator; relative to init, not data."""
    value = value.strip()
    if not value:
        return call_path.parent / default_name
    if os.name != "nt":
        if ntpath.splitdrive(value)[0]:
            raise ValueError("Windows drive or network path cannot be resolved on this system: "
                             "{!r}. Use a local mount path or relative path.".format(value))
        value = value.replace("\\", "/")
    path = Path(value).expanduser()
    return Path(os.path.abspath(str(path if path.is_absolute() else call_path.parent / path)))


def load_call_atlas_file(path, schemas):
    """Keep the selected init editable even when its linked files cannot load.

    Unavailable configuration/settings become unsaved drafts, with a warning.
    Missing settings are created only at a valid resolved destination. Existing
    linked files are never overwritten while loading.
    """
    path = Path(path).expanduser().resolve()
    init = EditorModel(schemas["call_atlas.ini"])
    init.load(path)
    models = {"call_atlas.ini": init}
    actions, missing, paths = ["Loaded call_atlas.ini"], [], {path}

    def draft(filename, key, target, error):
        model = EditorModel(schemas[filename])
        suggested = path.parent / filename
        try:
            loaded_paths = {item.path for item in models.values() if item.path is not None}
            if (target is not None and target.parent.is_dir() and not target.is_dir()
                    and target.resolve() not in loaded_paths):
                suggested = target
        except (OSError, ValueError, RuntimeError):
            pass
        model.suggested_save_path = suggested
        model.relink_on_save = True
        model.load_warning = (
            "Could not load {}: {}\n{}\n"
            "Correct {} in {} and select the init again, or fill in the empty "
            "{} draft and save it. Saving the draft at a different path updates "
            "the init entry in the editor; save the init too before running."
        ).format(filename, target if target is not None else init.values[("explicit_paths", key)],
                 error, key, path.name, filename)
        actions.append("Prepared empty " + filename + " draft")
        models[filename] = model

    for filename, key in (("config_file.ini", "atlas_configuration_file"),
                          ("settings_file.ini", "atlas_settings_file")):
        raw = init.values[("explicit_paths", key)]
        target = None
        try:
            target = linked_ini_path(path, raw, filename)
            canonical = target.resolve()
            if canonical in paths:
                raise ValueError("Each INI tab must use a different file.")
            paths.add(canonical)
            model = EditorModel(schemas[filename])
            if target.exists() or target.is_symlink():
                model.load(target)
                actions.append("Loaded " + filename)
            elif filename == "config_file.ini":
                raise ValueError("Configuration file not found.")
            elif not target.parent.is_dir():
                raise ValueError("The INI destination folder does not exist.")
            else:
                missing.append((filename, key, model, target))
            models[filename] = model
        except (OSError, ValueError, RuntimeError) as error:
            draft(filename, key, target, error)
    for filename, key, model, target in missing:
        try:
            try:
                with target.open("xb"):
                    pass
                actions.append("Created " + filename)
            except FileExistsError:
                actions.append("Loaded " + filename)
            model.load(target)
        except (OSError, ValueError, RuntimeError) as error:
            draft(filename, key, target, error)
    return models, actions
