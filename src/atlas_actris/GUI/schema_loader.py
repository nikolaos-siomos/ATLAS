"""Read ATLAS schemas and help text from source without importing ATLAS."""
import ast
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path


# GUI choices are intentionally independent snapshots of today's beginner
# template selections. Editing template_profiles.py must not change these.
RECOMMENDED_FIELDS = {
    "call_atlas.ini": frozenset({
        ("trimming_options", "slice_measurement"),
        ("trimming_options", "exclude_measurement"),
        ("parsing_options", "rsonde_station_name"),
        ("parsing_options", "rsonde_station_wmo_id"),
        ("parsing_options", "cloudnet_station_name"),
    }),
    "settings_file.ini": frozenset({
        ("quadrant_telecover", "normalization_region"),
        ("ring_telecover", "normalization_region"),
        ("polarization_calibration", "calibration_region"),
        ("polarization_calibration", "rayleigh_region"),
    }),
}


@dataclass(frozen=True)
class Unresolved:
    expression: str


class SourceError(ValueError):
    pass


class SourceReader:
    """Small evaluator for declarations and simple text helper functions.

    Imports and application code are never executed. Unsupported expressions
    remain explicit Unresolved values; structural schema errors stop loading.
    """
    def __init__(self, path, seeds=None):
        self.path = Path(path)
        tree = ast.parse(self.path.read_text(encoding="utf-8-sig"), filename=str(path))
        self.env = dict(seeds or {})
        self.definitions = {}
        self.functions = {}
        self.updates = []
        self.active = set()
        for node in tree.body:
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for target in targets:
                    if isinstance(target, ast.Name) and node.value is not None:
                        self.definitions[target.id] = node.value
                    elif isinstance(target, ast.Subscript):
                        self.updates.append((target, node.value))
            elif isinstance(node, ast.FunctionDef):
                self.functions[node.name] = node

    def get(self, name):
        if name in self.env:
            return self.env[name]
        if name not in self.definitions or name in self.active:
            return Unresolved(name)
        self.active.add(name)
        try:
            value = self.evaluate(self.definitions[name], {})
            self.env[name] = value
            for target, expression in self.updates:
                if isinstance(target.value, ast.Name) and target.value.id == name:
                    if isinstance(value, dict):
                        value[self.evaluate(target.slice, {})] = self.evaluate(expression, {})
            return value
        finally:
            self.active.remove(name)

    def evaluate(self, node, local):
        def ev(n):
            return self.evaluate(n, local)
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id in local:
                return local[node.id]
            if node.id in ("str", "int", "float", "bool"):
                return {"str": str, "int": int, "float": float, "bool": bool}[node.id]
            return self.get(node.id)
        if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
            items = [ev(n) for n in node.elts]
            return {ast.List: list, ast.Tuple: tuple, ast.Set: set}[type(node)](items)
        if isinstance(node, ast.Dict):
            result = {}
            for key, value in zip(node.keys, node.values):
                if key is None:
                    extra = ev(value)
                    if not isinstance(extra, dict):
                        return Unresolved(ast.unparse(node))
                    result.update(extra)
                else:
                    result[ev(key)] = ev(value)
            return result
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            value = ev(node.operand)
            if isinstance(value, (int, float)):
                return -value if isinstance(node.op, ast.USub) else value
        if isinstance(node, ast.BinOp):
            left, right = ev(node.left), ev(node.right)
            if not isinstance(left, Unresolved) and not isinstance(right, Unresolved):
                if isinstance(node.op, ast.Add) and isinstance(left, (str, list, tuple, int, float)) and isinstance(right, type(left)):
                    return left + right
                if isinstance(node.op, ast.BitOr) and isinstance(left, set) and isinstance(right, set):
                    return left | right
        if isinstance(node, ast.Subscript):
            value, key = ev(node.value), ev(node.slice)
            if isinstance(value, (dict, list, tuple)) and not isinstance(key, Unresolved):
                try:
                    return value[key]
                except (KeyError, IndexError, TypeError):
                    pass
        if isinstance(node, ast.JoinedStr):
            parts = []
            for part in node.values:
                value = ev(part.value) if isinstance(part, ast.FormattedValue) else ev(part)
                if isinstance(value, Unresolved):
                    return value
                parts.append(str(value))
            return "".join(parts)
        if isinstance(node, ast.Call):
            args = [ev(n) for n in node.args]
            if any(isinstance(a, Unresolved) for a in args) or node.keywords:
                return Unresolved(ast.unparse(node))
            if isinstance(node.func, ast.Attribute) and node.func.attr == "keys" and not args:
                value = ev(node.func.value)
                if isinstance(value, dict):
                    return list(value.keys())
            if isinstance(node.func, ast.Name):
                name = node.func.id
                if name in ("sorted", "set", "list", "tuple", "str"):
                    return {"sorted": sorted, "set": set, "list": list, "tuple": tuple, "str": str}[name](*args)
                function = self.functions.get(name)
                if function is not None and len(function.args.args) == len(args):
                    # Only straight-line text helpers; loops, I/O and arbitrary
                    # expression statements are unsupported.
                    scope = dict(zip((a.arg for a in function.args.args), args))
                    for statement in function.body:
                        if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Constant):
                            continue
                        if isinstance(statement, ast.Assign) and all(isinstance(t, ast.Name) for t in statement.targets):
                            for target in statement.targets:
                                scope[target.id] = self.evaluate(statement.value, scope)
                        elif isinstance(statement, ast.Return):
                            return self.evaluate(statement.value, scope)
                        else:
                            break
        return Unresolved(ast.unparse(node))


def display(value):
    if isinstance(value, Unresolved):
        return "Unavailable from supplied files: " + value.expression
    if value is None or value == [] or value == "":
        return "automatic / not specified"
    if isinstance(value, (list, tuple)):
        return ", ".join("None" if v is None else str(v) for v in value)
    return str(value)


@dataclass(frozen=True)
class ValueSyntax:
    separator: str = ","
    drop_empty: bool = True
    blank_tokens: tuple = ()


def value_syntax(parser):
    """Recognize the supplied list formats without executing parser functions."""
    function = parser.functions.get("_split_list")
    if function is None or not isinstance(function.body[-1], ast.Return):
        raise SourceError("Cannot resolve list syntax in " + str(parser.path))
    expression = ast.dump(function.body[-1].value, include_attributes=False)
    blank_tokens = parser.get("blank_tokens") if "blank_tokens" in parser.definitions else []
    if not isinstance(blank_tokens, (list, tuple)) or not all(isinstance(t, str) for t in blank_tokens):
        raise SourceError("Cannot resolve blank tokens in " + str(parser.path))
    # Both supplied parser generations use these literal forms. Unrecognized
    # changes fail explicitly instead of applying stale validation rules.
    for separator in (",", ", "):
        for drop_empty in (True, False):
            candidate = '[x.strip() for x in s.replace(";", {0!r}).split({0!r}){1}]'.format(
                separator, ' if x.strip() != ""' if drop_empty else '')
            if expression == ast.dump(ast.parse(candidate, mode="eval").body, include_attributes=False):
                return ValueSyntax(separator, drop_empty, tuple(blank_tokens))
    raise SourceError("Unsupported list syntax in " + str(parser.path))


@dataclass
class Field:
    section: str
    key: str
    meta: dict
    flavor: dict
    source: str
    value_syntax: ValueSyntax = ValueSyntax()
    gui_category: str = ""

    @property
    def identity(self):
        return self.section, self.key

    @property
    def category(self):
        return self.gui_category or self.meta.get("category", "optional")

    def visible(self, recommended=True, optional=True):
        return self.category == "mandatory" or (self.category == "recommended" and recommended) or (self.category == "optional" and optional)

    def help_text(self):
        lines = [self.flavor.get("description", "No description provided in the flavor file."), "",
                 "Category: " + self.category,
                 "Type: " + self.meta["dtype"].__name__ + (" list" if self.meta.get("is_list") else ""),
                 "Default: " + display(self.meta.get("default"))]
        if self.gui_category:
            lines.append("Category assigned by the GUI's independent recommended-field list.")
        elif "category" not in self.meta:
            lines.append("Category is not declared; this GUI treats the setting as optional.")
        for key, label in (("allowed", "Allowed values"), ("min", "Minimum"), ("max", "Maximum"), ("size", "List size"), ("check_path", "Path type")):
            if key in self.meta:
                lines.append(label + ": " + display(self.meta[key]))
        if self.flavor.get("example"):
            lines.extend(["", "Example: " + display(self.flavor["example"])])
        legacy = self.flavor.get("legacy", {})
        if isinstance(legacy, dict):
            details = [name.replace("_", " ").capitalize() + ": " + display(value)
                       for name, value in legacy.items() if value]
            if details:
                lines.extend(["", "Version notes:"] + details)
        lines.extend(["", "Schema source: " + self.source])
        return "\n".join(str(line) for line in lines)


@dataclass
class FileSchema:
    filename: str
    sections: dict

    @property
    def fields(self):
        return [item for fields in self.sections.values() for item in fields]


def source_root(folder):
    folder = Path(folder).expanduser().resolve()
    for candidate in (folder, folder / "src" / "atlas_actris", folder / "atlas_actris"):
        if (candidate / "utils" / "parse_init_file.py").is_file():
            return candidate
    # Also accept a flat folder containing the six uploaded files.
    if (folder / "parse_init_file.py").is_file():
        return folder
    raise SourceError("Select the ATLAS project folder, src/atlas_actris, or a folder containing the six source files.")


def load_schemas(folder):
    root = source_root(folder)
    def reader(subdir, name, seeds=None):
        path = root / subdir / name if (root / subdir).is_dir() else root / name
        if not path.is_file():
            raise SourceError("Missing source file: " + str(path))
        return SourceReader(path, seeds)
    init = reader("utils", "parse_init_file.py")
    config = reader("utils", "parse_config_file.py")
    settings = reader("utils", "parse_settings_file.py")
    init_flavor = reader("templates", "init_template_flavor.py", {
        name: init.get(name) for name in ("INIT_FILE_SECTIONS", "slice_exclude_allowed_keys", "version_warning_nrm", "version_warning_pcb")})
    config_flavor = reader("templates", "config_template_flavor.py")
    settings_flavor = reader("templates", "settings_template_flavor.py")
    schemas = OrderedDict()

    def assemble(filename, layout, technical, flavor, parser):
        sections = OrderedDict()
        syntax = value_syntax(parser)
        if not isinstance(layout, dict):
            raise SourceError("Cannot resolve the section layout in " + str(parser.path))
        for section, keys in layout.items():
            items = []
            for key in keys:
                meta = technical.get((section, key))
                if not isinstance(meta, dict) or meta.get("dtype") not in (str, int, float, bool):
                    raise SourceError("Unsupported field schema: [{}] {}".format(section, key))
                for name, value in meta.items():
                    if isinstance(value, Unresolved) and name != "allowed":
                        raise SourceError("Cannot resolve {} for {}: {}".format(name, key, value.expression))
                if meta.get("category", "optional") not in ("mandatory", "recommended", "optional"):
                    raise SourceError("Unknown category for " + key)
                help_data = flavor.get((section, key), {})
                if not isinstance(help_data, dict) or any(isinstance(v, Unresolved) for v in help_data.values()):
                    raise SourceError("Cannot resolve help text for " + key)
                category = ("recommended" if meta.get("category") != "mandatory"
                            and (section, key) in RECOMMENDED_FIELDS.get(filename, ()) else "")
                items.append(Field(section, key, meta, help_data, str(parser.path), syntax, category))
            sections[section] = items
        schemas[filename] = FileSchema(filename, sections)

    init_meta, init_layout, init_text = init.get("SCHEMA"), init.get("INIT_FILE_SECTIONS"), init_flavor.get("INIT_FLAVOR")
    if not all(isinstance(v, dict) for v in (init_meta, init_layout, init_text)):
        raise SourceError("Cannot resolve initialization schema.")
    if {k for keys in init_layout.values() for k in keys} != set(init_meta):
        raise SourceError("Initialization section layout does not match SCHEMA.")
    assemble("call_atlas.ini", init_layout, {(s, k): init_meta[k] for s, keys in init_layout.items() for k in keys},
             {(s, k): init_text.get(k, {}) for s, keys in init_layout.items() for k in keys}, init)

    config_meta, config_text = config.get("SCHEMA"), config_flavor.get("CONFIG_FLAVOR")
    if not all(isinstance(value, dict) for value in (config_meta, config_text)):
        raise SourceError("Cannot resolve configuration schema or help text.")
    # Recent flavor files only provide descriptions. The parser supplies the
    # section membership and SCHEMA insertion order; older layouts still work.
    template_layout = (config_flavor.get("CONFIG_TEMPLATE_SECTIONS")
                       if "CONFIG_TEMPLATE_SECTIONS" in config_flavor.definitions else {})
    if not isinstance(template_layout, dict) or any(
            not isinstance(keys, (list, tuple)) for keys in template_layout.values()):
        raise SourceError("Unsupported configuration template layout.")
    groups = ("SYSTEM_KEYS", "CHANNEL_KEYS", "POL_CAL_KEYS", "GL_KEYS", "WV_KEYS", "TEMP_KEYS")
    section_names = config.get("recognized_sections")
    if not isinstance(section_names, list) or len(section_names) != len(groups):
        raise SourceError("Unsupported configuration section layout.")
    config_layout = OrderedDict()
    for section, group in zip(section_names, groups):
        members = config.get(group)
        if not isinstance(members, set):
            raise SourceError("Cannot resolve configuration group " + group)
        ordered = [k for k in template_layout.get(section, []) if k in members]
        ordered += [k for k in config_meta if k in members and k not in ordered]
        config_layout[section] = ordered
    if {k for keys in config_layout.values() for k in keys} != set(config_meta):
        raise SourceError("Configuration groups do not match SCHEMA.")
    assemble("config_file.ini", config_layout, {(s, k): config_meta[k] for s, keys in config_layout.items() for k in keys},
             {(s, k): config_text.get(k, {}) for s, keys in config_layout.items() for k in keys}, config)

    settings_meta, settings_text = settings.get("SCHEMA"), settings_flavor.get("SETTINGS_FLAVOR")
    names = settings.get("recognized_sections")
    layout = OrderedDict((section, list(settings_meta[alias])) for alias, section in names.items())
    assemble("settings_file.ini", layout,
             {(section, k): meta for alias, section in names.items() for k, meta in settings_meta[alias].items()},
             {(section, k): settings_text.get(alias, {}).get(k, {}) for alias, section in names.items() for k in settings_meta[alias]}, settings)
    return schemas
