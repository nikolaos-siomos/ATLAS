# `atlas-generate-templates`

Use `atlas-generate-templates` in the environment where ATLAS is installed. To write templates into a folder of your choice, run:

```bash
atlas-generate-templates -o ./my_templates
```

The folder is created automatically. Relative paths are resolved from your current working directory. You do not need to locate the ATLAS project folder. With `-o` and no explicit `--target`, only INI files are generated.

Generate compact bare templates:

```bash
atlas-generate-templates -o ./my_templates --profile bare
```

The equivalent long output option is `--output_folder`. Quote paths containing spaces:

```powershell
atlas-generate-templates -o "C:\ATLAS runs\templates" --profile bare
```

## Profiles and filenames

| Profile | Files generated |
| --- | --- |
| `full` | `call_atlas.ini`, `config_file.ini`, `settings_file.ini`, `intercomparison.ini` |
| `bare` | `call_atlas_bare.ini`, `config_file_bare.ini`, `settings_file_bare.ini`, `intercomparison_bare.ini` |
| `beginner` | `call_atlas_beginner.ini`, `config_file_beginner.ini`, `settings_file_beginner.ini` |
| `all` (default) | All the files above. |

Full templates contain parameter comments. Bare templates omit per-parameter comments and blank lines between variables, keeping one blank line between sections. Beginner templates provide the selected beginner parameters and their comments. There is no separate beginner intercomparison template.

Generation overwrites files with matching names in the destination folder. Use a separate folder for templates if you have already edited configuration files there.

## Output and target options

| Option | Behavior |
| --- | --- |
| `-o`, `--output_folder` | Optional destination for INI files. Without it, INI files go to the repository's `templates/` directory. |
| `--profile all/full/bare/beginner` | Select INI template profiles; default `all`. |
| `--target ini/docs/all` | Select output types. Default `ini` with `-o`, otherwise `all`. |
| `--repo-root` | Override the repository root used for default paths and generated documentation. |
| `--check` | Compare expected files with existing files without writing; fail if missing or different. |

For the existing developer workflow, run from the repository:

```bash
atlas-generate-templates
```

This generates INI files under `templates/` and reference pages under `docs/generated/`.

Generate only documentation:

```bash
atlas-generate-templates --target docs
```

`-o` controls only INI files, so combining it with `--target docs` is rejected. With `-o ./my_templates --target all`, INI files go into `./my_templates`, while documentation still goes into the repository's `docs/generated/` directory.

Check a custom bare-template folder without changing it:

```bash
atlas-generate-templates -o ./my_templates --profile bare --check
```

Check the repository outputs:

```bash
atlas-generate-templates --check
```

After changing parser schemas or template metadata, rerun the generator. Editable reinstallation is only required when console-script definitions change.
