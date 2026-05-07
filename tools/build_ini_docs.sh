#!/usr/bin/env bash
set -euo pipefail

python tools/gen_ini_docs.py \
  --ini call_atlas.ini --out docs/config/initialization.md --title "Initialization" \
  --ini config_file.ini --out docs/config/configuration.md --title "Configuration" \
  --ini settings_file.ini --out docs/config/settings.md --title "Settings"
