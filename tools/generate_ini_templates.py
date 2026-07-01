#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Developer wrapper for ATLAS template/documentation generation."""

from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from atlas_actris.templates.generator import main


if __name__ == "__main__":
    raise SystemExit(main())
