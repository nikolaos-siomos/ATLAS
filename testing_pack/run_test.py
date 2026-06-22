from __future__ import annotations

from pathlib import Path
import sys

from atlas_actris.testing.smoke_test import run_smoke_test


if __name__ == "__main__":
    pack_dir = Path(__file__).resolve().parent
    raise SystemExit(run_smoke_test(testing_pack=pack_dir))
