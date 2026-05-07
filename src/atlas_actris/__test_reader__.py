#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 19:56:33 2025

@author: nikos
"""


#!/usr/bin/env python3
"""
Smoke-test the reader functions: call them and report exceptions.
Edit the imports and the TESTS list to match your project.
"""

import unittest
import traceback
import time
from pathlib import Path
from readers.read_scc import read_dataset as reader_scc
from readers.read_licel import read_dataset as reader_licel
from readers.read_polly_xt import read_dataset as reader_polly_xt
from readers.read_licel_matlab import read_dataset as reader_licel_matlab
from readers.read_polly_xt_first import read_dataset as reader_polly_xt_first


PARAMS = [
    ("reader_scc", reader_scc, ("./sample_files/scc", "ray")),
    ("reader_licel", reader_licel, ("./sample_files/licel", "ray")),
    ("reader_polly_xt", reader_polly_xt, ("./sample_files/polly_xt", "ray")),
    ("reader_licel_matlab", reader_licel_matlab, ("./sample_files/licel_matlab", "ray")),
    ("reader_polly_xt_first", reader_polly_xt_first, ("./sample_files/polly_xt_first", "ray")),
    ("reader_licel_old2rack", reader_licel, ("./sample_files/licel_old2rack", "ray")),
]


def _missing_first_path(args) -> bool:
    """Return True if first arg looks like a filesystem path and doesn't exist."""
    if not args:
        return False
    first = args[0]
    if isinstance(first, (str, Path)):
        p = Path(first)
        if ("/" in str(p) or "\\" in str(p)) and not p.exists():
            return True
    return False


class TestReaders(unittest.TestCase):
    # store (status, name, seconds, info)
    results = []

    @classmethod
    def setUpClass(cls):
        print("=== Reader smoke tests ===")

    def test_readers_smoke(self):
        for name, func, args in PARAMS:
            with self.subTest(reader=name):
                if _missing_first_path(args):
                    msg = f"Missing path: {args[0]}"
                    print(f"[SKIP] {name:>24s}  {msg}")
                    # Record and mark subtest as skipped in unittest
                    self.results.append(("SKIP", name, 0.0, msg))
                    self.skipTest(msg)

                t0 = time.perf_counter()
                try:
                    print('')
                    print(f"Checking reader: {name}")
                    func(*args)  # just call; pass if no exception
                    dt = time.perf_counter() - t0
                    print(f"[PASS] {name:>24s}  ({dt:.3f}s)")
                    self.results.append(("PASS", name, dt, ""))
                except Exception as e:
                    dt = time.perf_counter() - t0
                    tb = traceback.format_exc()
                    print(f"[FAIL] {name:>24s}  ({dt:.3f}s)")
                    print(f"{e.__class__.__name__}: {e}\n{tb}")
                    self.results.append(("FAIL", name, dt, f"{e.__class__.__name__}: {e}\n{tb}"))
                    # Make unittest record the failure (but continue to next subTest)
                    self.fail(f"{name} raised {e.__class__.__name__}: {e}\n{tb}")

    @classmethod
    def tearDownClass(cls):
        # Print a compact summary at the end
        passed = sum(1 for s, *_ in cls.results if s == "PASS")
        failed = sum(1 for s, *_ in cls.results if s == "FAIL")
        skipped = sum(1 for s, *_ in cls.results if s == "SKIP")
        total = len(cls.results)
        print(f"\nSummary: {passed} passed, {failed} failed, {skipped} skipped (out of {total})")


if __name__ == "__main__":
    unittest.main(verbosity=2)