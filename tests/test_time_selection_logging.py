"""Tests of interval reporting and its wiring into the actual slicing function.

Signal-trimming functions are loaded in isolation to keep these coordinate-only
checks independent of xarray/Dask and the full processing installation.
"""
import ast
import contextlib
import io
import re
import sys
import unittest
import copy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

SOURCE = Path(__file__).resolve().parents[1]/"src"/"atlas_actris"
sys.path.insert(0, str(SOURCE))
from processor.time_selection_logging import print_time_selections
from utils.parse_init_file import _expand_slice_exclude_bundles


class Coordinate:
    def __init__(self, values):
        self.values = np.asarray(values)
        self.size = self.values.size
        self.data = self.values
    def __ge__(self, value): return Coordinate(self.values >= np.datetime64(value))
    def __le__(self, value): return Coordinate(self.values <= np.datetime64(value))
    def __gt__(self, value): return Coordinate(self.values > np.datetime64(value))
    def __lt__(self, value): return Coordinate(self.values < np.datetime64(value))
    def __or__(self, other): return Coordinate(self.values | other.values)
    def __and__(self, other): return Coordinate(self.values & other.values)
    def __invert__(self): return Coordinate(~self.values)


class TimeInfo:
    def __init__(self, times): self.time = Coordinate(np.array(times, dtype="datetime64[s]"))
    def __getitem__(self, key): return self.time
    def sel(self, key): return SimpleNamespace(values=self.time.values.astype(str))


def isolated_trimming_functions():
    names = {"hhmm_to_datetime", "_parse_slice_time", "make_interval",
             "_resolve_loading_map_key", "compute_slice_and_exclude"}
    tree = ast.parse((SOURCE/"processor"/"signal_trimming.py").read_text())
    funcs = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    module = ast.Module(body=[ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0)]+funcs, type_ignores=[])
    namespace = dict(np=np, pd=pd, re=re, copy=copy,
                     xr=SimpleNamespace(
                         ones_like=lambda arr, dtype: Coordinate(np.ones(arr.size, dtype=dtype)),
                         zeros_like=lambda arr, dtype: Coordinate(np.zeros(arr.size, dtype=dtype))),
                     print_time_selections=print_time_selections,
                     iso_to_datetimes=lambda values: pd.to_datetime(values),
                     CustomWarning=lambda text: print(text),
                     endpoint=lambda ID: (_ for _ in ()).throw(SystemExit(ID)))
    def capture_selection(keys, masks, output):
        output["captured_masks"] = {key: value.values for key, value in masks.items()}
        return output
    namespace["select_by_time_mask"] = capture_selection
    exec(compile(ast.fix_missing_locations(module), str(SOURCE/"processor"/"signal_trimming.py"), "exec"), namespace)
    return namespace


class TimeSelectionLoggingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.functions = isolated_trimming_functions()

    def run_rules(self, slices, exclusions, measurements, logs=None, loading_map=None):
        caller = dict(slice_measurement=slices, exclude_measurement=exclusions,
                      loading_map=loading_map or {})
        if logs is not None:
            caller["_time_selection_requests"] = logs
        data = {"time_info": {key: TimeInfo(times) for key, times in measurements.items()},
                "profile": {key: object() for key in measurements}}
        with contextlib.redirect_stdout(io.StringIO()) as output:
            result = self.functions["compute_slice_and_exclude"]({"caller_info": caller}, data)
        return result, output.getvalue().splitlines()

    def test_three_requested_rules_produce_three_rows_not_six(self):
        slices_log, exclude_log = [], []
        slices = _expand_slice_exclude_bundles(
            ["ray", "2000", "2100", "tlc", "2100", "2200"], request_log=slices_log)
        exclusions = _expand_slice_exclude_bundles(
            ["ray", "2010", "2020"], request_log=exclude_log)
        data = {"ray": ["2025-10-29T20:00", "2025-10-29T20:10", "2025-10-29T20:20", "2025-10-29T21:00"]}
        data.update({key: ["2025-10-29T21:00", "2025-10-29T22:00"]
                     for key in ("tlc_north", "tlc_south", "tlc_east", "tlc_west")})
        result, lines = self.run_rules(slices, exclusions, data,
                    {"slice_measurement": slices_log, "exclude_measurement": exclude_log})
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("-- Slice ray:"))
        self.assertTrue(lines[1].startswith("-- Slice tlc:"))
        self.assertIn("8/8 profiles", lines[1])
        self.assertIn("2 profiles removed", lines[2])
        self.assertEqual(result["captured_masks"]["ray"].tolist(), [True, False, False, True])

    def test_2245_to_2240_prints_next_day_and_no_change(self):
        result, lines = self.run_rules(["tlc_north", "2245", "2240"], [],
                   {"tlc_north": ["2025-10-29T22:45:11", "2025-10-29T23:06:20"]})
        self.assertEqual(len(lines), 1)
        self.assertIn("2025-10-29 22:45 -> 2025-10-30 22:40 UTC", lines[0])
        self.assertIn("2/2 profiles in interval (unchanged)", lines[0])
        self.assertTrue(result["captured_masks"]["tlc_north"].all())

    def test_repeated_slice_intervals_keep_union_and_report_each(self):
        result, lines = self.run_rules(["ray", "2000", "2010", "ray", "2030", "2040"], [],
                     {"ray": ["2025-10-29T20:00", "2025-10-29T20:20", "2025-10-29T20:40"]})
        self.assertEqual(len(lines), 2)
        self.assertEqual(result["captured_masks"]["ray"].tolist(), [True, False, True])

    def test_exclusion_counts_only_new_removals_after_slicing(self):
        result, lines = self.run_rules(["ray", "2000", "2030"],
                  ["ray", "2010", "2020", "ray", "2010", "2020", "ray", "2100", "2200"],
                  {"ray": ["2025-10-29T20:00", "2025-10-29T20:10", "2025-10-29T21:30"]})
        self.assertEqual(len(lines), 4)
        self.assertIn("1 profile removed", lines[1])
        self.assertIn("0 profiles removed (unchanged)", lines[2])
        self.assertIn("0 profiles removed (unchanged)", lines[3])

    def test_missing_measurement_has_one_short_skipped_row(self):
        _, lines = self.run_rules(["pcb_p45", "2000", "2100"], [], {"ray": ["2025-10-29T20:30"]})
        self.assertEqual(len(lines), 1)
        self.assertIn("skipped (no data)", lines[0])

    def test_explicit_sector_overrides_bundle_in_expansion_and_logging(self):
        log = []
        expanded = _expand_slice_exclude_bundles(
            ["tlc_rin", "2100", "2200", "tlc_inner", "2130", "2140"], request_log=log)
        self.assertEqual(log[0]["members"], ["tlc_outer"])
        _, lines = self.run_rules(expanded, [],
                  {key: ["2025-10-29T21:30"] for key in ("tlc_inner", "tlc_outer")},
                  {"slice_measurement": log, "exclude_measurement": []})
        self.assertEqual(len(lines), 2)
        self.assertIn("21:00 ->", lines[0])
        self.assertIn("21:30 ->", lines[1])

    def test_alias_source_is_shown(self):
        _, lines = self.run_rules(["drk_ray", "2000", "2100"], [],
                  {"drk": ["2025-10-29T20:30"]}, loading_map={"drk_ray": "drk"})
        self.assertIn("Slice drk_ray (drk):", lines[0])

    def test_explicit_dates_and_seconds_are_shown(self):
        _, lines = self.run_rules(["ray", "20251029_224511", "20251029_230620"], [],
                    {"ray": ["2025-10-29T22:45:11", "2025-10-29T23:06:20"]})
        self.assertIn("22:45:11 -> 2025-10-29 23:06:20", lines[0])

    def test_no_rules_produce_no_selection_messages(self):
        _, lines = self.run_rules([], [], {"ray": ["2025-10-29T20:30"]})
        self.assertEqual(lines, [])

    def test_distinct_sector_dates_are_not_hidden_by_grouping(self):
        log = []
        expanded = _expand_slice_exclude_bundles(["tlc_rin", "2245", "2240"], request_log=log)
        _, lines = self.run_rules(expanded, [], {"tlc_inner": ["2025-10-29T23:00"],
                    "tlc_outer": ["2025-10-30T00:30"]},
                    {"slice_measurement": log, "exclude_measurement": []})
        self.assertEqual(len(lines), 1)
        self.assertIn("tlc_inner: 2025-10-29 22:45", lines[0])
        self.assertIn("tlc_outer: 2025-10-30 22:45", lines[0])

    def test_stale_parser_log_falls_back_to_actual_runtime_rules(self):
        _, lines = self.run_rules(["ray", "2000", "2100"], [],
                  {"ray": ["2025-10-29T20:30"]},
                  {"slice_measurement": [dict(key="tlc", start="2200", stop="2300", members=["tlc_north"])]})
        self.assertTrue(lines[0].startswith("-- Slice ray:"))


if __name__ == "__main__":
    unittest.main()
