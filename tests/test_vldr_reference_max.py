"""NumPy regressions for automatic VLDR scaling, without full lidar I/O."""
import ast
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

SOURCE = Path(__file__).resolve().parents[1]/"src"/"atlas_actris"
sys.path.insert(0, str(SOURCE))
from utils.parse_settings_file import SCHEMA_QCK_VLDR


def load_functions(path, names, namespace):
    tree = ast.parse(path.read_text())
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)


class VLDRReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scope = {"np": np}
        load_functions(SOURCE/"visualizer"/"smoothing.py",
                       {"quicklook_normalization_max", "_as_odd_window", "_coord_step", "_limits_to_slice"}, cls.scope)
        load_functions(SOURCE/"visualizer"/"plot_vldr.py",
                       {"vldr_reference_max", "_manual_y_lims", "_safe_positive_vmin", "_safe_vldr_color_limits"}, cls.scope)

    def setUp(self):
        self.x = np.arange(101)*.1
        self.y = np.full((2, 101), .05)
        self.y[0, 20] = 1.05
        self.y[:, 80] = 1000.  # Outside y_max_zone: must not set the scale.
        self.settings = dict(y_lims=[], y_max_zone=[1., 3.], smoothing_window=.5,
                             smoothing_range=[0., 10.], smooth=False)

    def maximum(self, values=None, settings=None):
        return self.scope["vldr_reference_max"](self.y if values is None else values,
                                                self.x, self.settings if settings is None else settings)

    def limits(self, values=None, settings=None):
        return self.scope["_safe_vldr_color_limits"](self.y if values is None else values,
                                                     self.x, self.settings if settings is None else settings)

    def test_default_schema_is_automatic_again(self):
        self.assertEqual(SCHEMA_QCK_VLDR["y_lims"]["default"], [])

    def test_spike_is_spatially_smoothed_for_reference(self):
        self.assertAlmostEqual(self.maximum(), .25)
        self.assertEqual(self.y[0, 20], 1.05)
        self.assertAlmostEqual(self.limits()[1], .25)
        self.assertEqual(self.limits()[0], 0.)

    def test_reference_does_not_depend_on_display_smoothing_flag(self):
        smoothed = dict(self.settings, smooth=True)
        self.assertEqual(self.maximum(settings=smoothed), self.maximum())

    def test_variable_window_uses_the_shared_normal_quicklook_reference(self):
        settings = dict(self.settings, smoothing_window=[.3, .7])
        expected = self.scope["quicklook_normalization_max"](self.y, self.x, settings)
        self.assertEqual(self.maximum(settings=settings), expected)
        self.assertLess(expected, 1.05)

    def test_explicit_color_limits_override_automatic_reference(self):
        self.assertEqual(self.limits(settings=dict(self.settings, y_lims=[-.1, .5])), (-.1, .5))

    def test_cached_raw_reference_is_reused_without_resmoothing_display(self):
        settings = dict(self.settings, _vldr_reference_max=.25)
        with patch.dict(self.scope, {"vldr_reference_max": lambda *args: self.fail("Reference recomputed")}):
            self.assertEqual(self.limits(values=np.full_like(self.y, 99.), settings=settings), (0., .25))

    def test_invalid_cached_reference_goes_directly_to_one(self):
        for value in (0., -1., np.nan, np.inf, None, "invalid"):
            with self.subTest(value=value):
                self.assertEqual(self.limits(settings=dict(self.settings, _vldr_reference_max=value)), (0., 1.))

    def test_no_full_data_percentile_fallback(self):
        data = np.full_like(self.y, np.nan)
        data[:, 80:] = 1000.  # Large finite data outside reference zone.
        self.assertEqual(self.limits(values=data), (0., 1.))

    def test_nan_zero_and_negative_references_fall_back_to_one(self):
        for value in (np.nan, 0., -.1):
            with self.subTest(value=value):
                self.assertEqual(self.limits(values=np.full_like(self.y, value)), (0., 1.))

    def test_out_of_bounds_zone_falls_back_to_one(self):
        self.assertEqual(self.limits(settings=dict(self.settings, y_max_zone=[20., 30.])), (0., 1.))

    def test_invalid_window_or_coordinates_fall_back_to_one(self):
        self.assertEqual(self.limits(settings=dict(self.settings, smoothing_window=None)), (0., 1.))
        self.assertEqual(self.scope["vldr_reference_max"](np.ones((1,1)), np.array([1.]), self.settings), 1.)

    def test_log_scaling_keeps_positive_lower_limit(self):
        lower, upper = self.limits(settings=dict(self.settings, use_log_y_scale=True))
        self.assertGreater(lower, 0.)
        self.assertLess(lower, upper)
        self.assertAlmostEqual(upper, .25)


if __name__ == "__main__":
    unittest.main()
