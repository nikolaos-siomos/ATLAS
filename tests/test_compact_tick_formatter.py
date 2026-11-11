"""Agg rendering regressions for fixed-layout polarization tick formatting."""
import ast
import sys
import unittest
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.ticker import FixedLocator, MultipleLocator

SOURCE = Path(__file__).resolve().parents[1]/"src"/"atlas_actris"
sys.path.insert(0, str(SOURCE))
from visualizer.compact_tick_formatter import CompactScalarFormatter
from visualizer.generate_polarization_calibration_axis import (
    get_x_ticks, get_calibration_y_limits, get_rayleigh_y_limits,
)


class CompactTickFormatterTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def axis(self, ticks, reference):
        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_axes([.2, .2, .7, .7])
        ax.set_ylim(min(ticks), max(ticks))
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        formatter = CompactScalarFormatter(reference)
        ax.yaxis.set_major_formatter(formatter)
        fig.canvas.draw()
        return fig, ax, formatter

    def test_middle_short_decimal_labels_remain_plain(self):
        _, _, f = self.axis([20.1, 20.2, 20.3], "0.275")
        self.assertEqual(f.orderOfMagnitude, 0)
        _, _, f = self.axis([.02, .04, .06], "0.275")
        self.assertEqual(f.orderOfMagnitude, 0)

    def test_middle_screenshot_labels_fit_reference_width(self):
        _, ax, f = self.axis(np.arange(.1, .276, .025), "0.275")
        self.assertEqual(f.orderOfMagnitude, 0)
        self.assertEqual(ax.yaxis.get_offset_text().get_text(), "")
        self.assertEqual(f.offset, 0)

    def test_decimal_point_and_minus_have_their_actual_widths(self):
        _, ax, f = self.axis([.1, .2, .3], "0.275")
        font = ax.yaxis.get_major_ticks(1)[0].label1.get_fontproperties()
        self.assertLess(f._label_width(".", font), f._label_width("0", font))
        self.assertGreater(f._label_width("−", font), f._label_width(".", font))
        self.assertEqual(f._label_width("20.1", font), f._label_width("0.02", font))
        self.assertLess(f._label_width("1.1", font), f._label_width("111", font))

    def test_middle_wider_decimal_labels_use_scientific_format(self):
        _, _, f = self.axis([.1001, .1002, .1003], "0.275")
        self.assertNotEqual(f.orderOfMagnitude, 0)

    def test_right_minus_0001_stays_plain(self):
        _, _, f = self.axis([-.001, 0, .001], "-0.001")
        self.assertEqual(f.orderOfMagnitude, 0)

    def test_right_long_labels_switch_to_shared_exponent(self):
        _, _, f = self.axis([-.00025, 0, .00025], "-0.001")
        self.assertNotEqual(f.orderOfMagnitude, 0)

    def test_minus_sign_counts_towards_budget(self):
        _, _, positive = self.axis([.0012, .0014, .0016], "-0.001")
        _, _, negative = self.axis([-.0016, -.0014, -.0012], "-0.001")
        self.assertEqual(positive.orderOfMagnitude, 0)
        self.assertNotEqual(negative.orderOfMagnitude, 0)

    def test_redraw_rechecks_changed_limits(self):
        fig, ax, f = self.axis([-.00025, 0, .00025], "-0.001")
        self.assertNotEqual(f.orderOfMagnitude, 0)
        ax.yaxis.set_major_locator(FixedLocator([0, .02, .04]))
        ax.set_ylim(0, .04)
        fig.canvas.draw()
        self.assertEqual(f.orderOfMagnitude, 0)
        self.assertEqual(ax.yaxis.get_offset_text().get_text(), "")

    def test_out_of_view_ticks_do_not_force_scientific_format(self):
        fig, ax, f = self.axis([.001, .002], "-0.001")
        ax.yaxis.set_major_locator(FixedLocator([-.00025, .001, .002]))
        fig.canvas.draw()
        self.assertEqual(f.orderOfMagnitude, 0)

    def test_formatting_preserves_limits_ticks_and_panel_geometry(self):
        fig, ax, f = self.axis([-.00025, 0, .00025], "-0.001")
        original_position = ax.get_position().bounds
        original_ticks = ax.get_yticks().copy()
        original_limits = ax.get_ylim()
        fig.canvas.draw()
        self.assertEqual(ax.get_position().bounds, original_position)
        self.assertEqual(ax.get_ylim(), original_limits)
        np.testing.assert_array_equal(ax.get_yticks(), original_ticks)
        self.assertEqual(f.offset, 0)

    def test_neighbouring_panel_reduces_allowed_tick_width(self):
        fig = plt.figure(figsize=(16, 4))
        fig.add_axes([.375, .14, .285, .58])
        ax = fig.add_axes([.705, .14, .285, .58])
        ax.set_ylabel("Linear dep. ratio")
        ax.set_ylim(0., .0117)
        ax.yaxis.set_major_locator(FixedLocator([0., .002, .004, .006, .008, .010]))
        formatter = CompactScalarFormatter("-0.001")
        ax.yaxis.set_major_formatter(formatter)
        original = ax.get_position().bounds
        fig.canvas.draw()
        self.assertNotEqual(formatter.orderOfMagnitude, 0)
        self.assertEqual(ax.get_position().bounds, original)

    def test_precision_is_preserved_for_narrow_ranges(self):
        _, ax, _ = self.axis([20.123, 20.124, 20.125], "0.275")
        labels = [label.get_text() for label in ax.get_yticklabels()]
        self.assertEqual(len(set(labels)), len(labels))

    def test_actual_three_panel_plot_keeps_original_positions(self):
        # Isolate the plotting functions from unrelated file-export dependencies.
        path = SOURCE/"visualizer"/"plot_polarization_calibration.py"
        tree = ast.parse(path.read_text())
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        scope = dict(np=np, plt=plt, MultipleLocator=MultipleLocator,
                     CompactScalarFormatter=CompactScalarFormatter,
                     get_x_ticks=get_x_ticks,
                     get_calibration_y_limits=get_calibration_y_limits,
                     get_rayleigh_y_limits=get_rayleigh_y_limits,
                     get_vertical_axis_label=lambda kind: "Range from lidar [km]",
                     add_fitting_suptitle=lambda fig, title, **kwargs: fig.suptitle(title),
                     export_plot=lambda fig, args: fig)
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), scope)
        x = np.linspace(0, 20, 401)
        signals = {key: np.exp(-x/2.) for key in
                   ("ray_r", "ray_t", "pcb_p45_r", "pcb_p45_t", "pcb_m45_r", "pcb_m45_t")}
        calibration = {key: .15 + .025*np.sin(x) for key in
                       ("gain_ratio", "gain_ratio_p45", "gain_ratio_m45")}
        rayleigh = {key: .0002*np.sin(x) for key in ("calibrated_ratio", "vldr", "mldr")}
        args = dict(title="Polarization calibration: fixed panel geometry", vertical_scale="range",
                    x_lims_signals=[0, 8], x_lims_calibration=[0, 8], x_lims_rayleigh=[0, 20],
                    y_lims_calibration=[.1001, .1009], y_lims_rayleigh=[-.00025, .00025],
                    calibration_region=[2, 3], rayleigh_region=[6, 7],
                    calibration_text={}, rayleigh_text={})
        fig = scope["generate_plot"]({key: x for key in signals}, signals, {}, x,
                                    calibration, {}, x, rayleigh, {}, args)
        fig.axes[1].yaxis.set_major_locator(FixedLocator(np.arange(.1001, .10091, .0001)))
        fig.canvas.draw()
        expected = [[.045, .14, .28, .58], [.375, .14, .285, .58], [.705, .14, .285, .58]]
        for ax, position in zip(fig.axes, expected):
            np.testing.assert_allclose(ax.get_position().bounds, position)
        self.assertEqual(fig.axes[1].yaxis.get_major_formatter().max_label, "0.275")
        self.assertEqual(fig.axes[2].yaxis.get_major_formatter().max_label, "-0.001")
        self.assertNotEqual(fig.axes[1].yaxis.get_major_formatter().orderOfMagnitude, 0)
        self.assertNotEqual(fig.axes[2].yaxis.get_major_formatter().orderOfMagnitude, 0)
        renderer = fig.canvas.get_renderer()
        for previous, current in zip(fig.axes[:-1], fig.axes[1:]):
            self.assertGreater(current.yaxis.label.get_window_extent(renderer).x0,
                               previous.get_window_extent(renderer).x1)
        # Reproduce the reported 532 and 355 nm tick ranges. Both fit the
        # reference string but need scientific notation to clear the y-title.
        middle, right = fig.axes[1:]
        middle.set_ylim(10., 43.)
        middle.yaxis.set_major_locator(FixedLocator([15., 20., 25., 30., 35., 40.]))
        for ticks, limits in (([0., .002, .004, .006, .008, .010], [0., .0117]),
                              (np.arange(0., .0176, .0025), [0., .019])):
            right.set_ylim(limits)
            right.yaxis.set_major_locator(FixedLocator(ticks))
            fig.canvas.draw()
            self.assertNotEqual(right.yaxis.get_major_formatter().orderOfMagnitude, 0)
            self.assertEqual(middle.yaxis.get_major_formatter().orderOfMagnitude, 0)
            for ax, position in zip(fig.axes, expected):
                np.testing.assert_allclose(ax.get_position().bounds, position)
            renderer = fig.canvas.get_renderer()
            self.assertGreater(right.yaxis.label.get_window_extent(renderer).x0,
                               middle.get_window_extent(renderer).x1 + 2.)
        # Optional local preview for visual inspection during development.
        import os
        if os.environ.get("ATLAS_FORMATTER_PREVIEW"):
            fig.savefig(os.environ["ATLAS_FORMATTER_PREVIEW"], dpi=150)


if __name__ == "__main__":
    unittest.main()
