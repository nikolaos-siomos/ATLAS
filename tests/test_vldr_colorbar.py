"""Render VLDR colorbars without importing unrelated data-processing modules."""
import ast
import sys
import unittest
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.ticker import FixedLocator, LogFormatterSciNotation

SOURCE = Path(__file__).resolve().parents[1]/"src"/"atlas_actris"
sys.path.insert(0,str(SOURCE))
from visualizer.compact_tick_formatter import CompactScalarFormatter
from utils.parse_settings_file import SCHEMA_QCK_VLDR


class VLDRColorbarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = SOURCE/"visualizer"/"plot_vldr.py"
        tree = ast.parse(path.read_text())
        function = next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=="add_colorbar")
        cls.scope = dict(CompactScalarFormatter=CompactScalarFormatter, LogNorm=LogNorm,
                         LogFormatterSciNotation=LogFormatterSciNotation)
        exec(compile(ast.Module(body=[function],type_ignores=[]),str(path),'exec'),cls.scope)

    def tearDown(self):
        plt.close('all')

    def render(self, lower, upper, ticks=None, logarithmic=False):
        fig, ax = plt.subplots(figsize=(10,5))
        norm = LogNorm(lower,upper) if logarithmic else Normalize(lower,upper)
        mesh = ax.pcolormesh(np.ones((3,3)), norm=norm)
        self.scope['add_colorbar'](fig,mesh)
        cax = fig.axes[-1]
        if ticks is not None:
            cax.yaxis.set_major_locator(FixedLocator(ticks))
        fig.canvas.draw()
        return cax.yaxis.get_major_formatter()

    def test_default_reference_zone_is_01_to_1_km(self):
        self.assertEqual(SCHEMA_QCK_VLDR['y_max_zone']['default'],[.1,1.])
        self.assertEqual(SCHEMA_QCK_VLDR['y_lims']['default'],[])

    def test_short_labels_stay_decimal_without_a_gutter_penalty(self):
        formatter=self.render(0.,.02,[0.,.005,.01,.015,.02])
        self.assertEqual(formatter.max_label,'-0.001')
        self.assertFalse(formatter.check_panel_gap)
        self.assertEqual(formatter.orderOfMagnitude,0)

    def test_long_labels_use_shared_scientific_scale(self):
        formatter=self.render(-.00025,.00025,[-.00025,0.,.00025])
        self.assertNotEqual(formatter.orderOfMagnitude,0)

    def test_minus_sign_can_trigger_scientific_scale(self):
        positive=self.render(.0012,.0016,[.0012,.0014,.0016])
        negative=self.render(-.0016,-.0012,[-.0016,-.0014,-.0012])
        self.assertEqual(positive.orderOfMagnitude,0)
        self.assertNotEqual(negative.orderOfMagnitude,0)

    def test_1000_forces_scientific_even_when_the_text_fits(self):
        formatter=self.render(0.,1000.,[0.,250.,500.,750.,1000.])
        self.assertEqual(formatter.orderOfMagnitude,3)

    def test_very_large_and_negative_values_use_scientific_scale(self):
        for lower,upper,ticks in ((0.,1.e8,[0.,5.e7,1.e8]),
                                  (-1000.,0.,[-1000.,-500.,0.])):
            with self.subTest(limits=(lower,upper)):
                formatter=self.render(lower,upper,ticks)
                self.assertNotEqual(formatter.orderOfMagnitude,0)

    def test_logarithmic_colorbar_keeps_log_scientific_formatter(self):
        formatter=self.render(1.e-5,.1,logarithmic=True)
        self.assertIsInstance(formatter,LogFormatterSciNotation)


if __name__=='__main__':
    unittest.main()
