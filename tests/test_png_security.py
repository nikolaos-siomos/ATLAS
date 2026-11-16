"""Real-image regressions for PNG-only plot processing, without lidar I/O."""
import ast
import contextlib
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

import numpy as np
from PIL import Image, PngImagePlugin, UnidentifiedImageError

SOURCE = Path(__file__).resolve().parents[1] / "src" / "atlas_actris"
sys.path.insert(0, str(SOURCE))
from visualizer import export_html as reports
from utils.png_images import open_png


def _plot_image_functions():
    # These Pillow operations do not need xarray or the full processing stack.
    # Compile their actual source, following the lightweight plotting tests.
    path = SOURCE / "visualizer" / "plot_utils.py"
    names = {"add_plot_metadata", "perform_color_reduction"}
    functions = [node for node in ast.parse(path.read_text(encoding="utf-8")).body
                 if isinstance(node, ast.FunctionDef) and node.name in names]
    namespace = {"Image": Image, "PngImagePlugin": PngImagePlugin, "np": np,
                 "open_png": open_png}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


class PngSecurityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.functions = _plot_image_functions()

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas-png-test-")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def png(self, name="plot.png", **metadata):
        path = self.folder / name
        info = PngImagePlugin.PngInfo()
        for key, value in metadata.items():
            info.add_text(key, value)
        with Image.new("RGB", (8, 6), (30, 60, 90)) as image:
            image.save(path, "PNG", pnginfo=info)
        return path

    def foreign_images(self, prefix="foreign"):
        paths = []
        for format_name in ("GIF", "BMP", "JPEG"):
            path = self.folder / (prefix + "_" + format_name + ".png")
            with Image.new("RGB", (2, 2), (30, 60, 90)) as image:
                image.save(path, format_name)
            paths.append(path)
        # Minimal, benign 1x1 RGB PSD with uncompressed pixels. No exploit data.
        path = self.folder / (prefix + "_PSD.png")
        path.write_bytes(struct.pack(">4sH6sHIIHH", b"8BPS", 1, b"\0" * 6, 3, 1, 1, 8, 3)
                         + b"\0" * 12 + b"\0\0" + bytes((30, 60, 90)))
        paths.append(path)
        return paths

    @contextlib.contextmanager
    def foreign_decoder_guards(self):
        Image.init()
        guards = {name: Mock(side_effect=AssertionError("Non-PNG decoder called"))
                  for name in ("GIF", "BMP", "JPEG", "PSD")}
        replacements = {name: (guard, Image.OPEN[name][1]) for name, guard in guards.items()}
        with patch.dict(Image.OPEN, replacements):
            yield
        for guard in guards.values():
            guard.assert_not_called()

    def test_metadata_update_preserves_pixels_and_text(self):
        path = self.png()
        self.functions["add_plot_metadata"](
            path, {"station": "München", "QA_test_ID": "qck_ray", "ignored": None},
            {"atlas_channel_id": "0355xta", "bins": np.array([1, 2])})
        with Image.open(path, formats=["PNG"]) as image:
            self.assertEqual(image.size, (8, 6))
            self.assertEqual(image.getpixel((0, 0)), (30, 60, 90))
            self.assertEqual(image.text["station"], "München")
            self.assertEqual(image.text["atlas_channel_id"], "0355xta")
            self.assertNotIn("ignored", image.text)
        path.rename(self.folder / "renamed.png")

    def test_colour_reduction_keeps_valid_png_and_allows_metadata_update(self):
        path = self.png()
        self.functions["perform_color_reduction"](True, path)
        self.functions["add_plot_metadata"](path, {"station": "test"})
        with Image.open(path, formats=["PNG"]) as image:
            self.assertEqual(image.mode, "P")
            self.assertEqual(image.convert("RGB").getpixel((0, 0)), (30, 60, 90))
            self.assertEqual(image.text["station"], "test")
        path.unlink()

    def test_updates_reject_other_formats_renamed_png_without_touching_them(self):
        for path in self.foreign_images():
            before = path.read_bytes()
            for name, args in (("add_plot_metadata", (path, {"station": "test"})),
                               ("perform_color_reduction", (True, path))):
                with self.subTest(file=path.name, function=name), self.foreign_decoder_guards():
                    with self.assertRaises(UnidentifiedImageError):
                        self.functions[name](*args)
                self.assertEqual(path.read_bytes(), before)

    def test_disabled_colour_reduction_does_not_open_or_modify_file(self):
        path = self.folder / "unused.png"
        path.write_bytes(b"Not an image")
        self.functions["perform_color_reduction"](False, path)
        self.assertEqual(path.read_bytes(), b"Not an image")

    def test_animated_png_is_rejected_before_pillow_is_called(self):
        path = self.folder / "animation.png"
        with Image.new("RGB", (2, 2), "red") as first, Image.new("RGB", (2, 2), "blue") as second:
            first.save(path, "PNG", save_all=True, append_images=[second], duration=100)
        with patch.object(Image, "open", side_effect=AssertionError("Pillow must not open APNG")) as decoder:
            with self.assertRaisesRegex(UnidentifiedImageError, "static PNG"):
                with open_png(path):
                    self.fail("Animated PNG accepted")
            decoder.assert_not_called()

    def test_invalid_png_framing_is_rejected_before_pillow_is_called(self):
        signature = b"\x89PNG\r\n\x1a\n"
        cases = (b"", signature, signature + struct.pack(">I4s", 13, b"IHDR"),
                 signature + struct.pack(">I4s", 0, b"IEND") + b"\0" * 4,
                 signature + struct.pack(">I4s", 0x80000000, b"IHDR"))
        for index, content in enumerate(cases):
            path = self.folder / (str(index) + ".png")
            path.write_bytes(content)
            with self.subTest(case=index), patch.object(Image, "open") as decoder:
                with self.assertRaises(UnidentifiedImageError):
                    with open_png(path):
                        self.fail("Malformed PNG accepted")
                decoder.assert_not_called()

    def test_png_image_and_stream_close_after_processing_error(self):
        path = self.png()
        with self.assertRaisesRegex(RuntimeError, "processing stopped"):
            with open_png(path) as image:
                stream = image.fp
                raise RuntimeError("processing stopped")
        self.assertTrue(stream.closed)
        self.assertIsNone(image.fp)

    def test_report_index_skips_other_formats_without_calling_their_decoders(self):
        valid = self.png(QA_test_ID="qck_ray", atlas_channel_id="0355xta")
        self.foreign_images()
        with self.foreign_decoder_guards():
            indexed = reports._index_plot_metadata(str(self.folder))
        self.assertEqual([Path(path) for path, _ in indexed["qck_ray"]], [valid])

    def test_legacy_quicklook_selector_uses_only_png_decoder(self):
        valid = self.png("valid_qck_ray_0355xta.png", atlas_channel_id="0355xta")
        self.foreign_images(prefix="invalid_qck_ray_0355xta")
        with self.foreign_decoder_guards():
            selected = reports._select_quicklook_images(str(self.folder), "ray")
        self.assertEqual([Path(path) for path in selected], [valid])

    def test_combined_report_plot_preserves_png_metadata(self):
        quicklook = self.png("qck.png", QA_test_ID="qck_ray", atlas_channel_id="0355xta")
        background = self.png("bgd.png", QA_test_ID="bgd_ray", atlas_channel_id="0355xta")
        target = self.folder / "combined" / "result.png"
        combined = reports._combine_qck_and_bgd_images(
            {"path": str(quicklook)}, {"path": str(background)}, str(target))
        self.assertTrue(combined["combined_qck_bgd"])
        self.assertTrue(combined["data_uri"].startswith("data:image/png;base64,"))
        with Image.open(target, formats=["PNG"]) as image:
            self.assertEqual(image.size, (28, 6))
            self.assertEqual(image.text["QA_test_ID"], "qck_ray")
            self.assertEqual(image.text["bgd_QA_test_ID"], "bgd_ray")

    def test_combining_rejects_other_formats_in_either_panel(self):
        valid = self.png()
        target = self.folder / "combined" / "rejected.png"
        for foreign in self.foreign_images():
            for left, right in ((foreign, valid), (valid, foreign)):
                with self.subTest(left=left.name, right=right.name), self.foreign_decoder_guards():
                    with self.assertRaises(UnidentifiedImageError):
                        reports._combine_qck_and_bgd_images(
                            {"path": str(left)}, {"path": str(right)}, str(target))
                self.assertFalse(target.exists())

    def test_matplotlib_png_still_produces_html_report(self):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        figure = Figure(figsize=(2, 1), dpi=80)
        FigureCanvasAgg(figure)
        figure.subplots().plot([0, 1, 2], [0, 1, 0])
        path = self.folder / "matplotlib.png"
        figure.savefig(path)
        self.functions["perform_color_reduction"](True, path)
        self.functions["add_plot_metadata"](
            path, {"QA_test_ID": "qck_ray", "atlas_channel_id": "0355xta"})
        output = self.folder / "report.html"
        reports.QA_report(str(self.folder), str(output), export_docx=False)
        self.assertIn("data:image/png;base64,", output.read_text())


if __name__ == "__main__":
    unittest.main()
