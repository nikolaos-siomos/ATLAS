"""VLDR PNGs generated with default metadata must reach both report formats."""
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile
from xml.etree import ElementTree

from PIL import Image, PngImagePlugin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src/atlas_actris"))
from visualizer import export_html as reports


class VldrReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas-vldr-report-")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.pairs = ("0355xva", "1640xvp")
        for index, pair in enumerate(self.pairs):
            info = PngImagePlugin.PngInfo()
            # Same fields as qa_test_quicklook_vldr: no atlas_channel_id.
            for key, value in {"QA_test_ID": "qck_vldr", "vldr_id": pair,
                               "atlas_channel_id_r": pair + "r",
                               "atlas_channel_id_t": pair + "t"}.items():
                info.add_text(key, value)
            with Image.new("RGB", (8, 6), (index * 80, 60, 90)) as image:
                image.save(self.folder / (pair + ".png"), pnginfo=info)

    def test_index_accepts_exact_vldr_id_and_collector_keeps_all_pairs(self):
        info = PngImagePlugin.PngInfo()
        info.add_text("QA_test_ID", "qck_vldr_unrecognized")
        info.add_text("vldr_id", "invalid")
        with Image.new("RGB", (2, 2)) as image:
            image.save(self.folder / "unknown.png", pnginfo=info)
        indexed = reports._index_plot_metadata(str(self.folder))
        self.assertEqual(len(indexed["qck_vldr"]), 2)
        self.assertNotIn("qck_vldr_unrecognized", indexed)
        data = reports._collect_report_data(str(self.folder))
        self.assertEqual(set(data["vldr"]), set(self.pairs))

    def test_html_contains_vldr_section_and_every_pair(self):
        target = self.folder / "report.html"
        reports.QA_report(str(self.folder), str(target), export_docx=False)
        output = target.read_text(encoding="utf-8")
        self.assertIn("<h1>VLDR Quicklooks</h1>", output)
        for pair in self.pairs:
            self.assertIn("<h2>" + pair.upper() + "</h2>", output)
        self.assertEqual(output.count('class="report-plot"'), 2)

    def test_docx_contains_vldr_section_and_embedded_plots(self):
        target = self.folder / "report.html"
        reports.QA_report(str(self.folder), str(target), export_docx=True)
        with zipfile.ZipFile(target.with_suffix(".docx")) as archive:
            document = ElementTree.fromstring(archive.read("word/document.xml"))
            namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
            texts = [node.text for node in document.findall(".//w:t", namespace)]
            self.assertIn("VLDR Quicklooks", texts)
            for pair in self.pairs:
                self.assertIn(pair.upper(), texts)
            self.assertEqual(len([n for n in archive.namelist() if n.startswith("word/media/")]), 2)


if __name__ == "__main__":
    unittest.main()
