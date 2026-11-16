"""Report metadata must remain visible text, never HTML markup."""
import contextlib
import html
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from PIL import Image, PngImagePlugin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src/atlas_actris"))
from visualizer import export_html as reports


class ReportHtmlSecurityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="atlas-html-test-")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.payload = '0355xta<script>window.ATLAS_TEST=true</script>&"'

    def test_real_png_metadata_is_escaped_in_generated_report(self):
        metadata = PngImagePlugin.PngInfo()
        metadata.add_text("QA_test_ID", "qck_ray")
        metadata.add_text("atlas_channel_id", self.payload)
        with Image.new("RGB", (2, 2)) as image:
            image.save(self.folder / "plot.png", pnginfo=metadata)
        output = self.folder / "ATLAS & test.html"
        reports.QA_report(str(self.folder), str(output), export_docx=False)
        result = output.read_text(encoding="utf-8")
        self.assertIn(html.escape(self.payload), result)
        self.assertIn("ATLAS &amp; test.html", result)
        self.assertNotIn("<script>", result)
        self.assertIn("data:image/png;base64,", result)

    def test_every_report_heading_branch_escapes_supplied_channel_text(self):
        image = {"data_uri": "data:image/png;base64,audit"}
        data = {key: {} for key in ("qck_ray", "qck_ray_pcb", "qck_tlc", "qck_tlc_rin",
                                    "qck_pcb", "qck_drk", "drk", "vldr", "ray", "ray_pcb",
                                    "tlc", "tlc_rin", "pcb", "pcb_aux")}
        for key in ("qck_ray", "drk", "vldr", "ray", "tlc"):
            data[key] = {self.payload: dict(image)}
        data["pcb"] = {"pair": dict(image, atlas_channel_id_r=self.payload,
                                     atlas_channel_id_t=self.payload + "T")}
        output = self.folder / "report.html"
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(reports, "_collect_report_data", return_value=data))
            for name in ("channel_limit_table", "channel_background_table", "channel_polarization"):
                stack.enter_context(patch.object(reports, name))
            reports.QA_report(str(self.folder), str(output), export_docx=False)
        result = output.read_text(encoding="utf-8")
        self.assertNotIn("<script>", result)
        self.assertNotIn("<SCRIPT>", result)
        self.assertGreaterEqual(result.count(html.escape(self.payload)), 6)
        self.assertIn(html.escape(self.payload.upper()), result)
        self.assertIn(html.escape(self.payload) + " to " + html.escape(self.payload + "T"), result)


if __name__ == "__main__":
    unittest.main()
