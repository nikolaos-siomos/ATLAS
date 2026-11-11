"""Offline tests of elastic H consistency checks and SCC export retention."""
import contextlib
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "atlas_actris"))
# requests is used only to download SCC data. These tests use in-memory CSV
# data, so importing the exporter must not require networking dependencies.
try:
    import requests
except ImportError:
    sys.modules["requests"] = types.ModuleType("requests")
    try:
        from utils import get_scc_config as scc
    finally:
        del sys.modules["requests"]
else:
    from utils import get_scc_config as scc


class SCCChannelTypeTests(unittest.TestCase):
    def classify(self, label, H):
        default = "t" if label.startswith("elT") else "p" if label.startswith("elPP") else "c"
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = scc._elastic_type_from_H(label, H, "1996", default, atlas_channel_id="0532f?pt")
        return result, output.getvalue()

    def test_finite_thresholds_and_label_variants(self):
        rules = {
            "elT": [(-1., ""), (-.3, ""), (-.299, "t"), (0., "t"), (.299, "t"), (.3, ""), (1., "")],
            "elPP": [(-1., ""), (-.3, ""), (-.299, ""), (0., ""), (.299, ""), (.3, "p"), (1., "p")],
            "elCP": [(-1., "c"), (-.301, "c"), (-.3, ""), (-.299, ""), (0., ""), (.3, ""), (1., "")],
        }
        for family, cases in rules.items():
            suffixes = ("", "nr", "fr", "unr") if family == "elT" else ("", "nr", "fr")
            for suffix in suffixes:
                for H, expected in cases:
                    with self.subTest(label=family+suffix, H=H):
                        actual, warning = self.classify(family+suffix, str(H))
                        self.assertEqual(actual, expected)
                        self.assertEqual(bool(warning), expected == "")
                        if warning:
                            self.assertIn("SCC channel 1996 [0532f?pt]:", warning)
                            expected = {"elT": "abs(H) < 0.3", "elPP": "H >= 0.3", "elCP": "H < -0.3"}[family]
                            self.assertIn(expected, warning)
                            self.assertNotIn("finite", warning)
                            self.assertTrue(warning.rstrip().endswith("Please revise the SCC configuration."))

    def test_missing_invalid_and_nonfinite_H(self):
        for H in (None, "", "_", "bad", "nan", "inf", "-inf", float("nan"), float("inf")):
            for label in ("elT", "elPP", "elCP"):
                with self.subTest(label=label, H=H):
                    actual, warning = self.classify(label, H)
                    self.assertEqual(actual, "t" if label == "elT" else "")
                    self.assertEqual(bool(warning), label != "elT")

    def test_Raman_types_are_unchanged(self):
        for label, default in (("vrRN2", "v"), ("vrRN2fr", "r"), ("vrRH2O", "v"),
                               ("pRRlow", "r"), ("pRRhighnr", "r")):
            for H in (None, -1., 0., 1.):
                with self.subTest(label=label, H=H):
                    with contextlib.redirect_stdout(io.StringIO()) as output:
                        self.assertEqual(scc._elastic_type_from_H(label, H, "1", default), default)
                    self.assertEqual(output.getvalue(), "")

    @staticmethod
    def data():
        return {
            "HoiChannels": {"id": ["1", "2", "3", "4", "5", "6"],
                            "_scat_mechanism_id_id": ["elT", "elPPnr", "elCPfr", "elPP", "vrRN2", "pRRlow"],
                            "emission_wavelength": ["532"]*6,
                            "if_center": ["532", "532", "532", "532", "607", "530"]},
            "PolarizationCrosstalkParameter": {
                "channel_id": ["1", "2", "3", "4", "5", "6"],
                "g": ["1"]*6, "h": ["0.8", "0.02", "-0.9", "-0.8", "0.7", "0.9"]},
            "PRODUCSTCHANNELS": {"channelID": ["1", "2", "3", "4"], "signalTypeID": ["99", "7", "6", "7"]},
            "atlas_channels": {"scc_channel_id": ["1", "2", "3", "4", "5", "6"],
                               "telescope_type": ["f"]*6, "acquisition_mode": ["p"]*6},
            "atlas_system": {"station_id": ["test"]},
        }

    def test_channel_classification_and_subtypes_together(self):
        data = self.data()
        with contextlib.redirect_stdout(io.StringIO()):
            scc._get_channel_types(data)  # Also loads H for standalone callers.
        self.assertEqual(data["atlas_channels"]["channel_type"], ["", "", "c", "", "v", "r"])
        self.assertEqual(data["atlas_channels"]["channel_subtype"], ["x", "t", "r", "t", "n", "l"])
        self.assertEqual(data["atlas_channels"]["H"], ["0.8", "0.02", "-0.9", "-0.8", "0.7", "0.9"])

    def test_missing_polarization_section(self):
        for parameters in (None, {}, {"channel_id": []}):
            data = self.data()
            data["PolarizationCrosstalkParameter"] = parameters
            with contextlib.redirect_stdout(io.StringIO()):
                scc._get_channel_types(data)
            self.assertEqual(data["atlas_channels"]["channel_type"], ["t", "", "", "", "v", "r"])
            self.assertEqual(data["atlas_channels"]["H"], ["_"]*6)

    def test_last_crosstalk_record_wins_and_missing_rows_are_aligned(self):
        data = self.data()
        data["PolarizationCrosstalkParameter"] = {
            "channel_id": [" 2 ", "2", "3"], "h": ["-0.9", "0E-8", "-0.9"], "g": ["1", "2", "3"]}
        scc._get_GH(data)
        self.assertEqual(data["atlas_channels"]["H"], ["_", "0.", "-0.9", "_", "_", "_"])
        self.assertEqual(data["atlas_channels"]["G"], ["_", "2", "3", "_", "_", "_"])

    def test_export_pipeline_loads_H_before_classification_and_ID_generation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"configuration.ini"
            with contextlib.ExitStack() as stack:
                stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                stack.enter_context(patch.object(scc, "_get_scc_file", return_value="in-memory.csv"))
                stack.enter_context(patch.object(scc, "read_sectioned_csv", return_value=self.data()))
                stack.enter_context(patch.object(scc, "_copy_atlas_keys", side_effect=lambda data, **kwargs: data))
                for name in ("_get_zero_bin", "_get_telescope_type", "_get_acquisition_mode",
                             "_get_laser_parameters", "_get_dead_time"):
                    stack.enter_context(patch.object(scc, name, side_effect=lambda data: data))
                result = scc.export_scc_config("713", str(path), "1", directory,
                                               scc_compatible_format=True)
            channels = result["atlas_channels"]
            types_ = dict(zip(channels["scc_channel_id"], channels["channel_type"]))
            self.assertEqual(types_["2"], "")
            self.assertEqual(types_["3"], "c")
            self.assertEqual(types_["1"], "")
            ids = dict(zip(channels["scc_channel_id"], channels["atlas_channel_id"]))
            self.assertEqual(len(ids["2"]), 8)
            self.assertEqual(ids["2"][5], "?")
            self.assertEqual(ids["3"][5], "c")
            self.assertTrue(path.exists())

    def test_warnings_use_exact_incomplete_dropdown_ID_and_expected_condition(self):
        data = self.data()
        with contextlib.redirect_stdout(io.StringIO()) as output:
            scc._get_channel_types(data)
            self.assertEqual(output.getvalue(), "")  # IDs not constructed yet.
            scc._get_atlas_channel_id(data)
        warnings = output.getvalue()
        ids = dict(zip(data["HoiChannels"]["id"], data["atlas_channels"]["atlas_channel_id"]))
        self.assertIn(f"SCC channel 1 [{ids['1']}]:", warnings)
        self.assertIn(f"SCC channel 2 [{ids['2']}]:", warnings)
        self.assertIn("expected abs(H) < 0.3 for total channels", warnings)
        self.assertIn("expected H >= 0.3 for parallel channels", warnings)
        self.assertNotIn("_pending_channel_warnings", data)
        with contextlib.redirect_stdout(io.StringIO()) as again:
            scc._get_atlas_channel_id(data)
        self.assertEqual(again.getvalue(), "")

    def test_ambiguous_subtype_warning_also_includes_ID_and_expectation(self):
        data = self.data()
        data["PRODUCSTCHANNELS"] = {"channelID": ["3", "3"], "signalTypeID": ["6", "7"]}
        with contextlib.redirect_stdout(io.StringIO()) as output:
            scc._get_channel_types(data)
            scc._get_atlas_channel_id(data)
        atlas_id = data["atlas_channels"]["atlas_channel_id"][2]
        matching = [line for line in output.getvalue().splitlines() if "ambiguous subtype" in line]
        self.assertEqual(len(matching), 1)
        self.assertIn(f"SCC channel 3 [{atlas_id}]:", matching[0])
        self.assertIn("expected one reflected/transmitted assignment", matching[0])
        self.assertTrue(matching[0].endswith("Please revise the SCC configuration."))

    def test_ID_placeholders_preserve_positions_and_empty_parameters(self):
        data = self.data()
        with contextlib.redirect_stdout(io.StringIO()):
            scc._get_channel_types(data)
            data["atlas_channels"]["telescope_type"][0] = ""
            data["atlas_channels"]["acquisition_mode"][0] = None
            data["atlas_channels"]["channel_subtype"][0] = ""
            scc._get_atlas_channel_id(data)
        self.assertEqual(data["atlas_channels"]["atlas_channel_id"][0], "0532????")
        self.assertEqual(data["atlas_channels"]["channel_type"][0], "")
        self.assertTrue(all(len(value) == 8 for value in data["atlas_channels"]["atlas_channel_id"]))

    def test_underscore_displays_as_question_mark_but_keeps_exclusion_policy(self):
        data = self.data()
        with contextlib.redirect_stdout(io.StringIO()):
            scc._get_channel_types(data)
            data["atlas_channels"]["channel_subtype"][0] = "_"
            scc._get_atlas_channel_id(data)
            self.assertEqual(data["atlas_channels"]["atlas_channel_id"][0], "0532f?p?")
            scc._exclude_undifined_scc_channels(data)
        self.assertNotIn("1", data["atlas_channels"]["scc_channel_id"])

    def prompted_data(self, answers):
        data = self.data()
        with contextlib.redirect_stdout(io.StringIO()) as output:
            scc._get_channel_types(data)
            scc._get_atlas_channel_id(data)
            with patch("builtins.input", side_effect=answers):
                scc._get_recorder_channel_id(data)
        return data, output.getvalue()

    def test_partial_recorder_IDs_retain_all_channels_and_export_empty_slots(self):
        answers = ["1", "", "3", "", "5", "6"]
        data, output = self.prompted_data(answers)
        channels = data["atlas_channels"]
        self.assertEqual(channels["recorder_channel_id"], answers)
        self.assertEqual(len(channels["scc_channel_id"]), 6)
        self.assertTrue(all(len(value) == 6 for value in channels.values() if isinstance(value, list)))
        self.assertIn("Press Enter to skip", output)
        self.assertIn("partially filled", output)
        self.assertIn("complete empty entries manually", output)
        self.assertNotIn("Excluded!", output)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"configuration.ini"
            with contextlib.redirect_stdout(io.StringIO()):
                scc.export_selected_to_ini(data, str(path), scc.sections_map)
            import configparser
            config = configparser.ConfigParser()
            config.read(path)
            exported = [v.strip() for v in config["Channels"]["recorder_channel_id"].split(",")]
            self.assertEqual(exported, answers)

    def test_all_recorder_IDs_can_be_skipped(self):
        data, output = self.prompted_data([""]*6)
        self.assertEqual(data["atlas_channels"]["recorder_channel_id"], [""]*6)
        self.assertEqual(len(data["atlas_channels"]["scc_channel_id"]), 6)
        self.assertIn("recorder_channel_id is empty", output)

    def test_complete_recorder_IDs_have_no_partial_warning(self):
        data, output = self.prompted_data(["1", "2", "3", "4", "5", "6"])
        self.assertEqual(data["atlas_channels"]["recorder_channel_id"], ["1", "2", "3", "4", "5", "6"])
        self.assertNotIn("recorder_channel_id is partially filled", output)
        self.assertNotIn("recorder_channel_id is empty", output)

    def test_duplicate_filled_recorder_IDs_are_rejected_with_or_without_blanks(self):
        for answers in (["1", "1", "3", "4", "5", "6"], ["1", "", "1", "", "", ""]):
            with self.subTest(answers=answers):
                with self.assertRaisesRegex(scc.ConfigError, "Duplicates"):
                    self.prompted_data(answers)

    def test_whitespace_only_recorder_ID_is_a_skipped_entry(self):
        data, output = self.prompted_data([" 1 ", " ", "3", "4", "5", "6"])
        self.assertEqual(data["atlas_channels"]["recorder_channel_id"][:2], ["1", ""])
        self.assertIn("partially filled", output)

    def test_rejected_types_remain_empty_in_exported_configuration(self):
        data = self.data()
        with contextlib.redirect_stdout(io.StringIO()):
            scc._get_GH(data)
            scc._get_channel_types(data)
            scc._get_atlas_channel_id(data)
            scc._exclude_undifined_scc_channels(data)
        self.assertEqual(len(data["atlas_channels"]["scc_channel_id"]), 6)
        self.assertEqual(data["atlas_channels"]["channel_type"].count(""), 3)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"configuration.ini"
            scc.export_selected_to_ini(data, str(path), scc.sections_map)
            import configparser
            config = configparser.ConfigParser()
            config.read(path)
            ids = [v.strip() for v in config["Channels"]["scc_channel_id"].split(",")]
            types_ = [v.strip() for v in config["Channels"]["channel_type"].split(",")]
            mapping = dict(zip(ids, types_))
            self.assertEqual(mapping, {"1": "", "2": "", "3": "c", "4": "", "5": "v", "6": "r"})


if __name__ == "__main__":
    unittest.main()
