"""Coordinate-only validation tests; runnable without scientific dependencies."""
import contextlib
import io
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "atlas_actris"))
from processor.channel_validation import (
    calibration_pair_ids, drop_measurement, resolve_measurement,
    telecover_channels_match, telecover_test_allowed, validate_channels, warn_once,
)


class Values(list):
    def tolist(self):
        return list(self)


def coordinate(values):
    return SimpleNamespace(values=Values(values))


def profile(*channels):
    return SimpleNamespace(channel=coordinate(channels))


class PairTable:
    """Small labelled table exposing the coordinate API used by validation."""
    dims = ("parameters", "pair")
    parameters = coordinate(["ch_r", "ch_t"])

    def __init__(self, rows):
        self.rows = rows
        self.pair = coordinate(rows)

    def sel(self, *, parameters, pair):
        value = self.rows[pair][0 if parameters == "ch_r" else 1]
        return SimpleNamespace(item=lambda: value)


class ChannelValidationTests(unittest.TestCase):
    def setUp(self):
        self.caller = {"process": ["pcb", "tlc"], "process_qck": ["pcb", "tlc"],
                       "process_bgd": ["pcb", "tlc"], "loading_map": {}}
        self.output = io.StringIO()
        self.redirect = contextlib.redirect_stdout(self.output)
        self.redirect.__enter__()
        self.addCleanup(self.redirect.__exit__, None, None, None)
        self.info = PairTable({5: ("r1", "t1"), 9: ("r2", "t2")})
        self.profiles = {key: profile("r1", "t1", "r2", "t2")
                         for key in ("ray", "pcb_p45", "pcb_m45")}
        self.metadata = {key: self.info for key in self.profiles}

    def active(self):
        return calibration_pair_ids(self.profiles, self.metadata, self.caller)

    def test_primary_does_not_require_auxiliary(self):
        self.assertEqual(self.active(), {"pcb": [5, 9], "pcb_aux": []})

    def test_missing_primary_channel_rejects_only_its_pair(self):
        self.profiles["pcb_m45"] = profile("r1", "t1", "r2")
        self.assertEqual(self.active()["pcb"], [5])
        self.assertIn("pcb_m45: t2", self.output.getvalue())

    def test_ray_is_required_and_ray_pcb_is_not_a_substitute(self):
        self.profiles["ray_pcb"] = self.profiles.pop("ray")
        self.assertEqual(self.active()["pcb"], [])

    def test_explicit_ray_alias_is_resolved(self):
        self.profiles["physical_ray"] = self.profiles.pop("ray")
        self.caller["loading_map"]["ray"] = "physical_ray"
        self.assertEqual(self.active()["pcb"], [5, 9])

    def test_auxiliary_missing_channel_does_not_reject_primary(self):
        self.profiles.update(pcb_aux_p45=profile("r1", "t1", "r2"),
                             pcb_aux_m45=profile("r1", "t1", "r2", "t2"))
        self.assertEqual(self.active(), {"pcb": [5, 9], "pcb_aux": [5]})

    def test_single_auxiliary_orientation_keeps_primary(self):
        self.profiles["pcb_aux_p45"] = profile("r1", "t1")
        self.assertEqual(self.active(), {"pcb": [5, 9], "pcb_aux": []})
        self.assertIn("missing pcb_aux_m45", self.output.getvalue())

    def test_auxiliary_cannot_activate_a_rejected_primary_pair(self):
        self.profiles["ray"] = profile("r1", "t1")
        self.profiles.update(pcb_aux_p45=profile("r1", "t1", "r2", "t2"),
                             pcb_aux_m45=profile("r1", "t1", "r2", "t2"))
        self.assertEqual(self.active(), {"pcb": [5], "pcb_aux": [5]})

    def test_telecover_order_does_not_matter(self):
        sectors = {"tlc_north": profile("a", "b"), "tlc_south": profile("b", "a")}
        self.assertTrue(telecover_channels_match("tlc", sectors, self.caller))

    def test_telecover_mismatch_does_not_disable_diagnostics(self):
        sectors = {"tlc_north": profile("a", "b"), "tlc_south": profile("a")}
        validate_channels(self.caller, sectors, {})
        self.assertFalse(telecover_test_allowed("tlc", sectors, self.caller))
        self.assertEqual(self.caller["process_qck"], ["pcb", "tlc"])
        self.assertEqual(self.caller["process_bgd"], ["pcb", "tlc"])

    def test_ring_check_is_independent(self):
        sectors = {"tlc_inner": profile("a"), "tlc_outer": profile("a", "b")}
        self.assertFalse(telecover_channels_match("tlc_rin", sectors, self.caller))

    def test_removed_sector_cannot_be_silently_ignored(self):
        sectors = {"tlc_north": profile("a"), "tlc_south": profile()}
        metadata = {name: {key: object() for key in sectors} for name in ("shots", "time_info", "channel_info")}
        drop_measurement("tlc_south", sectors, metadata, self.caller)
        self.assertTrue(all("tlc_south" not in group for group in metadata.values()))
        self.assertFalse(telecover_channels_match("tlc", sectors, self.caller))

    def test_sector_removed_later_in_processing_is_not_ignored(self):
        sectors = {"tlc_north": profile("a"), "tlc_south": profile("a")}
        validate_channels(self.caller, sectors, {})
        sectors.pop("tlc_south")
        self.assertFalse(telecover_test_allowed("tlc", sectors, self.caller))

    def test_alias_cycles_terminate(self):
        self.caller["loading_map"] = {"a": "b", "b": "a"}
        self.assertIsNone(resolve_measurement("a", {}, self.caller))

    def test_repeated_checks_emit_warning_once(self):
        warn_once(self.caller, "one warning")
        warn_once(self.caller, "one warning")
        self.assertEqual(self.output.getvalue().count("one warning"), 1)


if __name__ == "__main__":
    unittest.main()
