"""Numerical regressions for measurements with different channel sets.

Run in the ATLAS environment: python -m pytest tests/test_channel*.py
"""
import contextlib
import io
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "atlas_actris"))
try:
    import numpy as np
    import pandas as pd
    import xarray as xr
except ImportError:
    xr = None


@unittest.skipIf(xr is None, "Requires ATLAS numerical dependencies (numpy, pandas, xarray).")
class ChannelProcessingTests(unittest.TestCase):
    def setUp(self):
        self.redirect = contextlib.redirect_stdout(io.StringIO())
        self.redirect.__enter__()
        self.addCleanup(self.redirect.__exit__, None, None, None)

    @staticmethod
    def signal(channels, times=(0,), value=10.0):
        return xr.DataArray(np.full((len(times), len(channels), 3), value),
                            dims=("time", "channel", "bins"),
                            coords={"time": list(times), "channel": list(channels), "bins": [0, 1, 2]})

    def test_selection_and_id_mapping_use_labels_not_positions(self):
        from processor.modify import remove_unrecognised_channels, get_atlas_channel_id
        config = dict(recorder_channel_id=["2", "1", "3"],
                      detected_wavelength=[355., 532., 607.], telescope_type=["f"]*3,
                      channel_type=["c", "p", "v"], channel_subtype=["r", "t", "n"],
                      acquisition_mode=["p"]*3)
        profiles = {"ray": self.signal(["1", "3", "2"]),
                    "tlc_north": self.signal(["1", "2", "unknown"])}
        metadata = {"shots": {key: arr.isel(bins=0, drop=True) for key, arr in profiles.items()}}
        remove_unrecognised_channels({}, config, profiles, metadata)
        self.assertEqual(profiles["tlc_north"].channel.values.tolist(), ["2", "1"])
        get_atlas_channel_id(config, profiles, metadata)
        self.assertEqual(profiles["tlc_north"].channel.values.tolist(), ["0355fcpr", "0532fppt"])
        xr.testing.assert_equal(profiles["tlc_north"].channel, metadata["shots"]["tlc_north"].channel)

    def test_empty_measurement_cleanup_and_global_endpoint(self):
        from processor.modify import remove_unrecognised_channels
        profiles = {"ray": self.signal(["1"]), "tlc_north": self.signal(["x"])}
        metadata = {"shots": {key: arr.isel(bins=0, drop=True) for key, arr in profiles.items()},
                    "channel_info": {key: object() for key in profiles}}
        caller = {}
        remove_unrecognised_channels(caller, {"recorder_channel_id": ["1"]}, profiles, metadata)
        self.assertEqual(list(profiles), ["ray"])
        self.assertNotIn("tlc_north", metadata["channel_info"])
        with self.assertRaisesRegex(SystemExit, "Termination ID 2"):
            remove_unrecognised_channels(caller, {"recorder_channel_id": ["absent"]}, profiles, metadata)

    def test_metadata_loading_does_not_shift_channel_or_pair_parameters(self):
        from processor.modify import load_metadata
        config = {"recorder_channel_id": ["1", "2", "3"], "bins": [None]*3,
                  "station_id": "station", "K": [0.8, 0.9], "ch_r": ["r1", "r2"]}
        metadata = {"time_info": {"ray": object(), "tlc_north": object()},
                    "system_info": {key: pd.Series(dtype=object) for key in ("ray", "tlc_north")},
                    "channel_info": {
                        "ray": pd.DataFrame({"bins": [3000, 1000]}, index=["3", "1"]),
                        "tlc_north": pd.DataFrame({"bins": [2000]}, index=["2"]),
                    }}
        result = load_metadata(config, metadata)
        self.assertEqual(result["bins"], ["1000", "2000", "3000"])
        self.assertEqual(result["K"], [0.8, 0.9])
        self.assertEqual(result["ch_r"], ["r1", "r2"])
        self.assertEqual(config["bins"], [None]*3)

    def test_store_metadata_matches_each_measurements_channels(self):
        from processor.modify import store_updated_metadata
        config = {"atlas_channel_id": ["a", "b", "c"], "recorder_channel_id": ["1", "2", "3"],
                  "range_resolution": [1., 2., 3.], "station_id": "station", "K": [0.8, 0.9],
                  "ch_r": ["a", "b"], "ch_t": ["b", "c"]}
        profiles = {"ray": self.signal(["a", "b", "c"]), "tlc_north": self.signal(["c", "a"])}
        metadata = {"time_info": {key: pd.DataFrame({"filename": ["file"]}, index=[0]) for key in profiles},
                    "system_info": {}, "channel_info": {}}
        result = store_updated_metadata(config, metadata, profiles)
        self.assertEqual(result["channel_info"]["tlc_north"].channel.values.tolist(), ["c", "a"])
        self.assertEqual(result["channel_info"]["tlc_north"].sel(parameters="range_resolution").values.tolist(), [3., 1.])
        self.assertEqual(result["pol_cal_info"]["ray"].sel(parameters="K").values.tolist(), [0.8, 0.9])

    def test_packaging_keeps_unique_channels_and_missing_periods_for_all_groups(self):
        from processor.packaging import combine_QA_pack, packing_map
        for group, sources in packing_map.items():
            with self.subTest(group=group):
                first, second = sources[:2]
                a = self.signal(["shared", "first_only"], (0, 1), value=1.)
                b = self.signal(["shared", "second_only"], (2, 3), value=2.)
                metadata_a = xr.DataArray([[10., 20.]], dims=("parameters", "channel"),
                                          coords={"parameters": ["resolution"], "channel": a.channel})
                metadata_b = xr.DataArray([[10., 30.]], dims=("parameters", "channel"),
                                          coords={"parameters": ["resolution"], "channel": b.channel})
                data = {"profile": {first: a, second: b},
                        "background": {first: a.isel(bins=0, drop=True), second: b.isel(bins=0, drop=True)},
                        "channel_info": {first: metadata_a, second: metadata_b},
                        "range": {first: a.isel(time=0, drop=True), second: b.isel(time=0, drop=True)}}
                result = combine_QA_pack(data)
                combined = result["profile"][group]
                self.assertEqual(set(combined.channel.values), {"shared", "first_only", "second_only"})
                self.assertTrue(bool(combined.sel(channel="second_only", time=[0, 1]).isnull().all()))
                self.assertTrue(bool(combined.sel(channel="first_only", time=[2, 3]).isnull().all()))
                self.assertEqual(result["channel_info"][group].sel(parameters="resolution", channel="second_only").item(), 30.)
                self.assertEqual(set(result["range"][group].channel.values), set(combined.channel.values))
                xr.testing.assert_identical(data["profile"][first], a)
                self.assertNotIn(group, data["profile"])

    def test_dark_correction_preserves_missing_and_photon_channels(self):
        from processor.signal_processing import compute_dark_correction
        from processor.definitions import profile_instances
        sig = self.signal(["analog_present", "analog_absent", "photon"], (1, 2))
        dark = self.signal(["analog_present"], (0,), value=2.)
        data = {key: {"ray": sig.copy()} for key in profile_instances}
        data["profile_mean"]["drk"] = dark
        data["channel_info"] = {"ray": xr.DataArray([["a", "a", "p"]],
                                 dims=("parameters", "channel"),
                                 coords={"parameters": ["acquisition_mode"], "channel": sig.channel})}
        caller = {"loading_map": {"drk_ray": "drk"}}
        result = compute_dark_correction({"caller_info": caller}, data)
        corrected = result["profile"]["ray"]
        self.assertEqual(corrected.channel.values.tolist(), sig.channel.values.tolist())
        np.testing.assert_allclose(corrected.sel(channel="analog_present"), 8.)
        np.testing.assert_allclose(corrected.sel(channel="analog_absent"), 10.)
        np.testing.assert_allclose(corrected.sel(channel="photon"), 10.)
        xr.testing.assert_equal(data["profile"]["ray"], sig)

    def test_auxiliary_factor_is_never_applied_to_ray(self):
        from processor.pol_cal_calculations import _find_eta_for_pairs, _collect_eta_pairs_from_pcb, make_ratio_id
        r, t = "0532fcpr", "0532fppt"
        pair = make_ratio_id(r, t, "e")
        def info(value):
            return xr.DataArray([[r], [t], ["eta"], [value], [0.1]],
                                dims=("parameters", "pair"),
                                coords={"parameters": ["ch_r", "ch_t", "ratio_type", "mean", "sem"], "pair": [pair]})
        data = {"pol_cal_info": {"pcb": info(1.1), "pcb_aux": info(1.9)}}
        eta, _ = _find_eta_for_pairs(data, [r], [t], ["calibrated"])
        self.assertAlmostEqual(float(eta.item()), 1.1)
        data["pol_cal_info"].pop("pcb")
        self.assertEqual(_collect_eta_pairs_from_pcb(data)[0], [])
        with self.assertRaises(KeyError):
            _find_eta_for_pairs(data, [r], [t], ["calibrated"])

    def test_active_pair_selection_keeps_its_original_calibration_values(self):
        from processor.channel_validation import calibration_pair_ids
        info = xr.DataArray([["r1", "r2"], ["t1", "t2"], [0.8, 0.9]],
                            dims=("parameters", "pair"),
                            coords={"parameters": ["ch_r", "ch_t", "K"], "pair": [5, 9]})
        profiles = {key: self.signal(["r2", "t2"]) for key in ("pcb_p45", "pcb_m45", "ray")}
        active = calibration_pair_ids(profiles, {"pcb_p45": info}, {})
        selected = info.sel(pair=active["pcb"])
        self.assertEqual(selected.pair.values.tolist(), [9])
        self.assertAlmostEqual(float(selected.sel(parameters="K").item()), 0.9)

    def test_lazy_packaging_does_not_compute_signal_data(self):
        try:
            import dask.array as da
            from dask import delayed
        except ImportError:
            self.skipTest("Requires dask")
        from processor.packaging import combine_QA_pack
        @delayed
        def unavailable_data():
            raise AssertionError("Signal data were computed during packaging")
        lazy = da.from_delayed(unavailable_data(), shape=(1, 1, 3), dtype=float)
        arr = xr.DataArray(lazy, dims=("time", "channel", "bins"),
                           coords={"time": [0], "channel": ["first"], "bins": [0, 1, 2]})
        other = arr.assign_coords(time=[1], channel=["second"])
        result = combine_QA_pack({"profile": {"pcb_p45": arr, "pcb_m45": other}})
        self.assertIsInstance(result["profile"]["pcb"].data, da.Array)


if __name__ == "__main__":
    unittest.main()
