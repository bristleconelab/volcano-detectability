"""Synthetic checks with known answers. Run: python -I -m unittest discover -s tests -v"""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import vdlib as V  # noqa: E402


class TestMogi(unittest.TestCase):
    def test_peak_uplift(self):
        # 1e6 m^3 at 5 km: u_z(0) = (1 - nu) dV / (pi d^2) = 0.75e6 / (pi * 25e6) = 9.549 mm
        _, _, uz = V.mogi_greens(0.0, 0.0, 5.0)
        self.assertAlmostEqual(uz * 1e6 * 1e3, 9.5493, places=3)

    def test_horizontal_peak_at_d_over_sqrt2(self):
        r = np.linspace(0.1, 20, 2000)
        ue, _, _ = V.mogi_greens(r, 0 * r, 5.0)
        self.assertAlmostEqual(r[np.argmax(ue)], 5 / np.sqrt(2), delta=0.02)

    def test_uplift_volume_ratio(self):
        # integral of u_z over the surface = 2 (1 - nu) dV (= 1.5 dV for nu = 0.25); out to radius R it is
        # 1.5 * (1 - d / sqrt(R^2 + d^2))
        R, d = 400.0, 5.0
        r = np.linspace(0, R, 400001)
        _, _, uz = V.mogi_greens(r, 0 * r, d)
        vol = np.trapezoid(uz * 2 * np.pi * r * 1e3, r * 1e3)
        self.assertAlmostEqual(vol, 1.5 * (1 - d / np.hypot(R, d)), places=5)


class TestSill(unittest.TestCase):
    def test_peak_twice_mogi(self):
        _, _, us = V.sill_greens(0.0, 0.0, 5.0)
        _, _, um = V.mogi_greens(0.0, 0.0, 5.0)
        self.assertAlmostEqual(us / um, 1.5 / 0.75, places=9)

    def test_uplift_volume_equals_dv(self):
        R, d = 400.0, 5.0
        r = np.linspace(0, R, 400001)
        _, _, uz = V.sill_greens(r, 0 * r, d)
        vol = np.trapezoid(uz * 2 * np.pi * r * 1e3, r * 1e3)
        self.assertAlmostEqual(vol, 1 - (d / np.hypot(R, d)) ** 3, places=5)

    def test_sill_option_used(self):
        dvm, *_ = V.min_detectable_volume([0.0], [0.0], [[1.0, 1.0, 0.002]], 5.0)
        dvs, *_ = V.min_detectable_volume([0.0], [0.0], [[1.0, 1.0, 0.002]], 5.0, source="sill")
        self.assertAlmostEqual(dvm / dvs, 2.0, places=9)


class TestDetectability(unittest.TestCase):
    def test_single_station_overhead(self):
        # 2 mm vertical noise at the epicentre, K = 1: dV_min = 2e-3 / (0.75 / (pi * 25e6)) = 2.094e5 m^3
        dv, i, c = V.min_detectable_volume([0.0], [0.0], [[1.0, 1.0, 0.002]], 5.0, k=1)
        self.assertAlmostEqual(dv, 0.002 * np.pi * 25e6 / 0.75, delta=1)
        self.assertEqual((i, c), (0, 2))

    def test_k_scales_linearly(self):
        dv1, *_ = V.min_detectable_volume([3.0], [4.0], [[0.001, 0.001, 0.003]], 5.0, k=1)
        dv3, *_ = V.min_detectable_volume([3.0], [4.0], [[0.001, 0.001, 0.003]], 5.0, k=3)
        self.assertAlmostEqual(dv3 / dv1, 3.0, places=9)

    def test_best_station_wins(self):
        x, y = [50.0, 2.0], [0.0, 0.0]
        dv, i, _ = V.min_detectable_volume(x, y, [[1e-3] * 3, [1e-3] * 3], 5.0)
        self.assertEqual(i, 1)

    def test_network_single_station_vertical(self):
        # one station overhead, only the vertical informative: dV = snr * sigma / Gz
        _, _, gz = V.mogi_greens(0.0, 0.0, 5.0)
        dv = V.network_min_volume([0.0], [0.0], [[1e9, 1e9, 0.002]], 5.0, snr=5)
        self.assertAlmostEqual(dv, 5 * 0.002 / gz, delta=dv * 1e-9)

    def test_network_more_stations_better(self):
        one = V.network_min_volume([3.0], [0.0], [[1e-3] * 3], 5.0)
        four = V.network_min_volume([3.0, -3.0, 0.0, 0.0], [0.0, 0.0, 3.0, -3.0], [[1e-3] * 3] * 4, 5.0)
        self.assertAlmostEqual(one / four, 2.0, places=9)  # sqrt(4) identical stations

    def test_network_empty(self):
        self.assertTrue(np.isinf(V.network_min_volume([], [], np.zeros((0, 3)), 5.0)))

    def test_no_station_is_blind(self):
        dv, i, c = V.min_detectable_volume([], [], np.zeros((0, 3)), 5.0)
        self.assertTrue(np.isinf(dv))
        self.assertEqual(i, -1)


class TestInsar(unittest.TestCase):
    def test_summit_pixel(self):
        # coherent at the summit (f = 1), 5 km, 5.9 cm, LOS factor 1: dV = 0.059 pi 25e6 / 0.75 = 6.18e6 m^3
        dv = V.insar_min_volume(1.0, 5.0, 0.059, los=1.0)
        self.assertAlmostEqual(float(dv), 0.059 * np.pi * 25e6 / 0.75, delta=1)

    def test_consistent_with_mogi(self):
        # the returned volume must produce exactly the threshold LOS uplift at a pixel 3 km from the summit
        d, r, thr = 5.0, 3.0, 0.01
        f = (1 + (r / d) ** 2) ** -1.5
        dv = float(V.insar_min_volume(f, d, thr))
        _, _, uz = V.mogi_greens(r, 0.0, d)
        self.assertAlmostEqual(uz * dv * V.LOS_VERT, thr, places=12)

    def test_no_coherence_is_blind(self):
        self.assertTrue(np.isinf(V.insar_min_volume(0.0, 5.0, 0.01)))


class TestPrefilter(unittest.TestCase):
    def test_window_contains_100km_at_high_latitude(self):
        for lat in (0, 45, 70, 78, -78, 85):
            w = V.lon_window_deg(lat, 100.0)
            for la in (lat - 0.9, lat, lat + 0.9):
                la = max(min(la, 89.9), -89.9)
                _, _, d = V.GEOD.inv(0.0, lat, w, la)
                self.assertGreater(d / 1e3, 100.0, f"lat {lat}")


class TestNoise(unittest.TestCase):
    def test_white_noise_recovered(self):
        rng = np.random.default_rng(1)
        n = 730
        df = pd.DataFrame({"mjd": np.arange(60000, 60000 + n), "e": rng.normal(0, 0.0015, n),
                           "n": rng.normal(0, 0.0012, n), "u": 0.01 * np.arange(n) / 365 + rng.normal(0, 0.005, n)})
        df.loc[400:, "u"] += 0.05  # an offset must not inflate the estimate
        sig, npairs = V.white_noise_sigma(df)
        self.assertEqual(npairs, n - 1)
        np.testing.assert_allclose(sig, [0.0015, 0.0012, 0.005], rtol=0.08)

    def test_block_noise_white(self):
        # pure white noise: a 30-day mean has sigma / sqrt(30)
        rng = np.random.default_rng(2)
        n = 720
        df = pd.DataFrame({"mjd": np.arange(n), "e": rng.normal(0, 0.003, n), "n": rng.normal(0, 0.003, n),
                           "u": 0.02 * np.arange(n) / 365 + rng.normal(0, 0.003, n)})
        sig, npairs = V.block_noise_sigma(df, w=30)
        self.assertEqual(npairs, 23)
        np.testing.assert_allclose(sig, 0.003 / np.sqrt(30), rtol=0.45)

    def test_block_noise_removes_seasonal(self):
        # 4 mm annual vertical cycle on top of 2 mm white noise: the seasonal fit must keep it out of the 30-day sigma
        rng = np.random.default_rng(4)
        n = 730
        t = np.arange(n)
        u = 0.004 * np.sin(2 * np.pi * t / 365.25) + rng.normal(0, 0.002, n)
        df = pd.DataFrame({"mjd": 60000 + t, "e": rng.normal(0, 0.002, n), "n": rng.normal(0, 0.002, n), "u": u})
        with_fit, _ = V.block_noise_sigma(df, w=30)
        without, _ = V.block_noise_sigma(df, w=30, seasonal=False)
        self.assertLess(with_fit[2], 0.6 * without[2])
        self.assertLess(with_fit[2], 3 * 0.002 / np.sqrt(30))

    def test_block_noise_random_walk_larger(self):
        rng = np.random.default_rng(3)
        n = 720
        rw = np.cumsum(rng.normal(0, 0.0005, n))
        df = pd.DataFrame({"mjd": np.arange(n), "e": rw, "n": rw, "u": rw})
        sig, _ = V.block_noise_sigma(df, w=30)
        self.assertGreater(sig[0], 0.0005 / np.sqrt(30) * 3)

    def test_gaps_not_paired(self):
        df = pd.DataFrame({"mjd": [1, 2, 4, 5, 9], "e": [0, 1, 0, 1, 0], "n": 0, "u": 0})
        _, npairs = V.white_noise_sigma(df, min_pairs=1)
        self.assertEqual(npairs, 2)

    def test_parse_truncated_tail(self):
        text = ("ATION 99XXX 2024.0 junk\n"
                "P316 24SEP20 2024.7201 60573 2332 5 -124.1   1155  0.635662   4602670  0.175617   234  0.967697"
                "  0.0083 0.000906 0.001244 0.005015  0.115506 -0.110096 -0.440029  41.55 -124.08   234.96\n")
        df = V.parse_tenv3_tail(text)
        self.assertEqual(len(df), 1)
        self.assertAlmostEqual(df.loc[0, "e"], 1155.635662)
        self.assertAlmostEqual(df.loc[0, "ant"], 0.0083)
        self.assertEqual(df.loc[0, "mjd"], 60573)

    def test_block_noise_never_collapses_below_white_noise_floor(self):
        rng = np.random.default_rng(7)
        for nd in (300, 510):
            n_trials = 40
            ests = []
            for _ in range(n_trials):
                df = pd.DataFrame({
                    "mjd": 60000 + np.arange(nd),
                    "e": rng.normal(0, 0.001, nd),
                    "n": rng.normal(0, 0.001, nd),
                    "u": rng.normal(0, 0.004, nd),
                })
                wn, _ = V.white_noise_sigma(df, min_pairs=50)
                s30, _ = V.block_noise_sigma(df, w=30)
                self.assertTrue(np.all(s30 >= wn / np.sqrt(30) - 1e-12))
                ests.append(s30[2] / (0.004 / np.sqrt(30)))
            self.assertGreater(np.median(ests), 0.90)

    def test_block_noise_rejects_intrablock_antenna_toggle_and_spikes(self):
        rng = np.random.default_rng(11)
        n = 720
        u = rng.normal(0, 0.004, n)
        ant = np.zeros(n)
        u[60:64] -= 0.500
        ant[60:64] = 0.500
        u[150] += 0.300  # isolated single-day glitch
        df = pd.DataFrame({
            "mjd": 60000 + np.arange(n),
            "e": rng.normal(0, 0.001, n),
            "n": rng.normal(0, 0.001, n),
            "u": u,
            "ant": ant,
        })
        s30, _ = V.block_noise_sigma(df, w=30)
        self.assertLess(s30[2], 0.002)


class TestPipelineHelpers(unittest.TestCase):
    def test_insar_tile_name_and_antimeridian_intervals(self):
        import insar_coherence as ic

        self.assertEqual(ic.tile_name(0, 0), "N00E000")
        self.assertEqual(ic.tile_name(-1, -78), "S01W078")
        self.assertEqual(ic.tile_name(-16, -181), "S16E179")
        self.assertEqual(ic.tile_name(52, 180), "N52W180")
        self.assertEqual(ic._lon_intervals(-180.2, -179.8), [(179.8, 180.0), (-180.0, -179.8)])
        self.assertEqual(ic._lon_intervals(179.8, 180.2), [(179.8, 180.0), (-180.0, -179.8)])

    def test_mantel_haenszel_excludes_zero_variance_strata(self):
        import ranking as R

        df = pd.DataFrame([
            # Country A: genuine variation in outcome & exposure
            {"country": "A", "ae": True, "gnss_blind": True},
            {"country": "A", "ae": True, "gnss_blind": False},
            {"country": "A", "ae": False, "gnss_blind": True},
            {"country": "A", "ae": False, "gnss_blind": True},
            # Country B: has both ae and ~ae, but all blind -> v_i == 0 (non-informative)
            {"country": "B", "ae": True, "gnss_blind": True},
            {"country": "B", "ae": False, "gnss_blind": True},
            {"country": "B", "ae": False, "gnss_blind": True},
        ])
        res = R.mantel_haenszel(df, "ae", "gnss_blind", "country")
        self.assertEqual(res["informative_strata"], 1)

    def test_pop_bounds_separation(self):
        import ranking as R

        df = pd.DataFrame([
            {"lon": 0.0, "lat": 0.0, "pop30": 500_000.0},
            {"lon": 0.2, "lat": 0.0, "pop30": 400_000.0},  # ~22 km away -> overlaps
            {"lon": 2.0, "lat": 0.0, "pop30": 300_000.0},  # ~222 km away -> separate
        ])
        lo, up = R.pop_bounds(df, sep_km=60.0)
        self.assertEqual(up, 1_200_000.0)
        self.assertEqual(lo, 800_000.0)

    def test_export_web_num(self):
        import export_web as ew

        self.assertIsNone(ew.num(None))
        self.assertIsNone(ew.num(float("nan")))
        self.assertIsNone(ew.num(float("inf")))
        self.assertEqual(ew.num(12345.0, nd=3), 12300.0)

    def test_combine_insar_table_area_and_missing_season(self):
        from unittest.mock import patch
        import combine as C

        synth = pd.DataFrame([
            # Volcano 1: 4 ok seasons, but 2 seasons have frac_coh_5km < 0.10 -> insA worst becomes inf
            {"vnum": 1, "season": "winter", "status": "ok", "frac_coh_5km": 0.05, "nearest_coh_km": 0.1, "f_d2": 0.9, "f_d5": 0.99, "f_d10": 1.0},
            {"vnum": 1, "season": "spring", "status": "ok", "frac_coh_5km": 0.25, "nearest_coh_km": 0.2, "f_d2": 0.8, "f_d5": 0.95, "f_d10": 0.99},
            {"vnum": 1, "season": "summer", "status": "ok", "frac_coh_5km": 0.30, "nearest_coh_km": 0.3, "f_d2": 0.7, "f_d5": 0.90, "f_d10": 0.98},
            {"vnum": 1, "season": "fall",   "status": "ok", "frac_coh_5km": 0.08, "nearest_coh_km": 0.1, "f_d2": 0.9, "f_d5": 0.99, "f_d10": 1.0},
            # Volcano 2: only 3 ok seasons -> both ins & insA worst must be inf
            {"vnum": 2, "season": "winter", "status": "no_tile", "frac_coh_5km": np.nan, "nearest_coh_km": np.nan, "f_d2": np.nan, "f_d5": np.nan, "f_d10": np.nan},
            {"vnum": 2, "season": "spring", "status": "ok", "frac_coh_5km": 0.40, "nearest_coh_km": 0.5, "f_d2": 0.5, "f_d5": 0.8, "f_d10": 0.9},
            {"vnum": 2, "season": "summer", "status": "ok", "frac_coh_5km": 0.50, "nearest_coh_km": 0.4, "f_d2": 0.6, "f_d5": 0.85, "f_d10": 0.95},
            {"vnum": 2, "season": "fall",   "status": "ok", "frac_coh_5km": 0.60, "nearest_coh_km": 0.3, "f_d2": 0.7, "f_d5": 0.90, "f_d10": 0.97},
        ])
        with patch.object(C.pd, "read_csv", return_value=synth):
            out = C.insar_table().set_index("vnum")
        self.assertTrue(np.isfinite(out.loc[1, "insA_d5_t59_best"]))
        self.assertTrue(np.isinf(out.loc[1, "insA_d5_t59_worst"]))
        self.assertTrue(np.isfinite(out.loc[1, "ins_d5_t59_worst"]))
        self.assertTrue(np.isinf(out.loc[2, "ins_d5_t59_worst"]))
        self.assertTrue(np.isinf(out.loc[2, "insA_d5_t59_worst"]))

    def test_exported_volcanoes_geojson_matches_ranking(self):
        import json

        rk = pd.read_csv(V.DATA / "ranking.csv").set_index("vnum")
        vj = json.loads((V.ROOT / "docs" / "data" / "volcanoes.geojson").read_text())["features"]
        self.assertEqual(len(vj), len(rk))
        comb5_det = 0
        merapi_f5b = None
        for ft in vj:
            p = ft["properties"]
            a, b = p["g5"], p["ib5"]
            comb = b if a is None else (a if b is None else min(a, b))
            if comb is not None and comb <= 1e7:
                comb5_det += 1
            if p["n"] == "Merapi":
                merapi_f5b = p["f5b"]
        self.assertEqual(comb5_det, int((rk["main_d5"] <= 1e7).sum()))
        self.assertIsNotNone(merapi_f5b)
        self.assertLess(merapi_f5b, 0.10)


if __name__ == "__main__":
    unittest.main()
