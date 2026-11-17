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

    def test_no_station_is_blind(self):
        dv, i, c = V.min_detectable_volume([], [], np.zeros((0, 3)), 5.0)
        self.assertTrue(np.isinf(dv))
        self.assertEqual(i, -1)


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
        self.assertEqual(df.loc[0, "mjd"], 60573)


if __name__ == "__main__":
    unittest.main()
