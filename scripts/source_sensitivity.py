"""How much do the GNSS results depend on the source model?

For every subaerial volcano (5 km depth, daily noise, K = 1):
  dv_mogi      Mogi below the summit (the main case)
  dv_sill      horizontal point sill below the summit (Okada 1985, dip 0)
  dv_off_med   Mogi offset 5 km from the summit, median over 8 azimuths
  dv_off_max   ... worst azimuth
Then recounts the headline: active (since 1900) and exposed (>= 100k within 30 km) volcanoes that open GNSS cannot see
at 1e7 m^3. Output docs/data/source_sensitivity.csv and a printed summary.
    python scripts/source_sensitivity.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import vdlib as V
from detectability import MAX_KM, load_noise

DEPTH, OFFSET_KM, AZ = 5.0, 5.0, np.radians(np.arange(0, 360, 45))


def main():
    noise = load_noise()
    sig = noise[["sig_e_mm", "sig_n_mm", "sig_u_mm"]].to_numpy() / 1e3
    slat, slon = noise["lat"].to_numpy(), noise["lon"].to_numpy()
    rk = pd.read_csv(V.DATA / "ranking.csv")
    rows = []
    for v in rk.itertuples():
        near = (np.abs(slat - v.lat) < 1.2) & (np.abs(((slon - v.lon + 180) % 360) - 180) < 3.0)
        idx = np.flatnonzero(near)
        x, y = V.local_km(slon[idx], slat[idx], v.lon, v.lat) if idx.size else (np.zeros(0), np.zeros(0))
        row = dict(vnum=v.vnum)
        ok = np.hypot(x, y) <= MAX_KM
        row["dv_mogi"] = V.min_detectable_volume(x[ok], y[ok], sig[idx][ok], DEPTH)[0]
        row["dv_sill"] = V.min_detectable_volume(x[ok], y[ok], sig[idx][ok], DEPTH, source="sill")[0]
        offs = []
        for a in AZ:
            xs, ys = x - OFFSET_KM * np.sin(a), y - OFFSET_KM * np.cos(a)
            okk = np.hypot(xs, ys) <= MAX_KM
            offs.append(V.min_detectable_volume(xs[okk], ys[okk], sig[idx][okk], DEPTH)[0])
        row["dv_off_med"] = float(np.median(offs))
        row["dv_off_max"] = float(np.max(offs))
        rows.append(row)
    s = pd.DataFrame(rows)
    out = rk[["vnum", "name", "country", "active", "exposed", "pop30"]].merge(s, on="vnum")
    out.to_csv(V.DATA / "source_sensitivity.csv", index=False, float_format="%.7g")
    chk = np.allclose(out["dv_mogi"].replace(np.inf, -1), rk["dv_d5_k1"].replace(np.inf, -1), rtol=1e-3)
    print(f"Mogi column reproduces dv_d5_k1: {chk}")
    ae = out[out["active"] & out["exposed"]]
    for col in ("dv_mogi", "dv_sill", "dv_off_med", "dv_off_max"):
        b = ae[ae[col] > 1e7]
        fin = out[col][np.isfinite(out[col])]
        print(f"{col:11s} active+exposed blind at 1e7: {len(b):3d}/{len(ae)} ({b['pop30'].sum() / 1e6:5.1f} M people); "
              f"all subaerial <=1e7: {(out[col] <= 1e7).sum()}; median finite {np.median(fin):.3g}")
    r = out["dv_mogi"] / out["dv_sill"]
    r = r[np.isfinite(r)]
    print(f"Mogi/sill dV ratio: median {r.median():.2f}, 5-95 % {r.quantile(.05):.2f}-{r.quantile(.95):.2f}")


if __name__ == "__main__":
    main()
