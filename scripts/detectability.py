"""Minimum detectable Mogi volume change below every GVP Holocene volcano from open NGL GNSS stations.

Output docs/data/volcano_detectability.csv, one row per volcano:
  n_10km, n_20km, n_50km   active stations with a usable noise estimate within that distance
  dv_d{2,5,10}_k{1,3}      smallest detectable volume change (m^3) for a source at 2, 5, 10 km depth, K = 1 or 3 sigma
  best_site_d5, best_km_d5 station that decides dv_d5_k1 and its distance
Volcanoes with no station inside MAX_KM get inf ("geodetically blind" for open GNSS).
    python scripts/detectability.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import vdlib as V

MAX_KM = 100.0
DEPTHS = (2.0, 5.0, 10.0)
KS = (1.0, 3.0)


def load_noise() -> pd.DataFrame:
    n = pd.read_csv(V.DATA / "station_noise.csv")
    n = n[n["status"] == "ok"].drop_duplicates("site").reset_index(drop=True)
    return n


def per_volcano(vol: pd.DataFrame, noise: pd.DataFrame) -> pd.DataFrame:
    sig = noise[["sig_e_mm", "sig_n_mm", "sig_u_mm"]].to_numpy() / 1e3
    cols30 = ["sig30_e_mm", "sig30_n_mm", "sig30_u_mm"]
    sig30 = noise[cols30].to_numpy() / 1e3 if set(cols30) <= set(noise.columns) else None
    slat, slon = noise["lat"].to_numpy(), noise["lon"].to_numpy()
    rows = []
    for v in vol.itertuples():
        near = (np.abs(slat - v.lat) < 1.0) & (np.abs(((slon - v.lon + 180) % 360) - 180) < 2.5)
        idx = np.flatnonzero(near)
        row = dict(vnum=v.vnum, name=v.name, country=v.country, lat=v.lat, lon=v.lon, elev_m=v.elev_m,
                   last_eruption=v.last_eruption)
        if idx.size:
            x, y = V.local_km(slon[idx], slat[idx], v.lon, v.lat)
            dist = np.hypot(x, y)
            ok = dist <= MAX_KM
            idx, x, y, dist = idx[ok], x[ok], y[ok], dist[ok]
        else:
            x = y = dist = np.zeros(0)
        for r in (10, 20, 50):
            row[f"n_{r}km"] = int((dist <= r).sum())
        for d in DEPTHS:
            for k in KS:
                dv, i, _ = V.min_detectable_volume(x, y, sig[idx], d, k=k)
                row[f"dv_d{int(d)}_k{int(k)}"] = dv
                if d == 5.0 and k == 1.0:
                    row["best_site_d5"] = noise.loc[idx[i], "site"] if i >= 0 else ""
                    row["best_km_d5"] = round(float(dist[i]), 2) if i >= 0 else np.nan
            if sig30 is not None:
                s30 = sig30[idx]
                good = np.isfinite(s30).all(axis=1)
                for k in (1.0, 2.0):
                    row[f"dv30_d{int(d)}_k{int(k)}"] = V.min_detectable_volume(
                        x[good], y[good], s30[good], d, k=k)[0]
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    vol = V.load_volcanoes()
    noise = load_noise()
    out = per_volcano(vol, noise)
    out.to_csv(V.DATA / "volcano_detectability.csv", index=False, float_format="%.7g")
    sub = out[out["elev_m"] >= 0]
    print(f"{len(noise)} stations with noise; {len(out)} volcanoes ({len(sub)} with elevation >= 0 m)")
    for col in [c for c in ("dv_d5_k1", "dv_d5_k3", "dv30_d5_k1", "dv30_d5_k2") if c in sub]:
        dv = sub[col]
        print(f"{col}: no station within {MAX_KM:.0f} km {np.isinf(dv).sum()}; "
              f"<=1e5 {(dv <= 1e5).sum()}, <=1e6 {(dv <= 1e6).sum()}, <=1e7 {(dv <= 1e7).sum()}, "
              f"<=1e8 {(dv <= 1e8).sum()} m^3")


if __name__ == "__main__":
    main()
