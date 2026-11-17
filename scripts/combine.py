"""Combine the GNSS and Sentinel-1 layers into one detectability table and a world map.

InSAR thresholds (cumulative line-of-sight displacement):
  T = 5.9 cm  operational: global detection threshold of the automated COMET/LiCSAR system, 2015-2020
              (Biggs et al. 2022, Bull. Volcanol., doi:10.1007/s00445-022-01608-x)
  T = 1.0 cm  optimistic: careful time-series analysis (rate uncertainties ~0.1-0.6 cm/yr, Albino et al. 2022)
Per volcano the best and worst season are kept (snow and vegetation make coherence seasonal).
Output: docs/data/volcano_detectability_combined.csv, docs/img/world_d5.png
    python scripts/combine.py
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

import vdlib as V  # noqa: E402

THRESH = {"t59": 0.059, "t10": 0.010}
DEPTHS = (2, 5, 10)
LAND = V.RAW / "ne_110m_land.geojson"


def insar_table() -> pd.DataFrame:
    c = pd.read_csv(V.DATA / "insar_coherence.csv")
    c = c[c["status"] == "ok"].drop_duplicates(["vnum", "season"])
    rows = []
    for vnum, g in c.groupby("vnum"):
        row = dict(vnum=vnum, n_seasons=len(g), frac_coh_5km_best=g["frac_coh_5km"].max(),
                   frac_coh_5km_worst=g["frac_coh_5km"].min(), nearest_coh_km_best=g["nearest_coh_km"].min())
        for d in DEPTHS:
            f = g[f"f_d{d}"].to_numpy()
            for name, t in THRESH.items():
                dv = V.insar_min_volume(f, d, t)
                row[f"ins_d{d}_{name}_best"] = float(dv.min())
                row[f"ins_d{d}_{name}_worst"] = float(dv.max()) if len(g) == 4 else np.inf
        rows.append(row)
    return pd.DataFrame(rows)


def world_map(df: pd.DataFrame, col: str, path):
    fig, ax = plt.subplots(figsize=(11, 5.4))
    if LAND.exists():
        for feat in json.loads(LAND.read_text())["features"]:
            geom = feat["geometry"]
            polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
            for poly in polys:
                xy = np.asarray(poly[0])
                ax.fill(xy[:, 0], xy[:, 1], color="0.9", lw=0)
    blind = ~np.isfinite(df[col])
    sc = ax.scatter(df.loc[~blind, "lon"], df.loc[~blind, "lat"], c=df.loc[~blind, col], s=10,
                    norm=LogNorm(1e5, 1e9), cmap="magma_r", edgecolors="k", linewidths=0.2, zorder=3)
    ax.scatter(df.loc[blind, "lon"], df.loc[blind, "lat"], marker="x", s=10, c="tab:blue", lw=0.7, zorder=3,
               label=f"no open GNSS station within 100 km and no coherent Sentinel-1 pixel ({blind.sum()})")
    cb = fig.colorbar(sc, ax=ax, extend="both", pad=0.01)
    cb.set_label("Smallest detectable volume change (m$^3$)")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 75)
    ax.set_aspect("equal")
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title("Open-data geodetic detectability of Holocene volcanoes: Mogi source at 5 km, "
                 "best of GNSS (1$\\sigma$ daily) and Sentinel-1 (5.9 cm, best season)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200)


def main():
    g = pd.read_csv(V.DATA / "volcano_detectability.csv")
    ins = insar_table()
    df = g.merge(ins, on="vnum", how="left")
    for d in DEPTHS:
        for name in THRESH:
            for s in ("best", "worst"):
                c = f"ins_d{d}_{name}_{s}"
                df[c] = df[c].fillna(np.inf)
        df[f"comb_d{d}"] = np.minimum(df[f"dv_d{d}_k1"], df[f"ins_d{d}_t59_best"])
    df.to_csv(V.DATA / "volcano_detectability_combined.csv", index=False, float_format="%.7g")
    sub = df[df["elev_m"] >= 0]
    print(f"{len(sub)} subaerial volcanoes; InSAR rows {sub['n_seasons'].notna().sum()}")
    for col in ("dv_d5_k1", "ins_d5_t59_best", "ins_d5_t59_worst", "ins_d5_t10_best", "comb_d5"):
        v = sub[col]
        print(f"{col:18s} inf {np.isinf(v).sum():4d}  <=1e6 {(v <= 1e6).sum():4d}  <=1e7 {(v <= 1e7).sum():4d}  "
              f"<=1e8 {(v <= 1e8).sum():4d}  median {np.median(v[np.isfinite(v)]):.3g}")
    world_map(sub, "comb_d5", V.ROOT / "docs" / "img" / "world_d5.png")


if __name__ == "__main__":
    main()
