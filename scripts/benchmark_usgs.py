"""Benchmark against USGS SIR 2024-5062-D (Montgomery-Brown et al., 2024, Fig. D4) and SIR 2026-5117-A.

USGS: Mogi sources on a grid at 5 km depth are inflated until the predicted displacement exceeds the GNSS white noise
at one station; between Mount Shasta and Medicine Lake a large area allows > 0.001 km^3 (1e6 m^3) to intrude undetected.
Here the same map is rebuilt from open NGL stations and our own per-station noise.
Outputs: docs/data/benchmark_shasta_medicine.csv (grid), docs/img/benchmark_shasta_medicine.{png,pdf} (the PDF is
Figure S1 of the paper), printed summary.
    python scripts/benchmark_usgs.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

import vdlib as V  # noqa: E402
from detectability import load_noise  # noqa: E402

BOX = dict(lon0=-122.75, lon1=-121.0, lat0=41.0, lat1=42.0)
DEPTH_KM, K, STEP_KM, MAX_KM = 5.0, 1.0, 1.0, 100.0
VOLCS = ("Shasta", "Medicine Lake")


def main():
    noise = load_noise()
    vol = V.load_volcanoes()
    vs = vol[vol["name"].isin(VOLCS)]
    lat_c = (BOX["lat0"] + BOX["lat1"]) / 2
    dlat = STEP_KM / 111.0
    dlon = STEP_KM / (111.0 * np.cos(np.radians(lat_c)))
    lats = np.arange(BOX["lat0"], BOX["lat1"] + 1e-9, dlat)
    lons = np.arange(BOX["lon0"], BOX["lon1"] + 1e-9, dlon)
    near = noise[(noise["lat"].between(BOX["lat0"] - 1, BOX["lat1"] + 1))
                 & (noise["lon"].between(BOX["lon0"] - 1.3, BOX["lon1"] + 1.3))].reset_index(drop=True)
    sig = near[["sig_e_mm", "sig_n_mm", "sig_u_mm"]].to_numpy() / 1e3
    grid = np.full((lats.size, lons.size), np.inf)
    for i, la in enumerate(lats):
        for j, lo in enumerate(lons):
            x, y = V.local_km(near["lon"].to_numpy(), near["lat"].to_numpy(), lo, la)
            ok = np.hypot(x, y) <= MAX_KM
            grid[i, j] = V.min_detectable_volume(x[ok], y[ok], sig[ok], DEPTH_KM, k=K)[0]
    LO, LA = np.meshgrid(lons, lats)
    pd.DataFrame({"lon": LO.ravel(), "lat": LA.ravel(), "dv_m3": grid.ravel()}).to_csv(
        V.DATA / "benchmark_shasta_medicine.csv", index=False, float_format="%.7g")

    print(f"{len(near)} NGL stations with noise in/near the box; grid {lats.size} x {lons.size} at {STEP_KM} km")
    print(f"grid dv (m^3): min {grid.min():.3g}, median {np.median(grid):.3g}, max {grid.max():.3g}; "
          f"fraction > 1e6: {(grid > 1e6).mean():.2f}")
    for v in vs.itertuples():
        x, y = V.local_km(near["lon"].to_numpy(), near["lat"].to_numpy(), v.lon, v.lat)
        ok = np.hypot(x, y) <= MAX_KM
        dv, ii, c = V.min_detectable_volume(x[ok], y[ok], sig[ok], DEPTH_KM, k=K)
        site = near.loc[np.flatnonzero(ok)[ii], "site"]
        print(f"  {v.name}: dv_min {dv:.3g} m^3 (deciding station {site}, "
              f"{np.hypot(x[ok][ii], y[ok][ii]):.1f} km, component {'ENU'[c]})")
    # midpoint between the two edifices, the area USGS flags
    mlo, mla = vs["lon"].mean(), vs["lat"].mean()
    x, y = V.local_km(near["lon"].to_numpy(), near["lat"].to_numpy(), mlo, mla)
    ok = np.hypot(x, y) <= MAX_KM
    print(f"  midpoint ({mla:.3f}, {mlo:.3f}): dv_min {V.min_detectable_volume(x[ok], y[ok], sig[ok], DEPTH_KM)[0]:.3g}")

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    pc = ax.pcolormesh(LO, LA, grid, norm=LogNorm(1e5, 1e8), cmap="magma_r", shading="auto")
    ax.contour(LO, LA, grid, levels=[1e6], colors="k", linewidths=0.8)
    ins = near[near["lat"].between(BOX["lat0"], BOX["lat1"]) & near["lon"].between(BOX["lon0"], BOX["lon1"])]
    ax.plot(ins["lon"], ins["lat"], "^", ms=5, mfc="w", mec="k", label="NGL station (active)")
    ax.plot(vs["lon"], vs["lat"], "*", ms=12, mfc="tab:cyan", mec="k", label="GVP volcano")
    for v in vs.itertuples():
        ax.annotate(v.name, (v.lon, v.lat), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlim(BOX["lon0"], BOX["lon1"])
    ax.set_ylim(BOX["lat0"], BOX["lat1"])
    ax.set_aspect(1 / np.cos(np.radians(lat_c)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    cb = fig.colorbar(pc, ax=ax, extend="both")
    cb.set_label(f"Smallest detectable volume change (m$^3$), {DEPTH_KM:g} km, {K:g}$\\sigma$")
    ax.legend(loc="lower left", fontsize=7)
    ax.set_title("Open NGL GNSS detectability, Mount Shasta to Medicine Lake (contour: 10$^6$ m$^3$)", fontsize=9)
    fig.tight_layout()
    img = V.ROOT / "docs" / "img"
    fig.savefig(img / "benchmark_shasta_medicine.png", dpi=200)
    # vector copy for the Supporting Information (Figure S1); no creation date, so reruns give identical files
    fig.savefig(img / "benchmark_shasta_medicine.pdf", metadata={"CreationDate": None})


if __name__ == "__main__":
    main()
