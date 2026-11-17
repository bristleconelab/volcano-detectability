"""Shared helpers for volcano-detectability.

Detectability criterion (after Montgomery-Brown 2024, USGS SIR 2024-5062-D, and 2026, SIR 2026-5117-A): a Mogi point
source at depth d below a volcano is "detected" when its predicted displacement exceeds K times the daily white-noise
scatter on at least one component (east, north or up) of at least one GNSS station. Because displacement is linear in
the volume change dV, the smallest detectable volume has a closed form:

    dV_min = min over stations s and components k of  K * sigma_{s,k} / |G_{s,k}|

where G is the Mogi displacement per cubic metre. Inputs are open: Smithsonian GVP Holocene volcanoes and Nevada
Geodetic Laboratory (NGL) daily IGS20 .tenv3 series.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Geod

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DATA = ROOT / "docs" / "data"
NU = 0.25  # Poisson's ratio, the usual Mogi choice (Segall 2010)
GEOD = Geod(ellps="WGS84")
NGL_TENV3 = "https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/IGS20/{site}.tenv3"


# ------------------------------------------------------------------------------------------------ inputs
def load_volcanoes(path: Path = RAW / "gvp_holocene.json") -> pd.DataFrame:
    """GVP Holocene volcanoes from the WFS GeoJSON (geometry is lon, lat; properties repeat Latitude/Longitude)."""
    g = json.loads(path.read_text())
    rows = []
    for f in g["features"]:
        p = f["properties"]
        lon, lat = f["geometry"]["coordinates"]
        if abs(lat - p["Latitude"]) > 1e-6 or abs(lon - p["Longitude"]) > 1e-6:
            raise ValueError(f"coordinate order mismatch for {p['Volcano_Name']}")
        rows.append(dict(vnum=p["Volcano_Number"], name=p["Volcano_Name"], country=p["Country"], lat=lat, lon=lon,
                         elev_m=p["Elevation"], last_eruption=p["Last_Eruption_Year"],
                         vtype=p["Primary_Volcano_Type"], region=p["Region"]))
    return pd.DataFrame(rows)


def load_holdings(path: Path = RAW / "DataHoldings.txt") -> pd.DataFrame:
    """NGL station list: site, lat, lon (-180..180), first and last day with a solution, number of solutions."""
    rows = []
    for line in path.read_text().splitlines()[1:]:
        t = line.split()
        if len(t) < 11:
            continue
        lon = float(t[2])
        rows.append(dict(site=t[0], lat=float(t[1]), lon=lon - 360 if lon > 180 else lon,
                         start=t[7], end=t[8], nsol=int(t[10])))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ geometry
def local_km(lon, lat, lon0: float, lat0: float):
    """Azimuthal-equidistant (x east, y north) km about (lon0, lat0), from WGS84 geodesics."""
    lon, lat = np.broadcast_arrays(np.asarray(lon, float), np.asarray(lat, float))
    az, _, d = GEOD.inv(np.full(lon.shape, float(lon0)), np.full(lat.shape, float(lat0)), lon, lat)
    az, d = np.radians(az), np.asarray(d) / 1e3
    return d * np.sin(az), d * np.cos(az)


def mogi_greens(x_km, y_km, depth_km: float, nu: float = NU):
    """Surface displacement (m) per 1 m^3 of Mogi volume change at (0, 0, -depth): (ue, un, uu).

    u = (1 - nu) dV / pi * (x, y, d) / R^3 (Segall 2010, eq. 7.14)."""
    x, y, d = np.asarray(x_km) * 1e3, np.asarray(y_km) * 1e3, depth_km * 1e3
    r3 = (x ** 2 + y ** 2 + d ** 2) ** 1.5
    c = (1 - nu) / np.pi
    return c * x / r3, c * y / r3, c * d / r3


def min_detectable_volume(x_km, y_km, sig_m, depth_km: float, k: float = 1.0, nu: float = NU):
    """Smallest |dV| (m^3) detected by >= 1 component of >= 1 station; also index of the deciding station.

    x_km, y_km: station positions relative to the source (arrays, n). sig_m: (n, 3) white-noise sigma in metres
    for east, north, up. Returns (dv_min, i_station, i_component); inf and -1 when there are no stations."""
    x_km, y_km = np.atleast_1d(x_km), np.atleast_1d(y_km)
    if x_km.size == 0:
        return np.inf, -1, -1
    g = np.abs(np.column_stack(mogi_greens(x_km, y_km, depth_km, nu)))  # (n, 3) m per m^3
    with np.errstate(divide="ignore"):
        dv = np.where(g > 0, k * np.asarray(sig_m, float) / g, np.inf)
    i = int(np.argmin(dv))
    return float(dv.flat[i]), i // 3, i % 3


LOS_VERT = float(np.cos(np.radians(39.0)))  # Sentinel-1 IW mid-swath incidence ~39 deg; vertical-only LOS


def insar_min_volume(f, depth_km: float, threshold_m: float, los: float = LOS_VERT, nu: float = NU):
    """Smallest Mogi volume change (m^3) whose LOS uplift reaches threshold_m at the best coherent pixel.

    Peak uplift per m^3 is (1 - nu) / (pi d^2); a pixel at distance r sees the fraction f = (1 + (r/d)^2)^-1.5 of it,
    and the line of sight sees los times the vertical (horizontal motion ignored). f = 0 (no coherent pixel) -> inf."""
    f = np.asarray(f, float)
    peak = (1 - nu) / (np.pi * (depth_km * 1e3) ** 2)
    with np.errstate(divide="ignore"):
        return np.where(f > 0, threshold_m / (peak * los * f), np.inf)


# ------------------------------------------------------------------------------------------------ noise
def parse_tenv3_tail(text: str) -> pd.DataFrame:
    """Parse the rows of a (possibly byte-range-truncated) NGL .tenv3 file; displacements in metres."""
    rows = []
    for line in text.splitlines():
        t = line.split()
        if len(t) < 23 or t[0] == "site":
            continue
        try:
            rows.append((float(t[3]), float(t[7]) + float(t[8]), float(t[9]) + float(t[10]),
                         float(t[11]) + float(t[12])))
        except ValueError:
            continue  # first line of a byte range is usually cut
    df = pd.DataFrame(rows, columns=["mjd", "e", "n", "u"]).drop_duplicates("mjd").sort_values("mjd")
    return df.reset_index(drop=True)


def block_noise_sigma(df: pd.DataFrame, w: int = 30, min_cover: float = 0.7, min_pairs: int = 8):
    """Scatter (m) of a w-day mean position, from differences of consecutive non-overlapping w-day block means.

    Unlike the daily estimate this includes coloured (flicker, random-walk) noise at the w-day scale, which is what
    limits detection of slow inflation. Blocks need >= min_cover of their days; the median difference (a constant
    velocity) is removed by the robust scale. Returns (sigma[3], number of block pairs); NaN when too few."""
    mjd = df["mjd"].to_numpy()
    blk = ((mjd - mjd.min()) // w).astype(int)
    g = df.assign(blk=blk).groupby("blk")
    means = g[["e", "n", "u"]].mean()
    cnt = g.size()
    means = means[cnt >= min_cover * w]
    idx = means.index.to_numpy()
    ok = np.diff(idx) == 1
    npairs = int(ok.sum())
    if npairs < min_pairs:
        return np.full(3, np.nan), npairs
    out = []
    for c in ("e", "n", "u"):
        d = np.diff(means[c].to_numpy())[ok]
        out.append(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))
    return np.array(out), npairs


def white_noise_sigma(df: pd.DataFrame, min_pairs: int = 100) -> tuple[np.ndarray, int]:
    """Daily white-noise sigma (m) for e, n, u from consecutive-day differences.

    For white noise of sigma w, the difference of two consecutive days has sigma w*sqrt(2); flicker and random-walk
    noise and slow signals add little at a one-day lag. A robust scale (1.4826 * MAD) ignores offsets and outliers.
    Returns (sigma[3], number of consecutive-day pairs); NaN when fewer than min_pairs pairs."""
    mjd = df["mjd"].to_numpy()
    ok = np.diff(mjd) == 1
    npairs = int(ok.sum())
    if npairs < min_pairs:
        return np.full(3, np.nan), npairs
    out = []
    for c in ("e", "n", "u"):
        d = np.diff(df[c].to_numpy())[ok]
        mad = np.median(np.abs(d - np.median(d)))
        out.append(1.4826 * mad / np.sqrt(2))
    return np.array(out), npairs
