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


def lon_window_deg(lat: float, km: float) -> float:
    """Longitude half-width (deg) that contains every point within km of latitude lat (prefilter, generous)."""
    return min(180.0, km / (111.0 * max(np.cos(np.radians(abs(lat) + km / 111.0)), 0.01)) + 0.1)


def mogi_greens(x_km, y_km, depth_km: float, nu: float = NU):
    """Surface displacement (m) per 1 m^3 of Mogi volume change at (0, 0, -depth): (ue, un, uu).

    u = (1 - nu) dV / pi * (x, y, d) / R^3 (Segall 2010, eq. 7.14)."""
    x, y, d = np.asarray(x_km) * 1e3, np.asarray(y_km) * 1e3, depth_km * 1e3
    r3 = (x ** 2 + y ** 2 + d ** 2) ** 1.5
    c = (1 - nu) / np.pi
    return c * x / r3, c * y / r3, c * d / r3


def sill_greens(x_km, y_km, depth_km: float):
    """Surface displacement (m) per 1 m^3 of opening of a horizontal point crack (sill) at depth d.

    Okada (1985) point tensile source with dip 0: the nu-dependent I-terms carry sin^2(dip) and vanish, leaving
    u = 3 dV d^2 / (2 pi R^5) * (x, y, d). Peak uplift is twice the Mogi value and the uplift volume equals dV."""
    x, y, d = np.asarray(x_km) * 1e3, np.asarray(y_km) * 1e3, depth_km * 1e3
    r5 = (x ** 2 + y ** 2 + d ** 2) ** 2.5
    c = 3 * d ** 2 / (2 * np.pi)
    return c * x / r5, c * y / r5, c * d / r5


def min_detectable_volume(x_km, y_km, sig_m, depth_km: float, k: float = 1.0, nu: float = NU,
                          source: str = "mogi"):
    """Smallest |dV| (m^3) detected by >= 1 component of >= 1 station; also index of the deciding station.

    x_km, y_km: station positions relative to the source (arrays, n). sig_m: (n, 3) white-noise sigma in metres
    for east, north, up. Returns (dv_min, i_station, i_component); inf and -1 when there are no stations."""
    x_km, y_km = np.atleast_1d(x_km), np.atleast_1d(y_km)
    if x_km.size == 0:
        return np.inf, -1, -1
    greens = sill_greens(x_km, y_km, depth_km) if source == "sill" else mogi_greens(x_km, y_km, depth_km, nu)
    g = np.abs(np.column_stack(greens))  # (n, 3) m per m^3
    with np.errstate(divide="ignore"):
        dv = np.where(g > 0, k * np.asarray(sig_m, float) / g, np.inf)
    i = int(np.argmin(dv))
    return float(dv.flat[i]), i // 3, i % 3


LOS_VERT = float(np.cos(np.radians(39.0)))  # Sentinel-1 IW mid-swath incidence ~39 deg; vertical-only LOS


def network_min_volume(x_km, y_km, sig_m, depth_km: float, snr: float = 5.0, nu: float = NU,
                       source: str = "mogi") -> float:
    """Smallest |dV| (m^3) that a matched filter over all stations and components detects at the given SNR.

    The template is the known source; with independent Gaussian noise sigma_d on each displacement estimate the
    optimal statistic has SNR = |dV| * sqrt(sum (G / sigma_d)^2), so dV_min = snr / sqrt(sum (G / sigma_d)^2).
    Pass sig_m as the noise of a displacement (difference of two epoch estimates, i.e. sqrt(2) x the epoch sigma).
    SNR 5 corresponds to a one-sided Gaussian false-alarm probability of about 3e-7 per test. inf without stations."""
    x_km, y_km = np.atleast_1d(x_km), np.atleast_1d(y_km)
    if x_km.size == 0:
        return np.inf
    greens = sill_greens(x_km, y_km, depth_km) if source == "sill" else mogi_greens(x_km, y_km, depth_km, nu)
    g = np.column_stack(greens)
    s = np.asarray(sig_m, float)
    ok = np.isfinite(s) & (s > 0)
    q = np.sum((g[ok] / s[ok]) ** 2)
    return float(snr / np.sqrt(q)) if q > 0 else np.inf


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
                         float(t[11]) + float(t[12]), float(t[13])))
        except ValueError:
            continue  # first line of a byte range is usually cut
    df = pd.DataFrame(rows, columns=["mjd", "e", "n", "u", "ant"]).drop_duplicates("mjd").sort_values("mjd")
    return df.reset_index(drop=True)


def _clean_block_rows(sub: pd.DataFrame, min_days: float, has_ant: bool) -> pd.DataFrame | None:
    """Drop intra-block antenna-height toggles and isolated daily spikes before averaging a block."""
    if len(sub) < min_days:
        return None
    if has_ant:
        ants = sub["ant"].to_numpy(float)
        if np.ptp(ants) > 1e-4:
            vals, counts = np.unique(np.round(ants, 4), return_counts=True)
            best = int(np.argmax(counts))
            if counts[best] < min_days:
                return None  # mixed antenna heights across >30% of the block
            sub = sub[np.round(ants, 4) == vals[best]]
    keep = np.ones(len(sub), dtype=bool)
    for c in ("e", "n", "u"):
        vals = sub[c].to_numpy(float)
        med = np.median(vals)
        mad = np.median(np.abs(vals - med))
        keep &= np.abs(vals - med) <= 5.0 * max(1.4826 * mad, 0.002)
    if int(keep.sum()) < min_days:
        return None
    return sub.iloc[keep]


def block_noise_sigma(df: pd.DataFrame, w: int = 30, min_cover: float = 0.7, min_pairs: int = 8,
                      seasonal: bool = True):
    """Scatter (m) of a w-day mean position, from differences of consecutive non-overlapping w-day block means.

    Unlike the daily estimate this includes coloured (flicker, random-walk) noise at the w-day scale, which is what
    limits detection of slow inflation. Blocks need >= min_cover of their days (after rejecting intra-block antenna
    height switches and >5-MAD daily spikes); differences regress out velocity and seasonal terms with an inlier guard,
    degrees-of-freedom correction, and the physical white-noise lower bound sig_daily / sqrt(w).
    Returns (sigma[3], number of block pairs); NaN when too few."""
    if df.empty:
        return np.full(3, np.nan), 0
    mjd = df["mjd"].to_numpy()
    blk = ((mjd - mjd.min()) // w).astype(int)
    has_ant = "ant" in df.columns
    min_days = min_cover * w
    blocks = []
    for b, sub in df.assign(blk=blk).groupby("blk"):
        clean = _clean_block_rows(sub, min_days, has_ant)
        if clean is not None:
            blocks.append((int(b), float(clean["mjd"].mean()), float(clean["e"].mean()),
                           float(clean["n"].mean()), float(clean["u"].mean())))
    if len(blocks) < 2:
        return np.full(3, np.nan), 0
    bdf = pd.DataFrame(blocks, columns=["blk", "tmid", "e", "n", "u"])
    ok = np.diff(bdf["blk"].to_numpy()) == 1
    npairs = int(ok.sum())
    if npairs < min_pairs:
        return np.full(3, np.nan), npairs
    # regress block-mean differences on a constant (velocity) and on the differences of annual and semi-annual
    # sinusoids, so that seasonal motion does not count as noise; differencing keeps the fit insensitive to steps
    tmid = bdf["tmid"].to_numpy() / 365.25
    t0, t1 = tmid[:-1][ok], tmid[1:][ok]
    cols = [np.ones_like(t0)]
    if npairs >= 2 * min_pairs and seasonal:
        for f in (1.0, 2.0):
            w_ = 2 * np.pi * f
            cols += [np.sin(w_ * t1) - np.sin(w_ * t0), np.cos(w_ * t1) - np.cos(w_ * t0)]
    A = np.column_stack(cols)
    p = A.shape[1]
    wn, _ = white_noise_sigma(df, min_pairs=min(50, max(10, len(df) // 4)))
    out = []
    for c_idx, c in enumerate(("e", "n", "u")):
        d = np.diff(bdf[c].to_numpy())[ok]
        keep = np.ones(d.size, bool)
        min_keep = min(d.size, max(p + 4, int(np.ceil(0.6 * d.size))))
        if d.size > p:
            for _ in range(3):  # least squares with guarded 3-MAD outlier rejection (steps, transients)
                coef, *_ = np.linalg.lstsq(A[keep], d[keep], rcond=None)
                r = d - A @ coef
                rk = r[keep]
                med_r = np.median(rk)
                mad = np.median(np.abs(rk - med_r))
                new_keep = np.abs(r - med_r) <= 3 * 1.4826 * max(mad, 1e-9)
                if int(new_keep.sum()) < min_keep:
                    break
                keep = new_keep
            coef, *_ = np.linalg.lstsq(A[keep], d[keep], rcond=None)
            r = d - A @ coef
        else:
            r = d
        rk = r[keep]
        k_eff = int(rk.size)
        dof_corr = np.sqrt(k_eff / max(k_eff - p, 1))
        est = 1.4826 * np.median(np.abs(rk - np.median(rk))) * dof_corr / np.sqrt(2)
        if np.isfinite(wn[c_idx]):
            est = max(est, wn[c_idx] / np.sqrt(w))
        out.append(float(est))
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
