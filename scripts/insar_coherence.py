"""Sentinel-1 coherence around every subaerial GVP Holocene volcano, per season.

Data: Global Seasonal Sentinel-1 Interferometric Coherence (Kellndorfer et al. 2022, Sci. Data 9:73,
doi:10.1038/s41597-022-01189-6), Dec 2019 - Nov 2020, 3 arcsec, 12-day VV coherence (x100, uint8), public on AWS.
Only a window of +-BOX_KM around each volcano is read (cloud-optimised GeoTIFF over HTTP), never whole tiles.

For each volcano and season the script stores, for sources at 2, 5 and 10 km depth, the best "visibility factor"
f = max over coherent pixels of (1 + (r/d)^2)^-1.5, the fraction of the Mogi peak uplift that the most favourable
coherent pixel would see (r = distance from the summit). A pixel is coherent when the 5 x 5 mean coherence >= COH_MIN.
Output: docs/data/insar_coherence.csv (resumable).
    python scripts/insar_coherence.py [--workers 8]
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import rasterio
from rasterio.merge import merge
from scipy.ndimage import uniform_filter

import vdlib as V

URL = ("/vsicurl/https://sentinel-1-global-coherence-earthbigdata.s3.us-west-2.amazonaws.com/data/tiles/"
       "{tile}/{tile}_{season}_{pol}_COH12.tif")
SEASONS = ("winter", "spring", "summer", "fall")
POLS = ("vv", "hh")
BOX_KM = 20.0
COH_MIN = 0.35
DEPTHS = (2.0, 5.0, 10.0)
OUT = V.DATA / "insar_coherence.csv"
FIELDS = ["vnum", "season", "status", "pol", "frac_coh_5km", "nearest_coh_km"] + [f"f_d{int(d)}" for d in DEPTHS]


def tile_name(lat_top: int, lon_left: int) -> str:
    """Tiles are named by their upper-left corner, e.g. N42W123 spans 41-42 N, 123-122 W."""
    lon_left = ((int(lon_left) + 180) % 360) - 180
    return f"{'N' if lat_top >= 0 else 'S'}{abs(lat_top):02d}{'E' if lon_left >= 0 else 'W'}{abs(lon_left):03d}"


def _lon_intervals(w: float, e: float) -> list[tuple[float, float]]:
    """Split a longitude interval [w, e] across the antimeridian into sub-intervals inside [-180, 180]."""
    if w < -180.0:
        return [(w + 360.0, 180.0), (-180.0, e)]
    if e > 180.0:
        return [(w, 180.0), (-180.0, e - 360.0)]
    return [(w, e)]


def _read_interval(w: float, e: float, s: float, n: float, lat0: float, season: str, pol: str):
    names = {tile_name(math.ceil(la), math.floor(lo))
             for la in (s, n, lat0) for lo in (w, max(w, e - 1e-6), 0.5 * (w + e))}
    srcs = []
    for t in sorted(names):
        try:
            srcs.append(rasterio.open(URL.format(tile=t, season=season, pol=pol)))
        except rasterio.errors.RasterioIOError:
            continue
    if not srcs:
        return None
    try:
        arr, tr = merge(srcs, bounds=(w, s, e, n), nodata=0)
    finally:
        for src in srcs:
            src.close()
    rows, cols = np.indices(arr[0].shape)
    lon = tr.c + (cols + 0.5) * tr.a
    lat = tr.f + (rows + 0.5) * tr.e
    return arr[0].astype(float) / 100.0, lon, lat


def read_window(lon0, lat0, season):
    dlat = BOX_KM / 111.0
    dlon = min(179.0, BOX_KM / (111.0 * max(math.cos(math.radians(lat0)), 0.05)))
    w, e, s, n = lon0 - dlon, lon0 + dlon, max(-89.9, lat0 - dlat), min(89.9, lat0 + dlat)
    intervals = _lon_intervals(w, e)
    for pol in POLS:
        parts = [res for iw, ie in intervals if (res := _read_interval(iw, ie, s, n, lat0, season, pol)) is not None]
        if not parts:
            continue
        if len(parts) == 1:
            coh, lon, lat = parts[0]
        else:
            h = min(p[0].shape[0] for p in parts)
            coh = np.hstack([p[0][:h] for p in parts])
            lon = np.hstack([p[1][:h] for p in parts])
            lat = np.hstack([p[2][:h] for p in parts])
        return coh, lon, lat, pol
    return None, None, None, ""


def analyse(v, season):
    coh, lon, lat, pol = read_window(v.lon, v.lat, season)
    if coh is None:
        return dict(vnum=v.vnum, season=season, status="no_tile", pol="")
    x, y = V.local_km(lon, lat, v.lon, v.lat)
    r = np.hypot(x, y)
    valid = (coh > 0) & (coh <= 1.0)
    smooth = uniform_filter(np.where(valid, coh, 0.0), 5) / np.maximum(uniform_filter(valid.astype(float), 5), 1e-9)
    good = valid & (smooth >= COH_MIN) & (r <= BOX_KM)
    out = dict(vnum=v.vnum, season=season, status="ok" if valid.any() else "empty", pol=pol)
    in5 = valid & (r <= 5)
    out["frac_coh_5km"] = round(float(good[r <= 5].sum() / max(in5.sum(), 1)), 3)
    out["nearest_coh_km"] = round(float(r[good].min()), 3) if good.any() else math.inf
    for d in DEPTHS:
        out[f"f_d{int(d)}"] = round(float(((1 + (r[good] / d) ** 2) ** -1.5).max()), 5) if good.any() else 0.0
    return out


def _load_existing() -> list[dict]:
    if not OUT.exists():
        return []
    with OUT.open() as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        if "pol" not in r or not r["pol"]:
            r["pol"] = "vv" if r["status"] == "ok" else ""
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--recheck-missing", action="store_true",
                    help="Re-run non-ok rows (e.g. Antarctic HH coverage) and antimeridian volcanoes (|lon|>179).")
    a = ap.parse_args()
    vol = V.load_volcanoes()
    vol = vol[vol["elev_m"] >= 0]
    existing = _load_existing()
    old_map = {(int(r["vnum"]), r["season"]): r for r in existing}
    if a.recheck_missing and existing:
        edge_vnums = set(vol.loc[vol["lon"].abs() > 179.0, "vnum"].astype(int))
        kept = [r for r in existing if r["status"] == "ok" and int(r["vnum"]) not in edge_vnums]
        done = {(int(r["vnum"]), r["season"]) for r in kept}
    else:
        kept = existing
        done = {(int(r["vnum"]), r["season"]) for r in existing}
    jobs = [(v, s) for v in vol.itertuples() for s in SEASONS if (v.vnum, s) not in done]
    if a.limit:
        jobs = jobs[: a.limit]
    print(f"{len(vol)} volcanoes x {len(SEASONS)} seasons; {len(done)} kept, {len(jobs)} to do", flush=True)
    results = list(kept)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                      GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="2"):
        with ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(analyse, v, s): (v.vnum, s) for v, s in jobs}
            for i, fu in enumerate(as_completed(futs), 1):
                try:
                    results.append(fu.result())
                except Exception as e:
                    key = futs[fu]
                    if key in old_map:
                        results.append(old_map[key])
                    print(f"  {key}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
                if i % 200 == 0:
                    print(f"  {i}/{len(jobs)}", flush=True)
    if a.recheck_missing or not existing:
        dedup = {(int(r["vnum"]), r["season"]): r for r in results}
        results = [dedup[k] for k in sorted(dedup, key=lambda x: (x[0], SEASONS.index(x[1])))]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(results)


if __name__ == "__main__":
    main()
