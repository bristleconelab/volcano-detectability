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
       "{tile}/{tile}_{season}_vv_COH12.tif")
SEASONS = ("winter", "spring", "summer", "fall")
BOX_KM = 20.0
COH_MIN = 0.35
DEPTHS = (2.0, 5.0, 10.0)
OUT = V.DATA / "insar_coherence.csv"
FIELDS = ["vnum", "season", "status", "frac_coh_5km", "nearest_coh_km"] + [f"f_d{int(d)}" for d in DEPTHS]


def tile_name(lat_top: int, lon_left: int) -> str:
    """Tiles are named by their upper-left corner, e.g. N42W123 spans 41-42 N, 123-122 W."""
    return f"{'N' if lat_top >= 0 else 'S'}{abs(lat_top):02d}{'E' if lon_left >= 0 else 'W'}{abs(lon_left):03d}"


def read_window(lon0, lat0, season):
    dlat = BOX_KM / 111.0
    dlon = BOX_KM / (111.0 * max(math.cos(math.radians(lat0)), 0.05))
    w, e, s, n = lon0 - dlon, lon0 + dlon, lat0 - dlat, lat0 + dlat
    names = {tile_name(math.ceil(la), math.floor(lo)) for la in (s, n, lat0) for lo in (w, e, lon0)}
    srcs = []
    for t in sorted(names):
        try:
            srcs.append(rasterio.open(URL.format(tile=t, season=season)))
        except rasterio.errors.RasterioIOError:
            continue  # no tile: ocean or outside coverage
    if not srcs:
        return None, None
    try:
        arr, tr = merge(srcs, bounds=(w, s, e, n), nodata=0)
    finally:
        for src in srcs:
            src.close()
    return arr[0].astype(float) / 100.0, tr


def analyse(v, season):
    coh, tr = read_window(v.lon, v.lat, season)
    if coh is None:
        return dict(vnum=v.vnum, season=season, status="no_tile")
    rows, cols = np.indices(coh.shape)
    lon = tr.c + (cols + 0.5) * tr.a
    lat = tr.f + (rows + 0.5) * tr.e
    x, y = V.local_km(lon, lat, v.lon, v.lat)
    r = np.hypot(x, y)
    valid = coh > 0
    smooth = uniform_filter(np.where(valid, coh, 0.0), 5) / np.maximum(uniform_filter(valid.astype(float), 5), 1e-9)
    good = valid & (smooth >= COH_MIN) & (r <= BOX_KM)
    out = dict(vnum=v.vnum, season=season, status="ok" if valid.any() else "empty")
    in5 = valid & (r <= 5)
    out["frac_coh_5km"] = round(float(good[r <= 5].sum() / max(in5.sum(), 1)), 3)
    out["nearest_coh_km"] = round(float(r[good].min()), 3) if good.any() else math.inf
    for d in DEPTHS:
        out[f"f_d{int(d)}"] = round(float(((1 + (r[good] / d) ** 2) ** -1.5).max()), 5) if good.any() else 0.0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    vol = V.load_volcanoes()
    vol = vol[vol["elev_m"] >= 0]
    done = set()
    if OUT.exists():
        with OUT.open() as f:
            done = {(int(r["vnum"]), r["season"]) for r in csv.DictReader(f)}
    jobs = [(v, s) for v in vol.itertuples() for s in SEASONS if (v.vnum, s) not in done]
    if a.limit:
        jobs = jobs[: a.limit]
    print(f"{len(vol)} volcanoes x {len(SEASONS)} seasons; {len(done)} done, {len(jobs)} to do", flush=True)
    new = not OUT.exists()
    with OUT.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                          GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="2"):
            with ThreadPoolExecutor(a.workers) as ex:
                futs = {ex.submit(analyse, v, s): (v.vnum, s) for v, s in jobs}
                for i, fu in enumerate(as_completed(futs), 1):
                    try:
                        w.writerow(fu.result())
                    except Exception as e:  # retried on the next run
                        print(f"  {futs[fu]}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
                    if i % 200 == 0:
                        f.flush()
                        print(f"  {i}/{len(jobs)}", flush=True)


if __name__ == "__main__":
    main()
