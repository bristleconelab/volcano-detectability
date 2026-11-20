"""Daily white-noise sigma for every active NGL station within MAX_KM of a GVP Holocene volcano.

Only the last ~260 kB of each .tenv3 file is fetched (HTTP Range), which holds about three years of daily solutions;
the newest 730 days are used. Two noise levels: daily white noise and the scatter of 30-day means. Results are appended to docs/data/station_noise.csv, so the run can be resumed.
    python scripts/station_noise.py [--workers 6]
"""
from __future__ import annotations

import argparse
import csv
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import requests

import vdlib as V

MAX_KM = 100.0          # Mogi signals beyond ~100 km are negligible for sources shallower than 10 km
ACTIVE_SINCE = "2026-08-01"  # station must have solutions in the last few months of the NGL archive
TAIL_BYTES = 260_000
WINDOW_DAYS = 730
OUT = V.DATA / "station_noise.csv"
FIELDS = ["site", "lat", "lon", "start", "end", "first_mjd", "last_mjd", "ndays", "npairs",
          "sig_e_mm", "sig_n_mm", "sig_u_mm", "npairs30", "sig30_e_mm", "sig30_n_mm", "sig30_u_mm",
          "status"]
_TLS = threading.local()


def _session() -> requests.Session:
    if not hasattr(_TLS, "session"):
        _TLS.session = requests.Session()
    return _TLS.session


def candidate_stations(active_since: str = ACTIVE_SINCE):
    vol = V.load_volcanoes()
    st = V.load_holdings()
    st = st[st["end"] >= active_since].reset_index(drop=True)
    keep = np.zeros(len(st), bool)
    for _, v in vol.iterrows():
        near = (np.abs(st["lat"] - v.lat) < 1.0) & (np.abs(((st["lon"] - v.lon + 180) % 360) - 180) < V.lon_window_deg(v.lat, MAX_KM))
        if not near.any():
            continue
        idx = np.flatnonzero(near.to_numpy())
        _, _, d = V.GEOD.inv(np.full(idx.size, v.lon), np.full(idx.size, v.lat),
                             st.loc[idx, "lon"].to_numpy(), st.loc[idx, "lat"].to_numpy())
        keep[idx[np.asarray(d) / 1e3 <= MAX_KM]] = True
    return st[keep].reset_index(drop=True)


def fetch(site: str, session: requests.Session | None = None) -> str:
    sess = session or _session()
    for attempt in range(4):
        try:
            r = sess.get(V.NGL_TENV3.format(site=site), headers={"Range": f"bytes=-{TAIL_BYTES}"}, timeout=60)
            if r.status_code in (200, 206):
                return r.text
            if r.status_code == 404:
                return ""
        except requests.RequestException:
            pass
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"{site}: download failed")


def process(row, session: requests.Session | None = None):
    text = fetch(row.site, session)
    base = dict(site=row.site, lat=row.lat, lon=row.lon, start=row.start, end=row.end)
    if not text:
        return {**base, "status": "missing"}
    df = V.parse_tenv3_tail(text)
    if df.empty:
        return {**base, "status": "unparsed"}
    df = df[df["mjd"] > df["mjd"].max() - WINDOW_DAYS]
    sig, npairs = V.white_noise_sigma(df)
    s30, n30 = V.block_noise_sigma(df, w=30)
    return {**base, "first_mjd": int(df["mjd"].min()), "last_mjd": int(df["mjd"].max()), "ndays": len(df),
            "npairs": npairs, "sig_e_mm": round(sig[0] * 1e3, 3), "sig_n_mm": round(sig[1] * 1e3, 3),
            "sig_u_mm": round(sig[2] * 1e3, 3), "npairs30": n30, "sig30_e_mm": round(s30[0] * 1e3, 3),
            "sig30_n_mm": round(s30[1] * 1e3, 3), "sig30_u_mm": round(s30[2] * 1e3, 3), "status": "ok" if np.isfinite(sig).all() else "short"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--active-since", default=ACTIVE_SINCE)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    st = candidate_stations(a.active_since)
    done = set()
    if OUT.exists() and not a.overwrite:
        with OUT.open() as f:
            done = {r["site"] for r in csv.DictReader(f)}
    todo = st[~st["site"].isin(done)]
    if a.limit:
        todo = todo.head(a.limit)
    print(f"{len(st)} candidate stations, {len(done)} done, {len(todo)} to fetch", flush=True)
    target = OUT.with_suffix(".csv.tmp") if a.overwrite else OUT
    new = a.overwrite or not OUT.exists()
    failed = []
    with target.open("w" if a.overwrite else "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        with ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(process, r): r.site for r in todo.itertuples()}
            for i, fu in enumerate(as_completed(futs), 1):
                try:
                    w.writerow(fu.result())
                except Exception as e:  # keep going; retried on next run or checked below
                    failed.append(futs[fu])
                    print(f"  {futs[fu]}: {e}", file=sys.stderr, flush=True)
                if i % 100 == 0:
                    f.flush()
                    print(f"  {i}/{len(todo)}", flush=True)
    if a.overwrite:
        if len(failed) > max(1, int(0.01 * len(todo))):
            raise RuntimeError(f"Aborted overwrite: {len(failed)}/{len(todo)} stations failed")
        target.replace(OUT)
        if failed:
            print(f"Retrying {len(failed)} transient failures...", flush=True)
            rem = st[st["site"].isin(failed)]
            with OUT.open("a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=FIELDS)
                for r in rem.itertuples():
                    w.writerow(process(r))


if __name__ == "__main__":
    main()
