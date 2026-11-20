"""Export compact GeoJSON for the web atlas (docs/atlas/): volcanoes with detectability & context fields, NGL stations.
    python scripts/export_web.py
"""
from __future__ import annotations

import json
import math

import pandas as pd

import vdlib as V

OUT = V.ROOT / "docs" / "atlas" / "data"


def num(x, nd=4):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    return float(f"{x:.{nd}g}")


def rnum(x, dp=2):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    return round(float(x), dp)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(V.DATA / "ranking.csv")
    nn = pd.read_csv(V.DATA / "national_networks.csv").set_index("vnum")["category"].to_dict()
    st = pd.read_csv(V.DATA / "station_noise.csv")
    st = st[st["status"] == "ok"]
    st_coords = {str(s.site): (round(float(s.lon), 4), round(float(s.lat), 4)) for s in st.itertuples()}

    feats = []
    for r in df.itertuples():
        ns = r.nearest_site if isinstance(r.nearest_site, str) and r.nearest_site else None
        ns_lon, ns_lat = st_coords.get(ns, (None, None)) if ns else (None, None)
        elev = None if (r.elev_m is None or (isinstance(r.elev_m, float) and math.isnan(r.elev_m))) else int(round(r.elev_m))
        p = {
            "n": r.name, "c": r.country, "v": int(r.vnum),
            "el": elev, "vt": r.vtype if isinstance(r.vtype, str) else "",
            "p30": num(r.pop30, 4), "e19": int(r.n_erupt_1900), "vei": num(r.vei_max, 2),
            "s20": int(r.n30_20km), "s50": int(r.n30_50km),
            "sa20": int(r.n_20km), "sa50": int(r.n_50km),
            "ns": ns, "nk": num(r.nearest_km, 3), "ns_lon": ns_lon, "ns_lat": ns_lat,
            "nbs": r.net_best_site_d5 if isinstance(r.net_best_site_d5, str) and r.net_best_site_d5 else None,
            "nbk": num(r.net_best_km_d5, 3),
            "f5b": rnum(r.frac_coh_5km_best, 3),
            "f5w": rnum(r.frac_coh_5km_worst, 3),
            "act": bool(r.active), "exp": bool(r.exposed), "fld": bool(r.field),
            "nn": nn.get(int(r.vnum)),
        }
        for d in (2, 5, 10):
            p[f"g{d}"] = num(getattr(r, f"net_d{d}"), 4)          # main: network SNR 5, 30-day means
            p[f"gs{d}"] = num(getattr(r, f"netsill_d{d}"), 4)
            p[f"g1s{d}"] = num(getattr(r, f"dv_d{d}_k1"), 4)       # first pass: USGS-style 1 sigma daily
            p[f"ib{d}"] = num(getattr(r, f"insA_d{d}_t59_best"), 4)  # Sentinel-1 with coherent-area requirement
            p[f"iw{d}"] = num(getattr(r, f"insA_d{d}_t59_worst"), 4)
        feats.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(r.lon, 4), round(r.lat, 4)]},
                      "properties": p})
    (OUT / "volcanoes.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                       separators=(",", ":")))
    sf = [{"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(s.lon, 4), round(s.lat, 4)]},
           "properties": {
               "s": s.site,
               "e": rnum(s.sig_e_mm, 2), "n": rnum(s.sig_n_mm, 2), "u": rnum(s.sig_u_mm, 2),
               "e30": rnum(s.sig30_e_mm, 2), "n30": rnum(s.sig30_n_mm, 2), "u30": rnum(s.sig30_u_mm, 2),
           }} for s in st.itertuples()]
    (OUT / "stations.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": sf},
                                                      separators=(",", ":")))
    print({p.name: p.stat().st_size for p in OUT.iterdir()})


if __name__ == "__main__":
    main()
