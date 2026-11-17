"""Export compact GeoJSON for the web atlas (docs/atlas/): volcanoes with detectability fields, NGL stations, land.
    python scripts/export_web.py
"""
from __future__ import annotations

import json
import math

import pandas as pd
import requests

import vdlib as V

OUT = V.ROOT / "docs" / "atlas" / "data"
NE = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/{name}.geojson"


def num(x, nd=3):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    return float(f"{x:.{nd}g}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(V.DATA / "ranking.csv")
    nn = pd.read_csv(V.DATA / "national_networks.csv").set_index("vnum")["category"].to_dict()
    feats = []
    for r in df.itertuples():
        p = {"n": r.name, "c": r.country, "v": int(r.vnum), "p30": num(r.pop30), "p100": num(r.pop100),
             "e19": int(r.n_erupt_1900), "vei": num(r.vei_max, 2), "s20": int(r.n_20km), "s50": int(r.n_50km),
             "bs": r.best_site_d5 if isinstance(r.best_site_d5, str) else None, "bk": num(r.best_km_d5),
             "act": bool(r.active), "exp": bool(r.exposed), "nn": nn.get(int(r.vnum))}
        for d in (2, 5, 10):
            p[f"g{d}"] = num(getattr(r, f"dv_d{d}_k1"))
            p[f"gs{d}"] = num(getattr(r, f"dvsill_d{d}_k1"))
            p[f"ib{d}"] = num(getattr(r, f"ins_d{d}_t59_best"))
            p[f"iw{d}"] = num(getattr(r, f"ins_d{d}_t59_worst"))
            p[f"it{d}"] = num(getattr(r, f"ins_d{d}_t10_best"))
        feats.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(r.lon, 4), round(r.lat, 4)]},
                      "properties": p})
    (OUT / "volcanoes.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                       separators=(",", ":")))
    st = pd.read_csv(V.DATA / "station_noise.csv")
    st = st[st["status"] == "ok"]
    sf = [{"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(s.lon, 4), round(s.lat, 4)]},
           "properties": {"s": s.site, "e": s.sig_e_mm, "n": s.sig_n_mm, "u": s.sig_u_mm}} for s in st.itertuples()]
    (OUT / "stations.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": sf},
                                                      separators=(",", ":")))
    for name in ("ne_50m_land", "ne_50m_admin_0_boundary_lines_land"):
        path = OUT / f"{name}.geojson"
        if not path.exists():
            g = requests.get(NE.format(name=name), timeout=120).json()
            for f in g["features"]:
                f["properties"] = {}
            path.write_text(json.dumps(g, separators=(",", ":")))
    print({p.name: p.stat().st_size for p in OUT.iterdir()})


if __name__ == "__main__":
    main()
