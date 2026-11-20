"""Fetch GVMID (Widiwijayanti et al. 2024) deformation stations for the active, exposed, GNSS-blind volcanoes.

GVMID web app JSON APIs (no login; non-commercial terms, cite the paper):
  https://wovodat.org/gvmid/php/download_api.php?type[]=deformation   all deformation stations
  https://wovodat.org/gvmid/php/api.php?vd_num=<GVP volcano number>   per-volcano instruments
Writes data/raw/gvmid/{dl_def.json, vol/<vnum>.json, the70.csv}. Then:
  python scripts/national_networks/classify.py data/raw/gvmid data/raw/gvmid/classified.json
  python scripts/national_networks/build.py data/raw/gvmid docs/data/national_networks.csv docs/data/station_noise.csv
build.py also holds the hand-checked literature evidence (each item tagged page-read or search-summary).
"""
import json
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import vdlib as V  # noqa: E402

D = V.RAW / "gvmid"
(D / "vol").mkdir(parents=True, exist_ok=True)
rk = pd.read_csv(V.DATA / "ranking.csv")
blind = rk[rk["active"] & rk["exposed"] & rk["gnss_blind"]][["vnum", "name", "country", "lat", "lon"]]
blind.to_csv(D / "blind_ae.csv", index=False)
s = requests.Session()
r = s.get("https://wovodat.org/gvmid/php/download_api.php", params={"type[]": "deformation"}, timeout=120)
(D / "dl_def.json").write_text(r.text)
for v in blind["vnum"]:
    p = D / "vol" / f"{v}.json"
    if not p.exists():
        p.write_text(s.get("https://wovodat.org/gvmid/php/api.php", params={"vd_num": v}, timeout=60).text)
print(len(blind), "volcanoes")
