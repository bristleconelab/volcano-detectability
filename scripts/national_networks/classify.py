import collections
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import vdlib as V  # noqa: E402

D = Path(sys.argv[1])
OUT = Path(sys.argv[2])
OPEN = {"NGL", "UNAVCO", "PBO", "USGS"}

vol_csv = D / "blind_ae.csv"
if not vol_csv.exists():
    vol_csv = D / "the70.csv"
rows = list(csv.DictReader(vol_csv.open()))
world = json.loads((D / "dl_def.json").read_text())["deformation"]

# instrument types by station code gathered from all per-volcano files
itype: dict[str, set[str]] = collections.defaultdict(set)
for r in rows:
    j = json.loads((D / "vol" / f"{r['vnum']}.json").read_text())
    de = j["deformation"]
    code = {s["station_id"]: s["station_code"].strip() for s in de["station"]}
    for g in ("general", "tilt"):
        for inst in de["instrument"][g]:
            if inst["station_id"] in code:
                itype[code[inst["station_id"]]].add(inst["inst_type_abbre"])


def key(code: str) -> str:
    c = code.strip()
    m = re.match(r"NGL_(.{4})_\d+$", c) or re.match(r"UNAVCO_(.{4})$", c)
    if m:
        return m.group(1).rstrip("_")
    m = re.match(r"(.+)_GPS$", c)
    if m:
        return m.group(1)
    return c


def isgnss(code: str, ty: set[str]) -> bool:
    return bool(re.match(r"(NGL|UNAVCO)_", code)) or code.endswith("_GPS") or any("GPS" in t or "GNSS" in t for t in ty)


out = []
for r in rows:
    la, lo = float(r["lat"]), float(r["lon"])
    j = json.loads((D / "vol" / f"{r['vnum']}.json").read_text())
    cand = []
    for s in j["deformation"]["station"]:
        src = {s.get(f"datasource_{k}") for k in (1, 2, 3)} - {None}
        cand.append((s["station_code"].strip(), float(s["station_lat"]), float(s["station_lon"]), src))
    for w in world:
        src = {w["datasource1"], w["datasource2"], w["datasource3"]} - {None}
        cand.append((w["station_code"].strip(), float(w["station_lat"]), float(w["station_lon"]), src))
    g_open, g_nat, oth = (collections.defaultdict(set) for _ in range(3))
    for code, a, b, src in cand:
        _, _, d_m = V.GEOD.inv(lo, la, b, a)
        if d_m / 1e3 > 30.0:
            continue
        ty = itype.get(code, set())
        k = key(code)
        if isgnss(code, ty):
            (g_open if src <= OPEN or (src & {"NGL", "UNAVCO"} and not (src - OPEN - {"SNET"})) else g_nat)[k] |= src
        else:
            oth[k] |= src | {f"type:{'/'.join(sorted(ty))}" if ty else "type:unrecorded"}
    out.append(dict(
        vnum=r["vnum"], name=r["name"], country=r["country"], lat=la, lon=lo,
        gnss_open=sorted(g_open), gnss_nat=sorted(g_nat),
        nat_src=sorted(set().union(*g_nat.values())) if g_nat else [],
        other=len(oth),
        other_types=sorted({x for v in oth.values() for x in v if x.startswith("type:")}),
        other_src=sorted({x for v in oth.values() for x in v if not x.startswith("type:")}),
        seis=len(j["seismic"]["station"]),
    ))
OUT.write_text(json.dumps(out, indent=0))
for o in out:
    print(o["vnum"], o["name"][:16].ljust(16), "open", o["gnss_open"], "nat", o["gnss_nat"],
          o["nat_src"], "oth", o["other"], o["other_types"], o["other_src"])
