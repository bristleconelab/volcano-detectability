"""High-threat, geodetically blind volcanoes: activity and exposure (GVP) against open-data detectability.

Inputs: docs/data/volcano_detectability_combined.csv, GVP E3 layer (population within 5/10/30/100 km, from GVP) and
the GVP Holocene eruptions layer (start year, VEI). A volcano is
  active      confirmed eruption that started in or after 1900
  exposed     population within 30 km >= 100,000
  GNSS-blind  no open GNSS station sees a 1e7 m^3 Mogi change at 5 km (1 sigma, daily)
  fully blind neither GNSS nor Sentinel-1 (5.9 cm, worst season) sees 1e7 m^3 at 5 km
Outputs: docs/data/ranking.csv (all subaerial volcanoes with the fields above), printed summary.
    python scripts/ranking.py
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import vdlib as V

DV = 1e7
POP30 = 100_000


def gvp_extra() -> pd.DataFrame:
    e3 = json.loads((V.RAW / "gvp_e3_volcanoes.json").read_text())["features"]
    pop = pd.DataFrame([{"vnum": f["properties"]["VolcanoNumber"],
                         **{f"pop{r}": f["properties"][f"Within_{r}km"] for r in (5, 10, 30, 100)}} for f in e3])
    er = pd.DataFrame([f["properties"] for f in
                       json.loads((V.RAW / "gvp_eruptions.json").read_text())["features"]])
    er = er[er["Activity_Type"] == "Confirmed Eruption"]
    g = er.groupby("Volcano_Number")
    stats = pd.DataFrame({
        "n_erupt_holocene": g.size(),
        "n_erupt_1900": g["StartDateYear"].apply(lambda s: int((s >= 1900).sum())),
        "vei_max": g["ExplosivityIndexMax"].max(),
        "last_start": g["StartDateYear"].max(),
    }).reset_index().rename(columns={"Volcano_Number": "vnum"})
    return pop.merge(stats, on="vnum", how="left")


def main():
    df = pd.read_csv(V.DATA / "volcano_detectability_combined.csv")
    df = df[df["elev_m"] >= 0].merge(gvp_extra(), on="vnum", how="left")
    df["n_erupt_1900"] = df["n_erupt_1900"].fillna(0).astype(int)
    df["active"] = df["n_erupt_1900"] > 0
    df["exposed"] = df["pop30"].fillna(0) >= POP30
    df["gnss_blind"] = df["dv_d5_k1"] > DV
    df["fully_blind"] = df["gnss_blind"] & (df["ins_d5_t59_worst"] > DV)
    df.to_csv(V.DATA / "ranking.csv", index=False, float_format="%.4g")

    pop_missing = df["pop30"].isna().sum()
    print(f"{len(df)} subaerial volcanoes; population missing for {pop_missing}")
    a = df[df["active"]]
    ae = a[a["exposed"]]
    print(f"active (eruption since 1900): {len(a)}; of these exposed (>= {POP30:,} within 30 km): {len(ae)}")
    for col in ("gnss_blind", "fully_blind"):
        b = ae[ae[col]]
        print(f"  {col}: {len(b)} of {len(ae)}; people within 30 km of them: {b['pop30'].sum() / 1e6:.1f} M "
              f"(all active-exposed: {ae['pop30'].sum() / 1e6:.1f} M)")
    top = ae[ae["gnss_blind"]].sort_values("pop30", ascending=False)
    cols = ["name", "country", "pop30", "pop100", "n_erupt_1900", "vei_max", "n_20km", "dv_d5_k1",
            "ins_d5_t59_best", "ins_d5_t59_worst"]
    print(top[cols].head(25).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\nby country (active, exposed, GNSS-blind):")
    print(top.groupby("country").agg(n=("name", "size"), pop30_M=("pop30", lambda s: round(s.sum() / 1e6, 1)))
          .sort_values("pop30_M", ascending=False).head(12).to_string())


if __name__ == "__main__":
    main()
