"""Threat, exposure and open-data blindness, with the statistics asked for in review.

Inputs: docs/data/volcano_detectability_combined.csv, GVP E3 layer (population within 5/10/30/100 km) and the GVP
Holocene eruptions layer. Definitions:
  active        confirmed eruption that started in or after 1900
  exposed       population within 30 km >= 100,000
  gnss_blind    open GNSS network (SNR 5, 30-day means) cannot detect a 1e7 m^3 Mogi change at 5 km   [main]
  gnss_blind_1s first-pass USGS-style criterion (1 sigma daily, one component of one station)
  s1_sees_*     Sentinel-1 (5.9 cm, coherent area >= 10 % within 5 km) sees 1e7 m^3 in the best / worst season
  both_blind    gnss_blind and not s1_sees_worst
  field         GVP primary type is a volcanic field (no single summit reservoir; point source less meaningful)
Statistics: Fisher exact test of the blind share (active & exposed vs the rest), Mantel-Haenszel odds ratio and
test stratified by country, the same without Japan; population of the blind set as an upper bound (sum over
volcanoes) and a lower bound (largest set of volcanoes more than 60 km apart, so their 30 km circles cannot overlap).
Output: docs/data/ranking.csv, docs/data/ranking_stats.json, printed summary.
    python scripts/ranking.py
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

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
    st = pd.DataFrame({
        "n_erupt_holocene": g.size(),
        "n_erupt_1900": g["StartDateYear"].apply(lambda s: int((s >= 1900).sum())),
        "vei_max": g["ExplosivityIndexMax"].max(),
        "last_start": g["StartDateYear"].max(),
    }).reset_index().rename(columns={"Volcano_Number": "vnum"})
    return pop.merge(st, on="vnum", how="left")


def mantel_haenszel(df: pd.DataFrame, group: str, outcome: str, strata: str) -> dict:
    """Mantel-Haenszel common odds ratio and CMH chi-square (with continuity correction) over 2x2 strata."""
    num = den = a_sum = e_sum = v_sum = 0.0
    inf_strata = 0
    for _, s in df.groupby(strata):
        a = int((s[group] & s[outcome]).sum())
        b = int((s[group] & ~s[outcome]).sum())
        c = int((~s[group] & s[outcome]).sum())
        d = int((~s[group] & ~s[outcome]).sum())
        n = a + b + c + d
        if n < 2 or (a + b) == 0 or (c + d) == 0:
            continue
        r1, c1 = a + b, a + c
        v_i = r1 * (n - r1) * c1 * (n - c1) / (n ** 2 * (n - 1))
        if v_i > 0:
            inf_strata += 1
        num += a * d / n
        den += b * c / n
        a_sum += a
        e_sum += r1 * c1 / n
        v_sum += v_i
    chi = (abs(a_sum - e_sum) - 0.5) ** 2 / v_sum if v_sum > 0 else np.nan
    return {"or_mh": num / den if den > 0 else np.nan, "chi2": chi, "p": float(stats.chi2.sf(chi, 1)),
            "informative_strata": inf_strata}


def pop_bounds(sub: pd.DataFrame, sep_km: float = 60.0) -> tuple[float, float]:
    """Upper bound: sum of pop30. Lower bound: greedy set (largest first) of volcanoes more than sep_km apart."""
    upper = float(sub["pop30"].fillna(0).sum())
    s = sub.sort_values("pop30", ascending=False)
    chosen, lower = [], 0.0
    for r in s.itertuples():
        if all(V.GEOD.inv(r.lon, r.lat, c[0], c[1])[2] / 1e3 > sep_km for c in chosen):
            chosen.append((r.lon, r.lat))
            lower += 0 if np.isnan(r.pop30) else r.pop30
    return lower, upper


def share(d: pd.DataFrame, col: str) -> str:
    k, n = int(d[col].sum()), len(d)
    lo, hi = stats.binomtest(k, n).proportion_ci(method="wilson") if n else (np.nan, np.nan)
    return f"{k}/{n} = {100 * k / n:.1f} % (95 % CI {100 * lo:.1f}-{100 * hi:.1f})"


def main():
    df = pd.read_csv(V.DATA / "volcano_detectability_combined.csv")
    df = df[df["elev_m"] >= 0].merge(gvp_extra(), on="vnum", how="left")
    vt = V.load_volcanoes()[["vnum", "vtype"]]
    df = df.merge(vt, on="vnum", how="left")
    df["n_erupt_1900"] = df["n_erupt_1900"].fillna(0).astype(int)
    df["active"] = df["n_erupt_1900"] > 0
    df["exposed"] = df["pop30"].fillna(0) >= POP30
    df["ae"] = df["active"] & df["exposed"]
    df["gnss_blind"] = df["net_d5"].fillna(np.inf) > DV
    df["gnss_blind_1s"] = df["dv_d5_k1"].fillna(np.inf) > DV
    df["s1_sees_best"] = df["insA_d5_t59_best"] <= DV
    df["s1_sees_worst"] = df["insA_d5_t59_worst"] <= DV
    df["both_blind"] = df["gnss_blind"] & ~df["s1_sees_worst"]
    df["field"] = df["vtype"].fillna("").str.contains("field", case=False)
    df.to_csv(V.DATA / "ranking.csv", index=False, float_format="%.7g")

    ae = df[df["ae"]]
    out = {}
    print(f"{len(df)} subaerial volcanoes; population missing for {df['pop30'].isna().sum()}")
    print(f"active {df['active'].sum()}; exposed {df['exposed'].sum()}; active & exposed {len(ae)}")
    for col in ("gnss_blind", "gnss_blind_1s", "both_blind"):
        b = ae[ae[col]]
        lo, up = pop_bounds(b)
        print(f"  {col:13s} {share(ae, col)}; people within 30 km {lo / 1e6:.1f}-{up / 1e6:.1f} M")
        out[col] = dict(k=int(len(b)), n=int(len(ae)), pop_lower=lo, pop_upper=up)
    b = ae[ae["gnss_blind"]]
    print(f"  of the GNSS-blind: Sentinel-1 sees 1e7 in best season {int(b['s1_sees_best'].sum())}, "
          f"worst {int(b['s1_sees_worst'].sum())}; volcanic fields {int(b['field'].sum())}")
    print(f"  without volcanic fields: {share(ae[~ae['field']], 'gnss_blind')}")

    print("\nIs the blind share different for active & exposed volcanoes? (main GNSS criterion)")
    for name, d in (("all countries", df), ("without Japan", df[df["country"] != "Japan"])):
        tab = pd.crosstab(d["ae"], d["gnss_blind"]).reindex(index=[True, False], columns=[True, False], fill_value=0)
        orr, p = stats.fisher_exact(tab.to_numpy())
        mh = mantel_haenszel(d, "ae", "gnss_blind", "country")
        print(f"  {name:14s} active&exposed {share(d[d['ae']], 'gnss_blind')}; others {share(d[~d['ae']], 'gnss_blind')}; "
              f"Fisher OR {orr:.2f} p {p:.3f}; within-country MH OR {mh['or_mh']:.2f} p {mh['p']:.4f}")
        out[f"test_{name}"] = dict(fisher_or=orr, fisher_p=p, **mh)
    for name, grp in (("pop30 >= 100k (all)", df["exposed"]), ("eruption since 1900 (all)", df["active"])):
        mh = mantel_haenszel(df.assign(g=grp), "g", "gnss_blind", "country")
        print(f"  within-country MH, {name}: OR {mh['or_mh']:.2f} p {mh['p']:.4f}")
        out[f"mh_{name}"] = mh

    top = ae[ae["gnss_blind"]].sort_values("pop30", ascending=False)
    cols = ["name", "country", "pop30", "n_erupt_1900", "vei_max", "field", "n_20km", "net_d5", "dv_d5_k1",
            "insA_d5_t59_best", "insA_d5_t59_worst"]
    print("\n" + top[cols].head(25).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\nby country (active, exposed, GNSS-blind):")
    print(top.groupby("country").agg(n=("name", "size"), pop30_M=("pop30", lambda s: round(s.sum() / 1e6, 1)))
          .sort_values("pop30_M", ascending=False).head(12).to_string())
    eq = df["pop30"].notna() & df["pop100"].notna() & (df["pop30"] == df["pop100"])
    eq_pos = eq & (df["pop30"] > 0)
    out["pop_eq_summary"] = {
        "total": int(eq.sum()),
        "zero": int((eq & (df["pop30"] == 0)).sum()),
        "positive": int(eq_pos.sum()),
        "positive_volcanic_fields": int((eq_pos & df["field"]).sum()),
        "positive_islands_or_other": int((eq_pos & ~df["field"]).sum()),
    }
    (V.DATA / "ranking_stats.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
