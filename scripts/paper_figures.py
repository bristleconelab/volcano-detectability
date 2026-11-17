"""Figures and the numbers quoted in the manuscript draft, from the tables in docs/data/ only (plus the GVP and NGL
inputs for the Shasta benchmark check).

  Fig. 2  cumulative number of subaerial volcanoes whose smallest detectable volume change is <= dV (5 km, Mogi):
          open GNSS (daily and 30-day noise) against Sentinel-1 (5.9 cm, best and worst season);
          (a) all subaerial volcanoes, (b) active (eruption since 1900) and exposed (>= 100,000 people within 30 km)
  Fig. 3  threat x blindness: (a) population within 30 km against GNSS dV_min for active volcanoes, coloured by
          eruptions since 1900; (b) GNSS-blind active and exposed volcanoes by country, people within 30 km
Outputs: docs/img/paper/fig2_cumulative.{png,pdf}, docs/img/paper/fig3_threat_blindness.{png,pdf}; every number used
in the text is printed.
    python scripts/paper_figures.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap, LogNorm  # noqa: E402

import vdlib as V  # noqa: E402

OUT = V.ROOT / "docs" / "img" / "paper"
DV = 1e7
GNSS_C, SAR_C, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
FLOOR_5KM = V.insar_min_volume(1.0, 5.0, 0.059).item()  # Sentinel-1 5.9 cm, coherent pixel at the summit

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "savefig.dpi": 300, "pdf.fonttype": 42})


def counts(v: pd.Series) -> str:
    v = v.to_numpy(float)
    fin = v[np.isfinite(v)]
    return (f"none within reach {np.isinf(v).sum():4d}  <=1e6 {(v <= 1e6).sum():4d}  <=1e7 {(v <= 1e7).sum():4d}  "
            f"<=1e8 {(v <= 1e8).sum():4d}  median(finite) {np.median(fin):.3g}")


def cumulative(ax, v, n, **kw):
    v = np.sort(v[np.isfinite(v)])
    ax.step(np.r_[1e4, v, 1e10], np.r_[0, np.arange(1, v.size + 1), v.size], where="post", **kw)


def fig2(rk: pd.DataFrame, ae: pd.DataFrame):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
    for ax, d, title in ((axs[0], rk, f"(a) All subaerial Holocene volcanoes (n = {len(rk):,})"),
                         (axs[1], ae, f"(b) Active and exposed volcanoes (n = {len(ae)})")):
        cumulative(ax, d["dv_d5_k1"].to_numpy(), len(d), color=GNSS_C, lw=1.6, label="GNSS, daily noise, 1$\\sigma$")
        cumulative(ax, d["dv30_d5_k1"].to_numpy(), len(d), color=GNSS_C, lw=1.0, ls=":",
                   label="GNSS, 30-day noise, 1$\\sigma$")
        cumulative(ax, d["ins_d5_t59_best"].to_numpy(), len(d), color=SAR_C, lw=1.6,
                   label="Sentinel-1, 5.9 cm, best season")
        cumulative(ax, d["ins_d5_t59_worst"].to_numpy(), len(d), color=SAR_C, lw=1.0, ls="--",
                   label="Sentinel-1, 5.9 cm, worst season")
        ax.axvline(DV, color=MUTED, lw=0.6, zorder=0)
        ax.axvspan(1e4, FLOOR_5KM, color="0.93", lw=0, zorder=-1)
        ax.axhline(len(d), color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=0)
        ax.set_xscale("log")
        ax.set_xlim(1e5, 1e9)
        ax.set_ylim(0, len(d) * 1.04)
        ax.set_xlabel("Volume change at 5 km depth, $\\Delta V$ (m$^3$)")
        ax.set_title(title, loc="left")
        ax.grid(axis="y", color="0.9", lw=0.5)
    axs[0].set_ylabel("Volcanoes with $\\Delta V_{min} \\leq \\Delta V$")
    axs[0].text(1.15e5, len(rk) * 0.97, "below the Sentinel-1\nfloor: GNSS only", va="top", fontsize=6.5, color=MUTED)
    axs[0].text(1.05e7, 25, "10$^7$ m$^3$", fontsize=6.5, color=MUTED)
    axs[0].legend(loc="lower right", frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig2_cumulative.{ext}")
    plt.close(fig)


def fig3(a: pd.DataFrame, blind: pd.DataFrame):
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw=dict(width_ratios=[1.45, 1]))
    ax = axs[0]
    x = a["pop30"].fillna(0).clip(lower=1e3).to_numpy()
    top = 3e9  # row for "no open station within 100 km"
    y = np.where(np.isfinite(a["dv_d5_k1"]), a["dv_d5_k1"], top)
    c = a["n_erupt_1900"].to_numpy()
    o = np.argsort(c)
    cmap = ListedColormap(plt.get_cmap("Oranges")(np.linspace(0.25, 1.0, 256)))  # keep 1 eruption visible
    sc = ax.scatter(x[o], y[o], c=c[o], s=11, cmap=cmap, norm=LogNorm(1, 100), edgecolors="0.25",
                    linewidths=0.3, zorder=3)
    ax.axhline(DV, color=MUTED, lw=0.6, zorder=1)
    ax.axhspan(DV, 1.5e9, color="0.95", lw=0, zorder=0)
    ax.axvline(1e5, color=MUTED, lw=0.6, ls="--", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e5, 6e9)
    ax.set_xlim(7e2, 2e7)
    ax.set_yticks([1e5, 1e6, 1e7, 1e8, 1e9, top])
    ax.set_yticklabels(["10$^5$", "10$^6$", "10$^7$", "10$^8$", "10$^9$", "none\n<100 km"])
    ax.set_xticks([1e3, 1e4, 1e5, 1e6, 1e7])
    ax.set_xticklabels(["$\\leq$10$^3$", "10$^4$", "10$^5$", "10$^6$", "10$^7$"])
    ax.set_xlabel("Population within 30 km (GVP)")
    ax.set_ylabel("GNSS $\\Delta V_{min}$ at 5 km (m$^3$)")
    ax.set_title(f"(a) Active volcanoes (n = {len(a)}); shaded: GNSS-blind at 10$^7$ m$^3$", loc="left")
    ax.text(1.1e5, 1.3e5, "exposed: $\\geq$100,000", fontsize=6.5, color=MUTED)
    cb = fig.colorbar(sc, ax=ax, pad=0.01, extend="max")
    cb.set_label("Confirmed eruptions since 1900")
    cb.outline.set_visible(False)

    ax = axs[1]
    g = blind.groupby("country").agg(n=("name", "size"), pop=("pop30", "sum")).sort_values("pop", ascending=False)
    small = g["pop"] < 1e6
    if small.any():
        g = pd.concat([g[~small], pd.DataFrame({"n": [g.loc[small, "n"].sum()], "pop": [g.loc[small, "pop"].sum()]},
                                                index=[f"{small.sum()} other countries"])])
    yy = np.arange(len(g))[::-1]
    ax.barh(yy, g["pop"] / 1e6, color=GNSS_C, height=0.7)
    for yi, (name, r) in zip(yy, g.iterrows()):
        ax.text(r["pop"] / 1e6 + 0.6, yi, f"{int(r['n'])}", va="center", fontsize=6.5, color=INK)
    ax.set_yticks(yy)
    ax.set_yticklabels(g.index)
    ax.set_xlabel("People within 30 km (millions)")
    ax.set_title(f"(b) GNSS-blind, active, exposed\n(n = {len(blind)}; label: volcanoes)", loc="left")
    ax.set_xlim(0, g["pop"].max() / 1e6 * 1.12)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig3_threat_blindness.{ext}")
    plt.close(fig)
    return g


def benchmark_check():
    from detectability import MAX_KM, load_noise
    noise = load_noise()
    vol = V.load_volcanoes()
    vs = vol[vol["name"].isin(("Shasta", "Medicine Lake"))]
    near = noise[noise["lat"].between(40.0, 43.0) & noise["lon"].between(-124.05, -119.7)]
    sig = near[["sig_e_mm", "sig_n_mm", "sig_u_mm"]].to_numpy() / 1e3
    pts = [(v.name, v.lon, v.lat) for v in vs.itertuples()] + [("midpoint", vs["lon"].mean(), vs["lat"].mean())]
    for name, lo, la in pts:
        x, y = V.local_km(near["lon"].to_numpy(), near["lat"].to_numpy(), lo, la)
        ok = np.hypot(x, y) <= MAX_KM
        dv, i, c = V.min_detectable_volume(x[ok], y[ok], sig[ok], 5.0)
        print(f"  {name:13s} dv_min {dv:.3g} m^3; deciding station {near['site'].to_numpy()[ok][i]} at "
              f"{np.hypot(x[ok][i], y[ok][i]):.1f} km, component {'ENU'[c]}; stations within 100 km {ok.sum()}")
    grid = pd.read_csv(V.DATA / "benchmark_shasta_medicine.csv")["dv_m3"]
    print(f"  grid ({len(grid)} cells): min {grid.min():.3g}, median {grid.median():.3g}, max {grid.max():.3g}; "
          f"share > 1e6: {(grid > 1e6).mean():.2f}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rk = pd.read_csv(V.DATA / "ranking.csv")
    noise = pd.read_csv(V.DATA / "station_noise.csv")
    src = pd.read_csv(V.DATA / "source_sensitivity.csv")
    ok = noise[noise["status"] == "ok"].drop_duplicates("site")
    a = rk[rk["active"]]
    ae = a[a["exposed"]]
    blind = ae[ae["gnss_blind"]]

    print("== Inputs")
    print(f"GVP Holocene subaerial volcanoes (elevation >= 0 m): {len(rk)}; NGL stations with noise: {len(ok)} "
          f"(30-day noise available: {ok['sig30_u_mm'].notna().sum()})")
    for k in ("e", "n", "u"):
        print(f"  median noise {k}: daily {ok[f'sig_{k}_mm'].median():.2f} mm (10-90 %: "
              f"{ok[f'sig_{k}_mm'].quantile(.1):.2f}-{ok[f'sig_{k}_mm'].quantile(.9):.2f}), "
              f"30-day {ok[f'sig30_{k}_mm'].median():.2f} mm")
    print(f"  days per station: median {ok['ndays'].median():.0f}")
    print(f"stations within 10/20/50 km: >=1 at {(rk.n_10km >= 1).sum()}/{(rk.n_20km >= 1).sum()}/"
          f"{(rk.n_50km >= 1).sum()} volcanoes; >=4 at {(rk.n_10km >= 4).sum()}/{(rk.n_20km >= 4).sum()}/"
          f"{(rk.n_50km >= 4).sum()}")

    print("== Global, 5 km unless stated (counts of subaerial volcanoes)")
    for col, lab in (("dv_d5_k1", "GNSS daily 1 sigma"), ("dv_d5_k3", "GNSS daily 3 sigma"),
                     ("dv30_d5_k1", "GNSS 30-day 1 sigma"), ("dv30_d5_k2", "GNSS 30-day 2 sigma"),
                     ("dv_d2_k1", "GNSS daily 1 sigma, 2 km"), ("dv_d10_k1", "GNSS daily 1 sigma, 10 km"),
                     ("ins_d5_t59_best", "S1 5.9 cm best season"), ("ins_d5_t59_worst", "S1 5.9 cm worst season"),
                     ("ins_d5_t10_best", "S1 1 cm best season"), ("ins_d5_t10_worst", "S1 1 cm worst season"),
                     ("ins_d2_t59_best", "S1 5.9 cm best, 2 km"), ("ins_d10_t59_best", "S1 5.9 cm best, 10 km"),
                     ("comb_d5", "GNSS or S1 best")):
        print(f"  {lab:28s} {counts(rk[col])}")
    worst_comb = np.minimum(rk["dv_d5_k1"], rk["ins_d5_t59_worst"])
    print(f"  {'GNSS or S1 worst':28s} {counts(worst_comb)}")
    print(f"S1 floor (5.9 cm, summit pixel, LOS cos 39 deg): 2 km {V.insar_min_volume(1.0, 2.0, 0.059).item():.3g}, "
          f"5 km {FLOOR_5KM:.3g}, 10 km {V.insar_min_volume(1.0, 10.0, 0.059).item():.3g}; "
          f"1 cm at 5 km {V.insar_min_volume(1.0, 5.0, 0.01).item():.3g} m^3")
    print(f"  smallest S1 value in the table at 5 km: {rk['ins_d5_t59_best'].min():.4g}")
    g, ib = rk["dv_d5_k1"], rk["ins_d5_t59_best"]
    print(f"  GNSS dV_min below the S1 floor: {(g < FLOOR_5KM).sum()}; seen at 1e7 by GNSS only "
          f"{((g <= DV) & (ib > DV)).sum()}, by S1 only {((g > DV) & (ib <= DV)).sum()}, by both "
          f"{((g <= DV) & (ib <= DV)).sum()}")
    nt = rk["n_seasons"].isna().sum()
    print(f"  no coherent pixel within 20 km in any season: {np.isinf(ib).sum()} (of which no coherence tile: {nt})")
    near = rk["nearest_coh_km_best"].replace(np.inf, np.nan)
    print(f"  nearest coherent pixel to the summit (best season): median {near.median() * 1e3:.0f} m; share of pixels "
          f"within 5 km that are coherent: median {rk.frac_coh_5km_best.median():.1%} (best season), "
          f"{rk.frac_coh_5km_worst.median():.1%} (worst season)")

    print("== Threat and exposure")
    print(f"active (confirmed eruption since 1900): {len(a)}; exposed (>=100k within 30 km): {rk.exposed.sum()}; "
          f"active and exposed: {len(ae)}; population missing for {rk.pop30.isna().sum()}")
    print(f"GNSS-blind at 1e7: all subaerial {rk.gnss_blind.sum()} ({rk.gnss_blind.mean():.1%}); active "
          f"{a.gnss_blind.sum()} ({a.gnss_blind.mean():.1%}); active and exposed {len(blind)} of {len(ae)} "
          f"({len(blind) / len(ae):.1%}); active with no station within 100 km {np.isinf(a.dv_d5_k1).sum()}")
    print(f"  people within 30 km (summed): active and exposed {ae.pop30.sum() / 1e6:.1f} M; "
          f"GNSS-blind ones {blind.pop30.sum() / 1e6:.1f} M")
    print(f"  GNSS-blind active and exposed: no station within 20 km {(blind.n_20km == 0).sum()}, within 50 km "
          f"{(blind.n_50km == 0).sum()}, within 100 km {np.isinf(blind.dv_d5_k1).sum()}")
    print(f"  ... Sentinel-1 sees 1e7 in best season {(blind.ins_d5_t59_best <= DV).sum()}, worst "
          f"{(blind.ins_d5_t59_worst <= DV).sum()}; median S1 dV_min (best) "
          f"{np.median(blind.ins_d5_t59_best[np.isfinite(blind.ins_d5_t59_best)]):.3g}; fully blind "
          f"{blind.fully_blind.sum()} ({blind[blind.fully_blind].pop30.sum() / 1e6:.1f} M); coherent share within 5 km: "
          f"median {blind.frac_coh_5km_best.median():.1%} best, {blind.frac_coh_5km_worst.median():.1%} worst season")
    print(f"  eruptions since 1900: at all subaerial {rk.n_erupt_1900.sum()}, at GNSS-blind volcanoes "
          f"{rk[rk.gnss_blind].n_erupt_1900.sum()}; at active-exposed {ae.n_erupt_1900.sum()}, at the blind ones "
          f"{blind.n_erupt_1900.sum()}")
    jp = ae[ae.country != "Japan"]
    print(f"  without Japan: {jp.gnss_blind.sum()} of {len(jp)} ({jp.gnss_blind.mean():.1%}); Japan "
          f"{(ae.country == 'Japan').sum()} active-exposed, blind {ae[ae.country == 'Japan'].gnss_blind.sum()}")
    ae_seen = ae[~ae.gnss_blind]
    print(f"  active-exposed seen by GNSS at 1e6: {(ae.dv_d5_k1 <= 1e6).sum()}; below the S1 floor: "
          f"{(ae.dv_d5_k1 < FLOOR_5KM).sum()}; seen ones by country: "
          f"{ae_seen.groupby('country').size().sort_values(ascending=False).head(5).to_dict()}")
    bins = pd.cut(rk.n_erupt_1900, [-1, 0, 3, 9, 10_000], labels=["0", "1-3", "4-9", ">=10"])
    tab = rk.groupby(bins, observed=True).agg(n=("vnum", "size"), blind=("gnss_blind", "mean"))
    print("  GNSS-blind share by eruptions since 1900: " +
          "; ".join(f"{i}: {r.blind:.0%} of {r.n}" for i, r in tab.iterrows()))
    pb = pd.cut(rk.pop30.fillna(0), [-1, 1e4, 1e5, 1e6, 1e10], labels=["<10k", "10k-100k", "100k-1M", ">=1M"])
    tab = rk.groupby(pb, observed=True).agg(n=("vnum", "size"), blind=("gnss_blind", "mean"))
    print("  GNSS-blind share by population within 30 km: " +
          "; ".join(f"{i}: {r.blind:.0%} of {r.n}" for i, r in tab.iterrows()))
    print("  top 10 GNSS-blind active-exposed by population:")
    for r in blind.sort_values("pop30", ascending=False).head(10).itertuples():
        print(f"    {r.name} ({r.country}): {r.pop30 / 1e6:.2f} M, {r.n_erupt_1900} eruptions since 1900, "
              f"GNSS {r.dv_d5_k1:.3g}, S1 best {r.ins_d5_t59_best:.3g}, worst {r.ins_d5_t59_worst:.3g}")
    print("  by country:")
    gb = blind.groupby("country").agg(n=("name", "size"), pop=("pop30", "sum"), er=("n_erupt_1900", "sum"))
    for c, r in gb.sort_values("pop", ascending=False).iterrows():
        print(f"    {c}: {r.n} volcanoes, {r['pop'] / 1e6:.2f} M, {r.er} eruptions since 1900")

    print("== Robustness (active and exposed, blind at 1e7 m^3, 5 km)")
    print(f"  30-day noise, 1 sigma: {(ae.dv30_d5_k1 > DV).sum()}; 30-day 2 sigma: {(ae.dv30_d5_k2 > DV).sum()}; "
          f"daily 3 sigma: {(ae.dv_d5_k3 > DV).sum()}; 2 km: {(ae.dv_d2_k1 > DV).sum()}; 10 km: "
          f"{(ae.dv_d10_k1 > DV).sum()}")
    sae = src[src.active & src.exposed]
    for col, lab in (("dv_mogi", "Mogi below summit"), ("dv_off_med", "Mogi 5 km offset, median azimuth"),
                     ("dv_off_max", "Mogi 5 km offset, worst azimuth"), ("dv_sill", "horizontal point sill")):
        b = sae[sae[col] > DV]
        print(f"  {lab:34s} blind {len(b):3d} of {len(sae)} ({b.pop30.sum() / 1e6:.1f} M); all subaerial seen at "
              f"1e7 {(src[col] <= DV).sum()}, at 1e6 {(src[col] <= 1e6).sum()}")
    r = (src.dv_mogi / src.dv_sill)
    r = r[np.isfinite(r)]
    print(f"  Mogi/sill dV ratio: median {r.median():.3f}, 5-95 % {r.quantile(.05):.3f}-{r.quantile(.95):.2f} "
          f"(n {len(r)})")

    print("== Benchmark, Mount Shasta to Medicine Lake (5 km, 1 sigma)")
    benchmark_check()

    fig2(rk, ae)
    g = fig3(a, blind)
    print("== Figures written to", OUT, "; Fig. 3b rows:", g.assign(pop=(g["pop"] / 1e6).round(2)).to_dict("index"))
    print(f"  Fig. 3a: active volcanoes with pop30 < 1e3 plotted at 1e3: {(a.pop30.fillna(0) < 1e3).sum()}")


if __name__ == "__main__":
    main()
