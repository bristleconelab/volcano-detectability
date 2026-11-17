"""Figures and every number quoted in the manuscript draft, from the tables in docs/data/ (plus the GVP eruption list
for alternative "active" definitions, the NGL station list for the Japan check, and the GVP/NGL inputs for the Shasta
benchmark check).

Main criteria (after the internal reviews):
  GNSS        network matched filter over all open NGL stations within 100 km, SNR 5, displacement = difference of two
              seasonally cleaned 30-day means (noise sqrt(2) x sigma30)                       column net_d{2,5,10}
  Sentinel-1  5.9 cm LOS at a coherent pixel, a season counts only if >= 10 % of the pixels within 5 km of the
              summit are coherent                                                             insA_d5_t59_{best,worst}
  Comparison  USGS-style first pass: 1 sigma daily noise on one component of one station      dv_d5_k1
              Sentinel-1 single coherent pixel                                                ins_d5_t59_{best,worst}

  Fig. 1  world map of what open data detect at 1e7 m^3, 5 km (categories), with zooms on Mexico/Central America and
          Indonesia/Philippines; active and exposed volcanoes outlined
  Fig. 2  cumulative number of volcanoes whose smallest detectable volume change at 5 km is <= dV
          (a) all subaerial volcanoes, (b) active (eruption since 1900) and exposed (>= 100,000 people within 30 km)
  Fig. 3  (a) GNSS-blind share (95 % Wilson CI) by country group, active and exposed vs other volcanoes;
          (b) GNSS-blind active and exposed volcanoes by country with population bounds
Outputs: docs/img/paper/fig{1,2,3}_*.{png,pdf}; the numbers are printed (redirect to notes/paper_numbers_v2.txt).
    python scripts/paper_figures.py > notes/paper_numbers_v2.txt
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy import stats  # noqa: E402

import vdlib as V  # noqa: E402
from ranking import mantel_haenszel, pop_bounds  # noqa: E402

OUT = V.ROOT / "docs" / "img" / "paper"
DV = 1e7
GNSS_C, GNSS_L, SAR_C, SAR_L, NONE_C = "#1f4e9c", "#4f97d9", "#b5520b", "#e8962e", "#3d3d3d"
OTHER_C, INK, MUTED = "#8a8a8a", "#0b0b0b", "#52514e"
FLOOR_5KM = V.insar_min_volume(1.0, 5.0, 0.059).item()  # Sentinel-1 5.9 cm, coherent pixel at the summit
CENTRAL_AM = {"Mexico", "Mexico-Guatemala", "Guatemala", "Guatemala-El Salvador", "El Salvador", "Honduras",
              "Nicaragua", "Costa Rica", "Panama"}
GROUPS = ("Japan", "Indonesia", "Philippines", "Mexico and\nC. America", "Rest of\nworld")

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "savefig.dpi": 300, "pdf.fonttype": 42})


# ------------------------------------------------------------------------------------------------ helpers
def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    ci = stats.binomtest(int(k), int(n)).proportion_ci(method="wilson")
    return ci.low, ci.high


def share(mask: pd.Series) -> str:
    k, n = int(mask.sum()), len(mask)
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {100 * k / n:.1f} % (95 % CI {100 * lo:.1f}-{100 * hi:.1f})"


def mh_ci(df: pd.DataFrame, group: str, outcome: str, strata: str) -> dict:
    """Mantel-Haenszel OR (same strata rules as ranking.mantel_haenszel) with the Robins-Breslow-Greenland 95 % CI."""
    sr = ss = spr = spsqr = sqs = 0.0
    used = 0
    for _, s in df.groupby(strata):
        a = int((s[group] & s[outcome]).sum())
        b = int((s[group] & ~s[outcome]).sum())
        c = int((~s[group] & s[outcome]).sum())
        d = int((~s[group] & ~s[outcome]).sum())
        n = a + b + c + d
        if n < 2 or (a + b) == 0 or (c + d) == 0:
            continue
        used += 1
        r, q = a * d / n, b * c / n
        p_, q_ = (a + d) / n, (b + c) / n
        sr += r
        ss += q
        spr += p_ * r
        spsqr += p_ * q + q_ * r
        sqs += q_ * q
    orr = sr / ss
    var = spr / (2 * sr ** 2) + spsqr / (2 * sr * ss) + sqs / (2 * ss ** 2)
    se = np.sqrt(var)
    base = mantel_haenszel(df, group, outcome, strata)
    assert abs(base["or_mh"] - orr) < 1e-9, (base["or_mh"], orr)
    return dict(or_mh=orr, lo=float(np.exp(np.log(orr) - 1.96 * se)), hi=float(np.exp(np.log(orr) + 1.96 * se)),
                p=base["p"], strata=used)


def counts(v: pd.Series) -> str:
    v = v.fillna(np.inf).to_numpy(float)
    fin = v[np.isfinite(v)]
    return (f"no detection at any volume {np.isinf(v).sum():4d}  <=1e6 {(v <= 1e6).sum():4d}  "
            f"<=1e7 {(v <= 1e7).sum():4d}  <=1e8 {(v <= 1e8).sum():4d}  median(finite) {np.median(fin):.3g}")


def group_of(country: str) -> str:
    if country == "Japan":
        return GROUPS[0]
    if country == "Indonesia":
        return GROUPS[1]
    if country == "Philippines":
        return GROUPS[2]
    if country in CENTRAL_AM:
        return GROUPS[3]
    return GROUPS[4]


def mpop(x: float) -> str:
    return f"{x / 1e6:.1f}"


def bounds_str(sub: pd.DataFrame) -> str:
    lo, up = pop_bounds(sub)
    return f"{mpop(lo)}-{mpop(up)} M"


# ------------------------------------------------------------------------------------------------ figures
def land(ax):
    path = V.RAW / "ne_110m_land.geojson"
    if not path.exists():
        return
    for feat in json.loads(path.read_text())["features"]:
        geom = feat["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polys:
            xy = np.asarray(poly[0])
            ax.fill(xy[:, 0], xy[:, 1], color="0.91", lw=0, zorder=0)


def categories(rk: pd.DataFrame) -> pd.Series:
    net, sb, sw = rk["net_d5"].fillna(np.inf), rk["insA_d5_t59_best"], rk["insA_d5_t59_worst"]
    cat = np.select([net <= 1e6, net <= DV, sw <= DV, sb <= DV], [0, 1, 2, 3], default=4)
    return pd.Series(cat, index=rk.index)


CAT_STYLE = [
    ("Open GNSS: 10$^6$ m$^3$", GNSS_C, "o"),
    ("Open GNSS: 10$^7$ m$^3$", GNSS_L, "o"),
    ("Sentinel-1 only, every season", SAR_C, "s"),
    ("Sentinel-1 only, best season", SAR_L, "s"),
    ("Neither at 10$^7$ m$^3$", NONE_C, "x"),
]


def scatter_cats(ax, rk, cat, size):
    for i, (_, col, mk) in enumerate(CAT_STYLE):
        d = rk[cat == i]
        kw = dict(c=col, lw=0.9) if mk == "x" else dict(c=col, edgecolors="white", linewidths=0.25)
        ax.scatter(d["lon"], d["lat"], s=size, marker=mk, zorder=3 + (i == 4), **kw)
    a = rk[rk["ae"]]
    ax.scatter(a["lon"], a["lat"], s=size * 3.2, facecolors="none", edgecolors=INK, linewidths=0.6, zorder=5)


def fig1(rk: pd.DataFrame):
    cat = categories(rk)
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.18, wspace=0.12)
    ax = fig.add_subplot(gs[0, :])
    land(ax)
    scatter_cats(ax, rk, cat, 6)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 78)
    ax.set_aspect("equal")
    boxes = {"b": (-106.5, -82.0, 8.0, 22.5), "c": (94.0, 129.0, -11.5, 19.5)}
    for key, (x0, x1, y0, y1) in boxes.items():
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=MUTED, lw=0.6, ls="--", zorder=6))
        ax.text(x0, y1 + 1.5, f"({key})", fontsize=7, color=MUTED, zorder=6)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_xticklabels([f"{abs(x)}°{'W' if x < 0 else ('E' if x > 0 else '')}" for x in range(-180, 181, 60)])
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_yticklabels([f"{abs(y)}°{'S' if y < 0 else ('N' if y > 0 else '')}" for y in (-60, -30, 0, 30, 60)])
    ax.text(-178, 72, "(a)", fontsize=8, fontweight="bold")
    for j, (key, title) in enumerate((("b", "Mexico and Central America"), ("c", "Indonesia and Philippines"))):
        x0, x1, y0, y1 = boxes[key]
        a2 = fig.add_subplot(gs[1, j])
        land(a2)
        inbox = rk["lon"].between(x0, x1) & rk["lat"].between(y0, y1)
        scatter_cats(a2, rk[inbox], cat[inbox], 16)
        a2.set_xlim(x0, x1)
        a2.set_ylim(y0, y1)
        a2.set_aspect("equal")
        a2.set_title(f"({key}) {title}", loc="left")
        a2.tick_params(labelsize=6.5)
        a2.xaxis.set_major_locator(plt.MultipleLocator(10 if key == "c" else 5))
        a2.xaxis.set_major_formatter(lambda v, _: f"{abs(v):.0f}°{'W' if v < 0 else 'E'}")
        a2.yaxis.set_major_formatter(lambda v, _: f"{abs(v):.0f}°{'S' if v < 0 else ('N' if v > 0 else '')}")
        for s in ("top", "right"):
            a2.spines[s].set_visible(True)
    handles = [Line2D([], [], marker=mk, ls="none", color=col, markersize=5, markeredgewidth=0.9 if mk == "x" else 0.3,
                      markeredgecolor=col if mk == "x" else "white", label=lab) for lab, col, mk in CAT_STYLE]
    handles.append(Line2D([], [], marker="o", ls="none", markerfacecolor="none", markeredgecolor=INK, markersize=7,
                          label="Active and exposed (ring)"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7, columnspacing=1.2)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.98, bottom=0.12)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig1_map.{ext}")
    plt.close(fig)
    return cat


def cumulative(ax, v, **kw):
    v = np.sort(v[np.isfinite(v)])
    ax.step(np.r_[1e4, v, 1e10], np.r_[0, np.arange(1, v.size + 1), v.size], where="post", **kw)


def fig2(rk: pd.DataFrame, ae: pd.DataFrame):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), sharex=True)
    for ax, d, title in ((axs[0], rk, f"(a) All subaerial Holocene volcanoes (n = {len(rk):,})"),
                         (axs[1], ae, f"(b) Active and exposed volcanoes (n = {len(ae)})")):
        cumulative(ax, d["net_d5"].fillna(np.inf).to_numpy(), color=GNSS_C, lw=1.8,
                   label="Open GNSS network, SNR 5, 30-day means (main)")
        cumulative(ax, d["dv_d5_k1"].to_numpy(), color=GNSS_C, lw=1.0, ls=":",
                   label="Open GNSS, 1$\\sigma$ at one station (USGS-style)")
        cumulative(ax, d["insA_d5_t59_best"].to_numpy(), color=SAR_C, lw=1.8,
                   label="Sentinel-1, 5.9 cm, coherent area, best season")
        cumulative(ax, d["insA_d5_t59_worst"].to_numpy(), color=SAR_C, lw=1.0, ls="--",
                   label="Sentinel-1, 5.9 cm, coherent area, worst season")
        ax.axvline(DV, color=MUTED, lw=0.6, zorder=0)
        ax.axvspan(1e4, FLOOR_5KM, color="0.93", lw=0, zorder=-1)
        ax.axhline(len(d), color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=0)
        ax.set_xscale("log")
        ax.set_xlim(1e5, 1e9)
        ax.set_ylim(0, len(d) * 1.05)
        ax.set_xlabel("Volume change at 5 km depth, $\\Delta V$ (m$^3$)")
        ax.set_title(title, loc="left")
        ax.grid(axis="y", color="0.9", lw=0.5)
        ax.text(1.08e7, len(d) * 0.03, "10$^7$ m$^3$", fontsize=6.5, color=MUTED)
        ax.text(1.12e5, len(d) * 1.0, f"all {len(d):,}", fontsize=6.5, color=MUTED, va="bottom")
    axs[0].set_ylabel("Volcanoes with $\\Delta V_{\\mathrm{min}} \\leq \\Delta V$")
    axs[0].text(1.12e5, len(rk) * 0.90, "below the Sentinel-1\nfloor: GNSS only", va="top", fontsize=6.5,
                color=MUTED)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig2_cumulative.{ext}")
    plt.close(fig)


def fig3(rk: pd.DataFrame, blind: pd.DataFrame, mh: dict):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.5), gridspec_kw=dict(width_ratios=[1.15, 1]))
    ax = axs[0]
    rows = []
    for i, g in enumerate(GROUPS):
        d = rk[rk["cgroup"] == g]
        for j, (lab, m, col) in enumerate((("Active and exposed", d["ae"], GNSS_C),
                                           ("Other volcanoes", ~d["ae"], OTHER_C))):
            s = d[m]
            k, n = int(s["gnss_blind"].sum()), len(s)
            lo, hi = wilson(k, n)
            x = i + (-0.2 if j == 0 else 0.2)
            ax.errorbar(x, 100 * k / n, yerr=[[100 * (k / n - lo)], [100 * (hi - k / n)]], fmt="o" if j == 0 else "s",
                        color=col, ms=4.5, elinewidth=1.1, capsize=0, label=lab if i == 0 else None)
            above = j == 1 and hi < 0.9
            ax.text(x, 100 * hi + 3 if above else 100 * lo - 3, f"{k}/{n}", ha="center",
                    va="bottom" if above else "top", fontsize=5.8, color=MUTED)
            rows.append((g.replace("\n", " "), lab, k, n, lo, hi))
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(GROUPS, fontsize=6.8)
    ax.set_ylim(-14, 104)
    ax.set_yticks(range(0, 101, 20))
    ax.set_xlim(-0.5, len(GROUPS) - 0.5)
    ax.set_ylabel("Blind to open GNSS at 10$^7$ m$^3$ (%, 95 % CI)")
    ax.set_title("(a) Blind share within country groups", loc="left")
    ax.grid(axis="y", color="0.9", lw=0.5)
    ax.legend(loc="lower right", frameon=False, bbox_to_anchor=(1.0, 0.02))
    ax.text(0.99, 0.985, f"Within-country odds ratio\n{mh['or_mh']:.2f} (95 % CI {mh['lo']:.2f}-{mh['hi']:.2f})",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color=INK)

    ax = axs[1]
    gb = []
    for c, s in blind.groupby("country"):
        lo, up = pop_bounds(s)
        gb.append(dict(country=c, n=len(s), lo=lo, up=up))
    g = pd.DataFrame(gb).set_index("country").sort_values("up", ascending=False)
    small = g["up"] < 1.5e6
    if small.any():
        other = pd.DataFrame({"n": [g.loc[small, "n"].sum()], "lo": [g.loc[small, "lo"].sum()],
                              "up": [g.loc[small, "up"].sum()]}, index=[f"{small.sum()} other countries"])
        g = pd.concat([g[~small], other])
    yy = np.arange(len(g))[::-1]
    ax.barh(yy, g["up"] / 1e6, color=GNSS_L, height=0.72, label="Upper bound (sum over volcanoes)")
    ax.barh(yy, g["lo"] / 1e6, color=GNSS_C, height=0.72, label="Lower bound (volcanoes > 60 km apart)")
    for yi, (_, r) in zip(yy, g.iterrows()):
        ax.text(r["up"] / 1e6 + 0.5, yi, f"{int(r['n'])}", va="center", fontsize=6.5, color=INK)
    ax.set_yticks(yy)
    ax.set_yticklabels(g.index)
    ax.set_xlabel("People within 30 km (millions)")
    ax.set_title("(b) Blind in open GNSS data, by country", loc="left")
    ax.set_xlim(0, g["up"].max() / 1e6 * 1.1)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="center right", frameon=False, fontsize=6.5, bbox_to_anchor=(1.0, 0.42))
    ax.text(0.99, 0.62, "number = volcanoes", transform=ax.transAxes, ha="right", fontsize=6.5, color=MUTED)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig3_country.{ext}")
    plt.close(fig)
    return rows, g


# ------------------------------------------------------------------------------------------------ checks
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


def eruption_table() -> pd.DataFrame:
    er = pd.DataFrame([f["properties"] for f in json.loads((V.RAW / "gvp_eruptions.json").read_text())["features"]])
    return er[er["Activity_Type"] == "Confirmed Eruption"]


def japan_station_check(rk: pd.DataFrame, ok: pd.DataFrame):
    """How many active NGL stations within 100 km of the Japanese active and exposed volcanoes carry a 4-digit
    GEONET-style original name in DataHoldings.txt."""
    names = {}
    for line in (V.RAW / "DataHoldings.txt").read_text().splitlines()[1:]:
        t = line.split()
        if len(t) >= 11:
            names[t[0]] = t[11] if len(t) > 11 else ""
    jp = rk[rk["ae"] & (rk["country"] == "Japan")]
    near = set()
    for v in jp.itertuples():
        x, y = V.local_km(ok["lon"].to_numpy(), ok["lat"].to_numpy(), v.lon, v.lat)
        near |= set(ok["site"].to_numpy()[np.hypot(x, y) <= 100])
    geonet = [s for s in near if names.get(s, "").isdigit() and len(names.get(s, "")) == 4]
    print(f"  Japan: active NGL stations within 100 km of the 21 active and exposed volcanoes {len(near)}; with a "
          f"4-digit GEONET-style original name {len(geonet)} ({len(geonet) / max(len(near), 1):.0%})")


# ------------------------------------------------------------------------------------------------ main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rk = pd.read_csv(V.DATA / "ranking.csv")
    rk["net_d5"] = rk["net_d5"].fillna(np.inf)
    rk["cgroup"] = rk["country"].map(group_of)
    noise = pd.read_csv(V.DATA / "station_noise.csv")
    src = pd.read_csv(V.DATA / "source_sensitivity.csv")
    nn = pd.read_csv(V.DATA / "national_networks.csv")
    stj = json.loads((V.DATA / "ranking_stats.json").read_text())
    ok = noise[noise["status"] == "ok"].drop_duplicates("site")
    a = rk[rk["active"]]
    ae = rk[rk["ae"]]
    blind = ae[ae["gnss_blind"]]

    print("== 1 Inputs and noise")
    print(f"GVP Holocene subaerial volcanoes (elevation >= 0 m): {len(rk)}; active NGL stations with noise: {len(ok)}; "
          f"with seasonally cleaned 30-day noise: {ok['sig30_u_mm'].notna().sum()}")
    for k in ("e", "n", "u"):
        print(f"  median noise {k}: daily {ok[f'sig_{k}_mm'].median():.2f} mm; 30-day (seasonal terms removed) "
              f"{ok[f'sig30_{k}_mm'].median():.2f} mm (10-90 %: {ok[f'sig30_{k}_mm'].quantile(.1):.2f}-"
              f"{ok[f'sig30_{k}_mm'].quantile(.9):.2f})")
    npr = ok["npairs30"].median()
    print(f"  days per station: median {ok['ndays'].median():.0f}; 30-day block-difference pairs per station: median "
          f"{npr:.0f}; approx. relative standard error of a MAD-based sigma from n pairs, 1.17/sqrt(n): "
          f"{1.17 / np.sqrt(npr):.0%}")
    print(f"stations within 10/20/50 km: >=1 at {(rk.n_10km >= 1).sum()}/{(rk.n_20km >= 1).sum()}/"
          f"{(rk.n_50km >= 1).sum()} volcanoes; >=4 at {(rk.n_10km >= 4).sum()}/{(rk.n_20km >= 4).sum()}/"
          f"{(rk.n_50km >= 4).sum()}")
    pfa = stats.norm.sf(5.0)
    print(f"False alarms: SNR 5, one-sided Gaussian {pfa:.2g} per test; monthly tests at {len(rk)} volcanoes for one year "
          f"-> expected {pfa * 12 * len(rk):.2g} false alarms (white Gaussian noise assumed). First pass: 1 sigma is "
          f"exceeded by noise on one component {2 * stats.norm.sf(1):.0%} of the time, on at least one of 3 components "
          f"{1 - (1 - 2 * stats.norm.sf(1)) ** 3:.0%}")

    print("== 2 Global counts, 5 km unless stated")
    worst_main = np.minimum(rk["net_d5"], rk["insA_d5_t59_worst"])
    for v, lab in ((rk["net_d5"], "GNSS network SNR 5 (main)"), (rk["dv_d5_k1"], "GNSS 1 sigma one station"),
                   (rk["net_d2"].fillna(np.inf), "GNSS network, 2 km"), (rk["net_d10"].fillna(np.inf), "GNSS network, 10 km"),
                   (rk["insA_d5_t59_best"], "S1 area best season (main)"), (rk["insA_d5_t59_worst"], "S1 area worst season"),
                   (rk["ins_d5_t59_best"], "S1 single pixel best"), (rk["ins_d5_t59_worst"], "S1 single pixel worst"),
                   (rk["insA_d5_t10_best"], "S1 area 1 cm best"), (rk["insA_d5_t10_worst"], "S1 area 1 cm worst"),
                   (rk["main_d5"], "GNSS or S1 area, best season"), (worst_main, "GNSS or S1 area, worst season")):
        print(f"  {lab:32s} {counts(v)}")
    no_reach = np.isinf(rk["net_d5"])
    print(f"  no active open station within 100 km: {no_reach.sum()} (first-pass daily-noise set: "
          f"{np.isinf(rk['dv_d5_k1']).sum()})")
    fin = rk.loc[~no_reach, "net_d5"]
    print(f"  median GNSS network dV_min where a station exists: {fin.median():.3g} m^3")
    agree = ((rk["net_d5"] > DV) == (rk["n_50km"] == 0)).mean()
    print(f"  GNSS-blind at 1e7 identical to 'no active station within 50 km' for {agree:.1%}; within 20 km "
          f"{((rk['net_d5'] > DV) == (rk['n_20km'] == 0)).mean():.1%}")
    print(f"S1 floor (5.9 cm, summit pixel, LOS cos 39 deg): 2 km {V.insar_min_volume(1.0, 2.0, 0.059).item():.3g}, "
          f"5 km {FLOOR_5KM:.3g}, 10 km {V.insar_min_volume(1.0, 10.0, 0.059).item():.3g}; at 5 km with 1 cm "
          f"{V.insar_min_volume(1.0, 5.0, 0.01).item():.3g}, 3.9 cm {V.insar_min_volume(1.0, 5.0, 0.039).item():.3g}, "
          f"6.3 cm {V.insar_min_volume(1.0, 5.0, 0.063).item():.3g}")
    f30, f46 = (V.insar_min_volume(1.0, 5.0, 0.059, los=np.cos(np.radians(t))).item() for t in (30, 46))
    print(f"  floor range for incidence 30-46 deg: {f30:.3g}-{f46:.3g} m^3")
    print(f"  smallest S1 area value at 5 km: {rk['insA_d5_t59_best'].min():.4g}")
    print(f"  GNSS network dV_min below the S1 floor: {(rk['net_d5'] < FLOOR_5KM).sum()}; at 1e7 seen by GNSS only "
          f"{((rk.net_d5 <= DV) & (rk.insA_d5_t59_best > DV)).sum()}, by S1 area (best) only "
          f"{((rk.net_d5 > DV) & (rk.insA_d5_t59_best <= DV)).sum()}, by both "
          f"{((rk.net_d5 <= DV) & (rk.insA_d5_t59_best <= DV)).sum()}")
    print(f"  coherent share of pixels within 5 km: median {rk.frac_coh_5km_best.median():.1%} best season, "
          f"{rk.frac_coh_5km_worst.median():.1%} worst; seasons usable (>= 10 %) best/worst: "
          f"{(rk.frac_coh_5km_best >= 0.1).sum()}/{(rk.frac_coh_5km_worst >= 0.1).sum()}")

    print("== 3 Active and exposed")
    print(f"active (confirmed eruption starting 1900 or later): {len(a)}; exposed (>= 100k within 30 km): "
          f"{rk.exposed.sum()}; active and exposed: {len(ae)}; population missing: {rk.pop30.isna().sum()} "
          f"(active among them {a.pop30.isna().sum()})")
    print(f"  people within 30 km of the 118: {bounds_str(ae)} (lower-upper)")
    for col, lab in (("gnss_blind", "GNSS-blind (main)"), ("both_blind", "GNSS-blind and S1-area-blind worst season"),
                     ("gnss_blind_1s", "GNSS-blind first pass 1 sigma")):
        print(f"  {lab:45s} {share(ae[col])}; people {bounds_str(ae[ae[col]])}; "
              f"ranking_stats.json {stj[col]['k']} ({mpop(stj[col]['pop_lower'])}-{mpop(stj[col]['pop_upper'])} M)")
    bb = ae[ae["gnss_blind"] & ~ae["s1_sees_best"]]
    print(f"  GNSS-blind and S1-area-blind in the best season too: {share(ae['gnss_blind'] & ~ae['s1_sees_best'])}; "
          f"people {bounds_str(bb)}")
    print(f"  all 70 first-pass blind are also main-blind: {bool((ae.gnss_blind_1s <= ae.gnss_blind).all())}; added "
          f"by the main criterion: {', '.join(ae[ae.gnss_blind & ~ae.gnss_blind_1s].name)}")
    print(f"  of the {len(blind)}: S1 area sees 1e7 best {int(blind.s1_sees_best.sum())}, worst "
          f"{int(blind.s1_sees_worst.sum())}; S1 single pixel best {(blind.ins_d5_t59_best <= DV).sum()}, worst "
          f"{(blind.ins_d5_t59_worst <= DV).sum()}; coherent share within 5 km median "
          f"{blind.frac_coh_5km_best.median():.1%} best, {blind.frac_coh_5km_worst.median():.1%} worst")
    print(f"  of the {len(blind)}: no station within 20 km {(blind.n_20km == 0).sum()}, 50 km {(blind.n_50km == 0).sum()},"
          f" 100 km {np.isinf(blind.net_d5).sum()}")
    print(f"  volcanic fields among the {len(blind)}: {int(blind.field.sum())} ({', '.join(blind[blind.field].name)}); "
          f"without fields {share(ae[~ae.field]['gnss_blind'])}; people {bounds_str(blind[~blind.field])}")
    print(f"  active and exposed seen by GNSS at 1e6: {(ae.net_d5 <= 1e6).sum()}; at 1e8 blind: "
          f"{(ae.net_d5 > 1e8).sum()}; at 1e6 blind: {(ae.net_d5 > 1e6).sum()}")
    print("  both-blind (GNSS main and S1 area worst season) by volume and depth:")
    for d in (2, 5, 10):
        net = ae[f"net_d{d}"].fillna(np.inf)
        sw = ae[f"insA_d{d}_t59_worst"]
        sb = ae[f"insA_d{d}_t59_best"]
        cells = ", ".join(f"{v:.0e}: {((net > v) & (sw > v)).sum()} (best season {((net > v) & (sb > v)).sum()})"
                          for v in (3e6, 1e7, 3e7, 1e8))
        print(f"    depth {d:2d} km, GNSS-blind at 1e7 {(net > DV).sum()}; both-blind {cells}")
    print(f"  eruptions since 1900: at the 118 {ae.n_erupt_1900.sum()}, at the {len(blind)} {blind.n_erupt_1900.sum()}")
    print(f"  GVP primary type of the {len(blind)}: {blind.vtype.value_counts().to_dict()}")
    jp = ae[ae.country != "Japan"]
    print(f"  Japan: {(ae.country == 'Japan').sum()} active and exposed, blind {ae[ae.country == 'Japan'].gnss_blind.sum()};"
          f" without Japan {share(jp['gnss_blind'])}")
    japan_station_check(rk, ok)
    print("  top 8 GNSS-blind active and exposed by population within 30 km:")
    for r in blind.sort_values("pop30", ascending=False).head(8).itertuples():
        print(f"    {r.name} ({r.country}{', volcanic field' if r.field else ''}): {r.pop30 / 1e6:.2f} M, "
              f"{r.n_erupt_1900} eruptions since 1900, GNSS {r.net_d5:.3g}, S1 area best {r.insA_d5_t59_best:.3g}, "
              f"worst {r.insA_d5_t59_worst:.3g}")

    print("== 4 Statistics (main GNSS criterion, blind at 1e7 m^3, 5 km)")
    rest = rk[~rk.ae]
    for name, d in (("all countries", rk), ("without Japan", rk[rk.country != "Japan"])):
        tab = pd.crosstab(d["ae"], d["gnss_blind"]).reindex(index=[True, False], columns=[True, False], fill_value=0)
        orr, p = stats.fisher_exact(tab.to_numpy())
        print(f"  pooled, {name}: active and exposed {share(d[d.ae]['gnss_blind'])}; others "
              f"{share(d[~d.ae]['gnss_blind'])}; Fisher OR {orr:.2f}, p {p:.3f}")
    for name, d in (("all countries", rk), ("without Japan", rk[rk.country != "Japan"])):
        m = mh_ci(d, "ae", "gnss_blind", "country")
        print(f"  within-country MH (strata = GVP country), active and exposed vs others, {name}: OR {m['or_mh']:.2f} "
              f"(95 % CI {m['lo']:.2f}-{m['hi']:.2f}), CMH p {m['p']:.1e}, informative strata {m['strata']}")
    mh_main = mh_ci(rk, "ae", "gnss_blind", "country")
    for g, lab in (("exposed", "exposed (>= 100k within 30 km) vs not"), ("active", "active (eruption since 1900) vs not")):
        m = mh_ci(rk, g, "gnss_blind", "country")
        print(f"  within-country MH, {lab}: OR {m['or_mh']:.2f} (95 % CI {m['lo']:.2f}-{m['hi']:.2f}), p {m['p']:.2g}")
    rk["freq4"] = rk.n_erupt_1900 >= 4
    m = mh_ci(rk, "freq4", "gnss_blind", "country")
    print(f"  within-country MH, >= 4 eruptions since 1900 vs fewer: OR {m['or_mh']:.2f} (95 % CI {m['lo']:.2f}-"
          f"{m['hi']:.2f}), p {m['p']:.2g}")
    m = mh_ci(rk, "ae", "gnss_blind_1s", "country")
    print(f"  same test, first-pass 1 sigma criterion: OR {m['or_mh']:.2f} (95 % CI {m['lo']:.2f}-{m['hi']:.2f}), "
          f"p {m['p']:.2g}")
    print(f"  rest of subaerial (not active and exposed): {share(rest['gnss_blind'])}")
    bins = pd.cut(rk.n_erupt_1900, [-1, 0, 3, 9, 10_000], labels=["0", "1-3", "4-9", ">=10"])
    print("  blind share by eruptions since 1900: " + "; ".join(
        f"{i}: {share(s['gnss_blind'])}" for i, s in rk.groupby(bins, observed=True)))
    pb = pd.cut(rk.pop30.fillna(0), [-1, 1e4, 1e5, 1e6, 1e10], labels=["<10k (incl. missing)", "10k-100k", "100k-1M", ">=1M"])
    print("  blind share by population within 30 km: " + "; ".join(
        f"{i}: {share(s['gnss_blind'])}" for i, s in rk.groupby(pb, observed=True)))
    print("  alternative definitions (blind share, main criterion):")
    er = eruption_table()
    vn = rk.set_index("vnum")
    pop_ok = rk.pop30.fillna(0) >= 1e5
    for lab, yr in (("eruption since 1500, pop30 >= 100k", 1500), ("eruption since 1950, pop30 >= 100k", 1950),
                    ("eruption since 2000, pop30 >= 100k", 2000)):
        act = rk.vnum.isin(er[er.StartDateYear >= yr].Volcano_Number)
        print(f"    {lab:40s} {share(rk[act & pop_ok]['gnss_blind'])}")
    print(f"    {'active and exposed, Holocene VEI max >= 4':40s} {share(ae[ae.vei_max >= 4]['gnss_blind'])}")
    print(f"    {'eruption since 1900, pop10 >= 10k':40s} {share(a[a.pop10.fillna(0) >= 1e4]['gnss_blind'])}")
    print(f"    {'eruption since 1900, pop100 >= 1M':40s} {share(a[a.pop100.fillna(0) >= 1e6]['gnss_blind'])}")

    print("== 5 Robustness table (active and exposed, n = 118, blind at 1e7 m^3)")
    sae = src[src.active & src.exposed].merge(rk[["vnum", "lat", "lon"]], on="vnum")
    rows = [("Main: network SNR 5, 30-day, Mogi 5 km below summit", ae, ae.net_d5, rk.net_d5),
            ("Network SNR 3", ae, ae.net_d5 * 3 / 5, rk.net_d5 * 3 / 5),
            ("Depth 2 km", ae, ae.net_d2.fillna(np.inf), rk.net_d2.fillna(np.inf)),
            ("Depth 10 km", ae, ae.net_d10.fillna(np.inf), rk.net_d10.fillna(np.inf)),
            ("Mogi 5 km off summit, median of 8 azimuths", sae, sae.dv_off_med, src.dv_off_med),
            ("Mogi 5 km off summit, worst azimuth", sae, sae.dv_off_max, src.dv_off_max),
            ("Horizontal point sill", sae, sae.dv_sill, src.dv_sill),
            ("First pass: 1 sigma daily noise, one station", ae, ae.dv_d5_k1, rk.dv_d5_k1),
            ("Without the 2 volcanic fields", ae[~ae.field], ae[~ae.field].net_d5, rk[~rk.field].net_d5)]
    for lab, d, v, allv in rows:
        b = d[v.to_numpy() > DV]
        k, n = len(b), len(d)
        lo, hi = wilson(k, n)
        print(f"  {lab:50s} {k:3d}/{n} ({100 * k / n:.0f} %, CI {100 * lo:.0f}-{100 * hi:.0f}); people "
              f"{bounds_str(b)}; all subaerial with <= 1e7 {(allv <= DV).sum()}")
    print(f"  sill column equals netsill_d5: {np.allclose(src.dv_sill.replace(np.inf, -1), rk.netsill_d5.fillna(np.inf).replace(np.inf, -1), rtol=1e-3)}")
    for lab, v in (("both-blind, S1 worst season (main)", ae.both_blind),
                   ("both-blind, S1 best season", ae.gnss_blind & ~ae.s1_sees_best),
                   ("both-blind, S1 1 cm worst season", ae.gnss_blind & ~(ae.insA_d5_t10_worst <= DV)),
                   ("both-blind, S1 single pixel worst season", ae.gnss_blind & ~(ae.ins_d5_t59_worst <= DV))):
        print(f"  {lab:50s} {share(v)}; people {bounds_str(ae[v])}")
    r = src.dv_mogi / src.dv_sill
    r = r[np.isfinite(r)]
    print(f"  Mogi/sill network dV ratio: median {r.median():.3f}, 5-95 % {r.quantile(.05):.3f}-{r.quantile(.95):.2f} "
          f"(n {len(r)})")
    up = ok["sig30_u_mm"].median() / 1e3 * np.sqrt(2)
    print(f"  one station directly above a 5 km source, median 30-day vertical noise x sqrt 2 ({up * 1e3:.2f} mm), "
          f"SNR 5 on the vertical alone: {5 * up / V.mogi_greens(0.0, 0.0, 5.0)[2]:.3g} m^3; peak uplift of 1e6 m^3 at "
          f"5 km {V.mogi_greens(0.0, 0.0, 5.0)[2] * 1e6 * 1e3:.2f} mm")

    print(f"== 6 National networks at the {len(nn)} GNSS-blind active and exposed volcanoes (national_networks.csv)")
    nn2 = nn.merge(rk[["vnum", "pop30", "lat", "lon"]], on="vnum", how="left")
    for c, s in nn2.groupby("category"):
        print(f"  {c:28s} {len(s):3d} volcanoes, people {bounds_str(s)}")
    print(f"  GVMID GNSS site within 30 km: {(nn.gvmid_gnss_n > 0).sum()}; Indonesian with no GVMID deformation "
          f"station: {((nn.country == 'Indonesia') & (nn.gvmid_gnss_n == 0) & (nn.gvmid_other_deformation_n == 0)).sum()}"
          f" of {(nn.country == 'Indonesia').sum()}; none_found in Indonesia: "
          f"{((nn.country == 'Indonesia') & (nn.category == 'none_found')).sum()}; rows resting on search summaries: "
          f"{nn.evidence.str.contains('search-summary').sum()}")
    newv = set(ae[ae.gnss_blind & ~ae.gnss_blind_1s].vnum)
    print(f"  added in the revision (GVMID only): {nn[nn.vnum.isin(newv)][['name', 'category']].values.tolist()}")

    print("== 7 Benchmark, Mount Shasta to Medicine Lake (first-pass criterion, 5 km, 1 sigma)")
    benchmark_check()

    print("== 8 Figures")
    cat = fig1(rk)
    print(f"  Fig. 1 categories (all subaerial): " + ", ".join(f"{CAT_STYLE[i][0]}: {(cat == i).sum()}" for i in range(5)))
    print(f"  Fig. 1 categories (active and exposed): " +
          ", ".join(f"{CAT_STYLE[i][0]}: {(cat[rk.ae] == i).sum()}" for i in range(5)))
    fig2(rk, ae)
    print(f"  Fig. 2 curve ends (finite values): GNSS main {np.isfinite(rk.net_d5).sum()}, GNSS 1 sigma "
          f"{np.isfinite(rk.dv_d5_k1).sum()}, S1 area best {np.isfinite(rk.insA_d5_t59_best).sum()}, worst "
          f"{np.isfinite(rk.insA_d5_t59_worst).sum()}; for the 118: {np.isfinite(ae.net_d5).sum()}, "
          f"{np.isfinite(ae.dv_d5_k1).sum()}, {np.isfinite(ae.insA_d5_t59_best).sum()}, "
          f"{np.isfinite(ae.insA_d5_t59_worst).sum()}")
    rows3, g = fig3(rk, blind, mh_main)
    print("  Fig. 3a (group, set, blind k/n, Wilson 95 % CI):")
    for gname, lab, k, n, lo, hi in rows3:
        print(f"    {gname:28s} {lab:20s} {k:3d}/{n:<4d} {100 * k / n:5.1f} % ({100 * lo:.1f}-{100 * hi:.1f})")
    print("  Fig. 3b rows (volcanoes, people lower-upper M):")
    for c, r in g.iterrows():
        print(f"    {c:28s} {int(r.n):3d}  {mpop(r.lo)}-{mpop(r.up)}")
    print("  written to", OUT)


if __name__ == "__main__":
    main()
