# Corrections

- 2026-11-17: output tables were written with `float_format="%.4g"`, which rounded volcano and station coordinates to
  four significant digits (up to 0.05 deg, about 5 km at longitudes above 100 deg). The detectability values were
  computed from the exact GVP and NGL coordinates and did not change; the CSV coordinates and the web atlas positions
  were wrong. All tables are now written with 7 significant digits; `scripts/source_sensitivity.py` checks that a
  recomputation from the tables reproduces `net_d5`, and the maximum coordinate difference from GVP is 0.
- 2026-11-17: README used GVP's last-eruption-year field for "erupted since 1900" (383 volcanoes, 256 blind, 171 without
  a station), while the ranking uses confirmed eruptions with start year >= 1900 (380, 253, 168). README now uses the
  ranking definition. README also said Sentinel-1 sees nothing below "about 6e6 m3"; the floor with the line-of-sight
  factor is 7.95e6 m3 (6.18e6 is the value without it). Both found by an independent re-computation of every number.
- 2026-11-17: the station prefilter used a fixed +-2.5 deg longitude window, narrower than 100 km poleward of about 69 deg.
  Three Antarctic NGL stations (BACK, MIN0, TOMO) were never fetched, and Hudson Mountains (Antarctica) was counted as
  having no station within 100 km (it has BACK at 92.6 km; dV at 5 km ~2e9 m3). Fixed with a latitude-dependent window
  (`vdlib.lon_window_deg`, unit-tested). Changes: stations 2,477 -> 2,480; volcanoes without a station 484 -> 483;
  fully blind (GNSS and Sentinel-1) 69 -> 68. No change to 173, 376, 118, 70, 72.3 M or the sill case. Found by the
  independent verification (notes/verification-2026-11-17.md, PASS on all other items).
- 2026-11-17: `docs/data/national_networks.csv` listed GVMID station codes in its evidence column. The WOVOdat data
  policy does not allow redistribution of GVMID data, so the table now gives counts and classes only (the category of
  every volcano is unchanged). Earlier commits in the history still contain the codes (no history rewrite).
- 2026-11-20: full codebase, scientific, and web-atlas audit fixes:
  - `scripts/vdlib.py`: `parse_tenv3_tail` now reads `_ant(m)` (`t[13]`), and `block_noise_sigma` rejects intra-block
    antenna toggles and >5-MAD daily spikes, enforces an OLS inlier-count guard (`>= max(p + 4, ceil(0.6 * n))`),
    applies a `sqrt(k / (k - p))` degrees-of-freedom correction, and enforces the physical white-noise lower bound
    `sig_daily / sqrt(30)`. Eliminated all 143 sub-white-noise component collapses (e.g. `CANI` 0.008 -> 0.69 mm) and
    intra-block `P656` corruption (40.8 -> 1.90 mm). Active & exposed GNSS-blind set updated from 75 to 76/118
    (Nyamulagira became detectable at 7.1e6 m3; Masaya and Cumbal moved just above 1e7 m3).
  - `scripts/detectability.py`: added explicit `n_100km`, `nearest_site`/`nearest_km`, `n30_*km`, and `net_best_site_d5`,
    separating volcanoes with no active station ≤100 km (483; 44 of 76 AE blind) from those lacking an ≥8-pair 30-day
    estimate (501; 48 of 76 AE blind).
  - `scripts/insar_coherence.py`: wrapped antimeridian tile names (`-180..+180`) and stitched windows across ±180°, plus
    automatic `HH` polarization fallback where `VV` tiles are absent (covering 16 of 18 Antarctic volcanoes, 10 usable).
  - `scripts/national_networks/`: prevented open-archive inactive sites (El Salvador) from falling into
    `observatory_cgnss_not_open`, and added `open_gnss_insufficient` for volcanoes with active open NGL stations ≤30 km
    (Popocatepetl, Masaya, Cumbal).
  - `docs/atlas/`: replaced flat Natural Earth polygon fills with real zoomable Esri Satellite, Ocean Bathymetry, and
    Terrain basemaps with city/place overlays; simplified to a minimalist search-first HUD with 1-click preset pills,
    a 30 km exposure circle and dashed GNSS link on selection, a plain-language volcano inspector card with collapsible
    2/5/10 km math breakdown, lazy-loaded GNSS station noise popups, XSS escaping, responsive mobile layout, and
    4-digit/3-dp GeoJSON precision (matching exact `736` headline detectability and trimming initial GeoJSON payload
    from ~2.98 MB to ~649 KB).
