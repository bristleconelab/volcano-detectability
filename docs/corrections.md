# Corrections

- 2026-11-17: output tables were written with `float_format="%.4g"`, which rounded volcano and station coordinates to
  four significant digits (up to 0.05 deg, about 5 km at longitudes above 100 deg). The detectability values were
  computed from the exact GVP and NGL coordinates and did not change; the CSV coordinates and the web atlas positions
  were wrong. All tables are now written with 7 significant digits; `scripts/source_sensitivity.py` checks that a
  recomputation from the tables reproduces `dv_d5_k1`, and the maximum coordinate difference from GVP is 0.
- 2026-11-17: README used GVP's last-eruption-year field for "erupted since 1900" (383 volcanoes, 256 blind, 171 without
  a station), while the ranking uses confirmed eruptions with start year >= 1900 (380, 253, 168). README now uses the
  ranking definition. README also said Sentinel-1 sees nothing below "about 6e6 m3"; the floor with the line-of-sight
  factor is 7.95e6 m3 (6.18e6 is the value without it). Both found by an independent re-computation of every number.
