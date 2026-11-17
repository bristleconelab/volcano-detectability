# Corrections

- 2026-11-17: output tables were written with `float_format="%.4g"`, which rounded volcano and station coordinates to
  four significant digits (up to 0.05 deg, about 5 km at longitudes above 100 deg). The detectability values were
  computed from the exact GVP and NGL coordinates and did not change; the CSV coordinates and the web atlas positions
  were wrong. All tables are now written with 7 significant digits; `scripts/source_sensitivity.py` checks that a
  recomputation from the tables reproduces `dv_d5_k1`, and the maximum coordinate difference from GVP is 0.
