# Prior work

Checked 2026-11-17 (searches: OpenAlex 2023-2026, EGU 2024-2026 via OpenAlex/web, AGU24-25 meeting search, journal sites).
A gap below means "not found in these searches", not proof of absence.

| Work | What it does | What it does not do |
|---|---|---|
| Montgomery-Brown et al. 2024, USGS SIR 2024-5062-D | Mogi grid at 5 km, inflate until displacement exceeds GNSS white noise at one station; one map (Mount Shasta to Medicine Lake), >1e6 m3 could intrude undetected between them | One area; no code or data release; no InSAR |
| Montgomery-Brown 2026, USGS SIR 2026-5117-A (doi:10.3133/sir20265117A) | Same method as a reusable methodology: 1 km grid, 2/5/8 km depths, per-station white noise on 3 components, seasonal loss; Lassen and Shasta examples; aims at all US volcanoes | US only; InSAR qualitative; no ranking |
| Widiwijayanti et al. 2024, GVMID (Frontiers) | Deformation instrument inventory, 571 volcanoes, 8,924 stations | Counts only; no sensitivity model |
| Brown et al. 2015 (GVM) | Monitoring levels 0-5 by instrument presence | Not quantitative; not updated |
| Biggs et al. 2022 (Bull. Volcanol.) | One global Sentinel-1 detection threshold (5.9 cm; 1.2 cm/yr) | Not per volcano; no volume |
| Espinosa-Ortega & Taisne 2024 (JVGR) | Seismic network detection capability at >100 volcanoes (GVMID stations) | Seismic, not geodetic |
| Lee 2023, NTU undergraduate project (supervisor Taisne) | Continuous-GPS network sensitivity at 184 volcanoes | Unpublished; full text not read (CAPTCHA) |
| Poland et al. 2026 (Bull. Volcanol. 88:74), GVEWERS | Satellite monitoring programme history | No per-volcano detectability (full text searched) |
| Pritchard et al. 2026 (Bull. Volcanol. 88:75) | Where daily high-resolution satellite imaging is needed, by threat and volcano type | Tasking priority, not geodetic detectability |
| Meredith et al. 2025 (NHESS 25:2731), 2026 (npj Nat. Hazards) | Urban exposure per Holocene volcano | No monitoring |

This repository differs by: global scope (all GVP Holocene volcanoes), open NGL stations with per-station noise estimated
from the data, a Sentinel-1 layer (planned), and a threat x geodetic-blindness ranking (planned).
