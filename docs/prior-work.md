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
| Biggs et al. 2026 (RSE 339:115377), Sentinel-1 decade review | States that detection thresholds differ by volcano and that there is no simple mapping from threshold to rate; uses ~5 cm as "mid-range" | No per-volcano threshold table or map (full text checked) |
| Ebmeier et al. 2013 (JGR, doi:10.1002/jgrb.50195) | Per-volcano minimum detectable InSAR rates at 20 Central American volcanoes (L-band, mean 2.4 cm/yr) | Regional; rates, not volumes |
| Anantrasirichai et al. 2019 (GRL) | Detection thresholds for automated InSAR: about 3.9 cm (signal only) to 6.3 cm (with atmosphere) | Global numbers, not per volcano |
| Espín Bedón et al. 2026 (SSRN, doi:10.2139/ssrn.7052099) | Coherence-guided Sentinel-1 time series at ~230 volcanic systems, 2014-2026 (abstract via Crossref) | No thresholds or volumes in the abstract; full text not read (403) |
| Hooper et al. EGU26-18724; Shkreli, Zebker et al. AGU25 G13C-0277 (NASA MEaSUREs) | Global InSAR processing at ~1,300 / ~1,600 volcanoes; MEaSUREs ~2 cm RMS vs GPS at 10 test volcanoes | No per-volcano detection limits announced; their measured noise could replace our a-priori constants once released |
| Meredith et al. 2025 (NHESS 25:2731), 2026 (npj Nat. Hazards) | Urban exposure per Holocene volcano | No monitoring |

This repository differs by: global scope (all GVP Holocene volcanoes), open NGL stations with per-station noise estimated
from the data, an a-priori Sentinel-1 layer expressed as volume on the same scale as GNSS, and a threat x geodetic-blindness ranking.
