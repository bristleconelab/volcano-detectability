# Internal reviews and responses (2026-11-17)

Two independent reviews of the first draft (a volcano-geodesy referee and a risk / statistics referee; full reports
kept privately) both recommended major revision. The arithmetic reproduced; the conclusions and criteria did not hold.

| Issue raised | Response |
|---|---|
| 1 sigma at one station has no false-alarm control; displacement noise is sqrt(2) sigma | Main criterion is now a network matched filter at SNR 5 on differences of 30-day means (sqrt(2) included). The USGS-style criterion is kept only for the benchmark. |
| 30-day noise inflated by seasonal motion | Velocity and annual / semi-annual terms are fitted to block-mean differences before the robust scatter (unit-tested). Median 30-day noise fell from 0.96 / 0.90 / 2.50 mm to 0.56 / 0.65 / 1.64 mm (E / N / U). |
| Sentinel-1 "detection" from one coherent pixel | A season counts only with >= 10 % coherent pixels within 5 km. Sentinel-1 now sees 1e7 m3 at 602 (best) / 276 (worst season) volcanoes instead of 849 / 470. |
| "Threat does not bring coverage" untested; Simpson's paradox by country | Fisher and country-stratified Mantel-Haenszel tests added. Within countries, populous active volcanoes are better covered (OR 0.05); between countries, data policy decides (Japan vs Indonesia). The claim was replaced. |
| Population double counting | Reported as a range: lower bound from volcanoes more than 60 km apart, upper bound from the sum. |
| Title contradicted by the Sentinel-1 layer | Both-technique blindness reported separately (53 of 118 with the area criterion, worst season). |
| Volcanic fields as point sources | Flagged (`field`); results given with and without them. |
| Sill and offset sources | Kept as sensitivity; reported on the main criterion (sill 93 of 118). |
| Licensing of derived data | Computed values CC BY 4.0; copied GVP / GVMID columns under their terms. |
| Missing references (NVEWS: Ewert 2005/2018, Moran 2008; Chaussard & Amelung; GNSS noise literature; compressibility) | To be added in the manuscript revision. |
| Open: per-volcano validation of Sentinel-1 detections, injection tests on real GNSS series, a union population grid | Not done yet. |
