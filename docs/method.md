# Method and limitations

## Question

For each Smithsonian GVP Holocene volcano: what is the smallest volume change of a magma body below the summit that
open geodetic data would notice? Small values mean well watched; large values or none mean blind spots.

## Source model

A Mogi point source (Segall 2010, eq. 7.14) in an elastic half-space, Poisson's ratio 0.25, below the GVP summit
coordinate at 2, 5 or 10 km depth. Surface displacement is linear in the volume change dV, so every threshold below
has a closed form. Topography, layering, other source shapes (sills, dikes, prolate bodies) and source offsets from the
summit are not modelled; for a sill or an elongated source the visible signal differs by a factor of order one.

## GNSS layer

- Stations: Nevada Geodetic Laboratory (NGL) IGS20 daily solutions; active = solutions after 2026-08-01; within 100 km.
- Noise per station and component (east, north, up), from the newest 730 days:
  - daily white noise = 1.4826 MAD of consecutive-day differences / sqrt(2) (robust to offsets and outliers);
  - 30-day noise = same estimator on consecutive non-overlapping 30-day means (includes coloured noise).
- Detection: predicted displacement > K sigma on at least one component of at least one station, after USGS SIR
  2024-5062-D and SIR 2026-5117-A. dV_min = min over stations and components of K sigma / |G|.
- Main numbers use daily noise and K = 1, as USGS; K = 3 and 30-day noise with K = 1, 2 are given for sensitivity.

A single station passing 1 sigma would not trigger an alarm in practice; the criterion is a lower bound on what is
needed, used because it matches the USGS maps and is easy to reproduce.

## Main GNSS criterion (after review)

- Noise per station and component from the newest 730 days: the 30-day noise is the robust scatter of differences of
  consecutive non-overlapping 30-day means after a least-squares fit (with 3-MAD rejection) of a constant velocity and
  annual and semi-annual terms to those differences, so seasonal motion is not counted as noise.
- A displacement is the difference of two 30-day means, so its noise is sqrt(2) x the 30-day sigma.
- Detection: a matched filter for the known source over all stations within 100 km and all components reaches
  signal-to-noise 5: dV_min = 5 / sqrt(sum (G / sigma_d)^2). SNR 5 corresponds to a one-sided Gaussian false-alarm
  probability of about 3e-7 per test, and to roughly one false alarm in several thousand station-years of monthly tests.
- Limits: the noise model is white at the 30-day scale after the seasonal fit (flicker noise at longer periods would
  make slow, multi-year accumulation harder to see); the template assumes the source location and depth are known.

## Sentinel-1 layer

- Coherence: global seasonal 12-day VV coherence, Dec 2019 to Nov 2020 (Kellndorfer et al. 2022), 3 arcsec; a pixel is
  usable when its 5 x 5 mean coherence is at least 0.35. Seasons are kept separately (snow, vegetation).
- The best usable pixel within 20 km of the summit sees the fraction f = (1 + (r/d)^2)^-1.5 of the peak uplift.
- Main criterion (after review): a season counts only when at least 10 % of the pixels within 5 km of the summit are
  coherent, so that a deformation pattern, not one pixel, can be seen. The single-pixel version is kept for comparison.
- Revisit: Sentinel-1C (launched Dec 2024) restores 6-day repeat over much of the globe; the 2019-2020 coherence
  product (with Sentinel-1A and 1B) is representative of 12-day pairs.
- Line of sight = cos(39 deg) x vertical (mid-swath incidence; horizontal motion ignored).
- Threshold on cumulative LOS displacement: 5.9 cm (global threshold of the automated COMET/LiCSAR system, Biggs et al.
  2022) or 1 cm (careful time-series analysis; optimistic). Published automated thresholds range from about 3.9 cm
  (signal only) to 6.3 cm (with atmosphere; Anantrasirichai et al. 2019), and Biggs et al. (2026) call ~5 cm a
  "mid-range" value: the true threshold differs by volcano. Volumes scale linearly with the threshold.
- This layer is an a-priori estimate from coherence and geometry. Per-volcano measured noise from global time-series
  efforts (COMET/LiCSAR; NASA MEaSUREs) can replace the constant threshold once released.
- Not modelled: atmospheric noise per site (stratified delay at tall cones can mimic a summit-centred signal), the need
  for a stable reference area, layover and shadow, and coherence change since 2020. The Sentinel-1 numbers are
  therefore optimistic.

## Threat and exposure

GVP population within 30 and 100 km (GVP E3 layer) and the GVP eruption list (confirmed eruptions, start year, VEI).
"Active" = confirmed eruption that started in 1900 or later; "exposed" = at least 100,000 people within 30 km.

## What "blind" means

Only open data are counted. Many observatories (for example in Indonesia, the Philippines and Central America) run GNSS
networks whose data are not in NGL. The atlas shows where independent, open checking of magma accumulation is not
possible, not where volcanoes are unmonitored. Retrospective sensitivity only; no eruption forecasts.

## Benchmark

USGS SIR 2024-5062-D (Fig. D4) reports that more than 1e6 m3 could intrude undetected between Mount Shasta and Medicine
Lake at 5 km depth. With 53 open NGL stations this repository gives 5.1e5 m3 at Shasta, 3.0e5 m3 at Medicine Lake and
2.3e6 m3 at the midpoint (`scripts/benchmark_usgs.py`). USGS grid values are not published, so the comparison is of
pattern and order of magnitude, not cell by cell.
