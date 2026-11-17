# Volcano geodetic detectability (open data)

How much magma can collect beneath each of the world's Holocene volcanoes before open GNSS data would notice?
For every Smithsonian GVP Holocene volcano, this repository computes the smallest Mogi volume change at 2, 5 and 10 km
depth that would exceed the daily white noise of at least one open Nevada Geodetic Laboratory (NGL) GNSS station on
at least one component. The criterion follows USGS SIR 2024-5062-D and SIR 2026-5117-A, applied globally.

> Work in progress (first pass, 2026-11-17). This measures what **open** data can see, not what volcano observatories
> see: many national networks are not in NGL.

## First numbers (5 km depth, 1 sigma, daily noise)

- 1,100 GVP Holocene volcanoes at or above sea level; 2,477 active NGL stations within 100 km of one.
- 484 have no active open station within 100 km.
- A 1e6 m3 change is detectable at 173; a 1e7 m3 change at 376.
- Of 383 that erupted since 1900, 256 cannot detect 1e7 m3 and 171 have no station within 100 km.
- Benchmark: the USGS Mount Shasta to Medicine Lake map is reproduced in pattern (2-5e5 m3 at the edifices, about
  2e6 m3 between them; USGS reports more than 1e6 m3 could go undetected there).

![Benchmark](docs/img/benchmark_shasta_medicine.png)

## Run

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python -I -m unittest discover -s tests -v
curl -o data/raw/gvp_holocene.json "https://webservices.volcano.si.edu/geoserver/GVP-VOTW/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=GVP-VOTW:Smithsonian_VOTW_Holocene_Volcanoes&outputFormat=application/json"
curl -o data/raw/DataHoldings.txt https://geodesy.unr.edu/NGLStationPages/DataHoldings.txt
.venv/bin/python scripts/station_noise.py      # ~250 kB per station via HTTP Range, resumable
.venv/bin/python scripts/detectability.py
.venv/bin/python scripts/benchmark_usgs.py
```

Prior work and what is new: [docs/prior-work.md](docs/prior-work.md).

## Data and licences

GVP Volcanoes of the World (Smithsonian, non-commercial terms), NGL GNSS (Blewitt, Hammond & Kreemer 2018, Eos 99).
Derived tables: CC BY-NC 4.0 because of the GVP terms. Code: MIT.
