# Topic: Seasonal Solar Signal Investigation

**Status: Initial Case Planning**

## Intent

The bimodal equinox pattern is not straightforwardly explained by the dominant mechanism in the literature (hydrological/snow loading). Hydrological forcing peaks at summer or winter solstice depending on hemisphere and regional climate — not at the equinoxes. The equinox signal is more consistent with a direct solar-geometric forcing (e.g., solar declination rate of change, Earth-Sun distance, or tidal stress geometry) or with a globally symmetric effect that cancels hemisphere-specific hydrological loading and leaves a residual equinox peak. This distinction is the core scientific question.

## Data Sources

All cases use the ISC-GEM catalog (9,210 events, M ≥ 6.0, 1950–2021). Prior results (Case 3A context) used the ComCat catalog (9,802 events). See Cases A0 and A0b in topic-adhoc for catalog differences.

**Base catalog** (`data/global-sets/`):
- `iscgem_global_events.csv` — full ISC-GEM catalog with ephemeris columns

**Declustered catalogs** (`data/iscgem/declustering-algorithm/`):
- `mainshocks_G-K_global.csv` / `aftershocks_G-K_global.csv` — Gardner-Knopoff method (5,883 / 3,327)
- `mainshocks_reas_global.csv` / `aftershocks_reas_global.csv` — Reasenberg method (8,265 / 945)
- `mainshocks_a1b_global.csv` / `aftershocks_a1b_global.csv` — A1b-informed custom window (7,137 / 2,073)

**Enrichment files** (`data/iscgem/`):
- `solar_geometry_global.csv` — solar declination, declination rate, and Earth-Sun distance per event
- `focal_join_global.csv` — GCMT focal mechanism join; mechanism, rake, strike, dip (52.9% match rate)
- `ocean_class_gshhg_global.csv` — ocean/continent classification via GSHHG coastline (primary)
- `ocean_class_ne_global.csv` — ocean/continent classification via Natural Earth coastline (secondary)
- `ocean_class_pb2002_global.csv` — ocean/continent classification via PB2002 proximity (coarse proxy)

## Case Table

| Priority | Case | Status  | Title                                                          |
| -------- | ---- | ------- | -------------------------------------------------------------- |
| 1        | A4   | Pending | Declustering Sensitivity Analysis                              |
| 2        | B6   | Pending | Rolling Window Stationarity Test                               |
| 3        | A1   | Pending | Schuster Spectrum and MFPA Periodicity Analysis                |
| 4        | B1   | Pending | Hemisphere Stratification — Phase Symmetry Test                |
| 5        | A3   | Pending | Magnitude Stratification of the Solar Signal                   |
| 6        | B2   | Pending | Ocean vs. Continent Location — Hydrological Loading Discrimination |
| 7        | B4   | Pending | Depth Stratification — Surface Loading Penetration Test        |
| 8        | A2   | Pending | b-Value Seasonal Variation                                     |
| 9        | B3   | Pending | Tectonic Regime Stratification                                 |
| 10       | B5   | Pending | Solar Declination Rate-of-Change vs. Position Test             |

**Eager Load Planning:** [warning user if topic status changes] Remove eager-loading when planning is complete: @topic-a2/.claude/docs/planning-initial.md

<!-- --resume "seismic-solar-phase-analysis"  -->