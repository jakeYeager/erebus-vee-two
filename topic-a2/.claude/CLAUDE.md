# Topic: Seasonal Solar Signal Investigation

**STATUS: Active**

## Intent

The bimodal equinox pattern is not straightforwardly explained by the dominant mechanism in the literature (hydrological/snow loading). Hydrological forcing peaks at summer or winter solstice depending on hemisphere and regional climate — not at the equinoxes. The equinox signal is more consistent with a direct solar-geometric forcing (e.g., solar declination rate of change, Earth-Sun distance, or tidal stress geometry) or with a globally symmetric effect that cancels hemisphere-specific hydrological loading and leaves a residual equinox peak. This distinction is the core scientific question.

## Data Sources

All cases use the ISC-GEM catalog (9,210 events, M ≥ 6.0, 1950–2021). Prior results (Case 3A context) used the ComCat catalog (9,802 events). See Cases A0 and A0b in topic-adhoc for catalog differences.

**Base catalog** (`data/iscgem/`):
- `iscgem_global_6-9_1950-2021.csv` — full ISC-GEM catalog with ephemeris columns

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
| 1        | A4   | Complete | Declustering Sensitivity Analysis                             |
| 2        | B6   | Complete | Rolling Window Stationarity Test                              |
| 3        | A1   | Complete | Schuster Spectrum and MFPA Periodicity Analysis               |
| 4        | B1   | Complete | Hemisphere Stratification — Phase Symmetry Test               |
| 5        | A3   | Pending | Magnitude Stratification of the Solar Signal                   |
| 6        | B2   | Pending | Ocean vs. Continent Location — Hydrological Loading Discrimination |
| 7        | B4   | Pending | Depth Stratification — Surface Loading Penetration Test        |
| 8        | A2   | Pending | b-Value Seasonal Variation                                     |
| 9        | B3   | Pending | Tectonic Regime Stratification                                 |
| 10       | B5   | Pending | Solar Declination Rate-of-Change vs. Position Test             |

## Analysis Framework

Only read this file if you need the full description of "Case A4: Declustering Sensitivity Analysis": `topic-a2/.claude/spec/case-a4-spec.md`

Only read this file if you need the full description of "Case B6: Rolling Window Stationarity Test": `topic-a2/.claude/spec/case-b6-spec.md`

Only read this file if you need the full description of "Case A1: Schuster Spectrum and MFPA Periodicity Analysis": `topic-a2/.claude/spec/case-a1-spec.md`

Only read this file if you need the full description of "Case B1: Hemisphere Stratification — Phase Symmetry Test": `topic-a2/.claude/spec/case-b1-spec.md`

Only read this file if you need the full description of "Case A3: Magnitude Stratification of the Solar Signal": `topic-a2/.claude/spec/case-a3-spec.md`
**Pre-run:** *User Confirmed* - apply Julian constant (31,557,600 s) should be used uniformly for year phase normalization.

Only read this file if you need the full description of "Case B2: Ocean vs. Continent Location — Hydrological Loading Discrimination": `topic-a2/.claude/spec/case-b2-spec.md`
**Pre-run:** (1) *User Confirmed* - apply Julian constant (31,557,600 s) should be used uniformly for year phase normalization. (2) Verify cartopy or basemap availability in the project Python environment for global map figures; if neither is available, use lat/lon scatter with Natural Earth coastline path via geopandas.

Only read this file if you need the full description of "Case B4: Depth Stratification — Surface Loading Penetration Test": `topic-a2/.claude/spec/case-b4-spec.md`
**Pre-run:** *User Confirmed* - apply Julian constant (31,557,600 s) should be used uniformly for year phase normalization.

Only read this file if you need the full description of "Case A2: b-Value Seasonal Variation": `topic-a2/.claude/spec/case-a2-spec.md`
**Pre-run:** *User Confirmed* - apply Julian constant (31,557,600 s) should be used uniformly for year phase normalization.

Only read this file if you need the full description of "Case B3: Tectonic Regime Stratification": `topic-a2/.claude/spec/case-b3-spec.md`
**Pre-run:** (1) *User Confirmed* - apply Julian constant (31,557,600 s) should be used uniformly for year phase normalization. (2) Verify that `focal_join_global.csv` includes `solar_secs` as a base-catalog-prefix column, or whether the raw catalog must be joined separately to supply `solar_secs` for matched events. *User Confirmed* - join exists already from data pipeline treatment.

Only read this file if you need the full description of "Case B5: Solar Declination Rate-of-Change vs. Position Test": `topic-a2/.claude/spec/case-b5-spec.md`
**Pre-run:** Verify actual min/max ranges of `solar_declination`, `declination_rate`, and `earth_sun_distance` in `solar_geometry_global.csv` against the REQ-3 spec ranges (−23.5 to +23.5°, ~−0.40 to +0.40 deg/day, ~0.983 to ~1.017 AU); update bin-computation formula if actual ranges differ.
