# Declustered Catalog Analysis (legacy Approach Four)

## Intent

This topic tests whether the solar_secs clustering signal identified in previous topics persists in background seismicity after removing aftershock sequences using the Gardner-Knopoff (1974) declustering algorithm. If the signal survives declustering, it strengthens the case for a global gravitational forcing mechanism rather than an artifact of aftershock clustering.

## Data Description

- Primary dataset: `data/declustering/mainshocks.csv` (declustered population)
- Secondary reference dataset: `data/declustering/aftershocks.csv` (G-K designated aftershock events)
- Source: USGS Earthquake Catalog; global post-1950 events M>=6.0, declustered via Gardner-Knopoff (1974)
- CSV column headers: `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `longitude`, `latitude`, `depth`
- Astronomical metrics recomputed with verified ephemeris calculations
  - Python library used: Skyfield 1.54 -- uses JPL DE421 ephemeris
  - Collected with code from https://github.com/jakeYeager/nornir-urd-two


## Analysis Framework: Cases D0, D1, D2, D3A/D3B, D4A/D4B

- Case D0: Declustered Population Description
- Case D1: Distribution Uniformity Testing
- Case D2: Inter-Event Interval Analysis
- Case D3A/B: Clustering Patterns
- Case D4A/B: Energy-Weighted Validation

Only read this file if you need to know the final summary of this topic: `topic-l4/final-summary.md`
