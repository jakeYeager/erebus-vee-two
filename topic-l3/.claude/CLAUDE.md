# Topic: Analysis of Global Dataset Description (legacy "Approach Three")

## Intent

This approach in the project reviews the statistical testing suites from Approaches One and Two, and leverages a fresh global dataset with post-1950 artifacts removed, as well as ephemeris-verified calculations for astronomical metrics. This approach removes the blind study isolation constraints applied to Approach Two and uses non-anonymized column labels, allowing full interpretation of results against the known seismic and gravitational forcing context.

## Data Description
- File location: `data/global-sets/global_events.csv`
- Total records: 9,802
- CSV columns headers: `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `longitude`
- Source: USGS Earthquake Catalog; location non-specific (global) post-1950 events M>=6.0
- Astronomical metrics recomputed with verified ephemeris calculations
  - Python library used: Skyfield 1.54 -- uses JPL DE421 ephemeris
  - Collected with code from https://github.com/jakeYeager/nornir-urd-two

## Analysis Framework: Cases 0, 1/1B, 2/2B, 3A/3B, 4A/4B

- Case 0: Population Description
- Case 1: Distribution Uniformity Testing
- Case 1B: Lunar Bin 16 Deficit Investigation
- Case 2,2B: Inter-Event Interval Analysis
- Case 3A/B: Clustering Patterns
- Case 4A/B: Energy-Weighted Validation

Only read this file if you need to know the full summary of this topic: `topic-global-description/final-summary.md`