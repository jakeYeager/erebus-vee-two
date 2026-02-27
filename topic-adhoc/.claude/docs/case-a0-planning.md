# Topic Adhoc - Case A0 Planning

## Intent

Two data files of events populations from different USGS catalogs: ComCat and ISC-GEM. The ComCat catalog was used in previous research and has known data quality limitations. This comparison documents the structural and population differences between the two catalogs to serve as a reference when interpreting legacy ComCat results alongside future ISC-GEM-based analyses.

All future case studies in this project use the ISC-GEM catalog.

---

## Rationale

**Rationale for comparison**

Extensive use of the ComCat population in prior analysis exposes those results to potential data quality issues. A documented record of the differences between the two catalogs is necessary for interpreting legacy results and understanding how the catalog change may affect cross-study comparisons.

**Rationale for catalog change**

The rationale for migrating the data pipeline to ISC-GEM can be found here: https://raw.githubusercontent.com/jakeYeager/nornir-urd-two/refs/heads/main/review/isc-gem_catalog_migration.md

**Data files**

Both files were collected using the data pipeline application https://github.com/jakeYeager/nornir-urd-two

- ComCat: `data/global-sets/comcat_global_6-9_1949-2021.csv`
- ISC-GEM: `data/global-sets/iscgem_global_events.csv`

---

## Case A0 - Catalog Comparison Reference Report

### Purpose

Generate a single descriptive comparison report documenting the structural and population differences between the ComCat and ISC-GEM catalogs. This report serves as the reference document for interpreting legacy ComCat analysis results in the context of future ISC-GEM results — not as a validation or pass/fail exercise.

### Key Findings to Document

The following differences were identified during planning and should be described and quantified in the report:

**1. Hybrid catalog composition (ComCat)**

ComCat is not a single-source catalog. Its event IDs reveal it sources from multiple contributing catalogs:

| ID prefix | Count | % of population |
| --- | --- | --- |
| `us...` (USGS native) | 6,716 | 68.5% |
| `iscgem...` (ISC-GEM sourced) | 2,973 | 30.3% |
| Other (regional networks) | 113 | 1.2% |

Approximately 30% of ComCat records are sourced directly from ISC-GEM, meaning the two catalogs are not independent. ISC-GEM is both a standalone catalog and a contributing source within ComCat. The report should characterize the temporal distribution of the `iscgem`-prefixed records within ComCat — specifically whether they are concentrated in the pre-1976 historical period where USGS native coverage is sparse.

**2. Population counts**

| Catalog | Events | Year range | Mag range |
| --- | --- | --- | --- |
| ComCat | 9,802 | 1949–2021 | M 6.0–9.50 |
| ISC-GEM | 9,210 | 1949–2021 | M 6.0–9.55 |

ComCat contains 592 more events (+6.4%). The temporal and magnitude range is essentially identical across both populations.

**3. Magnitude precision and methodology**

| Catalog | 1-decimal mags | 2-decimal mags |
| --- | --- | --- |
| ComCat | 7,601 (77.5%) | 2,201 (22.5%) |
| ISC-GEM | 1,575 (17.1%) | 7,635 (82.9%) |

ComCat is predominantly 1-decimal precision, reflecting aggregation of legacy magnitude scales (Ms, mb) that were historically reported at lower resolution. ISC-GEM uses unified moment magnitude (Mw) converted to continuous precision throughout. This difference is consequential for magnitude-based analyses (b-value, magnitude stratification) particularly near the M 6.0 completeness threshold.

**4. Schema**

Both files share an identical 10-column schema: `usgs_id, usgs_mag, event_at, solaration_year, solar_secs, lunar_secs, midnight_secs, latitude, longitude, depth`. No structural incompatibility exists between the two populations.

### Report Outputs

- Population summary table for both catalogs
- ID prefix breakdown with temporal distribution of hybrid ComCat records
- Magnitude distribution comparison (overlaid histograms or binned table)
- Magnitude precision breakdown table
- Summary narrative suitable for cross-referencing in future case whitepapers

### Out of Scope

This case does not perform event-level matching, cross-catalog deduplication, or statistical significance testing. It is a descriptive reference document only.
