# Case A0 Spec: Catalog Comparison Reference Report

## Status

Ready for implementation.

---

## Overview

Produce a descriptive reference report documenting the structural and population differences between the ComCat and ISC-GEM earthquake catalogs. This is a one-time reference document, not a recurring analysis. No event-level matching, cross-catalog deduplication, or statistical significance testing is performed.

---

## Data Inputs

| File | Path |
| --- | --- |
| ComCat | `data/global-sets/comcat_global_6-9_1949-2021.csv` |
| ISC-GEM | `data/global-sets/iscgem_global_events.csv` |

**Schema (both files):** `usgs_id, usgs_mag, event_at, solaration_year, solar_secs, lunar_secs, midnight_secs, latitude, longitude, depth`

Paths are relative to the project root (`/Users/jake/Projects/Code/prod/erebus-vee-two/`).

---

## Deliverables

All outputs go under `topic-adhoc/`.

| Type | Path |
| --- | --- |
| Analysis script | `src/case-a0-analysis.py` |
| Visualization script | `src/visualization-case-a0.py` |
| Test suite | `tests/test-case-a0.py` |
| Results JSON | `output/case-a0-results.json` |
| Whitepaper | `output/case-a0-whitepaper.md` |
| Image: mag distributions | `output/case-a0-magnitude-distributions.png` |
| Image: ComCat ID prefix temporal | `output/case-a0-comcat-prefix-temporal.png` |

---

## Analysis Script: `case-a0-analysis.py`

Loads both CSVs, computes all statistics, and writes `case-a0-results.json`. No visualizations in this script.

### Computations

**1. Population summary (both catalogs)**
- Event count (rows minus header)
- Year range: min and max of `solaration_year`
- Magnitude range: min and max of `usgs_mag`

**2. ComCat ID prefix breakdown**

Classify each `usgs_id` into one of three groups:
- `us_native`: starts with `"us"` but NOT `"uscgem"` or `"iscgem"` — i.e., matches `^us[^c]` or `^us$` — specifically: starts with `"us"` and does NOT start with `"iscgem"`
  - Clarification: prefix logic should be: if `usgs_id.startswith("iscgem")` → `iscgem`; elif `usgs_id.startswith("us")` → `us_native`; else → `other`
- `iscgem`: starts with `"iscgem"`
- `other`: anything else

For each group: count and percentage of total ComCat population.

**3. Temporal distribution of `iscgem`-prefixed ComCat records**

For ComCat rows where `usgs_id` starts with `"iscgem"`, compute counts by decade using `solaration_year`:
- Decades: 1940s, 1950s, 1960s, 1970s, 1980s, 1990s, 2000s, 2010s, 2020s
- Record the count per decade and the percentage of all `iscgem`-prefix ComCat records

**4. Magnitude precision breakdown (both catalogs)**

For each catalog, classify each `usgs_mag` value:
- `one_decimal`: magnitude value has exactly 1 decimal place (e.g., 6.2, 7.0)
  - Detection: `round(mag, 1) == mag and round(mag, 2) != mag` — more robustly: check if `(mag * 10) % 1 == 0`
  - Practical approach: convert to string, check decimal places after the decimal point — but floating point can be unreliable. Use: `len(str(round(mag, 2)).rstrip('0').split('.')[-1]) == 1` or equivalently check `(round(mag * 100) % 10) == 0` and `(round(mag * 10) % 10) != 0 or mag == round(mag, 1)`
  - Simplest reliable method: multiply by 100, round to integer, check if divisible by 10 → 1 decimal; if not → 2 decimal
- `two_decimal`: all others (2 decimal places)

For each catalog: count and percentage for each precision class.

**5. Magnitude distribution bins**

For binned comparison in visualization: compute event counts per 0.1-magnitude bin from 6.0 to 9.6 for both catalogs.

### Output: `case-a0-results.json`

```json
{
  "generated": "<ISO timestamp>",
  "catalogs": {
    "comcat": {
      "file": "data/global-sets/comcat_global_6-9_1949-2021.csv",
      "event_count": 9802,
      "year_range": {"min": 1949, "max": 2021},
      "mag_range": {"min": 6.0, "max": 9.5}
    },
    "iscgem": {
      "file": "data/global-sets/iscgem_global_events.csv",
      "event_count": 9210,
      "year_range": {"min": 1949, "max": 2021},
      "mag_range": {"min": 6.0, "max": 9.55}
    }
  },
  "comcat_prefix_breakdown": {
    "us_native": {"count": 6716, "pct": 68.5},
    "iscgem": {"count": 2973, "pct": 30.3},
    "other": {"count": 113, "pct": 1.2}
  },
  "comcat_iscgem_prefix_temporal": {
    "1940s": {"count": 0, "pct": 0.0},
    "1950s": {"count": ..., "pct": ...},
    ...
  },
  "magnitude_precision": {
    "comcat": {
      "one_decimal": {"count": 7601, "pct": 77.5},
      "two_decimal": {"count": 2201, "pct": 22.5}
    },
    "iscgem": {
      "one_decimal": {"count": 1575, "pct": 17.1},
      "two_decimal": {"count": 7635, "pct": 82.9}
    }
  },
  "magnitude_bins": {
    "bin_width": 0.1,
    "bins": [6.0, 6.1, ...],
    "comcat_counts": [...],
    "iscgem_counts": [...]
  }
}
```

Note: The values above (counts, percentages) are from the planning document and serve as expected reference values. The script computes them from the data; the test suite validates they match.

---

## Visualization Script: `visualization-case-a0.py`

Reads `output/case-a0-results.json` and produces two PNG images.

### Image 1: `case-a0-magnitude-distributions.png`

Overlaid magnitude distribution comparison.

- Type: overlapping step histograms (or bar charts with side-by-side bins)
- X-axis: Magnitude (6.0–9.6, bins of 0.1)
- Y-axis: Event count
- ComCat bars: steelblue, semi-transparent (alpha 0.6)
- ISC-GEM bars: a contrasting color (e.g., coral/salmon), semi-transparent (alpha 0.6)
- Legend identifying each catalog with event count in label (e.g., "ComCat (n=9,802)")
- Title: "Magnitude Distribution: ComCat vs ISC-GEM"
- Annotation: note the magnitude precision difference (e.g., "ComCat: 77.5% 1-decimal | ISC-GEM: 82.9% 2-decimal")
- 300 DPI, saved to `output/`

### Image 2: `case-a0-comcat-prefix-temporal.png`

Stacked or grouped bar chart showing ComCat event counts by decade, colored by ID prefix group.

- X-axis: Decade (1950s through 2020s)
- Y-axis: Event count
- Three groups: `us_native` (steelblue), `iscgem` (coral), `other` (gray)
- Stacked bar chart preferred
- Title: "ComCat Events by Decade and ID Prefix Source"
- Legend
- 300 DPI, saved to `output/`

---

## Test Suite: `test-case-a0.py`

All tests must pass. Tests run against the computed results JSON.

### Test cases

1. **Schema validation**: both CSV files contain exactly the 10 expected columns in any order
2. **ComCat event count**: equals 9,802
3. **ISC-GEM event count**: equals 9,210
4. **ComCat prefix counts sum**: `us_native + iscgem + other == 9802`
5. **ComCat `iscgem` prefix count**: equals 2,973 (±1 tolerance for rounding edge cases)
6. **Magnitude precision sums**: for each catalog, `one_decimal + two_decimal == event_count`
7. **Magnitude bin coverage**: all events accounted for in bins (sum of bin counts == event count, per catalog)
8. **Results JSON exists and is valid**: file is present, parses as JSON, contains all required top-level keys

---

## Whitepaper: `case-a0-whitepaper.md`

Follow the report writing rules (header/footer template, numbered sections).

**Title:** `Case A0: ComCat and ISC-GEM Catalog Comparison Reference`

**Sections:**

1. **Abstract** — Brief statement of purpose: this is a descriptive reference, not a validation exercise. States the key structural finding (hybrid composition of ComCat, magnitude precision difference).

2. **Data Sources** — Table of both catalogs with file paths, event counts, year range, mag range. Note that both files share an identical schema.

3. **Methodology** — Describe the three analyses performed: (a) population summary, (b) ID prefix classification and temporal distribution, (c) magnitude precision classification.

4. **Results**

   4.1 Population Summary — Reproduce the population table. Note ComCat has 592 more events (+6.4%).

   4.2 ComCat Hybrid Composition — Reproduce the ID prefix breakdown table. Embed `case-a0-comcat-prefix-temporal.png`. Describe the temporal pattern of `iscgem`-prefixed records (whether concentrated pre-1976 or distributed throughout).

   4.3 Magnitude Precision — Reproduce the precision table. Embed `case-a0-magnitude-distributions.png`. Describe the consequence for magnitude-based analyses near the M 6.0 threshold.

   4.4 Schema Compatibility — Note that both files use identical column structure; no structural incompatibility.

5. **Interpretation** — Objective summary: the two catalogs are not independent (30% overlap via ISC-GEM sourcing in ComCat). Legacy ComCat results are subject to lower magnitude precision. Do not editorialize beyond what the data shows.

6. **Limitations** — No event-level matching was performed; overlap percentage is estimated from ID prefix, not confirmed by coordinate/time matching. Magnitude precision classification relies on the stored float representation.

7. **References** — None required unless citing catalog documentation.

---

## Acceptance Criteria

- All tests in `test-case-a0.py` pass
- Both PNG images exist in `output/` at 300 DPI
- `case-a0-results.json` present and valid
- `case-a0-whitepaper.md` present with all 6 sections populated and both images embedded inline
