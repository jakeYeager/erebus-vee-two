# Case A0b Spec: Duplicate Detection and Cross-Catalog Event Accounting

## Status

Ready for implementation.

---

## Context

Case A0 established that ComCat contains 592 more events than ISC-GEM (+6.4%), that 30.3% of ComCat records carry an `iscgem` ID prefix (concentrated 85.3% in the 1950s–1960s), and that the M 6.0 bin is inflated in ComCat by +52% relative to ISC-GEM. This case investigates whether any of that excess reflects duplicate records within ComCat or systematic catalog-level differences, and translates findings into concrete data pipeline refinements.

**Prerequisite output:** `topic-adhoc/output/case-a0-results.json` must exist (produced by Case A0).

---

## Data Inputs

| File | Path |
| --- | --- |
| ComCat | `data/global-sets/comcat_global_6-9_1949-2021.csv` |
| ISC-GEM | `data/global-sets/iscgem_global_events.csv` |

Paths are relative to project root (`/Users/jake/Projects/Code/prod/erebus-vee-two/`).

Schema (both files): `usgs_id, usgs_mag, event_at, solaration_year, solar_secs, lunar_secs, midnight_secs, latitude, longitude, depth`

---

## Deliverables

All outputs go under `topic-adhoc/`.

| Type | Path |
| --- | --- |
| Analysis script | `src/case-a0b-analysis.py` |
| Visualization script | `src/visualization-case-a0b.py` |
| Test suite | `tests/test-case-a0b.py` |
| Results JSON | `output/case-a0b-results.json` |
| Whitepaper | `output/case-a0b-whitepaper.md` |
| Image: ID overlap | `output/case-a0b-id-overlap.png` |
| Image: cross-catalog match accounting | `output/case-a0b-event-accounting.png` |
| Image: unmatched events temporal | `output/case-a0b-unmatched-temporal.png` |

---

## Analysis Script: `case-a0b-analysis.py`

Performs three sequential analyses and writes all results to `case-a0b-results.json`.

### Analysis 1: Direct ID Cross-Reference

The 2,973 `iscgem`-prefixed ComCat records use the same ID scheme as the standalone ISC-GEM catalog. Compare them by exact string match on `usgs_id`.

Compute:
- Count of `iscgem`-prefixed ComCat IDs that appear verbatim in ISC-GEM (`id_matched`)
- Count that do not appear in ISC-GEM (`id_unmatched_in_iscgem`) — these are ComCat records attributed to ISC-GEM sourcing but not present in the current ISC-GEM file
- Count of ISC-GEM IDs that do NOT appear anywhere in ComCat (`iscgem_not_in_comcat`) — ISC-GEM events ComCat did not ingest
- Express each as count and percentage of the relevant population

### Analysis 2: Within-ComCat Duplicate Detection

Check whether any `us_native`-prefixed and `iscgem`-prefixed ComCat records represent the same physical event (different ID, same earthquake).

**Method:**
1. Split ComCat into two subsets: `us_native` (6,716 records) and `iscgem`-prefixed (2,973 records)
2. For each record in the smaller `iscgem` subset, find candidate matches in `us_native` using:
   - **Time window:** parse `event_at` as UTC datetime; absolute difference ≤ 60 seconds
   - **Distance:** Haversine distance ≤ 50 km (use `latitude`, `longitude`)
   - **Magnitude:** `|usgs_mag_a - usgs_mag_b|` ≤ 0.3
3. A record pair satisfying all three criteria is a **duplicate candidate**
4. Each `iscgem` record may match at most one `us_native` record (take the closest temporal match if multiple candidates exist)

**Haversine implementation:** implement in the script using the standard formula; do not rely on external geo libraries beyond `numpy` or `math`.

**Efficiency note:** Sort both subsets by `event_at` before matching. For each `iscgem` record, use binary search to find the time window in `us_native`, then apply spatial and magnitude filters only within that window. With 2,973 × 6,716 records and a 60-second window, the candidate set per record will typically be very small.

Compute:
- Total duplicate candidate pairs found
- As a percentage of the `iscgem`-prefixed subset
- Temporal distribution of duplicate candidates by decade (using `solaration_year`)
- Net ComCat effective event count after removing duplicates: `9802 - duplicate_count`

### Analysis 3: Full Cross-Catalog Proximity Matching

Match the full ComCat population (9,802) against the full ISC-GEM population (9,210) to categorize the 592-event difference.

**Method:** Same proximity criteria as Analysis 2 (60s, 50km, 0.3 mag). Use the same sorted time-window approach for efficiency.

Matching is one-to-one: once an ISC-GEM event is matched to a ComCat event, it cannot match another ComCat event (use greedy matching by closest temporal distance).

**Run at two tolerance levels** and report both:
- **Primary:** ±60s, 50km, 0.3 mag
- **Strict:** ±30s, 25km, 0.2 mag

For each tolerance level, compute:
- `matched`: events with a counterpart in the other catalog
- `comcat_only`: ComCat events with no ISC-GEM match
- `iscgem_only`: ISC-GEM events with no ComCat match
- Verify: `matched_comcat + comcat_only == 9802` and `matched_iscgem + iscgem_only == 9210`
- `comcat_only - iscgem_only` should approximate the 592-event gap net of any catalog-specific inclusions

Compute temporal distribution (by decade) and magnitude distribution (0.5-mag bins: 6.0–6.4, 6.5–6.9, 7.0–7.4, 7.5+) for both `comcat_only` and `iscgem_only` subsets. This identifies whether the gap is concentrated in a particular period or magnitude range.

---

## Output: `case-a0b-results.json`

```json
{
  "generated": "<ISO timestamp>",
  "prerequisites": {
    "case_a0_results": "topic-adhoc/output/case-a0-results.json"
  },
  "id_cross_reference": {
    "comcat_iscgem_prefix_count": 2973,
    "id_matched_in_iscgem": {"count": null, "pct_of_comcat_iscgem": null},
    "id_unmatched_in_iscgem": {"count": null, "pct_of_comcat_iscgem": null},
    "iscgem_ids_not_in_comcat": {"count": null, "pct_of_iscgem": null}
  },
  "within_comcat_duplicates": {
    "tolerances": {"time_seconds": 60, "distance_km": 50, "magnitude": 0.3},
    "duplicate_candidate_count": null,
    "pct_of_iscgem_prefix_subset": null,
    "net_effective_comcat_count": null,
    "temporal_distribution": {
      "1950s": null,
      "1960s": null,
      "1970s": null,
      "1980s_plus": null
    }
  },
  "cross_catalog_matching": {
    "primary": {
      "tolerances": {"time_seconds": 60, "distance_km": 50, "magnitude": 0.3},
      "matched": null,
      "comcat_only": null,
      "iscgem_only": null,
      "comcat_only_temporal": {},
      "iscgem_only_temporal": {},
      "comcat_only_by_mag_band": {},
      "iscgem_only_by_mag_band": {}
    },
    "strict": {
      "tolerances": {"time_seconds": 30, "distance_km": 25, "magnitude": 0.2},
      "matched": null,
      "comcat_only": null,
      "iscgem_only": null
    }
  }
}
```

---

## Visualization Script: `visualization-case-a0b.py`

Reads `output/case-a0b-results.json` and produces three PNG images.

### Image 1: `case-a0b-id-overlap.png`

Stacked bar chart showing ID-level accounting.

- Three bars: "ComCat iscgem-prefix" (total 2,973), "ISC-GEM full population" (9,210)
- Within each bar, show matched vs unmatched segments
- Colors: matched segment steelblue, unmatched segment light gray
- Title: "Direct ID Overlap: ComCat iscgem-prefixed Records vs ISC-GEM Catalog"
- Annotate each segment with count
- 300 DPI

### Image 2: `case-a0b-event-accounting.png`

Grouped bar chart showing event accounting at primary and strict tolerances.

- X groups: "ComCat Only", "Matched (both)", "ISC-GEM Only"
- Two bars per group: primary tolerance (steelblue), strict tolerance (lighter blue)
- Y-axis: event count
- Title: "Cross-Catalog Event Accounting: ComCat vs ISC-GEM"
- Annotate bars with counts
- 300 DPI

### Image 3: `case-a0b-unmatched-temporal.png`

Side-by-side bar chart showing temporal distribution of unmatched events (primary tolerance).

- X-axis: decade
- Two bar groups per decade: "ComCat Only" (steelblue), "ISC-GEM Only" (coral)
- Y-axis: event count
- Title: "Unmatched Events by Decade (Primary Tolerance)"
- Legend
- 300 DPI

---

## Test Suite: `test-case-a0b.py`

All tests must pass.

1. **Results JSON exists and is valid**: file present, parses as JSON, all top-level keys present
2. **ID accounting**: `id_matched + id_unmatched == 2973` (ComCat iscgem-prefix total)
3. **Within-ComCat duplicate accounting**: `duplicate_candidate_count` is a non-negative integer ≤ 2973
4. **Net effective count**: `net_effective_comcat_count == 9802 - duplicate_candidate_count`
5. **Cross-catalog primary accounting (ComCat)**: `matched + comcat_only == 9802`
6. **Cross-catalog primary accounting (ISC-GEM)**: `matched + iscgem_only == 9210`
7. **Strict tolerance ≤ primary**: `strict.matched <= primary.matched` (stricter criteria cannot match more events)
8. **All three PNG images exist** at expected paths

---

## Whitepaper: `case-a0b-whitepaper.md`

Follow report writing rules (header/footer template, numbered sections).

**Title:** `Case A0b: Duplicate Detection and Cross-Catalog Event Accounting`

**Sections:**

1. **Abstract** — States that this case extends A0 by performing event-level analysis to quantify how much of the 592-event difference between catalogs reflects genuine catalog divergence vs within-ComCat duplication. Summarizes whether significant internal duplication was found. Includes pipeline refinement recommendations.

2. **Data Sources** — Same as A0; reference Case A0 for detailed catalog description. Note that `case-a0-results.json` was used as an input.

3. **Methodology** — Describe the three analyses: (a) direct ID string matching, (b) within-ComCat spatial-temporal duplicate detection with specified tolerances, (c) full cross-catalog proximity matching at two tolerance levels. State the chosen tolerance values and rationale (seismological location uncertainty for M 6+ events is typically < 20 km for well-instrumented periods; ±60s is conservative; ±0.3 mag accounts for inter-scale conversion differences).

4. **Results**

   4.1 ID Cross-Reference — Reproduce ID overlap table. Embed `case-a0b-id-overlap.png`. Interpret what fraction of the `iscgem`-prefixed ComCat records are traceable to the current ISC-GEM file.

   4.2 Within-ComCat Duplicate Candidates — Report duplicate candidate count and percentage. Embed nothing here (temporal distribution shown in image 3). State whether duplicates are concentrated in the same pre-1980 period as the backfill.

   4.3 Cross-Catalog Event Accounting — Reproduce the accounting table for primary and strict tolerances. Embed `case-a0b-event-accounting.png`. Report `comcat_only` and `iscgem_only` counts. Embed `case-a0b-unmatched-temporal.png`. Describe whether the unmatched events are concentrated in a period or magnitude band.

5. **Interpretation** — Objective assessment only. Describe what the accounting reveals about the nature of the 592-event gap: Is it driven by ComCat capturing more events (broader completeness), ISC-GEM capturing different events (different selection criteria), duplicates, or some combination? Note any magnitude band where one catalog is systematically richer. Reference the M 6.0 bin inflation from A0 and whether the `comcat_only` events are disproportionately M 6.0.

   **5.1 Implications for Legacy G-K Declustering Results (Topics L4 and L5)**

   Topics L4 and L5 applied the Gardner-Knopoff (1974) declustering algorithm to this same ComCat population (n=9,802), designating 3,580 events (36.52%) as aftershocks and retaining 6,222 mainshocks. The A0b findings introduce two specific mechanisms by which ComCat data quality issues may have inflated the aftershock designation count, and a third mechanism by which they may have degraded the mainshock classification confidence. Each should be assessed against the A0b quantified results and stated objectively.

   **Mechanism 1: Within-catalog duplicate pairs as false aftershocks.**
   If ComCat contains duplicate records of the same physical event under two different IDs (one `us_native`, one `iscgem`-prefixed), G-K would process both records as independent earthquakes. Because the two records would share nearly identical coordinates and an event time separation of seconds to minutes — well within any G-K temporal window (minimum ~295 days for M 6.0) — G-K would classify the record with the lower magnitude or later timestamp as an aftershock of the other. Each confirmed within-ComCat duplicate pair represents one potential false-positive aftershock designation. State the Analysis 2 duplicate candidate count and express it as a fraction of the 3,580 designated aftershocks. If duplicates exceed ~1% of the aftershock population (~36 events), this constitutes a non-trivial source of inflation.

   **Mechanism 2: Near-threshold rounding artifacts at M 6.0 as spurious aftershock candidates.**
   Case A0 established that ComCat contains 643 more events at the M 6.0 bin than ISC-GEM (1,876 vs. 1,233), consistent with 1-decimal rounding capturing near-threshold events (~M 5.95–6.04 stored as M 6.0). These events, had they been measured at 2-decimal precision, may have fallen below the M 6.0 inclusion threshold entirely. When such a near-threshold event occurs in temporal and spatial proximity to a larger mainshock — as sub-threshold seismicity commonly does — G-K correctly classifies it as an aftershock relative to ComCat's inclusion criteria. However, if the event is a rounding artifact that ISC-GEM would exclude, its aftershock designation is cataloguing an event that should not be present. Cross-reference the A0b `comcat_only` magnitude distribution: if near-threshold `comcat_only` events (M 6.0 bin) are disproportionately concentrated in the post-mainshock temporal windows of larger events, this pattern would support Mechanism 2.

   **Mechanism 3: Backfill duplicates inflating the E1 classification confidence borderline fraction.**
   The L5 Case E1 analysis found that 46.61% of G-K mainshocks (2,900 of 6,222) fall within another mainshock's G-K exclusion window, yielding a classification confidence of 0.534 — substantially below the Nandan et al. (2024) reference threshold of ~0.67. The `iscgem`-prefixed backfill records are concentrated in the 1950s–1970s (85.3% of the 2,973 records, per A0 results). If any of these backfill records are near-duplicates of `us_native` records for the same events, they would appear as tight spatiotemporal clusters in the pre-1980 period. G-K, encountering such clusters, would designate one record per pair as a mainshock and the other as an aftershock, but the retained "mainshock" would necessarily be in close proximity to its near-duplicate — artificially elevating the borderline fraction. If the temporal distribution of within-ComCat duplicate candidates from Analysis 2 is concentrated pre-1980 (as expected given the backfill pattern), this is consistent with backfill duplication as a partial driver of the anomalously low E1 classification confidence.

   **Overall assessment:** Quantify what fraction of the 3,580 designated aftershocks could plausibly be affected by Mechanisms 1 and 2 combined. State whether the magnitude of the effect is material relative to the full aftershock population. Do not claim the solar signal suppression (60.2%, L5 Case E2) is definitively explained by false-positive aftershocks — the E5 finding that suppression is proportional to other lawful properties is not contradicted by a moderate inflation of the aftershock count. The appropriate framing is that the ComCat-based aftershock partition carries a quantifiable data quality uncertainty that was not accounted for in L4/L5, and that repeating the G-K declustering on the ISC-GEM catalog is the necessary validation step for any conclusion that depends on the mainshock/aftershock partition.

6. **Data Pipeline Refinement Recommendations** — Actionable suggestions derived from A0 and A0b findings. Address each finding explicitly:

   **R1 — Magnitude completeness threshold (from A0, M 6.0 bin inflation)**
   The M 6.0 bin in ComCat is inflated +52% relative to ISC-GEM (1,876 vs 1,233 events), consistent with 1-decimal rounding capturing near-threshold events. For analyses sensitive to the completeness threshold (b-value estimation, rate comparisons), apply a working threshold of M 6.1 when using ComCat to avoid contamination from rounding artifacts.

   **R2 — Magnitude precision flag (from A0)**
   Add a derived column `mag_precision` during ingestion (values: `"1d"`, `"2d"`) using the multiply-by-100 modulo method. This enables per-record precision awareness without removing data, and allows sensitivity analyses to be run on 2-decimal subsets only.

   **R3 — Source provenance column (from A0)**
   Add a derived column `id_source` (values: `"us_native"`, `"iscgem"`, `"other"`) derived from `usgs_id` prefix during ComCat ingestion. This allows any analysis to be stratified or filtered by source without re-deriving prefix logic ad hoc.

   **R4 — Within-catalog deduplication step (from A0b, quantify finding)**
   If duplicate candidates exceed 1% of the `iscgem`-prefixed subset (i.e., > ~30 events), add a deduplication step to the ComCat pipeline that removes `iscgem`-prefixed records that match a `us_native` record within primary tolerances. Retain the `us_native` record as authoritative for the post-1980 period. If duplicates are negligible (< 1%), document the finding and skip the deduplication step.

   **R5 — Catalog selection validation (from A0b)**
   The migration to ISC-GEM (already completed) is reinforced by these findings: ComCat is a composite source that introduces magnitude precision heterogeneity and potential duplication. Quantify how many `comcat_only` events (no ISC-GEM counterpart) fall outside the ISC-GEM selection criteria — if they are predominantly M 6.0 (rounding artifacts) or historical backfill anomalies, the migration carries no material loss of real seismic signal.

   **R6 — Year range metadata assertion (from A0)**
   Both files report a minimum year of 1950 despite the ComCat filename claiming 1949 coverage. Add a pipeline validation step that asserts `solaration_year.min()` matches the declared coverage range. Log a warning (not a failure) if the actual range differs from the filename metadata.

7. **Limitations** — Proximity matching cannot definitively confirm event identity — spatially coincident events within 60 seconds and 50 km could theoretically be distinct earthquakes (e.g., a large foreshock-mainshock pair). Duplicate detection is probabilistic, not deterministic. The ID cross-reference reflects the current state of both catalog files; catalog versions may differ from those used to build ComCat.

8. **References** — Reference Case A0 whitepaper.

---

## Acceptance Criteria

- All tests in `test-case-a0b.py` pass
- All three PNG images exist in `output/` at 300 DPI
- `case-a0b-results.json` present, valid, with no `null` values remaining
- `case-a0b-whitepaper.md` present with all 8 sections populated and all three images embedded inline
- Pipeline recommendations R1–R6 all addressed with quantified findings where applicable
