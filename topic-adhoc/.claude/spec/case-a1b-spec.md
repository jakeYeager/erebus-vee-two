# Case A1b Spec: Elevated Bin Event Characterization and Declustering Implications

## Status

Ready for implementation.

---

## Context

Case A1 established that the `solar_secs` chi-square signal is robust across bin counts (Bonferroni-significant at k=16, 24, 32). This case drills into the events driving those elevated bins to determine whether they show temporal and/or spatial clustering — which would be consistent with aftershock sequence residuals that G-K failed to capture — or whether they are dispersed, consistent with independent events genuinely preferring those solar phases.

The G-K windows for M≥6.0 events at global scale may be too permissive: spatial windows of ~49–493 km and temporal windows of ~295–2,117 days span ranges that mix genuine aftershock sequences with spatially and temporally proximate independent events. If elevated-bin events show clustering tighter than these windows, the G-K parameterization is not the limiting factor; if they are dispersed beyond these windows, G-K cannot capture them by design and a different algorithm geometry is needed.

---

## Data Inputs

| File | Path |
| --- | --- |
| ISC-GEM | `data/global-sets/iscgem_global_events.csv` |
| Plate boundaries | `lib/pb2002_boundaries.dig` |

Paths relative to project root (`/Users/jake/Projects/Code/prod/erebus-vee-two/`).

**`pb2002_boundaries.dig` format:**
- Header line: `PLATE_A-PLATE_B  source` (text)
- Coordinate lines: `+longitude, +latitude` (decimal degrees, scientific notation) — note lon-first order
- Segment separator: `*** end of line segment ***`

**Prerequisite output:** `topic-adhoc/output/case-a1-results.json` must exist.

---

## Deliverables

All outputs go under `topic-adhoc/`.

| Type | Path |
| --- | --- |
| Analysis script | `src/case-a1b-analysis.py` |
| Visualization script | `src/visualization-case-a1b.py` |
| Test suite | `tests/test-case-a1b.py` |
| Results JSON | `output/case-a1b-results.json` |
| Whitepaper | `output/case-a1b-whitepaper.md` |
| Image: phase coherence | `output/case-a1b-phase-coherence.png` |
| Image: temporal and magnitude | `output/case-a1b-temporal-magnitude.png` |
| Image: spatial | `output/case-a1b-spatial.png` |

---

## Analysis Script: `case-a1b-analysis.py`

### Phase Normalization (reuse from A1)

Apply the same phase normalization as Case A1 to all ISC-GEM events:
- `solar_secs`: divide by actual calendar year length (leap/non-leap) for each `solaration_year`
- All events should produce phase values in `[0, 1)`. Clamp any values ≥ 1.0 to 0.999999 (boundary effect, consistent with A1).

### Analysis 1: Phase Coherence and Combined Elevated Phase Set

For each k ∈ {16, 24, 32}:
1. Bin all 9,210 events by `floor(phase * k)`
2. Compute expected count = 9210 / k
3. Compute deviation = observed − expected for each bin
4. Select the top 3 bins by positive deviation
5. Record the phase range for each: `[i/k, (i+1)/k)`

Build the **combined elevated phase set**:
- Collect all phase ranges that appear in the top-3 for ≥2 of the 3 bin counts
- Merge overlapping ranges into contiguous intervals
- Extract all events whose normalized phase falls within any interval in this set
- De-duplicate by `usgs_id` — each event appears at most once regardless of how many bin counts it was elevated in

Record:
- The top-3 bins and phase ranges for each k
- The combined phase intervals (merged)
- The count of events in the combined elevated set (`n_elevated`)
- The count of full catalog events in each combined interval (for a baseline comparison: if a phase range covers 12% of the year, it should contain ~12% of all events under the null)

### Analysis 2: Temporal Distribution

**2a. Year distribution**

For elevated-bin events and for the full catalog:
- Count events per year (1950–2021)
- Compute the elevated set's yearly fraction relative to full catalog (elevated_year_count / total_year_count)
- If this fraction is roughly constant across years: temporally uniform. If concentrated in specific years: investigate those years.

**2b. Inter-event intervals (IEI) within elevated-bin population**

Sort elevated-bin events by `event_at` (parse as UTC datetime). Compute consecutive time differences in days.

For a baseline comparison, draw 1,000 random samples of size `n_elevated` from the full catalog, compute IEI for each, and record the 2.5th, 50th, 97.5th percentile of the median IEI across samples. This gives a null distribution for the expected median IEI of a random subset of this size.

Record:
- Median, mean, 10th and 90th percentile IEI for elevated-bin events
- The null distribution 2.5th/50th/97.5th percentile of median IEI
- Whether the elevated-set median IEI falls outside the null 95% interval

### Analysis 3: Spatial Distribution

**3a. Nearest-neighbor distance**

For each elevated-bin event, compute the Haversine distance to its nearest other elevated-bin event. Record the distribution (median, 25th, 75th percentile).

For a baseline comparison, draw 1,000 random samples of size `n_elevated` from the full catalog and compute the same nearest-neighbor statistics. Record the null distribution 2.5th/50th/97.5th percentile of median NN distance.

**3b. Plate boundary proximity**

Parse `pb2002_boundaries.dig`:
- Skip header lines (non-numeric) and separator lines
- Parse coordinate lines as `(longitude, latitude)` pairs
- Collect all coordinates as a flat array of (lon, lat) points representing boundary vertices

For each event in the **full catalog** and each event in the **elevated-bin set**, compute the minimum Haversine distance to any boundary vertex in the parsed array. (Vertex-based proximity is an approximation of segment distance; it is sufficient at the scale of this analysis given boundary vertex density.)

Classify each event by proximity:
- `near_boundary`: distance ≤ 100 km
- `transitional`: 100 km < distance ≤ 300 km
- `intraplate`: distance > 300 km

Compare the proximity distribution (% in each class) between elevated-bin events and the full catalog. A significant difference would indicate spatial structure in the elevated-bin population.

**3c. Latitude banding**

Count events per 30° latitude band (90°S–60°S, 60°S–30°S, 30°S–0°, 0°–30°N, 30°N–60°N, 60°N–90°N) for both elevated-bin and full catalog. Express as percentage of each population in each band.

### Analysis 4: Magnitude Distribution

Compare magnitude distributions using 0.5-magnitude bins: 6.0–6.4, 6.5–6.9, 7.0–7.4, 7.5+.

For each bin: count and percentage for elevated-bin events vs full catalog. Express the elevated set's magnitude distribution relative to the full catalog to identify any magnitude skew.

### Analysis 5: Declustering Window Estimation

Summarize the observed clustering footprint and compare to G-K standard windows.

**G-K reference windows for M6.0 (Gardner & Knopoff 1974):**
- Spatial: ~49 km
- Temporal: ~295 days

Compute from the elevated-bin population:
- Spatial footprint: 25th, 50th, 75th percentile of nearest-neighbor distance
- Temporal footprint: 25th, 50th, 75th percentile of IEI

For each G-K reference magnitude level (M6.0, M7.0, M8.0), state whether the observed elevated-bin clustering falls within or beyond that window.

Propose a data-informed window estimate: a spatial radius (km) and temporal window (days) that would capture the clustering observed in the elevated-bin population. Clearly label this as a quantified reference point derived from observed metrics, not a validated declustering algorithm.

---

## Output: `case-a1b-results.json`

```json
{
  "generated": "<ISO timestamp>",
  "data_source": "data/global-sets/iscgem_global_events.csv",
  "boundary_file": "lib/pb2002_boundaries.dig",
  "event_count": 9210,
  "phase_coherence": {
    "k16": {
      "top3_bins": [],
      "top3_phase_ranges": [],
      "top3_deviations": []
    },
    "k24": { ... },
    "k32": { ... },
    "combined_elevated_intervals": [],
    "n_elevated": null,
    "elevated_pct_of_catalog": null,
    "expected_pct_under_null": null
  },
  "temporal": {
    "elevated_iei_days": {
      "median": null, "mean": null, "p10": null, "p90": null
    },
    "null_iei_median_ci": {
      "p2_5": null, "p50": null, "p97_5": null
    },
    "elevated_median_outside_null_ci": null
  },
  "spatial": {
    "elevated_nn_km": {
      "p25": null, "p50": null, "p75": null
    },
    "null_nn_median_ci": {
      "p2_5": null, "p50": null, "p97_5": null
    },
    "elevated_nn_outside_null_ci": null,
    "boundary_proximity": {
      "elevated": {
        "near_boundary_pct": null,
        "transitional_pct": null,
        "intraplate_pct": null
      },
      "full_catalog": {
        "near_boundary_pct": null,
        "transitional_pct": null,
        "intraplate_pct": null
      }
    },
    "latitude_bands": {
      "elevated": {},
      "full_catalog": {}
    }
  },
  "magnitude": {
    "elevated": { "6.0-6.4": null, "6.5-6.9": null, "7.0-7.4": null, "7.5+": null },
    "full_catalog": { "6.0-6.4": null, "6.5-6.9": null, "7.0-7.4": null, "7.5+": null }
  },
  "declustering_window_estimate": {
    "observed_spatial_p50_km": null,
    "observed_spatial_p75_km": null,
    "observed_temporal_p50_days": null,
    "observed_temporal_p75_days": null,
    "gk_reference_m6": { "spatial_km": 49, "temporal_days": 295 },
    "gk_reference_m7": { "spatial_km": 156, "temporal_days": 790 },
    "gk_reference_m8": { "spatial_km": 493, "temporal_days": 2117 },
    "proposed_window": {
      "spatial_km": null,
      "temporal_days": null,
      "basis": null
    }
  }
}
```

---

## Visualization Script: `visualization-case-a1b.py`

Reads `output/case-a1b-results.json` and produces three PNG images.

### Image 1: `case-a1b-phase-coherence.png`

Three-row panel showing elevated bin positions at k=16, 24, 32 on a phase axis [0, 1).

- Each row is a bin count. Draw a horizontal axis from 0 to 1.
- Mark each bin as a rectangle; fill elevated (top-3) bins steelblue, all others light gray
- Add a fourth row below showing the combined elevated phase intervals (merged, filled red)
- Annotate approximate calendar month positions on the x-axis (Jan=0.0, Apr=0.25, Jul=0.5, Oct=0.75)
- Title: "Solar Phase Elevated Bins by Bin Count"
- 300 DPI

### Image 2: `case-a1b-temporal-magnitude.png`

Two-panel figure (side by side).

**Left panel: IEI distribution**
- Overlapping step histograms: elevated-bin events (steelblue) vs a single representative random sample of equal size (gray)
- X-axis: IEI in days (log scale recommended given the wide range)
- Y-axis: count
- Annotate median IEI for each distribution
- Mark the null 95% CI band for the median as a shaded region
- Title: "Inter-Event Intervals: Elevated Bins vs Random Baseline"

**Right panel: Magnitude distribution**
- Side-by-side bar chart (0.5-magnitude bins)
- Elevated-bin events: steelblue; full catalog: light gray
- Express as percentage of each population (not raw count) so populations of different sizes are comparable
- Title: "Magnitude Distribution: Elevated Bins vs Full Catalog"

300 DPI

### Image 3: `case-a1b-spatial.png`

Three-panel figure.

**Left panel: Geographic scatter**
- All ISC-GEM events: light gray dots, size 1
- Elevated-bin events: steelblue dots, size 3, alpha 0.6
- Simple equirectangular projection (latitude vs longitude)
- Title: "Elevated Bin Events: Geographic Distribution"

**Center panel: Nearest-neighbor distance distribution**
- Overlapping step histograms: elevated-bin NN distances (steelblue) vs null distribution median (gray dashed line) with 95% CI band shaded
- X-axis: distance in km
- Annotate the 50th percentile for elevated-bin events
- Mark G-K M6.0 spatial window (49 km) as a vertical dashed red line
- Title: "Nearest-Neighbor Distance: Elevated Bins vs Null"

**Right panel: Plate boundary proximity**
- Stacked horizontal bar chart: three proximity classes (near_boundary, transitional, intraplate)
- Two bars: elevated-bin events (steelblue shades) vs full catalog (gray shades)
- X-axis: percentage of population
- Title: "Boundary Proximity: Elevated Bins vs Full Catalog"

300 DPI

---

## Test Suite: `test-case-a1b.py`

All tests must pass.

1. **Results JSON exists and is valid**: file present, all top-level keys present, no `null` values
2. **Phase coherence — event count**: `n_elevated` is a positive integer < 9210
3. **Phase coherence — de-duplication**: recompute elevated events from scratch; verify count matches `n_elevated` and no `usgs_id` appears more than once
4. **Phase coherence — combined interval coverage**: `expected_pct_under_null` equals the total fraction of [0,1) covered by the combined elevated intervals (within floating point tolerance)
5. **Temporal — IEI sum**: the number of IEI values equals `n_elevated - 1` (one fewer than the event count)
6. **Spatial — NN count**: the number of nearest-neighbor distances equals `n_elevated`
7. **Spatial — boundary proximity sums**: for both elevated and full catalog, `near_boundary_pct + transitional_pct + intraplate_pct ≈ 100.0` (within 0.1)
8. **Magnitude sums**: for both elevated and full catalog, sum of magnitude band counts equals respective population size
9. **Declustering proposal populated**: `proposed_window.spatial_km` and `proposed_window.temporal_days` are positive numbers and `basis` is a non-empty string
10. **All three PNG images exist** at expected paths

---

## Whitepaper: `case-a1b-whitepaper.md`

Follow report writing rules (header/footer template, numbered sections).

**Title:** `Case A1b: Elevated Bin Event Characterization and Declustering Implications`

**Sections:**

1. **Abstract** — States the purpose: characterize the events driving the `solar_secs` elevated bins identified in A1 to determine whether they are temporally/spatially clustered (sequence residuals) or dispersed (genuine phase-preferring events). Summarizes the key findings and the declustering window estimate.

2. **Data Sources** — ISC-GEM catalog (n=9,210). PB2002 plate boundary file (`lib/pb2002_boundaries.dig`) used for boundary proximity classification. Note that boundary proximity uses vertex-based distance as an approximation of segment distance.

3. **Methodology**

   3.1 Phase Coherence and Combined Elevated Set — Describe top-3 bin selection, phase range mapping, and the ≥2-of-3 union criterion. State the combined elevated interval(s) and total event count.

   3.2 Temporal Analysis — Describe IEI computation within elevated-bin population and the null comparison via 1,000 random size-matched samples.

   3.3 Spatial Analysis — Describe nearest-neighbor distance computation, null comparison method, boundary proximity classification (thresholds: 100 km, 300 km), and latitude banding.

   3.4 Magnitude Analysis — Describe the 0.5-magnitude bin comparison.

4. **Results**

   4.1 Phase Coherence — Embed `case-a1b-phase-coherence.png`. State the combined elevated phase intervals and their approximate calendar date equivalents. State whether the top-3 bins are coherent across k=16, 24, 32.

   4.2 Temporal — Embed `case-a1b-temporal-magnitude.png` (left panel discussion here). Report the elevated-bin IEI statistics and whether the median IEI falls outside the null 95% CI.

   4.3 Spatial — Embed `case-a1b-spatial.png`. Report nearest-neighbor statistics and null comparison. Report boundary proximity percentages for both populations. Report latitude band distribution.

   4.4 Magnitude — Discuss the right panel of `case-a1b-temporal-magnitude.png`. State whether the elevated-bin population shows magnitude skew relative to the full catalog.

5. **Interpretation** — Objective assessment of all four analyses together. Address the core question: do the elevated-bin events look like sequence residuals (clustered in time and space, lower magnitude) or like independent events (dispersed, magnitude-uniform)? State explicitly which evidence points in each direction without reconciling into a single definitive conclusion if the evidence is mixed.

6. **Declustering Implications** — State the observed spatial and temporal clustering metrics. Compare to G-K windows at M6.0, M7.0, M8.0. State the proposed window parameter estimate and its basis. Explicitly note this is a data-informed reference point derived from the observed clustering footprint of n_elevated events, not a validated replacement for G-K — further validation against known aftershock sequences would be required before applying as a declustering criterion.

7. **Limitations** — Boundary proximity uses boundary vertex distance as a proxy for segment distance; accuracy degrades for coarser boundary segments. The combined elevated phase set is defined by a ≥2-of-3 coherence criterion; events coherent at only one bin count are excluded. IEI analysis captures internal temporal clustering only and does not test whether elevated-bin events cluster with non-elevated events.

8. **References** — Case A1 (binning sensitivity baseline). Bird, P. (2003). An updated digital model of plate boundaries. *Geochemistry, Geophysics, Geosystems*, 4(3). Gardner, J.K. & Knopoff, L. (1974). Is the sequence of earthquakes in Southern California, with aftershocks removed, Poissonian? *BSSA* 64(5), 1363–1367.

---

## Acceptance Criteria

- All 10 tests in `test-case-a1b.py` pass
- All three PNG images exist in `output/` at 300 DPI
- `case-a1b-results.json` present, valid, no `null` values
- `case-a1b-whitepaper.md` present with all 8 sections populated and all three images embedded inline
- The proposed declustering window in the results JSON has a populated `basis` field explaining how the values were derived from the observed metrics
