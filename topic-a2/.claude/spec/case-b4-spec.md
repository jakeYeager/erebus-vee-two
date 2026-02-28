# Case B4: Depth Stratification — Surface Loading Penetration Test

**Status:** Pending

**Intent statement:**
Case B4 stratifies the ISC-GEM catalog by focal depth into four bands — shallow (0–20 km), mid-crustal (20–70 km), intermediate (70–300 km), and deep (>300 km if sample permits) — and independently computes solar-phase bin statistics for each band. The original motivation was that hydrological loading primarily affects the shallow crust, so a signal persisting at depth would rule out surface loading as sufficient. However, Zhan & Shearer (2015) found that 70% of deep-focus M>7 earthquakes (depth >500 km) occurred in April–September, a suggestive seasonal signal with no known surface mechanism, which directly inverts the original prediction. The literature now hints that large deep events may show more seasonality rather than less, making depth stratification a multi-directional test: if the signal decreases monotonically with depth, surface loading is implicated; if the signal is present at all depths including >300 km, a geometric or deep forcing is implied; if deep events show a different phase than shallow events, the Zhan & Shearer pattern may be present. The deep band (>300 km) is included if sample size permits meaningful statistics (n ≥ 100 as minimum threshold).

**Relationship to prior topics:**
Topic L4 confirmed that the `depth` column is present in the ISC-GEM catalog schema. Topic L5 showed that declustering removes events at all depth ranges, but its effect on depth-stratified signal is unknown and left to this case. B4 does not require declustered catalogs — the planning doc uses the raw ISC-GEM catalog — but the A4 results should be cross-referenced in the whitepaper. The Zhan & Shearer (2015) finding of a deep-event seasonal pattern specifically covers M>7, so the B4 deep-band result should note the magnitude constraint when comparing.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/global-sets/iscgem_global_events.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `latitude`, `longitude`, `depth` |

Depth bands derived at analysis time from `depth` (units: km):
- Band 1 shallow: `0 <= depth < 20`
- Band 2 mid-crustal: `20 <= depth < 70`
- Band 3 intermediate: `70 <= depth < 300`
- Band 4 deep: `depth >= 300` (include in analysis only if n ≥ 100; log a warning if below threshold)

Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year [confirm before running: consistent with prior cases — verify uniform Julian constant vs per-year values].

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-b4-analysis.py` — main analysis: depth stratification, per-band bin statistics, monotonicity test, Zhan & Shearer comparison, writes results JSON
- `src/visualization-case-b4.py` — generates all PNG figures
- `tests/test-case-b4.py` — test suite
- `output/case-b4-results.json` — per-band statistics, trend analysis, prediction matching
- `output/case-b4-whitepaper.md` — methodology, results, cross-topic comparison
- `output/case-b4-binplots.png` — multi-panel bin distributions (one per depth band) at k=24
- `output/case-b4-depth-trend.png` — Cramér's V and Rayleigh R vs depth band (trend plot)
- `output/case-b4-deep-events-phase.png` — detailed phase distribution for deep (>300 km) events with April–September highlight

---

## 1. Environment and data loading

In `src/case-b4-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define `RAW_PATH = BASE_DIR.parent / "data" / "global-sets" / "iscgem_global_events.csv"`
- Load catalog; assert n=9210
- Define depth band boundaries:
  ```python
  DEPTH_BANDS = [
      {"label": "shallow_0-20km",    "min": 0,   "max": 20},
      {"label": "midcrustal_20-70km","min": 20,  "max": 70},
      {"label": "intermediate_70-300km", "min": 70,  "max": 300},
      {"label": "deep_300km+",       "min": 300, "max": 9999},
  ]
  ```
- Create subsets; log band sizes
- Check for NaN in `depth` column: log count of NaN events; exclude from depth-band analysis; record `n_depth_null`
- Assert sum of band sizes + n_depth_null == 9210
- Log warning if deep band n < 100; still include in analysis but flag in results
- Compute phase for each subset using Julian year constant

---

## 2. Per-band bin statistics

For each depth band at each of k=16, 24, 32:

1. Compute observed bin counts, expected counts, chi-square, p-value, Cramér's V, Rayleigh R, Rayleigh p, mean phase (same procedure as A4)
2. Elevated bins at 1-SD threshold; elevated phase intervals (same as A4)
3. Compare elevated intervals to A1b baseline (same matching logic as A4 Sub-analysis B)

Output structure under key `"band_stats"` in results JSON:
```json
{
  "band_stats": {
    "shallow_0-20km": {
      "n": int,
      "sufficient_n": bool,
      "k16": {"chi2": float, "p_chi2": float, "cramer_v": float, "rayleigh_R": float,
               "p_rayleigh": float, "mean_phase": float, "bin_counts": [...], "elevated_intervals": [...]},
      "k24": {...}, "k32": {...}
    },
    "midcrustal_20-70km": {...},
    "intermediate_70-300km": {...},
    "deep_300km+": {...}
  },
  "n_depth_null": int
}
```

---

## 3. Depth trend analysis

Using results at k=24:

1. Extract Cramér's V for each band (ordered shallow → deep), excluding bands with insufficient n
2. Compute Spearman rank correlation of Cramér's V vs depth band index (1–4); record rho and p-value
3. Classify monotonicity:
   - "decreasing with depth": Spearman rho < -0.5 (signal weakens going deeper)
   - "increasing with depth": Spearman rho > 0.5 (signal strengthens going deeper)
   - "non-monotonic": |rho| ≤ 0.5 (no clear direction)
4. Identify which bands are individually significant (p_chi2 < 0.05 at k=24)
5. Check whether deep band (>300 km) is significant; record as separate flag

**Prediction matching:**
- **Surface loading hypothesis**: signal strongest at 0–20 km; absent at >70 km → evaluate: shallow significant AND intermediate not significant
- **Geometric/deep forcing hypothesis**: signal present at all depths including >300 km → evaluate: all bands significant including deep
- **Zhan & Shearer pattern**: deep events show seasonality; check whether deep mean phase is in April–September range (phase 0.23–0.67 for a year starting Jan 1) → evaluate: deep band significant AND mean_phase in [0.23, 0.67]

Output under key `"trend_analysis"` in results JSON:
```json
{
  "trend_analysis": {
    "cramer_v_by_band": {"shallow_0-20km": float, "midcrustal_20-70km": float, "intermediate_70-300km": float, "deep_300km+": float},
    "spearman_rho": float,
    "spearman_p": float,
    "monotonicity_classification": "decreasing | increasing | non-monotonic",
    "significant_bands": [...],
    "deep_band_significant": bool,
    "deep_mean_phase": float,
    "deep_mean_phase_in_apr_sep": bool,
    "prediction_support": {
      "surface_loading": "supported | partially | not supported",
      "geometric_deep_forcing": "supported | partially | not supported",
      "zhan_shearer_pattern": "supported | partially | not supported"
    }
  }
}
```

---

## 4. Visualizations

In `src/visualization-case-b4.py`:

**Figure 1 — Multi-panel bin distributions** (`output/case-b4-binplots.png`):
- 4-panel layout (2×2), one per depth band, all at k=24
- Each panel: horizontal steelblue bar chart; dashed expected-count line; 1-SD threshold dashed line
- Gray shaded bands for A1b baseline intervals (phases 0.1875–0.25, 0.625–0.656, 0.875–0.917)
- Annotate each panel with band label, n, χ², p-value, Cramér's V; if n < 100, add "(low n)" annotation
- Shared x-axis: phase fraction 0–1; label equinox and solstice positions
- 300 DPI

**Figure 2 — Depth trend plot** (`output/case-b4-depth-trend.png`):
- Same format as A3 Figure 2: x-axis = depth band (shallow to deep as ordinal positions 1–4); dual y-axes for Cramér's V (steelblue) and Rayleigh R (orange)
- Mark significance (p < 0.05) with filled circles; non-significant with open circles; bands with n < 100 marked with 'x' symbols
- Annotate Spearman rho and monotonicity classification
- 300 DPI

**Figure 3 — Deep events phase detail** (`output/case-b4-deep-events-phase.png`):
- Single panel for the deep band (>300 km events) at k=24
- Shade the April–September range (phase 0.23–0.67) in light orange to reference the Zhan & Shearer (2015) finding
- Same bin-distribution format; annotate n, χ², p, mean phase; add note "Zhan & Shearer (2015): 70% of deep M>7 events in April–September"
- 300 DPI

---

## 5. Test suite

In `tests/test-case-b4.py`:

- `test_catalog_load`: assert n=9210
- `test_depth_band_partition`: assert sum of all band sizes + n_depth_null == 9210
- `test_depth_band_sizes_positive`: assert all four bands have n > 0; log warning if deep band < 100
- `test_no_overlap_between_bands`: assert no `usgs_id` appears in more than one depth band
- `test_phase_range`: assert all phases in [0.0, 1.0)
- `test_chi_square_all_bands`: assert chi2 and p_chi2 are finite for all bands at k=16, 24, 32
- `test_cramer_v_range`: assert Cramér's V in [0.0, 1.0] for all bands
- `test_spearman_rho_range`: assert Spearman rho in [-1.0, 1.0]
- `test_monotonicity_classification_valid`: assert monotonicity_classification is "decreasing", "increasing", or "non-monotonic"
- `test_zhan_shearer_phase_range`: assert `deep_mean_phase_in_apr_sep` is bool; assert April–September phase range check uses [0.23, 0.67] (Jan 1 = phase 0.0)
- `test_prediction_support_valid`: assert all prediction_support values are one of "supported", "partially", "not supported"
- `test_results_json_keys`: load `output/case-b4-results.json`; assert keys "band_stats" and "trend_analysis" present

---

## 6. Whitepaper

In `output/case-b4-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (does the solar signal persist at depth?), four depth bands, monotonicity result, comparison to Zhan & Shearer pattern, and mechanistic implication
2. **Data Source** — ISC-GEM catalog (n=9,210); report band sizes; note depth null count; note ISC-GEM depth precision
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1, `rules/data-handling.md`)
   - 3.2 Depth band definitions and physical rationale
   - 3.3 Chi-square, Rayleigh, Cramér's V per band
   - 3.4 Depth trend analysis: Spearman rank correlation and monotonicity classification
   - 3.5 Prediction framework: surface loading, geometric/deep forcing, Zhan & Shearer
4. **Results**
   - 4.1 Per-band distributions: embed Figure 1; tabulate chi2, p, Cramér's V, Rayleigh R for all four bands at k=24
   - 4.2 Depth trend: embed Figure 2; state Spearman rho, monotonicity classification, and which bands are individually significant
   - 4.3 Deep events: embed Figure 3; compare deep mean phase to April–September range; report Zhan & Shearer pattern match
   - 4.4 Prediction matching: evaluate all three predictions
5. **Cross-Topic Comparison** — compare to Zhan & Shearer (2015) finding (deep M>7 events 70% in April–September); compare to Dutilleul et al. (2021) Parkfield depth-signal analysis if relevant; note depth distribution differences between ISC-GEM and catalogs used in tidal literature
6. **Interpretation** — state mechanism implication from monotonicity and deep-band results; acknowledge the inverted prediction (deep signal is not ruled out by literature); maintain objectivity
7. **Limitations** — deep band has small sample (n < 100 likely); ISC-GEM depth estimates have varying precision, particularly for historical events; declustering not applied
8. **References** — Zhan & Shearer (2015), Dutilleul et al. (2021), Johnson et al. (2017), Adhoc A4, Case 3A

---

## 7. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case B4: Depth Stratification — Surface Loading Penetration Test
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [Cramér's V per depth band at k=24; monotonicity classification; deep band n and significance; Zhan & Shearer match; prediction support summary]
  ```
- Update Case B4 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
