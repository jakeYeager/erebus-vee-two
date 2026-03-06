# Case A3.B3: Ocean/Coast Sequential Threshold Sensitivity

**Status:** Complete

**Source case:** A2.B2
**Reference case:** A2.B3 (GCMT join reference)

**Intent statement:**
Case A3.B3 addresses the core sensitivity finding of A2.B2: the oceanic subset (GSHHG primary) just misses significance at p=0.061, yet its Cramér's V (0.0276) exceeds the continental V (0.0245) — the significance gap is entirely a sample-size artifact. The transitional zone (50–200 km offshore) is the most significant class, and the PB2002 broader oceanic definition flips the oceanic result to significant (p=5.24×10⁻³). The A2.B2 conclusion rests on a single boundary drawn at 200 km offshore that has not been systematically tested.

A3.B3 sequentially tightens the transitional-to-oceanic boundary from 200 km down to 25 km in 25 km increments (8 threshold steps), recomputing class sizes and chi-square statistics at each step. The goal is to identify the specific threshold at which the solar-phase signal migrates from the transitional class into the oceanic class — or fails to — directly quantifying the geographic extent of the mechanism.

Additionally, A3.B3 computes a per-event **distance-to-subduction-zone** metric using `PB2002_steps.dat` (SUB-type boundary segments), enabling a cross-tabulation of ocean classification with subduction proximity. This tests whether the transitional zone's signal is driven by subduction zone geometry specifically, or by coastal proximity in general. The GCMT mechanism column provides independent validation of the subduction proximity metric.

**Relationship to prior cases:**
A2.B2 established the baseline three-class chi-square results (oceanic p=0.061, transitional p=2.67×10⁻⁴, continental p=4.61×10⁻⁴) under GSHHG at k=24 with inner threshold 50 km and outer threshold 200 km. A3.B3 holds the inner threshold fixed at 50 km and sweeps the outer threshold. A2.B3 provided the GCMT focal mechanism join used here for subduction proxy validation.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |
| GSHHG classification (primary) | `data/iscgem/plate-location/ocean_class_gshhg_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |
| PB2002 classification (secondary) | `data/iscgem/plate-location/ocean_class_pb2002_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |
| GCMT focal mechanism join | `data/iscgem/focal-mechanism/focal_join_global.csv` | 9,210 | `usgs_id`, `mechanism`, `rake`, `match_confidence` |
| PB2002 boundary steps | `lib/PB2002_steps.dat` | 5,819 rows | row_idx, plate_pair, lon1, lat1, lon2, lat2, ..., boundary_type |
| PB2002 boundary polylines | `lib/pb2002_boundaries.dig` | 229 segments | plate_pair header + lon/lat coordinate rows |

**GSHHG baseline counts (T_outer=200 km, T_inner=50 km fixed):**
- continental (dist ≤ 50 km): 3,799
- transitional (50 < dist ≤ 200 km): 3,459
- oceanic (dist > 200 km): 1,952

**A1b baseline intervals (k=24 bin mapping):**
- Interval 1: bins 4–5 (phase [0.1667, 0.2500)) — March equinox (~0.19–0.25)
- Interval 2: bin 15 (phase [0.6250, 0.6667)) — mid-August (~0.625–0.656)
- Interval 3: bin 21 (phase [0.8750, 0.9167)) — late-November (~0.875–0.917)

Phase normalization: `phase = (solar_secs / 31_557_600.0) % 1.0` (Julian year constant).

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a3/`
- `LIB_DIR = BASE_DIR.parent / "lib"`
- All output paths: `BASE_DIR / "output" / ...`
- All data paths: `BASE_DIR.parent / "data" / "iscgem" / ...`

**Planned outputs:**
- `src/case-a3-b3-analysis.py` — threshold sweep, subduction proximity computation, cross-tabulations, GCMT validation, writes results JSON
- `src/visualization-case-a3-b3.py` — generates all PNG figures (uses cartopy)
- `tests/test-case-a3-b3.py` — test suite (all tests must pass)
- `output/case-a3-b3-results.json` — per-threshold stats, subduction proximity metrics, cross-tabulations
- `output/case-a3-b3-whitepaper.md` — methodology, results, interpretation
- `output/case-a3-b3-threshold-sweep.png` — chi-square p + Cramér's V + class n across 8 threshold steps
- `output/case-a3-b3-global-map.png` — Cartopy global map with event classification + SUB boundary overlay
- `output/case-a3-b3-subduction-crosstab.png` — cross-tabulation of ocean class × subduction proximity at each threshold
- `output/case-a3-b3-binplots.png` — bin distributions at baseline (T=200) and the key threshold where oceanic becomes significant
- `output/case-a3-b3-region-maps.png` — 3-panel regional detail maps (Japan, Philippines, New Zealand) showing coastline resolution and event classification at baseline

---

## 1. Environment and data loading

In `src/case-a3-b3-analysis.py`:

- Imports: `pandas`, `numpy`, `scipy.stats`, `scipy.spatial`, `pyproj`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`; `LIB_DIR = BASE_DIR.parent / "lib"`
- Define paths:
  ```python
  RAW_PATH      = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
  GSHHG_PATH    = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
  PB2002_PATH   = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_pb2002_global.csv"
  FOCAL_PATH    = BASE_DIR.parent / "data" / "iscgem" / "focal-mechanism" / "focal_join_global.csv"
  STEPS_PATH    = LIB_DIR / "PB2002_steps.dat"
  ```
- Load raw catalog; parse `event_at` as UTC; assert n=9210
- Load GSHHG and PB2002 classification files; assert n=9210 each
- Load focal join file; assert n=9210
- Merge all onto raw catalog by `usgs_id`; assert no NaN in `dist_to_coast_km` columns after merge
- Compute `phase = (solar_secs / 31_557_600.0) % 1.0`
- Define constants:
  ```python
  K_BINS = 24
  T_INNER = 50.0          # fixed inner threshold (continental/transitional boundary), km
  T_OUTER_STEPS = [200, 175, 150, 125, 100, 75, 50, 25]  # 8 threshold test points
  SUB_PROXIMITY_THRESHOLD_KM = 200.0   # near-subduction classification threshold
  INTERVAL_BINS = {"interval_1": [4, 5], "interval_2": [15], "interval_3": [21]}
  ```

---

## 2. Subduction zone proximity computation

```python
def parse_sub_boundaries(steps_path: Path) -> np.ndarray:
```
- Read `PB2002_steps.dat` line by line
- Split each line on whitespace; skip blank lines
- Strip `:` prefix from the last field (boundary type) and first two fields (row index, plate pair)
- Filter to rows where stripped boundary type is `"SUB"` or `"OCB"` (subduction + oceanic convergent)
- Extract columns: `lon1=col[2]`, `lat1=col[3]`, `lon2=col[4]`, `lat2=col[5]` (float)
- For each segment, generate three sample points: (lon1, lat1), (lon2, lat2), midpoint ((lon1+lon2)/2, (lat1+lat2)/2)
- Return as `np.ndarray` of shape `(N, 2)` with columns [lon, lat]; log N (expected ~850 for SUB+OCB)

```python
def latlon_to_unit_sphere(lon_deg: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
```
- Convert degrees to radians; return `(N, 3)` array of unit-sphere Cartesian coordinates:
  - `x = cos(lat) * cos(lon)`, `y = cos(lat) * sin(lon)`, `z = sin(lat)`

```python
def compute_subduction_distances(
    event_lons: np.ndarray,
    event_lats: np.ndarray,
    sub_boundary_lonlat: np.ndarray,
) -> np.ndarray:
```
- Convert both event and boundary points to unit-sphere coordinates
- Build `scipy.spatial.cKDTree` from boundary unit-sphere points
- For each event, query k=5 nearest neighbors (to ensure the true nearest is found despite Cartesian approximation near the antimeridian)
- For each event, compute geodesic distance to each of the 5 candidate neighbors using `pyproj.Geod(ellps="WGS84").inv()`:
  - `_, _, dist_m = geod.inv(event_lon, event_lat, cand_lon, cand_lat)`
  - `dist_km = dist_m / 1000.0`
- Take the minimum as `dist_to_subduction_km` for that event
- Return `np.ndarray` of shape `(n_events,)` — one distance per event; log mean, median, and fraction within 200 km

Add `dist_to_subduction_km` as a column to the main merged dataframe. Classify:
```python
df["near_subduction"] = df["dist_to_subduction_km"] <= SUB_PROXIMITY_THRESHOLD_KM
```

---

## 3. Threshold sweep analysis

```python
def classify_at_threshold(dist_series: pd.Series, t_outer: float, t_inner: float = 50.0) -> pd.Series:
```
- Returns per-event class label:
  - `"continental"` if `dist <= t_inner`
  - `"oceanic"` if `dist > t_outer`
  - `"transitional"` if `t_inner < dist <= t_outer` (empty when `t_outer <= t_inner`)

```python
def compute_chi2_stats(phases: np.ndarray, k: int = 24) -> dict:
```
- `n = len(phases)`; if n < k: return null stats dict with `n` recorded
- `observed = np.bincount(np.floor(phases * k).astype(int) % k, minlength=k)`
- `chi2_stat, p_chi2 = scipy.stats.chisquare(observed, np.full(k, n / k))`
- `cramers_v = np.sqrt(chi2_stat / (n * (k - 1)))`
- Per A1b interval z-scores (same as A3.B1):
  - `z_i = (observed[bins].sum() - len(bins) * n / k) / np.sqrt(len(bins) * n / k)`
- Return: `{"n": int, "chi2_k24": float, "p_chi2_k24": float, "cramers_v": float, "bin_counts": list, "interval_1_z": float, "interval_2_z": float, "interval_3_z": float}`

For each `t_outer` in `T_OUTER_STEPS`:
1. Classify all events via `classify_at_threshold(df["dist_km_gshhg"], t_outer)` → `class_col`
2. For each class label in `["oceanic", "transitional", "continental"]`:
   - Filter subset; compute stats via `compute_chi2_stats`
3. Also compute subduction cross-tabulation for this threshold:
   - For each class: `n_near_sub = (class_events["near_subduction"]).sum()`; `pct_near_sub = n_near_sub / n_class`
4. Record step dict:
```json
{
  "t_outer_km": float,
  "t_inner_km": 50.0,
  "oceanic":     {"n": int, "chi2_k24": float, "p_chi2_k24": float, "cramers_v": float,
                  "interval_1_z": float, "interval_2_z": float, "interval_3_z": float,
                  "n_near_sub": int, "pct_near_sub": float},
  "transitional": {same},
  "continental":  {same}
}
```
5. Flag `"oceanic_significant"` (p_chi2_k24 < 0.05) and `"transitional_significant"` for each step

Identify:
- `first_oceanic_significant_threshold`: the largest T at which oceanic p < 0.05 (working inward from T=200)
- `signal_migration_step`: whether the oceanic signal appears as T decreases (oceanic significance increases with smaller T as more transitional events are reclassified as oceanic)

---

## 4. GCMT validation of subduction proximity proxy

```python
def validate_subduction_proxy(df: pd.DataFrame) -> dict:
```
Using events with `match_confidence == "proximity"` (GCMT-matched):
- `pct_thrust_near_sub = (df.loc[near_sub_mask & matched_mask, "mechanism"] == "thrust").mean()`
- `pct_thrust_far_sub = (df.loc[~near_sub_mask & matched_mask, "mechanism"] == "thrust").mean()`
- `thrust_enrichment = pct_thrust_near_sub / max(pct_thrust_far_sub, 0.001)`
- Flag `proxy_validated = (pct_thrust_near_sub > pct_thrust_far_sub)` — SUB boundary events should preferentially show thrust mechanism

Record under `"gcmt_validation"` in results JSON:
```json
{
  "gcmt_validation": {
    "n_matched": int,
    "n_near_sub_matched": int,
    "pct_thrust_near_sub": float,
    "pct_thrust_far_sub": float,
    "thrust_enrichment_ratio": float,
    "proxy_validated": bool
  }
}
```

---

## 5. Results JSON structure

```json
{
  "case": "A3.B3",
  "title": "Ocean/Coast Sequential Threshold Sensitivity",
  "parameters": {
    "n_catalog": 9210,
    "k_bins": 24,
    "t_inner_km": 50.0,
    "t_outer_steps": [200, 175, 150, 125, 100, 75, 50, 25],
    "sub_proximity_threshold_km": 200.0,
    "julian_year_secs": 31557600.0
  },
  "subduction_proximity": {
    "n_sub_boundary_points": int,
    "n_near_subduction": int,
    "pct_near_subduction": float,
    "mean_dist_to_sub_km": float,
    "median_dist_to_sub_km": float
  },
  "gcmt_validation": {...},
  "threshold_sweep": [
    /* list of 8 step dicts, one per T_OUTER value */
  ],
  "summary": {
    "baseline_t200_oceanic_p": float,
    "first_oceanic_significant_threshold_km": float,
    "signal_migration_observed": bool,
    "transitional_pct_near_sub_at_baseline": float
  }
}
```

---

## 6. Visualizations

In `src/visualization-case-a3-b3.py`:

Import: `matplotlib`, `cartopy`, `cartopy.crs`, `cartopy.feature`, `numpy`, `json`, `pathlib`, `logging`

---

**Figure 1 — Threshold sweep trajectory** (`output/case-a3-b3-threshold-sweep.png`):
- 3-row stacked subplot sharing x-axis (T_outer in km, x-axis reversed so 200 is on left and 25 on right — moving from B2 baseline toward coastline)
- Row 1: chi-square p-value (log scale); three lines: oceanic=steelblue, transitional=green, continental=red; horizontal dashed line at p=0.05 (gray)
- Row 2: Cramér's V (linear); same three colored lines
- Row 3: class n (linear); stacked area or three lines; shows how class sizes change as threshold sweeps
- Vertical dashed line at T=200 (B2 baseline) labeled "A2.B2 baseline"
- Mark oceanic significance threshold crossing (if any) with a vertical dotted line
- x-axis: tick at each threshold step (200, 175, ..., 25), reversed; label "Outer threshold (km from coast)"
- Title: "Threshold Sensitivity — GSHHG Classification"
- 300 DPI, publication quality

---

**Figure 2 — Global map with subduction boundaries** (`output/case-a3-b3-global-map.png`):
- Cartopy `Robinson` projection, global extent
- Plot all 9,210 ISC-GEM events colored by GSHHG baseline classification (T=200):
  - oceanic: steelblue, alpha=0.4
  - transitional: orange, alpha=0.4
  - continental: red, alpha=0.4
  - Point size proportional to magnitude: `size = (usgs_mag - 5.5) ** 2 * 1.5`
- Overlay SUB-type PB2002 boundary segments from `lib/PB2002_steps.dat` as black lines (linewidth=1.0), drawn as individual segments from (lon1, lat1) to (lon2, lat2)
- Add `cartopy.feature.COASTLINE` (resolution='110m', linewidth=0.5, color='gray')
- Legend: event classes with n counts; note "Black lines = PB2002 subduction boundaries"
- Title: "ISC-GEM Events — GSHHG Classification at T=200 km (A2.B2 baseline)"
- 300 DPI

---

**Figure 3 — Subduction proximity cross-tabulation** (`output/case-a3-b3-subduction-crosstab.png`):
- 2-row stacked subplot
- Row 1: Stacked bar chart; x-axis = T_outer threshold steps; y-axis = n per class; bars split by `near_subduction` (dark shade = near subduction, light shade = far from subduction); one bar group per threshold per class, or grouped differently for readability
- Row 2: Line plot; x-axis = T_outer; y-axis = fraction of oceanic events that are near subduction (`pct_near_sub` for oceanic class); steelblue line; horizontal dashed line at the fraction for the full catalog
- x-axis reversed (200 → 25, matching Figure 1 orientation)
- Title: "Subduction Zone Proximity vs. Threshold Reclassification"
- 300 DPI

---

**Figure 4 — Bin distributions at baseline and key threshold** (`output/case-a3-b3-binplots.png`):
- 2 rows × 3 columns: rows = baseline (T=200) and key threshold (the smallest T where oceanic p < 0.05, or T=100 if no threshold crosses); columns = oceanic, transitional, continental
- Each panel: horizontal steelblue bar chart at k=24; dashed expected-count line; 1-SD threshold; gray shaded bands for A1b intervals
- Annotate: n, χ², p, Cramér's V
- Row labels: "T=200 km (A2.B2 baseline)", "T=[key threshold] km"
- 300 DPI

---

**Figure 5 — Regional detail maps** (`output/case-a3-b3-region-maps.png`):
- 1 row × 3 columns; each panel is a Cartopy regional inset using `PlateCarree` projection
- Event classification at GSHHG baseline (T=200): oceanic=steelblue, transitional=orange, continental=red; alpha=0.5; point size proportional to magnitude as in Figure 2
- Overlay SUB-type PB2002 boundary segments as black lines (linewidth=1.2)
- Add `cartopy.feature.COASTLINE` at resolution `'50m'` (higher detail than global map), linewidth=0.6, color='dimgray'
- Add `cartopy.feature.LAND` at resolution `'50m'`, facecolor='whitesmoke', zorder=0 (land as background)
- **Panel A — Japan**: extent lon [128, 148], lat [30, 46]; rationale: dense seismicity with well-defined Pacific and Philippine plate subduction arcs; complex island coastline tests the transitional zone boundary precisely
- **Panel B — Philippines**: extent lon [114, 130], lat [4, 22]; rationale: multiple active subduction zones (Manila Trench, Philippine Trench, Cotabato Trench); island-arc geometry makes the coastal/ocean boundary definition non-trivial
- **Panel C — New Zealand**: extent lon [164, 180], lat [-48, -32]; rationale: Hikurangi subduction zone (NE coast) transitions to Alpine strike-slip fault (SW); provides a contrasting fault regime within a single region to assess whether classification correlates with tectonic type
- Each panel annotated with: region name, n_events in extent, n_transitional and n_oceanic counts
- Shared legend below panels: event class colors and n; note "Black lines = PB2002 subduction boundaries"
- Title: "Regional Coastline Resolution — GSHHG Classification at T=200 km (A2.B2 baseline)"
- 300 DPI, figsize=(18, 7)

---

## 7. Test suite

In `tests/test-case-a3-b3.py`:

- `test_catalog_load`: assert n=9210; assert `event_at` parses without NaT; assert phase in [0.0, 1.0)
- `test_classification_joins`: assert GSHHG and PB2002 joins produce 9210 rows each with no NaN in `dist_to_coast_km`
- `test_baseline_class_counts`: assert GSHHG at T=200: oceanic≈1952, transitional≈3459, continental≈3799 (±5 tolerance for merge handling)
- `test_threshold_step_count`: assert threshold sweep produces exactly 8 step records
- `test_class_partition_all_steps`: for each threshold step, assert oceanic + transitional + continental == 9210
- `test_n_monotonic_oceanic`: as T_outer decreases from 200 to 25, assert oceanic n is non-decreasing (more events reclassified as oceanic at tighter thresholds)
- `test_chi2_bounds`: for all threshold steps and all classes, assert chi2_k24 ≥ 0 and p_chi2_k24 in [0.0, 1.0]
- `test_cramers_v_bounds`: assert cramers_v ≥ 0 for all
- `test_baseline_matches_b2`: assert T=200 oceanic p_chi2_k24 is within 0.01 of 0.061 (A2.B2 reference); assert transitional p_chi2_k24 < 0.01
- `test_sub_boundary_point_count`: assert n_sub_boundary_points is between 700 and 1100 (expected ~850 for SUB+OCB at 3 points per segment)
- `test_dist_to_sub_nonneg`: assert all dist_to_subduction_km values ≥ 0
- `test_dist_to_sub_plausible`: assert mean dist_to_subduction_km is between 100 and 2000 km (plausible global range)
- `test_gcmt_validation_structure`: load results JSON; assert `"gcmt_validation"` present with `"proxy_validated"` bool key
- `test_results_json_completeness`: assert `"threshold_sweep"` list has length 8; assert each entry contains all three class keys with `chi2_k24` and `p_chi2_k24`
- `test_sub_parse_type_filter`: confirm that parsing `PB2002_steps.dat` and filtering type `SUB` produces no rows with type `OTF`, `OSR`, `CTF`, `CRB` (spot-check that filter is working)
- `test_region_maps_output_exists`: assert `output/case-a3-b3-region-maps.png` exists and file size > 100 KB (confirms all three panels rendered)

---

## 8. Whitepaper

In `output/case-a3-b3-whitepaper.md`:

Standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Code [model name]).

### Sections:

1. **Abstract** (150–200 words): state the question (does the solar-phase signal migrate from the transitional zone into the oceanic class as the boundary is tightened?), summarize the threshold sweep design and subduction proximity extension, state the key outcome, and state the implication for the mechanism discrimination question.

2. **Data Source**: describe raw catalog (n=9,210); GSHHG classification file with pre-computed `dist_to_coast_km`; PB2002 classification file (secondary); GCMT focal join (52.9% match, used for proxy validation only); `PB2002_steps.dat` and its SUB-type boundary segments.

3. **Methodology**
   - 3.1 Phase-normalized binning: Julian year constant; cite `data-handling.md`
   - 3.2 Baseline classification (A2.B2 reference): GSHHG three-class system at T_inner=50 km, T_outer=200 km; state baseline counts
   - 3.3 Threshold sweep design: 8 test points from T=200 to T=25 km in 25 km steps; T_inner held fixed at 50 km; note that at T_outer ≤ T_inner the transitional zone collapses
   - 3.4 Chi-square and Cramér's V per class per threshold step
   - 3.5 Subduction zone proximity: PB2002_steps.dat SUB+OCB type parsing; 3-point-per-segment sampling; unit-sphere cKDTree nearest-neighbor + pyproj WGS84 geodesic refinement; 200 km proximity threshold
   - 3.6 GCMT mechanism validation of subduction proxy: thrust enrichment ratio in near-vs-far subduction events
   - 3.7 Cross-tabulation: ocean class × subduction proximity at each threshold step

4. **Results**
   - 4.1 Threshold sweep trajectory: embed Figure 1; describe at which threshold (if any) oceanic p crosses 0.05; describe Cramér's V trajectory for oceanic vs. transitional; state whether signal migrates or remains confined to transitional class
   - 4.2 Subduction proximity: embed Figure 3; report fraction of transitional and oceanic events near SUB boundaries at baseline; describe whether reclassified events (as T tightens) are disproportionately near subduction zones
   - 4.3 Global map: embed Figure 2; describe geographic distribution of classes at baseline; note alignment of transitional zone events with subduction arcs visible in the SUB boundary overlay
   - 4.4 Bin distributions: embed Figure 4; compare oceanic/transitional/continental bin structure at baseline vs. key threshold; note interval-level changes
   - 4.5 GCMT validation: report thrust enrichment ratio; state whether subduction proximity proxy is validated

5. **Cross-Topic Comparison**: compare to A2.B2 (baseline result and PB2002 flip); compare to A2.B3 (tectonic regime — thrust events should preferentially appear in transitional zone if subduction drives the signal); compare to A2.B4 (depth stratification — mid-crustal band dominance; subduction zone events are mid-crustal)

6. **Interpretation**: state whether the transitional zone's signal is driven by subduction zone geometry or coastal proximity in general; assess whether the signal is broad enough to be considered "oceanic" at any tested threshold; guard against both confirming and dismissing the hydrological vs. geometric framing based on this single case.

7. **Limitations**: the 200 km subduction proximity threshold is arbitrary; PB2002 subduction segments do not capture all active subduction zones (some minor arcs may be absent); cKDTree nearest-neighbor uses 3-point-per-segment sampling which undersamples long boundary arcs; GCMT validation covers only 52.9% of events; the inner threshold (50 km, continental/transitional boundary) is not swept in this case.

8. **References**:
   - Yeager, J. (2026). A2.B2: Ocean vs. Continent Location — Hydrological Loading Discrimination. erebus-vee-two internal report.
   - Yeager, J. (2026). A2.B3: Tectonic Regime Stratification. erebus-vee-two internal report.
   - Yeager, J. (2026). A2.B4: Depth Stratification — Surface Loading Penetration Test. erebus-vee-two internal report.
   - Bird, P. (2003). An updated digital model of plate boundaries. *Geochemistry, Geophysics, Geosystems*, 4(3).
   - Scholz, J. R. et al. (2019). Volcanic and tidal stresses triggering earthquakes. *Nature Communications*.

---

## 9. Update context docs

After all outputs are generated and tests pass:

- Append to `topic-a3/docs/topic-summary.md`:
  ```
  ## Case A3.B3: Ocean/Coast Sequential Threshold Sensitivity
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [threshold at which oceanic p crosses 0.05 (or "never" if it doesn't); Cramér's V for oceanic at baseline vs. tightest threshold; pct of transitional events near subduction at baseline; GCMT proxy validated T/F; signal migration observed T/F]
  ```
- Update Case B3 status from `Planning` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a3/CLAUDE.md` Case Table
