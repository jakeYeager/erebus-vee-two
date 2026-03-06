# Case A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation

**Status:** Complete

**Source case:** A2.B4
**Reference case:** A3.B3 (geospatial implementation)

**Intent statement:**
A2.B4 established that only the mid-crustal band (20–70 km, n=4,561, χ²=85.48, p=4.02×10⁻⁹, V=0.0285) carries a significant solar-phase signal. All other depth bands (shallow 0–20 km, intermediate 70–300 km, deep >300 km) are non-significant. The console summary flagged a magnitude–depth confound: larger events tend to be deeper, so the mid-crustal dominance may partly reflect the magnitude trend from A2.A3. A3.B3 introduced a further confound: the mid-crustal band is the depth range of slab interface seismicity (thrust events cluster at 20–50 km), and A3.B3 found that 65.8% of transitional-zone events are within 200 km of a subduction boundary (GCMT thrust enrichment ratio 1.97). The A2.B4 mid-crustal signal may therefore be a proxy for subduction geometry rather than crustal loading depth.

A3.B4 is a three-component case designed to disentangle these confounds:

1. **Two-way stratification sub-test:** Partition the 9,210 events into a depth-band × magnitude-band matrix and compute chi-square independently in each cell, holding one variable constant while varying the other. This tests whether the mid-crustal signal persists within the M6.0–6.9 majority (ruling out the magnitude confound) and whether signal appears in any other depth band under magnitude control.

2. **Moho isolation sub-test:** Using CRUST1.0 (Laske et al. 2013, 1°×1° global crustal model) to assign per-event Moho depth, test a narrow proximity band around the local crust-mantle boundary (`|focal_depth − moho_depth_km| ≤ delta`) to determine whether the signal is specifically concentrated at the Moho transition. Three delta values are tested: 5, 10, and 15 km. This replaces the fixed-depth-proxy approach with spatially adaptive Moho assignment.

3. **Subduction proximity cross-tabulation sub-test:** Within the significant depth×magnitude cells (expected: mid-crustal × M6.0–6.9 and mid-crustal × M7.0–7.9), cross-tabulate events by subduction proximity (reusing the B3 PB2002 cKDTree pipeline) and compare chi-square statistics between near-subduction and far-subduction subsets. This tests whether the mid-crustal signal survives geographical partitioning consistent with its subduction-geometry alternative explanation.

**Relationship to prior cases:**
A2.B4 provided the depth stratification baseline. A2.A3 established the magnitude–signal relationship. A3.B3 introduced and validated the subduction proximity pipeline now reused here. A3.B4 feeds directly into A3.C1 (subduction zone subset test), which will further refine the geographic partitioning.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |
| GSHHG classification | `data/iscgem/plate-location/ocean_class_gshhg_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |
| GCMT focal mechanism join | `data/iscgem/focal-mechanism/focal_join_global.csv` | 9,210 | `usgs_id`, `mechanism`, `rake`, `match_confidence` |
| PB2002 boundary steps | `lib/PB2002_steps.dat` | 5,819 rows | row_idx, plate_pair, lon1, lat1, lon2, lat2, ..., boundary_type |
| CRUST1.0 boundaries | `lib/crust1.bnds` | 64,800 rows (180×360) | 9 layer boundary depths (km, negative = below sea level); col 8 = Moho |

**A2.B4 baseline results at k=24 (regression anchors):**
- shallow 0–20 km: n=3,063, χ²=24.66, p=0.368, V=0.0187
- mid-crustal 20–70 km: n=4,561, χ²=85.48, p=4.02×10⁻⁹, V=0.0285
- intermediate 70–300 km: n=1,030, χ²=17.05, p=0.807, V=0.0268
- deep >300 km: n=556, χ²=20.26, p=0.626, V=0.0398

**Projected depth×magnitude cell sizes** (approximate, based on ISC-GEM distribution ~73%/22%/5% across M6-7/7-8/8+):

| Depth band | M6.0–6.9 (~73%) | M7.0–7.9 (~22%) | M8.0+ (~5%) |
|-----------|----------------|----------------|------------|
| shallow (3,063) | ~2,236 | ~674 | ~153 |
| mid-crustal (4,561) | ~3,329 | ~1,003 | ~228 |
| intermediate (1,030) | ~752 | ~227 | ~52 |
| deep (556) | ~406 | ~122 | ~28 |

**Adaptive k rule:** k=24 if n≥500; k=16 if 200≤n<500; k=12 if 100≤n<200; flag as low-n (report only) if n<100.

**A1b baseline intervals (k=24 bin mapping):**
- Interval 1: bins 4–5 (phase [0.1667, 0.2500)) — March equinox region
- Interval 2: bin 15 (phase [0.6250, 0.6667)) — mid-August region
- Interval 3: bin 21 (phase [0.8750, 0.9167)) — late-November region

Phase normalization: `phase = (solar_secs / 31_557_600.0) % 1.0` (Julian year constant).

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a3/`
- `LIB_DIR = BASE_DIR.parent / "lib"`
- All output paths: `BASE_DIR / "output" / ...`
- All data paths: `BASE_DIR.parent / "data" / "iscgem" / ...`

**Planned outputs:**
- `src/case-a3-b4-analysis.py` — CRUST1.0 load + Moho assignment, subduction proximity, two-way stratification, Moho isolation, subduction cross-tabulation, writes results JSON
- `src/visualization-case-a3-b4.py` — generates all PNG figures
- `tests/test-case-a3-b4.py` — test suite (all tests must pass)
- `output/case-a3-b4-results.json` — full stratification matrix, Moho isolation stats, subduction cross-tabulation
- `output/case-a3-b4-whitepaper.md` — methodology, results, interpretation
- `output/case-a3-b4-stratification-matrix.png` — heatmap of p-values across depth×magnitude cells
- `output/case-a3-b4-cramer-trends.png` — Cramér's V per depth band, one line per magnitude band
- `output/case-a3-b4-moho-isolation.png` — bin distributions and significance metrics at delta=5, 10, 15 km
- `output/case-a3-b4-subduction-crosstab.png` — near-sub vs far-sub chi-square within mid-crustal magnitude cells
- `output/case-a3-b4-moho-map.png` — Cartopy global map: CRUST1.0 Moho depth as background contour + events colored by Moho-proximal status

---

## 1. Environment and data loading

In `src/case-a3-b4-analysis.py`:

- Imports: `pandas`, `numpy`, `scipy.stats`, `scipy.spatial`, `scipy.interpolate`, `pyproj`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`; `LIB_DIR = BASE_DIR.parent / "lib"`
- Define paths:
  ```python
  RAW_PATH     = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
  GSHHG_PATH   = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
  FOCAL_PATH   = BASE_DIR.parent / "data" / "iscgem" / "focal-mechanism" / "focal_join_global.csv"
  STEPS_PATH   = LIB_DIR / "PB2002_steps.dat"
  CRUST1_PATH  = LIB_DIR / "crust1.bnds"
  ```
- Load raw catalog; parse `event_at` as UTC; assert n=9210
- Load GSHHG and focal join; merge onto raw by `usgs_id`; assert no NaN in `dist_to_coast_km` after merge
- Compute `phase = (solar_secs / 31_557_600.0) % 1.0`
- Check NaN in `depth` column: log count; exclude from all depth analyses; record `n_depth_null`
- Define constants:
  ```python
  K_BINS = 24
  MAG_BANDS = [
      {"label": "m6_6.9", "min": 6.0, "max": 7.0},
      {"label": "m7_7.9", "min": 7.0, "max": 8.0},
      {"label": "m8_plus", "min": 8.0, "max": 99.0},
  ]
  DEPTH_BANDS = [
      {"label": "shallow_0-20km",       "min": 0,   "max": 20},
      {"label": "midcrustal_20-70km",   "min": 20,  "max": 70},
      {"label": "intermediate_70-300km","min": 70,  "max": 300},
      {"label": "deep_300km+",          "min": 300, "max": 9999},
  ]
  MOHO_DELTAS_KM = [5.0, 10.0, 15.0]
  SUB_PROXIMITY_THRESHOLD_KM = 200.0
  INTERVAL_BINS = {"interval_1": [4, 5], "interval_2": [15], "interval_3": [21]}
  ```
- Adaptive k function:
  ```python
  def select_k(n: int) -> int | None:
      if n >= 500: return 24
      if n >= 200: return 16
      if n >= 100: return 12
      return None  # low-n; flag but still compute at k=12 for record
  ```

---

## 2. CRUST1.0 Moho depth assignment

```python
def load_crust1_moho(crust1_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
```
- Load `crust1.bnds` via `np.loadtxt(crust1_path)` → shape (64800, 9)
- Reshape to (180, 360, 9): `bnds = raw.reshape(180, 360, 9)`
- Extract Moho layer: `moho_grid = -bnds[:, :, 8]` (negate: values are negative depths, result is positive km below sea level)
- Define grid axes:
  - `lats = np.arange(89.5, -90.5, -1.0)` — 89.5°N to -89.5°N (180 values, cell centers)
  - `lons = np.arange(-179.5, 180.5, 1.0)` — -179.5°E to 179.5°E (360 values, cell centers)
- Build bilinear interpolation surface:
  ```python
  from scipy.interpolate import RectBivariateSpline
  moho_interp = RectBivariateSpline(lats[::-1], lons, moho_grid[::-1], kx=1, ky=1)
  ```
  (flip lat axis so input is ascending for RectBivariateSpline)
- Return `(moho_grid, lats, lons)` and the `moho_interp` object
- Log summary: min, max, mean Moho depth across grid; expected range ~5–80 km

```python
def assign_moho_depth(
    event_lats: np.ndarray,
    event_lons: np.ndarray,
    moho_interp: RectBivariateSpline,
) -> np.ndarray:
```
- For each event, evaluate `moho_interp(lat, lon)` as a scalar; clip to [3.0, 90.0] km to handle edge cells
- Return `np.ndarray` of shape `(n_events,)` — per-event Moho depth in km
- Log mean, median, and std of assigned Moho depths

Add `moho_depth_km` column to main DataFrame. Compute Moho-proximity flags for each delta:
```python
for delta in MOHO_DELTAS_KM:
    col = f"moho_proximal_{int(delta)}km"
    df[col] = np.abs(df["depth"] - df["moho_depth_km"]) <= delta
```
Log n_proximal for each delta.

**Attribution:** CRUST1.0 boundary data sourced from https://igppweb.ucsd.edu/~gabi/crust1.html (Laske et al. 2013).

---

## 3. Subduction zone proximity computation

Reuse the B3 pipeline exactly:

```python
def parse_sub_boundaries(steps_path: Path) -> np.ndarray:
```
- Read `PB2002_steps.dat` line by line; split on whitespace; skip blank lines
- Strip `:` prefix from last field (boundary type); filter to type `"SUB"` or `"OCB"`
- Extract `lon1=col[2]`, `lat1=col[3]`, `lon2=col[4]`, `lat2=col[5]` (float)
- Generate 3 sample points per segment: both endpoints + midpoint
- Return `np.ndarray` shape `(N, 2)` [lon, lat]; log N

```python
def latlon_to_unit_sphere(lon_deg: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
```
- Convert to radians; return `(N, 3)` unit-sphere Cartesian: `x=cos(lat)cos(lon)`, `y=cos(lat)sin(lon)`, `z=sin(lat)`

```python
def compute_subduction_distances(
    event_lons: np.ndarray,
    event_lats: np.ndarray,
    sub_boundary_lonlat: np.ndarray,
) -> np.ndarray:
```
- Build `scipy.spatial.cKDTree` from boundary unit-sphere points
- Query k=5 nearest neighbors per event; refine with `pyproj.Geod(ellps="WGS84").inv()`
- Return minimum geodesic distance per event in km

Add `dist_to_subduction_km` and `near_subduction` columns to main DataFrame.

---

## 4. Chi-square statistics

```python
def compute_chi2_stats(phases: np.ndarray, k: int) -> dict:
```
- `n = len(phases)`; `observed = np.bincount(np.floor(phases * k).astype(int) % k, minlength=k)`
- `chi2_stat, p_chi2 = scipy.stats.chisquare(observed, np.full(k, n / k))`
- `cramers_v = np.sqrt(chi2_stat / (n * (k - 1)))`
- Per A1b interval z-scores (at k=24 only; skip for k<24): `z_i = (observed[bins].sum() - len(bins) * n / k) / np.sqrt(len(bins) * n / k)`
- Return: `{"n": int, "k": int, "low_n": bool, "chi2": float, "p_chi2": float, "cramers_v": float, "bin_counts": list, "interval_1_z": float|None, "interval_2_z": float|None, "interval_3_z": float|None}`

---

## 5. Two-way stratification analysis

For each `depth_band` in `DEPTH_BANDS`:
  For each `mag_band` in `MAG_BANDS`:
  1. Subset: `cell_df = df[(df["depth"] >= d_min) & (df["depth"] < d_max) & (df["usgs_mag"] >= m_min) & (df["usgs_mag"] < m_max)]`
  2. `k = select_k(len(cell_df))`; set `low_n = (k is None)`; if `low_n`, use k=12 for record but flag
  3. Compute stats via `compute_chi2_stats(cell_df["phase"].values, k or 12)`
  4. Record cell dict (see JSON structure below)

Also compute depth-band totals (all magnitudes combined) to serve as regression anchor vs A2.B4.

---

## 6. Moho isolation analysis

For each `delta` in `MOHO_DELTAS_KM`:
  1. Subset: `moho_df = df[df[f"moho_proximal_{int(delta)}km"]]`; log n
  2. `k = select_k(len(moho_df))`
  3. Compute overall stats via `compute_chi2_stats`
  4. Split by `ocean_class`:
     - `continental_moho = moho_df[moho_df["ocean_class"] == "continental"]` — continental Moho (~25–35 km)
     - `oceanic_moho = moho_df[moho_df["ocean_class"] == "oceanic"]` — oceanic Moho (~7–12 km)
     - `transitional_moho = moho_df[moho_df["ocean_class"] == "transitional"]` — ambiguous, report separately
  5. Compute stats per split; record all under `moho_isolation`

---

## 7. Subduction proximity cross-tabulation

For the mid-crustal band only, for each magnitude band:
1. Split into `near_sub` and `far_sub` subsets using `near_subduction` flag
2. Compute chi-square stats for each subset independently
3. Record near-sub and far-sub stats side by side

Also compute for the full mid-crustal band (all magnitudes) for comparison to A3.B3's baseline.

---

## 8. Results JSON structure

```json
{
  "case": "A3.B4",
  "title": "Depth × Magnitude Two-Way Stratification with Moho Isolation",
  "parameters": {
    "n_catalog": 9210,
    "n_depth_null": int,
    "k_bins_primary": 24,
    "adaptive_k_thresholds": {"k24": 500, "k16": 200, "k12": 100},
    "moho_deltas_km": [5.0, 10.0, 15.0],
    "sub_proximity_threshold_km": 200.0,
    "julian_year_secs": 31557600.0,
    "crust1_source": "https://igppweb.ucsd.edu/~gabi/crust1.html"
  },
  "moho_assignment": {
    "n_events_assigned": int,
    "mean_moho_depth_km": float,
    "median_moho_depth_km": float,
    "std_moho_depth_km": float,
    "n_proximal_by_delta": {"5": int, "10": int, "15": int}
  },
  "subduction_proximity": {
    "n_sub_boundary_points": int,
    "n_near_subduction": int,
    "pct_near_subduction": float,
    "mean_dist_to_sub_km": float
  },
  "depth_band_totals": {
    "shallow_0-20km":        {"n": int, "k": 24, "chi2": float, "p_chi2": float, "cramers_v": float},
    "midcrustal_20-70km":    {"n": int, "k": 24, "chi2": float, "p_chi2": float, "cramers_v": float},
    "intermediate_70-300km": {"n": int, "k": 24, "chi2": float, "p_chi2": float, "cramers_v": float},
    "deep_300km+":           {"n": int, "k": 24, "chi2": float, "p_chi2": float, "cramers_v": float}
  },
  "stratification_matrix": {
    "shallow_0-20km": {
      "m6_6.9":   {"n": int, "k": int, "low_n": bool, "chi2": float, "p_chi2": float, "cramers_v": float, "interval_1_z": float|null, "interval_2_z": float|null, "interval_3_z": float|null},
      "m7_7.9":   {same},
      "m8_plus":  {same}
    },
    "midcrustal_20-70km":    {same 3 mag bands},
    "intermediate_70-300km": {same 3 mag bands},
    "deep_300km+":           {same 3 mag bands}
  },
  "moho_isolation": {
    "delta_5km":  {"n_total": int, "k": int, "chi2": float, "p_chi2": float, "cramers_v": float,
                   "continental": {same keys}, "oceanic": {same keys}, "transitional": {same keys}},
    "delta_10km": {same},
    "delta_15km": {same}
  },
  "subduction_crosstab_midcrustal": {
    "all_magnitudes": {
      "near_sub": {"n": int, "k": int, "chi2": float, "p_chi2": float, "cramers_v": float},
      "far_sub":  {same}
    },
    "m6_6.9":  {"near_sub": {same}, "far_sub": {same}},
    "m7_7.9":  {"near_sub": {same}, "far_sub": {same}},
    "m8_plus": {"near_sub": {same}, "far_sub": {same}}
  }
}
```

---

## 9. Visualizations

In `src/visualization-case-a3-b4.py`:

Import: `matplotlib`, `cartopy`, `cartopy.crs`, `cartopy.feature`, `numpy`, `json`, `pathlib`, `logging`

---

**Figure 1 — Two-way stratification matrix** (`output/case-a3-b4-stratification-matrix.png`):
- 4 rows × 3 columns; rows = depth bands (shallow → deep), columns = magnitude bands (M6–7, M7–8, M8+)
- Each cell: color-coded by p-value using a white (p=1.0) → red (p=0) scale; cells with p < 0.05 additionally marked with a bold border
- Annotate each cell with: n (top line), p-value (middle line, 2 sig figs scientific notation), k used (bottom line, gray), "LOW-N" watermark if n < 100
- Row label on left: depth band name; column label at top: magnitude band name
- Add a reference row at the right edge showing A2.B4 full-band p-values (all magnitudes combined) for regression comparison
- Colorbar: "Chi-square p-value"; title: "Solar Phase Chi-Square P-Value by Depth × Magnitude"
- 300 DPI, figsize=(12, 10)

---

**Figure 2 — Cramér's V trends** (`output/case-a3-b4-cramer-trends.png`):
- Single panel; x-axis = depth band (ordinal 1–4, shallow to deep); y-axis = Cramér's V
- One line per magnitude band: M6–7 (steelblue solid), M7–8 (steelblue dashed), M8+ (steelblue dotted)
- Additional line for depth-band totals (all magnitudes): red solid
- Filled circle markers for p < 0.05; open circle markers for p ≥ 0.05; 'x' marker for low-n cells
- Horizontal dashed line at the full-catalog global V (reference)
- Annotate each marker with n count (small gray text)
- x-tick labels: depth band names; title: "Cramér's V by Depth Band and Magnitude Band"
- 300 DPI

---

**Figure 3 — Moho isolation bin distributions** (`output/case-a3-b4-moho-isolation.png`):
- 3 columns (one per delta: 5 km, 10 km, 15 km) × 2 rows (top: all tectonic settings; bottom: continental vs. oceanic side-by-side bar chart)
- Top row each panel: horizontal steelblue bar chart at k=12 (or k=16/24 per adaptive rule); dashed expected-count line; annotate n, χ², p, Cramér's V; gray shaded bands for A1b intervals
- Bottom row each panel: grouped bar chart comparing continental (red) vs. oceanic (steelblue) Moho-proximal bin counts at k=12; legend; annotate n per group
- Column titles: "Δ=5 km", "Δ=10 km", "Δ=15 km"
- Note in each panel: continental Moho assigned via CRUST1.0 (1°×1°)
- 300 DPI, figsize=(18, 8)

---

**Figure 4 — Subduction proximity cross-tabulation** (`output/case-a3-b4-subduction-crosstab.png`):
- 2-row stacked subplot; x-axis = magnitude band (M6–7, M7–8, M8+, All) within mid-crustal depth band
- Row 1: grouped bar chart of chi-square p-value for near-subduction (orange) vs. far-subduction (steelblue) subsets; y-axis log scale; horizontal dashed line at p=0.05; annotate n above each bar
- Row 2: grouped bar chart of Cramér's V; same grouping and colors; horizontal dashed line at full-catalog V (reference from A3.B3)
- Title: "Mid-Crustal Solar Signal: Near-Subduction vs. Far-Subduction (PB2002 SUB+OCB, ≤200 km)"
- 300 DPI

---

**Figure 5 — Moho depth global map with event overlay** (`output/case-a3-b4-moho-map.png`):
- Cartopy `Robinson` projection, global extent
- Background: CRUST1.0 Moho depth grid as a filled contour or pcolormesh using white (shallow Moho, ~5 km oceanic) → blue (deep Moho, ~60 km+ continental) color scale; plot every 4th row/column of the 1°×1° grid for performance
- Overlay all 9,210 ISC-GEM events at `delta=10 km` classification:
  - Moho-proximal events: orange filled circles, size proportional to magnitude (`size = (mag - 5.5)**2 * 2.0`), alpha=0.7
  - Non-proximal events: gray filled circles, size proportional to magnitude, alpha=0.2
- Add `cartopy.feature.COASTLINE` (resolution='110m', linewidth=0.4, color='black')
- Legend: "Moho-proximal (Δ≤10 km)" and "Other events"; colorbar for Moho depth labeled "CRUST1.0 Moho depth (km)"
- Title: "CRUST1.0 Moho Depth and Moho-Proximal Events (Δ≤10 km)"
- 300 DPI

---

## 10. Test suite

In `tests/test-case-a3-b4.py`:

- `test_catalog_load`: assert n=9210; assert phase in [0.0, 1.0); assert `depth` has no unexpected negatives (log count of depth < 0)
- `test_depth_band_partition`: assert sum of all band sizes + n_depth_null == 9210; assert no overlap between bands
- `test_magnitude_band_partition`: for each depth band, assert sum across 3 magnitude bands equals depth band total (no events lost or double-counted)
- `test_adaptive_k`: for each cell, assert k assignment matches rule: k=24 if n≥500, k=16 if 200≤n<500, k=12 otherwise
- `test_chi2_bounds`: for all cells in stratification_matrix, assert chi2 ≥ 0 and p_chi2 ∈ [0.0, 1.0]
- `test_cramer_v_bounds`: for all cells, assert cramers_v ≥ 0
- `test_crust1_load`: load `crust1.bnds` with `np.loadtxt`; assert shape after reshape is (180, 360, 9); assert Moho values (layer 8) are all negative (depth below sea level reference); assert range is within [-90, 0] km
- `test_moho_assignment_range`: assert all `moho_depth_km` values are in [3.0, 90.0] km (plausible global range after negation and clipping)
- `test_moho_proximal_counts`: assert n_proximal at delta=10 > 0 and < 9210; assert n_proximal increases monotonically with delta (5 < 10 < 15 km)
- `test_sub_boundary_points`: assert n_sub_boundary_points between 700 and 1100 (matching B3 test)
- `test_sub_distances_nonneg`: assert all dist_to_subduction_km ≥ 0
- `test_a2b4_regression`: load results JSON; assert depth_band_totals["midcrustal_20-70km"]["chi2"] is within ±2.0 of 85.48 (A2.B4 reference); assert depth_band_totals["midcrustal_20-70km"]["p_chi2"] < 1e-6
- `test_midcrustal_m6_significant`: assert stratification_matrix["midcrustal_20-70km"]["m6_6.9"]["p_chi2"] < 0.05 (mid-crustal signal expected to persist in M6–7 majority; if this fails, it is a scientifically significant null result — test should be flagged, not silenced)
- `test_results_json_completeness`: assert all four depth bands present in stratification_matrix; assert all three mag bands present per depth band; assert moho_isolation has keys "delta_5km", "delta_10km", "delta_15km"; assert subduction_crosstab_midcrustal present with expected keys

---

## 11. Whitepaper

In `output/case-a3-b4-whitepaper.md`:

Standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Code [model name]).

### Sections:

1. **Abstract** (150–200 words): state the three confounds (magnitude, Moho/depth, subduction geometry); describe the three-component design; state the key outcome on whether the mid-crustal signal survives magnitude control; note the CRUST1.0 Moho isolation result; state the subduction proximity finding and its implication for mechanism discrimination.

2. **Data Source**: raw catalog (n=9,210); GSHHG classification for tectonic setting (ocean_class); GCMT focal join (52.9% match); PB2002_steps.dat for subduction proximity; CRUST1.0 `crust1.bnds` (1°×1°, Laske et al. 2013) for per-event Moho depth — note attribution to https://igppweb.ucsd.edu/~gabi/crust1.html.

3. **Methodology**
   - 3.1 Phase-normalized binning: Julian year constant; cite `data-handling.md`
   - 3.2 Depth bands: same four definitions as A2.B4 (0–20, 20–70, 70–300, >300 km); note depth null count
   - 3.3 Magnitude bands: M6.0–6.9, M7.0–7.9, M8.0+; rationale (three roughly even-logarithmic subdivisions of the catalog range)
   - 3.4 Adaptive k selection: state thresholds (k=24/16/12 by n); explain why variable k is necessary for sparse cells
   - 3.5 CRUST1.0 Moho depth assignment: describe `crust1.bnds` format (64,800 cells, 9 layers); Moho = layer 8 (negative km); bilinear interpolation via `scipy.interpolate.RectBivariateSpline`; attribution to Laske et al. 2013 and https://igppweb.ucsd.edu/~gabi/crust1.html
   - 3.6 Moho isolation proximity filter: `|focal_depth − moho_depth_km| ≤ delta`; three delta values tested (5, 10, 15 km); continental vs. oceanic Moho split using GSHHG `ocean_class`
   - 3.7 Subduction proximity pipeline: PB2002_steps.dat SUB+OCB filter; 3-point-per-segment sampling; unit-sphere cKDTree + pyproj WGS84 geodesic refinement; 200 km threshold (same as A3.B3)

4. **Results**
   - 4.1 Two-way stratification matrix: embed Figure 1; report which cells are significant; focus on mid-crustal column — does signal persist in M6–7 subset independently?; describe whether any other depth band becomes significant under magnitude stratification
   - 4.2 Cramér's V trends: embed Figure 2; describe trajectory across depth bands for each magnitude band; compare to A2.B4 depth-band totals; note whether mid-crustal peak is specific to a particular magnitude subset
   - 4.3 Moho isolation: embed Figure 3; report n_proximal at each delta; report chi-square and Cramér's V for Moho-proximal events at delta=10 km (primary); compare to the mid-crustal baseline; distinguish continental vs. oceanic Moho proximal results
   - 4.4 Subduction proximity cross-tabulation: embed Figure 4; report near-sub vs. far-sub chi-square and V in mid-crustal band across magnitude bands; state whether signal persists in far-subduction mid-crustal events
   - 4.5 Global Moho map: embed Figure 5; describe geographic distribution of Moho-proximal events; note overlap with subduction arcs and A3.B3 transitional zone

5. **Cross-Topic Comparison**:
   - **Depth Stratification — Surface Loading Penetration Test (A2.B4):** Regression to A2.B4 mid-crustal χ²=85.48 confirms reproducibility. The depth-band total comparison tests whether the three-component stratification changes the headline result.
   - **Magnitude Stratification of the Solar Signal (A2.A3):** A2.A3 established the magnitude–signal relationship in the full catalog. If mid-crustal signal persists in M6–7 alone, it is not the magnitude effect from (A2.A3) expressed through depth.
   - **Ocean/Coast Sequential Threshold Sensitivity (A3.B3):** (A3.B3) showed 65.8% of transitional-zone events are subduction-proximal with GCMT thrust enrichment ratio 1.97. A3.B4's subduction cross-tabulation tests the depth dimension of this same mechanism.

6. **Interpretation**: state which confound (magnitude, Moho proximity, subduction geometry) is most supported as the driver of mid-crustal dominance; note implications for mechanism discrimination (surface loading, crustal resonance, or subduction geometry); maintain objectivity — the three sub-tests are designed to distinguish, not confirm, a preferred interpretation.

7. **Limitations**:
   - CRUST1.0 1°×1° resolution may inadequately resolve Moho variability in regions with rapid crustal thickness changes (e.g., subduction zones, Tibetan Plateau margins); ISC-GEM focal depths have varying precision, particularly for pre-1970 events
   - Magnitude bands are defined by integer thresholds; the ~5% M8.0+ population produces low-n cells in all depth bands except mid-crustal
   - The subduction proximity threshold (200 km) and Moho delta values (5–15 km) are arbitrary; sensitivity to these choices is not systematically explored in A3.B4
   - A2.B4's non-monotonic Spearman classification ("increasing with depth" rho=0.80) reflects only four band-level V values with one outlier; this classification is not reproduced or extended in A3.B4
   - GCMT mechanism data covers only 52.9% of events; the GSHHG ocean_class proxy for Moho type (continental/oceanic) has a ~50 km uncertainty in the transitional zone

8. **References**:
   - Yeager, J. (2026). A2.B4: Depth Stratification — Surface Loading Penetration Test. erebus-vee-two internal report.
   - Yeager, J. (2026). A2.A3: Magnitude Stratification of the Solar Signal. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B3: Ocean/Coast Sequential Threshold Sensitivity. erebus-vee-two internal report.
   - Laske, G., Masters, G., Ma, Z., and Pasyanos, M. (2013). Update on CRUST1.0 — A 1-degree Global Model of Earth's Crust. *Geophysical Research Abstracts*, 15, Abstract EGU2013-2658. Data: https://igppweb.ucsd.edu/~gabi/crust1.html
   - Bird, P. (2003). An updated digital model of plate boundaries. *Geochemistry, Geophysics, Geosystems*, 4(3).
   - Zhan, Z. and Shearer, P.M. (2015). Possible seasonality in large deep-focus earthquakes. *Geophysical Research Letters*, 42(18), 7366–7373.

---

## 12. Update context docs

After all outputs are generated and tests pass:

- Append to `topic-a3/docs/topic-summary.md`:
  ```
  ## Case A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [whether mid-crustal signal survives M6–7 magnitude control (p<0.05 T/F); V range across significant cells; Moho isolation result at delta=10 km (p, n); near-sub vs far-sub mid-crustal p comparison; CRUST1.0 mean Moho depth assigned]
  ```
- Update Case B4 status from `Planning` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a3/CLAUDE.md` Case Table
