# Case A3.B5: Corrected Null-Distribution Geometric Variable Test

**Status:** Planning

**Source case:** A2.B5
**Informed by:** A3.B2 (phase alignment), A3.B3 (tectonic class), A3.B4 (depth localization)

**Intent statement:**
Case A2.B5 tested three solar geometric variables — `solar_declination`, `declination_rate`, `earth_sun_distance` — against the base cyclic variable `solar_secs`, comparing chi-square statistics to rank which variable best explains the observed seismic clustering. However, A2.B5 applied a uniform expected-count null to non-cyclic variables, which is incorrect: the Sun does not spend equal time at each value of declination, declination rate, or Earth-Sun distance. For example, solar declination changes slowly near solstices and rapidly near equinoxes, so the Earth spends more time at high-magnitude declination values than at near-zero values. A chi-square test against a uniform null will therefore register spurious non-uniformity for any dataset, regardless of the physical mechanism.

A3.B5 corrects this by computing **time-weighted expected distributions** for each non-cyclic variable using a dense analytic synthetic time series (daily resolution, 1950–2021, ~26,000 points). The corrected expected counts replace the uniform null in the chi-square test, yielding a statistic that isolates genuine excess clustering beyond what the variable's own temporal distribution would predict.

Additionally, A3.B5 extends A2.B5 by stratifying the variable comparison by tectonic class (continental, from A3.B3) and depth band (mid-crustal 20–70 km, from A3.B4). The signal concentrates in these subpopulations; if a geometric variable drives the mechanism, it should exhibit its strongest corrected chi-square in the stratum where the signal is most robust.

**Relationship to prior cases:**
A2.B5 established the four-variable comparison framework and reported uncorrected chi-square statistics. A3.B5 supersedes A2.B5 as the definitive geometric variable test by correcting the null. A3.B2 found that the NH peaks near the March equinox (phase ≈0.23) and SH near mid-August (phase ≈0.65), both equinox-adjacent — the declination rate peaks at equinoxes, making it mechanistically motivated as the candidate variable. A3.B3 and A3.B4 provide the tectonic and depth stratification targets.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| Solar geometry catalog | `data/iscgem/celestial-geometry/solar_geometry_global.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth`, `solar_declination`, `declination_rate`, `earth_sun_distance` |
| GSHHG classification | `data/iscgem/plate-location/ocean_class_gshhg_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |

**Confirmed variable ranges (from solar_geometry_global.csv):**
- `solar_declination`: [−23.4452, +23.4459] degrees
- `declination_rate`: [−0.3898, +0.3957] degrees/day
- `earth_sun_distance`: [0.9832, 1.0168] AU

**Phase normalization:** `solar_phase = (solar_secs / 31_557_600.0) % 1.0` (Julian year constant; uniform null is correct for this cyclic variable).

**Tectonic class boundaries (GSHHG, matching A3.B3 baseline):**
- `continental`: `dist_to_coast_km ≤ 50`
- `transitional`: `50 < dist_to_coast_km ≤ 200`
- `oceanic`: `dist_to_coast_km > 200`

**Mid-crustal depth band (from A3.B4):** `20 ≤ depth < 70` km.

**Adaptive k rule:** k=24 if n≥500; k=16 if 200≤n<500; k=12 if 100≤n<200; flag n<100.

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a3/`
- All output paths: `BASE_DIR / "output" / ...`
- All data paths: `BASE_DIR.parent / "data" / "iscgem" / ...`

**Planned outputs:**
- `src/case-a3-b5-analysis.py` — null generation, corrected chi-square for all variables × strata, writes results JSON
- `src/visualization-case-a3-b5.py` — all PNG figures
- `tests/test-case-a3-b5.py` — test suite (all tests must pass)
- `output/case-a3-b5-results.json`
- `output/case-a3-b5-whitepaper.md`
- `output/case-a3-b5-null-distributions.png` — Figure 1: time-weighted vs. uniform null for all three non-cyclic variables
- `output/case-a3-b5-binplots.png` — Figure 2: corrected observed vs. expected distributions for all four variables
- `output/case-a3-b5-variable-ranking.png` — Figure 3: corrected Cramér's V by variable × stratum
- `output/case-a3-b5-correction-delta.png` — Figure 4: uncorrected vs. corrected chi-square comparison
- `output/case-a3-b5-a1b-alignment.png` — Figure 5: A1b interval physical-value alignment

---

## 1. Environment and data loading

In `src/case-a3-b5-analysis.py`:

- Imports: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define paths:
  ```python
  SOLAR_PATH = BASE_DIR.parent / "data" / "iscgem" / "celestial-geometry" / "solar_geometry_global.csv"
  GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
  ```
- Load solar geometry catalog; assert n=9210; log column names
- Assert required columns: `solar_secs`, `solar_declination`, `declination_rate`, `earth_sun_distance`, `depth`, `latitude`
- Log actual min/max for each variable; compare to confirmed ranges; warn if outside bounds
- Load GSHHG classification; assert n=9210; merge onto solar catalog by `usgs_id`
- Assert no NaN in `dist_to_coast_km` after merge
- Compute `solar_phase = (solar_secs / 31_557_600.0) % 1.0`
- Define constants:
  ```python
  JULIAN_YEAR_SECS = 31_557_600.0
  MIDCRUSTAL_MIN_KM = 20.0
  MIDCRUSTAL_MAX_KM = 70.0
  VARIABLES = ["solar_phase", "solar_declination", "declination_rate", "earth_sun_distance"]
  STRATA = {
      "full": None,                       # no filter
      "continental": "ocean_class == 'continental'",
      "midcrustal": "20 <= depth < 70",
      "continental_midcrustal": "ocean_class == 'continental' and 20 <= depth < 70",
  }
  ```

---

## 2. Analytic null distribution generation

```python
def generate_analytic_null(year_start: int, year_end: int, k: int) -> dict:
    """
    Generate time-weighted expected bin fractions for solar_declination,
    declination_rate, and earth_sun_distance using a dense daily analytic model.

    Parameters
    ----------
    year_start : int
        First year of catalog coverage (inclusive).
    year_end : int
        Last year of catalog coverage (inclusive).
    k : int
        Number of bins.

    Returns
    -------
    dict with keys "solar_declination", "declination_rate", "earth_sun_distance",
    each mapping to a np.ndarray of shape (k,) representing the expected fraction
    of time in each bin (sums to 1.0).
    """
```

Implementation:
- Generate daily time steps: `days = np.arange(0, (year_end - year_start + 1) * 365.25)` — approximately 26,300 steps
- Compute J2000 offset D: `D = days + (year_start - 2000) * 365.25 - 0.5`
- Solar mean longitude: `L = np.radians((280.46 + 0.9856474 * D) % 360)`
- Mean anomaly: `g = np.radians((357.528 + 0.9856003 * D) % 360)`
- Ecliptic longitude (first-order correction): `lam = L + np.radians(1.915 * np.sin(g) + 0.020 * np.sin(2 * g))`
- Declination: `dec = np.degrees(np.arcsin(np.sin(np.radians(23.439)) * np.sin(lam)))`
- Earth-Sun distance: `dist = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)`
- Declination rate: finite difference: `rate = np.gradient(dec, days)` — units: degrees/day

For each variable, bin the synthetic series into k equal-width bins over the variable's observed range (use the actual data min/max as bin edges to match the per-event binning exactly):
  ```python
  edges = np.linspace(data_min, data_max, k + 1)
  null_counts, _ = np.histogram(synthetic_values, bins=edges)
  null_fractions = null_counts / null_counts.sum()
  ```

Return `{"solar_declination": null_fractions_dec, "declination_rate": null_fractions_rate, "earth_sun_distance": null_fractions_dist}`.

Also record the synthetic series stats for QC:
- `n_synthetic_points`: total synthetic time steps
- `dec_synthetic_range`: [min, max] of synthetic declination
- `dist_synthetic_range`: [min, max] of synthetic distance

The `solar_phase` variable uses a uniform null (each bin expected fraction = 1/k) — no synthetic series needed.

---

## 3. Per-event binning

```python
def assign_bins(df: pd.DataFrame, k: int, var_ranges: dict) -> pd.DataFrame:
    """
    Assign bin indices for all four variables using equal-width bins over observed ranges.

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with solar_phase, solar_declination, declination_rate, earth_sun_distance.
    k : int
        Number of bins.
    var_ranges : dict
        Maps variable name to (min, max) tuple from actual data.

    Returns
    -------
    DataFrame with added columns bin_solar_phase, bin_solar_declination,
    bin_declination_rate, bin_earth_sun_distance.
    """
```

- `bin_solar_phase = np.floor(df["solar_phase"] * k).astype(int).clip(0, k - 1)`
- For each non-cyclic variable `v`: `bin_v = np.floor((df[v] - vmin) / (vmax - vmin) * k).astype(int).clip(0, k - 1)`
- Record `var_ranges` under results JSON parameters for traceability

---

## 4. Corrected chi-square statistics

```python
def compute_corrected_chi2(
    bin_col: np.ndarray,
    k: int,
    null_fractions: np.ndarray,
    n: int,
) -> dict:
    """
    Chi-square test using time-weighted expected counts.

    Parameters
    ----------
    bin_col : np.ndarray
        Integer bin assignments for each event.
    k : int
        Number of bins.
    null_fractions : np.ndarray
        Expected fraction of time in each bin (shape k, sums to 1.0).
    n : int
        Total events in the subset.

    Returns
    -------
    dict with: n, k, chi2_corrected, p_corrected, cramers_v_corrected,
               chi2_uniform, p_uniform, cramers_v_uniform,
               bin_counts, expected_corrected, expected_uniform.
    """
```

- `observed = np.bincount(bin_col, minlength=k)`
- `expected_corrected = null_fractions * n`
- `expected_uniform = np.full(k, n / k)`
- `chi2_c, p_c = scipy.stats.chisquare(observed, expected_corrected)`
- `chi2_u, p_u = scipy.stats.chisquare(observed, expected_uniform)`
- `cramers_v_c = np.sqrt(chi2_c / (n * (k - 1)))`
- `cramers_v_u = np.sqrt(chi2_u / (n * (k - 1)))`
- Return both corrected and uniform statistics for direct comparison

For `solar_phase` (cyclic, uniform null): `null_fractions = np.full(k, 1.0 / k)`; corrected and uniform are identical — record once.

---

## 5. Full analysis loop

For each stratum in `STRATA`:
1. Filter catalog to stratum subset; apply adaptive k; if n < 100 skip and record null entry
2. Assign bins at adaptive k
3. Generate corrected null at adaptive k: `generate_analytic_null(1950, 2021, k=adaptive_k)`
4. For each variable, call `compute_corrected_chi2` with appropriate null
5. Record all four variable results under stratum key

Also record for the full catalog at k=24, k=16, k=32 (for multi-k comparison, matching A2.B5 structure).

---

## 6. Variable ranking

For the full catalog at k=24:
- Rank four variables by `cramers_v_corrected` (descending)
- Record `most_significant_variable_corrected`: variable with lowest `p_corrected`
- Record `most_significant_variable_uniform`: variable with lowest `p_uniform` (for comparison with A2.B5 result)
- For each stratum, record which variable has the highest `cramers_v_corrected`

A1b interval physical-value alignment (at k=24):
- Interval 1 center phase 0.22 → DOY ≈ 80 → compute expected variable values:
  - `solar_declination` ≈ 0.0° (March equinox crossing)
  - `declination_rate` ≈ +0.39 deg/day (maximum positive rate)
  - `earth_sun_distance` ≈ 0.9990 AU (near mean, approaching aphelion)
- Interval 2 center phase 0.64 → DOY ≈ 234 → compute expected variable values:
  - `solar_declination` ≈ +12.0° (post-summer-solstice, declining)
  - `declination_rate` ≈ −0.26 deg/day (declining toward autumn equinox)
  - `earth_sun_distance` ≈ 1.0143 AU (near aphelion)
- Interval 3 center phase 0.90 → DOY ≈ 329 → compute expected variable values:
  - `solar_declination` ≈ −21.0° (approaching winter solstice)
  - `declination_rate` ≈ −0.12 deg/day (near solstice minimum)
  - `earth_sun_distance` ≈ 0.9877 AU (approaching perihelion)
- For each interval, identify which bin each expected value falls in; check whether that bin is elevated (observed > corrected_expected + sqrt(corrected_expected)) in the full-catalog corrected distribution

Store under `"variable_ranking"` and `"a1b_alignment"` in results JSON.

---

## 7. Results JSON structure

```json
{
  "case": "A3.B5",
  "title": "Corrected Null-Distribution Geometric Variable Test",
  "parameters": {
    "n_catalog": 9210,
    "julian_year_secs": 31557600.0,
    "midcrustal_depth_range_km": [20.0, 70.0],
    "tectonic_class_boundaries_km": {"continental_max": 50, "transitional_max": 200},
    "analytic_null_year_range": [1950, 2021],
    "var_ranges": {
      "solar_declination": [float, float],
      "declination_rate": [float, float],
      "earth_sun_distance": [float, float]
    }
  },
  "null_generation": {
    "n_synthetic_points": int,
    "dec_synthetic_range": [float, float],
    "dist_synthetic_range": [float, float],
    "note": "Analytic solar model: mean longitude + first-order aberration + Kepler distance; daily resolution"
  },
  "strata": {
    "full":                    {"n": int, "k": int, "solar_phase": {...}, "solar_declination": {...}, "declination_rate": {...}, "earth_sun_distance": {...}},
    "continental":             {"n": int, "k": int, "solar_phase": {...}, ...},
    "midcrustal":              {"n": int, "k": int, "solar_phase": {...}, ...},
    "continental_midcrustal":  {"n": int, "k": int, "solar_phase": {...}, ...}
  },
  "full_multibin": {
    "k16": {"solar_phase": {...}, "solar_declination": {...}, "declination_rate": {...}, "earth_sun_distance": {...}},
    "k24": {/* same */},
    "k32": {/* same */}
  },
  "variable_ranking": {
    "by_cramers_v_corrected": [{"variable": str, "cramers_v_corrected": float, "p_corrected": float, "cramers_v_uniform": float, "p_uniform": float}],
    "most_significant_variable_corrected": str,
    "most_significant_variable_uniform": str,
    "per_stratum_top_variable": {"full": str, "continental": str, "midcrustal": str, "continental_midcrustal": str}
  },
  "a1b_alignment": {
    "interval_1": {"phase_center": 0.22, "doy_approx": 80, "solar_declination_expected": float, "declination_rate_expected": float, "earth_sun_distance_expected": float, "elevated_variables": [str]},
    "interval_2": {"phase_center": 0.64, "doy_approx": 234, ...},
    "interval_3": {"phase_center": 0.90, "doy_approx": 329, ...}
  }
}
```

Each variable stats dict within a stratum:
```json
{
  "n": int, "k": int, "low_n": bool,
  "chi2_corrected": float, "p_corrected": float, "cramers_v_corrected": float,
  "chi2_uniform": float,  "p_uniform": float,  "cramers_v_uniform": float,
  "bin_counts": [int, ...],
  "expected_corrected": [float, ...],
  "peak_bin": int, "peak_phase_fraction": float
}
```

---

## 8. Visualizations

In `src/visualization-case-a3-b5.py`:

Import: `matplotlib`, `matplotlib.pyplot`, `numpy`, `json`, `pathlib`, `logging`

---

**Figure 1 — Null distribution correction** (`output/case-a3-b5-null-distributions.png`):
- 1-row × 3-column grid; one panel per non-cyclic variable (solar_declination, declination_rate, earth_sun_distance)
- Each panel: x-axis = variable value (physical units); y-axis = fraction of time / probability density
  - Steelblue bars: time-weighted corrected null distribution (from analytic synthetic series)
  - Gray dashed line: uniform null (1/k for each bin)
  - Difference shaded: area where corrected null deviates from uniform (green = corrected > uniform, red = corrected < uniform)
- x-axis ticks at physical values (degrees, deg/day, AU); label bin boundaries
- Annotate: variable name, total synthetic points, max deviation from uniform (%)
- Title: "Time-Weighted Null Distributions — Corrected vs. Uniform (A3.B5)"
- 300 DPI, figsize=(14, 5)

---

**Figure 2 — Corrected bin distributions** (`output/case-a3-b5-binplots.png`):
- 2-row × 2-column grid: solar_phase (top-left), solar_declination (top-right), declination_rate (bottom-left), earth_sun_distance (bottom-right); all at k=24 full catalog
- Each panel: steelblue bars = observed; solid line = corrected expected; dashed line = uniform expected; 1-SD band around corrected expected as light fill
- For solar_phase panel: add A1b interval gray shading (bins 4–5, 15, 21) and equinox/solstice markers
- For non-cyclic panels: x-axis labeled with physical units; mark A1b interval expected values with vertical dotted red lines
- Annotate each panel: χ²_corrected, p_corrected, Cramér's V_corrected; and χ²_uniform, p_uniform in smaller text for comparison
- Title: "Solar Variable Distributions vs. Corrected Null (Full Catalog, k=24)"
- 300 DPI, figsize=(14, 10)

---

**Figure 3 — Variable ranking by stratum** (`output/case-a3-b5-variable-ranking.png`):
- 2-row × 2-column grid: one panel per stratum (full, continental, mid-crustal, continental × mid-crustal)
- Each panel: grouped horizontal bar chart; y-axis = four variable names; x-axis = Cramér's V_corrected
  - Steelblue bars = corrected; light gray bars = uniform (paired for each variable)
  - Significance markers: ** if p_corrected < 0.01, * if p_corrected < 0.05, ns otherwise
- Panel subtitle: stratum name and n
- Title: "Corrected Cramér's V by Variable and Stratum (A3.B5)"
- 300 DPI, figsize=(13, 9)

---

**Figure 4 — Correction impact** (`output/case-a3-b5-correction-delta.png`):
- Scatter plot; one point per (variable × stratum) combination; x-axis = χ²_uniform, y-axis = χ²_corrected
- Identity line (y=x) as gray diagonal; points above line = correction increases chi-square; below = decreases
- Color by variable: solar_declination=steelblue, declination_rate=darkorange, earth_sun_distance=green, solar_phase=red (should lie exactly on y=x)
- Marker style by stratum: full=circle, continental=square, midcrustal=triangle, continental_midcrustal=diamond
- Annotate any point where |χ²_corrected − χ²_uniform| / χ²_uniform > 0.20 with its label
- Title: "Chi-Square Before vs. After Null Correction (A3.B5)"
- 300 DPI, figsize=(8, 7)

---

**Figure 5 — A1b interval alignment** (`output/case-a3-b5-a1b-alignment.png`):
- 3-row × 3-column grid: rows = A1b intervals 1, 2, 3; columns = solar_declination, declination_rate, earth_sun_distance
- Each cell: steelblue bars = corrected-null-normalized observed distribution (observed / expected_corrected, so 1.0 = no excess); vertical red dotted line at the expected physical value for that A1b interval's phase center; horizontal dashed line at 1.0 (no excess reference)
- Row labels: "Interval 1 (Phase ≈0.22, ~Mar 22)", "Interval 2 (Phase ≈0.64, ~Aug 22)", "Interval 3 (Phase ≈0.90, ~Nov 24)"
- Column labels: variable name with units
- Annotate each cell with whether the interval's expected-value bin is elevated (obs/exp_corrected > 1 + 1/√exp_corrected)
- Title: "A1b Interval Alignment with Solar Variables — Corrected Null (A3.B5)"
- 300 DPI, figsize=(13, 10)

---

## 9. Test suite

In `tests/test-case-a3-b5.py`:

- `test_solar_geo_load`: assert n=9210; assert all required columns present; assert no NaN in key columns
- `test_column_ranges`: assert `solar_declination` range within [−24, +24]; `declination_rate` within [−0.45, +0.45]; `earth_sun_distance` within [0.975, 1.025]
- `test_phase_range`: assert all `solar_phase` values in [0.0, 1.0)
- `test_gshhg_merge`: assert no NaN in `dist_to_coast_km` after merge
- `test_null_generation_shape`: assert `generate_analytic_null(1950, 2021, k=24)` returns dict with three keys each mapping to array of shape (24,); each sums to 1.0 within floating-point tolerance
- `test_null_generation_nonuniform`: assert that at least one bin in the corrected null for `solar_declination` deviates from 1/24 by more than 5% (confirms the correction is non-trivial)
- `test_null_generation_k_parametric`: assert `generate_analytic_null(1950, 2021, k=16)` returns arrays of shape (16,)
- `test_bins_in_range`: assert all bin assignments for all four variables are in [0, k-1] for k=24
- `test_corrected_chi2_solar_phase_equals_uniform`: assert that for `solar_phase`, `chi2_corrected == chi2_uniform` (within 1e-6 tolerance) in the full stratum — the cyclic variable uses the same null
- `test_all_strata_present`: load results JSON; assert `strata` key contains all four stratum keys with non-null `solar_phase` entries
- `test_variable_ranking_complete`: assert `variable_ranking → by_cramers_v_corrected` has exactly 4 entries covering all four variable names
- `test_most_significant_valid`: assert `most_significant_variable_corrected` is one of the four variable names
- `test_a1b_alignment_keys`: assert `a1b_alignment` has keys `interval_1`, `interval_2`, `interval_3` each with `elevated_variables` list
- `test_solar_phase_full_significant`: assert full stratum `solar_phase` `p_corrected < 0.05` (regression anchor — the primary signal must survive) — log a warning but do not fail if p≥0.05
- `test_output_figures_exist`: assert all 5 PNG files exist and size > 50 KB

---

## 10. Whitepaper

In `output/case-a3-b5-whitepaper.md`:

Standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Code [model name]).

### Sections:

1. **Abstract** (150–200 words): state that A2.B5 applied an incorrect uniform null to non-cyclic variables; describe the correction (time-weighted analytic null); state the four variables tested and four strata; report the ranking result and whether any geometric variable outperforms `solar_secs` after correction; state the A1b interval alignment finding.

2. **Data Source**: solar geometry catalog (n=9,210; columns: `solar_declination`, `declination_rate`, `earth_sun_distance`, `solar_secs`; computed via Skyfield 1.54 + JPL DE421 in an external pipeline). GSHHG tectonic classification for stratification. State confirmed variable ranges.

3. **Methodology**
   - 3.1 Phase normalization: Julian year constant for `solar_secs`; uniform null appropriate for cyclic variable
   - 3.2 Non-cyclic binning: equal-width bins over observed range; formula for bin assignment; note that bin edges are anchored to observed data min/max to ensure consistency between per-event binning and null generation
   - 3.3 Analytic null distribution: describe the synthetic time series approach (daily resolution, 1950–2021, ~26,000 points); solar mean longitude + first-order aberration for declination; Kepler distance formula; finite-difference declination rate; note that the analytic model has ~0.3° accuracy for declination, sufficient for k=24 null estimation
   - 3.4 Corrected chi-square: describe how corrected expected counts replace uniform expected counts; record both corrected and uniform statistics for comparison with A2.B5
   - 3.5 Stratification: four strata (full, continental, mid-crustal, continental × mid-crustal); adaptive k rule; cite A3.B3 and A3.B4 for stratum definitions

4. **Results**
   - 4.1 Null corrections: embed Figure 1; describe how much the corrected null deviates from uniform for each variable; identify which variable has the largest correction magnitude
   - 4.2 Corrected distributions: embed Figure 2; report χ²_corrected, p_corrected, Cramér's V_corrected for all four variables at k=24 in the full catalog; compare to χ²_uniform for each to quantify correction impact
   - 4.3 Variable ranking by stratum: embed Figure 3; state which variable has the highest corrected Cramér's V in each stratum; note whether `declination_rate` outperforms `solar_secs` in any stratum
   - 4.4 Correction impact: embed Figure 4; describe which variables are most affected by the null correction; note any cases where the correction reverses a significance call (significant under uniform, non-significant under corrected, or vice versa)
   - 4.5 A1b interval alignment: embed Figure 5; for each A1b interval, state whether the corresponding physical variable bin is elevated in the corrected distribution; summarize which intervals have single-variable explanations and which do not

5. **Cross-Topic Comparison**
   - **Solar Declination Rate-of-Change vs. Position Test (A2.B5):** A2.B5 applied a uniform null to all variables. A3.B5 corrects this; any change in ranking relative to A2.B5 is attributable to the null correction rather than data differences.
   - **Hemisphere Stratification Refinement (A3.B2):** A3.B2 found that NH populations peak near the March equinox (phase ≈0.23) and SH near mid-August (phase ≈0.65), both equinox-adjacent. The declination rate peaks at equinoxes; if `declination_rate` achieves higher corrected Cramér's V than `solar_secs`, this provides variable-level support for an equinox-timing mechanism.
   - **Ocean/Coast Sequential Threshold Sensitivity (A3.B3):** A3.B3 found the continental class most robustly significant. The continental stratum in A3.B5 tests whether the variable preference differs in the most signal-stable population.
   - **Depth × Magnitude Stratification with Moho Isolation (A3.B4):** A3.B4 localized the signal to the mid-crustal band. If the geometric variable signal is also most concentrated there, that reinforces the depth-specific mechanism interpretation.

6. **Interpretation**: state objectively which variable (if any) best explains the seismic clustering after null correction; note that `solar_secs` represents the annual cycle as a whole and may outperform any single geometric decomposition if the mechanism is multi-component; guard against concluding a null result across all geometric variables as definitive, given sample-size limitations in the stratified subsets.

7. **Limitations**: the analytic null model has ~0.3° accuracy for solar declination (sufficient for k=24, marginal for k=32); declination rate is computed by finite difference of the analytic model, introducing a second-order approximation; the null distribution is computed over the full 1950–2021 span uniformly, whereas the catalog has non-uniform temporal coverage in its earlier years; declustered catalogs are not used, so aftershock contamination is present; the continental × mid-crustal stratum may have marginal sample size at k=24.

8. **References**
   - Yeager, J. (2026). A2.B5: Solar Declination Rate-of-Change vs. Position Test. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B2: Hemisphere Stratification Refinement. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B3: Ocean/Coast Sequential Threshold Sensitivity. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation. erebus-vee-two internal report.

---

## 11. Update context docs

After all outputs are generated and tests pass:

- Append to `topic-a3/docs/topic-summary.md`:
  ```
  ## Case A3.B5: Corrected Null-Distribution Geometric Variable Test
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [corrected Cramér's V ranking of four variables at k=24 full catalog; most significant variable under corrected null; whether correction changed ranking vs. uniform null; which variable dominates in continental and mid-crustal strata; A1b interval alignment summary]
  ```
- Update Case B5 status from `Planning` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a3/CLAUDE.md` Case Table
