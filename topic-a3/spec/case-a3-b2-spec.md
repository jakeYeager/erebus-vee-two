# Case A3.B2: Hemisphere Stratification Refinement

**Status:** Complete

**Source case:** A2.B1
**Informed by:** A3.B3 (tectonic composition findings), A3.B4 (depth signal localization)

**Intent statement:**
Case A3.B2 tests whether the NH/SH solar-phase signal asymmetry observed in A2.B1 is a tectonic-composition artifact or a genuine hemispheric phase difference. A3.B3 established that the solar-phase signal concentrates in continental and transitional (subduction-proximal) tectonic classes; A3.B4 localized the signal to the mid-crustal band (20–70 km). Because the NH holds a disproportionate share of continental landmass and active subduction systems relative to the SH, a raw hemisphere split partially recovers this tectonic imbalance rather than isolating a solar-forcing difference.

A3.B2 addresses this with four sub-tests: (1) tectonic-matched hemisphere comparison — NH/SH χ² within each tectonic class; (2) mid-crustal hemisphere split — hemisphere comparison locked to the signal-bearing depth band; (3) phase alignment comparison — whether significant NH and SH cells peak at the same solar phase; (4) revised Interval 1 SH threshold sensitivity — does the A2.B1 Interval 1 SH absence persist when restricted to continental-only SH events, or does it dissolve when tectonic composition is controlled? Declustering (G-K, Reasenberg) is incorporated as a sensitivity layer on sub-test (1) rather than as standalone primary tests.

**Relationship to prior cases:**
A2.B1 established the hemisphere split framework, three-interval structure (from Adhoc A1b), and four symmetry tests. A3.B1 confirmed non-stationarity across rolling windows. A3.B3 established the tectonic classification (GSHHG baseline T_outer=200 km) and confirmed transitional+continental significance; A3.B4 established mid-crustal dominance and magnitude-independence. A3.B2 uses the GSHHG tectonic classification from B3 and the mid-crustal depth band from B4 as primary stratification dimensions.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |
| GSHHG classification | `data/iscgem/plate-location/ocean_class_gshhg_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |
| G-K mainshocks | `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv` | 5,883 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |
| Reasenberg mainshocks | `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv` | 8,265 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |

**Phase normalization:** `phase = (solar_secs / 31_557_600.0) % 1.0` (Julian year constant).

**Hemisphere split:** `NH = latitude > 0`; `SH = latitude < 0`; `equatorial = latitude == 0` (excluded, typically n≈0).

**Tectonic class boundaries (GSHHG baseline, matching A3.B3 T_outer=200 km):**
- `continental`: `dist_to_coast_km ≤ 50`
- `transitional`: `50 < dist_to_coast_km ≤ 200`
- `oceanic`: `dist_to_coast_km > 200`

**Mid-crustal depth band (from A3.B4):** `20 ≤ depth < 70` km.

**Adaptive k rule (matching A3.B4):** k=24 if n≥500; k=16 if 200≤n<500; k=12 if 100≤n<200; flag n<100 and skip χ² (record null stats).

**A1b baseline intervals at k=24:**
- Interval 1: bins 4–5 (phase [0.1667, 0.2500)) — March equinox
- Interval 2: bin 15 (phase [0.6250, 0.6667)) — mid-August
- Interval 3: bin 21 (phase [0.8750, 0.9167)) — late-November

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a3/`
- All output paths: `BASE_DIR / "output" / ...`
- All data paths: `BASE_DIR.parent / "data" / "iscgem" / ...`

**Planned outputs:**
- `src/case-a3-b2-analysis.py` — all four sub-tests; writes results JSON
- `src/visualization-case-a3-b2.py` — all PNG figures
- `tests/test-case-a3-b2.py` — test suite (all tests must pass)
- `output/case-a3-b2-results.json`
- `output/case-a3-b2-whitepaper.md`
- `output/case-a3-b2-tectonic-heatmap.png` — sub-test 1 signal strength grid (Figure 1)
- `output/case-a3-b2-midcrustal-binplots.png` — sub-test 2 NH/SH mid-crustal bin distributions (Figure 2)
- `output/case-a3-b2-phase-alignment.png` — sub-test 3 peak-phase comparison (Figure 3)
- `output/case-a3-b2-declustering-sensitivity.png` — sub-test 1 declustering sensitivity (Figure 4)
- `output/case-a3-b2-interval1-threshold.png` — sub-test 4 Interval 1 SH threshold sweep (Figure 5)

---

## 1. Environment and data loading

In `src/case-a3-b2-analysis.py`:

- Imports: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define paths:
  ```python
  RAW_PATH   = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
  GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
  GK_PATH    = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv"
  REAS_PATH  = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv"
  ```
- Load full catalog; assert n=9210; parse `event_at` as UTC; log row count
- Load GSHHG classification; assert n=9210; no NaN in `dist_to_coast_km`
- Load G-K mainshocks; assert n=5883; log
- Load Reasenberg mainshocks; assert n=8265; log
- Merge GSHHG classification (`ocean_class`, `dist_to_coast_km`) onto all four DataFrames by `usgs_id`; assert no NaN in `dist_to_coast_km` after merge for each
- Compute phase for all four DataFrames:
  ```python
  JULIAN_YEAR_SECS = 31_557_600.0
  df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
  ```
- Define constants:
  ```python
  ADAPTIVE_K_THRESHOLDS = {"k24": 500, "k16": 200, "k12": 100}
  TECTONIC_CLASSES = ["continental", "transitional", "oceanic"]
  MIDCRUSTAL_DEPTH_MIN = 20.0
  MIDCRUSTAL_DEPTH_MAX = 70.0
  INTERVAL_BINS = {"interval_1": [4, 5], "interval_2": [15], "interval_3": [21]}
  INTERVAL1_THRESHOLDS = [0.33, 0.40, 0.45, 0.50]
  CATALOGS = {"full": df_full, "gk": df_gk, "reas": df_reas}
  ```
- Apply hemisphere split to each catalog:
  ```python
  def split_hemisphere(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
      """Return (df_nh, df_sh, n_equatorial)."""
      n_equatorial = int((df["latitude"] == 0).sum())
      return df[df["latitude"] > 0].copy(), df[df["latitude"] < 0].copy(), n_equatorial
  ```
- Apply tectonic class split from `ocean_class` column (values: `"continental"`, `"transitional"`, `"oceanic"`)

---

## 2. Core statistics helper

```python
def adaptive_k(n: int) -> int | None:
    """Return bin count per adaptive-k rule; None if n < 100."""
    if n >= 500:
        return 24
    elif n >= 200:
        return 16
    elif n >= 100:
        return 12
    return None
```

```python
def compute_chi2_stats(phases: np.ndarray) -> dict:
    """
    Compute chi-square uniformity test on solar phase distribution.

    Parameters
    ----------
    phases : np.ndarray
        Array of phase values in [0, 1).

    Returns
    -------
    dict with keys: n, k, low_n, chi2, p_chi2, cramers_v, bin_counts,
                    peak_bin, peak_phase, interval_1_z, interval_2_z, interval_3_z.
    """
    n = len(phases)
    k = adaptive_k(n)
    if k is None:
        return {"n": n, "k": None, "low_n": True, "chi2": None, "p_chi2": None,
                "cramers_v": None, "bin_counts": None, "peak_bin": None,
                "peak_phase": None, "interval_1_z": None, "interval_2_z": None,
                "interval_3_z": None}
    bin_counts = np.bincount(np.floor(phases * k).astype(int) % k, minlength=k)
    expected = n / k
    chi2_stat, p_chi2 = scipy.stats.chisquare(bin_counts, np.full(k, expected))
    cramers_v = np.sqrt(chi2_stat / (n * (k - 1)))
    peak_bin = int(np.argmax(bin_counts))
    peak_phase = (peak_bin + 0.5) / k
    interval_z = {}
    for iname, bins in INTERVAL_BINS.items():
        obs_i = sum(bin_counts[b] for b in bins if b < k)
        exp_i = len([b for b in bins if b < k]) * expected
        interval_z[iname] = float((obs_i - exp_i) / np.sqrt(exp_i)) if exp_i > 0 else None
    return {
        "n": n, "k": k, "low_n": False,
        "chi2": float(chi2_stat), "p_chi2": float(p_chi2),
        "cramers_v": float(cramers_v),
        "bin_counts": bin_counts.tolist(),
        "peak_bin": peak_bin,
        "peak_phase": float(peak_phase),
        "interval_1_z": interval_z["interval_1"],
        "interval_2_z": interval_z["interval_2"],
        "interval_3_z": interval_z["interval_3"],
    }
```

---

## 3. Sub-test 1 — Tectonic-matched hemisphere comparison

```python
def run_tectonic_hemisphere_comparison(
    catalogs: dict[str, pd.DataFrame],
    tectonic_classes: list[str],
) -> dict:
```

For each catalog name in `["full", "gk", "reas"]`:
  For each tectonic class in `["continental", "transitional", "oceanic"]`:
    For each hemisphere in `["nh", "sh"]`:
      - Filter: `df[(df["ocean_class"] == tclass) & (df["latitude"] [>/< ] 0)]`
      - Call `compute_chi2_stats(subset["phase"].values)`
      - Record result dict under key `catalog → tectonic_class → hemisphere`

Return nested dict with 18 cells. Example structure:
```json
{
  "full": {
    "continental": {
      "nh": {"n": int, "k": int, "chi2": float, "p_chi2": float, "cramers_v": float, ...},
      "sh": {...}
    },
    "transitional": {"nh": {...}, "sh": {...}},
    "oceanic":      {"nh": {...}, "sh": {...}}
  },
  "gk":   {...},
  "reas": {...}
}
```

After computing all cells, for each catalog annotate each tectonic-class row with:
- `"nh_significant"`: bool (p_chi2 < 0.05)
- `"sh_significant"`: bool (p_chi2 < 0.05)
- `"both_significant"`: bool
- `"neither_significant"`: bool

Store under results key `"subtest_1_tectonic_hemisphere"`.

---

## 4. Sub-test 2 — Mid-crustal hemisphere split

```python
def run_midcrustal_hemisphere_split(catalogs: dict[str, pd.DataFrame]) -> dict:
```

For each catalog in `["full", "gk", "reas"]`:
  - Filter to mid-crustal band: `(df["depth"] >= 20.0) & (df["depth"] < 70.0)`
  - For each hemisphere in `["nh", "sh"]`:
    - Filter by hemisphere and call `compute_chi2_stats`
  - Also compute global (NH+SH combined) mid-crustal stats as regression anchor — full catalog mid-crustal global result should reproduce A3.B4's χ²≈85.48, p≈4.0×10⁻⁹ (±2.0 tolerance)
  - Record under catalog key

Return structure:
```json
{
  "full": {
    "global": {"n": int, "chi2": float, "p_chi2": float, "cramers_v": float, ...},
    "nh":     {"n": int, "chi2": float, "p_chi2": float, "cramers_v": float, ...},
    "sh":     {"n": int, "chi2": float, "p_chi2": float, "cramers_v": float, ...}
  },
  "gk":   {...},
  "reas": {...}
}
```

Store under results key `"subtest_2_midcrustal_hemisphere"`.

---

## 5. Sub-test 3 — Phase alignment comparison

```python
def run_phase_alignment(
    subtest1: dict,
    subtest2: dict,
) -> dict:
```

Collect all significant NH/SH pairs from sub-tests 1 and 2. A pair is formed from one NH cell and its corresponding SH cell (same catalog, same tectonic class or same mid-crustal band).

For each pair where **at least one** of NH or SH is significant (p_chi2 < 0.05) and has a valid `peak_phase`:
- Compute wrapped phase offset:
  ```python
  delta = (nh_peak_phase - sh_peak_phase + 0.5) % 1.0 - 0.5
  # Result in [-0.5, 0.5]; 0 = in-phase, ±0.5 = anti-phase
  ```
- Classify alignment:
  - `|delta| < 0.083` (2 bins at k=24): `"in_phase"` (consistent with solar-geometric mechanism)
  - `|delta| > 0.417` (within 2 bins of anti-phase): `"anti_phase"` (consistent with loading mechanism)
  - Otherwise: `"offset"`
- Record: `{"source": str, "catalog": str, "nh_peak_phase": float, "sh_peak_phase": float, "delta_phase": float, "alignment": str, "nh_significant": bool, "sh_significant": bool}`

Summarize:
```json
{
  "n_pairs_evaluated": int,
  "n_in_phase": int,
  "n_anti_phase": int,
  "n_offset": int,
  "dominant_alignment": "in_phase | anti_phase | offset | mixed",
  "pairs": [...]
}
```

Store under results key `"subtest_3_phase_alignment"`.

---

## 6. Sub-test 4 — Revised Interval 1 SH threshold sensitivity

```python
def run_interval1_threshold_sensitivity(df_full: pd.DataFrame) -> dict:
```

Using the full catalog only. Two SH populations:
- `sh_all`: all SH events (`latitude < 0`)
- `sh_continental`: SH events with `ocean_class == "continental"`

For each population × each threshold t in `[0.33, 0.40, 0.45, 0.50]`:
1. Compute k=24 bin counts for SH phases
2. Compute the fraction of Interval 1 bins that are elevated (observed > expected + sqrt(expected)):
   - `expected = n_sh / 24`
   - `elevated_bins = [b for b in [4, 5] if bin_counts[b] > expected + np.sqrt(expected)]`
   - `overlap_fraction = len(elevated_bins) / 2` (Interval 1 spans 2 bins at k=24)
3. Classify: `"present"` if `overlap_fraction >= t`, else `"absent"`

Record for each population × threshold:
```json
{
  "threshold": float,
  "n_sh": int,
  "overlap_fraction": float,
  "interval_1_classification": "present | absent",
  "bin_4_obs": int, "bin_5_obs": int, "expected": float
}
```

Also note the threshold at which classification flips (if any), separately for `sh_all` and `sh_continental`. A flip point ≤ 0.40 indicates Interval 1 is marginal; no flip (absent at all thresholds) indicates robust SH absence; present at all thresholds indicates the A2.B1 result was driven by composition.

Store under results key `"subtest_4_interval1_threshold"`.

---

## 7. Results JSON structure

```json
{
  "case": "A3.B2",
  "title": "Hemisphere Stratification Refinement",
  "parameters": {
    "n_catalog": 9210,
    "n_gk": 5883,
    "n_reas": 8265,
    "julian_year_secs": 31557600.0,
    "adaptive_k_thresholds": {"k24": 500, "k16": 200, "k12": 100},
    "tectonic_class_boundaries_km": {"continental_max": 50, "transitional_max": 200},
    "midcrustal_depth_range_km": [20.0, 70.0],
    "interval1_threshold_sweep": [0.33, 0.40, 0.45, 0.50],
    "phase_alignment_bin_tolerance": 2
  },
  "catalog_sizes": {
    "full": {"n_nh": int, "n_sh": int, "n_equatorial": int},
    "gk":   {"n_nh": int, "n_sh": int, "n_equatorial": int},
    "reas": {"n_nh": int, "n_sh": int, "n_equatorial": int}
  },
  "subtest_1_tectonic_hemisphere": { /* 3 catalogs × 3 classes × 2 hemispheres */ },
  "subtest_2_midcrustal_hemisphere": { /* 3 catalogs × global+NH+SH */ },
  "subtest_3_phase_alignment": { /* pairs summary */ },
  "subtest_4_interval1_threshold": {
    "sh_all":          [/* 4 threshold records */],
    "sh_continental":  [/* 4 threshold records */],
    "sh_all_flip_threshold": float | null,
    "sh_continental_flip_threshold": float | null
  }
}
```

---

## 8. Visualizations

In `src/visualization-case-a3-b2.py`:

Import: `matplotlib`, `matplotlib.pyplot`, `numpy`, `json`, `pathlib`, `logging`

---

**Figure 1 — Tectonic × hemisphere signal heatmap** (`output/case-a3-b2-tectonic-heatmap.png`):
- 3-row × 6-column grid: rows = tectonic classes (continental, transitional, oceanic); columns = NH and SH for each of three catalogs (full, G-K, Reasenberg) arranged as triplets: [Full-NH, Full-SH | GK-NH, GK-SH | Reas-NH, Reas-SH]
- Each cell: Cramér's V as fill color (white=0, red=max); annotate with χ², p-value, and n; bold border for p<0.05 cells; gray diagonal striping for low-n (n<100) cells
- Column headers: "Full NH", "Full SH", "G-K NH", "G-K SH", "Reas NH", "Reas SH"
- Row labels: tectonic class names
- Color scale: white-to-red gradient; shared scale across all cells; include colorbar labeled "Cramér's V"
- Title: "Tectonic-Matched Hemisphere χ² Signal Strength (A3.B2)"
- 300 DPI, figsize=(16, 7)

---

**Figure 2 — Mid-crustal hemisphere bin distributions** (`output/case-a3-b2-midcrustal-binplots.png`):
- 2-row × 3-column grid: rows = NH and SH; columns = full catalog, G-K, Reasenberg
- Each panel: horizontal steelblue bar chart at k=24 (or adaptive k if smaller); dashed expected-count line; 1-SD threshold dashed line; gray shaded bands for A1b Interval 1 (bins 4-5), Interval 2 (bin 15), Interval 3 (bin 21)
- Annotate each panel: n, χ², p, Cramér's V, k used; peak bin marker (vertical dotted line)
- Row labels: "Northern Hemisphere (20–70 km)", "Southern Hemisphere (20–70 km)"
- Column labels: "Full catalog", "G-K mainshocks", "Reasenberg mainshocks"
- Title: "Mid-Crustal Solar Phase Distribution by Hemisphere (A3.B2)"
- 300 DPI, figsize=(15, 9)

---

**Figure 3 — Phase alignment comparison** (`output/case-a3-b2-phase-alignment.png`):
- Scatter plot: x-axis = NH peak phase [0,1], y-axis = SH peak phase [0,1]
- One point per evaluated pair; color by source: tectonic class pairs (one color per class) and mid-crustal pairs (distinct color)
- Marker style: filled circle if both NH and SH significant; half-filled (open) if only one significant
- Reference lines: y=x diagonal (in-phase) as solid gray; y=x+0.5 and y=x-0.5 (anti-phase) as dashed gray
- Annotate each point with short label (e.g., "Cont-Full", "Mid-GK")
- Include shading: ±2-bin (0.083 phase unit) tolerance band around y=x in light blue; ±2-bin band around anti-phase lines in light orange
- Legend: pair type colors, significance styles, reference line meanings
- Title: "NH vs. SH Peak Phase — Alignment Test (A3.B2)"
- 300 DPI, figsize=(8, 7)

---

**Figure 4 — Declustering sensitivity on tectonic-matched comparison** (`output/case-a3-b2-declustering-sensitivity.png`):
- 2-row × 3-column grid: rows = NH and SH; columns = continental, transitional, oceanic
- Each panel: grouped bar chart; x-axis = catalog (Full, G-K, Reasenberg); y-axis = Cramér's V; steelblue bars; error bars not required; horizontal dashed line at Cramér's V for full catalog result; annotate bars with p-value (if p<0.05 mark with *, if p<0.01 mark with **)
- Title: "Declustering Sensitivity — Tectonic × Hemisphere (A3.B2)"
- 300 DPI, figsize=(14, 8)

---

**Figure 5 — Interval 1 SH threshold sensitivity** (`output/case-a3-b2-interval1-threshold.png`):
- 1-row × 2-column grid: columns = "All SH events" and "Continental SH events only"
- Each panel:
  - x-axis: overlap threshold (0.33, 0.40, 0.45, 0.50)
  - y-axis: Interval 1 overlap fraction (0.0 to 1.0); horizontal dashed line at each threshold value in matching color
  - Single horizontal line showing the Interval 1 overlap fraction (constant — it does not vary with threshold; only the classification boundary moves)
  - Fill below the fraction line: steelblue if fraction ≥ threshold ("present"), light gray if fraction < threshold ("absent")
  - Annotate with "present" / "absent" classification at each threshold tick
  - Include n for the population in panel subtitle
- Title: "Interval 1 SH Presence — Threshold Sensitivity (A3.B2)"
- 300 DPI, figsize=(10, 5)

---

## 9. Test suite

In `tests/test-case-a3-b2.py`:

- `test_catalog_loads`: assert full n=9210, GK n=5883, Reasenberg n=8265; assert all have `phase` column with values in [0.0, 1.0)
- `test_gshhg_merge`: assert no NaN in `dist_to_coast_km` after merge for full, GK, Reasenberg catalogs
- `test_hemisphere_split_full`: assert n_nh + n_sh + n_equatorial == 9210 for full catalog; assert n_nh > 0 and n_sh > 0
- `test_tectonic_partition_full`: assert continental + transitional + oceanic == 9210 for full catalog
- `test_adaptive_k_logic`: assert `adaptive_k(500)` == 24; `adaptive_k(300)` == 16; `adaptive_k(150)` == 12; `adaptive_k(50)` is None
- `test_subtest1_cell_count`: assert results JSON `subtest_1_tectonic_hemisphere` contains exactly 3 catalogs, each with 3 tectonic classes, each with 2 hemispheres (18 leaf cells total)
- `test_subtest1_continental_nh_significant`: assert `full → continental → nh → p_chi2 < 0.05` (scientific assertion: continental NH should be significant given A3.B3 finding that continental class is significant globally)
- `test_subtest2_midcrustal_global_regression`: assert `subtest_2_midcrustal_hemisphere → full → global → chi2` is within 2.0 of 85.48 (A3.B4 anchor); assert p_chi2 < 1e-7
- `test_subtest2_cell_count`: assert `subtest_2_midcrustal_hemisphere` contains 3 catalogs, each with `global`, `nh`, `sh` keys
- `test_phase_alignment_delta_range`: for all pairs in `subtest_3_phase_alignment`, assert `|delta_phase| <= 0.5` (wrapped correctly)
- `test_phase_alignment_classification_valid`: assert all `alignment` values are one of `"in_phase"`, `"anti_phase"`, `"offset"`
- `test_interval1_threshold_count`: assert `subtest_4_interval1_threshold → sh_all` has 4 records (one per threshold); same for `sh_continental`
- `test_interval1_classification_valid`: assert all `interval_1_classification` values are `"present"` or `"absent"`
- `test_results_json_completeness`: load results JSON; assert all four `subtest_*` keys present; assert `catalog_sizes` present with full, gk, reas entries; assert `parameters` present
- `test_output_figures_exist`: assert all 5 PNG files exist and size > 50 KB

---

## 10. Whitepaper

In `output/case-a3-b2-whitepaper.md`:

Standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Code [model name]).

### Sections:

1. **Abstract** (150–200 words): state that A2.B1's raw NH/SH comparison is confounded by tectonic composition; describe the four-sub-test design that controls for tectonic class (sub-test 1) and depth (sub-test 2); state whether the signal asymmetry persists after controlling for these variables; state the mechanistic implication of the phase alignment result (sub-test 3); and state the Interval 1 SH threshold sensitivity outcome (sub-test 4).

2. **Data Source**: ISC-GEM full catalog (n=9,210); GSHHG tectonic classification (three classes at T_outer=200 km); G-K mainshocks (n=5,883) and Reasenberg mainshocks (n=8,265) as declustering sensitivity layers. State hemisphere sizes for the full catalog (n_nh, n_sh).

3. **Methodology**
   - 3.1 Phase-normalized binning: Julian year constant; cite `data-handling.md`
   - 3.2 Tectonic classification: GSHHG baseline boundaries (continental ≤50 km, transitional 50–200 km, oceanic >200 km); cite A3.B3
   - 3.3 Adaptive k rule: describe the n-threshold logic; cite A3.B4 for prior application
   - 3.4 Sub-test 1 — tectonic-matched hemisphere comparison: describe design (18 cells); note oceanic is borderline-confidence per A3.B3
   - 3.5 Sub-test 2 — mid-crustal hemisphere split: describe depth filter (20–70 km from A3.B4); describe regression anchor check against A3.B4 global result
   - 3.6 Sub-test 3 — phase alignment: describe wrapped delta-phase computation; describe in-phase/anti-phase classification thresholds (2-bin tolerance)
   - 3.7 Sub-test 4 — Interval 1 SH threshold sensitivity: describe four thresholds (33%–50%); describe all-SH vs. continental-SH comparison design

4. **Results**
   - 4.1 Tectonic-matched hemisphere comparison: embed Figure 1; tabulate which cells achieve p<0.05; describe whether NH/SH asymmetry persists within tectonic classes; compare to A3.B3 global-class result
   - 4.2 Mid-crustal hemisphere split: embed Figure 2; state NH and SH mid-crustal χ², p, Cramér's V for full catalog; compare signal strength between hemispheres; state declustering suppression rate
   - 4.3 Phase alignment: embed Figure 3; describe whether significant NH and SH cells are in-phase or anti-phase; state dominant alignment classification; note any exceptions
   - 4.4 Declustering sensitivity: embed Figure 4; describe whether tectonic-class results are stable across G-K and Reasenberg; compare suppression patterns to A3.B1 findings
   - 4.5 Interval 1 SH threshold sensitivity: embed Figure 5; state overlap fractions for all-SH and continental-SH at each threshold; state whether classification flips; state mechanistic interpretation

5. **Cross-Topic Comparison**
   - **Hemisphere Stratification — Phase Symmetry Test (A2.B1):** The A2.B1 raw NH/SH comparison established the three-interval structure and the Interval 1 SH threshold-sensitive absence. A3.B2 extends this by stratifying within tectonic class, showing whether the A2.B1 hemispheric asymmetry persists after tectonic-composition control.
   - **Ocean/Coast Sequential Threshold Sensitivity (A3.B3):** A3.B3 established that the solar-phase signal is concentrated in continental and transitional classes (p=0.0005 and p=0.0003 respectively). A3.B2 inherits these boundaries and tests whether the same tectonic-class signal is hemispherically symmetric or skewed.
   - **Depth × Magnitude Stratification with Moho Isolation (A3.B4):** A3.B4 localized the signal to the mid-crustal band (20–70 km, p=4×10⁻⁹). A3.B2 sub-test 2 applies this depth window to the hemisphere split to test whether the dominant depth signal is globally symmetric.
   - **Declustering Sensitivity Analysis (A2.A4):** A2.A4 established that G-K declustering substantially suppresses the signal. A3.B2 sub-test 1's declustering sensitivity layer confirms whether tectonic-class hemisphere results are equally suppressed.

6. **Interpretation**: state whether the NH/SH asymmetry in the solar-phase signal is explained by tectonic composition (if NH-only significance disappears within matched tectonic classes), or whether it persists (if NH/SH difference remains even within continental or mid-crustal subsets). State the mechanistic implication of the phase alignment result. Guard against both confirming a hemisphere-specific loading mechanism and dismissing one — the sample sizes in the stratified cells are substantially smaller than the full catalog and may reduce statistical power.

7. **Limitations**: stratified cell sample sizes are reduced (particularly SH continental and SH transitional, given the SH's more oceanic character); adaptive k reduces to k=16 or k=12 in some cells, limiting phase resolution; the GSHHG 200 km outer threshold is the A3.B3 baseline and has not been swept in this case; declustering catalogs (G-K and Reasenberg) are included as sensitivity but not primary; the phase alignment comparison uses peak bin as a coarse proxy for interval center.

8. **References**
   - Yeager, J. (2026). A2.B1: Hemisphere Stratification — Phase Symmetry Test. erebus-vee-two internal report.
   - Yeager, J. (2026). A2.A4: Aftershock Phase-Preference Analysis. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B1: Rolling-Window Chi-Square Repeat. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B3: Ocean/Coast Sequential Threshold Sensitivity. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation. erebus-vee-two internal report.

---

## 11. Update context docs

After all outputs are generated and tests pass:

- Append to `topic-a3/docs/topic-summary.md`:
  ```
  ## Case A3.B2: Hemisphere Stratification Refinement
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [n_nh, n_sh for full catalog; tectonic-matched hemisphere: which cells significant (catalog × class × hemisphere); mid-crustal NH p and SH p for full catalog; phase alignment dominant classification; Interval 1 SH threshold flip point for all-SH and continental-SH]
  ```
- Update Case B2 status from `Planning` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a3/CLAUDE.md` Case Table
