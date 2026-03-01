# Case A3: Magnitude Stratification of the Solar Signal

**Status:** Pending

**Intent statement:**
Case A3 stratifies the ISC-GEM catalog into four magnitude bands (M 6.0–6.4, M 6.5–6.9, M 7.0–7.4, M ≥ 7.5) and independently computes solar-phase statistics for each band. The goal is to determine whether the Case 3A solar signal is evenly distributed across magnitude or concentrated in a specific band, which directly constrains the physical mechanism. The tidal literature provides competing predictions: Métivier (2009) found tidal triggering inversely dependent on magnitude (strongest for small events); Ide et al. (2016) found M < 8.2 indistinguishable from random for tidal phase; Hao et al. (2018) found diurnal periodicity stronger at larger magnitudes in Japan. These findings point in opposite directions at different timescales, making any positive result at global M ≥ 6 scale a novel finding. Adhoc A1b provides preliminary directional evidence: the elevated-bin events (n=1,438, 15.6% of catalog) showed no magnitude skew relative to the full catalog, suggesting magnitude-independence. Case A3 formalizes this with per-band chi-square tests and Cramér's V. The ISC-GEM catalog's 82.9% 2-decimal magnitude precision (vs ComCat's 77.5%) makes per-band stratification more reliable because magnitude rounding artifacts that would inflate the M 6.0 bin in ComCat are substantially reduced.

**Relationship to prior topics:**
Topic L3 established magnitude column `usgs_mag` in the ISC-GEM catalog. Topic L5 (Declustering Degradation) documented how declustering affects chi-square signal, providing context for interpreting band-specific results. Adhoc A1b's no-magnitude-skew finding in elevated-bin events is the primary prior-topic result directly informing this case's expected outcome. Adhoc A0 established ISC-GEM's magnitude precision advantage (82.9% 2-decimal vs ComCat's 77.5% 1-decimal), which specifically motivated placing b-value and magnitude-stratification tests in this topic rather than replicating them on ComCat.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `latitude`, `longitude`, `depth` |

Magnitude bands derived at analysis time from `usgs_mag`:
- Band 1: `6.0 <= usgs_mag < 6.5` (M 6.0–6.4)
- Band 2: `6.5 <= usgs_mag < 7.0` (M 6.5–6.9)
- Band 3: `7.0 <= usgs_mag < 7.5` (M 7.0–7.4)
- Band 4: `usgs_mag >= 7.5` (M ≥ 7.5)

Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year [confirm before running: consistent with A4/B6/A1/B1 — verify uniform Julian constant vs per-year values].

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-a3-analysis.py` — main analysis: magnitude stratification, per-band chi-square/Rayleigh/Cramér's V, trend analysis, writes results JSON
- `src/visualization-case-a3.py` — generates all PNG figures
- `tests/test-case-a3.py` — test suite
- `output/case-a3-results.json` — per-band statistics, effect-size trend, prediction matching
- `output/case-a3-whitepaper.md` — methodology, results, cross-topic comparison
- `output/case-a3-binplots.png` — four-panel bin distribution plot, one per magnitude band at k=24
- `output/case-a3-effect-trend.png` — Cramér's V and Rayleigh R vs magnitude band (trend plot)

---

## 1. Environment and data loading

In `src/case-a3-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define `RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"`
- Load catalog; assert n=9210; log row count
- Define magnitude band boundaries and labels:
  ```python
  BANDS = [
      {"label": "M6.0-6.4", "min": 6.0, "max": 6.5},
      {"label": "M6.5-6.9", "min": 6.5, "max": 7.0},
      {"label": "M7.0-7.4", "min": 7.0, "max": 7.5},
      {"label": "M7.5+",    "min": 7.5, "max": 99.0},
  ]
  ```
- For each band, create subset: `df_band = df[(df["usgs_mag"] >= band["min"]) & (df["usgs_mag"] < band["max"])]`
- Log band sizes; assert sum of band sizes == 9210 (all events assigned to exactly one band)
- Compute phase for each subset using Julian year constant

---

## 2. Per-band statistics

For each of four magnitude bands at each of three bin counts (k=16, 24, 32):

1. Compute phase and bin indices (same as A4)
2. Observed counts `O`, expected counts `E = n_band / k`
3. Chi-square: `chi2, p_chi2 = scipy.stats.chisquare(O, E)`
4. Cramér's V: `V = sqrt(chi2 / (n * (k-1)))`
5. Rayleigh R and p-value (same as A4)
6. Mean phase angle fraction (same as A4)
7. Elevated bins at 1-SD threshold; elevated phase intervals

Output structure under key `"band_stats"` in results JSON:
```json
{
  "band_stats": {
    "M6.0-6.4": {
      "n": int,
      "k16": {"chi2": float, "p_chi2": float, "cramer_v": float, "rayleigh_R": float, "p_rayleigh": float, "mean_phase": float, "bin_counts": [...], "elevated_intervals": [...]},
      "k24": {...},
      "k32": {...}
    },
    "M6.5-6.9": {...},
    "M7.0-7.4": {...},
    "M7.5+": {...}
  }
}
```

---

## 3. Effect-size trend analysis

Using results at k=24 as primary:

1. Extract Cramér's V and Rayleigh R for each band in ascending magnitude order
2. Compute Spearman rank correlation of Cramér's V vs band index (1–4); record rho and p-value
3. Compute Spearman rank correlation of Rayleigh R vs band index; record rho and p-value
4. Classify trend direction:
   - "increasing with magnitude": Spearman rho > 0.5 (positive correlation, effect grows with M)
   - "decreasing with magnitude": Spearman rho < -0.5 (negative correlation, effect shrinks with M)
   - "flat": |rho| <= 0.5 (no clear directional trend)
5. Check whether any band has p_chi2 < 0.05; record which bands are significant
6. Check whether M ≥ 7.5 band has p_chi2 < 0.05 (key discriminator for hydrological hypothesis)

Prediction matching:
- **Hydrological loading**: effect weakens or disappears at M ≥ 7.5 → evaluate: M7.5+ band p_chi2 >= 0.05 AND trend "decreasing"
- **Solar-geometric forcing**: effect flat or increases with magnitude → evaluate: trend "flat" or "increasing" AND M7.5+ band significant
- **Consistent with tidal literature pattern**: effect weakest at M ≥ 6 base → evaluate: M6.0–6.4 shows strongest signal; higher bands weaker or absent
- **Magnitude-independent (A1b preliminary)**: Cramér's V approximately uniform across bands → evaluate: max V / min V < 1.5

Output under key `"trend_analysis"` in results JSON:
```json
{
  "trend_analysis": {
    "cramer_v_by_band": {"M6.0-6.4": float, "M6.5-6.9": float, "M7.0-7.4": float, "M7.5+": float},
    "rayleigh_R_by_band": {...},
    "spearman_rho_cramer_v": float,
    "spearman_p_cramer_v": float,
    "spearman_rho_rayleigh": float,
    "spearman_p_rayleigh": float,
    "trend_classification": "increasing | decreasing | flat",
    "significant_bands": ["M6.0-6.4", ...],
    "m75_significant": bool,
    "prediction_support": {
      "hydrological_loading": "supported | partially | not supported",
      "solar_geometric": "supported | partially | not supported",
      "tidal_literature_pattern": "supported | partially | not supported",
      "magnitude_independent": "supported | partially | not supported"
    }
  }
}
```

---

## 4. Visualizations

In `src/visualization-case-a3.py`:

**Figure 1 — Four-panel bin distributions** (`output/case-a3-binplots.png`):
- 2×2 panel layout, one panel per magnitude band (M6.0–6.4, M6.5–6.9, M7.0–7.4, M≥7.5), all at k=24
- Each panel: horizontal bar chart (steelblue); dashed expected-count line; 1-SD threshold dashed line
- Mark A1b baseline intervals (phases 0.1875–0.25, 0.625–0.656, 0.875–0.917) with gray shaded bands
- Annotate each panel with band label, n, χ², p-value, Cramér's V
- Shared x-axis: phase fraction 0–1; label equinox (0.19, 0.69) and solstice (0.44, 0.94) positions
- 300 DPI

**Figure 2 — Effect-size trend plot** (`output/case-a3-effect-trend.png`):
- Dual-axis line chart: x-axis = magnitude band (M6.0–6.4, M6.5–6.9, M7.0–7.4, M≥7.5 as ordinal positions 1–4)
- Left y-axis: Cramér's V (steelblue line with points); right y-axis: Rayleigh R (orange line with points)
- Error bars on Cramér's V: use 95% CI from bootstrap (1000 bootstrap resamplings of each band); bootstrap is done in analysis script
- Mark significance (p < 0.05) with filled circles; non-significant with open circles
- Annotate Spearman rho and p-value for Cramér's V trend
- 300 DPI

---

## 5. Bootstrap confidence intervals for Cramér's V

Add to `src/case-a3-analysis.py`:

For each magnitude band at k=24:
- Set `np.random.seed(42)`
- Generate 1000 bootstrap resamples of the band: `resample = np.random.choice(len(band), size=len(band), replace=True)`
- For each resample, compute chi-square and Cramér's V
- Record 2.5th and 97.5th percentile of bootstrap Cramér's V distribution as 95% CI
- Add to `"band_stats"` per band: `"cramer_v_ci95_lower": float, "cramer_v_ci95_upper": float`

---

## 6. Test suite

In `tests/test-case-a3.py`:

- `test_catalog_load`: assert n=9210 rows
- `test_band_partition`: assert sum of band sizes == 9210 (no events unassigned or double-counted)
- `test_band_sizes_positive`: assert all four bands have n > 0; log warn if any band < 100
- `test_phase_range`: assert all phases in [0.0, 1.0)
- `test_chi_square_all_bands`: assert chi2 and p_chi2 are finite floats for all four bands at k=16, 24, 32
- `test_cramer_v_range`: assert Cramér's V in [0.0, 1.0] for all bands
- `test_bootstrap_ci_order`: assert cramer_v_ci95_lower <= cramer_v <= cramer_v_ci95_upper for all bands
- `test_spearman_rho_range`: assert Spearman rho in [-1.0, 1.0]
- `test_trend_classification_valid`: assert trend_classification is one of "increasing", "decreasing", "flat"
- `test_prediction_support_valid`: assert all prediction_support values are one of "supported", "partially", "not supported"
- `test_results_json_keys`: load `output/case-a3-results.json`; assert keys "band_stats" and "trend_analysis" present; assert band_stats has all four band keys
- `test_m75_band_n`: assert M7.5+ band n > 200 (sanity check: ISC-GEM should have ~400+ M≥7.5 events)

---

## 7. Whitepaper

In `output/case-a3-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (does solar signal vary with magnitude?), four bands tested, trend direction found, and primary mechanistic implication
2. **Data Source** — ISC-GEM catalog (n=9,210); note ISC-GEM magnitude precision advantage (82.9% 2-decimal vs ComCat 77.5%); report band sizes
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1, `rules/data-handling.md`)
   - 3.2 Magnitude band definitions and rationale (M 6.0–6.4, 6.5–6.9, 7.0–7.4, ≥7.5)
   - 3.3 Chi-square, Rayleigh, Cramér's V computation per band
   - 3.4 Bootstrap Cramér's V CI procedure (1000 resamples, seed 42)
   - 3.5 Effect-size trend: Spearman rank correlation and classification
   - 3.6 Prediction framework for mechanism discrimination
4. **Results**
   - 4.1 Per-band distributions: embed Figure 1; tabulate chi2, p, Cramér's V, and Rayleigh R for all four bands at k=24; note which bands are significant
   - 4.2 Effect-size trend: embed Figure 2; state Spearman rho and trend classification
   - 4.3 Prediction matching: evaluate each of the four competing predictions against the results
5. **Cross-Topic Comparison** — compare to Adhoc A1b no-magnitude-skew finding in elevated-bin events; compare to Métivier (2009) inverse tidal magnitude dependence; compare to Ide et al. (2016) no M≥6 tidal signal; compare to Hao et al. (2018) larger-magnitude diurnal signal
6. **Interpretation** — state mechanism implication from trend direction; note any unexpected results (e.g., M7.5+ significance would be novel); maintain objectivity
7. **Limitations** — M7.5+ band has smaller sample (~400 events) limiting statistical power; chi-square sensitivity decreases with smaller n; declustering not applied (uses raw catalog); phase normalization uses Julian year approximation
8. **References** — Hao et al. (2018), Ide et al. (2016), Métivier et al. (2009), Johnson et al. (2017), Adhoc A1b, Adhoc A0

---

## 8. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case A3: Magnitude Stratification of the Solar Signal
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [Cramér's V per band at k=24; Spearman rho trend; trend classification; which bands are significant; prediction support summary]
  ```
- Update Case A3 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
