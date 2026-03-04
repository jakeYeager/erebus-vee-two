# Case A3.B1: Rolling-Window Chi-Square Repeat

**Status:** Complete

**Source case:** A2.B6

**Intent statement:**
Case A3.B1 repeats the A2.B6 rolling-window stationarity analysis with chi-square (k=24) promoted to the primary statistic and Rayleigh demoted to secondary. In A2.B6, chi-square reached significance in 71% of rolling windows versus Rayleigh's 38.7% — a 32-percentage-point divergence flagged in the B6 console summary as evidence of "multi-modal within-window structure that the unimodal Rayleigh test misses." Because the A2.B6 stationarity conclusion (classification: non-stationary) rested on Rayleigh, the divergence was not resolved, and its implications for all downstream mechanistic cases were deferred. A3.B1 resolves that deferral: by making chi-square the primary criterion, it determines whether the phase-distribution non-uniformity is temporally consistent and whether the apparent non-stationarity was a consequence of applying a unimodal statistic to a multi-modal distribution.

Additionally, A3.B1 introduces two extensions not present in B6: (1) interval-level tracking within each window — which of the three A1b baseline intervals (Interval 1 ~0.19–0.25, Interval 2 ~0.625–0.656, Interval 3 ~0.875–0.917) are elevated — to identify whether non-stationarity is global or interval-specific; and (2) sequence density analysis across the three declustered catalogs, checking whether windows with elevated chi-square also show elevated aftershock sequence counts (a direct diagnostic for 2004-Sumatra-style contamination).

**Note on window parameters:** The A3.B1 case details document references "5-year window, 1-year stride, 62 windows" in the open-questions block. This is a transcription error. A2.B6 used a 10-year window, 1-year stride, 62 windows (range(1950, 2012) → 62 start years). A3.B1 retains the A2.B6 window parameters unchanged: **10-year window, 1-year stride, 62 windows.**

**Relationship to prior topics:**
A2.B6 established the 62-window rolling trajectory and classified the signal as non-stationary based on Rayleigh (38.7% windows significant, circular SD = 84.4°). The chi-square divergence in B6 is the direct motivation for this case. A2.A4 established that aftershock sequences are phase-preferring and that the 2003–2014 windows are the most significant in B6 — contemporaneous with the 2004 Sumatra M9.1 aftershock sequence. A3.B1's sequence density extension directly tests whether those window elevations are sequence-driven.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |
| G-K mainshocks (enriched) | `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv` | 5,883 | same + `aftershock_count`, `foreshock_count`, `window_secs`, `window_km` |
| Reasenberg mainshocks (enriched) | `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv` | 8,265 | same + `aftershock_count`, `foreshock_count`, `window_secs`, `window_km` |
| A1b mainshocks (enriched) | `data/iscgem/declustering-algorithm/mainshocks_a1b-seq_global.csv` | 7,137 | same + `aftershock_count`, `foreshock_count`, `window_secs`, `window_km` |

Parse `event_at` as UTC datetime. Phase normalization: `phase = (solar_secs / 31_557_600.0) % 1.0` (Julian year constant, same as A2.B6).

**A1b baseline intervals (k=24 bin mapping):**
- Interval 1: bins 4–5 (phase [0.1667, 0.2500)) — March equinox region (~0.19–0.25)
- Interval 2: bin 15 (phase [0.6250, 0.6667)) — mid-August region (~0.625–0.656)
- Interval 3: bin 21 (phase [0.8750, 0.9167)) — late-November region (~0.875–0.917)

Bin index (0-based) for k=24: `bin_i = floor(phase * 24)`.

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a3/`
- All output paths: `BASE_DIR / "output" / ...`
- All data paths: `BASE_DIR.parent / "data" / "iscgem" / ...`

**Planned outputs:**
- `src/case-a3-b1-analysis.py` — rolling window computation for all 4 catalogs; chi-square primary, Rayleigh secondary; interval tracking; sequence density; writes results JSON
- `src/visualization-case-a3-b1.py` — generates all PNG figures
- `tests/test-case-a3-b1.py` — test suite
- `output/case-a3-b1-results.json` — all per-window statistics, per-catalog stationarity classifications, interval elevation summaries
- `output/case-a3-b1-whitepaper.md` — methodology, results, stationarity reinterpretation
- `output/case-a3-b1-trajectory.png` — chi-square primary + Rayleigh secondary + mean phase time trajectory (raw catalog)
- `output/case-a3-b1-catalog-comparison.png` — chi-square p-value trajectory overlay for all 4 catalogs
- `output/case-a3-b1-interval-heatmap.png` — per-window z-score heatmap for three A1b intervals
- `output/case-a3-b1-phase-stability.png` — circular plot of mean phase angle per window (raw catalog)

---

## 1. Environment and data loading

In `src/case-a3-b1-analysis.py`:

- Imports: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define catalog paths:
  ```python
  CATALOGS = {
      "raw":        BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv",
      "gk":         BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv",
      "reasenberg": BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv",
      "a1b":        BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_a1b-seq_global.csv",
  }
  EXPECTED_N = {"raw": 9210, "gk": 5883, "reasenberg": 8265, "a1b": 7137}
  ```
- For each catalog:
  - Load CSV; parse `event_at` as `pd.Timestamp` (UTC); log row count; assert n == `EXPECTED_N[key]`
  - Compute `event_year = event_at.dt.year` (integer)
  - Compute `phase = (solar_secs / 31_557_600.0) % 1.0`
- Define: `WINDOW_YEARS = 10`, `STEP_YEARS = 1`, `K_BINS = 24`
- Define window start years: `list(range(1950, 2021 - WINDOW_YEARS + 1))` → 62 windows
- Define 1970s window flag: `window_start` in range(1970, 1980)
- Define A1b interval bin ranges (0-based, k=24):
  ```python
  INTERVAL_BINS = {
      "interval_1": [4, 5],   # phase [0.1667, 0.2500)
      "interval_2": [15],     # phase [0.6250, 0.6667)
      "interval_3": [21],     # phase [0.8750, 0.9167)
  }
  ```

---

## 2. Rolling window statistics

For each catalog and each window start year `y` (covering `[y, y + WINDOW_YEARS)`):

1. Filter catalog to `event_year >= y` and `event_year < y + WINDOW_YEARS`
2. Record `n_window = len(subset)`
3. Compute phase-normalized bin counts (k=24):
   - `bin_indices = np.floor(subset["phase"] * K_BINS).astype(int) % K_BINS`
   - `observed = np.bincount(bin_indices, minlength=K_BINS)` (length-24 array)
   - `expected = np.full(K_BINS, n_window / K_BINS)`
4. **Chi-square (primary):**
   - `chi2_stat, p_chi2 = scipy.stats.chisquare(observed, expected)`
5. **Rayleigh (secondary):**
   - `angles = 2 * np.pi * subset["phase"]`
   - `mean_cos = np.mean(np.cos(angles))`
   - `mean_sin = np.mean(np.sin(angles))`
   - `R = np.sqrt(mean_cos**2 + mean_sin**2)`
   - `rayleigh_z = n_window * R**2`
   - `p_rayleigh = np.exp(-rayleigh_z)`
   - `mean_angle_rad = np.arctan2(mean_sin, mean_cos)`
   - `mean_phase = (mean_angle_rad / (2 * np.pi)) % 1.0`
6. **Interval-level tracking:**
   For each interval in `INTERVAL_BINS`:
   - `interval_count = observed[bins].sum()` where `bins` is the bin index list
   - `interval_expected = len(bins) * (n_window / K_BINS)`
   - `interval_z = (interval_count - interval_expected) / np.sqrt(interval_expected)`
   - `interval_elevated = bool(interval_z > 1.96)`
   Record `interval_1_count`, `interval_1_z`, `interval_1_elevated`; same for intervals 2 and 3.
7. **Sequence density** (declustered catalogs only; skip for raw):
   - `mean_aftershock_count = subset["aftershock_count"].mean()` (float)
   - `total_aftershock_count = subset["aftershock_count"].sum()` (int)
   - `pct_isolated = (subset["aftershock_count"] == 0).mean()` — fraction of mainshocks with no claimed events
8. Record per-window dict:
   ```json
   {
     "window_start": int, "window_end": int, "n": int,
     "chi2_k24": float, "p_chi2_k24": float,
     "rayleigh_R": float, "rayleigh_z": float, "p_rayleigh": float,
     "mean_phase": float,
     "interval_1_count": int, "interval_1_z": float, "interval_1_elevated": bool,
     "interval_2_count": int, "interval_2_z": float, "interval_2_elevated": bool,
     "interval_3_count": int, "interval_3_z": float, "interval_3_elevated": bool,
     "mean_aftershock_count": float,   // null for raw
     "total_aftershock_count": int,    // null for raw
     "pct_isolated": float,            // null for raw
     "is_1970s_window": bool
   }
   ```

---

## 3. Stationarity classification

Performed separately for each catalog. Chi-square is the **primary statistic**; Rayleigh is secondary.

1. `n_significant_chi2_p05 = count(p_chi2_k24 < 0.05)`
2. `n_bonferroni_significant_chi2 = count(p_chi2_k24 < 0.05 / 62)`
3. `n_significant_rayleigh_p05 = count(p_rayleigh < 0.05)` (secondary, for comparison)
4. Compute circular phase stability (from Rayleigh mean phases):
   - `angles_w = [2π * w["mean_phase"] for w in windows]`
   - `mean_cos_w = mean(cos(angles_w))`, `mean_sin_w = mean(sin(angles_w))`
   - `R_of_means = sqrt(mean_cos_w^2 + mean_sin_w^2)`
   - `circ_var = 1 - R_of_means`
   - `circ_std_deg = sqrt(-2 * log(1 - circ_var)) * 180 / pi` (if circ_var < 1; else set to 180°)
5. **Classify stationarity** (chi-square primary):
   - `pct_chi2_sig = n_significant_chi2_p05 / 62`
   - "stationary": `pct_chi2_sig >= 0.70` AND `circ_std_deg < 20`
   - "partially stationary": `0.30 <= pct_chi2_sig < 0.70` OR `20 <= circ_std_deg <= 40`
   - "non-stationary": `pct_chi2_sig < 0.30` OR `circ_std_deg > 40`
   - When criteria conflict (e.g., chi2 pct ≥ 70% but circ_std > 40°), apply the more conservative (lower) classification and record both conflicting metrics.
6. **Interval stationarity summary** (per catalog):
   - `n_windows_interval_1_elevated = count(interval_1_elevated == True)`; same for 2 and 3
   - `pct_interval_1_elevated = n_windows_interval_1_elevated / 62`; same for 2 and 3
   - Classify each interval: "globally elevated" (≥70%), "partially elevated" (30–70%), "absent" (<30%)
7. **1970s anomaly check** (same as B6):
   - `mean_R_1970s = mean(rayleigh_R for is_1970s_window)` (using Rayleigh R for comparability with B6)
   - `mean_R_non1970s = mean(rayleigh_R for not is_1970s_window)`
   - `ratio = mean_R_1970s / mean_R_non1970s`; flag if ratio > 1.5

Output under key `"catalogs"` → per-catalog key → `"stationarity"`:
```json
{
  "catalogs": {
    "raw": {
      "stationarity": {
        "n_windows": 62,
        "n_significant_chi2_p05": int,
        "pct_significant_chi2_p05": float,
        "n_bonferroni_significant_chi2": int,
        "n_significant_rayleigh_p05": int,
        "pct_significant_rayleigh_p05": float,
        "circular_std_deg": float,
        "classification": "stationary | partially stationary | non-stationary",
        "classification_basis": "string explaining which criterion drove classification",
        "interval_1_pct_elevated": float,
        "interval_2_pct_elevated": float,
        "interval_3_pct_elevated": float,
        "interval_1_classification": "globally elevated | partially elevated | absent",
        "interval_2_classification": "globally elevated | partially elevated | absent",
        "interval_3_classification": "globally elevated | partially elevated | absent",
        "1970s_mean_R": float,
        "non_1970s_mean_R": float,
        "1970s_anomaly_ratio": float,
        "1970s_anomaly_flagged": bool
      },
      "windows": [...]
    },
    "gk": { ... },
    "reasenberg": { ... },
    "a1b": { ... }
  }
}
```

---

## 4. Sequence density correlation

For each declustered catalog (gk, reasenberg, a1b):

1. Compute Pearson correlation between `p_chi2_k24` and `mean_aftershock_count` across all 62 windows. Record `r_chi2_vs_aftershock`, `p_chi2_vs_aftershock`.
2. Identify the 10 windows with lowest `p_chi2_k24` (most significant chi-square). Record their `mean_aftershock_count` relative to the catalog-wide window mean.
3. Report `seq_density_elevated_in_high_chi2_windows`: bool — True if the top-10-chi2 windows have mean aftershock count > 1.5× the overall window mean aftershock count.

Output under `"catalogs"` → per-catalog key → `"sequence_density"`:
```json
{
  "sequence_density": {
    "r_chi2_vs_aftershock": float,
    "p_chi2_vs_aftershock": float,
    "mean_aftershock_count_all_windows": float,
    "mean_aftershock_count_top10_chi2_windows": float,
    "seq_density_elevated_in_high_chi2_windows": bool
  }
}
```

---

## 5. Visualizations

In `src/visualization-case-a3-b1.py`:

Load `output/case-a3-b1-results.json` and render all figures.

---

**Figure 1 — Chi-square + Rayleigh trajectory, raw catalog** (`output/case-a3-b1-trajectory.png`):
- 3-row stacked subplot sharing x-axis (window center year = `window_start + 5`)
- Row 1: Chi-square p-value (log scale, y-range [1e-12, 1.0]); steelblue line; horizontal dashed lines at p=0.05 (gray) and Bonferroni threshold p=0.05/62≈0.000806 (orange); label "Chi-square p (k=24) — primary"
- Row 2: Rayleigh p-value (log scale, y-range [0.001, 1.0]); red line; horizontal dashed line at p=0.05; label "Rayleigh p — secondary"
- Row 3: Mean phase fraction (0–1); black line with points; horizontal dashed lines at A1b interval centers (0.208 for Interval 1, 0.646 for Interval 2, 0.896 for Interval 3)
- Mark 1970s windows with vertical shaded band (light yellow)
- x-axis ticks every 5 years; label "Window center year"
- Title: "Raw Catalog — Chi-Square Primary (k=24)"
- 300 DPI, publication quality

---

**Figure 2 — Multi-catalog chi-square comparison** (`output/case-a3-b1-catalog-comparison.png`):
- Single panel; x-axis = window center year; y-axis = chi-square p-value (log scale, [1e-12, 1.0])
- Four lines, one per catalog:
  - Raw: steelblue (solid)
  - G-K: red (dashed)
  - Reasenberg: green (dotted)
  - A1b: purple (dash-dot)
- Horizontal dashed line at p=0.05 (gray)
- Mark 1970s windows with vertical shaded band (light yellow)
- Legend; title "Chi-Square p-Value Trajectory by Catalog (k=24)"
- 300 DPI

---

**Figure 3 — Interval-level elevation heatmap** (`output/case-a3-b1-interval-heatmap.png`):
- x-axis = window center year (62 positions); y-axis = 3 rows (Interval 1, Interval 2, Interval 3)
- Color = per-window interval z-score; use a white-to-red gradient (white=0, red=max) with contour lines at z=1.96 and z=3.0
- Separate subplot panel per catalog (4 panels stacked vertically, one per catalog), sharing x-axis
- Title each panel with catalog label and its interval stationarity classifications (e.g., "Raw — Int1: globally elevated, Int2: globally elevated, Int3: partially elevated")
- x-axis ticks every 5 years; label "Window center year"
- 300 DPI

---

**Figure 4 — Phase stability circular plot, raw catalog** (`output/case-a3-b1-phase-stability.png`):
- Same design as A2.B6 Figure 2: polar plot of mean phase angle per window, point radius ∝ Rayleigh R
- Points colored by decade: 1950s=blue, 1960s=steelblue, 1970s=orange, 1980s=green, 1990s=purple, 2000s=red, 2010s=black
- Gray shaded wedges on the circle perimeter marking the three A1b baseline intervals:
  - Interval 1: angle range [2π×0.1667, 2π×0.2500]
  - Interval 2: angle range [2π×0.6250, 2π×0.6667]
  - Interval 3: angle range [2π×0.8750, 2π×0.9167]
- Title "Mean Phase Angle per 10-Year Window — Raw Catalog"
- 300 DPI

---

## 6. Test suite

In `tests/test-case-a3-b1.py`:

- `test_catalog_load_raw`: assert n=9210, `event_at` parses without NaT
- `test_catalog_load_declustered`: assert n=5883 (gk), n=8265 (reasenberg), n=7137 (a1b); assert `aftershock_count` column present in all three
- `test_phase_range`: for all 4 catalogs, assert all computed phases in [0.0, 1.0)
- `test_window_count`: for all 4 catalogs, assert exactly 62 windows
- `test_window_n_positive`: for all 4 catalogs, assert all window n > 0; warn (log) if any window n < 100
- `test_chi2_nonneg`: for all 4 catalogs, assert all `chi2_k24` values ≥ 0
- `test_chi2_p_bounds`: for all 4 catalogs, assert all `p_chi2_k24` in [0.0, 1.0]
- `test_rayleigh_R_bounds`: for all 4 catalogs, assert all `rayleigh_R` in [0.0, 1.0]
- `test_mean_phase_bounds`: for all 4 catalogs, assert all `mean_phase` in [0.0, 1.0)
- `test_interval_z_type`: for all 4 catalogs, assert `interval_1_z`, `interval_2_z`, `interval_3_z` are floats
- `test_interval_elevated_type`: for all 4 catalogs, assert `interval_1_elevated`, `interval_2_elevated`, `interval_3_elevated` are bools
- `test_1970s_flagging`: assert windows with `window_start` in 1970–1979 are `is_1970s_window=True`; assert windows with `window_start` in 1980–2012 are `False`
- `test_stationarity_classification_all_catalogs`: load `output/case-a3-b1-results.json`; for each catalog key, assert `"classification"` is one of "stationary", "partially stationary", "non-stationary"
- `test_sequence_density_declustered_only`: for raw catalog, assert `mean_aftershock_count` is null/None in all window records; for gk/reasenberg/a1b, assert it is a non-negative float
- `test_chi2_uniform`: generate 1000 uniform random phases; compute chi-square at k=24; assert p > 0.05 in at least 90% of 100 bootstrap trials
- `test_chi2_concentrated`: generate 1000 phases concentrated in bins 15–16 (phase 0.625–0.667, σ=0.01); assert chi-square p < 0.001
- `test_b6_raw_chi2_pct_match`: load results JSON; assert raw catalog `pct_significant_chi2_p05` is between 0.65 and 0.80 (expected ~0.71 matching A2.B6's documented 71% chi-square significance rate)

---

## 7. Whitepaper

In `output/case-a3-b1-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and standard footer (Generated with Claude Code [model name]).

### Sections:

1. **Abstract** (150–200 words): state the question (does chi-square as primary change the stationarity classification from A2.B6?), summarize the rolling-window approach with interval tracking, state the key findings for raw vs. declustered catalogs, and state the implication for the interpretation of A2.B6's non-stationarity conclusion.

2. **Data Source**: describe all four catalogs (raw ISC-GEM n=9,210; G-K mainshocks n=5,883; Reasenberg n=8,265; A1b n=7,137); note the 1970s density spike (Adhoc A0b) context. Note that the declustered datasets are sequence-enriched, providing `aftershock_count` per mainshock row.

3. **Methodology**
   - 3.1 Phase-normalized binning: state Julian year constant used; cite A2.A1 and `data-handling.md` rules
   - 3.2 Sliding window design: 10-year window, 1-year step, 62 windows (1950–1959 through 2011–2020); note window-parameter retention from A2.B6 and correction of the "5-year" transcription error in the A3.B1 case details document
   - 3.3 Chi-square test (k=24) as primary: formula, degrees of freedom (k−1=23), significance threshold p<0.05 and Bonferroni-corrected threshold p<0.05/62
   - 3.4 Rayleigh statistic as secondary: formula for R, z, and mean phase; rationale for retaining as secondary (phase stability metric requires the angular mean)
   - 3.5 A1b interval-level tracking: define three intervals and their k=24 bin mappings; define elevation z-score formula; threshold z > 1.96
   - 3.6 Sequence density analysis: define `mean_aftershock_count` per window; describe Pearson correlation with chi-square p-value
   - 3.7 Stationarity classification criteria: state the 70%/30% thresholds and circular SD thresholds; describe conflict resolution rule when criteria disagree

4. **Results**
   - 4.1 Raw catalog trajectory: embed Figure 1; report `n_significant_chi2_p05`, `pct_significant_chi2_p05`; describe trajectory shape noting 2003–2014 elevated windows
   - 4.2 Stationarity comparison to A2.B6: table comparing chi-square % and Rayleigh % per catalog; state final chi-square-primary classifications; explicitly compare to A2.B6's Rayleigh-primary classification of "non-stationary"
   - 4.3 Interval-level stationarity: embed Figure 3; report per-interval classification (globally elevated / partially elevated / absent) for each catalog; describe which intervals drive chi-square significance across windows
   - 4.4 Multi-catalog comparison: embed Figure 2; describe how declustering suppresses or preserves chi-square significance across windows
   - 4.5 Sequence density results: report `r_chi2_vs_aftershock` and `seq_density_elevated_in_high_chi2_windows` per declustered catalog; describe whether 2003–2014 chi-square peaks are associated with elevated sequence density
   - 4.6 Phase stability: embed Figure 4; report `circular_std_deg`; describe whether mean phase drifts even when chi-square significance is consistent
   - 4.7 1970s anomaly: report ratios for all 4 catalogs; state whether flagged

5. **Cross-Topic Comparison**: compare to Bradley & Hubbard (2024) pre/post-2000 replication failure; compare to A2.B6 finding and explain why chi-square diverges from Rayleigh (multi-modal vs. unimodal); note A2.A4's aftershock phase-preference finding and its connection to sequence density results here

6. **Interpretation**: state whether switching to chi-square changes the stationarity narrative from A2.B6; address whether multi-modal distribution non-uniformity is temporally consistent even when unimodal Rayleigh is not; assess whether interval-level elevations are globally consistent or temporally concentrated; maintain objectivity against confirmation bias

7. **Limitations**: chi-square is sensitive to bin-edge effects (mitigated by phase-normalized binning); declustered catalog window n values smaller than raw, reducing chi-square power for some windows; interval z-score assumes Poisson counts (valid for large n); circular SD is derived from Rayleigh mean phase, which is a poor central tendency estimate for multi-modal distributions

8. **References**: Bradley & Hubbard (2024), Dutilleul et al. (2021), Mardia & Jupp (2000), Adhoc A0b, Adhoc A1, A2.B6, A2.A4

---

## 8. Update context docs

After all outputs are generated:

- Append to `topic-a3/CLAUDE.md` Case Table: update A3.B1 status from `Planning` to `Complete` (or `Blocked`/`Abandoned`)
- Create or append to `topic-a3/docs/topic-summary.md`:
  ```
  ## Case A3.B1: Rolling-Window Chi-Square Repeat
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [chi2 primary stationarity classification per catalog; pct chi2 windows significant; interval classifications (I1/I2/I3); sequence density correlation finding; comparison to A2.B6 non-stationary Rayleigh result]
  ```
