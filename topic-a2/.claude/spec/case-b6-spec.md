# Case B6: Rolling Window Stationarity Test

**Status:** Pending

**Intent statement:**
Case B6 tests whether the solar-phase signal observed in the full ISC-GEM catalog is stationary across the 72-year record (1950–2021) or concentrated in a subset of years. Bradley & Hubbard (2024) documented that tidal correlations detected in pre-2000 data consistently failed replication in post-2000 data — the canonical demonstration of non-stationarity as an artifact signature. Dutilleul et al. (2021) found that Parkfield's semiannual periodicity shifted phase after the 2004 M6.0 earthquake. A non-stationary solar signal would substantially weaken the interpretation of the Case 3A result and would deprioritize the mechanistic investigations in downstream cases. The test uses a sliding 10-year window to compute the Rayleigh statistic, mean vector length (R), and mean phase angle of the solar-phase distribution at each window position. The 1970s windows should be flagged for anomalous behavior given the Adhoc A0b finding that ISC-GEM has a 414-event density spike in that decade relative to ComCat, indicating a period-specific difference in catalog construction methodology.

**Relationship to prior topics:**
Topic L3 established the global ISC-GEM catalog (n=9,210, M≥6.0, 1950–2021) and its ephemeris columns. Topic L5 demonstrated that signal strength varies substantially after declustering on ComCat (60.2% chi-square reduction), motivating concern about whether the signal is similarly temporally concentrated. B6 is elevated to second in execution order precisely because a non-stationary result would call into question all downstream mechanistic cases — paralleling the role A4 plays for the aftershock contamination question.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `latitude`, `longitude`, `depth` |

Parse `event_at` as UTC datetime for windowing. Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year (31,557,600 s) unless per-year values are available [confirm before running: same as A4 — verify whether per-year solar year lengths are pre-computed or whether the Julian constant is used uniformly].

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-b6-analysis.py` — rolling window computation: Rayleigh R, p-value, mean phase angle per window; flags 1970s windows; writes results JSON
- `src/visualization-case-b6.py` — generates all PNG figures
- `tests/test-case-b6.py` — test suite
- `output/case-b6-results.json` — all window statistics, full trajectory, stationarity classification
- `output/case-b6-whitepaper.md` — methodology, results, stationarity interpretation
- `output/case-b6-trajectory.png` — time trajectory of Rayleigh R, p-value, and mean phase across all windows
- `output/case-b6-phase-stability.png` — circular plot of mean phase angle per window, colored by decade

---

## 1. Environment and data loading

In `src/case-b6-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define `RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"`
- Load CSV; parse `event_at` as `pd.Timestamp` (UTC); log row count; assert n=9210
- Compute `event_year` from `event_at` as integer calendar year
- Compute phase: `phase = (solar_secs / 31_557_600.0) % 1.0` (Julian year constant; see confirm note above)
- Define `WINDOW_YEARS = 10` and `STEP_YEARS = 1` for the sliding window
- Define window start years: `range(1950, 2021 - WINDOW_YEARS + 1)` → 62 windows (1950–1959 through 2012–2021)
- Flag 1970s windows: any window whose range overlaps [1970, 1979] (i.e., start years 1961–1979 have some 1970s years; flag start years 1970–1979 as "1970s window")

---

## 2. Rolling window statistics

For each window start year `y` (covering `[y, y + WINDOW_YEARS)`):

1. Filter catalog to `event_year >= y` and `event_year < y + WINDOW_YEARS`
2. Record `n_window = len(subset)`
3. Compute Rayleigh statistic:
   - `angles = 2 * np.pi * subset["phase"]`
   - `mean_cos = np.mean(np.cos(angles))`
   - `mean_sin = np.mean(np.sin(angles))`
   - `R = np.sqrt(mean_cos**2 + mean_sin**2)` (mean vector length)
   - `rayleigh_z = n_window * R**2`
   - `p_rayleigh = np.exp(-rayleigh_z)` (Rayleigh p-value approximation; valid for n > 10)
4. Compute mean phase angle: `mean_angle_rad = np.arctan2(mean_sin, mean_cos)` → convert to phase fraction: `mean_phase = (mean_angle_rad / (2 * np.pi)) % 1.0`
5. Compute supplementary chi-square at k=24 (apply phase-normalized binning): `chi2_stat, p_chi2 = scipy.stats.chisquare(observed_counts, expected_counts)`
6. Record `{"window_start": y, "window_end": y + WINDOW_YEARS - 1, "n": n_window, "rayleigh_R": R, "rayleigh_z": rayleigh_z, "p_rayleigh": p_rayleigh, "mean_phase": mean_phase, "chi2_k24": chi2_stat, "p_chi2_k24": p_chi2, "is_1970s_window": bool}`

---

## 3. Stationarity classification

After computing all windows:

1. Count windows where `p_rayleigh < 0.05` → `n_significant_windows`
2. Count windows where `p_rayleigh < 0.05 / 62` (Bonferroni threshold) → `n_bonferroni_significant`
3. Compute phase stability: std of `mean_phase` across all windows (circular standard deviation):
   - `circ_var = 1 - R_mean_of_window_means` where `R_mean_of_window_means` is the mean resultant length of the window mean phase angles
   - `circ_std_deg = sqrt(-2 * log(1 - circ_var)) * 180 / pi`
4. Classify stationarity:
   - "stationary": >= 70% of windows significant and circ_std_deg < 20°
   - "partially stationary": 30–70% of windows significant, or circ_std_deg 20–40°
   - "non-stationary": < 30% of windows significant or circ_std_deg > 40°
5. Check for 1970s anomaly: compute mean R for windows flagged as 1970s vs non-1970s; compute ratio; flag if ratio > 1.5

Output under key `"stationarity"` in results JSON:
```json
{
  "stationarity": {
    "n_windows": 62,
    "n_significant_rayleigh_p05": int,
    "n_bonferroni_significant": int,
    "circular_std_deg": float,
    "classification": "stationary | partially stationary | non-stationary",
    "1970s_mean_R": float,
    "non_1970s_mean_R": float,
    "1970s_anomaly_ratio": float,
    "1970s_anomaly_flagged": bool
  }
}
```

Output per-window records under key `"windows"` in results JSON:
```json
{
  "windows": [
    {"window_start": int, "window_end": int, "n": int, "rayleigh_R": float, "rayleigh_z": float,
     "p_rayleigh": float, "mean_phase": float, "chi2_k24": float, "p_chi2_k24": float, "is_1970s_window": bool},
    ...
  ]
}
```

---

## 4. Visualizations

In `src/visualization-case-b6.py`:

**Figure 1 — Time trajectory** (`output/case-b6-trajectory.png`):
- 3-row stacked line plot sharing the x-axis (window center year)
- Row 1: Rayleigh R (mean vector length), y-range [0.0, 0.1]; steelblue line
- Row 2: Rayleigh p-value (log scale, y-range [0.001, 1.0]); red line; horizontal dashed lines at p=0.05 and p=0.001
- Row 3: mean phase fraction (0–1); black line with points; horizontal dashed lines at the three A1b baseline phase centers (0.22, 0.64, 0.90)
- Mark 1970s windows with vertical shaded band (light yellow)
- x-axis ticks: every 5 years (1955, 1960, ..., 2015); label "Window center year"
- 300 DPI, publication quality

**Figure 2 — Phase stability circular plot** (`output/case-b6-phase-stability.png`):
- Circular/polar plot where each window is a point on the unit circle at angle `2π * mean_phase`
- Point radius proportional to Rayleigh R (scaled to visible range)
- Points colored by decade: 1950s=blue, 1960s=steelblue, 1970s=orange, 1980s=green, 1990s=purple, 2000s=red, 2010s=black
- Label the three A1b baseline intervals as arcs on the circle perimeter (gray shaded wedges)
- Legend for decades; title "Mean Phase Angle per 10-Year Window"
- 300 DPI

---

## 5. Test suite

In `tests/test-case-b6.py`:

- `test_catalog_load`: assert n=9210 rows loaded; assert `event_at` column parses without NaT values
- `test_phase_range`: assert all computed phases in [0.0, 1.0)
- `test_window_count`: assert exactly 62 windows are produced (1950–1959 through 2012–2021)
- `test_window_n_positive`: assert all window n values > 0; assert no window has n < 50 (warn if any window < 100)
- `test_rayleigh_uniform`: generate 500 uniform random phases; assert Rayleigh p > 0.05 in 95% of 100 bootstrap trials
- `test_rayleigh_concentrated`: generate 500 phases concentrated at 0.2 (Gaussian noise σ=0.02); assert Rayleigh p < 0.001
- `test_rayleigh_R_bounds`: assert all window Rayleigh R values in [0.0, 1.0]
- `test_mean_phase_bounds`: assert all window mean_phase values in [0.0, 1.0)
- `test_1970s_flagging`: assert windows with `window_start` in 1970–1979 are flagged as `is_1970s_window=True`; assert windows with `window_start` in 1980–2012 are flagged `False`
- `test_stationarity_classification_present`: load `output/case-b6-results.json`; assert key `"stationarity"` is present and `"classification"` is one of the three valid values
- `test_chi2_k24_all_windows`: assert all `chi2_k24` values are non-negative floats

---

## 6. Whitepaper

In `output/case-b6-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (is the 72-year solar-phase signal stationary?), sliding window approach, stationarity classification result, and implication for downstream cases
2. **Data Source** — ISC-GEM catalog description (n=9,210, M≥6.0, 1950–2021); note the Adhoc A0b 1970s density spike finding and its relevance to this case
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1 and `rules/data-handling.md`)
   - 3.2 Sliding window design: 10-year window, 1-year step, 62 windows from 1950–2021; choice of window size rationale
   - 3.3 Rayleigh statistic and mean phase angle computation
   - 3.4 Chi-square supplementary test at k=24
   - 3.5 Stationarity classification criteria
   - 3.6 1970s window flagging rationale (cite Adhoc A0b)
4. **Results**
   - 4.1 Trajectory results: embed Figure 1; tabulate % of windows with p < 0.05; describe R and phase trajectories
   - 4.2 Phase stability: embed Figure 2; report circular std of mean phase; describe whether phase drifts across decades
   - 4.3 Stationarity classification: state final classification and supporting metrics
   - 4.4 1970s anomaly check: report ratio of 1970s vs non-1970s mean R; state whether flagged
5. **Cross-Topic Comparison** — compare to Bradley & Hubbard (2024) finding of pre/post-2000 replication failure; compare to Dutilleul et al. (2021) Parkfield phase shift; state how B6 results confirm or contrast these patterns
6. **Interpretation** — state whether the signal is temporally robust or concentrated; note implications for the downstream program; maintain objectivity regarding confirmation and hypothesis bias
7. **Limitations** — 10-year window has limited statistical power for low-n decades; early decades (1950s–1960s) have fewer events; ISC-GEM 1970s density spike could inflate signal for those windows
8. **References** — Bradley & Hubbard (2024), Dutilleul et al. (2021), Mardia & Jupp (2000), Adhoc A0b, Adhoc A1

---

## 7. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case B6: Rolling Window Stationarity Test
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [stationarity classification; % windows significant at p<0.05; circular std of mean phase; 1970s anomaly ratio; R trajectory summary]
  ```
- Update Case B6 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
