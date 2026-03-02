# Case A2: b-Value Seasonal Variation

**Status:** Pending

**Intent statement:**
Case A2 divides the ISC-GEM catalog into solar-cycle phase bins and computes the Gutenberg-Richter b-value independently for each bin using maximum likelihood estimation (MLE), testing whether b-value varies significantly across the annual solar cycle. Colledge et al. (2025) found annual b-value variation in Nepal (~0.1 per kPa) peaking during the monsoon loading phase, with no corresponding robust seasonal rate variation. Ide et al. (2016) established that b-value decreases as tidal shear stress amplitude increases for very large events — confirming that b-value responds to stress perturbations even when seismicity rate does not. If rate and b-value vary inversely in phase with the solar cycle (b-value peaks near solstice when rate is low, troughs near equinox when rate is high), the equinox excess in rate is partly a magnitude-distribution effect: more events capable of reaching the M 6.0 threshold during unclamping, not only a nucleation rate increase. This would be a new finding at global M ≥ 6 scale. The ISC-GEM catalog's 82.9% 2-decimal magnitude precision (established in Adhoc A0) makes b-value MLE substantially more reliable than ComCat (77.5% 1-decimal). Bin count is selected from Adhoc Case A1's finding that solar signal is strongest at k=24 and k=32.

**Relationship to prior topics:**
Adhoc A0 established the ISC-GEM magnitude precision advantage that specifically motivates this case. Adhoc A1 identified k=24 and k=32 as the bin counts with strongest solar signal in chi-square, directly informing the bin-count selection for the b-value phase binning. Case 3A found the overall solar signal using rate (event count); A2 tests the complementary b-value dimension. Topic L3 established magnitude completeness for the ISC-GEM catalog at the global M ≥ 6.0 threshold used here.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `latitude`, `longitude`, `depth` |

Magnitude column: `usgs_mag` — used directly for b-value MLE. No magnitude band pre-filtering; full M ≥ 6.0 catalog used. Minimum magnitude for b-value MLE: `Mc = 6.0` (catalog threshold).

Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year [confirm before running: consistent with prior cases — verify uniform Julian constant vs per-year values].

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-a2-analysis.py` — main analysis: phase binning, b-value MLE per bin, bootstrap CI, phase-variation significance test, writes results JSON
- `src/visualization-case-a2.py` — generates all PNG figures
- `tests/test-case-a2.py` — test suite
- `output/case-a2-results.json` — per-bin b-values, CIs, phase-variation test results
- `output/case-a2-whitepaper.md` — methodology, results, cross-topic comparison
- `output/case-a2-bvalue-phase.png` — b-value vs phase bin with CI bars at k=24
- `output/case-a2-bvalue-k32.png` — same at k=32
- `output/case-a2-rate-bvalue-overlay.png` — overlay of rate signal (event count per bin) and b-value signal at k=24

---

## 1. Environment and data loading

In `src/case-a2-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define `RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"`
- Load catalog; assert n=9210; assert no NaN in `usgs_mag`; log row count
- Compute phase using Julian year constant
- Set `BIN_COUNTS = [24, 32]` (primary=24; k=32 as secondary per Adhoc A1 recommendation)
- Assign bin index: `bin_idx = np.floor(phase * k).astype(int)` for each k

---

## 2. b-Value MLE per phase bin

**Maximum likelihood b-value estimation (Aki 1965 formula):**
- For a set of magnitudes `{m_i}` with completeness magnitude `Mc`:
  - `b_mle = log10(e) / (mean(m_i) - Mc)`
  - Standard error: `sigma_b = b_mle / sqrt(n)` (Shi & Bolt 1982 approximation)
  - This is the standard MLE for the Gutenberg-Richter relation

For each bin count k (24, 32) and each bin index j (0 to k-1):

1. Extract magnitudes for events in bin j: `mags_j = df[df["bin_idx_k"] == j]["usgs_mag"].values`
2. Record `n_j = len(mags_j)`
3. If `n_j < 20`: flag bin as low-n; still compute b-value but mark with `"low_n_flag": true`
4. Compute b-value MLE: `b_j = np.log10(np.e) / (np.mean(mags_j) - 6.0)`
5. Compute standard error: `se_j = b_j / np.sqrt(n_j)`
6. Bootstrap 95% CI for b_j:
   - Set `np.random.seed(42)` on first use
   - Generate 1000 bootstrap resamples of `mags_j`
   - Compute b-value MLE for each resample
   - Record 2.5th and 97.5th percentiles as `b_ci95_lower`, `b_ci95_upper`

Record per-bin: `{"bin_idx": j, "n": n_j, "b_mle": float, "b_se": float, "b_ci95_lower": float, "b_ci95_upper": float, "low_n_flag": bool}`

---

## 3. Phase-variation significance test

After computing b-values for all bins:

1. Compute `b_mean = np.mean([b_j for all bins])` and `b_std = np.std([b_j for all bins])`
2. Peak and trough identification:
   - `b_max_bin = argmax(b_j)` — bin with highest b-value; record its phase center `phase_center = (b_max_bin + 0.5) / k`
   - `b_min_bin = argmin(b_j)` — bin with lowest b-value; record its phase center
   - `b_range = b_max - b_min` (peak-to-trough range)
3. Test whether b-value varies across phase bins using one-way ANOVA:
   - Groups = magnitudes in each phase bin
   - `F_stat, p_anova = scipy.stats.f_oneway(*[mags_in_bin_j for j in range(k)])`
4. Test inverse-phase relationship with rate signal:
   - Load or recompute event count per bin (`count_j = len(df[df["bin_idx"] == j])` for each j)
   - Compute Pearson correlation between `count_j` and `b_j`: `r_rate_b, p_rate_b = scipy.stats.pearsonr(counts, b_values)`
   - A negative correlation indicates inverse phase relationship (high rate ↔ low b-value)
5. Test whether b-value peak aligns with solstice and trough with equinox:
   - Equinox phases: ~0.19 (March) and ~0.69 (September); Solstice phases: ~0.44 (June) and ~0.94 (December)
   - Classify b_max_bin phase: `"near_solstice"` if within 0.1 of 0.44 or 0.94; `"near_equinox"` if within 0.1 of 0.19 or 0.69; `"other"` otherwise
   - Classify b_min_bin phase analogously

Output structure under key `"phase_variation"` in results JSON (per k):
```json
{
  "phase_variation": {
    "k24": {
      "bins": [{"bin_idx": int, "phase_center": float, "n": int, "b_mle": float, "b_se": float,
                "b_ci95_lower": float, "b_ci95_upper": float, "low_n_flag": bool}, ...],
      "b_mean": float, "b_std": float, "b_range": float,
      "b_max_phase_center": float, "b_max_phase_class": "near_solstice | near_equinox | other",
      "b_min_phase_center": float, "b_min_phase_class": "near_solstice | near_equinox | other",
      "f_stat": float, "p_anova": float,
      "r_rate_b": float, "p_rate_b": float,
      "inverse_phase_relationship": bool
    },
    "k32": {...}
  }
}
```

---

## 4. Visualizations

In `src/visualization-case-a2.py`:

**Figure 1 — b-Value vs phase at k=24** (`output/case-a2-bvalue-phase.png`):
- Line chart: x-axis = phase bin center (0 to 1, 24 equally spaced points); y-axis = b-value
- Steelblue line connecting bin b-values; error bars showing 95% CI from bootstrap
- Dashed horizontal line at b_mean
- Gray shaded bands at A1b baseline intervals (phases 0.1875–0.25, 0.625–0.656, 0.875–0.917)
- Vertical dotted lines at equinox positions (0.19, 0.69) and solstice positions (0.44, 0.94)
- Annotate: ANOVA F-stat and p-value; Pearson r with rate signal; label b_max and b_min bins with their phase class
- Low-n bins marked with open circles (vs filled); include "(low n < 20)" in legend
- 300 DPI

**Figure 2 — b-Value vs phase at k=32** (`output/case-a2-bvalue-k32.png`):
- Same format as Figure 1 but at k=32 (32 bins)
- 300 DPI

**Figure 3 — Rate and b-value overlay** (`output/case-a2-rate-bvalue-overlay.png`):
- Dual-axis chart at k=24
- Left y-axis (steelblue): event count per bin (bar chart, steelblue bars)
- Right y-axis (orange): b-value per bin (line with points, orange)
- Both x-axes aligned at phase 0–1
- Annotate Pearson r between rate and b-value; note direction (positive = in-phase, negative = inverse)
- 300 DPI

---

## 5. Test suite

In `tests/test-case-a2.py`:

- `test_catalog_load`: assert n=9210; assert no NaN in usgs_mag
- `test_phase_range`: assert all phases in [0.0, 1.0)
- `test_bin_partition_k24`: assert sum of per-bin counts at k=24 == 9210
- `test_bin_partition_k32`: assert sum of per-bin counts at k=32 == 9210
- `test_b_mle_formula`: for a synthetic magnitude set [6.0, 6.1, 6.2, 6.3, 6.4] with Mc=6.0, assert b_mle ≈ log10(e) / 0.2 ≈ 2.172 (within 0.01)
- `test_b_mle_range`: assert all computed b_mle values in [0.5, 3.0] (physically reasonable range for global M≥6 catalog)
- `test_bootstrap_ci_contains_mle`: for each bin, assert b_ci95_lower <= b_mle <= b_ci95_upper
- `test_low_n_flag_trigger`: assert bins with n < 20 have low_n_flag=True; assert bins with n >= 20 have low_n_flag=False
- `test_anova_p_nonnegative`: assert p_anova in [0.0, 1.0] for both k=24 and k=32
- `test_pearson_r_range`: assert r_rate_b in [-1.0, 1.0]
- `test_phase_class_valid`: assert b_max_phase_class and b_min_phase_class are one of "near_solstice", "near_equinox", "other"
- `test_results_json_keys`: load `output/case-a2-results.json`; assert key "phase_variation" present with sub-keys "k24" and "k32"

---

## 6. Whitepaper

In `output/case-a2-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (does b-value vary with solar phase?), MLE method, bins tested, ANOVA and correlation results, and diagnostic implication for mechanism
2. **Data Source** — ISC-GEM catalog (n=9,210, M≥6.0); note ISC-GEM magnitude precision advantage (82.9% 2-decimal vs ComCat 77.5%) and its relevance to MLE reliability (cite Adhoc A0)
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1, `rules/data-handling.md`); bin count selection rationale (k=24, 32 from Adhoc A1 signal peak)
   - 3.2 Gutenberg-Richter b-value MLE: Aki (1965) formula; completeness threshold Mc=6.0; Shi & Bolt (1982) standard error
   - 3.3 Bootstrap 95% CI: 1000 resamples, seed 42
   - 3.4 ANOVA across bins: test for significant phase-dependent b-value variation
   - 3.5 Rate-b correlation: Pearson r between per-bin event count and per-bin b-value; inverse relationship interpretation
   - 3.6 Phase alignment classification: solstice/equinox labels for b_max and b_min bins
4. **Results**
   - 4.1 b-value at k=24: embed Figure 1; report b_mean, b_range, b_max phase class, b_min phase class; note low-n bins
   - 4.2 b-value at k=32: embed Figure 2; compare to k=24 results
   - 4.3 Rate-b overlay: embed Figure 3; report Pearson r and direction; state whether inverse-phase relationship is supported
   - 4.4 Statistical tests: report ANOVA F-stat and p-value for both k values; state significance conclusion
5. **Cross-Topic Comparison** — compare to Colledge et al. (2025) Nepal b-value seasonal variation (~0.1 per kPa); compare to Ide et al. (2016) tidal b-value response in large events; note that global M≥6 result is distinct from regional Nepal or tidal-period findings
6. **Interpretation** — if inverse rate-b relationship confirmed: interpret as magnitude-distribution effect (more events crossing M6.0 threshold at unclamping phase); if not confirmed: discuss implications; maintain objectivity regarding hypothesis bias
7. **Limitations** — b-value MLE requires sufficient n per bin (flag low-n bins); completeness magnitude Mc=6.0 assumes uniform completeness across phases and decades (not tested here); declustering not applied
8. **References** — Aki (1965), Colledge et al. (2025), Bettinelli et al. (2008), Ide et al. (2016), Shi & Bolt (1982), Adhoc A0, Adhoc A1

---

## 7. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case A2: b-Value Seasonal Variation
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [b_mean, b_range at k=24; ANOVA p-value; Pearson r between rate and b-value; b_max and b_min phase classes; inverse-phase relationship supported/not]
  ```
- Update Case A2 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
