# Case A4: Declustering Sensitivity Analysis

**Status:** Pending

**Intent statement:**
Case A4 is the mandatory first-execution case for Topic A2. It tests whether the solar-phase signal observed in the full ISC-GEM catalog (Case 3A: χ²=45.61, p=6.13×10⁻⁵) survives aftershock removal under three independent declustering methods: Gardner-Knopoff (1974), Reasenberg (1985), and a data-informed custom window derived from Adhoc Case A1b (83.2 km / 95.6 days). Bradley & Hubbard (2024) demonstrated that aftershock temporal clustering is large enough to produce artificial periodic signal in Schuster-type tests; Park et al. (2021) showed standard Schuster requires prior declustering. This case applies that validation directly to the three-interval solar phase structure identified in Adhoc A1b. Three sub-analyses are performed: (A) scalar signal survival on mainshock catalogs, (B) post-declustering interval structure comparison, and (C) aftershock population phase preference. The three-way method comparison specifically tests whether G-K over-suppresses genuine signal by using windows (~295 days at M6.0) approximately 3× longer than the observed clustering footprint. Because this topic uses ISC-GEM, the ComCat-specific M 6.0 rounding artifact concern (Adhoc A0b) is mitigated, and ISC-GEM-based results are not expected to match the ComCat-based L4/L5 reference (60.2% suppression, χ²=45.61→18.13).

**Relationship to prior topics:**
Topic L4 established G-K declustering on ComCat (6,222 mainshocks / 3,580 aftershocks) and Topic L5 measured 60.2% chi-square degradation from declustering on ComCat. A4 replicates this measurement on the ISC-GEM catalog with two additional methods (Reasenberg and the A1b-informed window), providing ISC-GEM-specific suppression baselines. The A1b-informed window (83.2 km / 95.6 days) was derived from elevated-bin clustering metrics in Adhoc Case A1b, directly incorporating that topic's findings into the declustering design. The L4/L5 ComCat results serve as a cross-catalog reference point but should not be treated as a prediction for ISC-GEM outcomes.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `latitude`, `longitude`, `depth` |
| G-K mainshocks | `data/iscgem/declustering-algorithm/mainshocks_G-K_global.csv` | 5,883 | same schema as raw |
| G-K aftershocks | `data/iscgem/declustering-algorithm/aftershocks_G-K_global.csv` | 3,327 | same schema as raw |
| Reasenberg mainshocks | `data/iscgem/declustering-algorithm/mainshocks_reas_global.csv` | 8,265 | same schema as raw |
| Reasenberg aftershocks | `data/iscgem/declustering-algorithm/aftershocks_reas_global.csv` | 945 | same schema as raw |
| A1b mainshocks | `data/iscgem/declustering-algorithm/mainshocks_a1b_global.csv` | 7,137 | same schema as raw |
| A1b aftershocks | `data/iscgem/declustering-algorithm/aftershocks_a1b_global.csv` | 2,073 | same schema as raw |

Phase normalization: compute `phase = (solar_secs / solar_year_secs_for_event) % 1`, then `bin_index = floor(phase * k)` for bin counts k=16, 24, 32. The solar year length per event is the actual year length in seconds for the event's `solaration_year`.

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- Cross-topic output paths: `BASE_DIR.parent / "topic-l4" / "output" / ...` (if referencing prior topic outputs)
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-a4-analysis.py` — main analysis: loads all seven catalogs, runs sub-analyses A, B, C, writes results JSON
- `src/case-a4-sub-a.py` — Sub-analysis A: chi-square and Rayleigh on each mainshock catalog at k=16,24,32
- `src/case-a4-sub-b.py` — Sub-analysis B: elevated-bin interval identification on each mainshock catalog; interval structure comparison to A1b baseline
- `src/case-a4-sub-c.py` — Sub-analysis C: chi-square on each aftershock catalog; phase preference test
- `src/visualization-case-a4.py` — generates all PNG figures
- `tests/test-case-a4.py` — test suite for computed values
- `output/case-a4-results.json` — all statistics, sub-analysis results, interval comparisons
- `output/case-a4-whitepaper.md` — methodology, results, cross-topic comparison
- `output/case-a4-sub-a-binplot.png` — phase bin distributions for all four catalogs (raw + 3 mainshock), side-by-side at k=24
- `output/case-a4-sub-b-intervals.png` — elevated phase interval maps for each mainshock catalog vs A1b baseline
- `output/case-a4-sub-c-aftershock.png` — aftershock phase bin distributions for all three methods

---

## 1. Environment and data loading

In `src/case-a4-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define data paths using BASE_DIR:
  - `RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"`
  - `DECLUSTER_DIR = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm"`
  - Mainshock paths: `mainshocks_G-K_global.csv`, `mainshocks_reas_global.csv`, `mainshocks_a1b_global.csv`
  - Aftershock paths: `aftershocks_G-K_global.csv`, `aftershocks_reas_global.csv`, `aftershocks_a1b_global.csv`
- Load all seven CSVs; log row counts after load; assert counts match expected values (9210, 5883, 3327, 8265, 945, 7137, 2073) and log any discrepancy
- Assert partition integrity: for each method, `len(mainshocks) + len(aftershocks) == 9210` and zero `usgs_id` overlap between mainshock and aftershock files; raise a logged warning (not an error) if counts differ
- Define `SOLAR_YEAR_SECS_BY_YEAR: dict[int, float]` — use the Julian year approximation (365.25 × 86400 = 31,557,600 s) applied uniformly, unless `solaration_year`-specific values are available in the data [confirm before running: verify whether `solaration_year` column contains the actual calendar year integer and whether per-year solar year lengths are pre-computed anywhere in the project, or whether the Julian constant should be used uniformly]
- Phase normalization function: `def compute_phase(solar_secs: pd.Series, year_lengths: pd.Series) -> pd.Series` — returns `(solar_secs / year_lengths) % 1.0`

---

## 2. Sub-analysis A — Scalar signal survival (mainshock catalogs)

In `src/case-a4-sub-a.py` (callable from main analysis):

For each of four catalogs (raw, G-K mainshocks, Reasenberg mainshocks, A1b mainshocks), at each of three bin counts (k=16, 24, 32):

1. Compute phase for each event using the phase normalization function
2. Assign bin indices: `bin_idx = np.floor(phase * k).astype(int)`
3. Compute observed bin counts `O` (array of length k)
4. Expected counts: `E = np.full(k, len(catalog) / k)`
5. Chi-square statistic: `chi2_stat, p_val = scipy.stats.chisquare(O, E)`
6. Rayleigh statistic: compute mean resultant length `R = |sum(exp(2πi * phase))| / n`; Rayleigh p-value via `p_rayleigh = exp(-n * R²)` (Mardia & Jupp 2000 approximation)
7. Cramér's V: `V = sqrt(chi2_stat / (n * (k - 1)))`
8. Mean phase angle: `mean_angle = atan2(mean(sin(2π*phase)), mean(cos(2π*phase)))` in radians; convert to fraction of year (0–1)
9. Record for each catalog+k combination: `{catalog_name, n, k, chi2, p_chi2, cramer_v, rayleigh_R, p_rayleigh, mean_phase_fraction}`

Suppression calculation (Sub-analysis A summary):
- For each method and each k, compute `chi2_suppression_pct = (chi2_raw - chi2_mainshock) / chi2_raw * 100`
- Record alongside the per-catalog results

Output structure in `output/case-a4-results.json` under key `"sub_a"`:
```json
{
  "sub_a": {
    "raw": {"k16": {...}, "k24": {...}, "k32": {...}},
    "gk_mainshocks": {"k16": {...}, "k24": {...}, "k32": {...}},
    "reas_mainshocks": {"k16": {...}, "k24": {...}, "k32": {...}},
    "a1b_mainshocks": {"k16": {...}, "k24": {...}, "k32": {...}},
    "suppression_summary": [
      {"method": "gk", "k16": {...}, "k24": {...}, "k32": {...}},
      ...
    ]
  }
}
```

Each per-k record: `{"n": int, "k": int, "chi2": float, "p_chi2": float, "cramer_v": float, "rayleigh_R": float, "p_rayleigh": float, "mean_phase_fraction": float, "bin_counts": [int, ...]}`

---

## 3. Sub-analysis B — Post-declustering interval structure (mainshock catalogs)

In `src/case-a4-sub-b.py` (callable from main analysis):

**A1b baseline intervals (reference, from Adhoc Case A1b):**
- Interval 1: phase [0.1875, 0.25] — ~March 10 to April 1, centered on March equinox
- Interval 2: phase [0.625, 0.656] — ~August 16 to August 28, ~1 month before September equinox
- Interval 3: phase [0.875, 0.917] — ~November 16 to December 1, ~1 month before December solstice

For each of three mainshock catalogs (G-K, Reasenberg, A1b) at each of k=16, 24, 32:

1. Compute phase and bin counts as in Sub-analysis A
2. Identify elevated bins: bins whose observed count exceeds expected count by more than 1 standard deviation (`threshold = E + sqrt(E)`)
3. Merge adjacent elevated bins into contiguous phase intervals; record `[phase_start, phase_end]` for each merged interval
4. Classify each recovered interval against the A1b baseline:
   - "matches interval 1": recovered interval overlaps phase [0.1875, 0.25] by > 50% of its width
   - "matches interval 2": overlaps phase [0.625, 0.656] by > 50%
   - "matches interval 3": overlaps phase [0.875, 0.917] by > 50%
   - "new interval": does not overlap any A1b baseline interval
5. For each baseline interval, determine survival status across the three catalogs: "survives all", "survives some", "absent"

Compute phase coherence for each recovered interval:
- Mean resultant length R of phases within the interval: `R_interval = |sum(exp(2πi * phase[in_interval]))| / n_interval`

Output structure in results JSON under key `"sub_b"`:
```json
{
  "sub_b": {
    "a1b_baseline_intervals": [
      {"id": 1, "phase_start": 0.1875, "phase_end": 0.25, "calendar": "~Mar 10 – Apr 1"},
      {"id": 2, "phase_start": 0.625, "phase_end": 0.656, "calendar": "~Aug 16 – Aug 28"},
      {"id": 3, "phase_start": 0.875, "phase_end": 0.917, "calendar": "~Nov 16 – Dec 1"}
    ],
    "gk_mainshocks": {
      "k16": {"recovered_intervals": [...], "baseline_survival": {...}},
      "k24": {...},
      "k32": {...}
    },
    "reas_mainshocks": {...},
    "a1b_mainshocks": {...},
    "interval_survival_summary": {
      "interval_1": {"gk": "survives/absent", "reas": "...", "a1b": "..."},
      "interval_2": {...},
      "interval_3": {...}
    }
  }
}
```

---

## 4. Sub-analysis C — Aftershock population phase preference (aftershock catalogs)

In `src/case-a4-sub-c.py` (callable from main analysis):

For each of three aftershock catalogs (G-K, Reasenberg, A1b) at each of k=16, 24, 32:

1. Compute phase and bin counts as in Sub-analysis A
2. Compute chi-square, p-value, Cramér's V, Rayleigh R, and mean phase angle (same method as Sub-analysis A)
3. If p_chi2 < 0.05: identify which bins are elevated (same 1-SD threshold as Sub-analysis B); record the elevated phase intervals
4. Compare elevated phase intervals in aftershock population against: (a) A1b baseline intervals, (b) elevated intervals found in the corresponding mainshock catalog from Sub-analysis B

Classify aftershock phase preference result for each method:
- "no preference": p_chi2 >= 0.05 at all k
- "same intervals as mainshocks": elevated intervals overlap mainshock elevated intervals
- "different intervals from mainshocks": significant but non-overlapping intervals
- "significant, intervals match A1b baseline": significant and overlapping A1b baseline

Output structure in results JSON under key `"sub_c"`:
```json
{
  "sub_c": {
    "gk_aftershocks": {
      "n": 3327,
      "k16": {"chi2": float, "p_chi2": float, "cramer_v": float, "rayleigh_R": float, "p_rayleigh": float, "elevated_intervals": [...]},
      "k24": {...},
      "k32": {...},
      "classification": "no preference | same as mainshocks | different | matches A1b baseline"
    },
    "reas_aftershocks": {...},
    "a1b_aftershocks": {...}
  }
}
```

---

## 5. Visualizations

In `src/visualization-case-a4.py`:

**Figure 1 — Sub-analysis A bin distributions** (`output/case-a4-sub-a-binplot.png`):
- 4-panel plot (raw, G-K mainshocks, Reasenberg mainshocks, A1b mainshocks), all at k=24
- Each panel: horizontal bar chart of observed bin counts, with a dashed horizontal line at expected count (n/k)
- Bars colored steelblue; bars exceeding 1-SD threshold colored orange
- Annotate each panel with: n, χ², p-value, Cramér's V
- X-axis: phase fraction (0 to 1); label equinox positions (0.19, 0.69) and solstice positions (0.44, 0.94) with vertical dotted lines
- 300 DPI, titles, axis labels, legend

**Figure 2 — Sub-analysis B interval maps** (`output/case-a4-sub-b-intervals.png`):
- 3-row x 3-column grid: rows = catalogs (G-K, Reasenberg, A1b mainshocks); columns = bin counts (k=16, 24, 32)
- Each cell: phase axis (0–1), with A1b baseline intervals shown as shaded gray bands; recovered elevated intervals shown as colored bands (green if matching baseline, red if new)
- Title each cell with catalog name and k value
- 300 DPI

**Figure 3 — Sub-analysis C aftershock distributions** (`output/case-a4-sub-c-aftershock.png`):
- 3-panel plot (G-K, Reasenberg, A1b aftershocks), all at k=24
- Same format as Figure 1 panels
- Annotate with n, χ², p-value; add classification label
- 300 DPI

---

## 6. Test suite

In `tests/test-case-a4.py`:

- `test_catalog_counts`: assert loaded row counts match expected values (9210, 5883, 3327, 8265, 945, 7137, 2073); assert partition integrity (mainshocks + aftershocks = 9210 per method; zero usgs_id overlap)
- `test_phase_normalization`: create a synthetic 100-event catalog with known solar_secs values spanning 0–31,557,600; assert all computed phases are in [0, 1); assert phase for solar_secs=0 is 0.0; assert phase for solar_secs equal to year length is 0.0
- `test_chi_square_uniform`: generate a perfectly uniform bin distribution; assert chi2 ≈ 0.0 and p ≈ 1.0
- `test_chi_square_spike`: generate a distribution with a 3-SD spike in one bin; assert p < 0.01
- `test_rayleigh_uniform`: generate 1000 uniform random phases; assert Rayleigh p > 0.05
- `test_cramer_v_range`: assert Cramér's V for all catalogs at all k is in [0.0, 1.0]
- `test_raw_chi2_range`: assert the raw catalog chi2 at k=24 is > 20.0 (consistency with prior result); assert p < 0.01
- `test_suppression_direction`: assert that for G-K, chi2_mainshock < chi2_raw at all k (declustering reduces signal); assert A1b suppression < G-K suppression at all k (A1b removes fewer events)
- `test_interval_classification_logic`: create synthetic recovered intervals; assert interval overlapping [0.19, 0.25] is classified as "matches interval 1"; assert non-overlapping interval is classified as "new interval"
- `test_results_json_keys`: load `output/case-a4-results.json`; assert keys "sub_a", "sub_b", "sub_c" are present; assert "sub_a" contains "raw", "gk_mainshocks", "reas_mainshocks", "a1b_mainshocks"

---

## 7. Whitepaper

In `output/case-a4-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state the question (does the solar-phase signal survive declustering?), methods (three declustering methods, three sub-analyses), and key findings
2. **Data Source** — describe ISC-GEM catalog (n=9,210, M≥6.0, 1950–2021); note all six declustered files and their counts; note ISC-GEM vs ComCat distinction
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1 and `rules/data-handling.md`; describe the normalization formula)
   - 3.2 Declustering methods: describe G-K (cite Gardner & Knopoff 1974; note unverified table caveat from data-requirements.md), Reasenberg (parameters), A1b-informed window (83.2 km / 95.6 days; cite Adhoc A1b)
   - 3.3 Sub-analysis A: chi-square, Rayleigh, Cramér's V procedure
   - 3.4 Sub-analysis B: elevated-bin identification and interval matching procedure; list A1b baseline intervals
   - 3.5 Sub-analysis C: aftershock phase preference procedure
4. **Results**
   - 4.1 Sub-analysis A: table of chi2, p, Cramér's V, and suppression % for all four catalogs at k=24 (embed Figure 1)
   - 4.2 Sub-analysis B: interval survival table; embed Figure 2; state which baseline intervals survived under which methods
   - 4.3 Sub-analysis C: aftershock chi2, p, and classification for each method; embed Figure 3
5. **Cross-Topic Comparison** — compare ISC-GEM G-K suppression to L4/L5 ComCat baseline (60.2% suppression, χ²=45.61→18.13); note that values are not expected to match due to catalog differences; describe what the comparison implies
6. **Interpretation** — apply the outcome classification framework from the planning doc: which of the five sub-analysis B outcome scenarios occurred? Which of the three sub-analysis C outcome scenarios? State conclusions with appropriate objectivity
7. **Limitations** — G-K window table unverified; 52.9% GCMT match rate not applicable here but note the general data quality; phase normalization uses Julian year approximation (if applicable)
8. **References** — Bradley & Hubbard (2024), Gardner & Knopoff (1974), Park et al. (2021), Reasenberg (1985), Adhoc A1, Adhoc A1b, L4/L5 topic reference

---

## 8. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case A4: Declustering Sensitivity Analysis
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [chi2 raw, chi2 G-K mainshocks, chi2 Reasenberg mainshocks, chi2 A1b mainshocks at k=24; interval survival summary; aftershock phase preference classification per method]
  ```
- Update Case A4 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
