# Case B1: Hemisphere Stratification — Phase Symmetry Test

**Status:** Pending

**Intent statement:**
Case B1 splits the ISC-GEM catalog by hemisphere and independently computes solar-phase bin distributions for Northern Hemisphere (NH, latitude > 0°) and Southern Hemisphere (SH, latitude < 0°) events. The original bimodal equinox framing — testing whether NH and SH peaks are in-phase (geometric forcing) or anti-phased by 6 months (hydrological loading) — has been replaced by a richer four-question test based on the three-interval structure identified in Adhoc A1b. The test asks: (1) do all three elevated intervals appear in both hemispheres, (2) is interval 1 (March equinox, phase 0.1875–0.25) symmetric in both hemispheres, (3) are intervals 2 (phase 0.625–0.656) and 3 (phase 0.875–0.917) hemisphere-specific, and (4) are any intervals 0.5 cycles apart between hemispheres. A symmetric result (all three intervals in both hemispheres at the same phase) is the strongest evidence for geometric solar forcing. A hemisphere-specific result for intervals 2 and 3 is consistent with hydrological or seasonal forcing at those phases. No published study has applied this specific three-interval symmetry test to a global M≥6.0 catalog.

**Relationship to prior topics:**
The hemisphere split variable (`latitude`) was added to the ISC-GEM catalog in Topic L4 (Declustering), where latitude and depth columns were confirmed present. Adhoc A1b identified the three elevated phase intervals (0.1875–0.25, 0.625–0.656, 0.875–0.917) that define the specific questions this case tests. B1 draws directly on A4's declustering results: if intervals 2 and 3 were found to be partially aftershock-driven in A4, then the hemisphere-specific analysis of those intervals requires careful interpretation, since aftershock sequences are spatially clustered and may not be evenly distributed between hemispheres.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/global-sets/iscgem_global_events.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `latitude`, `longitude`, `depth` |

Hemisphere split: `NH = df[df["latitude"] > 0]`; `SH = df[df["latitude"] < 0]`; events at latitude=0 are excluded (equatorial events, expected n≈0). Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year [confirm before running: consistent with A4/B6/A1 — verify uniform Julian constant vs per-year values].

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-b1-analysis.py` — main analysis: hemisphere split, per-hemisphere bin distributions, four symmetry tests, interval survival classification, writes results JSON
- `src/visualization-case-b1.py` — generates all PNG figures
- `tests/test-case-b1.py` — test suite
- `output/case-b1-results.json` — hemisphere sizes, per-hemisphere statistics, four symmetry test results, interval classification
- `output/case-b1-whitepaper.md` — methodology, results, interpretation
- `output/case-b1-binplot-nh.png` — NH phase bin distribution at k=24
- `output/case-b1-binplot-sh.png` — SH phase bin distribution at k=24
- `output/case-b1-interval-comparison.png` — side-by-side NH vs SH elevated interval maps across k=16, 24, 32

---

## 1. Environment and data loading

In `src/case-b1-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define `RAW_PATH = BASE_DIR.parent / "data" / "global-sets" / "iscgem_global_events.csv"`
- Load catalog; assert n=9210; log row count
- Create hemisphere subsets:
  - `df_nh = df[df["latitude"] > 0].copy()` — log n_nh
  - `df_sh = df[df["latitude"] < 0].copy()` — log n_sh
  - `n_equatorial = len(df[df["latitude"] == 0])` — log and record; typically 0 or near 0
  - Assert n_nh + n_sh + n_equatorial == 9210
- Compute phase for each subset using phase normalization (Julian year constant)

---

## 2. Per-hemisphere bin statistics

For each of two subsets (NH, SH) at each of three bin counts (k=16, 24, 32):

1. Compute observed bin counts `O` (array of length k)
2. Expected counts: `E = np.full(k, n_hemisphere / k)`
3. Chi-square: `chi2, p_chi2 = scipy.stats.chisquare(O, E)`
4. Cramér's V: `V = sqrt(chi2 / (n * (k-1)))`
5. Rayleigh R and p-value (same formula as A4/B6)
6. Mean phase angle (same formula as A4/B6)
7. Elevated bins: bins where `O > E + sqrt(E)` (1-SD threshold)
8. Elevated phase intervals: merge adjacent elevated bins into contiguous intervals; record `[phase_start, phase_end]`

Record all per-hemisphere, per-k statistics in results JSON under key `"hemisphere_stats"`:
```json
{
  "hemisphere_stats": {
    "n_nh": int, "n_sh": int, "n_equatorial": int,
    "nh": {
      "k16": {"chi2": float, "p_chi2": float, "cramer_v": float, "rayleigh_R": float, "p_rayleigh": float,
               "mean_phase": float, "elevated_intervals": [...]},
      "k24": {...}, "k32": {...}
    },
    "sh": {"k16": {...}, "k24": {...}, "k32": {...}}
  }
}
```

---

## 3. Four symmetry tests

Perform the four symmetry tests defined in the planning doc. Use results from k=24 as primary; note if k=16 or k=32 differs substantially.

**Test 1 — Global symmetry:** Do all three A1b baseline intervals appear in both NH and SH?
- For each interval (1, 2, 3): check whether the NH elevated intervals include a match (overlap > 50%) and whether the SH elevated intervals include a match
- Record per-interval: `{"in_nh": bool, "in_sh": bool, "symmetric": bool}` (symmetric = appears in both)
- Classify: "fully symmetric" if all three appear in both; "partially symmetric" if some appear in both; "asymmetric" if none appear in both

**Test 2 — Interval 1 symmetry (March equinox, phase 0.1875–0.25):**
- Specifically test whether the March-equinox interval is present in both hemispheres at the same phase position (within 1 bin width at k=24, i.e., ±0.042 phase units)
- Record: `{"interval_1_nh": bool, "interval_1_sh": bool, "phase_offset": float | null}` (phase_offset = NH mean interval phase − SH mean interval phase)

**Test 3 — Intervals 2 and 3 hemisphere-specificity:**
- For interval 2 (phase 0.625–0.656): does it appear in NH only, SH only, both, or neither?
- For interval 3 (phase 0.875–0.917): does it appear in NH only, SH only, both, or neither?
- Record classification for each interval

**Test 4 — Half-cycle offset test (hydrological hypothesis):**
- For each elevated interval in NH, compute the expected SH counterpart phase: `counterpart_phase = (nh_interval_center + 0.5) % 1.0`
- Check whether any SH elevated interval is within 1 bin width (0.042 at k=24) of that counterpart phase
- Record: `{"any_half_cycle_offset_found": bool, "details": [...]}`

Output symmetry test results under key `"symmetry_tests"` in results JSON:
```json
{
  "symmetry_tests": {
    "test_1_global": {"interval_1": {...}, "interval_2": {...}, "interval_3": {...}, "classification": "..."},
    "test_2_interval_1": {"interval_1_nh": bool, "interval_1_sh": bool, "phase_offset": float},
    "test_3_interval_23_specificity": {
      "interval_2": {"hemisphere": "nh_only | sh_only | both | neither"},
      "interval_3": {"hemisphere": "nh_only | sh_only | both | neither"}
    },
    "test_4_half_cycle_offset": {"any_half_cycle_offset_found": bool, "details": [...]}
  }
}
```

---

## 4. Prediction matching

Evaluate which of the three mechanistic predictions is best supported:

- **Geometric hypothesis** (predicted by planning doc): all three intervals in both hemispheres at the same phase
  - Evaluate: all Test 1 intervals symmetric AND Test 4 no half-cycle offset found
- **Hydrological hypothesis**: NH and SH peaks offset by ~0.5 cycle; intervals 2 and 3 hemisphere-specific
  - Evaluate: Test 4 half-cycle offset found AND Test 3 shows interval 2 or 3 hemisphere-specific
- **Mixed hypothesis**: Interval 1 symmetric (genuine), intervals 2 and 3 hemisphere-specific or absent
  - Evaluate: Test 2 interval 1 symmetric AND Test 3 interval 2 and/or 3 hemisphere-specific

Record prediction scores as a summary:
```json
{
  "prediction_support": {
    "geometric": "supported | partially supported | not supported",
    "hydrological": "supported | partially supported | not supported",
    "mixed": "supported | partially supported | not supported",
    "primary_conclusion": "geometric | hydrological | mixed | ambiguous"
  }
}
```

---

## 5. Visualizations

In `src/visualization-case-b1.py`:

**Figure 1 — NH bin distribution** (`output/case-b1-binplot-nh.png`):
- Horizontal bar chart at k=24; steelblue bars; dashed expected-count line; 1-SD threshold dashed line
- Mark A1b baseline intervals (phases 0.1875–0.25, 0.625–0.656, 0.875–0.917) with gray shaded bands
- Annotate: n_nh, χ², p-value, Cramér's V; title "Northern Hemisphere (n=...)"
- 300 DPI

**Figure 2 — SH bin distribution** (`output/case-b1-binplot-sh.png`):
- Same format as Figure 1; title "Southern Hemisphere (n=...)"
- 300 DPI

**Figure 3 — Side-by-side interval comparison** (`output/case-b1-interval-comparison.png`):
- 2-column x 3-row grid: columns = NH, SH; rows = k=16, 24, 32
- Each cell: phase axis 0–1; elevated bins shown as steelblue vertical bars; A1b baseline intervals as gray shaded bands
- Annotate each cell with hemisphere, k, and whether each A1b interval was recovered
- 300 DPI

---

## 6. Test suite

In `tests/test-case-b1.py`:

- `test_hemisphere_split`: assert n_nh + n_sh + n_equatorial == 9210; assert n_nh > 0 and n_sh > 0
- `test_phase_range`: assert all phases in [0.0, 1.0)
- `test_chi_square_both_hemispheres`: assert chi2 and p_chi2 are finite floats for both hemispheres at k=16, 24, 32
- `test_cramer_v_range`: assert Cramér's V for both hemispheres at all k is in [0.0, 1.0]
- `test_interval_classification_keys`: load results JSON; assert `"symmetry_tests"` key present with sub-keys `"test_1_global"`, `"test_2_interval_1"`, `"test_3_interval_23_specificity"`, `"test_4_half_cycle_offset"`
- `test_prediction_support_valid_values`: assert `"prediction_support"` values are one of "supported", "partially supported", "not supported"; assert `"primary_conclusion"` is one of "geometric", "hydrological", "mixed", "ambiguous"
- `test_elevated_intervals_phase_range`: assert all recovered elevated interval phase_start and phase_end values in [0.0, 1.0]
- `test_half_cycle_offset_logic`: create synthetic NH interval at center 0.22; assert expected SH counterpart is computed as 0.72 (within tolerance)
- `test_all_k_computed`: assert results JSON contains k16, k24, k32 entries for both "nh" and "sh" in `"hemisphere_stats"`

---

## 7. Whitepaper

In `output/case-b1-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (are elevated solar-phase intervals symmetric between hemispheres?), three-interval framework from A1b, four symmetry tests, primary finding, and mechanistic implication
2. **Data Source** — ISC-GEM catalog split: NH (n=...), SH (n=...); note hemisphere sizes and any equatorial events excluded
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1 and `rules/data-handling.md`)
   - 3.2 A1b three-interval baseline: describe the three intervals (phase, calendar dates) that define the symmetry questions
   - 3.3 Four symmetry tests: describe each test in the order in which they are evaluated
   - 3.4 Mechanistic prediction framework: describe the geometric, hydrological, and mixed predictions
4. **Results**
   - 4.1 NH distribution: embed Figure 1; state chi2, p, Cramér's V; identify recovered intervals
   - 4.2 SH distribution: embed Figure 2; same
   - 4.3 Symmetry test outcomes: embed Figure 3; tabulate results of all four tests
   - 4.4 Prediction matching: state which predictions are supported/not supported; state primary conclusion
5. **Cross-Topic Comparison** — compare to the prior bimodal equinox framing (Case 3A); compare to Adhoc A1b three-interval characterization; note which intervals survived which declustering methods in A4 and how that informs the hemisphere test interpretation
6. **Interpretation** — state mechanistic conclusion with appropriate objectivity; note that A4 declustering results should be considered when interpreting hemisphere-specific intervals 2 and 3
7. **Limitations** — NH events substantially outnumber SH events (seismicity concentrated in NH due to continental geometry); SH sample size limits statistical power; phase normalization uses Julian year approximation
8. **References** — Adhoc A1b, Colledge et al. (2025) [if relevant], Case 3A, B6, A4

---

## 8. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case B1: Hemisphere Stratification — Phase Symmetry Test
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [n_nh, n_sh; chi2/p for each hemisphere at k=24; four symmetry test outcomes; primary conclusion (geometric/hydrological/mixed/ambiguous)]
  ```
- Update Case B1 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
