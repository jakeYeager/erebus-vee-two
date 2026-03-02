# Case B5: Solar Declination Rate-of-Change vs. Position Test

**Status:** Pending

**Intent statement:**
Case B5 re-encodes ISC-GEM events using three distinct solar geometric variables — solar declination angle, rate of change of solar declination, and Earth-Sun distance — and computes bin distributions for each independently, testing which variable produces the most coherent and statistically significant clustering. Case 3A encoded events by their position in the solar annual cycle (`solar_secs`), which conflates these distinct physical variables. The declination rate peaks at equinoxes, the Earth-Sun distance has a single minimum (January perihelion) and maximum (July aphelion), and the declination angle has zero-crossings at equinoxes. However, Adhoc A1b's three-interval structure (March equinox phase 0.19, mid-August phase 0.64, late November phase 0.90) does not cleanly match any single solar variable: only interval 1 aligns with declination rate maximum at the March equinox; interval 2 (~mid-August) falls between perihelion and autumnal equinox, near the Earth-Sun distance midpoint; interval 3 (~late November) is near the December solstice, where declination rate is at a minimum. This case is therefore more likely to produce a negative result for all three variables individually, which is itself informative: it would suggest the mechanism is more complex than any single solar parameter, or that one or more of the three A1b intervals reflects sequence contamination rather than a physical forcing response. This case should be run after A4 (declustering) and B1 (hemisphere symmetry) to inform interpretation.

**Relationship to prior topics:**
The solar geometry enrichment file (`solar_geometry_global.csv`) was produced for this topic (data-requirements.md REQ-3) by extending the existing Skyfield 1.54 + JPL DE421 ephemeris pipeline from topic L3. The file contains `solar_declination`, `declination_rate`, and `earth_sun_distance` as new columns, plus all base catalog columns. Topic L3 established the ephemeris computation framework. Case 3A's `solar_secs` variable (the cyclic annual position) was the basis of the original finding; B5 tests whether any of the three new geometric variables provides a better physical explanation for that signal.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| Solar geometry enriched catalog | `data/iscgem/celestial-geometry/solar_geometry_global.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `latitude`, `longitude`, `depth`, `solar_declination`, `declination_rate`, `earth_sun_distance` |

Column units and ranges (from data-requirements.md REQ-3):
- `solar_declination`: degrees, range −23.5 to +23.5
- `declination_rate`: degrees/day, range ~−0.40 to +0.40
- `earth_sun_distance`: AU, range ~0.983 to ~1.017

Phase normalization for `solar_secs`: as in all prior cases (Julian year constant).

Phase binning for declination and distance: these are continuous, non-cyclic variables. They must be discretized differently than `solar_secs`:
- `solar_declination`: range [−23.5, +23.5], total range 47°; bin as: `declination_bin = floor((solar_declination - (-23.5)) / 47.0 * k)`, clip to [0, k-1]
- `declination_rate`: range [−0.40, +0.40], total range 0.80 deg/day; bin as: `rate_bin = floor((declination_rate - (-0.40)) / 0.80 * k)`, clip to [0, k-1]
- `earth_sun_distance`: range [0.983, 1.017], total range 0.034 AU; bin as: `dist_bin = floor((earth_sun_distance - 0.983) / 0.034 * k)`, clip to [0, k-1]

Note: these are non-cyclic bins. The chi-square test still applies (expected uniform distribution over k bins under the null). The comparison to `solar_secs` (cyclic) is a comparison of which encoding produces stronger non-uniformity.

[confirm before running: verify that the full range of solar_declination, declination_rate, and earth_sun_distance in the actual delivered file matches the REQ-3 specified ranges (−23.5 to +23.5, ~−0.40 to +0.40, ~0.983 to ~1.017); if actual ranges differ, update bin computation accordingly]

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-b5-analysis.py` — main analysis: load solar geometry file, compute bins for all four variables, chi-square/Rayleigh per variable, variable comparison ranking, writes results JSON
- `src/visualization-case-b5.py` — generates all PNG figures
- `tests/test-case-b5.py` — test suite
- `output/case-b5-results.json` — per-variable statistics, variable ranking, A1b interval alignment analysis
- `output/case-b5-whitepaper.md` — methodology, results, variable comparison
- `output/case-b5-binplots.png` — 4-panel bin distributions (solar_secs, declination, rate, distance) at k=24
- `output/case-b5-variable-ranking.png` — bar chart of Cramér's V by variable
- `output/case-b5-a1b-alignment.png` — overlay showing which physical variable values correspond to each A1b interval

---

## 1. Environment and data loading

In `src/case-b5-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define `SOLAR_GEO_PATH = BASE_DIR.parent / "data" / "iscgem" / "celestial-geometry" / "solar_geometry_global.csv"`
- Load file; assert n=9210; log row count
- Assert required columns present: `solar_secs`, `solar_declination`, `declination_rate`, `earth_sun_distance`
- Log actual min/max of each variable; compare to expected ranges; warn if outside expected bounds
- Compute `solar_phase = (solar_secs / 31_557_600.0) % 1.0` (cyclic; Julian year)
- Compute non-cyclic bins for each variable at k=24 and k=32:
  - Use actual min/max from data for bin computation if they differ from REQ-3 spec ranges; log which range is used
  - Define `bin_solar_phase_k24 = floor(solar_phase * 24).clip(0, 23)`
  - Define `bin_declination_k24 = floor((solar_declination - actual_min_dec) / actual_range_dec * 24).clip(0, 23)`
  - Define `bin_rate_k24 = floor((declination_rate - actual_min_rate) / actual_range_rate * 24).clip(0, 23)`
  - Define `bin_distance_k24 = floor((earth_sun_distance - actual_min_dist) / actual_range_dist * 24).clip(0, 23)`
  - Repeat for k=32

---

## 2. Per-variable bin statistics

For each of four variables (`solar_phase`, `solar_declination`, `declination_rate`, `earth_sun_distance`) at each of k=16, 24, 32:

1. Compute observed bin counts, expected counts (n/k), chi-square, p-value, Cramér's V, Rayleigh R (Rayleigh applies to cyclic `solar_phase` only; for non-cyclic variables, substitute the Kolmogorov-Smirnov test against a uniform distribution), mean phase or mean bin fraction
2. For `solar_phase` (cyclic): compute Rayleigh R and p-value; compute mean phase angle
3. For non-cyclic variables: compute `ks_stat, ks_p = scipy.stats.kstest(normalized_values, 'uniform')` where `normalized_values = (variable - min) / range`; record KS stat and p-value

Output structure under key `"variable_stats"` in results JSON:
```json
{
  "variable_stats": {
    "solar_phase": {
      "k16": {"chi2": float, "p_chi2": float, "cramer_v": float, "rayleigh_R": float, "p_rayleigh": float, "bin_counts": [...]},
      "k24": {...}, "k32": {...}
    },
    "solar_declination": {
      "k16": {"chi2": float, "p_chi2": float, "cramer_v": float, "ks_stat": float, "ks_p": float, "bin_counts": [...]},
      "k24": {...}, "k32": {...}
    },
    "declination_rate": {"k16": {...}, "k24": {...}, "k32": {...}},
    "earth_sun_distance": {"k16": {...}, "k24": {...}, "k32": {...}}
  }
}
```

---

## 3. Variable ranking and A1b alignment

Using results at k=24:

1. Rank variables by Cramér's V (descending): record `cramer_v_ranking = [("solar_phase", V), ("declination_rate", V), ...]`
2. Identify which variable produces the lowest p_chi2 (most significant); record `most_significant_variable`
3. Compute physical values at the three A1b interval phase centers:
   - A1b interval 1 center: phase 0.22 → calendar date ~March 22; compute expected `solar_declination` ≈ 0° (equinox), `declination_rate` ≈ +0.40 deg/day (maximum), `earth_sun_distance` ≈ 1.00 AU (near mean)
   - A1b interval 2 center: phase 0.64 → calendar date ~August 22; compute expected `solar_declination` ≈ +12° (post-summer solstice), `declination_rate` ≈ −0.25 deg/day (declining toward autumn), `earth_sun_distance` ≈ 1.012 AU (near maximum)
   - A1b interval 3 center: phase 0.90 → calendar date ~November 24; compute expected `solar_declination` ≈ −21° (approaching winter solstice), `declination_rate` ≈ −0.12 deg/day (near minimum), `earth_sun_distance` ≈ 0.988 AU (approaching perihelion)
   - These are theoretical values; record for comparison with actual bin distributions

4. Record which A1b intervals fall in elevated bins for each variable at k=24; classify each interval's alignment with each variable

Output under key `"variable_ranking"` and `"a1b_alignment"`:
```json
{
  "variable_ranking": [
    {"variable": "solar_phase", "cramer_v_k24": float, "p_chi2_k24": float},
    ...
  ],
  "most_significant_variable": "solar_phase | solar_declination | declination_rate | earth_sun_distance",
  "a1b_alignment": {
    "interval_1_phase_0.22": {
      "solar_declination_expected": 0.0,
      "declination_rate_expected": 0.40,
      "earth_sun_distance_expected": 1.00,
      "variable_with_elevated_bin": ["declination_rate"]
    },
    "interval_2_phase_0.64": {...},
    "interval_3_phase_0.90": {...}
  }
}
```

---

## 4. Visualizations

In `src/visualization-case-b5.py`:

**Figure 1 — Four-panel bin distributions** (`output/case-b5-binplots.png`):
- 2×2 panel layout: solar_phase (top-left), solar_declination (top-right), declination_rate (bottom-left), earth_sun_distance (bottom-right); all at k=24
- Each panel: horizontal steelblue bar chart; dashed expected-count line; 1-SD threshold dashed line
- For solar_phase panel: add A1b baseline interval gray bands and equinox/solstice vertical lines
- For non-cyclic panels: x-axis labeled with physical units (degrees, deg/day, AU) at actual min/max ticks
- Annotate each panel with variable name, χ², p-value, Cramér's V
- 300 DPI

**Figure 2 — Variable ranking bar chart** (`output/case-b5-variable-ranking.png`):
- Grouped horizontal bar chart: y-axis = variable names; x-axis = Cramér's V at k=24
- Steelblue bars; mark significance (p < 0.05) with asterisk; mark non-significant with 'ns'
- Title "Cramér's V by Solar Variable at k=24"
- 300 DPI

**Figure 3 — A1b interval alignment** (`output/case-b5-a1b-alignment.png`):
- 3-row × 3-column grid: rows = A1b intervals 1, 2, 3; columns = solar_declination, declination_rate, earth_sun_distance
- Each cell: show the full distribution of that variable (histogram) with a vertical red line at the expected value for that A1b interval phase center
- Title each row with interval phase and calendar date; title each column with variable name and units
- 300 DPI

---

## 5. Test suite

In `tests/test-case-b5.py`:

- `test_solar_geo_load`: assert n=9210; assert all four required columns present
- `test_column_ranges_approx`: assert solar_declination min/max within ±2° of [−23.5, +23.5]; assert declination_rate within ±0.05 deg/day of [−0.40, +0.40]; assert earth_sun_distance within ±0.002 AU of [0.983, 1.017]
- `test_phase_range`: assert solar_phase in [0.0, 1.0)
- `test_non_cyclic_bins_in_range`: assert all bin_declination, bin_rate, bin_distance values in [0, k-1] for k=24
- `test_chi_square_all_variables`: assert chi2 and p_chi2 are finite for all four variables at k=16, 24, 32
- `test_cramer_v_range`: assert Cramér's V in [0.0, 1.0] for all variables
- `test_rayleigh_solar_phase_only`: assert results JSON has rayleigh_R field for solar_phase but not for other variables
- `test_ks_non_cyclic_only`: assert results JSON has ks_stat field for solar_declination, declination_rate, and earth_sun_distance but not solar_phase
- `test_variable_ranking_complete`: load results JSON; assert variable_ranking has exactly 4 entries covering all four variable names
- `test_most_significant_valid`: assert most_significant_variable is one of the four variable names
- `test_a1b_alignment_keys`: assert a1b_alignment has keys for all three intervals

---

## 6. Whitepaper

In `output/case-b5-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (which solar variable best explains the seismic clustering?), four variables tested, ranking result, A1b interval alignment finding, and interpretation of negative/positive results
2. **Data Source** — solar geometry enriched catalog (n=9,210); describe the three new columns and how they were computed (Skyfield 1.54, JPL DE421, cite data-requirements.md REQ-3); describe `solar_secs` as the base cyclic variable
3. **Methodology**
   - 3.1 Phase normalization for `solar_secs` (cite Adhoc A1, `rules/data-handling.md`)
   - 3.2 Non-cyclic binning for declination, rate, distance: describe range normalization and bin-assignment formula; note difference from cyclic phase
   - 3.3 Chi-square and Cramér's V for all four variables; Rayleigh for cyclic only; KS test for non-cyclic
   - 3.4 Variable ranking by Cramér's V
   - 3.5 A1b interval physical-value alignment: describe the cross-reference of phase positions to physical solar variable values
4. **Results**
   - 4.1 Bin distributions: embed Figure 1; report chi2, p, Cramér's V for all four variables at k=24
   - 4.2 Variable ranking: embed Figure 2; state most significant variable; state whether solar_phase or any geometric variable dominates
   - 4.3 A1b alignment: embed Figure 3; for each A1b interval, describe which solar variable (if any) shows elevated events at the corresponding physical value
5. **Cross-Topic Comparison** — compare to Case 3A `solar_secs`-based finding; compare to Adhoc A1b three-interval structure; note that if declination_rate is not the dominant variable, the bimodal equinox interpretation from Case 3A is not supported by the variable decomposition
6. **Interpretation** — state which variable (if any) best explains the signal; discuss the three-interval challenge for single-variable explanations; maintain objectivity; note that a fully negative result across all three geometric variables is a meaningful finding
7. **Limitations** — non-cyclic binning is not directly comparable to cyclic phase; the physical variable ranges vary across years but bin thresholds are computed from overall min/max; case is run on raw catalog (declustered catalogs not used for this variable test); sequenced after A4 and B1 but does not use those results directly
8. **References** — Ader & Avouac (2013), Case 3A, Adhoc A1b, A4, B1, data-requirements.md REQ-3

---

## 7. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case B5: Solar Declination Rate-of-Change vs. Position Test
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [Cramér's V ranking of four variables at k=24; most significant variable; A1b interval alignment summary; interpretation of negative/positive variable-specific results]
  ```
- Update Case B5 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
