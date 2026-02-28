# Case B3: Tectonic Regime Stratification

**Status:** Pending

**Intent statement:**
Case B3 classifies ISC-GEM events by focal mechanism type (thrust, normal, strike-slip) using the GCMT proximity join file, then independently computes solar-phase bin statistics for each tectonic class. The original motivation — different fault geometries respond differently to Coulomb stress changes from surface loading, so thrust and normal faults should show approximately anti-phased patterns if loading drives the signal — is strengthened by Métivier et al. (2009), who found slightly greater tidal triggering in normal and strike-slip regions than thrust regions at global scale. This provides a reference pattern: if the annual solar signal matches the Métivier tidal pattern (normal > strike-slip > thrust), a similar stress-change geometry is implied at annual timescales even if the mechanism differs. If the pattern is reversed or absent for thrust faults but present for normals, the loading hypothesis is supported. The GCMT focal mechanism join covers 4,874 of 9,210 events (52.9% match rate); mechanism is classified from the `rake` column using standard quadrant boundaries. Pre-1976 ISC-GEM events (approximately 26% of catalog) have no GCMT match. The 4,336 unmatched events are excluded from tectonic stratification but analyzed separately as a completeness check.

**Relationship to prior topics:**
The focal mechanism join file was produced for this topic (data-requirements.md REQ-4) using GCMT proximity matching at ±60s, 50 km, ±0.3 magnitude — the same tolerances used in Adhoc Case A0b for cross-catalog matching. Topic L3 established the ISC-GEM catalog schema. The PB2002 classification (also available from REQ-5) can serve as a coarse proxy for the 47% unmatched events, but the planning doc explicitly states PB2002 is not a substitute for actual focal mechanism data. The tectonic stratification is the most data-dependent case in the topic due to the GCMT coverage limitation.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/global-sets/iscgem_global_events.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `latitude`, `longitude`, `depth` |
| GCMT focal mechanism join | `data/iscgem/focal-mechanism/focal_join_global.csv` | 9,210 | `usgs_id`, `gcmt_id`, `mechanism`, `rake`, `strike`, `dip`, `scalar_moment`, `centroid_depth`, `match_confidence` |

GCMT join schema: `mechanism` column already classified (`thrust`, `normal`, `strike_slip`, `oblique`, or null). `match_confidence` is `proximity` for matched events, null for unmatched. The file is an enriched catalog (all base columns present as prefix).

Focal mechanism classification from `rake` (if `mechanism` column is null or needs re-derivation):
- `thrust`: rake in [45°, 135°]
- `normal`: rake in [−135°, −45°]
- `strike_slip`: rake in (−45°, 45°] or (135°, 180°] / [−180°, −135°)
- `oblique`: remaining cases

Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year [confirm before running: consistent with prior cases — verify uniform Julian constant vs per-year values].

[confirm before running: verify that `focal_join_global.csv` includes the `solar_secs` column (expected as part of base catalog prefix) or whether the raw catalog must be joined separately to obtain `solar_secs` for the 4,874 matched events]

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-b3-analysis.py` — main analysis: load focal join, classify mechanisms, per-class bin statistics, Métivier comparison, writes results JSON
- `src/visualization-case-b3.py` — generates all PNG figures
- `tests/test-case-b3.py` — test suite
- `output/case-b3-results.json` — match coverage stats, per-mechanism statistics, pattern comparison
- `output/case-b3-whitepaper.md` — methodology, results, cross-topic comparison
- `output/case-b3-binplots.png` — 3-panel bin distributions (thrust, normal, strike-slip) at k=24
- `output/case-b3-mechanism-comparison.png` — Cramér's V and Rayleigh R by mechanism type (bar chart)
- `output/case-b3-coverage-map.png` — global map of events colored by mechanism class (matched events) and gray (unmatched)

---

## 1. Environment and data loading

In `src/case-b3-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define paths using BASE_DIR (see Data context block)
- Load raw catalog (n=9210) and focal join file (n=9210)
- If `focal_join_global.csv` contains `solar_secs` (base catalog prefix): use directly; otherwise join raw catalog to focal join on `usgs_id` to obtain `solar_secs`
- Classify matched events by mechanism:
  - Primary: use `mechanism` column if not null and not "oblique"
  - Fallback: if `mechanism` is null but `rake` is not null, re-classify from rake using the quadrant boundaries above; log count of events re-classified
  - Exclude unmatched (`match_confidence` is null) and oblique events from per-class analysis; log counts
- Log classification summary:
  - `n_thrust`, `n_normal`, `n_strike_slip`, `n_oblique`, `n_unmatched`, `n_total_matched`
  - Assert `n_total_matched + n_unmatched == 9210`
- Compute phase for each class subset using Julian year constant

---

## 2. Per-mechanism bin statistics

For each of three mechanism classes (thrust, normal, strike-slip) at each of k=16, 24, 32:

1. Compute observed bin counts, expected counts, chi-square, p-value, Cramér's V, Rayleigh R, Rayleigh p, mean phase (same procedure as A4)
2. Elevated bins at 1-SD threshold; elevated phase intervals (same as A4)
3. Compare elevated intervals to A1b baseline (phases 0.1875–0.25, 0.625–0.656, 0.875–0.917)

Output structure under key `"mechanism_stats"` in results JSON:
```json
{
  "mechanism_stats": {
    "coverage": {
      "n_thrust": int, "n_normal": int, "n_strike_slip": int,
      "n_oblique": int, "n_unmatched": int, "n_total": 9210
    },
    "thrust": {
      "n": int,
      "k16": {"chi2": float, "p_chi2": float, "cramer_v": float, "rayleigh_R": float,
               "p_rayleigh": float, "mean_phase": float, "bin_counts": [...], "elevated_intervals": [...]},
      "k24": {...}, "k32": {...}
    },
    "normal": {...},
    "strike_slip": {...}
  }
}
```

---

## 3. Métivier pattern comparison and prediction evaluation

Using results at k=24:

1. Rank the three mechanism classes by Cramér's V: record `cramer_v_rank_order` (e.g., ["normal", "strike_slip", "thrust"] for Métivier pattern)
2. Classify against Métivier (2009) tidal pattern: if rank order matches "normal >= strike_slip >= thrust", flag `matches_metivier_pattern = True`
3. Check anti-phase relationship for loading hypothesis:
   - Loading predicts thrust and normal are anti-phased (mean phase offset ~0.5)
   - Compute `phase_offset_thrust_normal = |mean_phase_thrust - mean_phase_normal|`; if > 0.5, use `1 - phase_offset`
   - Classify: `"anti_phased"` if offset in [0.4, 0.6]; `"in_phase"` if offset < 0.2; `"other"` otherwise
4. Check all-classes equinox excess for geometric hypothesis:
   - All three classes show elevated intervals overlapping A1b interval 1 (March equinox, phase 0.1875–0.25)
   - Record `geometric_equinox_in_all = True/False`

Prediction matching:
- **Loading hypothesis**: thrust and normal anti-phased, thrust shows solstice deficit
- **Geometric hypothesis**: all three classes show similar phase preference near equinox
- **Métivier tidal pattern**: normal > strike-slip > thrust by Cramér's V

Output under key `"prediction_evaluation"`:
```json
{
  "prediction_evaluation": {
    "cramer_v_rank_order": [...],
    "matches_metivier_pattern": bool,
    "phase_offset_thrust_normal": float,
    "thrust_normal_phase_relationship": "anti_phased | in_phase | other",
    "geometric_equinox_in_all": bool,
    "loading_supported": bool,
    "geometric_supported": bool,
    "metivier_pattern_supported": bool,
    "primary_conclusion": "loading | geometric | metivier_consistent | ambiguous"
  }
}
```

---

## 4. Unmatched event sensitivity check

Compute chi-square and Cramér's V for the unmatched event subset (n=4,336):
- If unmatched events show similar or weaker signal than the full catalog, the GCMT subsample is representative
- If unmatched events show substantially different signal, the 47% coverage introduces selection bias

Record in results JSON under key `"unmatched_check"`:
```json
{
  "unmatched_check": {
    "n_unmatched": 4336,
    "k24_chi2": float, "k24_p": float, "k24_cramer_v": float,
    "similar_to_full_catalog": bool
  }
}
```

---

## 5. Visualizations

In `src/visualization-case-b3.py`:

**Figure 1 — Three-panel bin distributions** (`output/case-b3-binplots.png`):
- 3-panel layout (1×3 horizontal): thrust, normal, strike-slip at k=24
- Each panel: horizontal steelblue bar chart; dashed expected-count line; 1-SD threshold; A1b baseline interval gray bands
- Annotate: mechanism label, n, χ², p-value, Cramér's V
- 300 DPI

**Figure 2 — Mechanism comparison bar chart** (`output/case-b3-mechanism-comparison.png`):
- Dual-axis grouped bar chart: x-axis = mechanism (thrust, normal, strike-slip); left y-axis = Cramér's V (steelblue bars); right y-axis = Rayleigh R (orange bars)
- Mark significance (p < 0.05) with asterisk above bars
- Annotate whether Métivier pattern (normal ≥ strike-slip ≥ thrust) is matched
- 300 DPI

**Figure 3 — Coverage map** (`output/case-b3-coverage-map.png`):
- Global event scatter map: thrust=red, normal=blue, strike-slip=green, unmatched=gray
- Point size proportional to magnitude (as in B2 global maps)
- Legend with mechanism colors and n per class; title "ISC-GEM Events by Focal Mechanism (GCMT Join)"
- Use cartopy/basemap if available [confirm before running: same environment check as B2]; otherwise lat/lon scatter
- 300 DPI

---

## 6. Test suite

In `tests/test-case-b3.py`:

- `test_focal_join_load`: assert focal join n=9210
- `test_mechanism_coverage`: assert n_thrust + n_normal + n_strike_slip + n_oblique + n_unmatched == 9210
- `test_matched_count_range`: assert n_total_matched is between 4800 and 5000 (consistent with 52.9% reported match rate ≈ 4874 ± tolerance)
- `test_phase_range`: assert all phases in [0.0, 1.0)
- `test_chi_square_all_mechanisms`: assert chi2 and p_chi2 are finite for all three mechanisms at k=16, 24, 32
- `test_cramer_v_range`: assert Cramér's V in [0.0, 1.0] for all mechanisms
- `test_mechanism_class_labels_valid`: assert `cramer_v_rank_order` contains exactly the three labels "thrust", "normal", "strike_slip" in some order
- `test_rake_classification_logic`: for rake=90, assert classification="thrust"; for rake=-90, assert "normal"; for rake=0, assert "strike_slip"; for rake=170, assert "strike_slip"
- `test_phase_offset_range`: assert phase_offset_thrust_normal in [0.0, 0.5]
- `test_unmatched_check_present`: load results JSON; assert `"unmatched_check"` key present with `"n_unmatched"` field
- `test_primary_conclusion_valid`: assert primary_conclusion is one of "loading", "geometric", "metivier_consistent", "ambiguous"

---

## 7. Whitepaper

In `output/case-b3-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (does solar signal differ by focal mechanism?), GCMT coverage, three predictions tested, key finding, and mechanism implication
2. **Data Source** — ISC-GEM catalog; GCMT focal mechanism join (n=4,874 matched, 52.9% match rate); describe the join strategy and pre-1976 coverage gap; report class sizes
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1, `rules/data-handling.md`)
   - 3.2 Focal mechanism classification: rake quadrant boundaries for thrust/normal/strike-slip; oblique exclusion
   - 3.3 Chi-square, Rayleigh, Cramér's V per mechanism class
   - 3.4 Métivier (2009) reference pattern and ranking comparison
   - 3.5 Anti-phase test for loading hypothesis (thrust vs normal phase offset)
   - 3.6 Unmatched event sensitivity check
4. **Results**
   - 4.1 Coverage and classification: report n per class; describe unmatched event distribution
   - 4.2 Per-mechanism distributions: embed Figure 1; tabulate chi2, p, Cramér's V, Rayleigh R for all three mechanisms at k=24
   - 4.3 Mechanism comparison: embed Figure 2; state Cramér's V rank order; state whether Métivier pattern is matched; state thrust-normal phase offset and classification
   - 4.4 Coverage map: embed Figure 3; describe geographic distribution of each mechanism class
   - 4.5 Unmatched check: report Cramér's V for unmatched events; state representativeness assessment
5. **Cross-Topic Comparison** — compare to Métivier et al. (2009) tidal triggering by tectonic regime (normal > strike-slip > thrust); note that tidal period is different from annual but mechanism geometry is analogous; compare to Johnson et al. (2017b) tectonic stress seasonality
6. **Interpretation** — state which prediction is best supported; discuss GCMT coverage limitation; maintain objectivity
7. **Limitations** — 47.1% null rate in GCMT join; pre-1976 events entirely excluded from matched sample; oblique mechanisms excluded; PB2002 proxy not used as substitute; 52.9% coverage may introduce geographic or magnitude-range selection bias
8. **References** — Métivier et al. (2009), Johnson et al. (2017b), GCMT catalog reference, Adhoc A0b, A4

---

## 8. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case B3: Tectonic Regime Stratification
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [n per mechanism class; Cramér's V rank order; Métivier pattern matched/not; thrust-normal phase offset classification; primary conclusion]
  ```
- Update Case B3 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
