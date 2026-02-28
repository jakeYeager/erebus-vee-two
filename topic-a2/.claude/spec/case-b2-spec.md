# Case B2: Ocean vs. Continent Location — Hydrological Loading Discrimination

**Status:** Pending

**Intent statement:**
Case B2 tests whether the solar-phase signal persists in purely oceanic seismicity, which would rule out hydrological and snow loading as sufficient explanations because those are land-surface processes. Events are classified into three categories — oceanic, transitional, and continental — using a distance-to-coastline approach, and solar-phase bin distributions are computed independently for each subset. Three coastline classification files are available (GSHHG, Natural Earth, PB2002 proximity) and all three are used in order to assess sensitivity to the classification method. The tidal literature (Scholz et al. 2019) adds a complication: strong tidal modulation has been found at Axial Volcano (mid-ocean ridge), establishing that oceanic settings can show periodic signals via mechanisms other than continental hydrology, so a positive oceanic result requires careful mechanism discrimination. Visualizations include global maps and zoomed regional maps (Philippines, Japan, Chile, Java) for each of the three classification treatments. Color scheme: oceanic=blue, transitional=green, continental=red. No published study has used ocean/continent separation as a mechanism discriminator for annual solar seismicity at global M≥6.0 scale.

**Relationship to prior topics:**
The ocean/continent classification files were produced specifically for this topic (data-requirements.md REQ-2) using three methods: GSHHG coastline (primary, highest resolution), Natural Earth coastline (secondary), and PB2002 plate boundary proximity (coarse proxy). Topic L4 confirmed `latitude` and `longitude` are present in the ISC-GEM catalog schema. Adhoc A1b's PB2002 boundary proximity analysis on elevated-bin events already provides partial context for the geographic distribution of the signal, making this case's geographic analysis a formal extension of that exploratory work.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/global-sets/iscgem_global_events.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `latitude`, `longitude`, `depth` |
| GSHHG classification (primary) | `data/iscgem/plate-location/ocean_class_gshhg_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |
| Natural Earth classification (secondary) | `data/iscgem/plate-location/ocean_class_ne_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |
| PB2002 classification (coarse proxy) | `data/iscgem/plate-location/ocean_class_pb2002_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |

Classification labels: `oceanic`, `transitional`, `continental`
GSHHG/NE counts: continental=3,799 / transitional=3,459 / oceanic=1,952
PB2002 counts: continental=2,851 / transitional=3,677 / oceanic=2,682

Join: merge classification files onto raw catalog by `usgs_id`; assert 9210 matched rows for each classification file.

Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year [confirm before running: consistent with prior cases — verify uniform Julian constant vs per-year values].

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-b2-analysis.py` — main analysis: join classification files, per-class bin distributions for all three classification methods, chi-square tests, writes results JSON
- `src/visualization-case-b2.py` — generates all PNG figures
- `tests/test-case-b2.py` — test suite
- `output/case-b2-results.json` — per-method, per-class statistics and prediction evaluation
- `output/case-b2-whitepaper.md` — methodology, results, cross-topic comparison
- `output/case-b2-global-map-gshhg.png` — global map of events colored by GSHHG classification
- `output/case-b2-global-map-ne.png` — global map colored by NE classification
- `output/case-b2-global-map-pb2002.png` — global map colored by PB2002 classification
- `output/case-b2-regional-maps.png` — 4-row x 3-column grid: regions (Philippines, Japan, Chile, Java) x classification methods
- `output/case-b2-binplots-gshhg.png` — 3-panel bin distributions (oceanic, transitional, continental) using GSHHG at k=24
- `output/case-b2-binplots-ne.png` — same using Natural Earth
- `output/case-b2-binplots-pb2002.png` — same using PB2002

---

## 1. Environment and data loading

In `src/case-b2-analysis.py`:

- Import: `pandas`, `numpy`, `scipy.stats`, `matplotlib`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define paths using BASE_DIR (see Data context block)
- Load raw catalog; assert n=9210
- Load all three classification files; assert n=9210 for each
- Join each classification file onto raw catalog by `usgs_id`:
  - `df_gshhg = raw.merge(class_gshhg[["usgs_id","ocean_class","dist_to_coast_km"]], on="usgs_id", how="left")`
  - Rename `ocean_class` to `ocean_class_gshhg`, `dist_to_coast_km` to `dist_km_gshhg`
  - Repeat for NE and PB2002 with appropriate column name suffixes
  - Assert no NaN in the ocean_class columns after join
- Compute phase for full merged dataframe using Julian year constant
- Define classification methods list: `[("gshhg", "ocean_class_gshhg"), ("ne", "ocean_class_ne"), ("pb2002", "ocean_class_pb2002")]`
- Define location classes: `["oceanic", "transitional", "continental"]`

---

## 2. Per-class bin statistics (all three classification methods)

For each classification method (gshhg, ne, pb2002) and each class (oceanic, transitional, continental) at each of k=16, 24, 32:

1. Filter catalog to class: `subset = df[df[ocean_class_col] == class_label]`
2. Log n per class; warn if n < 200 for any class
3. Compute observed bin counts, expected counts, chi-square, p-value, Cramér's V, Rayleigh R, Rayleigh p, mean phase angle (same procedure as A4)
4. Elevated bins at 1-SD threshold; elevated phase intervals

Output structure under key `"class_stats"` in results JSON:
```json
{
  "class_stats": {
    "gshhg": {
      "oceanic": {
        "n": 1952,
        "k16": {...}, "k24": {...}, "k32": {...}
      },
      "transitional": {...},
      "continental": {...}
    },
    "ne": {...},
    "pb2002": {...}
  }
}
```

---

## 3. Classification method sensitivity comparison

For each class, compare results across the three methods at k=24:
- Tabulate chi2, p_chi2, and Cramér's V for each method/class combination in a 3×3 matrix (classes × methods)
- Compute average Cramér's V across methods for each class
- Flag any class where the three methods disagree on significance (some p<0.05, some p>=0.05)
- Record: `{"class": str, "gshhg_p": float, "ne_p": float, "pb2002_p": float, "agreement": bool}`

Output under key `"method_sensitivity"` in results JSON.

---

## 4. Prediction evaluation

Using GSHHG (primary) results at k=24:

- **Geometric hypothesis**: signal appears in both oceanic and continental subsets → evaluate: oceanic p_chi2 < 0.05 AND continental p_chi2 < 0.05
- **Hydrological hypothesis**: signal confined to or significantly stronger in continental → evaluate: continental Cramér's V > 2× oceanic Cramér's V, OR oceanic p_chi2 >= 0.05

Additional diagnostic:
- If oceanic signal is significant and matches A1b phase intervals → note possible tidal-like mechanism at mid-ocean ridges (Scholz et al. 2019 precedent)

Output under key `"prediction_evaluation"` in results JSON:
```json
{
  "prediction_evaluation": {
    "primary_method": "gshhg",
    "oceanic_significant": bool,
    "continental_significant": bool,
    "continental_stronger": bool,
    "oceanic_matches_a1b_intervals": bool,
    "geometric_supported": bool,
    "hydrological_supported": bool,
    "primary_conclusion": "geometric | hydrological | ambiguous"
  }
}
```

---

## 5. Visualizations

In `src/visualization-case-b2.py`:

**Figures 1–3 — Global maps** (`output/case-b2-global-map-gshhg.png`, `-ne.png`, `-pb2002.png`):
- World map using `matplotlib` with `cartopy` or `basemap` [confirm before running: verify cartopy or basemap availability in project Python environment; if neither is available, use scatter plot on lat/lon axes with coastline data from the shapely/geopandas/Natural Earth path]
- Each event plotted as a point: oceanic=blue, transitional=green, continental=red
- Point size proportional to `usgs_mag` (e.g., size = (usgs_mag - 5.5) ** 2 * 2)
- Title: "ISC-GEM Events — [Method] Classification"
- Legend for color classes; include n per class in legend
- 300 DPI

**Figure 4 — Regional maps** (`output/case-b2-regional-maps.png`):
- 4-row × 3-column grid; rows = Philippines (lat 5–22, lon 115–130), Japan (lat 28–47, lon 128–148), Chile (lat -58 – -15, lon -80 – -60), Java (lat -12 – -4, lon 100–115); columns = GSHHG, NE, PB2002
- Same color scheme and point sizing as global maps; no coastline background required (lat/lon scatter sufficient); add lat/lon grid lines
- Title each column by method name; title each row by region name
- 300 DPI

**Figures 5–7 — Bin distribution plots** (`output/case-b2-binplots-gshhg.png`, `-ne.png`, `-pb2002.png`):
- 3-panel layout: oceanic, transitional, continental at k=24
- Same format as A4/B1 bin plots: horizontal steelblue bars, dashed expected-count line, 1-SD threshold, A1b baseline interval gray bands
- Annotate each panel with class label, n, χ², p-value, Cramér's V
- 300 DPI

---

## 6. Test suite

In `tests/test-case-b2.py`:

- `test_catalog_load`: assert raw catalog n=9210
- `test_classification_joins`: assert all three classification joins produce 9210 rows; assert no NaN in ocean_class columns
- `test_class_counts_gshhg`: assert GSHHG: oceanic=1952, transitional=3459, continental=3799 (within ±1 for any row exclusions)
- `test_class_counts_pb2002`: assert PB2002: oceanic=2682, transitional=3677, continental=2851 (within ±1)
- `test_phase_range`: assert all phases in [0.0, 1.0)
- `test_class_partition`: for each method, assert oceanic + transitional + continental == 9210
- `test_chi_square_all_classes`: assert chi2 and p_chi2 are finite for all 3 methods × 3 classes × 3 bin counts
- `test_cramer_v_range`: assert Cramér's V in [0.0, 1.0] for all
- `test_method_sensitivity_structure`: load results JSON; assert `"method_sensitivity"` key present with correct structure
- `test_prediction_evaluation_keys`: assert `"prediction_evaluation"` key present; assert `"primary_conclusion"` is "geometric", "hydrological", or "ambiguous"
- `test_classification_labels_valid`: assert all ocean_class values are one of "oceanic", "transitional", "continental"

---

## 7. Whitepaper

In `output/case-b2-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (does oceanic seismicity show solar signal?), three classification methods, geographic visualization scope, primary finding, and mechanism implication
2. **Data Source** — ISC-GEM catalog (n=9,210); describe three classification files (GSHHG, NE, PB2002), their methods, and counts per class; report join success
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1, `rules/data-handling.md`)
   - 3.2 Ocean/continent classification methods: GSHHG, Natural Earth, PB2002 proximity; classification thresholds (>200 km = oceanic, <50 km = continental, 50–200 km = transitional)
   - 3.3 Chi-square, Rayleigh, Cramér's V per class per method
   - 3.4 Classification method sensitivity comparison
   - 3.5 Prediction evaluation criteria
   - 3.6 Visualization approach: global maps and regional zooms
4. **Results**
   - 4.1 Global maps: embed Figures 1–3; describe geographic distribution of classes
   - 4.2 Regional maps: embed Figure 4; describe regional classification pattern differences between methods
   - 4.3 Bin distributions (GSHHG primary): embed Figure 5; tabulate chi2, p, Cramér's V for oceanic, transitional, continental
   - 4.4 Method sensitivity: tabulate chi2, p, Cramér's V for NE and PB2002; note agreement or disagreement with GSHHG; embed Figures 6–7
   - 4.5 Prediction evaluation: state whether signal is confined to continental, present in oceanic, or both; state primary conclusion
5. **Cross-Topic Comparison** — compare to Adhoc A1b elevated-bin boundary proximity results; compare to Scholz et al. (2019) tidal signal at Axial Volcano; note implication of any oceanic signal finding
6. **Interpretation** — state mechanism conclusion; note sensitivity of result to classification method; maintain objectivity
7. **Limitations** — PB2002 proximity classification is a coarse proxy; transitional class is large (~37% of events) and mechanistically ambiguous; classification thresholds (50/200 km) are somewhat arbitrary; visual map resolution limited
8. **References** — Scholz et al. (2019), Bird (2003), Adhoc A1b, Case 3A

---

## 8. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case B2: Ocean vs. Continent Location — Hydrological Loading Discrimination
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [chi2/p/Cramér's V for oceanic/transitional/continental under GSHHG at k=24; method sensitivity agreement; primary conclusion (geometric/hydrological/ambiguous)]
  ```
- Update Case B2 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
