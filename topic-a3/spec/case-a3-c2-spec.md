# Case A3.C2: Targeted Major Sequence Phased Declustering Test

**Status:** Planning

**Reference cases:** A2.B6, A2.A4, A3.B1

**Intent statement:**
Case A3.C2 directly tests whether the global solar-phase signal in the ISC-GEM catalog is diffuse across events or disproportionately concentrated in a small number of major aftershock sequences. A2.B6 identified the 2003–2014 rolling-window cluster as the most statistically significant period, contemporaneous with the 2004 Sumatra M9.1 aftershock sequence. A2.A4 showed that aftershock populations carry stronger solar-phase signals than mainshocks. However, A3.B1 found that rolling-window chi-square significance is negatively correlated with per-window aftershock density (Reasenberg r = −0.574, p < 0.001), suggesting the signal is not simply concentrated in high-density aftershock windows.

A3.C2 resolves this apparent contradiction through a targeted approach: sequentially removing the M≥8.5 events in the raw catalog and their attributed aftershock sequences, recomputing global chi-square, Cramér's V, and A1b interval structure after each removal. If the signal collapses with the first one or two removals, it is sequence-concentrated despite A3.B1's window-level null. If the signal degrades progressively or survives all removals, it is genuinely diffuse. The test runs on four parallel tracks: raw catalog with G-K aftershock removal, raw catalog with Reasenberg aftershock removal, G-K mainshock catalog (mainshock-only removal), and Reasenberg mainshock catalog (mainshock-only removal).

**A3.B1 context for sequence handling:**
A3.B1 found that rolling-window chi-square significance is *negatively* correlated with aftershock count per window in all three declustered catalogs (r = −0.574 for Reasenberg). This means the most chi-square-significant time windows do not coincide with the highest aftershock density windows. Despite this, the raw catalog shows substantially higher chi-square significance than all declustered versions (71% vs. 40% Reasenberg vs. 19% G-K vs. 6.5% A1b). The discrepancy motivates C2: the window-level null does not rule out a few specific mega-sequences (with very large `aftershock_count`) driving the full-catalog signal even while not dominating any particular 10-year window's per-window density. A3.B1's negative correlation is consistent with sequences being distributed across multiple overlapping windows rather than concentrated in one.

**Magnitude threshold decision:**
Use M≥8.5 from the raw ISC-GEM catalog as the primary threshold. This is expected to yield approximately 10–15 events in the 1950–2021 record. Report results for all removal steps; the spec does not pre-commit to a specific "collapse point."

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |
| G-K mainshocks (enriched) | `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv` | 5,883 | + `foreshock_count`, `aftershock_count`, `window_secs`, `window_km` |
| G-K aftershocks (enriched) | `data/iscgem/declustering-algorithm/aftershocks_gk-seq_global.csv` | 3,327 | + `parent_id`, `parent_magnitude`, `delta_t_sec`, `delta_dist_km` |
| Reasenberg mainshocks (enriched) | `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv` | 8,265 | + `foreshock_count`, `aftershock_count`, `window_secs`, `window_km` |
| Reasenberg aftershocks (enriched) | `data/iscgem/declustering-algorithm/aftershocks_reas-seq_global.csv` | 945 | + `parent_id`, `parent_magnitude`, `delta_t_sec`, `delta_dist_km` |

Phase normalization: `phase = (solar_secs / 31_557_600.0) % 1.0` (Julian year constant).

**A1b baseline intervals (k=24 bin mapping):**
- Interval 1: bins 4–5 (phase [0.1667, 0.2500)) — March equinox region
- Interval 2: bin 15 (phase [0.6250, 0.6667)) — mid-August region
- Interval 3: bin 21 (phase [0.8750, 0.9167)) — late-November region

Bin index (0-based) for k=24: `bin_i = floor(phase * 24)`.

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a3/`
- All output paths: `BASE_DIR / "output" / ...`
- All data paths: `BASE_DIR.parent / "data" / "iscgem" / ...`

**Planned outputs:**
- `src/case-a3-c2-analysis.py` — major event identification, sequence metrics, phased removal computation, writes results JSON
- `src/visualization-case-a3-c2.py` — generates all PNG figures
- `tests/test-case-a3-c2.py` — test suite (all tests must pass)
- `output/case-a3-c2-results.json` — removal steps × 4 runs; sequence metrics; per-step chi2, V, interval z-scores
- `output/case-a3-c2-whitepaper.md` — methodology, sequence metrics breakout, results, interpretation
- `output/case-a3-c2-degradation.png` — chi-square p-value + Cramér's V trajectory across removal steps (4 runs)
- `output/case-a3-c2-interval-decay.png` — per-interval z-score change across removal steps
- `output/case-a3-c2-sequence-summary.png` — visual table of sequence metrics for major events

---

## 1. Environment and data loading

In `src/case-a3-c2-analysis.py`:

- Imports: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define paths:
  ```python
  RAW_PATH        = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
  GK_MS_PATH      = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv"
  GK_AS_PATH      = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_gk-seq_global.csv"
  REAS_MS_PATH    = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv"
  REAS_AS_PATH    = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_reas-seq_global.csv"
  ```
- Load all five files; parse `event_at` as `pd.Timestamp` (UTC); log row counts; assert:
  - raw: n=9210; gk_ms: n=5883; gk_as: n=3327; reas_ms: n=8265; reas_as: n=945
- Compute `phase = (solar_secs / 31_557_600.0) % 1.0` on all five DataFrames
- Define constants: `K_BINS = 24`, `MAG_THRESHOLD = 8.5`
- Define A1b interval bins:
  ```python
  INTERVAL_BINS = {
      "interval_1": [4, 5],
      "interval_2": [15],
      "interval_3": [21],
  }
  ```

---

## 2. Major event identification

```python
def identify_major_events(raw_df: pd.DataFrame, mag_threshold: float) -> pd.DataFrame:
```
- Filter raw catalog to `usgs_mag >= mag_threshold`
- Sort by `usgs_mag` descending (largest first — this defines the removal order)
- Return DataFrame with columns: `usgs_id`, `usgs_mag`, `event_at`, `latitude`, `longitude`, `depth`
- Log count of major events found; expected approximately 10–15 for M≥8.5

For each major event, determine its declustering classification in G-K and Reasenberg:
```python
def classify_event_in_catalog(usgs_id: str, mainshock_df: pd.DataFrame, aftershock_df: pd.DataFrame) -> dict:
```
- Returns: `{"is_mainshock": bool, "is_aftershock": bool, "parent_id": str|None}`
- Check if `usgs_id` appears in `mainshock_df["usgs_id"]` → is_mainshock
- Check if `usgs_id` appears in `aftershock_df["usgs_id"]` → is_aftershock (note: aftershock file uses `usgs_id` as the event identifier)
- Log any event that is neither (present in raw but classified as neither mainshock nor aftershock — edge case to report but not error)

---

## 3. Sequence metrics computation

For the whitepaper Report Requirement: compute per-major-event sequence metrics from both G-K and Reasenberg attributions.

```python
def compute_sequence_metrics(
    major_event: pd.Series,
    mainshock_df: pd.DataFrame,
    aftershock_df: pd.DataFrame,
    label: str,
) -> dict:
```

For each major event (`usgs_id = major_event["usgs_id"]`):

**If the event is a mainshock in this declustering:**
- Locate the mainshock row: `ms_row = mainshock_df[mainshock_df["usgs_id"] == usgs_id].iloc[0]`
- `foreshock_count = ms_row["foreshock_count"]`
- `aftershock_count = ms_row["aftershock_count"]`
- `total_sequence = foreshock_count + aftershock_count + 1`
- `window_secs = ms_row["window_secs"]`; `window_days = window_secs / 86400`
- Retrieve attributed aftershocks: `as_rows = aftershock_df[aftershock_df["parent_id"] == usgs_id]`
- **Half-life split**: partition `as_rows` by `delta_t_sec` relative to `window_secs / 2`:
  - `early_count = (as_rows["delta_t_sec"] <= window_secs / 2).sum()`
  - `late_count = (as_rows["delta_t_sec"] > window_secs / 2).sum()`
  - `early_pct = early_count / max(aftershock_count, 1)`

**If the event is an aftershock in this declustering:**
- `foreshock_count = 0`, `aftershock_count = 0`, `total_sequence = 1` (the event itself)
- `window_secs = 0`, `window_days = 0`, `early_count = 0`, `late_count = 0`, `early_pct = None`
- Log: event classified as aftershock in `{label}` — no sequence attributed

**If the event is absent from both mainshock and aftershock files:**
- `foreshock_count = 0`, `aftershock_count = 0`, `total_sequence = 1`, `window_secs = 0`, `window_km = 0`, `early_count = 0`, `late_count = 0`, `early_pct = None`
- Log: event absent from `{label}` — may have been below M≥6.0 threshold for this algorithm

Return dict:
```json
{
  "usgs_id": str,
  "usgs_mag": float,
  "event_at": str,
  "latitude": float,
  "longitude": float,
  "declustering": "gk | reasenberg",
  "classification": "mainshock | aftershock | absent",
  "foreshock_count": int,
  "aftershock_count": int,
  "total_sequence": int,
  "window_secs": float,
  "window_days": float,
  "window_km": float,
  "early_count": int,
  "late_count": int,
  "early_pct": float
}
```

---

## 4. Chi-square and interval statistics

```python
def compute_chi2_stats(phases: np.ndarray, k: int = 24) -> dict:
```
- `n = len(phases)`
- `bin_indices = np.floor(phases * k).astype(int) % k`
- `observed = np.bincount(bin_indices, minlength=k)`
- `expected = np.full(k, n / k)`
- `chi2_stat, p_chi2 = scipy.stats.chisquare(observed, expected)`
- `cramers_v = np.sqrt(chi2_stat / (n * (k - 1)))` if n > 0 else 0.0
- Per interval:
  - `interval_count = observed[bins].sum()`
  - `interval_expected = len(bins) * (n / k)`
  - `interval_z = (interval_count - interval_expected) / np.sqrt(interval_expected)` if `interval_expected > 0` else 0.0
- Return: `{"n": n, "chi2_k24": float, "p_chi2_k24": float, "cramers_v": float, "bin_counts": list[int], "interval_1_z": float, "interval_2_z": float, "interval_3_z": float}`

---

## 5. Phased removal analysis

Four parallel removal tracks:

| Run key | Base catalog | Removal set |
|---------|-------------|-------------|
| `raw_gk` | raw (9,210 events) | major event + G-K aftershocks (via `parent_id`) |
| `raw_reas` | raw (9,210 events) | major event + Reasenberg aftershocks (via `parent_id`) |
| `mainshock_gk` | G-K mainshocks (5,883 events) | major event mainshock row only (if present as mainshock) |
| `mainshock_reas` | Reasenberg mainshocks (8,265 events) | major event mainshock row only (if present as mainshock) |

```python
def run_phased_removal(
    base_df: pd.DataFrame,
    major_events: pd.DataFrame,
    aftershock_df: pd.DataFrame | None,
    run_key: str,
) -> list[dict]:
```

**Step 0 (baseline):** compute stats on full `base_df` with no removals. Record `{"step": 0, "event_removed": null, "ids_removed_cumulative": [], "n_remaining": n_full, "stats": {...}}`.

**Steps 1 through len(major_events):**
For step `i` (1-indexed), removing event `major_events.iloc[i-1]`:
1. Collect IDs to remove at this step:
   - Always add `major_event_usgs_id` to the removal set
   - For `raw_gk` and `raw_reas`: additionally add all `usgs_id` values from `aftershock_df` where `parent_id == major_event_usgs_id`
   - For `mainshock_gk` and `mainshock_reas`: no aftershock removal (mainshock-only catalogs)
2. Accumulate removal set across steps (cumulative removal)
3. Filter `base_df` to exclude all accumulated removed IDs
4. Compute stats via `compute_chi2_stats` on the phases of the remaining events
5. Record step dict:
```json
{
  "step": int,
  "event_removed": {
    "usgs_id": str,
    "usgs_mag": float,
    "event_at": str,
    "latitude": float,
    "longitude": float,
    "n_removed_this_step": int
  },
  "n_removed_cumulative": int,
  "ids_removed_cumulative_count": int,
  "n_remaining": int,
  "stats": {
    "chi2_k24": float,
    "p_chi2_k24": float,
    "cramers_v": float,
    "interval_1_z": float,
    "interval_2_z": float,
    "interval_3_z": float
  }
}
```

Note: if `n_remaining < 24` (degenerate case after extreme removals), record `chi2_k24 = null` and log a warning. This should not occur for M≥8.5 removals on a 9,210-event catalog.

---

## 6. Results JSON structure

```json
{
  "case": "A3.C2",
  "title": "Targeted Major Sequence Phased Declustering Test",
  "parameters": {
    "n_raw": 9210,
    "n_gk_mainshocks": 5883,
    "n_gk_aftershocks": 3327,
    "n_reas_mainshocks": 8265,
    "n_reas_aftershocks": 945,
    "mag_threshold": 8.5,
    "k_bins": 24,
    "julian_year_secs": 31557600.0
  },
  "major_events": [
    {"usgs_id": str, "usgs_mag": float, "event_at": str, "latitude": float, "longitude": float,
     "removal_order": int}
  ],
  "sequence_metrics": {
    "gk": [/* list of compute_sequence_metrics results for each major event */],
    "reasenberg": [/* same */]
  },
  "runs": {
    "raw_gk":       {"steps": [/* step dicts */]},
    "raw_reas":     {"steps": [/* step dicts */]},
    "mainshock_gk": {"steps": [/* step dicts */]},
    "mainshock_reas": {"steps": [/* step dicts */]}
  }
}
```

---

## 7. Visualizations

In `src/visualization-case-a3-c2.py`:

Load `output/case-a3-c2-results.json` and render all figures.

---

**Figure 1 — Chi-square degradation trajectory** (`output/case-a3-c2-degradation.png`):
- 2-row stacked subplot sharing x-axis (removal step, 0 = baseline)
- Row 1: chi-square p-value (log scale, y-range [1e-12, 1.0]); four lines, one per run:
  - `raw_gk`: steelblue solid
  - `raw_reas`: steelblue dashed
  - `mainshock_gk`: red solid
  - `mainshock_reas`: red dashed
  - Horizontal dashed line at p=0.05 (gray)
- Row 2: Cramér's V (linear); same four lines, same colors/styles
- x-axis: integer tick per removal step; label each tick with the magnitude and year of the event removed (e.g., "M9.5\n1960"); baseline labeled "Baseline"
- Legend; title "Phased Removal of M≥8.5 Sequences — Chi-Square and Cramér's V Degradation"
- 300 DPI, publication quality

---

**Figure 2 — Interval z-score decay** (`output/case-a3-c2-interval-decay.png`):
- 3-row stacked subplot sharing x-axis (removal step)
- Each row: one A1b interval's z-score trajectory across removal steps; four lines per row (one per run), same color scheme as Figure 1
- Horizontal dashed line at z=1.96 (α=0.05 one-tailed threshold) in each row
- y-axis minimum: 0 (floor at 0 for readability, since negative z-scores indicate suppression below expected)
- Row labels: "Interval 1 (0.19–0.25)", "Interval 2 (0.625–0.656)", "Interval 3 (0.875–0.917)"
- x-axis: same event-labeled ticks as Figure 1
- Legend; title "A1b Interval Elevation Z-Score Through Sequential Removal"
- 300 DPI

---

**Figure 3 — Sequence summary visual table** (`output/case-a3-c2-sequence-summary.png`):
- Matplotlib table rendered as a figure (not a bar chart); rows = major events (sorted by magnitude desc); columns for G-K and Reasenberg attributes:
  - Event (mag + date)
  - G-K: aftershock_count, window_days, early_pct, classification
  - Reasenberg: aftershock_count, window_days, early_pct, classification
- Color rows alternating white/light gray; highlight rows where classification differs between G-K and Reasenberg
- Font size appropriate for readability at 300 DPI; no axis frame
- Title "Sequence Metrics for M≥8.5 Events — G-K vs. Reasenberg Attribution"
- 300 DPI

---

## 8. Test suite

In `tests/test-case-a3-c2.py`:

- `test_catalog_loads`: assert n=9210 (raw), n=5883 (gk_ms), n=3327 (gk_as), n=8265 (reas_ms), n=945 (reas_as); assert `event_at` parses without NaT; assert `phase` in [0.0, 1.0) for all
- `test_major_event_count`: assert at least 5 events identified at M≥8.5; assert all identified events have `usgs_mag >= 8.5`; assert sorted by `usgs_mag` descending (removal order)
- `test_removal_monotonic`: for each run, assert `n_remaining` is non-increasing across steps (steps are cumulative); assert step 0 `n_remaining` equals the full catalog size for that run
- `test_step0_matches_full_catalog_stats`: for each run, assert step 0 chi2 and Cramér's V match direct computation on the unmodified base catalog (tolerance: 1e-6)
- `test_chi2_bounds`: for all steps and all runs, assert `chi2_k24 >= 0` and `p_chi2_k24` in [0.0, 1.0]
- `test_cramers_v_bounds`: for all steps and all runs, assert `cramers_v >= 0`
- `test_interval_z_type`: for all steps, assert `interval_1_z`, `interval_2_z`, `interval_3_z` are floats
- `test_parent_id_removal`: for `raw_gk` run step 1 (first event removed), assert that the removed event's G-K aftershocks (from aftershock file `parent_id`) are absent from the remaining catalog; assert the major event itself is absent
- `test_mainshock_only_removal`: for `mainshock_gk` run step 1, assert `n_remaining = n_gk_mainshocks - 1` (only one event removed per step in mainshock-only runs, for events present as mainshocks)
- `test_sequence_metrics_structure`: load results JSON; assert `sequence_metrics` contains keys "gk" and "reasenberg"; assert each list has length equal to `len(major_events)`; assert each entry has keys `foreshock_count`, `aftershock_count`, `total_sequence`, `window_days`, `early_count`, `late_count`
- `test_results_json_completeness`: load results JSON; assert all four run keys present; for each run, assert steps list length equals `len(major_events) + 1` (baseline + one per removal)
- `test_raw_baseline_chi2`: assert raw catalog step 0 `chi2_k24` approximately matches the known full-catalog value from A2.A4 (the global chi-square was 69.37 at k=24); tolerance ± 2.0 (some variation expected due to Julian year constant vs. per-year normalization differences)

---

## 9. Whitepaper

In `output/case-a3-c2-whitepaper.md`:

Standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Code [model name]).

### Sections:

1. **Abstract** (150–200 words): state the question (does the global solar-phase signal survive sequential removal of M≥8.5 event sequences?), describe the phased removal design, state the key outcome (signal collapses / degrades progressively / survives), and state the implication for the diffuse-vs.-concentrated framing.

2. **Data Source**: describe all five input files (raw n=9,210; G-K mainshocks n=5,883, aftershocks n=3,327; Reasenberg mainshocks n=8,265, aftershocks n=945); note that aftershock attribution uses `parent_id` for precise sequence membership. Note relationship to A2.A4 and A3.B1.

3. **Methodology**
   - 3.1 Phase-normalized binning: Julian year constant; cite `data-handling.md`
   - 3.2 Major event identification: M≥8.5 threshold rationale; removal order (largest first)
   - 3.3 Sequence attribution: G-K `parent_id` vs. Reasenberg `parent_id` as removal sets; classification check (mainshock vs. aftershock in each declustering)
   - 3.4 Four removal tracks: describe `raw_gk`, `raw_reas`, `mainshock_gk`, `mainshock_reas`
   - 3.5 Per-step statistics: chi-square (k=24), Cramér's V formula, A1b interval z-scores
   - 3.6 Sequence metrics: foreshock/aftershock counts, window duration, half-life split definition
   - 3.7 Relationship to A3.B1: explain why A3.B1's negative rolling-window correlation does not rule out targeted sequence concentration; describe how C2 tests the complementary hypothesis

4. **Results**
   - 4.1 Major events identified: embed Figure 3 (sequence summary table); for each event note the canonical seismological name, magnitude, date, and whether it appears as a mainshock or aftershock in G-K and Reasenberg
   - 4.2 Sequence metrics breakout (Report Requirement):
     - Present a markdown table: Event | Mag | Date | G-K foreshocks | G-K aftershocks | G-K window (days) | G-K early% | G-K classification | Reasenberg aftershocks | Reasenberg window (days) | Reasenberg early% | Reasenberg classification
     - Note which events are "early-loaded" (early_pct > 0.70) vs. "late-loaded" (early_pct < 0.30) — these designations will be referenced in downstream A3.A1
   - 4.3 Chi-square and Cramér's V degradation: embed Figure 1; describe whether the signal collapses early (step 1–2), degrades gradually, or is robust; compare raw-GK vs. raw-Reas removal paths
   - 4.4 Interval-level decay: embed Figure 2; describe which intervals decay first and which persist across all removals; relate to A3.B1's finding that only Interval 2 was partially elevated in any catalog
   - 4.5 Mainshock-only removal: compare `mainshock_gk` and `mainshock_reas` degradation trajectories to the raw-catalog tracks; assess whether declustered signal is similarly sequence-concentrated

5. **Cross-Topic Comparison**: compare to A2.A4 (aftershock phase-preference); compare to A2.B6 (2003–2014 window elevation); compare to A3.B1 (negative rolling-window sequence density correlation). State whether C2 results are consistent or in tension with A3.B1.

6. **Interpretation**: state whether the signal is diffuse or concentrated; note whether the finding changes the interpretation of the A2 main result; guard against both confirmatory framing (eagerly attributing any collapse to sequence artifact) and dismissive framing (assuming robustness without verifying the magnitude of per-step changes).

7. **Limitations**: removal is cumulative — later steps remove both the current event and all previously removed events; the ordering (largest first) is one reasonable choice but a random or chronological ordering could yield different degradation curves; G-K and Reasenberg aftershock attributions differ substantially in n (3,327 vs. 945 G-K vs. Reasenberg aftershocks), so `raw_gk` removes far more events per step than `raw_reas`.

8. **References**: Gardner & Knopoff (1974), Reasenberg (1985), A2.A4, A2.B6, A3.B1

---

## 10. Update context docs

After all outputs are generated and tests pass:

- Append to `topic-a3/docs/topic-summary.md`:
  ```
  ## Case A3.C2: Targeted Major Sequence Phased Declustering Test
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [signal classification (diffuse/concentrated/partially concentrated); step at which p crosses 0.05 for raw_gk and raw_reas; final step Cramér's V vs. baseline; interval persistence summary; early/late classification for major events]
  ```
- Update Case C2 status from `Planning` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a3/CLAUDE.md` Case Table
