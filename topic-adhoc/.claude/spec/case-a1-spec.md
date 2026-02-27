# Case A1 Spec: Binning Increment Sensitivity Analysis

## Status

Ready for implementation.

---

## Context

Legacy Approach Three used 16 bins for all three astronomical metrics (`solar_secs`, `lunar_secs`, `midnight_secs`). That bin count was base-2 convenient but methodologically arbitrary. Additionally, the `lunar_secs` 16-bin result was found to be a binning artifact (p=0.0179 → p=0.606 after phase normalization). This case tests whether the `solar_secs` chi-square signal is robust across bin counts (8, 16, 24, 32) or dependent on the 16-bin choice, using phase-normalized binning as a new project standard applied consistently to all metrics.

Note: the legacy Approach Three results quoted in the planning doc used the ComCat catalog. This case uses ISC-GEM. The 16-bin baseline here will differ from legacy values; this is expected and not a discrepancy.

---

## Data Input

| File | Path |
| --- | --- |
| ISC-GEM | `data/global-sets/iscgem_global_events.csv` |

Path is relative to project root (`/Users/jake/Projects/Code/prod/erebus-vee-two/`).

Schema: `usgs_id, usgs_mag, event_at, solaration_year, solar_secs, lunar_secs, midnight_secs, latitude, longitude, depth`

Expected event count: 9,210.

---

## Deliverables

All outputs go under `topic-adhoc/`.

| Type | Path |
| --- | --- |
| Analysis script | `src/case-a1-analysis.py` |
| Visualization script | `src/visualization-case-a1.py` |
| Test suite | `tests/test-case-a1.py` |
| Results JSON | `output/case-a1-results.json` |
| Whitepaper | `output/case-a1-whitepaper.md` |
| Image: p-value sweep | `output/case-a1-pvalue-sweep.png` |
| Image: distributions grid | `output/case-a1-distributions.png` |

---

## Analysis Script: `case-a1-analysis.py`

Performs phase-normalized chi-square tests for all 12 (metric × bin_count) combinations and writes `case-a1-results.json`.

### Phase Normalization

Apply to all three metrics before binning. Normalization maps each event's raw seconds value to a phase in `[0, 1)` by dividing by the actual cycle length for that event.

**`solar_secs`**

The cycle length is the actual calendar year length in seconds for the event's `solaration_year`:
- Leap year: 366 × 86,400 = 31,622,400 seconds
- Non-leap year: 365 × 86,400 = 31,536,000 seconds
- Leap year rule: divisible by 4, except centuries unless divisible by 400

`phase = solar_secs / year_length_seconds(solaration_year)`

**`lunar_secs`**

The cycle length is the mean synodic month: 29.53059 days × 86,400 = 2,551,442.976 seconds. Use this constant for all events.

`phase = lunar_secs / 2551442.976`

**`midnight_secs`**

The cycle length is a fixed day: 86,400 seconds.

`phase = midnight_secs / 86400`

All three should produce values strictly in `[0, 1)`. Assert this after computing. Any value ≥ 1.0 indicates a data issue and should raise an error with the offending row.

### Binning

For each phase value and bin count k ∈ {8, 16, 24, 32}:

`bin_index = floor(phase * k)`  →  integer in `[0, k-1]`

Count events per bin. Expected count under null = 9210 / k.

### Chi-Square Test

For each (metric, bin_count):
- `chi2, p = scipy.stats.chisquare(observed_counts)`  (uniform expected)
- Cramér's V = `sqrt(chi2 / (n * (k - 1)))` where n=9210, k=bin_count
- Bonferroni-corrected significance: `significant = p < (0.05 / 12)`

### Results JSON: `case-a1-results.json`

```json
{
  "generated": "<ISO timestamp>",
  "data_source": "data/global-sets/iscgem_global_events.csv",
  "event_count": 9210,
  "methodology": {
    "phase_normalization": true,
    "normalization_factors": {
      "solar_secs": "actual_calendar_year_length_by_solaration_year",
      "lunar_secs": "mean_synodic_month_2551442.976s",
      "midnight_secs": "fixed_day_86400s"
    },
    "bin_counts": [8, 16, 24, 32],
    "n_tests": 12,
    "alpha_uncorrected": 0.05,
    "alpha_bonferroni": 0.004167,
    "robustness_criterion": "significant at >= 3 of 4 bin counts after Bonferroni correction"
  },
  "results": {
    "solar_secs": {
      "8":  {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "16": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "24": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "32": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "n_significant": null,
      "robust": null
    },
    "lunar_secs": {
      "8":  {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "16": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "24": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "32": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "n_significant": null,
      "robust": null
    },
    "midnight_secs": {
      "8":  {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "16": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "24": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "32": {"chi2": null, "p_value": null, "cramers_v": null, "significant_bonferroni": null},
      "n_significant": null,
      "robust": null
    }
  }
}
```

---

## Visualization Script: `visualization-case-a1.py`

Reads `output/case-a1-results.json` and produces two PNG images.

### Image 1: `case-a1-pvalue-sweep.png`

**Primary robustness visualization.**

Line chart showing p-value vs bin count for all three metrics.

- X-axis: bin count (8, 16, 24, 32); treat as categorical with equal spacing
- Y-axis: p-value on a **log10 scale**, inverted so smaller p-values are higher
- One line per metric:
  - `solar_secs`: steelblue, solid line with circle markers
  - `lunar_secs`: coral, solid line with square markers
  - `midnight_secs`: gray, dashed line with triangle markers
- Horizontal reference lines:
  - α = 0.05 (uncorrected): thin gray dashed
  - α = 0.004167 (Bonferroni): red dashed, labeled "Bonferroni α"
- Annotate each data point with its p-value (2 significant figures)
- Title: "Chi-Square p-value by Bin Count (Phase-Normalized)"
- Legend
- 300 DPI

### Image 2: `case-a1-distributions.png`

**Distribution grid.**

3 rows (one per metric) × 4 columns (one per bin count) = 12 subplots, each showing the phase-normalized event distribution as a bar chart.

- Bar color: steelblue
- Each subplot includes:
  - A horizontal dashed red line at the expected uniform count (9210 / bin_count)
  - X-axis label: "Bin" (integers 1 to k)
  - Y-axis label: "Count" (on leftmost column only)
  - Subplot title: e.g., "`solar_secs` — 16 bins (p=0.0023)"
- Row labels on the right margin identifying the metric
- Overall figure title: "Phase-Normalized Distributions by Metric and Bin Count"
- Figure size should be large enough that individual subplots are readable (e.g., 16×10 inches)
- 300 DPI

---

## Test Suite: `test-case-a1.py`

All tests must pass.

1. **Results JSON exists and is valid**: file present, parses as JSON, all required keys present
2. **Event count**: `event_count == 9210`
3. **All 12 results populated**: no `null` values remaining in `results`
4. **Chi-square validity**: for each (metric, bin_count), expected count per bin = 9210 / bin_count ≥ 5 — assert all ≥ 5 (this is structural, will always pass for the chosen counts, but makes the guarantee explicit)
5. **Bonferroni threshold**: `alpha_bonferroni` in JSON ≈ 0.05 / 12 (within floating point tolerance)
6. **Bin sum integrity**: for each (metric, bin_count), recompute binning from scratch from the raw CSV and verify the sum of bin counts == 9210 (confirms no events were dropped during phase normalization)
7. **Phase bounds**: assert all phase values for all three metrics are in `[0, 1)` — verify by recomputing from raw CSV
8. **`midnight_secs` behavioral check**: assert all four `midnight_secs` p-values are > 0.004167 (Bonferroni threshold); if any are significant, raise an informative assertion error rather than a silent failure, noting that this warrants investigation of the test setup
9. **Robustness flag consistency**: for each metric, assert `n_significant` matches the count of `significant_bonferroni: true` entries, and `robust` is `true` iff `n_significant >= 3`
10. **Both PNG images exist** at expected paths

---

## Whitepaper: `case-a1-whitepaper.md`

Follow report writing rules (header/footer template, numbered sections).

**Title:** `Case A1: Binning Increment Sensitivity Analysis for Astronomical Metrics`

**Sections:**

1. **Abstract** — States the purpose: test whether the chi-square signal for `solar_secs` is robust across bin counts (8, 16, 24, 32) using phase-normalized binning as a project standard. Summarizes the robustness finding for each metric and whether the `midnight_secs` behavioral reference held.

2. **Data Source** — ISC-GEM catalog, n=9,210, M ≥ 6.0, 1950–2021. One-line note that this differs from the ComCat-based legacy Approach Three results; the 16-bin baseline here will not replicate legacy values and is not expected to.

3. **Methodology**

   3.1 Phase Normalization — State that phase-normalized binning is used for all metrics per project standard (cite `rules/data-handling.md` and the establishing case, Approach Three Case 1B.3.1). Document the normalization factor for each metric: actual calendar year length for `solar_secs`; mean synodic month (2,551,442.976 s) for `lunar_secs`; fixed 86,400 s for `midnight_secs`.

   3.2 Chi-Square Test — Uniform expected distribution. Cramér's V reported as effect size measure.

   3.3 Multiple Comparisons — 12 tests (3 metrics × 4 bin counts). Bonferroni correction: α = 0.05/12 ≈ 0.0042. Robustness criterion: significant at ≥3 of 4 bin counts.

   3.4 Chi-Square Validity — State expected counts per bin for all four k values and confirm all exceed the ≥5 threshold.

4. **Results**

   4.1 Summary Table — Reproduce the full 12-cell matrix of chi-square, p-value, Cramér's V, and Bonferroni significance for all three metrics and four bin counts. Include a "Robust?" column per metric.

   4.2 P-Value Sweep — Embed `case-a1-pvalue-sweep.png`. Describe the pattern across bin counts for each metric.

   4.3 Distribution Grid — Embed `case-a1-distributions.png`. Note any bins that show consistent elevation or depression across bin counts.

5. **Interpretation** — Objective assessment. Address each metric:

   - `solar_secs`: is the signal robust by the stated criterion? Does the pattern hold consistently across bin counts or show sensitivity to a particular choice?
   - `lunar_secs`: does the signal remain non-significant after phase normalization across all bin counts? Does any bin count recover a significant result?
   - `midnight_secs`: did it behave as observed at global scale (non-significant across all bin counts)? Note this reflects behavior at global scale and does not generalize to stratified or regional analyses.

   Do not editorialize on whether the `solar_secs` result implies a physical mechanism. State only what the bin-count sweep shows about signal stability.

6. **Limitations** — Phase normalization for `lunar_secs` uses the mean synodic month; individual month lengths vary (29.18–29.93 days) and per-event actual month length was not computed. The four bin counts tested are not exhaustive; intermediate or higher values (e.g., 48, 64) could reveal non-monotonic behavior not captured here.

7. **References** — Approach Three Case 1B.3.1 (lunar phase normalization). `rules/data-handling.md` (phase normalization standard).

---

## Acceptance Criteria

- All 10 tests in `test-case-a1.py` pass
- Both PNG images exist in `output/` at 300 DPI
- `case-a1-results.json` present, valid, no `null` values
- `case-a1-whitepaper.md` present with all 7 sections populated and both images embedded inline
