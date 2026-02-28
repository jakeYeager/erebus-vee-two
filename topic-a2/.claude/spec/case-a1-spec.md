# Case A1: Schuster Spectrum and MFPA Periodicity Analysis

**Status:** Pending

**Intent statement:**
Case A1 applies two complementary periodicity detection methods — the cluster-robust Schuster Spectrum Test (Park et al. 2021) and Modified Fourier Power Analysis (MFPA, Dutilleul et al. 2015) — to the ISC-GEM catalog in order to formally characterize the periodic structure of the solar seismic signal. The standard Schuster test is not used alone because Bradley & Hubbard (2024) specifically identified two failure modes relevant here: aftershock temporal clustering biases phase distributions, and a ~20% tidal asymmetry in astronomical forcing creates spurious signal. The cluster-robust variant (Park et al. 2021) is designed to address the first failure mode. MFPA scans the full range from 6 hours to 18 months, testing whether the annual (365.25-day) signal is accompanied by sub-annual harmonics (6-month, 4-month). Results are run on three catalog versions — raw, G-K declustered, and A1b-informed declustered — to test whether G-K over-suppresses the periodicity signal relative to the data-informed window. Sub-annual harmonics are cross-referenced against the three A1b phase intervals: a 6-month harmonic would predict peaks at phases ~0.19 and ~0.69, consistent with intervals 1 and potentially 2, but interval 3 (phase ~0.90) does not fit a 6-month pattern, suggesting residual sequence contamination or a more complex structure.

**Relationship to prior topics:**
Ader & Avouac (2013) applied the standard Schuster spectrum to Himalayan seismicity and detected annual periodicity with 40% amplitude. Dutilleul et al. (2015) applied MFPA to Parkfield and resolved 6-month and 4-month sub-annual harmonics. Case A1 replicates these methods at global scale using the ISC-GEM catalog. Adhoc A1 established that the chi-square signal is robust across k=16, 24, 32 bins (Cramér's V 0.016–0.018), providing the bin-count rationale for the MFPA supplementary analysis. The three-catalog comparison links directly to A4's declustering sensitivity findings — A1 must be run after A4, and the A4 suppression results inform interpretation of any periodicity weakening under G-K.

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM raw catalog | `data/global-sets/iscgem_global_events.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solaration_year`, `solar_secs`, `lunar_secs`, `midnight_secs`, `latitude`, `longitude`, `depth` |
| G-K mainshocks (after A4) | `data/iscgem/declustering-algorithm/mainshocks_G-K_global.csv` | 5,883 | same schema as raw |
| A1b mainshocks (after A4) | `data/iscgem/declustering-algorithm/mainshocks_a1b_global.csv` | 7,137 | same schema as raw |

Parse `event_at` as UTC datetime. Phase normalization: `phase = (solar_secs / year_length_secs) % 1.0` using Julian year (31,557,600 s) [confirm before running: same as A4/B6 — verify whether per-year solar year lengths are pre-computed or Julian constant used uniformly].

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a2/`
- Cross-topic output paths: not required (A4 results consumed within this topic)
- All output paths: `BASE_DIR / "output" / ...`

**Planned Outputs:**
- `src/case-a1-analysis.py` — main analysis: Schuster spectrum and MFPA on all three catalogs; writes results JSON
- `src/case-a1-schuster.py` — cluster-robust Schuster spectrum implementation (callable module)
- `src/case-a1-mfpa.py` — MFPA scan implementation (callable module)
- `src/visualization-case-a1.py` — generates all PNG figures
- `tests/test-case-a1.py` — test suite
- `output/case-a1-results.json` — all Schuster and MFPA results for all three catalogs
- `output/case-a1-whitepaper.md` — methodology, results, cross-topic comparison
- `output/case-a1-schuster-spectrum.png` — Schuster power spectrum plots for all three catalogs
- `output/case-a1-mfpa-scan.png` — MFPA periodogram from 6 hours to 18 months for all three catalogs
- `output/case-a1-harmonic-intervals.png` — A1b phase intervals vs detected harmonic phases overlay

---

## 1. Environment and data loading

In `src/case-a1-analysis.py`:

- Import: `pandas`, `numpy`, `scipy`, `scipy.signal`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define paths using BASE_DIR (see Data context block)
- Load three CSVs; parse `event_at` as UTC datetime; log row counts; assert n=9210, 5883, 7137
- Compute phase for each catalog using phase normalization function (same as A4)
- Convert `event_at` to `event_time_days` (decimal days since 1950-01-01 00:00:00 UTC) for use in spectral analysis

---

## 2. Cluster-robust Schuster Spectrum Test

In `src/case-a1-schuster.py`:

The cluster-robust Schuster test (Park et al. 2021) modifies the standard Schuster test to account for temporal clustering.

**Standard Schuster statistic:**
- For a set of event times `{t_i}` and test period `T`:
  - Phase angles: `φ_i = 2π * (t_i mod T) / T`
  - Schuster statistic: `D² = (1/n) * |Σ exp(i * φ_i)|²` = `(1/n) * [(Σ cos φ_i)² + (Σ sin φ_i)²]`
  - Standard p-value: `p_standard = exp(-D²)` (assumes independent events)

**Cluster-robust modification (Park et al. 2021):**
- Group events into temporal clusters using a minimum inter-event time threshold `dt_cluster` (use `dt_cluster = 1 day` as default)
- For each cluster, compute a single representative phase (mean phase of the cluster) rather than summing all events individually
- Compute effective n as number of clusters rather than number of events
- Compute D² and p-value using the reduced cluster set

Implementation steps:
1. Sort events by time
2. Compute inter-event times `IET[i] = t[i+1] - t[i]` in days
3. Assign cluster IDs: new cluster when `IET[i] >= dt_cluster`
4. Compute mean phase per cluster
5. Compute D² using cluster mean phases and `n_clusters`
6. Return `{"n_events": int, "n_clusters": int, "D2": float, "p_standard": float, "p_cluster_robust": float}` for each test period

**Period scan range:**
- Scan periods: [0.25 days (6h), 0.5, 1, 2, 7, 13.67, 14.77, 27.32, 29.53, 91.3, 182.6, 365.25, 547.9 (18 months)]
- Record full Schuster power spectrum across 200 log-spaced periods from 0.25 days to 548 days
- Tidal periods to explicitly test: 12h (0.5 day), 24h (1 day), 14.77 days (synodic fortnight), 27.32 days (sidereal month)
- Annual periods to test: 365.25 days (tropical year), 182.6 days (half-year), 121.8 days (third-year)

For each of three catalogs, run the full period scan and the explicit tidal/annual period tests. Record all results.

Output structure under key `"schuster"` in results JSON:
```json
{
  "schuster": {
    "raw": {
      "spectrum": [{"period_days": float, "D2": float, "p_standard": float, "p_cluster_robust": float}, ...],
      "explicit_tests": {
        "annual_365": {...}, "half_year_182": {...}, "third_year_122": {...},
        "tidal_12h": {...}, "tidal_24h": {...}, "tidal_14d": {...}, "tidal_27d": {...}
      }
    },
    "gk_mainshocks": {...},
    "a1b_mainshocks": {...}
  }
}
```

---

## 3. MFPA Periodicity Analysis

In `src/case-a1-mfpa.py`:

MFPA (Dutilleul et al. 2015) is a modified Fourier approach that accounts for uneven sampling (irregular event times) via a weighted periodogram.

**Implementation:**
- For each test period `T` (log-spaced from 0.25 days to 548 days, 300 periods):
  1. Compute phase angles `φ_i = 2π * (t_i mod T) / T` for all events
  2. Compute the Fourier sum `F = Σ exp(i * φ_i)` and MFPA power `P_T = |F|² / n`
  3. Compute the significance of `P_T` against a null distribution of uniformly distributed phases:
     - Generate 1,000 bootstrap catalogs of n uniform random phases
     - Compute `P_null` for each bootstrap; record the 95th and 99th percentiles
     - MFPA p-value: proportion of bootstrap P_null values >= observed P_T
- Identify all peaks where `P_T > P_null_95th` (marginally significant) and `P_T > P_null_99th` (strongly significant)
- For each significant peak, record peak period, power, and p-value

**Cross-reference with A1b intervals:**
- For each significant period T, compute the expected peak phases: `k / (n_harmonics * T)` for integer k
- Determine which A1b baseline intervals (0.1875–0.25, 0.625–0.656, 0.875–0.917 in phase fraction) are predicted by a peak at T
- Record whether T is "consistent with interval 1 only", "consistent with intervals 1+3 (6-month)", "consistent with all three intervals (4-month if applicable)", or "inconsistent with all A1b intervals"

Output structure under key `"mfpa"` in results JSON:
```json
{
  "mfpa": {
    "raw": {
      "significant_periods": [
        {"period_days": float, "power": float, "p_mfpa": float, "a1b_consistency": "..."}
      ],
      "spectrum": [{"period_days": float, "power": float, "p95_threshold": float, "p99_threshold": float}, ...]
    },
    "gk_mainshocks": {...},
    "a1b_mainshocks": {...}
  }
}
```

---

## 4. Visualizations

In `src/visualization-case-a1.py`:

**Figure 1 — Schuster spectrum** (`output/case-a1-schuster-spectrum.png`):
- 3-panel stacked plot, one panel per catalog (raw, G-K mainshocks, A1b mainshocks)
- x-axis: period in days, log scale from 0.25 to 548; label key periods (12h, 24h, 14d, 27d, 91d, 182d, 365d, 548d)
- y-axis: cluster-robust Schuster p-value (log scale); horizontal dashed line at p=0.05 and p=0.001
- Plot standard p-value as thin gray line; cluster-robust p-value as thick steelblue line
- Mark the annual, half-year, and third-year periods with vertical dotted lines (red=annual, orange=half-year, green=third-year); mark tidal periods with dotted gray lines
- Annotate n and n_clusters for each panel
- 300 DPI

**Figure 2 — MFPA periodogram** (`output/case-a1-mfpa-scan.png`):
- 3-panel stacked plot, same layout as Figure 1
- y-axis: MFPA power; include 95th percentile and 99th percentile significance thresholds as dashed horizontal lines at each period (or as envelope curves)
- Mark significant peaks with filled orange triangles; label each with its period
- 300 DPI

**Figure 3 — Harmonic-interval overlay** (`output/case-a1-harmonic-intervals.png`):
- Single panel: phase fraction axis (0–1)
- Show A1b baseline intervals as gray shaded bands labeled "Interval 1", "Interval 2", "Interval 3"
- For each significant MFPA period detected, show predicted peak phases as vertical colored lines (one color per detected period)
- Legend identifying each period's color
- 300 DPI

---

## 5. Test suite

In `tests/test-case-a1.py`:

- `test_catalog_loads`: assert all three catalogs load with correct row counts (9210, 5883, 7137)
- `test_event_time_days_monotonic`: assert `event_time_days` is non-decreasing after sort by `event_at`
- `test_schuster_uniform`: generate 1000 uniform random times over 72 years; assert standard Schuster p > 0.05 for annual period in > 90% of 100 bootstrap trials
- `test_schuster_annual_signal`: generate 1000 events with strong annual clustering (90% near phase 0.2); assert Schuster standard p < 0.001 at 365.25 days
- `test_cluster_robust_n_clusters`: assert `n_clusters <= n_events` for all catalogs at all periods
- `test_cluster_robust_p_ge_standard`: assert cluster-robust p >= standard p for all computed periods (clustering can only inflate standard test, not deflate it)
- `test_mfpa_spectrum_length`: assert MFPA spectrum has exactly 300 entries
- `test_mfpa_period_range`: assert min MFPA period >= 0.25 days and max <= 548 days
- `test_mfpa_power_positive`: assert all MFPA power values >= 0
- `test_a1b_crossref_format`: assert all MFPA spectrum entries have an `"a1b_consistency"` field
- `test_results_json_structure`: load `output/case-a1-results.json`; assert keys "schuster" and "mfpa" present; assert each contains "raw", "gk_mainshocks", "a1b_mainshocks"

---

## 6. Whitepaper

In `output/case-a1-whitepaper.md`:

Use the standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Sonnet 4.6).

Sections:
1. **Abstract** — 150–200 words: state question (is there formal periodic structure in the solar signal, and do sub-annual harmonics exist?), methods (cluster-robust Schuster, MFPA), catalogs tested, key findings
2. **Data Source** — ISC-GEM catalog (n=9,210); G-K mainshocks (n=5,883); A1b mainshocks (n=7,137); note dependency on A4 declustering results; cross-reference A4 suppression findings
3. **Methodology**
   - 3.1 Phase-normalized binning (cite Adhoc A1, `rules/data-handling.md`) — note use for bin-distribution supplement; MFPA and Schuster use event times directly
   - 3.2 Cluster-robust Schuster Spectrum Test (Park et al. 2021): describe modification, clustering threshold, cluster-count effective n
   - 3.3 MFPA (Dutilleul et al. 2015): describe weighted periodogram, bootstrap significance, 300-period scan range
   - 3.4 A1b phase interval cross-reference procedure
   - 3.5 Replication note: standard Schuster not used alone, per Bradley & Hubbard (2024) critique
4. **Results**
   - 4.1 Schuster spectrum: embed Figure 1; state which periods are significant under standard vs cluster-robust test for each catalog; tabulate annual and tidal period p-values
   - 4.2 MFPA: embed Figure 2; list all significant periods detected; state for each whether it matches A1b interval predictions
   - 4.3 Harmonic-interval cross-reference: embed Figure 3; state consistency of detected harmonics with three-interval structure
5. **Cross-Topic Comparison** — compare to Ader & Avouac (2013) annual periodicity at 40% amplitude; compare to Dutilleul et al. (2015) 6-month and 4-month sub-annual harmonics at Parkfield; note similarities and differences at global scale vs regional catalog
6. **Interpretation** — state whether annual signal is confirmed; whether sub-annual harmonics appear; whether the three-interval structure is consistent with any harmonic pattern; maintain objectivity
7. **Limitations** — MFPA bootstrap has 1,000 iterations (may have limited power for marginal signals); cluster-robust test uses 1-day threshold (sensitivity to this choice not tested); ISC-GEM catalog completeness varies by decade
8. **References** — Ader & Avouac (2013), Bradley & Hubbard (2024), Dutilleul et al. (2015), Dutilleul et al. (2021), Park et al. (2021), Adhoc A1, Adhoc A1b, A4

---

## 7. Update context docs

- Append to `topic-a2/.claude/docs/topic-summary.md`:
  ```
  ## Case A1: Schuster Spectrum and MFPA Periodicity Analysis
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [annual period Schuster p-value (standard vs cluster-robust) for raw/G-K/A1b catalogs; MFPA significant periods detected; consistency with A1b three-interval structure]
  ```
- Update Case A1 status from `Pending` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a2/.claude/CLAUDE.md` Case Table
