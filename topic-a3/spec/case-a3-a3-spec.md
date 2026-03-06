# Case A3.A3: Phase-Concentration Audit

**Status:** Complete

**Source case:** A2.A4
**Informed by:** A3.C2 (sequential removal precedent), A3.B3 (tectonic classification), A3.B4 (depth bands)

**Intent statement:**
A2.A4 established that declustering suppresses the solar-phase signal, and that *which* events are removed matters more than how many. A3.C2 subsequently showed that sequential removal of the 12 largest M≥8.5 sequences produces smooth, diffuse chi-square degradation — no single sequence dominates the signal. A3.A3 extends this to a full catalog-level audit: characterizing whether the signal is concentrated in a small number of identifiable events or is genuinely diffuse.

Critically, the audit is designed to be symmetric across the oscillation. Prior analyses focused on elevated phase bins (equinox peaks). If the mechanism is a suppression/release oscillation, the signal manifests equally in suppressed bins (solstice troughs) as in elevated bins. Auditing only the high-point events introduces high-point bias. A3.A3 characterizes both sides — elevated-bin events and suppressed-bin events — and asks whether removing either population causes the full distribution to regress toward uniformity symmetrically (consistent with a genuine oscillation) or asymmetrically (peaks collapse but troughs remain, or vice versa).

**Empirical bin structure (k=24, full catalog, n=9,210):**
- Expected per bin: 383.8 (1-SD = 19.6)
- **Elevated bins** (z > +1.0): 4, 5, 6, 7, 15, 19, 21 — total ~3,000 events
- **Suppressed bins** (z < −1.0): 2, 8, 10, 11, 12, 13, 16, 18, 22 — total ~3,182 events
- Strongest elevation: bin 15 (phase 0.646, z=+4.30), bin 21 (phase 0.896, z=+3.08)
- Strongest suppression: bin 13 (phase 0.562, z=−2.59), bin 16 (phase 0.688, z=−2.08)
- Note: deepest suppression falls around July–September, not exactly at the June solstice — this asymmetry is documented in the whitepaper

**Sequence catalog sizes:**
- G-K aftershocks: 3,327 events (parent_id column present)
- Reasenberg aftershocks: 945 events
- A1b aftershocks: 2,073 events
- G-K mainshocks: 5,883; Reasenberg mainshocks: 8,265; A1b mainshocks: 7,137

**Data context block:**

| File | Path | n | Key columns |
|------|------|---|-------------|
| ISC-GEM full catalog | `data/iscgem/iscgem_global_6-9_1950-2021.csv` | 9,210 | `usgs_id`, `usgs_mag`, `event_at`, `solar_secs`, `latitude`, `longitude`, `depth` |
| GSHHG classification | `data/iscgem/plate-location/ocean_class_gshhg_global.csv` | 9,210 | `usgs_id`, `ocean_class`, `dist_to_coast_km` |
| G-K mainshocks | `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv` | 5,883 | `usgs_id`, `aftershock_count`, `foreshock_count` |
| G-K aftershocks | `data/iscgem/declustering-algorithm/aftershocks_gk-seq_global.csv` | 3,327 | `usgs_id`, `parent_id`, `parent_magnitude` |
| Reasenberg mainshocks | `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv` | 8,265 | `usgs_id`, `aftershock_count` |
| Reasenberg aftershocks | `data/iscgem/declustering-algorithm/aftershocks_reas-seq_global.csv` | 945 | `usgs_id`, `parent_id`, `parent_magnitude` |
| A1b mainshocks | `data/iscgem/declustering-algorithm/mainshocks_a1b-seq_global.csv` | 7,137 | `usgs_id`, `aftershock_count` |
| A1b aftershocks | `data/iscgem/declustering-algorithm/aftershocks_a1b-seq_global.csv` | 2,073 | `usgs_id`, `parent_id` |

**Phase normalization:** `phase = (solar_secs / 31_557_600.0) % 1.0`

**Tectonic class boundaries (GSHHG, A3.B3 baseline):**
- `continental`: dist_to_coast_km ≤ 50
- `transitional`: 50 < dist_to_coast_km ≤ 200
- `oceanic`: dist_to_coast_km > 200

**Depth bands (A3.B4):** shallow (<20 km), mid-crustal (20–70 km), intermediate (70–300 km), deep (≥300 km)

**Script path conventions:**
- `BASE_DIR = Path(__file__).resolve().parent.parent` — resolves to `topic-a3/`
- All output paths: `BASE_DIR / "output" / ...`
- All data paths: `BASE_DIR.parent / "data" / "iscgem" / ...`

**Planned outputs:**
- `src/case-a3-a3-analysis.py` — influence computation, permutation baseline, sequential removal curves, representativeness tests; writes results JSON
- `src/visualization-case-a3-a3.py` — all PNG figures
- `tests/test-case-a3-a3.py` — test suite (all tests must pass)
- `output/case-a3-a3-results.json`
- `output/case-a3-a3-whitepaper.md`
- `output/case-a3-a3-influence-distribution.png` — Figure 1: signed influence histogram with sequence annotation
- `output/case-a3-a3-degradation-curves.png` — Figure 2: two parallel removal curves vs. baselines
- `output/case-a3-a3-permutation-tails.png` — Figure 3: per-bin counts with permutation confidence bands
- `output/case-a3-a3-representativeness.png` — Figure 4: tectonic × depth profile of top-influence event groups
- `output/case-a3-a3-sequence-annotation.png` — Figure 5: sequence membership breakdown by bin type

---

## 1. Environment and data loading

In `src/case-a3-a3-analysis.py`:

- Imports: `pandas`, `numpy`, `scipy.stats`, `pathlib`, `json`, `logging`
- Set `BASE_DIR = Path(__file__).resolve().parent.parent`
- Define paths for all 8 data files (see data context block)
- Load full catalog; assert n=9210; compute `phase = (solar_secs / 31_557_600.0) % 1.0`
- Load GSHHG classification; assert n=9210; merge onto full catalog; assert no NaN in `ocean_class`
- Build sequence membership lookup from all three algorithm aftershock files:
  - Collect all `usgs_id` values that appear as `parent_id` in any aftershock file → "mainshock_with_aftershocks" set
  - Collect all `usgs_id` values that appear as event rows in any aftershock file → "aftershock_member" set
  - For each event in full catalog: classify as:
    - `"aftershock"` if in any aftershock member set
    - `"mainshock_with_sequence"` if in mainshock_with_aftershocks set (and not already aftershock)
    - `"isolated_mainshock"` otherwise
  - Record counts of each class; log
- Add `sequence_role` column to full catalog DataFrame
- Add `tectonic_class` from `ocean_class` column
- Add `depth_band`: `"shallow"` (<20), `"mid_crustal"` (20–70), `"intermediate"` (70–300), `"deep"` (≥300)
- Define constants:
  ```python
  K = 24
  JULIAN_YEAR_SECS = 31_557_600.0
  N_PERMUTATIONS = 1000
  N_RANDOM_BASELINES = 500
  ELEVATED_BINS = [4, 5, 6, 7, 15, 19, 21]   # z > +1.0 in full catalog
  SUPPRESSED_BINS = [2, 8, 10, 11, 12, 13, 16, 18, 22]  # z < -1.0
  A1B_INTERVALS = {"interval_1": [4, 5], "interval_2": [15], "interval_3": [21]}
  TOP_N_REPRESENTATIVENESS = 100
  ```

---

## 2. Signed chi-square influence

```python
def compute_signed_influence(phases: np.ndarray, k: int = 24) -> np.ndarray:
    """
    Compute each event's signed chi-square influence: the change in total chi-square
    when that event is removed.

    Parameters
    ----------
    phases : np.ndarray
        Phase values in [0, 1) for all events (length n).
    k : int
        Number of bins.

    Returns
    -------
    np.ndarray of length n. Positive = event is in an elevated bin (removal
    reduces chi-square). Negative = event is in a suppressed bin (removal
    increases chi-square, deepening the trough).

    Notes
    -----
    All events in the same bin receive identical influence values.
    Influence for an event in bin j:
        influence_j = (2 * (obs_j - expected) - 1) / expected
    where expected = n / k. This is the exact first-order change in chi-square.
    """
```

Implementation:
```python
n = len(phases)
bin_indices = np.floor(phases * k).astype(int) % k
obs = np.bincount(bin_indices, minlength=k).astype(float)
expected = n / k
# Per-bin influence value
bin_influence = (2.0 * (obs - expected) - 1.0) / expected
# Map each event to its bin's influence
event_influence = bin_influence[bin_indices]
return event_influence
```

Record the full influence array. Attach to the main DataFrame as column `chi2_influence`.

Bin-level summary (record in results JSON under `"bin_summary"`):
- For each of k=24 bins: obs, expected, z-score, bin_influence, bin_type (`"elevated"` / `"suppressed"` / `"neutral"`), count of events per sequence_role in that bin

---

## 3. Permutation baseline

```python
def run_permutation_baseline(
    phases: np.ndarray,
    k: int = 24,
    n_permutations: int = 1000,
    rng_seed: int = 42,
) -> dict:
    """
    Permute solar phases and compute per-bin count distributions.

    Parameters
    ----------
    phases : np.ndarray
        Observed phase values (length n).
    k : int
        Number of bins.
    n_permutations : int
        Number of random phase shuffles.
    rng_seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict with:
        "perm_bin_counts": np.ndarray shape (n_permutations, k)
        "perm_chi2": np.ndarray shape (n_permutations,)
        "bin_p5": np.ndarray shape (k,) — 5th percentile per bin
        "bin_p95": np.ndarray shape (k,) — 95th percentile per bin
        "sig_elevated_bins": list of bin indices where actual obs > bin_p95
        "sig_suppressed_bins": list of bin indices where actual obs < bin_p5
    """
```

Implementation:
- `rng = np.random.default_rng(rng_seed)`
- For each permutation: shuffle phases, bin, record bin counts and chi-square
- Compute 5th and 95th percentiles across permutations for each bin
- Compare actual obs to permutation bands; classify significantly elevated/suppressed bins

This sub-test formally establishes which solstice bins are suppressed beyond chance expectation — the first direct statistical test of the suppression side in this project.

---

## 4. Sequential removal curves

```python
def sequential_removal_curve(
    df: pd.DataFrame,
    removal_order: list[str],
    k: int = 24,
    track_bins: list[int] | None = None,
) -> list[dict]:
    """
    Sequentially remove events and track chi-square and bin z-scores at each step.

    Parameters
    ----------
    df : pd.DataFrame
        Full catalog with 'phase' and 'usgs_id' columns.
    removal_order : list[str]
        usgs_id values in the order they should be removed.
    k : int
        Number of bins.
    track_bins : list[int] | None
        Specific bins to track z-scores for. If None, tracks ELEVATED_BINS and SUPPRESSED_BINS.

    Returns
    -------
    List of dicts, one per removal step, each containing:
        step, n_removed, n_remaining, chi2, p_chi2,
        elevated_bin_zscores (dict bin→z), suppressed_bin_zscores (dict bin→z),
        pct_catalog_removed.
    """
```

**Removal order construction:**

*Elevated-bin removal order:* Sort bins by initial influence (descending). For each bin in order, include all events in that bin (order within bin: by usgs_id for determinism). Bins in descending influence: 15, 21, 6, 7, 5, 19, 4 (based on actual z-scores; confirm from data at runtime).

*Suppressed-bin removal order:* Sort bins by |negative influence| (descending). Bins in descending |influence|: 13, 16, 2, 12, 8, 11, 22, 18, 10 (confirm at runtime).

**Baseline curves:**

*G-K sequential baseline:* Use G-K aftershock catalog ordered by `|delta_t_sec|` ascending (temporally closest to parent first). Remove up to min(len(gk_aftershocks), len(removal_order)) events.

*Random removal baseline:* Draw 500 random removal sequences of `len(removal_order)` events from the full catalog (without replacement, seeded). Run curve for each; compute mean chi-square at each step and 10th–90th percentile band.

Run both the elevated-bin and suppressed-bin curves to completion (all events in their respective bin groups removed). Record signal persistence count: first step where p ≥ 0.05.

---

## 5. Representativeness test

```python
def representativeness_test(
    df: pd.DataFrame,
    top_n: int = 100,
) -> dict:
    """
    Test whether top-influence events are anomalous relative to the full catalog
    and the signal-bearing stratum.

    Parameters
    ----------
    df : pd.DataFrame
        Full catalog with chi2_influence, tectonic_class, depth_band columns.
    top_n : int
        Number of top positive and top |negative| influence events to test.

    Returns
    -------
    dict with representativeness results for four groups:
        top_positive_50, top_positive_100, top_negative_50, top_negative_100.
    """
```

For each group (top-50 and top-100 by positive influence; top-50 and top-100 by |negative| influence):
1. Extract group; compute tectonic class distribution and depth band distribution
2. Chi-square goodness-of-fit vs. full catalog distribution: `chi2, p = scipy.stats.chisquare(group_counts, full_counts * (top_n / n_full))`
3. Chi-square goodness-of-fit vs. signal-bearing stratum (continental + mid-crustal subset of full catalog)
4. Record: is the group significantly different from full catalog? From signal-bearing stratum?
5. Also compute sequence_role breakdown for the group

A group that is NOT significantly different from the signal-bearing stratum (p > 0.05 vs. that stratum) is "representative of the signal population" — direct refutation of the rogue-events framing.

---

## 6. Results JSON structure

```json
{
  "case": "A3.A3",
  "title": "Phase-Concentration Audit",
  "parameters": {
    "n_catalog": 9210,
    "k": 24,
    "julian_year_secs": 31557600.0,
    "elevated_bins": [4, 5, 6, 7, 15, 19, 21],
    "suppressed_bins": [2, 8, 10, 11, 12, 13, 16, 18, 22],
    "n_permutations": 1000,
    "n_random_baselines": 500,
    "top_n_representativeness": 100,
    "rng_seed": 42
  },
  "sequence_roles": {
    "n_isolated_mainshock": int,
    "n_mainshock_with_sequence": int,
    "n_aftershock": int
  },
  "bin_summary": [
    {"bin": int, "obs": int, "expected": float, "z": float,
     "bin_influence": float, "bin_type": str,
     "n_isolated": int, "n_mainshock_seq": int, "n_aftershock": int}
  ],
  "permutation_baseline": {
    "n_permutations": 1000,
    "bin_p5": [float, ...],
    "bin_p95": [float, ...],
    "sig_elevated_bins_permutation": [int, ...],
    "sig_suppressed_bins_permutation": [int, ...],
    "perm_chi2_p95": float
  },
  "degradation_curves": {
    "elevated_removal": {
      "n_events_in_group": int,
      "signal_persistence_step": int,
      "signal_persistence_pct_catalog": float,
      "curve": [{"step": int, "n_removed": int, "pct_catalog_removed": float,
                 "chi2": float, "p_chi2": float, "elevated_bin_zscores": {}, "suppressed_bin_zscores": {}}]
    },
    "suppressed_removal": { /* same structure */ },
    "gk_sequential_baseline": { /* same structure, up to min removal count */ },
    "random_baseline": {
      "n_draws": 500,
      "curve_mean_chi2": [float, ...],
      "curve_p10_chi2": [float, ...],
      "curve_p90_chi2": [float, ...]
    }
  },
  "representativeness": {
    "top_positive_50":  {"n": 50, "tectonic_dist": {}, "depth_dist": {}, "sequence_role_dist": {}, "p_vs_full": float, "p_vs_signal_stratum": float, "representative_of_signal_stratum": bool},
    "top_positive_100": { /* same */ },
    "top_negative_50":  { /* same */ },
    "top_negative_100": { /* same */ }
  },
  "summary": {
    "signal_diffuse": bool,
    "elevated_persistence_pct": float,
    "suppressed_persistence_pct": float,
    "degradation_symmetric": bool,
    "top100_elevated_representative_of_signal_stratum": bool,
    "top100_suppressed_representative_of_signal_stratum": bool
  }
}
```

---

## 7. Visualizations

In `src/visualization-case-a3-a3.py`:

Import: `matplotlib`, `matplotlib.pyplot`, `numpy`, `json`, `pathlib`, `logging`

---

**Figure 1 — Signed influence distribution** (`output/case-a3-a3-influence-distribution.png`):
- 1-row × 3-column grid; one panel per sequence role (isolated mainshock, mainshock with sequence, aftershock)
- Each panel: histogram of chi2_influence values for that role's events; x-axis = influence value (positive right, negative left); steelblue bars; vertical dashed line at 0
- Color the bars: bins with influence > 0 in coral (elevated-bin events), bins with influence < 0 in steelblue (suppressed-bin events), bins near 0 in gray
- Annotate each panel: n events, fraction in elevated bins, fraction in suppressed bins
- Shared x-axis range across panels
- Title: "Signed Chi-Square Influence by Sequence Role (A3.A3)"
- 300 DPI, figsize=(14, 5)

---

**Figure 2 — Parallel degradation curves** (`output/case-a3-a3-degradation-curves.png`):
- 2-row × 2-column grid
- **Top row** — Elevated-bin removal (left) and Suppressed-bin removal (right):
  - Main line: chi-square as events removed (x-axis = % catalog removed)
  - Overlaid: G-K sequential baseline (dashed orange), random baseline mean ± 10–90th percentile band (light gray fill)
  - Horizontal dashed red line at critical chi-square value corresponding to p=0.05
  - Mark signal persistence step (first crossing of p=0.05) with vertical dotted line and annotation
- **Bottom row** — Z-score tracking for elevated-bin removal (left) and suppressed-bin removal (right):
  - Multiple lines: one per A1b interval bin (bins 4–5 aggregate, bin 15, bin 21) in warm colors; one line for strongest suppressed bin (bin 13) in cool color
  - Horizontal dashed line at z=+1.0 and z=−1.0
  - Shows whether elevated z-scores and suppressed z-scores converge toward zero simultaneously (symmetric) or independently (asymmetric)
- All panels share x-axis (% catalog removed, 0 to max removal)
- Title: "Phase-Concentration Audit — Sequential Removal (A3.A3)"
- 300 DPI, figsize=(14, 10)

---

**Figure 3 — Permutation confidence bands** (`output/case-a3-a3-permutation-tails.png`):
- Single panel; x-axis = bin index 0–23; y-axis = event count
- Steelblue bars: actual bin counts
- Gray shaded band: permutation 5th–95th percentile range per bin
- Mark bins where actual count exceeds p95 (significantly elevated) with red asterisk above bar
- Mark bins where actual count falls below p5 (significantly suppressed) with blue asterisk below bar
- Annotate with approximate calendar months on secondary x-axis (same mapping as A3.B2 phase-curves)
- Dashed horizontal line at expected count (n/k)
- Title: "Per-Bin Counts vs. Permutation Confidence Bands (A3.A3, n=1,000 permutations)"
- 300 DPI, figsize=(13, 5)

---

**Figure 4 — Representativeness heatmap** (`output/case-a3-a3-representativeness.png`):
- 2-row × 2-column grid: rows = tectonic class distribution and depth band distribution; columns = top-100 positive-influence group and top-100 negative-influence group
- Each panel: grouped bar chart comparing the group's distribution (steelblue) to the full catalog distribution (light gray) and signal-bearing stratum (continental + mid-crustal, dashed outline)
- Annotate each panel with chi-square goodness-of-fit p-value vs. full catalog and vs. signal stratum
- Label: "ns" (p>0.05), "*" (0.01<p≤0.05), "**" (p≤0.01)
- Title: "Tectonic and Depth Profile of Top-Influence Event Groups (A3.A3)"
- 300 DPI, figsize=(13, 9)

---

**Figure 5 — Sequence membership by bin type** (`output/case-a3-a3-sequence-annotation.png`):
- Single stacked horizontal bar chart; y-axis = bin type (Elevated A1b Interval 1, Elevated A1b Interval 2, Elevated A1b Interval 3, Other Elevated, Suppressed Strong [z<-2], Suppressed Moderate [-2≤z<-1], Neutral)
- x-axis = fraction of events; stacked segments = sequence role (isolated mainshock, mainshock with sequence, aftershock) — colors: isolated=steelblue, mainshock_seq=darkorange, aftershock=gray
- Annotate each bar with total event count (n=...)
- Reference vertical line at full-catalog sequence role fractions (so deviation from catalog-wide proportions is visible)
- Title: "Sequence Role by Bin Type — Phase-Concentration Audit (A3.A3)"
- 300 DPI, figsize=(11, 7)

---

## 8. Test suite

In `tests/test-case-a3-a3.py`:

- `test_catalog_load`: assert n=9210; assert phase in [0.0, 1.0); assert `sequence_role` column added with no NaN
- `test_gshhg_merge`: assert no NaN in `ocean_class` and `tectonic_class` after merge
- `test_sequence_role_partition`: assert n_isolated + n_mainshock_seq + n_aftershock == 9210
- `test_sequence_catalog_sizes`: assert GK aftershocks n=3327; Reasenberg aftershocks n=945; A1b aftershocks n=2073
- `test_influence_formula`: for a synthetic array with known bin counts, verify that `compute_signed_influence` returns (2*(obs_j - expected) - 1) / expected for each event in bin j; verify events in the same bin get identical values
- `test_influence_length`: assert `chi2_influence` column has length 9210
- `test_elevated_suppressed_bins`: assert ELEVATED_BINS and SUPPRESSED_BINS are non-empty; assert they are disjoint; assert all indices in [0, 23]
- `test_permutation_shape`: assert permutation produces arrays of shape (1000, 24) for bin counts and (1000,) for chi2
- `test_permutation_bands`: assert `bin_p5` and `bin_p95` are arrays of length 24; assert all p5 ≤ actual obs ≤ p95 is NOT required but assert p5 < p95 for all bins
- `test_permutation_sig_suppressed_identified`: assert at least one bin is classified as significantly suppressed in permutation test (given the actual distribution, bins 13 and 16 should fail the permutation threshold)
- `test_degradation_curves_present`: assert both `elevated_removal` and `suppressed_removal` curves present in results JSON with at least 1 step each
- `test_signal_persistence_recorded`: assert `signal_persistence_step` and `signal_persistence_pct_catalog` are present and non-null for both removal curves
- `test_representativeness_groups`: assert all four representativeness groups present; assert each has `p_vs_full` and `p_vs_signal_stratum` float fields and `representative_of_signal_stratum` bool
- `test_summary_keys`: assert `summary` dict present with all required keys; assert `signal_diffuse` is bool; assert `degradation_symmetric` is bool
- `test_output_figures_exist`: assert all 5 PNG files exist and size > 50 KB

---

## 9. Whitepaper

In `output/case-a3-a3-whitepaper.md`:

Standard header (Author: Jake Yeager, Version: 1.0, Date: current date) and footer (Generated with Claude Code [model name]).

### Sections:

1. **Abstract** (150–200 words): state that A2.A4 and A3.C2 established that the signal is not dominated by large sequences; A3.A3 extends this to a full catalog-level audit using signed chi-square influence; describe the symmetric design (both elevated and suppressed bins tested); state whether the signal is diffuse or concentrated; state whether degradation is symmetric or asymmetric; state the representativeness finding for top-influence events.

2. **Data Source**: ISC-GEM full catalog (n=9,210); GSHHG tectonic classification; G-K, Reasenberg, and A1b sequence catalogs for sequence role annotation (3,327 / 945 / 2,073 aftershocks respectively).

3. **Methodology**
   - 3.1 Phase normalization: Julian year constant
   - 3.2 Signed chi-square influence: describe the exact formula (2*(obs_j − expected) − 1) / expected; note that all events in the same bin receive identical influence values; explain positive (elevated-bin) vs. negative (suppressed-bin) sign
   - 3.3 Sequence role annotation: describe non-circular use of sequence membership (based on spatial-temporal proximity windows, not solar phase)
   - 3.4 Permutation baseline: 1,000 shuffles; per-bin 5th–95th percentile confidence bands; significance classification for both elevated and suppressed tails
   - 3.5 Sequential removal curves: describe two parallel curves (elevated-bin and suppressed-bin removal orders); describe G-K sequential and random baselines; define signal persistence count (events-to-collapse as % of catalog); define symmetric degradation criterion (both elevated and suppressed z-scores converge toward zero at similar removal rates)
   - 3.6 Representativeness test: chi-square goodness-of-fit for top-influence groups vs. full catalog and signal-bearing stratum; describe "representative_of_signal_stratum" classification

4. **Results**
   - 4.1 Bin structure and influence distribution: embed Figure 1; describe distribution of influence values across sequence roles; state whether aftershock members are over-represented in elevated bins relative to their catalog share
   - 4.2 Permutation confidence bands: embed Figure 3; state which bins are significantly elevated and significantly suppressed beyond permutation expectation; note the actual location of deepest suppression (bins 13–16, approximately July–September) relative to the June solstice
   - 4.3 Sequential removal: embed Figure 2; report signal persistence counts for elevated-bin and suppressed-bin removal separately; compare degradation rates to G-K sequential and random baselines; state whether degradation is symmetric (both bin-type z-scores converge together) or asymmetric
   - 4.4 Representativeness: embed Figure 4; state whether top-100 elevated-bin and top-100 suppressed-bin events are representative of the signal-bearing stratum (continental + mid-crustal); state whether they are anomalous relative to the full catalog
   - 4.5 Sequence annotation by bin type: embed Figure 5; compare sequence role fractions in elevated, suppressed, and neutral bins to the catalog-wide fractions

5. **Cross-Topic Comparison**
   - **Aftershock Phase-Preference Analysis (A2.A4):** A2.A4 found that aftershock populations show *stronger* solar-phase clustering than mainshock catalogs. A3.A3's sequence role annotation directly tests whether elevated-bin overrepresentation of aftershock members is consistent with this finding.
   - **Targeted Major Sequence Phased Declustering Test (A3.C2):** A3.C2 established diffuse degradation under sequential removal of the 12 largest sequences. A3.A3 generalizes this: if the signal is diffuse at the sequence level (C2), it should also be diffuse at the individual-event influence level (A3.A3). Signal persistence count >5% would confirm this generalization.
   - **Corrected Null-Distribution Geometric Variable Test (A3.B5):** A3.B5 found `declination_rate` as the top-ranked variable, with alignment at Interval 1 (March equinox). If the phase-concentrated events are geographically and tectonically representative of the signal-bearing stratum, this is consistent with a physical mechanism acting broadly — not a few anomalous events driving a spurious correlation.

6. **Interpretation**: state the primary finding on signal diffuseness objectively; address the symmetric vs. asymmetric degradation question and its implication for the suppression/release oscillation model; note that representative top-influence events (matching the signal-bearing stratum profile) strengthen the case for a genuine physical mechanism, while anomalous top-influence events would warrant further investigation; guard against over-interpreting the absence of rogue events as definitive proof of a physical mechanism.

7. **Limitations**: all events in the same bin receive identical influence values — the audit characterizes *bins* rather than truly individual events; the sequential removal order within a bin is deterministic by usgs_id, not by any meaningful criterion; the permutation baseline preserves catalog size but does not preserve geographic or temporal structure; the representativeness test uses a chi-square goodness-of-fit which may have reduced power for small top-N groups; the signal persistence count is sensitive to which specific events are in the highest-influence bins.

8. **References**
   - Yeager, J. (2026). A2.A4: Aftershock Phase-Preference Analysis. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.C2: Targeted Major Sequence Phased Declustering Test. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B3: Ocean/Coast Sequential Threshold Sensitivity. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation. erebus-vee-two internal report.
   - Yeager, J. (2026). A3.B5: Corrected Null-Distribution Geometric Variable Test. erebus-vee-two internal report.

---

## 10. Update context docs

After all outputs are generated and tests pass:

- Append to `topic-a3/docs/topic-summary.md`:
  ```
  ## Case A3.A3: Phase-Concentration Audit
  **Status:** [Complete | Blocked | Abandoned]
  **Key results:** [signal diffuse T/F; elevated-bin persistence pct; suppressed-bin persistence pct; degradation symmetric T/F; top-100 elevated events representative of signal stratum T/F; top-100 suppressed events representative T/F; which bins are permutation-significant on suppressed tail]
  ```
- Update Case A3 status from `Planning` to `Complete` (or `Blocked`/`Abandoned`) in `topic-a3/CLAUDE.md` Case Table
