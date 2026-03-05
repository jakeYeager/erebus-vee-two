> **Status: Approved**

# Topic A3 Planning Refinement

## Overview

This document proposes cases for Topic A3 based on a review of all Topic A2 console summaries and the initial planning review (`planning-initial.md`). Cases are organized in three series:

- **A-series** — Analysis refinements: new analytical approaches applied to the existing framework
- **B-series** — Stratification refinements: revisiting A2 stratification cases with corrected methods or additional dimensions
- **C-series** — Synthesis: cross-case tests that resolve confounds identified across multiple A2 cases

Each case includes its source case reference(s), the gap or concern it addresses, and open data/methodology questions to resolve before case spec writing.

---

## Case Table

| Series | Case ID | Title | Source Case(s) | Origin |
| ------ | ------- | ----- | -------------- | ------ |
| A | A3.A1 | Aftershock Phase-Preference Characterization | A2.A4 Sub-C | Evaluation gap |
| A | A3.A2 | Aftershock Periodicity Analysis (Schuster/MFPA) | A2.A1, A2.A4 | Evaluation gap + planning-initial |
| A | A3.A3 | Phase-Aware Declustering Methodology | A2.A4 note | Console summary gap |
| B | A3.B1 | Rolling-Window Chi-Square Repeat | A2.B6 | Evaluation gap |
| B | A3.B2 | Hemisphere Stratification Refinement | A2.B1 | Console summary gap + planning-initial |
| B | A3.B3 | Ocean/Coast Sequential Threshold Sensitivity | A2.B2 | Evaluation gap + planning-initial |
| B | A3.B4 | Depth × Magnitude Two-Way Stratification with Moho Isolation | A2.B4 | Console summary gap + planning-initial |
| B | A3.B5 | Corrected Null-Distribution Geometric Variable Test | A2.B5 | Evaluation gap |
| C | A3.C1 | Subduction Zone Subset Test | A2.B2, A2.B4 | Console summary gap |
| C | A3.C2 | Targeted Major Sequence Phased Declustering Test | A2.B6, A2.A4 | Evaluation gap |

---

## Run Order

Cases have upstream dependencies that determine execution order. Cases with no dependencies can run in parallel; synthesis cases require their upstream cases to complete first.

### Dependency Map

| Case | Feeds | Reason |
| ---- | ----- | ------ |
| A3.B1 | A3.C2 | A3.C2's major-sequence removal targets are sharpened by A3.B1's interval-level temporal tracking (which years and intervals are non-stationary) |
| A3.C2 | A3.A1 | A3.A1 benefits from A3.C2's sequence identification and handling results before characterizing aftershock phase preferences |
| A3.A1 | A3.A2 | A3.A2 needs A3.A1's aftershock phase distribution before MFPA results are interpretable |
| A3.B3 | A3.C1 | A3.C1's subduction zone boundary definition is informed by where the signal migrates in A3.B3's threshold sweep |
| A3.B4 | A3.C1 | A3.C1 needs to know whether the mid-crustal signal is independently real or magnitude-driven before interpreting a subduction zone subset |

### Execution Tiers

**Tier 1 — Foundational (run first)**
- **A3.B1** — Interval-level stationarity tracking; provides interpretive context for B2/B3/B4 and gates A3.C2

**Tier 2 — After Tier 1; independent cases may run alongside A3.C2**
- **A3.C2** — After A3.B1; sequence identification and handling results inform A3.A1
- **A3.B3** — Standalone; needed for A3.C1
- **A3.B4** — Standalone; needed for A3.C1
- **A3.B2** — Standalone; A3.B1 results provide interpretive context but do not block it
- **A3.B5** — Standalone methodological correction
- **A3.A3** — Exploratory probe; does not block any other case in the current list

**Tier 3 — After Tier 2 completes**
- **A3.A1** — After A3.C2; most novel A2 finding, run with sequence handling context from A3.C2; gates A3.A2
- **A3.C1** — After A3.B3 + A3.B4

**Tier 4 — After Tier 3 completes**
- **A3.A2** — After A3.A1

> **Note:** Unlike A2, where A2.A4 (declustering sensitivity) served as a single-gating prerequisite for all downstream cases, A3 has a sequential foundational chain: A3.B1 → A3.C2 → A3.A1, with independent stratification cases (B2–B5, A3) running in parallel during Tier 2.

---

## A-Series: Analysis Refinements

### A3.A1: Aftershock Phase-Preference Characterization

**Source:** A2.A4 Sub-C

**Gap or concern:**
A2.A4's Sub-analysis C is the most novel finding from all of Topic A2 — all three aftershock populations show *stronger* solar-phase clustering than their corresponding mainshock populations (χ²=66–149, all p<10⁻⁵; Reasenberg aftershocks χ²=107 vs. mainshocks χ²=40). This finding directly contradicts the premise that declustering removes solar-phase contamination; the aftershock sequences are themselves non-uniformly distributed in solar phase. This deserves a dedicated case, not the parenthetical treatment it received in A2.A4's framing.

**Intent:** Characterize the aftershock phase distributions in detail:
- What interval structure do aftershocks show at k=24? Does it replicate the A1b baseline three-interval structure or show a shifted structure?
- Are the elevated aftershock χ² values concentrated in specific large event sequences (e.g., 2004 Sumatra M9.1 — identified in A2.B6 as the likely driver of the 2003–2014 significant window cluster), or is the aftershock phase preference diffuse across sequences?
- Does the aftershock phase preference change with sequence age (early vs. late aftershocks within a defined window)?

**Data source requirements:** Sequence-enriched aftershock catalogs (G-K, Reasenberg, A1b) are **required** — `parent_id`, `parent_magnitude`, `delta_t_sec`, and `delta_dist_km` are now available on all three algorithm datasets. Corresponding sequence-enriched mainshock catalogs provide `aftershock_count` and `window_secs` per mainshock row, which are needed to filter sequences large enough for meaningful subdivision and to compute the early/late temporal midpoint.

**Open questions:**
- What temporal window defines "early" vs. "late" aftershocks within a sequence? *Decision: Split aftershocks into first/second half of sequence train based on temporal duration of train (`window_secs / 2` as midpoint), not total event count. Note: a single late straggler extends `window_secs` and shifts the midpoint outward — carry as a sensitivity check comparing duration-based vs. count-based splits.*

---

### A3.A2: Aftershock Periodicity Analysis (Schuster/MFPA)

**Source:** A2.A1, A2.A4; see also planning-initial A2.A1 comparitive questions

**Gap or concern:**
A2.A1 ran the Schuster spectrum and MFPA analysis exclusively on mainshock catalogs. Given A2.A4 Sub-C's finding that aftershock populations carry stronger solar-phase signals, running this analysis on aftershock populations is a direct cross-case validation that A2.A1 cannot currently provide. Two unresolved tensions also emerge from A2.A1:

1. **Interval 1 contradiction:** A2.A1 identifies the ~75.6-day quarter-year period as the strongest robust MFPA detection, and this aligns with Interval 1 (March equinox). But A2.A4.4.2 shows Interval 1 does not survive post-declustering at k=24 under any method. These two analyses need direct comparative treatment.

2. **Structural improvements pending from planning-initial:** More bootstrap replicates for better resolution; threshold sensitivity tests at 3-day and 7-day cluster correction windows; Bonferroni or FDR correction applied post-hoc.

**Intent:** Run the full A2.A1 framework (Schuster spectrum + MFPA + cluster correction) on aftershock populations from all three declustering methods. Compare resulting periodicity structure to the mainshock results to determine whether the 75.6-day quarter-year period is an aftershock artifact, a mainshock artifact, or genuinely present in both populations. Apply structural improvements (bootstrap replicates, threshold sensitivity, multiple comparison correction) to the full analysis.

**Data source requirements:** Sequence-enriched aftershock catalogs (G-K, Reasenberg, A1b) are **suggested** — `parent_id` enables correct stratification of the aftershock population by parent event, and `delta_t_sec` enables temporal position within sequence as an optional secondary variable (e.g., testing whether the periodicity signal shifts between early and late aftershocks within sequences).

**Open questions:**
- How many additional bootstrap replicates are warranted? (Current: assess from A2.A1 spec)
- Should FDR correction use Benjamini-Hochberg or a more conservative method?
- Does the 3-day vs. 7-day cluster correction threshold materially change the 329× inflation finding from A2.A1?

---

### A3.A3: Phase-Aware Declustering Methodology

**Source:** A2.A4 note (console summary)

**Gap or concern:**
A2.A4's console summary included a methodological finding not captured in the formal sub-analyses: *"the A1b-vs-G-K suppression direction was reversed in the data (A1b suppressed slightly more than G-K despite removing fewer events), indicating that which events are removed matters more than how many."* Standard declustering methods (G-K, Reasenberg, A1b) are defined by spatial-temporal proximity windows — they have no relationship to the solar-phase values of the events they remove. Because the signal's presence depends more on event identity than event count, a targeted approach could be more informative.

**Intent:** Design and apply a phase-aware declustering probe: identify events that are disproportionately contributing to solar-phase bin elevations (i.e., phase-clustered events relative to a uniform null) and measure the signal before and after their removal. This is not proposed as a declustering standard — it is an exploratory methodological test to determine how much of the signal can be attributed to a small number of phase-concentrating events.

**Data source requirements:** Sequence-enriched aftershock catalogs (G-K, Reasenberg, A1b) are **suggested** — `parent_id` is now available on all three algorithm datasets, enabling sequence membership to be used as a structuring variable when identifying which events are phase-clustering candidates, reducing (but not eliminating) circular reasoning risk. Permutation framework still required to define phase-clustered thresholds.

**Open questions:**
- How to define "phase-clustered" events without circular reasoning (i.e., without selecting events based on the outcome being tested)? A leave-one-out jackknife or permutation threshold may be required.
- What comparison framework makes the result interpretable against G-K/Reasenberg baselines?

---

## B-Series: Stratification Refinements

### A3.B1: Rolling-Window Chi-Square Repeat

**Source:** A2.B6

**Gap or concern:**
A2.B6 used the Rayleigh test (unimodal) as its primary stationarity statistic. Chi-square (k=24) was run as secondary and reached significance in 71% of windows vs. Rayleigh's 38.7% — a 32-percentage-point divergence. The B6 console summary explicitly flagged this as "multi-modal within-window structure that the unimodal Rayleigh test misses, worth carrying forward into downstream cases." The stationarity conclusion in A2.B6 is based on Rayleigh; re-running with chi-square as primary may materially change the stationarity classification and its implications for all downstream mechanistic cases.

**Intent:** Repeat the A2.B6 rolling-window stationarity analysis with chi-square (k=24) as the primary statistic and Rayleigh as secondary. Add interval-level tracking within each window — which of the three A1b baseline intervals (Interval 1 ~0.19–0.25, Interval 2 ~0.625–0.656, Interval 3 ~0.875–0.917) are elevated — to determine whether the non-stationarity is global or interval-specific.

**Data source requirements:** Full ISC-GEM catalog with solar phase values. Sequence-enriched mainshock catalogs (G-K, Reasenberg, A1b) are **suggested** — `aftershock_count` per mainshock row enables sequence density to be computed per rolling window, providing a direct diagnostic check: windows with elevated chi-square can be tested for elevated sequence density, which would directly quantify the contribution of major aftershock sequences (e.g., 2004 Sumatra) to window-level significance.

**Open questions:**
- Should window size and stride remain unchanged from A2.B6 (5-year window, 1-year stride, 62 windows)? *Decision: apply what was used A2.B6 unless a change can materially benefit the chi-square tests*
- Should the stationarity threshold (70%) be retained or recalibrated for chi-square as primary? *Decision: apply what was used A2.B6 unless a change can materially benefit the chi-square tests*
- Should this run on declustered catalogs in addition to the raw catalog? *Decision: yes*

---

### A3.B2: Hemisphere Stratification Refinement

**Source:** A2.B1; see also planning-initial A2.B1 question

**Gap or concern (merged from two sources):**
A2.B1 classified Interval 1 (March equinox, phase ~0.19–0.25) as absent in the SH at 33% interval overlap — just below the 50% threshold used as the primary classification criterion. This is a threshold-sensitive result. The planning-initial question asks only about declustered population testing; the console summary gap identified the 50% threshold as an analytic choice that has not been tested for sensitivity. These two concerns are combined here.

**Intent:** Two-component case:

1. **Threshold sensitivity sub-test:** Repeat the A2.B1 Interval 1 SH overlap analysis at four overlap thresholds (33%, 40%, 45%, 50%) to establish where the classification boundary falls. A hard result would show SH Interval 1 absence across all thresholds; a soft result would flip at ≤40%. The classification boundary has direct bearing on the bilateral vs. hemisphere-specific mechanism interpretation.

2. **Declustered population repeat:** Run the full A2.B1 framework (NH/SH χ² at k=24, three-interval structure, half-cycle offset test) on all three declustered mainshock catalogs. Given A2.A4's finding that declustering substantially suppresses the signal, the bilateral structure of Intervals 2 and 3 may not survive declustering, which would weaken the argument against hemisphere-specific mechanisms.

**Data source requirements:** Full ISC-GEM and declustered mainshock catalogs (G-K, Reasenberg, A1b) with hemisphere labels and solar phase values. Existing data sufficient.

---

### A3.B3: Ocean/Coast Sequential Threshold Sensitivity

**Source:** A2.B2; see also planning-initial A2.B2 questions

**Gap or concern:**
A2.B2's result is highly sensitive to the coastline classification boundary. The oceanic subset (GSHHG primary) just misses significance at p=0.061, yet its Cramér's V (0.0276) is *larger* than the continental V (0.0245) — the significance difference is entirely a sample-size artifact, not an effect-size finding. When PB2002's broader oceanic class (absorbing back-arc and marginal-basin events) is used, the oceanic result flips to significant (p=5.24×10⁻³). The transitional zone (50–200 km offshore) is the most significant class under the primary GSHHG classification. The conclusion rests entirely on a boundary drawn at 50 km offshore, and this boundary has not been systematically tested.

**Intent:** Sequentially tighten the transitional zone boundary from 200 km offshore down to 0 km in defined incremental steps, recomputing oceanic/transitional/continental class sizes and chi-square statistics at each step. Identify the specific boundary at which the signal migrates from the transitional class into the oceanic class — or fails to — which directly quantifies the geographic extent of the mechanism. Additionally, assess whether GCMT-matched events can provide a data-driven tectonic setting attribute to supplement or replace the GSHHG proximity metric for the transitional zone definition.

**Data source requirements:** GSHHG ocean classification file, PB2002 classification file, focal mechanism join file, event lat/lon. Will require distance recalculation per event at each boundary increment. Python mapping libraries may be needed if not already present in the environment.

**Open questions:**
- What increment step for threshold tightening? (Suggest 25 km steps from 200 km → 0 km, yielding 8 test points.)
- Does the GCMT focal mechanism dataset (52.9% match rate) have sufficient geographic coverage to define subduction zone proximity as a tectonic setting metric?

---

### A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation

**Source:** A2.B4; see also planning-initial A2.B4 questions

**Gap or concern (merged from two sources):**
A2.B4's console summary explicitly flagged a confound: *"Magnitude and depth are correlated (larger events tend to be deeper), so the depth pattern partly reflects the magnitude trend from A3."* The mid-crustal (20–70 km) band dominates the global signal (χ²=85.48, p=4.02×10⁻⁹), but this band also contains the highest concentration of large-magnitude events. The depth result may be largely the A2.A3 magnitude result expressed through depth. The planning-initial question asks about Moho isolation — this is incorporated as a second sub-analysis.

**Intent:** Two-component case:

1. **Two-way stratification sub-test:** Partition events by both depth band and magnitude band simultaneously. For each depth band, test whether the solar-phase signal persists after holding magnitude roughly constant (e.g., only M 6.0–6.9 events within each depth band). This disentangles the depth and magnitude effects and determines whether the mid-crustal signal is independently real or magnitude-driven.

2. **Moho isolation sub-test:** The Mohorovičić discontinuity (~20–35 km depth in continental crust, ~7–10 km in oceanic crust) separates crustal from mantle behavior. Test a narrow band around the continental Moho (~20–35 km) as a discrete depth class to determine whether the signal is specifically concentrated at the crust-mantle transition.

**Data source requirements:** Full ISC-GEM catalog with depth and magnitude. Moho depth per event would ideally use a tectonic setting or receiver function dataset to assign continental vs. oceanic Moho depth, but a fixed 25 km depth cutoff can serve as a proxy.

**Open questions:**
- Are two-way stratification cell sizes sufficient at k=24? Some depth × magnitude combinations may have very few events. Reducing to k=12 may be required for smaller cells.
- Does existing data include per-event tectonic setting for Moho depth assignment, or will a fixed depth band be used as proxy?

---

### A3.B5: Corrected Null-Distribution Geometric Variable Test

**Source:** A2.B5

**Gap or concern:**
A2.B5 found that geometric variable Cramér's V values (declination_rate V=0.1668, earth_sun_distance V=0.1279, solar_declination V=0.1249) are ~7–9× larger than solar_phase V (0.0181). The B5 console summary identifies this as "largely a null-distribution artifact": events naturally accumulate at geometric extrema (solstices, perihelion) because the Sun moves slowly there. A uniform-in-time seismicity catalog would still produce non-uniform geometric bins. The B5 V values cannot be compared to solar_phase V without correcting for this expected non-uniformity, and extending the B5 framework (e.g., with declustered populations) without first applying this correction would compound the interpretive problem.

**Intent:** Derive the expected geometric variable bin distribution under a null hypothesis of uniform seismicity in time. Two approaches:

1. **Monte Carlo simulation:** Generate N synthetic uniform-in-time catalogs with the same event count as the ISC-GEM catalog and compute chi-square and Cramér's V for each geometric variable. The resulting null distribution defines corrected significance thresholds.

2. **Analytical integration:** Integrate each geometric variable over the solar year to produce the expected proportion of time spent in each bin under uniform seismicity. Use these expected proportions as the null in the chi-square computation (replacing the implicit uniform null).

After null correction, re-evaluate whether the three A1b intervals align with corrected geometric variable elevations. Only after this correction is established does extending the framework to declustered populations yield interpretable results.

**Data source requirements:** Solar geometry file, event timestamps, geometric variable columns. Simulation requires no new data.

**Open questions:**
- Which null derivation approach is preferred: Monte Carlo simulation or analytical integration? Analytical is more elegant but requires numerical integration of the geometric variable's annual function; Monte Carlo is more straightforward and makes fewer assumptions.
- How many Monte Carlo draws are sufficient for stable null distribution estimates?

---

## C-Series: Synthesis

### A3.C1: Subduction Zone Subset Test

**Source:** A2.B2, A2.B4 (console summary gap)

**Gap or concern:**
Both A2.B2 and A2.B4 converge on the same geographic confound from different analytical directions. A2.B2 identifies the transitional zone (50–200 km offshore) as the most significant class and names "subduction zone geometry as a confounding variable." A2.B4 states the mid-crustal dominance "may trace to subduction zone seismicity — the same geometry implicated by B2's significant transitional zone." A2.B3 adds a third convergent line: the unmatched events (pre-1976, n=4,336) carry a stronger signal than any mechanism class, and the thrust class (most consistent with subduction) comes closest to significance among matched events (p=0.084). Three independent cases point to subduction zones; no case has tested them directly.

**Intent:** Define a subduction zone event subset using one or more classification approaches and test the solar-phase signal within that subset vs. the non-subduction remainder of the catalog. If the signal is concentrated in the subduction zone subset, the mechanistic question shifts: not "solar forcing vs. hydrological loading globally" but "why are subduction zones specifically susceptible to solar-phase modulation?"

Classification approach options (in order of preference):
1. Depth + transitional zone intersection as proxy (existing data, no new source required)
2. GCMT thrust mechanism as proxy (available at 52.9% coverage via focal join file)
3. Proximity to Slab2 or equivalent subduction zone boundary dataset (new data source required)

**Data source requirements:** Full ISC-GEM catalog with location, depth, and solar phase. GCMT focal join file. Potentially Slab2 dataset for highest-quality subduction classification.

**Open questions:**
- Which classification approach is most feasible given existing data coverage? The 52.9% GCMT match rate may be a material limitation for the thrust-mechanism proxy.
- If a new dataset (Slab2) is introduced, what are the pipeline implications for event-to-subduction-zone distance computation?

---

### A3.C2: Targeted Major Sequence Phased Declustering Test

**Source:** A2.B6, A2.A4

**Gap or concern:**
A2.B6 identified that the most statistically significant rolling windows cluster in the 2003–2014 start-year range, contemporaneous with aftershock sequences from the 2004 Sumatra M9.1. A2.A4 showed that aftershock populations carry stronger solar-phase signals than mainshocks. Together, these findings raise a direct question the A2 framework did not test: is the global solar-phase signal diffuse across the catalog, or is it disproportionately driven by a small number of major event sequences? If removing the top 3–5 major sequences causes the global signal to collapse, the phenomenon is sequence-concentrated. If the signal survives, it is genuinely diffuse.

**Intent:** Sequentially remove the largest events in the ISC-GEM catalog and their associated aftershock sequences, recomputing the global solar-phase chi-square and interval structure after each removal:

1. Identify the 5 largest events in the catalog by magnitude (candidates: 1960 Chile M9.5, 1964 Alaska M9.2, 2004 Sumatra M9.1, 2011 Tohoku M9.0, 2010 Chile M8.8 or similar).
2. Remove each event and its aftershock sequence using a defined spatial-temporal window.
3. Recompute global chi-square, interval structure, and Cramér's V after each removal.
4. Assess whether the global signal degrades progressively, collapses at a specific removal, or remains robust.

This test should run on both the raw catalog and the declustered mainshock-only catalog, since A2.A4 showed the mainshock signal is also suppressed. If even the mainshock signal is concentrated in a few sequences, the declustering framing needs further revision.

**Report Requirement:** In the whitepaper "Results" block include a simple breakout to describe the major events removed, and basic metrics of their sequence train including: how many fore/after/complete, sequence window duration - mainshock to last aftershock, pre/post count of "last event" half-life (will inform subsequent "early" and "late" designation).

**Data source requirements:** Sequence-enriched mainshock and aftershock catalogs (G-K, Reasenberg) are **required**. `parent_id` on aftershock rows provides precise sequence membership for targeted removal without relying on a fixed spatial/temporal buffer. `foreshock_count`, `aftershock_count`, `window_secs`, and `window_km` on mainshock rows directly supply the sequence train metrics required by the Report Requirement block above, without additional derivation.


**Open questions:**
- Should aftershock windows for removal be defined using the existing G-K/Reasenberg declustering windows, or a fixed spatial/temporal buffer applied to each large event independently? *Decision: use G-K/Reasenberg declustering windows as typical in other cases*
- What magnitude threshold defines a "major sequence" for inclusion in the removal set? (Suggest M≥8.5 as a starting point, yielding approximately 10–15 events in the ISC-GEM catalog.) *Decision: Wait for the results of A3.B1 to inform, or use suggest M≥8.5 as fallback.*
