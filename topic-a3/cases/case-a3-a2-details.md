# A3.A2: Stratified Schuster/MFPA Periodicity Audit

**Reference case:** A2.A1

*(Title updated from "Aftershock Periodicity Analysis (Schuster/MFPA)" — aftershock focus dropped; case is now a stratified periodicity audit informed by signal-location and oscillation findings from Tier 2.)*


## Instructions

1. Read the source case specifications file to understand how the previous case was constructed.
2. Read the source case results file to understand what results were gained from the analysis.
3. Read the expectations for this current case A3.A2 and use the knowledge from the source case to inform your approach to satisfy the expectations of this case.
4. Build a specification file for this case at topic-a3/spec/case-a3-a2-spec.md. It should use the same structure and pattern of deliverables as the source case specification file. Use the data sources for A3.A2 for analysis. The spec file should encapsulate everything required for an agent to execute the specifications outside of the main context thread.


## Reference case A2.A1 summary/files

**Results file**: `@topic-a2/output/case-a1-results.json`

Key A2.A1 findings relevant to this reframe:
- Schuster spectrum identified the ~75.6-day quarter-year period as the strongest robust MFPA detection
- The ~75.6-day period aligns with Interval 1 (~March equinox, phase ~0.19–0.25)
- Cluster correction inflated detections 329× before correction (cluster-window sensitivity unresolved)
- Analysis was run exclusively on the raw full catalog mainshocks — no stratification by tectonic class, depth, or hemisphere


## Reframing rationale

A3.A2 was originally planned as an aftershock population periodicity analysis, gated on A3.A1 (Aftershock Phase-Preference Characterization). A3.A1 was abandoned because A3.A3 and A3.C2 collectively resolved the concentration-vs-diffuse and aftershock phase representation questions that motivated it.

However, A3.A2's core framework — Schuster spectrum and MFPA — has unique value independent of the aftershock focus. The Tier 2 cases produced new stratification information that A2.A1 never exploited, and A3.A3 introduced the oscillation framing that generates testable periodicity predictions. Specifically:

- **A3.B3/B4** identified the signal-bearing stratum (continental × mid-crustal). A2.A1's detections may have been diluted by the non-signal-bearing oceanic and deep populations included in the full catalog.
- **A3.B2** found a ~5-month NH/SH phase offset rather than a clean 6-month anti-phase. This predicts that NH and SH strata should show different dominant phase angles for any half-year period detected in the Schuster spectrum.
- **A3.B5** identified `declination_rate` as the top-ranked geometric variable, with its peak aligning precisely with Interval 1 (~DOY 80). The MFPA framework has never tested whether the ~75.6-day detection is strongest in the declination-rate-aligned (signal-bearing) stratum.
- **A3.A3** formally characterized a symmetric oscillation: both equinox elevations and solstice suppressions are statistically significant beyond permutation expectation. This predicts that both the ~91.25-day quarter-year period (one equinox peak) and the ~182.5-day half-year period (full peak-to-trough oscillation) should be detectable — and the relative strength of these two periods directly tests whether the mechanism drives one equinox event per year or a full suppression/release cycle.

The Interval 1 contradiction from A2.A1 also remains unresolved: MFPA detects the ~75.6-day period, yet Interval 1 does not survive post-declustering at k=24. These two analyses were never directly reconciled. Comparing full-catalog vs. G-K mainshock Schuster spectra in the periodicity domain would determine whether Interval 1's declustering sensitivity is also visible in the periodicity domain.


## Expectations for A3.A2

**Intent:** Apply the Schuster spectrum and MFPA framework from A2.A1 to stratified subsets of the full catalog to test:
1. Whether the ~75.6-day quarter-year period is a global catalog feature or stratum-specific (diluted by non-signal-bearing populations)
2. Whether the ~182.5-day half-year period is detectable alongside the quarter-year period, and which is stronger — the oscillation framing from A3.A3 predicts both should be present
3. Whether NH and SH strata show different phase angles for the half-year period — a direct, testable prediction from A3.B2's ~5-month offset finding
4. Whether the Interval 1 contradiction is visible in the periodicity domain: does the ~75.6-day period survive G-K declustering in the Schuster spectrum?

This is not a replication of A2.A1. It uses A2.A1's methods as a tool to probe the stratification structure and oscillation geometry revealed by Tier 2 cases.


## Sub-tests for A3.A2

**Sub-test 1 — Baseline Schuster spectrum replication and structural improvements:**

Replicate the A2.A1 Schuster spectrum on the full catalog (n=9,210) with structural improvements applied:
- Increase bootstrap replicates to at least 10,000 for stable significance thresholds
- Apply FDR correction (Benjamini-Hochberg) post-hoc to all period-level p-values
- Test cluster-window sensitivity: run at 3-day and 7-day inter-event cluster correction windows; compare detection counts and the inflation factor against A2.A1's 329× result
- Record the top-10 detected periods, their Schuster p-values, FDR-corrected significance, and whether each aligns with: (a) the quarter-year period ~75.6 days, (b) the half-year period ~182.5 days, (c) the full-year period ~365.25 days

**Sub-test 2 — Stratified spectrum: signal-bearing vs. non-signal-bearing:**

Run the Schuster spectrum independently on:
- Signal-bearing stratum: continental (dist_to_coast_km ≤ 50) AND mid-crustal (20–70 km depth) — n approximately 1,500–2,000
- Non-signal-bearing remainder: all other events

Compare the top detected periods across strata. If the ~75.6-day period is stronger in the signal-bearing stratum than in the remainder, this confirms it is a genuine feature of the signal-carrying population rather than a catalog-wide artifact.

Apply adaptive k rule to bin counts if chi-square supplementary tests are run: k=24 if n≥500, k=16 if 200≤n<500, k=12 if 100≤n<200.

**Sub-test 3 — NH/SH phase angle comparison:**

Run the Schuster spectrum independently on NH (latitude ≥ 0) and SH (latitude < 0) strata. For the half-year period (~182.5 days), extract the dominant phase angle from each stratum's Schuster statistic. Compute the phase difference between NH and SH dominant angles.

A3.B2 found a ~5-month (~152-day) peak offset rather than the 6-month (182.5-day) anti-phase that hemisphere-specific hydrological loading would predict. In the periodicity domain, this predicts that NH and SH phase angles for the half-year period should differ by approximately 5/6 of a cycle (~300° rather than 180°). Report the observed phase difference and whether it is consistent with the ~5-month offset.

**Sub-test 4 — Interval 1 contradiction: declustering sensitivity in the periodicity domain:**

Run the Schuster spectrum on three catalogs: (a) full catalog, (b) G-K mainshock catalog, (c) G-K aftershock catalog. For each, record whether the ~75.6-day period is detected at FDR-corrected significance.

A2.A4 showed that Interval 1 (March equinox) does not survive G-K declustering in the chi-square domain. If the ~75.6-day Schuster detection also fails to survive declustering, this reconciles the contradiction: the period is present in the full catalog but is carried by aftershock members, not mainshocks. If the period survives declustering, the contradiction deepens and the mechanism must act independently of sequence membership.


## Design decisions (resolved)

- **No aftershock focus as primary driver**: the aftershock population comparison (original A3.A2 intent) is dropped; aftershocks appear only in Sub-test 4 as one of three comparison catalogs for the Interval 1 contradiction test.
- **Stratification informed by completed cases**: signal-bearing stratum (B3+B4), NH/SH split (B2), and oscillation period predictions (A3.A3) are all imported from Tier 2 results — no new stratification design decisions needed.
- **Both quarter-year and half-year periods tested**: the oscillation framing from A3.A3 generates the half-year period prediction. Both must be tested; relative strength is the primary interpretive question.
- **Structural improvements applied to all sub-tests**: FDR correction, 10,000 bootstrap replicates, and cluster-window sensitivity (3-day vs. 7-day) are applied uniformly.


## Open questions (resolved)

- ~~What is the aftershock focus?~~ → Dropped; case is stratified by signal-location and hemisphere, not by sequence membership.
- ~~How many bootstrap replicates?~~ → 10,000 for stable significance thresholds.
- ~~FDR method?~~ → Benjamini-Hochberg; applied post-hoc to all period-level p-values.
- ~~3-day vs. 7-day cluster correction sensitivity?~~ → Both tested in Sub-test 1; inflation factor compared to A2.A1's 329× result.
- ~~What period to test for the oscillation?~~ → Both ~91.25-day (quarter-year, one equinox peak) and ~182.5-day (half-year, full oscillation cycle); relative strength is the primary interpretive question.


## Data sources for A3.A2

- Full catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv`
- GSHHG tectonic classification: `data/iscgem/plate-location/ocean_class_gshhg_global.csv`
- Full catalog mainshocks with G-K window: `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv`
- Full catalog aftershocks with G-K window: `data/iscgem/declustering-algorithm/aftershocks_gk-seq_global.csv`

---

## Revision Note

**Date:** March 5, 2026
**Changes from original planning-refinement.md entry:**
- Title updated from "Aftershock Periodicity Analysis (Schuster/MFPA)" to "Stratified Schuster/MFPA Periodicity Audit"
- Aftershock population as primary focus dropped; A3.A1 dependency removed
- Four sub-tests defined: (1) baseline replication with structural improvements, (2) signal-bearing vs. non-signal-bearing stratification, (3) NH/SH phase angle comparison, (4) Interval 1 contradiction declustering sensitivity
- Oscillation period predictions from A3.A3 (half-year period alongside quarter-year) added as explicit sub-test targets
- NH/SH phase angle comparison added as direct test of A3.B2's ~5-month offset prediction
- Structural improvements (10,000 bootstrap replicates, FDR correction, cluster-window sensitivity) resolved from open questions
