# A3.B2: Hemisphere Stratification Refinement

**Source cases:** A2.B1


## Instructions

1. Read the reference case specifications file to understand how to structure a spec file
2. Read the reference case results file to understand how to structure a `json` results file
3. Read the expectations for this current case A3.B2, as well as the topic summary file
4. Create a specifications file for this case at topic-a3/spec/case-a3-b2-spec.md :
   1. It should satisfy the expectations of this case
   2. Use the data sources for A3.B2 for analysis
   3. It should incapsulate everthing required for an agent to execute the specifications outside of the main context thread


## Source case files

**Specifications file:** @topic-a2/spec/case-b1-spec.md
**Results file:** @topic-a2/output/case-b1-results.json
**Topic Summary:** @topic-a3/docs/topic-summary.md (informative results from A3.B3 (geographic) and A3.B4 (depth) for where signal is strongest)


## Expectations for A3.B2

**Gap or concern (original):**
A2.B1 classified Interval 1 (March equinox, phase ~0.19–0.25) as absent in the SH at 33% interval overlap — just below the 50% threshold used as the primary classification criterion. This is a threshold-sensitive result. The planning-initial question asks only about declustered population testing; the console summary gap identified the 50% threshold as an analytic choice that has not been tested for sensitivity.

**Extended concern (from A3.B3 and A3.B4):**
A raw NH vs. SH comparison is confounded by tectonic composition. A3.B3 found that the solar-phase signal is concentrated in continental (p=0.0005) and transitional (p=0.0003) tectonic classes, while oceanic events are borderline non-significant (p=0.061) until subduction-proximal events are excluded at T=175 km. A3.B4 found that the mid-crustal band (20–70 km) carries essentially the entire signal (p=4×10⁻⁹) while shallow (<20 km) and intermediate/deep bands are non-significant. The NH holds a disproportionate share of active subduction systems (Japan, Cascadia, Alaska, Aegean, Himalaya) and continental landmass, while the SH is more oceanic. A direct NH/SH χ² comparison therefore partially recovers this tectonic and depth imbalance rather than isolating a hemispheric solar-forcing difference. The revised intent addresses this confound directly.

**Intent:** Four-component case:

1. **Tectonic-matched hemisphere comparison:** Within each tectonic class from A3.B3 (continental, transitional), split events by hemisphere (NH: latitude ≥ 0°; SH: latitude < 0°) and run χ²(k=24) independently. Report for oceanic as well but flag as borderline-confidence given A3.B3 results. This tests whether any NH/SH signal asymmetry persists after controlling for tectonic composition, or whether it is an artifact of the hemispheres having different tectonic makeups.

2. **Mid-crustal hemisphere split:** Within the mid-crustal band (20–70 km) — the signal-bearing depth from A3.B4 — split by hemisphere and run χ²(k=24). Compare signal strength (χ², Cramér's V) and peak-bin phase between NH and SH. This directly asks whether the dominant depth-signal is hemispherically symmetric (consistent with a solar-geometric mechanism) or concentrated in one hemisphere (consistent with a hemisphere-specific loading effect). Use adaptive k: k=24 if n≥500, k=16 if 200≤n<500, k=12 if 100≤n<200; flag n<100 as low-n.

3. **Phase alignment comparison:** For all cells in sub-tests (1) and (2) that achieve χ² significance (p<0.05), extract the peak bin and compute the corresponding solar phase. Compare NH peak phase vs. SH peak phase. If both hemispheres peak near the equinox phases (0.19–0.25 and 0.69–0.75), the signal is globally symmetric — strengthening the case for a solar-geometric mechanism. A ~0.5 offset between hemispheres would indicate a loading-driven signal with opposite seasonal peaks.

4. **Revised Interval 1 SH threshold sensitivity:** Repeat the A2.B1 Interval 1 SH overlap analysis at four overlap thresholds (33%, 40%, 45%, 50%) on (a) all SH events and (b) SH continental-only events. If Interval 1 SH absence survives in the continental-only subset, it is a genuine hemispheric phase asymmetry. If it disappears (continental SH events *do* show Interval 1), then the A2.B1 result was driven by SH oceanic composition rather than a true hemispheric difference. This recontextualizes the A2.B1 threshold-sensitivity result through the tectonic lens of A3.B3.

**Dropped from original plan:** The standalone "declustered population repeat" is removed as a primary sub-test. A2.A4 already established that declustering suppresses the signal; repeating that on the NH/SH split would be predictably flat. Instead, declustering is incorporated as a sensitivity layer on sub-test (1): run the tectonic-matched hemisphere comparison on full catalog, G-K mainshocks, and Reasenberg mainshocks, and report suppression rates. This preserves the declustering context without making it the primary question.


## Data sources for A3.B2

- Full catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv`
- GSHHG ocean/tectonic classification: `data/iscgem/plate-location/ocean_class_gshhg_global.csv` — `ocean_class`, `dist_to_coast_km` (tectonic class join; matches A3.B3)
- Full catalog declustered with G-K window: `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv`
- Full catalog declustered with Reasenberg algorithm: `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv`

---

## Revision Note

**Date:** March 4, 2026
**Reason:** Original two-component plan (threshold sensitivity + declustered repeat) was assessed as informationally flat given A3.B3 and A3.B4 findings. The NH/SH split is confounded by the hemispheres having materially different tectonic compositions: the NH is more continental and subduction-rich, which is precisely where the signal concentrates. Revised to four components that control for tectonic class and depth before interpreting hemispheric differences. The phase alignment comparison (sub-test 3) is new and directly diagnostic for the solar-geometric vs. loading-mechanism distinction.
