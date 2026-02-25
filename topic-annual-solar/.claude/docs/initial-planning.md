# Initial Planning: Annual Solar Signal Investigation

## Context

### Prior Result (Case 3A)

A full-population analysis of 9,802 global M >= 6.0 earthquakes over a 72-year period found:

- **Solar cycle signal confirmed**: chi-square χ² = 45.61, p = 6.13 × 10⁻⁵, Cramér's V = 0.0176
- **Bimodal structure**: excess events near equinox-equivalent positions (bins 4 and 15), deficits near solstice positions
- **Effect magnitude**: approximately ±10–14% deviation from uniform expectation
- **Robustness**: zero of 1,000 synthetic uniform catalogs matched the observed extremity (>99.9% confidence)
- **Controls null**: lunar signal p = 0.606, local solar time signal p = 0.368

Reference: https://github.com/jakeYeager/erebus-valid-one/blob/main/approach-three/output/case_3a_whitepaper.md

### Central Puzzle

The bimodal equinox pattern is not straightforwardly explained by the dominant mechanism in the literature (hydrological/snow loading). Hydrological forcing peaks at summer or winter solstice depending on hemisphere and regional climate — not at the equinoxes. The equinox signal is more consistent with a direct solar-geometric forcing (e.g., solar declination rate of change, Earth-Sun distance, or tidal stress geometry) or with a globally symmetric effect that cancels hemisphere-specific hydrological loading and leaves a residual equinox peak. This distinction is the core scientific question.

### Additional Context from Tidal Literature Review

A second literature search on solid-earth tidal triggering and large earthquake periodicity establishes the following context for case evaluation:

- Tidal triggering is detectable globally but the effect **inverts with magnitude**: strongest for small, shallow events; absent or indistinguishable from random at M >= 6 (Métivier 2009, Ide et al. 2016).
- For M >= 8, no calendar-day or lunar-cycle preference has been found in 400+ years of records (Hough 2018).
- Published tidal correlations frequently fail replication on independent time windows (Bradley & Hubbard 2024); aftershock contamination and the ~20% tidal asymmetry bias are primary causes of false positives in Schuster-type tests.
- A suggestive seasonal signal (70% of deep-focus M > 7 events in April–September) exists in one narrow population (Zhan & Shearer 2015) but lacks a physical mechanism and hasn't been independently confirmed.
- The Case 3A finding of a significant equinox signal in a global M >= 6 catalog has **no established precedent** in the literature — which strengthens its novelty and raises the methodological bar.

---

## Proposed Test Cases

Cases are grouped by whether they primarily replicate existing literature or break new ground. Each case is evaluated against both literature reviews.

---

### Group A: Literature Replication and Validation

---

#### Case A1: Schuster Spectrum and MFPA Periodicity Analysis

> **Verdict: REFRAME**

**Original motivation:** Ader & Avouac (2013) applied the Schuster spectrum to Himalayan seismicity and detected annual periodicity with 40% amplitude while detecting no tidal-period signal. Dutilleul et al. (2015) showed MFPA resolved sub-annual harmonics (6-month, 4-month) at Parkfield alongside the annual signal.

**Reframe reason:** Bradley & Hubbard (2024) specifically identified that the standard Schuster test is susceptible to false positives from two sources relevant here: (1) aftershock temporal clustering biases tidal phase distributions, and (2) a ~20% tidal asymmetry in the astronomical forcing itself creates spurious signal. While these critiques target tidal-period tests specifically, the same aftershock contamination risk applies to any periodic test on an unclustered catalog. The standard Schuster spectrum should not be used alone.

**Revised test:** Apply the cluster-robust Schuster Spectrum Test (Park et al. 2021) rather than the standard Schuster. Use MFPA (Dutilleul et al. 2015) as the primary periodicity method, scanning 6 hours to 18 months. Run both on the raw catalog and on GK-declustered catalog (after A4) for comparison.

**Expected outcome if literature holds:** Annual period (365.25 days) reaches significance; tidal periods (12h, 24h, 14 days) do not. Sub-annual harmonics at 6 and 4 months would confirm the bimodal equinox structure as the fundamental signal.

**Key references:** Ader & Avouac (2013), Dutilleul et al. (2015), Park et al. (2021), Bradley & Hubbard (2024)

---

#### Case A2: b-Value Seasonal Variation

> **Verdict: KEEP**

**Motivation:** Colledge et al. (2025) found that in Nepal, a 20-year dataset showed no robust seasonal rate variation but did show annual variation in the Gutenberg-Richter b-value (~0.1 per kPa), peaking during the monsoon loading phase. Additionally, Ide et al. (2016) found the b-value decreases as tidal shear stress amplitude increases for very large events — establishing that b-value responds to stress perturbations even when rate does not.

**Test:** Divide the 72-year catalog into solar-cycle phase bins (matching Case 3A binning). Compute the Gutenberg-Richter b-value independently for each bin using maximum likelihood estimation. Test whether b-value varies significantly across the annual solar cycle. Compute confidence intervals via bootstrap resampling.

**Expected outcome if literature holds:** b-value peaks near solstice positions (clamping phase), troughs near equinox positions — inverse of the rate signal.

**Diagnostic value:** If rate and b-value vary inversely in phase, the equinox excess is partly a magnitude-distribution effect (more events capable of reaching M 6.0 threshold during unclamping), not only a nucleation rate effect. This would be a new finding at global M >= 6 scale.

**Key references:** Colledge et al. (2025), Bettinelli et al. (2008), Ide et al. (2016)

---

#### Case A3: Magnitude Stratification of the Solar Signal

> **Verdict: KEEP — elevated priority**

**Motivation (strengthened):** The tidal literature review substantially reinforces this case. Métivier (2009) found tidal triggering shows inverse magnitude dependence (stronger for smaller events). Ide et al. (2016) found M < 8.2 indistinguishable from random for tidal phase. Hao et al. (2018) found diurnal periodicity stronger at larger magnitudes in Japan. These findings in opposite directions at different timescales make the magnitude-dependence of the annual solar signal a key discriminator: if the Case 3A signal concentrates in a specific magnitude band, that directly constrains the mechanism.

**Test:** Stratify the catalog into magnitude bands (M 6.0–6.4, M 6.5–6.9, M 7.0–7.4, M >= 7.5). Compute chi-square, Rayleigh statistic, and Cramér's V for each band independently. Test whether effect size increases, decreases, or is flat with magnitude.

**Expected outcome if hydrological loading is mechanism:** Effect should weaken or disappear for M >= 7.5.

**Expected outcome if solar-geometric forcing is mechanism:** Effect is flat or increases with magnitude.

**Expected outcome if consistent with tidal literature pattern:** Effect should be weakest at M >= 6 (the opposite of what Case 3A found globally) — making any positive result here a novel finding.

**Key references:** Hao et al. (2018), Johnson et al. (2017), Métivier et al. (2009), Ide et al. (2016)

---

#### Case A4: Declustering Sensitivity Analysis

> **Verdict: KEEP — must execute first**

**Motivation (strengthened):** Bradley & Hubbard (2024) specifically demonstrated that aftershock temporal clustering produces artificial tidal phase correlations, and that this artifact is large enough to explain several published positive results. The same mechanism applies directly to any periodic signal test. Park et al. (2021) showed that the standard Schuster test requires prior declustering and that failing to remove aftershocks produces false positives. This is the single most important methodological validation for Case 3A.

**Test:** Apply two standard declustering algorithms (Gardner-Knopoff 1974 window method; Reasenberg 1985 spatial-temporal linking) to the catalog. Repeat the Case 3A chi-square and Rayleigh analysis on each declustered version. Compare p-values and bin distributions against the undeclustered result.

**Expected outcome:** If the equinox signal is genuine, it persists with comparable p-values across all three catalog versions.

**Diagnostic value:** A signal that disappears after declustering would indicate the equinox structure is driven by temporal clustering of aftershock sequences — an artifact. Given the Bradley & Hubbard findings, this case is a prerequisite for trusting any result downstream.

**Key references:** Park et al. (2021), Bradley & Hubbard (2024), Gardner & Knopoff (1974), Ader & Avouac (2013)

---

### Group B: Novel Investigations

---

#### Case B1: Hemisphere Stratification — Phase Symmetry Test

> **Verdict: KEEP — highest priority novel case**

**Motivation:** The literature attributes seasonal seismicity to hemisphere-specific hydrological loading. If hydrology drives the equinox signal, the Northern and Southern Hemispheres should be approximately 6 months out of phase with each other. If the signal reflects a solar-geometric mechanism symmetric across hemispheres (equinox geometry affects both hemispheres identically), both hemispheres should peak at the same solar phase positions.

**Test:** Split the catalog into Northern Hemisphere (latitude > 0°) and Southern Hemisphere (latitude < 0°) subsets. Compute solar-phase bin distributions independently for each. Test whether the equinox excess appears at the same phase positions in both hemispheres or at positions 0.5 cycles apart.

**Novel contribution:** No published study has applied this specific phase-symmetry test to a global M >= 6.0 catalog to discriminate between hydrological and geometric solar forcing.

**Prediction (geometric hypothesis):** Both hemispheres peak at the same two equinox positions.
**Prediction (hydrological hypothesis):** NH and SH peaks are offset by ~half a year.

---

#### Case B2: Ocean vs. Continent Location — Hydrological Loading Discrimination

> **Verdict: KEEP**

**Motivation:** Hydrological and snow loading mechanisms are land-surface processes. If the equinox signal persists in purely oceanic seismicity, hydrological loading is not a sufficient explanation. The tidal literature adds relevant context: Scholz et al. (2019) found strong tidal modulation at Axial Volcano (mid-ocean ridge), demonstrating that oceanic/volcanic settings can show periodic signals via different mechanisms than continental hydrology. A positive oceanic result would therefore require careful mechanism discrimination.

**Test:** Classify events using a distance-to-coastline threshold (e.g., > 200 km offshore = oceanic, < 50 km or on-continent = continental). Compute the solar-phase bin distribution for each subset. Test whether the equinox signal appears in both, only one, or neither.

**Novel contribution:** No published study has used ocean/continent separation as a mechanism discriminator for annual solar seismicity at global scale.

**Prediction (geometric hypothesis):** Equinox signal appears in both oceanic and continental subsets.
**Prediction (hydrological hypothesis):** Signal is confined to or significantly stronger in continental earthquakes.

---

#### Case B3: Tectonic Regime Stratification

> **Verdict: KEEP — reframed with additional motivation**

**Original motivation:** Different fault geometries respond differently to Coulomb stress changes from surface loading. Thrust faults are clamped; normal faults are unclamped; patterns should be anti-phased if loading drives the signal.

**Reframe addition:** Métivier et al. (2009) found slightly greater tidal triggering in normal and strike-slip regions than thrust regions at the global scale — establishing a fault-mechanism baseline for periodic forcing sensitivity. This provides a reference pattern: if the annual solar signal matches the Métivier tidal pattern (normal > strike-slip > thrust), it suggests a similar stress-change geometry even at annual timescales. If it differs, a distinct mechanism is implied.

**Test:** Classify events by focal mechanism type (thrust, normal, strike-slip) using the GCMT catalog. Compute solar-phase distributions for each class separately. Test for differences in phase, amplitude, and significance.

**Prediction (loading hypothesis):** Thrust faults show solstice deficits; normal faults show solstice excess; patterns approximately anti-phased.
**Prediction (geometric hypothesis):** All three classes show equinox excess with similar phase.
**Prediction (consistent with Métivier tidal pattern):** Normal and strike-slip show stronger signal than thrust.

**Key references:** Métivier et al. (2009), Johnson et al. (2017b)

---

#### Case B4: Depth Stratification — Surface Loading Penetration Test

> **Verdict: REFRAME — now has stronger positive motivation**

**Original motivation:** Hydrological loading primarily affects the shallow crust (< 20 km). Intermediate-depth events should be shielded. If the signal persists at depth, surface loading is insufficient.

**Reframe reason:** Zhan & Shearer (2015) found that 70% of deep-focus M > 7 earthquakes (depth > 500 km) occurred in April–September — a suggestive seasonal signal in deep events with *no known surface mechanism*. This directly inverts the original prediction: the literature now hints that large deep events may show *more* seasonality, not less. This makes depth stratification a much richer test than simply checking whether shallow signal exceeds deep signal.

**Revised test:** Stratify the catalog by focal depth: shallow (0–20 km), mid-crustal (20–70 km), intermediate (70–300 km), and deep (> 300 km if sample permits). Compute solar-phase distributions for each band. Test whether signal amplitude increases, decreases, or is non-monotonic with depth.

**Revised predictions:**
- *Surface loading hypothesis*: Signal strongest at 0–20 km; absent at > 70 km.
- *Geometric/deep forcing hypothesis*: Signal present at all depth ranges including > 300 km.
- *Zhan & Shearer pattern*: Deep events show a seasonality, but potentially with different phase than shallow events — a finding of independent interest.

**Key references:** Zhan & Shearer (2015), Dutilleul et al. (2021), Johnson et al. (2017)

---

#### Case B5: Solar Declination Rate-of-Change vs. Position Test

> **Verdict: KEEP**

**Motivation:** Case 3A encoded seismic events by their position in the solar annual cycle. The bimodal equinox structure could reflect distinct solar variables that both peak at the equinoxes: (a) the rate of change of solar declination is maximum at the equinoxes, and (b) the Earth-Sun distance is near mean value at equinoxes. These variables have different physical implications and have not been distinguished in any published seismicity study.

**Test:** Re-encode each event using three separate phase variables: (1) solar declination angle (ranging –23.5° to +23.5°), (2) rate of change of declination (degrees/day), (3) Earth-Sun distance (AU). Compute bin distributions for each. Identify which variable produces the most coherent, highest-significance clustering.

**Novel contribution:** No published seismicity study has explicitly tested solar declination rate-of-change as a potential trigger variable, nor compared it to solar position and Earth-Sun distance within the same catalog.

---

#### Case B6: Rolling Window Stationarity Test

> **Verdict: KEEP — elevated priority**

**Original motivation:** A 72-year signal could be driven by a concentrated subset of years (non-stationarity).

**Elevation reason:** Bradley & Hubbard (2024) found that tidal correlations detected in pre-2000 data consistently failed to replicate in post-2000 data — the canonical demonstration of non-stationarity as an artifact signature. Dutilleul et al. (2021) found that Parkfield's semiannual periodicity shifted phase after the 2004 M6.0 earthquake. Both findings support stationarity testing as a critical robustness check. A signal that holds stable phase and amplitude across independent decade-long windows is far more credible than one that concentrates in a single era.

**Test:** Compute the Rayleigh statistic and mean phase angle of the solar-cycle distribution using a sliding 10-year window across the 72-year record. Plot the trajectory of p-value, mean vector length, and mean phase. Test whether the equinox phase is stable across decades or drifts.

**Key references:** Bradley & Hubbard (2024), Dutilleul et al. (2021)

---

## Prioritization (Revised)

| Case                              | Verdict | Type        | Mechanism Discriminated         | Priority          |
| --------------------------------- | ------- | ----------- | ------------------------------- | ----------------- |
| A4: Declustering sensitivity      | KEEP    | Replication | Artifact vs. genuine            | **Execute first** |
| B6: Rolling window stationarity   | KEEP ↑  | Novel       | Temporal persistence            | High              |
| A1: Schuster / MFPA spectrum      | REFRAME | Replication | Periodicity confirmation        | High              |
| B1: Hemisphere phase symmetry     | KEEP    | Novel       | Hydrological vs. geometric      | High              |
| A3: Magnitude stratification      | KEEP ↑  | Replication | Loading vs. nucleation scale    | High              |
| B2: Ocean vs. continent           | KEEP    | Novel       | Hydrological vs. geometric      | Medium            |
| B4: Depth stratification          | REFRAME | Novel       | Surface vs. deep forcing        | Medium            |
| A2: b-value seasonality           | KEEP    | Replication | Rate vs. magnitude distribution | Medium            |
| B3: Tectonic regime               | KEEP    | Novel       | Stress geometry                 | Medium            |
| B5: Declination rate vs. position | KEEP    | Novel       | Solar variable identification   | Low               |

**Recommended execution order:** A4 → B6 → A1 → B1 → A3 → B2 → B4 → A2 → B3 → B5

**Rationale for changes from original order:**
- A4 remains first: the Bradley & Hubbard replication failures make declustering validation a prerequisite for all downstream cases.
- B6 moves to second: stationarity is now a near-prerequisite given the documented failure of tidal correlations across time windows; a non-stationary signal would deprioritize the entire remaining program.
- A3 elevated: the tidal literature's clear magnitude-dependence finding makes this a high-value discriminator for the mechanism question.
- B4 and B6 reframed but kept: new literature (Zhan & Shearer 2015; Dutilleul 2021) adds positive motivation rather than just null expectations.
- No cases discarded: the tidal literature strengthens the case for every proposed test by either adding motivation or tightening the predictions.
