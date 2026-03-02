> **Status: Superseded.** See `topic-a2/.claude/CLAUDE.md` for current case status.

# Topic A2 Initial Planning: Annual Solar Signal Investigation

## Context

### Prior Result (Case 3A)

> **Catalog note:** Case 3A used the ComCat catalog (9,802 events). All cases in this topic use the ISC-GEM catalog (`data/global-sets/iscgem_global_events.csv`, 9,210 events, 1950–2021). See Cases A0 and A0b in topic-adhoc for a full description of differences between the two catalogs.

A full-population analysis of 9,802 global M >= 6.0 earthquakes over a 72-year period found:

- **Solar cycle signal confirmed**: chi-square χ² = 45.61, p = 6.13 × 10⁻⁵, Cramér's V = 0.0176
- **Bimodal structure**: excess events near equinox-equivalent positions (bins 4 and 15), deficits near solstice positions
- **Effect magnitude**: approximately ±10–14% deviation from uniform expectation
- **Robustness**: zero of 1,000 synthetic uniform catalogs matched the observed extremity (>99.9% confidence)
- **Controls null**: lunar signal p = 0.606, local solar time signal p = 0.368

Reference: https://github.com/jakeYeager/erebus-valid-one/blob/main/approach-three/output/case_3a_whitepaper.md

### Signal Characterization Update (topic-adhoc Cases A1 and A1b)

The Case 3A "bimodal equinox" characterization has been refined by the topic-adhoc binning sensitivity and elevated-bin analyses:

- **Adhoc A1**: The `solar_secs` chi-square is robust across bin counts (Bonferroni-significant at k=16, 24, 32). Cramér's V stable at 0.016–0.018. `lunar_secs` is non-significant at all bin counts after phase normalization — the legacy result was a binning artifact. `midnight_secs` behaves as a negative control at global scale.

- **Adhoc A1b**: The elevated bins resolve into **three phase intervals**, not two:
  1. Phase [0.1875, 0.25] → **~March 10 – April 1** (centered on March equinox)
  2. Phase [0.625, 0.656] → **~August 16 – August 28** (~1 month before September equinox)
  3. Phase [0.875, 0.917] → **~November 16 – December 1** (~1 month before December solstice)

  Only interval 1 cleanly aligns with an equinox. Interval 2 leads the September equinox by ~1 month. Interval 3 is near the December solstice, not the September equinox. This three-interval structure replaces the bimodal equinox characterization from Case 3A.

  The elevated-bin events (n=1,438, 15.6% of catalog) show strong temporal clustering (median IEI 1.46 days vs null CI 11.4–12.8 days) consistent with residual aftershock sequences, mild spatial over-clustering (median NN 35.1 km vs null CI 37.8–45.3 km), and no magnitude skew relative to the full catalog. A1b proposed a data-informed declustering reference window of 83.2 km / 95.6 days — far tighter than G-K's 295-day temporal window at M6.0.

### Central Puzzle (Revised)

The three-interval solar phase pattern is not straightforwardly explained by the dominant mechanism in the literature (hydrological/snow loading). Hydrological forcing peaks at summer or winter solstice depending on hemisphere and regional climate — not at the observed phase positions. The pattern is also not a clean bimodal equinox signal: only one of three intervals aligns with an equinox. The signal structure is more consistent with a complex solar-geometric forcing, a globally symmetric effect that partially cancels hemisphere-specific loading, or a combination of genuine phase preference and residual aftershock sequence contamination at specific solar phases. Distinguishing between these possibilities is the core scientific question.

### Additional Context from Tidal Literature Review

A second literature search on solid-earth tidal triggering and large earthquake periodicity establishes the following context for case evaluation:

- Tidal triggering is detectable globally but the effect **inverts with magnitude**: strongest for small, shallow events; absent or indistinguishable from random at M >= 6 (Métivier 2009, Ide et al. 2016).
- For M >= 8, no calendar-day or lunar-cycle preference has been found in 400+ years of records (Hough 2018).
- Published tidal correlations frequently fail replication on independent time windows (Bradley & Hubbard 2024); aftershock contamination and the ~20% tidal asymmetry bias are primary causes of false positives in Schuster-type tests.
- A suggestive seasonal signal (70% of deep-focus M > 7 events in April–September) exists in one narrow population (Zhan & Shearer 2015) but lacks a physical mechanism and hasn't been independently confirmed.
- The Case 3A finding of a significant solar signal in a global M >= 6 catalog has **no established precedent** in the literature — which strengthens its novelty and raises the methodological bar.

### Data Quality Context (topic-adhoc Cases A0 and A0b)

- ComCat contains 592 more events than ISC-GEM (+6.4%), with 30.3% of ComCat records being ISC-GEM-sourced backfill concentrated pre-1980.
- ComCat has 77.5% 1-decimal magnitude precision vs ISC-GEM's 82.9% 2-decimal. The M 6.0 bin in ComCat is inflated +52% relative to ISC-GEM due to rounding artifacts.
- Within-ComCat duplication is negligible (1 candidate pair). The 592-event gap decomposes into 1,389 ComCat-only events (92.5% in M 6.0–6.4) minus 797 ISC-GEM-only events.
- G-K declustering on ComCat carries quantifiable data quality uncertainty via M 6.0 rounding artifacts as spurious aftershock candidates (A0b Mechanism 2). This concern is mitigated for this topic since all cases use ISC-GEM.

### Cross-Cutting Methodology Standard

All cases in this topic must use **phase-normalized binning** for any analysis involving binned astronomical metrics (`solar_secs`, `lunar_secs`, `midnight_secs`), per `rules/data-handling.md`. Phase normalization maps each event's metric value to [0, 1) of its respective cycle before binning, eliminating period-length edge effects. Established in Approach Three Case 1B.3.1; validated across bin counts in Adhoc Case A1.

---

## Proposed Test Cases

Cases are grouped by whether they primarily replicate existing literature or break new ground. Each case is evaluated against both literature reviews and the topic-adhoc findings.

---

### Group A: Literature Replication and Validation

---

#### Case A1: Schuster Spectrum and MFPA Periodicity Analysis

> **Verdict: REFRAME** (updated with Adhoc A1b cross-reference)

**Original motivation:** Ader & Avouac (2013) applied the Schuster spectrum to Himalayan seismicity and detected annual periodicity with 40% amplitude while detecting no tidal-period signal. Dutilleul et al. (2015) showed MFPA resolved sub-annual harmonics (6-month, 4-month) at Parkfield alongside the annual signal.

**Reframe reason:** Bradley & Hubbard (2024) specifically identified that the standard Schuster test is susceptible to false positives from two sources relevant here: (1) aftershock temporal clustering biases tidal phase distributions, and (2) a ~20% tidal asymmetry in the astronomical forcing itself creates spurious signal. While these critiques target tidal-period tests specifically, the same aftershock contamination risk applies to any periodic test on an unclustered catalog. The standard Schuster spectrum should not be used alone.

**Revised test:** Apply the cluster-robust Schuster Spectrum Test (Park et al. 2021) rather than the standard Schuster. Use MFPA (Dutilleul et al. 2015) as the primary periodicity method, scanning 6 hours to 18 months. Run on three catalog versions: (a) raw catalog, (b) G-K declustered (after A4), and (c) A1b-informed custom-window declustered (after A4). This three-way comparison tests whether G-K over-suppresses relative to the data-informed window.

**Expected outcome if literature holds:** Annual period (365.25 days) reaches significance; tidal periods (12h, 24h, 14 days) do not. Sub-annual harmonics should be cross-referenced against the three A1b phase intervals to test consistency: a 6-month harmonic would explain intervals 1 and 3 (~6 months apart at phases 0.19 and 0.90); the August interval (phase 0.64) does not fit a 6-month harmonic and may indicate a more complex periodicity or residual sequence contamination at that phase.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210)
- G-K mainshocks (after A4): `data/iscgem/declustering-algorithm/mainshocks_G-K_global.csv` (n=5,883)
- A1b mainshocks (after A4): `data/iscgem/declustering-algorithm/mainshocks_a1b_global.csv` (n=7,137)

**Key references:** Ader & Avouac (2013), Dutilleul et al. (2015), Park et al. (2021), Bradley & Hubbard (2024), Adhoc Case A1b

---

#### Case A2: b-Value Seasonal Variation

> **Verdict: KEEP**

**Motivation:** Colledge et al. (2025) found that in Nepal, a 20-year dataset showed no robust seasonal rate variation but did show annual variation in the Gutenberg-Richter b-value (~0.1 per kPa), peaking during the monsoon loading phase. Additionally, Ide et al. (2016) found the b-value decreases as tidal shear stress amplitude increases for very large events — establishing that b-value responds to stress perturbations even when rate does not.

**Test:** Divide the 72-year catalog into solar-cycle phase bins using phase-normalized binning. Reference Adhoc Case A1 for bin count selection: signal is strongest at k=24 and k=32. Compute the Gutenberg-Richter b-value independently for each bin using maximum likelihood estimation. Test whether b-value varies significantly across the annual solar cycle. Compute confidence intervals via bootstrap resampling.

**ISC-GEM advantage:** A0 established that ISC-GEM has 82.9% 2-decimal magnitude precision vs ComCat's 77.5% 1-decimal. This makes b-value MLE substantially more reliable on ISC-GEM.

**Expected outcome if literature holds:** b-value peaks near solstice positions (clamping phase), troughs near equinox positions — inverse of the rate signal.

**Diagnostic value:** If rate and b-value vary inversely in phase, the equinox excess is partly a magnitude-distribution effect (more events capable of reaching M 6.0 threshold during unclamping), not only a nucleation rate effect. This would be a new finding at global M >= 6 scale.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210; magnitude bands derived at analysis time from `usgs_mag`)

**Key references:** Colledge et al. (2025), Bettinelli et al. (2008), Ide et al. (2016)

---

#### Case A3: Magnitude Stratification of the Solar Signal

> **Verdict: KEEP — elevated priority**

**Motivation (strengthened):** The tidal literature review substantially reinforces this case. Métivier (2009) found tidal triggering shows inverse magnitude dependence (stronger for smaller events). Ide et al. (2016) found M < 8.2 indistinguishable from random for tidal phase. Hao et al. (2018) found diurnal periodicity stronger at larger magnitudes in Japan. These findings in opposite directions at different timescales make the magnitude-dependence of the annual solar signal a key discriminator: if the Case 3A signal concentrates in a specific magnitude band, that directly constrains the mechanism.

**Preliminary directional evidence (Adhoc A1b):** The elevated-bin event population (n=1,438) shows no magnitude skew relative to the full catalog. This provides early evidence of magnitude-independence, which A3 will formalize with per-band chi-square tests.

**Test:** Stratify the catalog into magnitude bands (M 6.0–6.4, M 6.5–6.9, M 7.0–7.4, M >= 7.5). Compute chi-square, Rayleigh statistic, and Cramér's V for each band independently using phase-normalized binning. Test whether effect size increases, decreases, or is flat with magnitude.

**Expected outcome if hydrological loading is mechanism:** Effect should weaken or disappear for M >= 7.5.

**Expected outcome if solar-geometric forcing is mechanism:** Effect is flat or increases with magnitude.

**Expected outcome if consistent with tidal literature pattern:** Effect should be weakest at M >= 6 (the opposite of what Case 3A found globally) — making any positive result here a novel finding.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210; magnitude bands M 6.0–6.4, 6.5–6.9, 7.0–7.4, ≥7.5 derived at analysis time from `usgs_mag`)

**Key references:** Hao et al. (2018), Johnson et al. (2017), Métivier et al. (2009), Ide et al. (2016)

---

#### Case A4: Declustering Sensitivity Analysis

> **Verdict: REFRAME — must execute first** (updated with Adhoc A1b declustering window)

**Motivation (strengthened):** Bradley & Hubbard (2024) specifically demonstrated that aftershock temporal clustering produces artificial tidal phase correlations, and that this artifact is large enough to explain several published positive results. The same mechanism applies directly to any periodic signal test. Park et al. (2021) showed that the standard Schuster test requires prior declustering and that failing to remove aftershocks produces false positives. This is the single most important methodological validation for Case 3A.

**Reframe reason (Adhoc findings):** Topic-adhoc Cases A1b and L4/L5 substantially change the declustering landscape:

- **A1b** found that G-K's temporal window (295 days at M6.0) is ~3× longer than the observed clustering footprint of the elevated-bin events (95.6 days p90 IEI). The proposed data-informed reference window is 83.2 km / 95.6 days.
- **L4/L5** established the G-K suppression baseline: 60.2% chi-square degradation on ComCat, with classification confidence of only 0.534 at global M≥6.0 (46.6% of mainshocks borderline).
- **A0b** confirmed that Mechanism 2 (M 6.0 rounding artifacts as spurious aftershock candidates) is the primary data quality concern for ComCat-based declustering. Since this topic uses ISC-GEM, this concern is mitigated.

A two-algorithm comparison (G-K vs Reasenberg) is insufficient given the evidence that G-K over-suppresses at this magnitude level. A three-way comparison is needed.

**Revised test:** Apply three declustering methods to the ISC-GEM catalog:
1. **Gardner-Knopoff (1974)** standard window method — established reference
2. **Reasenberg (1985)** spatial-temporal linking — independent method
3. **A1b-informed custom window** (83.2 km spatial / 95.6 days temporal) — data-informed from elevated-bin clustering footprint

**Note:** The values for the G-K window are derived from widely available formula but have not been directly verified against Gardner & Knopoff (1974) BSSA 64(5), pp. 1363–1367.

For each declustering method, produce a mainshock catalog (retained events) and an aftershock catalog (removed events). Then perform three sub-analyses:

**Sub-analysis A: Scalar signal survival (mainshock catalogs)**

Repeat the chi-square and Rayleigh analysis on each declustered mainshock catalog using phase-normalized binning. Compare p-values, bin distributions, and Cramér's V against the undeclustered result. The comparison between methods (1) and (3) directly tests whether G-K over-suppresses the solar signal by using windows that extend far beyond the observed clustering footprint.

**Sub-analysis B: Post-declustering interval structure (mainshock catalogs)**

Re-run the Adhoc A1b elevated-bin identification on each declustered mainshock catalog: identify the top-3 bins at k=16, 24, 32; compute phase coherence; extract the combined elevated phase intervals. Compare the resulting interval structure to the undeclustered A1b baseline (three intervals at phases ~0.19–0.25, ~0.63–0.66, ~0.88–0.92). This answers the question A4's chi-square alone cannot: do the specific phase intervals persist, collapse, or restructure after declustering?

Possible outcomes and their interpretations:
- All three intervals survive all three methods → genuine three-interval structure, not sequence-driven
- Intervals 2 and 3 disappear under G-K but survive under the A1b window → G-K over-suppresses; those intervals contain a mix of sequences and genuine signal
- Only interval 1 (March equinox) survives across all methods → intervals 2 and 3 were sequence contamination; the underlying signal is a single equinox peak
- Interval structure changes but remains significant → declustering reveals a different underlying pattern that was masked by sequence overlay
- No intervals survive → the entire phase structure was aftershock-driven

**Sub-analysis C: Aftershock population phase preference (aftershock catalogs)**

Compute the chi-square on the aftershock population (removed events) from each declustering method, using the same phase-normalized binning. Test whether the removed events show their own phase preference.

Possible outcomes:
- Aftershock population shows significant phase preference at the same intervals → aftershock sequences are themselves phase-preferring, which is a novel finding (aftershock occurrence timing is modulated by solar phase). This would mean the "aftershock contamination" explanation is incomplete: even if aftershocks inflate the signal, they are not randomly distributed in solar phase.
- Aftershock population shows no phase preference → the signal is entirely in the mainshock population, consistent with the standard interpretation that genuine independent events prefer certain solar phases.
- Aftershock population shows phase preference at *different* intervals than the mainshock population → the full-catalog three-interval pattern is a mixture of two distinct phase preferences (mainshock and aftershock), which merge in the undeclustered catalog.

**Expected outcome:** If the signal is genuine, it persists with comparable or improved p-values under the A1b-informed window while potentially weakening under G-K. If the signal is entirely aftershock-driven, all three methods should reduce it proportionally to the number of events removed. Sub-analysis B will reveal whether the three-interval structure is robust or partially sequence-driven — the most likely outcome given A1b's mixed evidence (strong temporal clustering but no magnitude skew) is that interval 1 persists while intervals 2 and/or 3 are partially reduced.

**Diagnostic value:** A signal that disappears after declustering would indicate the structure is driven by temporal clustering of aftershock sequences. A signal that weakens under G-K but persists under the A1b window would indicate G-K is over-aggressive — removing genuine signal along with aftershocks. Sub-analysis C adds a unique discriminator: if aftershocks themselves are phase-preferring, the traditional "declustering removes contamination" framing is too simple, and the mechanism acts on both mainshocks and aftershocks. The L4/L5 ISC-GEM values will differ from the ComCat-based reference (60.2% suppression, χ²=45.61→18.13) and should not be expected to match.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210; undeclustered baseline)
- G-K mainshocks: `data/iscgem/declustering-algorithm/mainshocks_G-K_global.csv` (n=5,883)
- G-K aftershocks: `data/iscgem/declustering-algorithm/aftershocks_G-K_global.csv` (n=3,327)
- Reasenberg mainshocks: `data/iscgem/declustering-algorithm/mainshocks_reas_global.csv` (n=8,265)
- Reasenberg aftershocks: `data/iscgem/declustering-algorithm/aftershocks_reas_global.csv` (n=945)
- A1b mainshocks: `data/iscgem/declustering-algorithm/mainshocks_a1b_global.csv` (n=7,137)
- A1b aftershocks: `data/iscgem/declustering-algorithm/aftershocks_a1b_global.csv` (n=2,073)

**Key references:** Park et al. (2021), Bradley & Hubbard (2024), Gardner & Knopoff (1974), Ader & Avouac (2013), Adhoc Cases A1 and A1b

---

### Group B: Novel Investigations

---

#### Case B1: Hemisphere Stratification — Phase Symmetry Test

> **Verdict: REFRAME — highest priority novel case** (updated for three-interval pattern)

**Motivation:** The literature attributes seasonal seismicity to hemisphere-specific hydrological loading. The original framing tested whether NH and SH peaks are in-phase (geometric) or anti-phased by 6 months (hydrological) — assuming a bimodal equinox pattern.

**Reframe reason (Adhoc A1b):** The signal resolves into three phase intervals, not two. Only interval 1 (~March equinox) cleanly aligns with an equinox. Interval 2 (~mid-August) leads the September equinox by ~1 month. Interval 3 (~late November) is near the December solstice. This three-interval structure requires a richer hemisphere test than the original bimodal in-phase/anti-phase question.

**Revised test:** Split the catalog into Northern Hemisphere (latitude > 0°) and Southern Hemisphere (latitude < 0°) subsets. Compute phase-normalized solar-phase bin distributions independently for each. Test the following:
1. Do all three elevated phase intervals appear in both hemispheres? If yes, the signal is globally symmetric regardless of hemisphere — consistent with geometric forcing.
2. Is the interval 1 (March equinox) peak present in both hemispheres? This is the cleanest equinox-aligned interval and the strongest test of hemisphere symmetry.
3. Are intervals 2 and 3 hemisphere-specific? If they appear in only one hemisphere (e.g., interval 2 in NH summer, interval 3 in NH autumn), they may reflect hemisphere-specific hydrological or seasonal effects rather than global geometry.
4. Are any intervals 0.5 cycles apart between hemispheres? This would be the strongest evidence for hydrological loading.

**Novel contribution:** No published study has applied this specific phase-symmetry test to a global M >= 6.0 catalog, and the three-interval structure provides finer discriminating power than a simple bimodal test.

**Prediction (geometric hypothesis):** All three intervals appear in both hemispheres at the same phase positions.
**Prediction (hydrological hypothesis):** NH and SH peaks are offset by ~half a year; intervals 2 and 3 may appear only in one hemisphere.
**Prediction (mixed — sequence contamination + genuine signal):** Interval 1 appears in both hemispheres (genuine); intervals 2 and 3 are hemisphere-specific or disappear after declustering.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210; hemisphere split derived at analysis time from `latitude`)

---

#### Case B2: Ocean vs. Continent Location — Hydrological Loading Discrimination

> **Verdict: KEEP**

**Motivation:** Hydrological and snow loading mechanisms are land-surface processes. If the equinox signal persists in purely oceanic seismicity, hydrological loading is not a sufficient explanation. The tidal literature adds relevant context: Scholz et al. (2019) found strong tidal modulation at Axial Volcano (mid-ocean ridge), demonstrating that oceanic/volcanic settings can show periodic signals via different mechanisms than continental hydrology. A positive oceanic result would therefore require careful mechanism discrimination.

**Resource note (Adhoc A1b):** The PB2002 plate boundary file (`lib/pb2002_boundaries.dig`) is now available in the project and can supplement the ocean/continent classification. The A1b boundary proximity analysis (elevated events vs full catalog) already provides partial context for the geographic distribution of the signal.

**Test:** Classify events using a distance-to-coastline threshold (e.g., > 200 km offshore = oceanic, < 50 km or on-continent = continental). Compute the solar-phase bin distribution for each subset using phase-normalized binning. Test whether the signal appears in both, only one, or neither.

**Novel contribution:** No published study has used ocean/continent separation as a mechanism discriminator for annual solar seismicity at global scale.

**Prediction (geometric hypothesis):** Signal appears in both oceanic and continental subsets.
**Prediction (hydrological hypothesis):** Signal is confined to or significantly stronger in continental earthquakes.

**Data Visualizations:** Outputs should include dot-plot or heatmap-like visualizations at the global scale, as well as tighter window views of the  seismically active areas coastlines (i.e. Phillipines, Japan, Chile, Java). These detail images should be done for each mapping treatments in order to review the effects of the mapping method. Color scheme: `oceanic`=blue, `transitional`=green, `continental`=red.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210)
- GSHHG classification (primary): `data/iscgem/plate-location/ocean_class_gshhg_global.csv` (continental=3,799 / transitional=3,459 / oceanic=1,952)
- Natural Earth classification (secondary): `data/iscgem/plate-location/ocean_class_ne_global.csv` (same distribution as GSHHG)
- PB2002 classification (coarse proxy): `data/iscgem/plate-location/ocean_class_pb2002_global.csv` (continental=2,851 / transitional=3,677 / oceanic=2,682)

---

#### Case B3: Tectonic Regime Stratification

> **Verdict: KEEP — reframed with additional motivation**

**Original motivation:** Different fault geometries respond differently to Coulomb stress changes from surface loading. Thrust faults are clamped; normal faults are unclamped; patterns should be anti-phased if loading drives the signal.

**Reframe addition:** Métivier et al. (2009) found slightly greater tidal triggering in normal and strike-slip regions than thrust regions at the global scale — establishing a fault-mechanism baseline for periodic forcing sensitivity. This provides a reference pattern: if the annual solar signal matches the Métivier tidal pattern (normal > strike-slip > thrust), it suggests a similar stress-change geometry even at annual timescales. If it differs, a distinct mechanism is implied.

**Data dependency:** Requires focal mechanism classification from the GCMT catalog, which is not currently in the project data directory. PB2002 plate boundary type (subduction, ridge, transform) could serve as a rough proxy for tectonic regime, but is not a substitute for actual focal mechanism data.

**Test:** Classify events by focal mechanism type (thrust, normal, strike-slip) using the GCMT catalog. Compute solar-phase distributions for each class separately using phase-normalized binning. Test for differences in phase, amplitude, and significance.

**Prediction (loading hypothesis):** Thrust faults show solstice deficits; normal faults show solstice excess; patterns approximately anti-phased.
**Prediction (geometric hypothesis):** All three classes show equinox excess with similar phase.
**Prediction (consistent with Métivier tidal pattern):** Normal and strike-slip show stronger signal than thrust.

**Data sources:**
- Focal mechanism join: `data/iscgem/focal-mechanism/focal_join_global.csv` (n=9,210; 4,874 proximity matches / 4,336 null; mechanism classified from `rake` column)
- PB2002 fallback (if needed): `data/iscgem/plate-location/ocean_class_pb2002_global.csv` — coarse proxy only; not required while GCMT data is available

**Key references:** Métivier et al. (2009), Johnson et al. (2017b)

---

#### Case B4: Depth Stratification — Surface Loading Penetration Test

> **Verdict: KEEP — now has stronger positive motivation**

**Original motivation:** Hydrological loading primarily affects the shallow crust (< 20 km). Intermediate-depth events should be shielded. If the signal persists at depth, surface loading is insufficient.

**Reframe reason:** Zhan & Shearer (2015) found that 70% of deep-focus M > 7 earthquakes (depth > 500 km) occurred in April–September — a suggestive seasonal signal in deep events with *no known surface mechanism*. This directly inverts the original prediction: the literature now hints that large deep events may show *more* seasonality, not less. This makes depth stratification a much richer test than simply checking whether shallow signal exceeds deep signal.

**Revised test:** Stratify the catalog by focal depth: shallow (0–20 km), mid-crustal (20–70 km), intermediate (70–300 km), and deep (> 300 km if sample permits). Compute solar-phase distributions for each band using phase-normalized binning. Test whether signal amplitude increases, decreases, or is non-monotonic with depth.

**Revised predictions:**
- *Surface loading hypothesis*: Signal strongest at 0–20 km; absent at > 70 km.
- *Geometric/deep forcing hypothesis*: Signal present at all depth ranges including > 300 km.
- *Zhan & Shearer pattern*: Deep events show a seasonality, but potentially with different phase than shallow events — a finding of independent interest.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210; depth bands 0–20 km, 20–70 km, 70–300 km, >300 km derived at analysis time from `depth`)

**Key references:** Zhan & Shearer (2015), Dutilleul et al. (2021), Johnson et al. (2017)

---

#### Case B5: Solar Declination Rate-of-Change vs. Position Test

> **Verdict: KEEP — low priority**

**Motivation:** Case 3A encoded seismic events by their position in the solar annual cycle. The bimodal equinox structure could reflect distinct solar variables that both peak at the equinoxes: (a) the rate of change of solar declination is maximum at the equinoxes, and (b) the Earth-Sun distance is near mean value at equinoxes. These variables have different physical implications and have not been distinguished in any published seismicity study.

**Complication (Adhoc A1b):** The three-interval pattern complicates the simple association with any single solar variable. Declination rate peaks at equinoxes — only interval 1 (~March equinox) aligns. Earth-Sun distance has a single annual minimum (January) and maximum (July) that does not match three peaks. Solar declination angle has two zero-crossings at equinoxes, also not matching three peaks. This case is more likely to produce a negative result for all three variables individually, which would itself be informative — suggesting the mechanism is more complex than any single solar parameter, or that one or more of the three intervals reflects sequence contamination rather than a forcing response.

**Test:** Re-encode each event using three separate phase variables: (1) solar declination angle (ranging –23.5° to +23.5°), (2) rate of change of declination (degrees/day), (3) Earth-Sun distance (AU). Compute bin distributions for each using phase-normalized binning. Identify which variable produces the most coherent, highest-significance clustering.

**Novel contribution:** No published seismicity study has explicitly tested solar declination rate-of-change as a potential trigger variable, nor compared it to solar position and Earth-Sun distance within the same catalog.

**Sequencing note:** Should be run after A4 (declustering) to ensure sequence contamination is controlled, and after B1 (hemisphere) to determine whether all three intervals represent genuine signal.

**Data sources:**
- Solar geometry enriched catalog: `data/iscgem/celestial-geometry/solar_geometry_global.csv` (n=9,210; includes `solar_declination`, `declination_rate`, `earth_sun_distance` plus all base catalog columns)

---

#### Case B6: Rolling Window Stationarity Test

> **Verdict: KEEP — elevated priority**

**Original motivation:** A 72-year signal could be driven by a concentrated subset of years (non-stationarity).

**Elevation reason:** Bradley & Hubbard (2024) found that tidal correlations detected in pre-2000 data consistently failed to replicate in post-2000 data — the canonical demonstration of non-stationarity as an artifact signature. Dutilleul et al. (2021) found that Parkfield's semiannual periodicity shifted phase after the 2004 M6.0 earthquake. Both findings support stationarity testing as a critical robustness check. A signal that holds stable phase and amplitude across independent decade-long windows is far more credible than one that concentrates in a single era.

**Data density note (Adhoc A0b):** ISC-GEM has a 1970s spike of 414 unmatched events relative to ComCat, indicating a period-specific difference in catalog construction. The rolling window approach inherently accounts for this by testing each window independently, but the 1970s windows should be flagged if they show anomalous behavior relative to adjacent decades.

**Test:** Compute the Rayleigh statistic and mean phase angle of the solar-cycle distribution using a sliding 10-year window across the 72-year record. Plot the trajectory of p-value, mean vector length, and mean phase. Test whether the phase structure is stable across decades or drifts. If any binned supplementary tests are used alongside Rayleigh, apply phase-normalized binning.

**Data sources:**
- Raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv` (n=9,210)

**Key references:** Bradley & Hubbard (2024), Dutilleul et al. (2021)

---

## Prioritization (Revised — updated with Adhoc findings)

| Case                              | Verdict          | Type        | Mechanism Discriminated         | Priority          |
| --------------------------------- | ---------------- | ----------- | ------------------------------- | ----------------- |
| A4: Declustering sensitivity      | REFRAME          | Replication | Artifact vs. genuine            | **Execute first** |
| B6: Rolling window stationarity   | KEEP ↑           | Novel       | Temporal persistence            | High              |
| A1: Schuster / MFPA spectrum      | REFRAME          | Replication | Periodicity confirmation        | High              |
| B1: Hemisphere phase symmetry     | REFRAME          | Novel       | Hydrological vs. geometric      | High              |
| A3: Magnitude stratification      | KEEP ↑           | Replication | Loading vs. nucleation scale    | High              |
| B2: Ocean vs. continent           | KEEP             | Novel       | Hydrological vs. geometric      | Medium            |
| B4: Depth stratification          | KEEP             | Novel       | Surface vs. deep forcing        | Medium            |
| A2: b-value seasonality           | KEEP             | Replication | Rate vs. magnitude distribution | Medium            |
| B3: Tectonic regime               | KEEP             | Novel       | Stress geometry                 | Medium            |
| B5: Declination rate vs. position | KEEP             | Novel       | Solar variable identification   | Low               |

**Recommended execution order:** A4 → B6 → A1 → B1 → A3 → B2 → B4 → A2 → B3 → B5

**Rationale for changes from original order:**
- A4 remains first: the Bradley & Hubbard replication failures make declustering validation a prerequisite for all downstream cases. Now includes the A1b-informed custom window as a third declustering variant.
- B6 moves to second: stationarity is now a near-prerequisite given the documented failure of tidal correlations across time windows; a non-stationary signal would deprioritize the entire remaining program.
- A3 elevated: the tidal literature's clear magnitude-dependence finding makes this a high-value discriminator for the mechanism question. A1b's no-magnitude-skew finding provides preliminary directional support.
- B1 reframed: the three-interval pattern from A1b requires testing all three intervals per hemisphere, not just a bimodal in-phase/anti-phase comparison.
- B4 verdict simplified to KEEP: no further changes needed from the prior reframe.
- B5 remains low priority: the three-interval pattern complicates the simple mechanism identification this case was designed for. Should follow A4 and B1.
- No cases discarded: the tidal literature and Adhoc findings strengthen the case for every proposed test by either adding motivation, tightening predictions, or providing cross-reference data.

**Cross-cutting dependencies on Adhoc results:**
- All cases: phase-normalized binning (Adhoc A1, `rules/data-handling.md`)
- A4, A1: A1b-informed declustering window (83.2 km / 95.6 days) as third variant
- B1: three-interval phase structure from A1b replaces bimodal equinox assumption
- A3: A1b no-magnitude-skew finding as preliminary directional evidence
- B6: A0b 1970s ISC-GEM data density note
- A2: A0 magnitude precision advantage of ISC-GEM over ComCat for b-value MLE
- B2: PB2002 plate boundary file available from A1b
