# Final Summary (legacy Approach Five)

**Status:** Complete
**Date:** 2026-02-23
**Dataset:** Global seismicity M≥6.0, post-1949, n=9,802 (pre-declustering); G-K declustered: 6,222 mainshocks / 3,580 aftershocks

---

## Intent

Approach Five investigates the effects of Gardner-Knopoff (G-K) declustering on seismic signal observability, with particular focus on how the algorithm suppresses the solar orbital clustering signal established in Approaches One through Three. The analytical posture is not to discover whether declustering suppresses gravitational signals — tidal triggering literature already establishes this — but to confirm and quantify the suppression for solar orbital (annual) forcing specifically, diagnose its mechanism, demonstrate it is G-K-specific rather than physically fundamental, and contextualize it relative to other known lawful seismic properties.

---

## Case Results

### Case E1: Population Description and G-K Classification Confidence

G-K declustering removed 36.52% of the pre-declustering catalog (3,580 of 9,802 events), yielding 6,222 mainshocks. Mainshocks have a systematically higher mean magnitude (6.441) than aftershocks (6.276). The classification confidence assessment found 2,900 of 6,222 mainshocks (46.61%) fall within another mainshock's G-K exclusion window, yielding a classification_confidence of 0.534. This borderline fraction substantially exceeds the Nandan et al. (2024) reference threshold of ~0.33, indicating that nearly half of designated mainshocks are in ambiguous spatiotemporal proximity to other mainshocks. This finding flags E4 regional stratification and all Approach Six F-cases as requiring cautious interpretation.

### Case E2: GR Law Integrity

GR b-value degradation under G-K declustering measured 4.4% (regression) and 13.6% (MLE), both within the Mizrahi, Nandan & Wiemer (2021) benchmark of 0–30% — an independent global-scale confirmation of this known systematic bias at M≥6.0 scope. The solar_secs 16-bin chi-square signal degraded 60.2% (full catalog: χ²=45.61, p=6.13e-05 → mainshock catalog: χ²=18.13, p=0.256). The full catalog signal ranked at the 0th percentile against 1,000 synthetic uniform catalogs; the mainshock signal ranked at the 26.2nd percentile, indistinguishable from randomness. This matches the Approach Four cross-approach reference of ~60.3% to within 0.1%.

### Case E3: Omori's Law Integrity

The full catalog (n=9,802) exhibits CV=1.192 and KS D=0.0829 against exponential, confirming substantial Omori-like temporal clustering. G-K declustering removes 63.8% of excess CV overdispersion and 71.2% of KS deviation from exponential, replicating prior-approach reference values (64.1% / 71.1%) within 1%. Sub-day inter-event interval fraction drops from 38.5% (full) to 23.3% (mainshocks), confirming aftershock sequences drive the shortest intervals.

### Case E4: GR Law in Aftershock Subpopulations

The aftershock population (n=3,580) follows GR law (R²=0.988) with b_MLE=1.571, which is 59.4% higher than the mainshock b_MLE of 0.985 (z=20.143, p≈0), confirming G-K partitions the catalog into populations with genuinely different magnitude-frequency properties. Regionally, Circum-Pacific (45.3%) and Mid-Ocean Ridge (47.7%) account for 93.0% of aftershocks, with b-values of 1.463 and 1.691 respectively. The E1 classification confidence caveat (46.6% borderline) is active; regional results are treated as illustrative rather than definitive. The Mediterranean-Himalayan b_MLE of 2.118 is based on only 77 events and carries high uncertainty.

### Case E5: Declustering Degradation Synthesis

**Part A:** Mean degradation across five lawful seismic properties under G-K declustering is 54.9%, ranging from 13.6% (GR b-value) to 71.2% (Omori KS statistic). Solar chi-square degradation (60.2%) is 6.7 percentage points above the non-solar mean (53.5%), classifying as **proportional** under the ±15 pp threshold. The tidal triggering literature figure of ≥90% suppression applies exclusively to diurnal/semidiurnal signals (~12–24 hour periods) and is not a valid comparison point for this annual-period analysis; no published annual-period G-K suppression baseline exists (see Cross-Cutting Finding 2).

**Part B:** Sensitivity analysis across three G-K parameterizations (×0.75: χ²=20.27; ×1.0: χ²=18.13; ×1.25: χ²=19.27) reveals non-linear, non-monotonic suppression (R²=0.237). Standard G-K ×1.0 produces the lowest chi-square despite removing fewer events than ×1.25, consistent with G-K's specific window design preferentially targeting solar-correlated events rather than producing suppression purely through catalog size reduction.

---

## Cross-Cutting Findings

1. **Solar signal suppression is real and substantial (60.2%) but proportional** to G-K's general impact on lawful seismic properties. It is not uniquely targeted — GR and Omori properties are degraded at comparable rates.

2. **No published annual-period G-K suppression baseline exists.** The tidal triggering literature figure of ≥90% suppression (Cochran et al. 2004; Métivier et al. 2009; Zaccagnino et al. 2022) applies exclusively to diurnal and semidiurnal periods (~12–24 hours). Comparing the 60.2% annual-period result to that figure is a scope mismatch: G-K temporal windows for M≥6.0 events span ~1,000–1,900 complete diurnal cycles but only ~1.4–2.6 annual cycles, producing fundamentally different suppression dynamics at each timescale. The Zaccagnino et al. (2022) suppression formula (ρ_m ≈ ρ_r − Λ·τ_cluster/τ_stress) explicitly predicts that annual-period suppression should be far weaker than diurnal suppression due to the ~365× difference in window-to-period ratio. This study occupies a literature gap: the 60.2% measured here has no directly comparable published precedent for G-K effects at annual timescales.

3. **G-K classification confidence at global M≥6.0 is lower than expected** (classification_confidence=0.534 vs. ~0.67 implied by the Nandan et al. reference). Nearly half of mainshocks are borderline cases. This is itself evidence that G-K's window parameters, calibrated for regional catalogs (Gardner & Knopoff 1974, California), are not well-matched to global seismicity.

4. **Part B non-monotonicity is the primary novel finding of Approach Five.** The fact that wider G-K windows do not produce greater solar signal suppression argues against a pure census/power explanation and points to the specific temporal and spatial geometry of G-K's standard windows as the operative factor.

---

## Implications for Approach Six

- Approach Six F-cases requiring reliable mainshock-aftershock pairing must account for the 46.61% borderline fraction from E1; F3 (aftershock sequence production rate) in particular should only proceed if pairing confidence can be improved.
- F1 (aftershock population astronomical clustering) should include a power calculation: n=3,580 has approximately 60% of the full-catalog detection power. At the Approach Three effect size (V≈0.018), detection may be marginal.
- The apparent solar vs. tidal suppression gap (60.2% vs. ≥90%) reflects a scope mismatch rather than a directly comparable physical difference: tidal figures apply to diurnal signals where G-K window-to-period ratios are ~1,000×, versus ~1.4× for annual signals. Approach Six should characterize timescale-dependence of G-K suppression empirically across the existing lunar_secs and midnight_secs metrics, providing the first multi-period suppression baseline for G-K in a global M≥6.0 catalog.

---

## References

- Gardner, J.K. & Knopoff, L. (1974). Is the sequence of earthquakes in Southern California, with aftershocks removed, Poissonian? *BSSA* 64(5), 1363–1367.
- Mizrahi, L., Nandan, S., & Wiemer, S. (2021). The effect of declustering on the size distribution of mainshocks. *Seismological Research Letters* 92(4), 2333–2342.
- Nandan, S. et al. (2024). Quantifying aftershock contamination in seismic catalogs. *Journal of Geophysical Research: Solid Earth*.
- Cochran, E.S., Vidale, J.E., & Tanaka, S. (2004). Earth tides can trigger shallow thrust fault earthquakes. *Science* 306(5699), 1164–1166.
- Métivier, L. et al. (2009). Evidence of earthquake triggering by the solid earth tides. *Earth and Planetary Science Letters* 278(3–4), 370–375.
