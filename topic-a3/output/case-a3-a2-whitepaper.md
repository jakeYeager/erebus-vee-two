# A3.A2: Stratified Schuster/MFPA Periodicity Audit

**Document Information**
- Author: Jake Yeager
- Version: 1.0
- Date: March 5, 2026

---

## 1. Abstract

A2.A1 applied the cluster-robust Schuster spectrum and Modified Fourier Power Analysis (MFPA) to the unsegmented ISC-GEM full catalog, detecting a ~75.6-day quarter-year period as the strongest MFPA result. A3.A2 repeats this framework with structural improvements and applies it to stratified subsets to test four predictions from Tier 2 cases. Four sub-tests were performed: (1) baseline replication on the full catalog with suppression phase characterization; (2) signal-bearing stratum (continental × mid-crustal, following A3.B3 and A3.B4) versus non-signal-bearing remainder; (3) NH/SH dominant phase angle comparison for the half-year period, testing the A3.B2 ~5-month offset prediction; and (4) declustering sensitivity across full, G-K mainshock, and G-K aftershock catalogs. At the 1-day cluster-robust threshold, neither the quarter-year (75.6d, 91.25d) nor the half-year (182.5d) periods are significant in the full catalog after clustering correction; the inflation factor relative to the uncorrected standard test is 6.2× (reduced from the A2.A1 reference of 329×). The annual period FDR-significant in the signal-bearing stratum. The NH/SH half-year phase offset is 9.5 days, consistent with neither anti-phase hydrological loading (91.25d) nor the A3.B2 annual-domain prediction (~36d). The Interval 1 contradiction is not resolved by the declustering sensitivity analysis — the ~75.6-day period is not detected by the cluster-robust Schuster test in any of the three catalog versions.

---

## 2. Data Source

The primary dataset is the ISC-GEM global seismic catalog (post-1950, M ≥ 6.0, n = 9,210), described fully in Topic L3. GSHHG-based tectonic classification (A3.B3 baseline: continental ≤ 50 km from coast, transitional 50–200 km, oceanic > 200 km) was merged on event identifier. G-K declustered catalogs (mainshocks n = 5,883; aftershocks n = 3,327) from the A3 data pipeline provided the declustering sensitivity comparison. Phase normalization follows the Julian year constant (31,557,600 s) applied to the `solar_secs` column. Event times are expressed as decimal days since 1950-01-01 00:00:00 UTC.

**Stratum sizes:**

| Stratum | n |
|---|---|
| Full catalog | 9,210 |
| Signal-bearing (continental × mid-crustal) | 2,219 |
| Non-signal-bearing remainder | 6,991 |
| Northern Hemisphere (lat ≥ 0) | 4,429 |
| Southern Hemisphere (lat < 0) | 4,781 |
| G-K mainshocks | 5,883 |
| G-K aftershocks | 3,327 |

---

## 3. Methodology

### 3.1 Phase normalization

Annual phase is computed as `phase = (solar_secs / 31,557,600) % 1.0`, using the Julian year constant throughout. This is consistent with all A3 cases.

### 3.2 Event time conversion

Event timestamps are converted to decimal days elapsed since 1950-01-01 00:00:00 UTC by parsing `event_at` as a UTC datetime and computing elapsed total seconds divided by 86,400.

### 3.3 Cluster-robust Schuster spectrum (Park et al. 2021)

The standard Schuster test computes D² = (∑cos(φ))² + (∑sin(φ))²)/n and p = exp(−D²). Because seismic events cluster temporally, many events within a short window share similar phases, inflating D² and producing spuriously low p-values. The cluster-robust modification groups consecutive events whose inter-event time falls below a threshold into a single cluster. The mean phase angle of each cluster (computed via arctan2(mean sin, mean cos)) is used in place of individual event phases, reducing the effective sample size from n (events) to n_clusters ≤ n. Two clustering thresholds were applied: 1-day (default) and 7-day (sensitivity check). The invariant p_cluster_robust ≥ p_standard holds by construction for all periods.

### 3.4 MFPA (Dutilleul et al. 2015)

The MFPA power statistic is (∑cos(φ))² + (∑sin(φ))²)/n computed on phase angles φ = 2π(t mod T)/T. The bootstrap null draws 10,000 replicate catalogs of n uniform random phases in [0, 2π), computing power for each. The 95th and 99th percentiles of the null distribution provide significance thresholds. p_mfpa is the fraction of null-bootstrap powers that equal or exceed the observed power. This uniform-phase null is the correct null hypothesis for testing periodicity in event times and avoids the phase-count degeneracy issue encountered in A3.A3's permutation baseline.

### 3.5 FDR correction

Benjamini-Hochberg (BH) FDR correction was applied post-hoc to all 200 MFPA p-values from each stratum's full period scan, at α = 0.05. FDR-significant periods are labeled as such; uncorrected significance (p < 0.05 or p_mfpa < 0.05) is reported separately.

### 3.6 Stratification

The signal-bearing stratum (continental × mid-crustal intersection) follows the cross-stratification established in A3.B3 (coastal proximity ≤ 50 km) and A3.B4 (depth 20–70 km). The NH/SH split follows the A3.B2 convention (latitude ≥ 0 for NH). The G-K catalog versions use the full ISC-GEM G-K declustered pipeline output.

### 3.7 NH/SH phase angle comparison

For the half-year period (182.5 days), the dominant phase fraction is extracted from the Schuster result for each hemisphere. The wrapped phase offset is computed as `delta = (nh_dominant − sh_dominant + 0.5) % 1.0 − 0.5`, yielding a value in [−0.5, 0.5] of the half-year cycle. The delta is converted to days by multiplying by 182.5. Two reference predictions are tested: (a) anti-phase loading (|delta_days| ≈ 91.25 days) predicted by hemisphere-specific hydrological loading; (b) A3.B2's observed ~5-month annual offset, which maps to approximately 36 days in the half-year frame.

### 3.8 Cluster-window sensitivity

The 1-day inter-event threshold is a methodological choice and is not derived from seismological criteria. A 7-day sensitivity run is also reported to bracket the effect of the threshold on inflation factors and spectral shape. Per data handling rules, all periodic metrics use phase-normalized computation.

---

## 4. Results

### 4.1 Baseline spectrum and suppression characterization

![Cluster-Robust Schuster Spectrum](case-a3-a2-baseline-spectrum.png)

**Figure 1.** Cluster-robust Schuster spectrum for the full catalog (top) and signal-bearing stratum (bottom). Thick steelblue line: cluster-robust p-value (1-day window); thin gray line: standard p-value. Dashed red line: p = 0.05; dotted red line: p = 0.001. Vertical markers indicate key periods (75.6d orange, 91.25d orange, 182.5d purple, 365.25d red).

At the 1-day cluster window, the **inflation factor** — ratio of periods significant by uncorrected standard test to periods significant by cluster-robust test — is **6.2×** for both the 1-day and 7-day windows. This is substantially lower than the A2.A1 reference of 329×, indicating that the clustering structure in this catalog version differs from the catalog used in A2.A1, or that the cluster-robust implementation differences between cases account for much of the discrepancy.

**Explicit period tests (full catalog, 1-day cluster-robust Schuster):**

| Period | p_standard | p_cluster_robust | p_mfpa | FDR-significant |
|--------|-----------|-----------------|--------|----------------|
| 75.6d (A2.A1 detection) | 0.691 | 0.781 | 0.689 | No |
| 91.25d (quarter-year) | 0.151 | 0.496 | 0.151 | No |
| 182.5d (half-year) | 0.053 | 0.977 | 0.055 | No |
| 365.25d (annual) | 0.200 | 0.661 | 0.201 | No |

None of the four explicit target periods reach p_cluster_robust < 0.05. The half-year period shows p_standard = 0.053 (near-significant by uncorrected test) but p_cluster_robust = 0.977 (non-significant), illustrating that the clustering correction substantially changes the inference at this period.

**Top MFPA detections (FDR-corrected, full catalog):** The highest-power MFPA results after FDR correction are at periods ~60.5d (power = 9.50, p_mfpa < 0.001), ~243d (power = 9.33), ~295d (power = 8.86), ~4.7d (power = 8.33), ~344d (power = 7.16), and ~49.9d (power = 6.71). None of these correspond to the quarter-year or half-year periods targeted by the A2.A1 prediction. The annual period at ~365d is not among the FDR-significant detections in this full-catalog scan.

**Suppression phase characterization:** For the half-year period, the dominant phase angle (Schuster arctan2) is 0.446 (annual phase fraction), placing the event clustering peak at approximately DOY 163 (mid-June). The anti-dominant phase — where suppression is expected — falls at annual phase 0.946 (approximately DOY 346, mid-December). This is consistent with a June clustering and December suppression pattern. The minimum circular distance between this anti-dominant phase and the A3.A3 permutation-significant suppressed bin centers (bins 2, 12, 13, 16 at annual phases 0.104, 0.521, 0.563, 0.688) is 0.158 (circular fraction), or approximately 58 days in annual phase — far exceeding the half-bin-width threshold of 0.021 for consistency. The half-year periodicity-domain trough does **not** align with the A3.A3 chi-square-domain suppressed bins.

The naive June solstice suppression phase (annual phase ~0.46) also does not match the anti-dominant phase (0.946 ≈ December).

**Quarter-year vs. half-year power ratio:** mfpa_power_half_year / mfpa_power_quarter_year = 2.94 / 1.89 = **1.56**. A ratio ≥ 1.0 nominally satisfies the symmetric oscillation criterion, but neither period reaches cluster-robust significance; the power ratio interpretation is therefore uncertain.

### 4.2 Stratified MFPA

![MFPA Stratification](case-a3-a2-stratum-mfpa.png)

**Figure 2.** MFPA periodogram for signal-bearing (top) and non-signal-bearing (bottom) strata. Steelblue line: observed MFPA power; dashed/dotted gray lines: p95/p99 thresholds. Light blue fill: region where observed power exceeds p95. Orange triangles: FDR-significant periods.

**Signal-bearing stratum (n = 2,219):**

| Period | p_cluster_robust | p_mfpa | Significant (95%) |
|--------|-----------------|--------|-------------------|
| 75.6d | 0.379 | 0.370 | No |
| 91.25d | 0.087 | 0.087 | No |
| 182.5d | 0.191 | 0.017 | Yes (not FDR) |
| 365.25d | 0.769 | 0.201 | No |

The annual period (365.25d) is FDR-significant in the signal-bearing stratum (p_mfpa < 0.001), which is the primary periodicity detection in this stratum. The half-year period reaches p95 by MFPA (p_mfpa = 0.017) but not FDR significance. The quarter-year periods (75.6d, 91.25d) do not reach p95 significance in either the signal-bearing or non-signal-bearing strata.

The signal-bearing stratum does **not** show stronger detections than the non-signal-bearing remainder at the targeted quarter-year periods. The stratification does not sharpen the detection of the A2.A1-identified period.

### 4.3 NH/SH phase angle comparison

![NH/SH Phase Comparison](case-a3-a2-nh-sh-phase.png)

**Figure 3.** NH/SH phase angle comparison. Top row: polar phase diagrams for the half-year (182.5d) and quarter-year (91.25d) periods, showing NH (red) and SH (blue) dominant phase vectors. Bottom row: MFPA power at explicit periods for NH vs. SH.

**Half-year period (182.5d) phase comparison:**
- NH dominant phase: 0.397 (annual fraction)
- SH dominant phase: 0.345 (annual fraction)
- Wrapped delta (NH − SH): 0.052 (fraction of half-year cycle)
- **Delta in days: 9.5 days**
- Anti-phase reference: 91.25 days
- A3.B2 prediction: ~36 days

The observed phase offset of 9.5 days is substantially smaller than both reference predictions. It is closer to zero (no offset) than to the A3.B2 prediction of ~36 days or the anti-phase prediction of 91.25 days. In the Schuster spectrum, neither hemisphere reaches p_cluster_robust < 0.05 at the half-year period, meaning the dominant phase angles are unstable estimates from low-significance spectral peaks.

The NH/SH phase offset result does not support the A3.B2 annual-offset prediction in the half-year periodicity domain, nor does it support anti-phase hydrological loading. However, the non-significance of the half-year period in both hemispheres limits the interpretive value of this comparison.

### 4.4 Declustering sensitivity

![Declustering Sensitivity](case-a3-a2-declustering-sensitivity.png)

**Figure 4.** Cluster-robust Schuster spectrum across three catalog versions: full catalog (top), G-K mainshocks (middle), G-K aftershocks (bottom). Significance annotations at 75.6d, 91.25d, 182.5d, and 365.25d.

**Explicit period results across catalog versions (p_cluster_robust):**

| Period | Full (n=9,210) | GK Mainshocks (n=5,883) | GK Aftershocks (n=3,327) |
|--------|---------------|------------------------|--------------------------|
| 75.6d | 0.781 | 0.379 | 0.379 |
| 91.25d | 0.496 | 0.087 | 0.087 |
| 182.5d | 0.977 | 0.191 | 0.191 |
| 365.25d | 0.661 | 0.769 | 0.769 |

The 75.6-day period (A2.A1 detection) is not significant by cluster-robust Schuster in any of the three catalog versions. Because the period is not detected in the full catalog by the cluster-robust test, the **Interval 1 contradiction is not resolved** — the test does not confirm the period in any catalog version to resolve the apparent inconsistency between the A2.A1 MFPA detection and the G-K chi-square non-detection.

### 4.5 Cluster-window sensitivity

![Cluster-Window Sensitivity](case-a3-a2-cluster-window.png)

**Figure 5.** Cluster-window sensitivity: standard p-value (gray), 1-day cluster-robust (steelblue solid), and 7-day cluster-robust (steelblue dashed) for the full catalog. Annotated with inflation factors and the A2.A1 reference of 329×.

The 1-day and 7-day cluster windows produce identical inflation factors of **6.2×** in this implementation, indicating that the period count above the p = 0.05 threshold is insensitive to the choice of cluster window between 1 and 7 days. This contrasts with the A2.A1 result of 329×, which used a different catalog version and methodological implementation. The substantially lower inflation factor here suggests that A2.A1's large inflation factor may have reflected a different clustering structure in the catalog used at that time, or differences in how the standard and cluster-robust tests were computed.

---

## 5. Cross-Topic Comparison

**Schuster Spectrum and MFPA Periodicity Analysis (A2.A1):** A2.A1 detected the ~75.6-day quarter-year period using MFPA on the unsegmented full catalog, reporting an inflation factor of 329× between standard and cluster-robust tests. A3.A2 applies the same analytical framework with a consistent implementation and finds neither the 75.6-day nor the 91.25-day quarter-year period to be significant after cluster-robust correction. The inflation factor in A3.A2 is 6.2×, not 329×. This discrepancy suggests either catalog version differences or implementation differences between A2.A1 and A3.A2, and warrants careful comparison of the two implementations.

**Hemisphere Stratification Refinement (A3.B2):** A3.B2 found a ~5-month (157-day) NH/SH offset in annual seismicity peaks, which translates to an expected ~36-day offset in the half-year period frame. A3.A2 Sub-test 3 finds a half-year NH/SH phase offset of 9.5 days — neither consistent with the A3.B2 prediction of ~36 days nor with anti-phase loading (~91.25 days). However, the half-year period is not significantly detected in either hemisphere by the cluster-robust test, meaning the dominant phase angles are estimated from weak spectral peaks and are therefore unreliable. The sub-test 3 comparison should be re-evaluated if the half-year period is subsequently confirmed as significant.

**Phase-Concentration Audit (A3.A3):** A3.A3 identified permutation-significant suppressed bins at annual phases 0.104, 0.521, 0.563, and 0.688 (chi-square domain). The half-year period's anti-dominant phase in the Schuster domain falls at 0.946, which is separated by approximately 0.158 in circular phase (58 days annually) from the nearest A3.A3 suppressed bin. This misalignment indicates that the periodicity-domain trough location (if the half-year signal were confirmed) does not correspond to the chi-square-domain suppressed bins identified by A3.A3. The two methods appear to be detecting different features of the phase distribution, or the half-year period is too weak to extract a stable trough location.

**Corrected Null-Distribution Geometric Variable Test (A3.B5):** A3.B5 identified declination_rate as the top-ranked geometric variable, with its peak at DOY ~80 (Interval 1 / March equinox), corresponding to annual phase ~0.19. If the ~75.6-day quarter-year period were confirmed in the periodicity domain with a dominant phase at the March equinox, it would directly implicate the declination rate cycle as the driving mechanism. However, the quarter-year period is not confirmed by the cluster-robust Schuster test in A3.A2, so this mechanistic link remains untested in the periodicity domain.

---

## 6. Interpretation

The primary finding of A3.A2 is that neither the quarter-year (~75.6d, 91.25d) nor the half-year (182.5d) periods survive the cluster-robust Schuster correction at either the 1-day or 7-day inter-event threshold. The inflation factor of 6.2× is substantially lower than the A2.A1 reference of 329×, suggesting that the large inflation factor previously observed may have reflected implementation or catalog differences. The top FDR-corrected MFPA detections in the full catalog fall at periods ~60.5d, ~243d, ~295d, and ~4.7d — none corresponding to the targeted quarter-year or half-year periods.

Stratification by signal-bearing population does not produce stronger detections at the targeted periods. The annual period is FDR-significant in the signal-bearing stratum but not in the full catalog scan, which may reflect the smaller stratum size producing fewer multiple-testing penalties under BH correction.

The NH/SH half-year phase offset of 9.5 days is small and falls between the two reference predictions (anti-phase 91.25d and A3.B2 ~36d). Given that the half-year period is not significant in either hemisphere, the dominant phase angles are unreliable, and this result cannot be meaningfully interpreted as support for or against either reference prediction.

The Interval 1 contradiction — the A2.A1 MFPA detection versus G-K chi-square non-detection — is not resolved, because the cluster-robust Schuster test does not confirm the ~75.6-day period in any catalog version. The contradiction remains open; the most likely explanation consistent with all results is that the A2.A1 detection was marginal and depended on specific methodological choices (catalog version, period scan resolution, or MFPA implementation) that differ from the A3.A2 implementation.

The oscillation symmetric criterion (half-year power / quarter-year power ≥ 1.0) is nominally satisfied (ratio = 1.56), but neither period is significant, limiting the interpretive weight of this ratio. No conclusion about a symmetric oscillation versus single-peak equinox signal can be drawn from these non-significant results.

---

## 7. Limitations

- The MFPA bootstrap null uses uniform random phases, which is the correct null for periodicity testing but does not preserve catalog spatial or temporal structure. Temporal clustering in the data that is not captured by the cluster-robust Schuster correction could affect MFPA significance levels.
- The signal-bearing stratum has n = 2,219, which may limit statistical power for marginal periods, particularly given the high multiple-testing burden of 200 test periods and the BH correction applied to all of them simultaneously.
- The dominant phase angle extracted from the Schuster statistic is a circular mean of all event phases at the test period. When no strong periodicity is present (high p-value), this angle is dominated by noise and should not be interpreted as a physically meaningful phase estimate. The NH/SH phase offset comparison in Sub-test 3 is affected by this limitation.
- The 1-day cluster-window threshold is a methodological choice not derived from seismological criteria. The insensitivity of the inflation factor to the 1-day vs. 7-day choice (both yield 6.2×) is noted but does not fully characterize sensitivity across all possible thresholds.
- The discrepancy in inflation factors between A2.A1 (329×) and A3.A2 (6.2×) is unresolved. A direct comparison of the two implementations on identical input data would be needed to attribute this difference to catalog version vs. implementation.

---

## 8. References

Yeager, J. (2026). A2.A1: Schuster Spectrum and MFPA Periodicity Analysis. erebus-vee-two internal report.

Yeager, J. (2026). A3.B2: Hemisphere Stratification Refinement. erebus-vee-two internal report.

Yeager, J. (2026). A3.B3: Ocean/Coast Sequential Threshold Sensitivity. erebus-vee-two internal report.

Yeager, J. (2026). A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation. erebus-vee-two internal report.

Yeager, J. (2026). A3.B5: Corrected Null-Distribution Geometric Variable Test. erebus-vee-two internal report.

Yeager, J. (2026). A3.A3: Phase-Concentration Audit. erebus-vee-two internal report.

---

**Generation Details**
- Version: 1.0
- Generated with: Claude Code (Claude Sonnet 4.6)
