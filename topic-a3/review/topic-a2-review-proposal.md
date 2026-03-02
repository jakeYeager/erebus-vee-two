**Topic A2: Research Review Proposal**

Annual Solar Signal in Global Seismicity

*Prepared: March 2026*

**1. Overview of Topic A2**

Topic A2 is a systematic investigation into whether global M ≥ 6.0
seismicity shows a statistically significant and physically meaningful
signal correlated with position in the annual solar cycle. The work
builds on a preliminary result from Case 3A (ComCat catalog, n=9,802)
and replicates, extends, and stress-tests that result using the ISC-GEM
catalog (n=9,210, 1950--2021) across ten analysis cases.

The Adhoc cases (A0, A0b, A1, A1b) provided the methodological
foundation: catalog characterization and comparison, phase-normalized
binning validation, and data-informed characterization of the
elevated-phase event population. These Adhoc results shaped the design
of all ten Topic A2 cases.

**2. Overall Success of the A2 Investigation**

**2.1 Core Signal: Confirmed but Diminished Under Declustering**

The raw ISC-GEM catalog reproduces the solar-phase chi-square signal
with high confidence (χ²=69.37, p=1.52×10⁻⁶, Cramér\'s V=0.0181 at
k=24). This confirms that the initial Case 3A result is not a
catalog-specific artifact. However, all three declustering methods (G-K,
Reasenberg, A1b-informed) substantially suppress the signal: G-K by
53.75%, A1b-informed by 56.16%, and Reasenberg by 42.43%. After G-K and
A1b declustering, the signal does not reach conventional significance at
any bin count tested. Only the Reasenberg mainshock catalog remains
marginally significant at k=24 (p=0.016).

This is the central interpretive tension of the investigation: the
signal exists in the raw catalog at high confidence, but its status in a
fully declustered catalog is ambiguous. Whether this represents genuine
aftershock contamination of a real signal, or an artifact driven
primarily by aftershock clustering, cannot be resolved from the current
case set alone.

**2.2 Phase Structure: Three Intervals, Not Bimodal**

Adhoc A1b refined the initial bimodal equinox characterization into
three phase intervals: Interval 1 (phase 0.19--0.25, \~March equinox),
Interval 2 (phase 0.63--0.66, \~mid-August), and Interval 3 (phase
0.88--0.92, \~late November). This three-interval structure is the
primary structural finding of the program and is the reference baseline
against which all stratification cases are evaluated.

Under declustering, Intervals 2 and 3 survive across all three methods;
Interval 1 does not survive at k=24 under any declustering approach,
suggesting it may be more influenced by aftershock contamination or may
require a larger sample to resolve.

**2.3 Stratification Results Summary**

Across the B-series stratification cases:

-   Hemisphere (B1): Signal is significant in both Northern
    (p=3.96×10⁻⁷) and Southern (p=2.03×10⁻³) hemispheres. Interval 1 is
    present in NH only; Intervals 2 and 3 appear in both. This partial
    asymmetry weakly disfavors a simple global-geometric mechanism.

-   Tectonic setting (B2): The continental/transitional class is
    significant; the oceanic class approaches but does not reach
    significance (p=0.0607). Consistent with a surface loading or
    hydrological mechanism but not conclusive.

-   Tectonic regime / focal mechanism (B3): No mechanism class reaches
    significance individually. Normal faults show the highest nominal V
    (0.0353) but the sample (n=720) lacks power. No single mechanistic
    prediction is clearly supported.

-   Magnitude (A3): Effect size increases monotonically with magnitude
    (Spearman ρ=1.000), peaking at M≥7.5 (V=0.0779). This argues against
    an explanation that relies purely on aftershock contamination, since
    large mainshocks are less dominated by aftershock bias.

-   Depth (B4): Signal is significant only in the mid-crustal band
    (20--70 km, p=4.02×10⁻⁹). Shallow events (0--20 km), where surface
    loading effects would be strongest, are not significant. This is
    inconsistent with a simple surface-loading explanation.

-   Temporal stationarity (B6): Non-stationary. The Rayleigh mean vector
    length is elevated in post-2000 windows and reduced in earlier
    decades. The most significant windows coincide with the 2004 Sumatra
    M9.1 aftershock sequence period, raising the concern that a single
    large aftershock sequence is driving much of the raw-catalog signal.

-   Periodicity structure (A1): No cluster-robust Schuster significance
    at the annual period. MFPA identifies a robust \~75.6-day period but
    this does not predict all three A1b intervals. The standard Schuster
    annual signal is an artifact of aftershock clustering.

-   Solar geometry decomposition (B5): Geometric variables (declination,
    declination rate, Earth-Sun distance) produce larger Cramér\'s V
    than the phase variable, but this largely reflects the non-linear
    mapping between calendar time and geometry rather than genuine
    physical forcing. No single geometric variable accounts for the
    three-interval structure.

-   b-value seasonality (A2): Tests whether the magnitude-frequency
    distribution shifts with solar phase, complementing the event-rate
    analysis. Results are in the case record.

**3. Limitations Assessment**

**3.1 Documented Limitations by Theme**

**Declustering (A4, A1)**

-   G-K window parameters are unverified against the original 1974
    paper; secondary-literature values for M6.0 (49 km / 295 days)
    differ from values used (54 km / 915 days). This is flagged as
    unresolved and should be verified before any resubmission.

-   The A1b-informed window (83.2 km / 95.6 days) is data-derived, not
    literature-validated. Its conservatism relative to G-K (56% vs 54%
    suppression) may understate true aftershock removal.

-   Reasenberg parameters use published defaults; sensitivity to r_fact
    and tau_min/max is not tested.

-   Addressability: Verifiable. Obtain and check original G-K paper
    parameters and compare against ZMAP implementation. Low effort, high
    priority.

**Catalog Quality (B6, A0b)**

-   ISC-GEM has a documented \~414-event density excess in the 1970s
    (Adhoc A0b). Although the 1970s windows in B6 do not show anomalous
    Rayleigh R, the mechanism behind the excess is unexplained.

-   Early decades (1950s--1960s) have lower completeness, creating
    temporal non-uniformity that could bias periodic tests.

-   Addressability: The density excess requires ISC-GEM documentation
    review or contact with the ISC-GEM team. Completeness variation is
    inherent to the catalog and cannot be resolved without using a
    completeness-trimmed subset.

**Periodicity Testing (A1)**

-   MFPA bootstrap uses 1,000 replicates; resolution is limited to
    \~0.001. Marginal detections near the 95th percentile are not
    reproducibly stable.

-   The cluster threshold for Schuster (1-day) is not
    sensitivity-tested. A 6-hour or 7-day threshold would yield
    different cluster correction factors.

-   No multiple-comparison correction applied across 300 scanned
    periods. With \~15 expected false positives at the 5% threshold,
    individual period-level claims at marginal significance are
    unreliable.

-   Addressability: MFPA bootstrap can be extended to 10,000 replicates
    (moderate compute cost). Cluster threshold sensitivity requires
    re-running with 2--3 threshold values. Bonferroni or FDR correction
    is straightforward to add post-hoc.

**Focal Mechanism Coverage (B3)**

-   GCMT match rate is 52.9%, with systematic pre-1976 exclusion. This
    means the mechanism-stratified sub-sample is not representative of
    the full catalog.

-   The unmatched events (n=4,336) show V=0.0260, higher than any
    mechanism class, suggesting the matched sample may be
    compositionally biased toward settings where the signal is weaker.

-   Addressability: Partial. GCMT coverage before 1976 cannot be
    extended. A sensitivity test using the matched-only raw catalog vs.
    the full catalog chi-square would quantify bias. Additionally, ISC
    focal mechanism data could supplement GCMT for pre-1976 events.

**Declustering and Raw Catalog Use in B-Series Cases (B4, B3, B2)**

-   B4, B3, and B2 all use the raw catalog. Given A4\'s finding that
    aftershock contamination substantially inflates the chi-square, the
    stratification results may reflect the aftershock distribution
    across strata rather than mainshock phase preferences.

-   Addressability: All B-series cases should be re-run on declustered
    catalogs. The Reasenberg mainshock catalog (which retains marginal
    significance) is the most appropriate declustered version for these
    tests. This is a significant analysis extension.

**Temporal Non-Stationarity (B6)**

-   The rolling-window analysis uses overlapping 10-year windows,
    introducing strong autocorrelation. The Bonferroni correction treats
    windows as independent, overstating multiple-comparison control.

-   The post-2000 signal concentration coincides with the 2004 Sumatra
    aftershock period but is not definitively linked to it.

-   Addressability: Non-overlapping decade windows would provide
    independent comparisons, though at the cost of resolution. A
    targeted analysis removing the Sumatra sequence (and other M≥8.5
    aftershock trains) from the raw catalog would test whether a single
    sequence drives the stationarity pattern.

**4. Cross-Topic Synthesis Opportunities**

**4.1 Key Inter-Case Relationships**

Several consistent findings emerge across multiple cases that warrant
synthesis:

-   Intervals 2 and 3 are the robust intervals. They survive
    declustering (A4), appear in both hemispheres (B1), appear in
    continental and transitional settings (B2), and are structurally
    consistent across cases. Interval 1 (March equinox) is the weakest
    and least stable interval under scrutiny.

-   Mid-crustal seismicity (20--70 km) drives the signal (B4). This
    depth band overlaps strongly with subduction zone seismicity, which
    is predominantly thrust in mechanism (B3) and predominantly
    transitional/near-continental in setting (B2). The B2, B3, and B4
    results point to the same population of events.

-   Large magnitude events show the strongest effect (A3, V=0.0779 at
    M≥7.5). B6 implicates post-2000 windows. These are consistent if the
    signal is driven by a small number of large-magnitude event clusters
    in the 2000--2014 window.

-   The non-annual periodicity (\~75.6 days) from A1 MFPA does not
    predict all three A1b intervals. This is a genuine unresolved
    tension: the chi-square detects a non-uniform annual distribution,
    but the periodicity test does not find an annual period.

**4.2 Potential Synthesis Cases**

Based on inter-case relationships, the following synthesis analyses are
proposed:

1.  **Aftershock-sequence contribution analysis (Synthesis S1):**
    Identify the top 5--10 largest post-2000 mainshock sequences (M≥8.0)
    and assess their contribution to the raw catalog chi-square
    statistic. Remove each sequence and retest. This directly tests the
    B6 hypothesis that post-2000 stationarity is driven by specific
    large sequences.

2.  **Mid-crustal / near-continental combined stratification (Synthesis
    S2):** Cross-stratify by depth band (20--70 km) AND tectonic setting
    (continental/transitional) to isolate the event sub-population that
    consistently drives the signal across B2 and B4. Test this joint
    sub-population on both raw and declustered catalogs.

3.  **Large-magnitude declustered phase analysis (Synthesis S3):** Focus
    on M≧7.0 events in the declustered catalog. The A3 magnitude trend
    and B6 large-event contribution suggest that large mainshocks may
    retain the phase signal more cleanly than the full catalog. With
    n≈700--900 events at M≧7.0, statistical power will be reduced but
    the result would be free of aftershock contamination.

4.  **Pre-/post-2000 stratified replication (Synthesis S4):** Split the
    catalog at year 2000 and run all primary tests (chi-square, interval
    recovery) on each sub-period independently. If the signal is
    concentrated in post-2000 data, the pre-2000 catalog should show
    near-null results. This is the cleanest test of temporal
    non-stationarity from B6.

**5. Interpretation Summary and Forward Statement**

**5.1 What the Body of Work Establishes**

The following statements are supportable from the completed case set:

-   The ISC-GEM global M≥6.0 catalog shows a statistically significant
    non-uniform distribution of events across the solar year in the raw
    catalog. This result is robust to bin count (k=16, 24, 32) and is
    not present in the lunar-phase or local-time distributions (p\>0.3
    for both).

-   The phase structure consists of three elevated intervals, two of
    which (August and November) are stable across declustering,
    hemispheres, and tectonic settings. The March equinox interval
    (Interval 1) is less stable.

-   The raw-catalog signal is substantially suppressed by all three
    declustering methods. Only the Reasenberg mainshock catalog retains
    marginal significance, and no catalog-wide Schuster annual signal
    survives cluster correction. The signal is therefore partially
    attributable to aftershock clustering.

-   The signal increases monotonically with magnitude band, is
    concentrated in mid-crustal depths (20--70 km), is stronger in
    near-continental settings, and is temporally concentrated in the
    post-2000 era. These findings point toward a small number of large,
    mid-crustal, near-continental earthquake sequences as the primary
    contributors to the signal.

-   No single physical mechanism is clearly identified. The depth result
    (shallow events not significant) disfavors a simple surface-loading
    explanation. The hemisphere asymmetry of Interval 1 is consistent
    with a hemisphere-specific mechanism for that interval. The
    post-2000 concentration raises the possibility that the signal is
    driven by a limited number of extraordinary seismic sequences rather
    than a persistent annual forcing.

**5.2 What Remains Unresolved**

-   Whether the Reasenberg-significant mainshock signal represents a
    genuine annual forcing or a residual aftershock artifact not removed
    by the Reasenberg algorithm.

-   Whether the three-interval phase structure is produced by a single
    underlying period (which MFPA does not identify at the annual scale)
    or by a multi-component, non-periodic process.

-   Whether the post-2000 signal concentration reflects a physical
    change in the annual seismicity pattern, a change in catalog
    completeness, or the coincidence of major aftershock sequences with
    particular solar-phase bins.

-   The mechanism: hydrological loading, geometric/tidal forcing, or a
    confound driven by aftershock clustering within specific sequences.

**5.3 Proposed Forward Statement for Topic A2**

> *Global M≥6.0 seismicity in the ISC-GEM catalog (1950--2021) shows a
> non-uniform distribution of event timing across the solar year,
> characterized by three elevated phase intervals near the March
> equinox, mid-August, and late November. This pattern is statistically
> robust in the raw catalog (χ²=69.37, p\<10⁻⁵) and survives hemisphere,
> ocean/continent, and magnitude stratification in its core structure
> (Intervals 2 and 3), but is substantially suppressed by declustering.
> The signal is concentrated in mid-crustal, near-continental seismicity
> and in the post-2000 observation window, consistent with a
> contribution from major subduction-zone aftershock sequences. No
> single physical mechanism is identified. The investigation
> demonstrates that the initial Case 3A solar-phase signal is replicable
> and partially generalizable, but that its physical origin requires a
> targeted analysis of large-magnitude mainshock sequences in
> declustered catalogs before a mechanism-level conclusion can be
> drawn.*

**5.4 New Topic or Extended Literature Review?**

Based on the current evidence, the body of work does not yet warrant a
fully independent new topic. The central ambiguity---whether the signal
survives in mainshock catalogs at the required significance level---is
not resolved. The proposed synthesis cases (S1--S4) should be completed
within the existing Topic A2 framework before a new topic is designated.

A new topic with a fresh literature review would be warranted if
synthesis cases establish either: (a) a significant annual signal in
large-magnitude declustered catalogs that clearly exceeds aftershock
contamination levels, or (b) a mechanistic identification that points to
a specific physical process absent from the current literature review
(e.g., far-field triggering by specific seasonal geophysical loads in
subduction zones). If neither condition is met, the appropriate outcome
is a refined null finding and a characterization of the aftershock
contribution.

**6. Additional Topics and Gaps Not Covered**

**6.1 Topics That Should Be Added to the Review**

5.  **Cross-case power analysis:** No case computes the minimum
    detectable effect size (MDES) for stratified sub-samples. Several
    negative results (B3 normal faults, B4 shallow events) may reflect
    insufficient statistical power rather than absence of signal. A
    post-hoc power calculation for each stratification sub-sample would
    clarify whether non-significant results are informative or simply
    underpowered.

6.  **Aftershock phase preference audit:** Case A4 noted that
    aftershocks preferentially align with the A1b solar intervals---this
    is the key mechanistic question. A dedicated case examining whether
    aftershock events are clustered in specific phase bins (and which
    primary sequences contribute most) would directly address whether
    aftershock contamination explains or extends the signal.

7.  **Spatially-resolved signal mapping:** All cases operate on global
    or broad-class aggregates. A grid-based or region-based chi-square
    map (e.g., 10°×10° cells) would identify which geographic regions
    drive the global signal. This could resolve whether the effect is
    concentrated in specific subduction zones (e.g., Sumatra-Andaman,
    Tonga, Japan) or is distributed globally.

8.  **Comparison with tidal triggering literature:** The literature
    review focuses on annual/seasonal forcing. The tidal triggering
    literature (Hao et al. 2018, Yan et al. 2022) is cited but not
    systematically reviewed in the context of sub-annual harmonic
    structure. A targeted comparison of the MFPA \~75.6-day period
    against tidal-cycle harmonics would determine whether the detected
    sub-annual periodicity overlaps with known tidal periods.

9.  **b-value seasonal analysis completeness check (A2 integration):**
    Case A2 (b-value seasonality) examines the magnitude-frequency
    distribution by solar phase. Its results should be explicitly
    integrated with A3 (magnitude stratification). If b-value shifts
    correlate with the elevated-phase intervals, this would indicate
    that the elevated-interval events are not simply more numerous but
    also statistically larger, which would be a significant mechanistic
    constraint.

**6.2 Clarifying Questions**

The following questions should be resolved before or during the review
process:

-   Declustering standard: Which declustering result is treated as the
    reference catalog for the B-series stratification cases going
    forward---raw, G-K, Reasenberg, or A1b-informed? The current cases
    use raw, which creates a known aftershock confound. The review
    should establish a consistent reference.

-   Significance threshold: Is the review using α=0.05 per individual
    test or a family-wise error rate across the full case set? With 10
    primary cases and multiple sub-analyses per case, uncorrected
    p\<0.05 thresholds likely produce false positives.

-   Treatment of Interval 1: Given that the March equinox interval does
    not survive declustering and is absent from the Southern Hemisphere,
    should it be retained as a primary result or demoted to a
    provisional finding requiring further validation?

-   B6 post-2000 concentration: Has the specific contribution of the
    2004 Sumatra and 2011 Tōhoku sequences to the raw-catalog chi-square
    been quantified? If not, this is the single most important
    outstanding question for interpreting the stationarity result.

-   G-K parameter verification: Has the G-K window table been verified
    against the original 1974 paper or ZMAP implementation? This is
    flagged as unresolved in the data requirements document and affects
    the validity of A4 and all A4-dependent cases.

-   Catalog scope: Is there a plan to test the findings on an
    independent catalog (e.g., USGS PDE, ISC reviewed catalog, or a
    regional catalog) to confirm that the ISC-GEM-specific completeness
    patterns are not driving the results?

-   Publication intent: Is the goal of this review to prepare the body
    of work for publication, to inform a new investigation, or to
    determine whether to close the A2 program? The answer affects which
    gaps are prioritized.

**7. Proposed Review Process**

The following sequence is recommended for conducting the Topic A2
review:

10. **Phase 1 --- Clarification (Pre-Review):** Resolve the clarifying
    questions in Section 6.2 before the review begins. Verify G-K
    parameters. Establish the reference declustering catalog. Define the
    significance framework.

11. **Phase 2 --- Case-by-Case Audit:** Review each case in execution
    order (A4 → B6 → A1 → B1 → A3 → B2 → B4 → A2 → B3 → B5). For each
    case: (a) confirm the Limitations block items are documented; (b)
    check which limitations are addressable vs. inherent; (c) note any
    unresolved cross-case dependencies.

12. **Phase 3 --- Cross-Case Synthesis:** Compile the Cross-Topic
    Comparison blocks across all cases. Assess whether the inter-case
    relationships support or undermine the proposed forward statement in
    Section 5.3. Identify which synthesis cases (S1--S4) are necessary
    to resolve the central ambiguity.

13. **Phase 4 --- Synthesis Case Execution:** Execute synthesis cases S1
    (aftershock-sequence contribution) and S4 (pre-/post-2000 split) as
    the highest priority, followed by S3 (large-magnitude declustered
    analysis) and S2 (mid-crustal/near-continental joint
    stratification).

14. **Phase 5 --- Forward Decision:** Based on synthesis case results,
    decide between: (a) close A2 with a refined null finding and
    write-up; (b) extend A2 with addressable limitation fixes (G-K
    verification, declustered B-series re-runs, power analysis); or (c)
    open a new topic with a focused literature review if the signal is
    confirmed in large-magnitude declustered catalogs.

*This proposal was generated as a structured review of the Topic A2
whitepaper case set. All statistics cited are drawn directly from the
individual case whitepapers.*
