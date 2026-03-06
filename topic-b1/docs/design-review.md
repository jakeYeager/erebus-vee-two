# Topic B1: Sequence Topic — Design Review

**Date:** March 5, 2026
**Source conversation:** A3.A2 post-analysis discussion

---

## Origin

This topic was motivated by a finding from A3.A2 (Stratified Schuster/MFPA Periodicity Audit). The MFPA period scan detected several FDR-significant inter-event recurrence periods in the ISC-GEM catalog that are not calendar-year harmonics and have no direct counterpart in the solar-phase chi-square analyses of Topics A2 and A3. These detections raise a question the A-series framework was not designed to answer: do earthquake sequences themselves have internal timing structure, and does that structure interact with the solar signal?

---

## Key empirical inputs from A3.A2

The Schuster/MFPA test operates on `event_time_days` — absolute time since 1950-01-01, monotonically increasing with no modular structure. A detection at period T means events tend to recur at a consistent phase within every T-day window. These are empirical statements about the catalog's internal inter-event timing, independent of the solar-year phase analysis.

**FDR-significant MFPA detections (full catalog, A3.A2 Table 2):**

| Period | Annual fraction | Power | Notes |
|--------|-----------------|-------|-------|
| ~4.7d  | 0.013 | 8.33 | Short-range clustering; likely intra-sequence inter-event grain |
| ~49.9d | 0.137 | 6.71 | Possible sequence activity window |
| ~60.5d | 0.166 | 9.50 | Strongest detection; possible sequence lifetime or inter-sequence rhythm |
| ~243d  | 0.665 | 9.33 | ~8-month recurrence; ~4× harmonic of ~60.5d |
| ~295d  | 0.808 | 8.86 | ~10-month recurrence; no obvious harmonic relationship |
| ~344d  | 0.942 | 7.16 | Near-annual sub-calendar recurrence; warrants investigation |

The ~60.5d and ~243d detections are approximately in a 1:4 harmonic ratio (60.5 × 4 = 242d), suggesting they may reflect the same underlying structure at different timescales. The ~4.7d detection describes the temporal grain at which consecutive events within sequences are coherently spaced.

---

## The core framing shift

Topics A2 and A3 treat the catalog as a **population of events** and ask: does solar phase predict event occurrence? Declustering is used to remove aftershocks as presumed contamination before this question is posed.

However, A2.A4 Sub-C established that aftershock populations carry *stronger* solar-phase signals than mainshock populations — meaning sequences are not noise relative to the solar signal, they are signal-bearing units. Removing them via declustering discards information rather than clarifying the picture.

A Sequence Topic would shift the unit of analysis from **individual events** to **sequences**, asking whether sequence-level properties are modulated by solar phase.

---

## Design principle: empirically grounded cluster windows

Current declustering methods (G-K, Reasenberg, A1b) define sequence membership using magnitude-dependent spatial-temporal formulas or adaptive algorithms. These windows are calibrated to identify aftershock contamination for seismic hazard purposes — not to describe the natural temporal grain of clustering in the catalog.

The A3.A2 detections provide an alternative, data-grounded basis for defining cluster windows:

- **~4.7d** as a short-range intra-sequence temporal threshold — events within ~5 days of each other are sufficiently coherent to produce a detectable Schuster signal. This could serve as a natural cluster-definition window, independent of the G-K magnitude formula.
- **~60.5d and ~243d** as medium-range windows — may describe the active lifetime of sequences at different scales, or a recurrence rhythm in sequence *initiation* times.
- **~344d** as a near-annual recurrence — if confirmed, suggests sequence initiation itself has a near-annual rhythm slightly shorter than the calendar year, which would be a direct connection back to the solar signal at the sequence level.

These empirically derived timescales should be tested as candidate cluster windows in a Sequence Topic rather than defaulting to G-K/Reasenberg formula windows.

---

## Research questions for a Sequence Topic

1. **Sequence initiation and solar phase:** Does the solar phase at the time of a mainshock predict whether it initiates a large sequence? Is there a phase preference for sequence *initiation* comparable to the catalog-wide chi-square signal?

2. **Sequence duration and solar phase:** Do sequences initiated near the equinox (elevated chi-square bins) differ in duration or total event count from those initiated near the solstice (suppressed bins)?

3. **Inter-sequence interval structure:** Does the time between sequence initiations follow the ~344d or ~60.5d recurrence patterns detected in A3.A2? If so, are those inter-sequence intervals themselves modulated by solar phase?

4. **Early vs. late aftershock solar-phase preference:** Within a sequence, do early aftershocks (first half of sequence lifetime by duration) show stronger solar-phase clustering than late aftershocks? This was the motivating question for A3.A1 (abandoned) and remains unanswered. It would distinguish a triggering-pulse mechanism (equinox initiates the sequence; aftershocks decay independently) from a sustained-modulation mechanism (solar forcing remains active throughout the sequence lifetime).

5. **Clustered catalog as primary unit:** Build a catalog where each "event" is a sequence (represented by its mainshock) with sequence-level attributes: initiation solar phase, duration, aftershock count, energy release, inter-sequence interval. Apply the full chi-square and Schuster framework to this sequence-level catalog.

---

## Connection to the solar signal body of work

The A2/A3 oscillation framing (equinox elevation + solstice suppression, established in A3.A3) should be carried into a Sequence Topic as a primary interpretive lens. Specifically:

- If sequences initiated at equinox phases are systematically larger or longer than those initiated at solstice phases, the solar modulation is not just a triggering effect — it shapes the *character* of the seismicity it produces.
- If inter-sequence intervals cluster near the ~344d period, the sequence catalog has a near-annual recurrence rhythm that could be the sequence-level manifestation of the annual solar signal seen in the event-level chi-square.
- The A3.B2 ~5-month NH/SH offset finding should be tested at the sequence level: do NH and SH sequences show different phase preferences for their initiation times, consistent with the hemisphere-asymmetric solar forcing identified in A3.B2?

---

## Relationship to abandoned A3.A1 and A3.A2

A3.A1 (Aftershock Phase-Preference Characterization) was abandoned because A3.A3 and A3.C2 addressed the concentration-vs-diffuse question. However, A3.A1's early/late temporal decay question — does solar-phase clustering within a sequence decay with sequence age? — was explicitly deferred to a future Sequence Topic (see `topic-a3/cases/case-a3-a1-details.md` abandonment rationale). This question maps directly to Research Question 4 above.

A3.A2's ~4.7d, ~60.5d, ~243d, and ~344d period detections provide the empirical timescales for defining what "early" and "late" within a sequence means — replacing the arbitrary `window_secs / 2` midpoint proposed in the original A3.A1 spec with data-grounded thresholds.

---

## Comparison Guidelines

The A2/A3 body of work provides a set of quantitative reference values against which Sequence Topic results can be evaluated. These benchmarks are not pass/fail thresholds — they are interpretive anchors. Results consistent with the benchmarks suggest the sequence data model is capturing the same underlying structure as the event-level analyses. Results that diverge indicate either a genuine data model improvement (signal clarified) or a loss of information (signal degraded or distorted).

### Solar-phase signal strength reference values

These are the primary chi-square and effect-size benchmarks from the event-level analyses at k=24, ISC-GEM full catalog (n=9,210):

| Population | χ² (k=24) | p-value | Cramér's V | Source |
|---|---|---|---|---|
| Full catalog (raw) | 69.37 | 1.52×10⁻⁶ | 0.0181 | A2.A1 / A3.B5 |
| G-K mainshocks only | ~32 | ~0.05 (marginal) | — | A2.A4 Sub-A |
| Reasenberg mainshocks only | ~33 | ~0.03 (marginal) | — | A2.A4 Sub-A |
| G-K aftershocks only | ~66–149 | highly significant | — | A2.A4 Sub-C |
| Mid-crustal band (20–70 km) | 85.48 | 4.02×10⁻⁹ | 0.0285 | A2.B4 / A3.B4 |
| M ≥ 7.5 only | — | 0.027 | 0.0779 | A2.A3 |
| Continental × mid-crustal | — | significant | — | A3.B4 |

**Interpretation key for a sequence-level catalog:**
- If the sequence catalog chi-square falls *below* the G-K mainshock-only level (~32), the data model is losing information relative to naive declustering — a failure mode.
- If it falls *above* the full catalog raw level (69.37), the model is concentrating signal — a strong positive result.
- The aftershock finding from A2.A4 Sub-C (aftershock populations carry stronger signal than mainshocks) sets an upper bound expectation: a sequence catalog that includes aftershock-level information within sequences should approach, but likely not exceed, those aftershock-level chi-square values.

### Cramér's V magnitude trend

A2.A3 found a monotonically increasing Cramér's V by magnitude band (Spearman ρ=1.000). In a sequence catalog, the unit of analysis is a sequence (represented by its mainshock magnitude). The magnitude trend should be preserved or steepened if the sequence framing is adding information. A flattening or reversal of the trend would indicate the sequence model is obscuring a magnitude-driven signal.

| Magnitude band | Cramér's V (event-level, A2.A3) | Expected direction in sequence catalog |
|---|---|---|
| M 6.0–6.4 | 0.0186 | Baseline |
| M 6.5–6.9 | ~0.024 | > M 6.0–6.4 |
| M 7.0–7.4 | ~0.030 | > M 6.5–6.9 |
| M ≥ 7.5 | 0.0779 | Highest |

### Cluster-robust Schuster inflation factor

A3.A2 measured a 6.2× inflation factor between the standard Schuster p-value and the cluster-robust p-value at the 1-day cluster window. This quantifies how much event-level temporal clustering inflates the nominal significance. In a sequence catalog — where each unit of analysis is already a sequence rather than an individual event — the inflation factor should contract toward 1.0× if the data model correctly captures the natural cluster unit. An inflation factor significantly greater than 6.2× would indicate the sequence catalog still contains within-cluster correlations (cluster windows too short); significantly less than 6.2× indicates the clustering is being correctly accounted for.

**Reference:** A3.A2 cluster-robust Schuster inflation factor = 6.2× (1-day window), 329× reported in A2.A1 (methodological discrepancy unresolved).

### MFPA period replication

The A3.A2 FDR-significant period detections are the primary candidate timescales for the Sequence Topic cluster windows. If a sequence catalog built around these windows reproduces the same periods in a subsequent MFPA scan of sequence initiation times or inter-sequence intervals, that is strong self-consistency evidence. The expected replication hierarchy:

| Period | Power (A3.A2) | Replication priority | Interpretation if reproduced |
|---|---|---|---|
| ~60.5d | 9.50 | High | Sequence lifetime or inter-sequence rhythm confirmed at sequence level |
| ~243d | 9.33 | High | 4× harmonic of ~60.5d; if both reproduced, harmonic structure is real |
| ~295d | 8.86 | Medium | No obvious harmonic; warrants independent confirmation |
| ~4.7d | 8.33 | High | Intra-sequence grain; if reproduced in inter-sequence intervals, suggests sequence-to-sequence entrainment |
| ~344d | 7.16 | High | Near-annual; if reproduced at sequence level, direct link to annual solar signal |
| ~49.9d | 6.71 | Low | Weaker detection; lower replication priority |

If the ~344d period is reproduced in the distribution of inter-sequence initiation times, that is the most direct possible connection between the sequence data model and the solar signal observed in A2/A3.

### Stationarity benchmark

A2.B6 and A3.B1 both classified the full event-level catalog as non-stationary:
- Rayleigh significance rate: 38.7% of 10-year rolling windows (threshold: 70% for stationary)
- Chi-square significance rate: 71% of windows (meets threshold, but overridden by high circular SD)
- Circular SD of mean phase angles: 84.4° (threshold: ≤ 40° for stationary)

A sequence catalog, if the sequence-level solar signal is more coherent than the event-level signal, should show improvement on at least one of these metrics — ideally a higher chi-square significance rate with a lower circular SD. If the sequence catalog shows *lower* stationarity than the event-level result, the data model is introducing noise.

### Phase structure reference (A1b intervals)

The three elevated solar-phase intervals from Adhoc A1b provide a phase-domain fingerprint for any new analysis to compare against:

| Interval | Phase range | Approximate calendar period |
|---|---|---|
| Interval 1 | ~0.1875–0.25 | March equinox (~DOY 69–91) |
| Interval 2 | ~0.625–0.656 | Mid-August (~DOY 228–240) |
| Interval 3 | ~0.875–0.917 | Late November (~DOY 319–335) |

A sequence catalog should show mainshock initiation phase preferences consistent with at least Intervals 2 and 3 (most robust across all A2/A3 stratifications). Interval 1 is less stable (absent from SH in A2.B1, inconsistent across declustering methods in A2.A4). Appearance of Interval 1 in a sequence-level analysis would be a notable positive result; absence is expected and not a failure.

### Hemisphere asymmetry benchmark

A3.B2 found that mid-crustal NH and SH signal strengths differ substantially:
- NH mid-crustal: p=2.0×10⁻⁹
- SH mid-crustal: p=0.0017

The NH/SH ratio of significance implies the NH signal is roughly 5–6 orders of magnitude stronger by p-value. A sequence catalog that removes tectonic-composition imbalance (as A3.B2 attempted) should narrow this gap if the asymmetry is compositional. If the NH/SH ratio persists at sequence level after tectonic matching, the asymmetry is real and not a sampling artifact — a meaningful finding for mechanism interpretation.

### Summary: what constitutes a well-functioning sequence data model

A sequence data model is performing well if, compared to the event-level baselines above:

1. Solar-phase chi-square signal is *at least as strong* as the declustered mainshock catalogs and ideally approaches the full raw catalog level.
2. Cramér's V magnitude trend is preserved (monotonically increasing with mainshock magnitude).
3. Cluster-robust Schuster inflation factor is *lower* than the 6.2× event-level value.
4. At least one of the A3.A2 MFPA periods (~60.5d, ~243d, ~344d) is reproduced in sequence initiation timing.
5. Phase preferences for sequence initiation are consistent with A1b Intervals 2 and 3.

Failure on all five criteria simultaneously would indicate the sequence data model is not capturing the solar signal structure seen in Topics A2 and A3 and would require fundamental redesign of cluster window definitions or sequence representation.
