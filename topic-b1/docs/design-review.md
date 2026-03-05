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
