# Topic Adhoc - Case A1b Planning

## Intent

Drill down into the elevated bins identified in Case A1 to test whether the events driving the `solar_secs` chi-square signal show temporal and/or spatial clustering. If they do, that clustering behavior can inform the construction of a more calibrated declustering window than the generic G-K parameters — specifically addressing the concern that the standard G-K windows are too permissive at global M≥6.0 scale.

**Data Source:** ISC-GEM catalog: `data/global-sets/iscgem_global_events.csv`

**Prerequisite:** Case A1 results and analysis script (`topic-adhoc/output/case-a1-results.json`, `topic-adhoc/src/case-a1-analysis.py`)

---

## Background

Case A1 established that `solar_secs` chi-square significance is robust across bin counts. The signal strengthens with resolution (χ²=42.5 at k=16 → χ²=78.2 at k=32) while Cramér's V remains stable (0.0175–0.0177). This flat effect size across bin counts means the signal is distributing across multiple coherent elevated bins at each resolution rather than concentrating in a single spike — a pattern consistent with a broad phase preference rather than a point artifact.

The question A1b addresses: **are the events driving those elevated bins independent seismic events, or are they temporally and spatially clustered in a way that suggests sequence contamination?** The answer has two possible implications:

- **Clustered** → consistent with aftershock or sequence residuals that G-K failed to remove; the elevated bins may be inflated by sequence events, and a tighter declustering window could reduce them
- **Dispersed** → consistent with independent events genuinely preferring those solar phases; the signal is real and G-K suppression (L4/L5) was removing signal, not noise

---

## Scope

Focus exclusively on `solar_secs`. `lunar_secs` and `midnight_secs` showed no significance in A1 and provide no elevated bins of interest.

Bin counts to analyze: **k=16, 24, 32** (the three Bonferroni-significant results). k=8 is excluded — it failed correction and at that coarseness each bin spans ~46 days, which is too broad to resolve useful spatial-temporal structure.

Elevated bins per count: **top 3 by positive deviation from expected** (observed − expected).

---

## Proposed Analyses

### Analysis 1: Phase Mapping and Cross-Resolution Coherence

Before any event-level work, map the top 3 elevated bins at each bin count to solar phase ranges `[phase_start, phase_end)` and convert to approximate calendar date ranges.

The key question: do the top 3 bins at k=16, 24, 32 point to the same or overlapping phase intervals? If the elevated phase ranges are coherent across all three bin counts, this confirms the signal is a stable annual feature rather than a resolution artifact. If they are incoherent, the "elevated bins" are noise fluctuations aggregated by chi-square.

Output: a phase alignment table showing the elevated phase ranges for each k, and whether they overlap. A single combined "elevated phase set" — the union of phase ranges that are elevated across ≥2 of the 3 bin counts — will be the primary event population for all subsequent analyses.

### Analysis 2: Temporal Distribution

For events in the combined elevated phase set, examine their distribution across the study period (1950–2021).

- Events per year: does the elevated-bin population distribute uniformly across the 72-year span, or is it concentrated in specific years or decades?
- If concentrated in early decades (1950s–1970s): possible catalog artifact (ISC-GEM backfill density) — compare to the full catalog year distribution to distinguish signal from data density
- Inter-event intervals (IEI): compute IEIs for elevated-bin events and compare to the full catalog IEI distribution. Short IEIs (days to weeks) within the elevated-bin population would indicate temporal clustering consistent with aftershock sequence contamination

### Analysis 3: Spatial Distribution

For events in the combined elevated phase set, examine geographic distribution.

- Latitude/longitude scatter: is the elevated-bin population spatially uniform or concentrated in specific regions?
- Mean nearest-neighbor distance: compare elevated-bin events to a random draw of equal size from the full catalog. If nearest-neighbor distance is significantly smaller, the population is spatially clustered
- Latitude banding: compute events per 30° latitude band for elevated-bin vs full catalog. Solar orbital effects would be expected to be globally symmetric; regionally concentrated elevation would suggest a different mechanism

### Analysis 4: Magnitude Distribution

Compare the magnitude distribution of elevated-bin events to the full catalog.

- If elevated bins skew toward M6.0–6.1: consistent with near-threshold rounding artifacts (as established in A0) or aftershock contamination (aftershocks are magnitude-suppressed relative to mainshocks)
- If magnitude-uniform: consistent with the magnitude-independence finding from Approach Three and a genuine global effect

### Analysis 5: Declustering Window Estimation

This is the applied output of the case. Based on Analyses 2–4, characterize the spatial and temporal footprint of the elevated-bin event population and compare it to the standard G-K window parameters for M≥6.0:

| G-K parameter | M=6.0 | M=7.0 | M=8.0 |
| --- | --- | --- | --- |
| Spatial window | ~49 km | ~156 km | ~493 km |
| Temporal window | ~295 days | ~790 days | ~2,117 days |

If the elevated-bin events show clustering tighter than these windows (e.g., within 20 km and 60 days), the G-K windows are not the limiting factor — something else is driving the clustering. If they are dispersed beyond these windows, G-K cannot capture them by design, and a custom algorithm with different spatial or temporal parameters would be needed.

Propose specific window parameter ranges (spatial radius in km, temporal window in days) informed by the observed clustering metrics — not as a final calibration, but as a quantified starting point for further algorithm development.

---

## Interpretive Framework

The analysis should resist two failure modes:

1. **Confirmation bias**: finding clustering and concluding the signal is "just aftershocks" without quantifying whether the clustering is materially different from the full catalog baseline
2. **Hypothesis bias**: finding dispersion and concluding the signal is "real" without checking whether the elevated-bin population differs from the full catalog in any other systematic way (e.g., geographic concentration in a seismically active region)

The declustering implication section should be framed as: *if the observed clustering metrics were used to calibrate a window, the window parameters would be X — this is a quantified reference point, not a recommendation to adopt those parameters without further validation.*

---

## Questions Before Spec

1. **Tectonic zone labels:** Spatial analysis would be substantially more informative with tectonic zone classification (subduction, intraplate, ridge, etc.). Do these exist in any project data file, or would they need to be derived from lat/lon? If unavailable, the analysis will use latitude banding as a proxy.

2. **Overlap handling:** An event that falls in an elevated bin at all three bin counts (k=16, 24, 32) would be counted three times if the elevated sets are unioned naively. The combined elevated phase set should use phase range union (de-duplicated by event), not set union of three event lists. Confirm this is the right approach.

3. **IEI scope:** Inter-event intervals should be computed within the elevated-bin population only (not against the full catalog chronologically). This measures internal temporal clustering. Alternatively, IEI could be computed globally and then masked to elevated-bin events to see if they preferentially occur in shorter-interval windows. Which is more informative for the declustering question?
