# Final Summary (legacy "Approach Four")

## Declustered Catalog Analysis (Gardner-Knopoff Method)

**Status:** *Complete*

**Intent:** Test whether the solar_secs clustering signal identified in previous topics persists in background seismicity after removing aftershock sequences using the Gardner-Knopoff (1974) declustering algorithm. If the signal survives declustering, it strengthens the case for a global gravitational forcing mechanism rather than an artifact of aftershock clustering.

**Analysis Cases:**
- Case D0: Population description of declustered and aftershock catalogs — **Complete**
  - 36.52% of events removed as aftershocks; solar_secs non-uniformity (bins 4, 15 elevated) persists after declustering
- Case D1: Distribution uniformity on declustered solar_secs, lunar_secs, midnight_secs — **Complete**
  - All three variables show uniform distributions after declustering (all p > 0.05)
  - solar_secs signal lost significance (Approach Three p = 5.18e-05 → D1 p = 0.2437); may reflect aftershock-driven artifact or reduced statistical power
- Case D2: Inter-event interval analysis on declustered catalog — **Complete**
  - CV dropped from 1.192 to 1.069 (random regime); KS statistic dropped 71% (0.083 → 0.024)
  - Declustered catalog substantially closer to Poisson process; residual deviation likely from non-stationarity
- Case D3A: Clustering patterns — full declustered population with synthetic validation — **Complete**
  - All three variables show no clustering (all p > 0.05, synthetic percentiles 24–85%)
  - solar_secs chi2 dropped 60% (45.61 → 18.13), synthetic percentile 0.0% → 24.2% (clustering → inconclusive)
  - 1,000 uniform random synthetic catalogs confirm observed statistics are consistent with null hypothesis
- Case D3B: Clustering patterns — stratified by magnitude quartiles — **Complete**
  - All variables inconclusive in most strata; solar_secs absent in all four (17%, 38%, 61%, 42%)
  - lunar_secs shows robust pattern in Q4 (M≥6.62, p = 0.028, 3rd pct) — consistent with Approach Three 3B Q4 finding
  - midnight_secs Q3 borderline (7th pct) but inconsistent with Q4 (100th pct); likely multiple-comparison artifact
  - solar_secs magnitude-independence from Approach Three does NOT persist after declustering
- Case D4A: Energy-weighted clustering — full declustered population — **Complete**
  - Chi-square uninformative: all energy chi2 ~ 10^14–10^15, p=0.0 for real and all 1,000 synthetic catalogs
  - Max-to-min energy ratio 177,828×; effective n ≈ 9.89 events — same Gutenberg-Richter saturation as Approach Three 4A
  - Energy-weighted Rayleigh tests: all p > 0.25 — no circular concentration of seismic energy in any variable
  - D3A count-based results remain authoritative
- Case D4B: Energy-weighted clustering — stratified by magnitude quartiles — **Abandoned**
  - Not performed. D4A confirmed full Gutenberg-Richter saturation (effective n ≈ 9.9 events, p=0.0 for all real and synthetic catalogs). Stratifying an already-saturated test into smaller sub-populations would reduce effective n further and produce no interpretable metrics. Count-based D3B results are the authoritative stratified analysis.

**Cross-Case Summary (TBD):** After all cases complete, a standalone summary whitepaper will synthesize findings across D0–D4A and compare to Approach Three results. Will address the key question and assess overall signal robustness after declustering.

**Key Question:** Does the solar_secs clustering signal persist after aftershock removal? If yes, the signal reflects a genuine global modulation of independent mainshock occurrence. If no, the signal may be driven by aftershock sequences occurring preferentially at certain solar orbital positions.