# Tpoic A3 Console Summary Reports

Case A3.B1 is complete. All 17 tests pass. Summary of key findings:

**Stationarity classification: Non-stationary (unchanged from A2.B6)**
The raw catalog hits 71.0% chi-square significance (just at the 70% threshold), but the circular SD of 84.4° triggers the conflict-resolution rule — conservatively yielding non-stationary. Switching primary statistics does not change the overall classification.                             

**Interval-level results — no persistent signal:**
None of the three A1b baseline intervals is globally elevated in any catalog. Only Interval 2 (mid-August, ~0.625–0.656) qualifies as partially elevated in the raw catalog (38.7% of windows).
Chi-square significance is driven by rotating bin elevations across time, not persistent signal at the A1b positions.

**Sequence density does not explain elevated windows:**
Pearson correlation between chi-square p-values and mean aftershock counts is negative (Reasenberg r = −0.574, p < 0.001) — the most significant chi-square windows do not correspond to elevated sequence density. Sumatra-style contamination is not the driver.

**Declustering suppression cascade:**               

  ┌────────────┬─────────────────────────────┐
  │  Catalog   │ Chi-sq. windows significant │
  ├────────────┼─────────────────────────────┤
  │ Raw        │ 71.0% (44/62)               │
  ├────────────┼─────────────────────────────┤
  │ Reasenberg │ 40.3%                       │
  ├────────────┼─────────────────────────────┤
  │ G-K        │ 19.4%                       │
  ├────────────┼─────────────────────────────┤
  │ A1b        │ 6.5%                        │
  └────────────┴─────────────────────────────┘

This ladder directly parallels A2.A4's finding that aftershock sequences are themselves phase-preferring — their removal progressively collapses the detectable signal.

---

