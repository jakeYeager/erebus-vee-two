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

Case A3.C2 is complete. All 12 tests pass. Summary of key findings:

**Signal is diffuse — survives all 13 M≥8.5 removals**
Removing all 13 M≥8.5 events with their G-K aftershock sequences (163 events, 1.8% of catalog) only degrades chi-square from p=2.3×10⁻⁶ → p=1.3×10⁻⁴. Never approaches p=0.05. Reasenberg removal is even weaker. The signal is not sequence-concentrated.

**Interval-level detail:**
- Interval 2 (mid-August, ~0.625) is the most persistent — z≈4.3 throughout all raw-catalog removal steps
- Interval 1 (March equinox) drops below z=1.96 after step 4 in the G-K track — slightly more sequence-sensitive
- Interval 3 (late-November) shows intermediate persistence

**Mainshock-only tracks:** Major event rows themselves carry no detectable concentration — the mainshock-only removal steps produce no meaningful signal degradation in either G-K or Reasenberg declustered catalogs.

**Notable anomalies flagged in whitepaper:**
1. early_pct > 1.0 for some events due to foreshocks being included in aftershock file retrieval — flagged for A3.A1 downstream use
2. iscgem879134 (M8.60, 1960 Valdivia-2) is already a G-K/Reasenberg aftershock of the M9.55 mainshock, so its removal step removes 0 new raw-catalog rows
3. Steps 1–3 show a slight chi-square increase before declining — the Valdivia/Sumatra sequences pull slightly against the signal direction (their phase positions don't align with signal peaks)

Relationship to A3.B1: Both cases now consistently point to a diffuse, distributed signal — B1 via negative rolling-window correlation, C2 via direct targeted removal.

---

