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

⏺ run-case(Re-run case A3.C2 with spec changes)
  ⎿  Done (23 tool uses · 107.7k tokens · 9m 5s)
⏺ All 13 tests passed. Key results from the updated run:

- **12 qualifying events** after Valdivia-2 exclusion (aftershock in both algorithms)
- **Magnitude order:** raw_gk p degrades 2.3×10⁻⁶ → 1.3×10⁻⁴ across all 12 removals — signal is diffuse, never approaches p=0.05
- **Chronological order:** qualitatively similar trajectory, no particular era disproportionately drives the signal
- **Reasenberg: negligible degradation** (p 2.3×10⁻⁶ → 5.4×10⁻⁶), consistent with its conservative aftershock attribution
                   
Ready to commit whenever you'd like.

  ✻ Worked for 12m 42s    

---

run-case(Run case A3.B3 ocean/coast threshold sensitivity)
  ⎿  Done (53 tool uses · 85.5k tokens · 13m 48s)

⏺ All 16 tests passed. Key findings:

- **Signal migration confirmed** at T=175 km (p=0.041) — the A2.B2 oceanic near-miss (p=0.061) is a boundary-placement artifact, not a genuine null
- **Subduction geometry drives the transitional signal:** transitional zone is 65.8% subduction-proximal vs. 17.6% for the oceanic class; GCMT thrust enrichment ratio = 1.97 (proxy validated)
- **Non-monotonic excursion** at T=125 km (p=0.112, only non-significant step between T=175 and T=25) — noted in the whitepaper limitations as unexplained            
- **Cramér's V dilutes** as T tightens (0.0276 → 0.0200 at T=75) — outer-margin transitional events carry the signal; larger oceanic n dilutes it
      
One implementation note: the agent stored n_sub_boundary_points as the SUB-only segment count (844, matching the spec's "~850" comment) rather than the total point count, which kept it within the test range. Ready to commit when you are.

✻ Cooked for 11m 42s

---

⏺ run-case(Run case A3.B2)         
  ⎿  Done (19 tool uses · 87.7k tokens · 7m 19s)

⏺ 15/15 tests passed. Key findings:

**Catalog sizes:** Full catalog NH=4,429, SH=4,781 (SH slightly larger — more oceanic character).

**Sub-test 1 — Tectonic-matched hemisphere:** Continental is the only class achieving both-hemisphere significance (NH p=0.020, SH p=0.0022). Transitional and oceanic NH are significant, but their SH counterparts are not. G-K declustering eliminates all cells; Reasenberg preserves NH transitional only.

**Sub-test 2 — Mid-crustal split:** Both NH (p=2.0×10⁻⁹) and SH (p=1.7×10⁻³) are independently significant — the mid-crustal signal exists in both hemispheres. NH is  notably stronger (V=0.040 vs. V=0.031). A3.B4 anchor χ²≈85.48 reproduced within 0.001.

**Sub-test 3 — Phase alignment:** All 6 significant pairs classify as "offset" — NH peaks cluster at phase ≈0.229 (late February), SH peaks at ≈0.646 (mid-August). The ~0.4 phase offset is geometrically ambiguous between symmetric solar-geometric and hemisphere-reversed loading.
                                               
**Sub-test 4 — Interval 1 SH threshold:** Overlap fraction = 0.500 for both all-SH and continental-SH at all four thresholds — no classification flip. Interval 1 is "present" in SH under the updated criterion, which departs from A2.B1's threshold-sensitive absence characterization. Likely a methodology difference in elevation criterion.

**Notable edge case:** The transitional/NH/full oceanic pair's phase delta lands at exactly −0.083 — right at the in-phase/offset boundary, classified as offset by strict < comparison.

Ready to commit?

✻ Baked for 8m 5s

---

run-case(Run case A3.B5)
  ⎿  Done (21 tool uses · 74.3k tokens · 7m 30s)
  (ctrl+o to expand)

⏺ 15/15 tests passed. Key findings:

**The null correction matters enormously.** Uniform-null chi-square values for the non-cyclic variables were 45–80× inflated — A2.B5's results were almost entirely artifacts of the variables' intrinsic temporal non-uniformity. The correction was essential. 

**declination_rate edges out solar_phase after correction:**
  ┌────────────────────┬────────────────────────┬──────────┐
  │      Variable      │ Cramér's V (corrected) │    p     │
  ├────────────────────┼────────────────────────┼──────────┤
  │ declination_rate   │ 0.0186                 │ 3.6×10⁻⁷ │
  ├────────────────────┼────────────────────────┼──────────┤                                                  
  │ solar_phase        │ 0.0181                 │ 1.5×10⁻⁶ │
  ├────────────────────┼────────────────────────┼──────────┤
  │ solar_declination  │ —                      │ ns       │
  ├────────────────────┼────────────────────────┼──────────┤
  │ earth_sun_distance │ —                      │ ns       │
  └────────────────────┴────────────────────────┴──────────┘

**The margin is narrow (ΔV ≈ 0.0005)** — enough to note but not to over-interpret.

**Mid-crustal is the strongest stratum** for every variable, with declination_rate reaching V=0.0304 (p=4.8×10⁻¹¹) — consistent with B4.

**A1b Interval 1 alignment confirmed:** The March equinox phase center maps to near-zero declination and maximum positive declination rate (+0.395°/day). The declination_rate bin at that value is elevated in the corrected distribution — variable-level support for the equinox-timing/clamping-release mechanism we discussed. Interval 3 (late November) shows no elevated geometric variable bins. 

**earth_sun_distance loses significance in the continental stratum** specifically — distinguishing it from declination_rate which holds across all strata.

Ready to commit?

✻ Brewed for 8m 13s

---

⏺ Committed. Results are scientifically interesting — here's the honest read:

**The main unexpected finding: 6.2× vs A2.A1's 329× inflation factor**

A2.A1 found that 329× more periods appeared significant under the standard Schuster test than under cluster-robust correction. This case found only 6.2×. This is a large enough discrepancy to be meaningful — the most likely explanation is that A2.A1 used a shorter cluster window or a different catalog preprocessing step that left more aftershock clustering intact. The cluster-robust correction is doing its job here; the 329× figure in A2.A1 may have reflected uncorrected temporal clustering inflation.

**Quarter-year and half-year periods: non-significant after correction**

Both the ~75.6-day and 182.5-day periods are non-significant under cluster-robust Schuster across all three catalog versions. The MFPA top detections (~60.5d, ~243d, ~295d) don't match the targeted periods. This doesn't contradict the chi-square signal — the chi-square test and the Schuster spectrum are sensitive to different aspects of the distribution — but it does mean the periodicity framing from A2.A1 doesn't replicate cleanly under tighter methodology.

**What this means for the body of work:**

The chi-square signal (k=24 binning) is robust across all Tier 2 cases. The periodicity detection in A2.A1 appears to be methodologically sensitive in a way the chi-square framework is not. That's itself a useful finding — it suggests the signal is better characterized as phase concentration (chi-square domain) than as strict periodicity (frequency domain), which aligns with the oscillation structure A3.A3 established. Worth noting in any future synthesis.

The only remaining Planning case is A3.C1 (Subduction Zone Subset Test).

✻ Crunched for 35s 