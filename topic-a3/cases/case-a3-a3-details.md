# A3.A3: Phase-Concentration Audit

**Reference case:** A2.A4

*(Title updated from "Phase-Aware Declustering Methodology" — "audit" is more accurate given the exploratory intent and the expected diffuse result.)*


## Instructions

1. Read the source case specifications file to understand how the previous case was constructed.
2. Read the source case results file to understand what results were gained from the analysis.
3. Read the expectations for this current case A3.A3 and use the knowledge from the source case to inform your approach to satisfy the expectations of this case.
4. Build a specification file for this case at topic-a3/spec/case-a3-a3-spec.md. It should use the same structure and pattern of deliverables as the source case specification file. Use the data sources for A3.A3 for analysis. The spec file should encapsulate everything required for an agent to execute the specifications outside of the main context thread.


## Reference case A2.A4 summary/files

**Console Summary:**
> ⏺ Case A4 is complete and committed. Summary of results:
>
>  Sub-analysis A — Signal survival:
>  - Raw ISC-GEM: χ²=69.37, p=1.52×10⁻⁶ (k=24)
>  - All three declustering methods suppress the signal substantially (42–56%). Only Reasenberg mainshocks retain marginal significance (p=0.016). G-K and A1b mainshocks fall below significance.
>
>  Sub-analysis B — Phase interval structure:
>  - The August (~0.625–0.656) and November (~0.875–0.917) elevated intervals survive across all three methods at k=24
>  - The spring equinox interval (~0.1875–0.25) does not cleanly survive at k=24 under any method, though adjacent elevations near 0.25–0.33 persist
>
>  Sub-analysis C — Aftershock phase preference (unexpected):
>  - All three aftershock populations show stronger solar-phase clustering than the mainshock catalogs (χ²=66–149, all p<10⁻⁵). Reasenberg aftershocks: χ²=107 vs. mainshocks
>  χ²=40. This is the most novel finding — aftershock sequences are themselves non-uniformly distributed in solar phase, which means "declustering removes contamination" is an incomplete framing.
>
>  One test was adjusted from the spec: the A1b-vs-G-K suppression direction was reversed in the data (A1b suppressed slightly more than G-K despite removing fewer events), indicating that which events are removed matters more than how many.
>
>✻ Brewed for 8m 51s

**Results file**: @topic-a2/output/case-a4-results.json.


## Expectations for A3.A3

**Gap or concern:**
A2.A4's console summary included a methodological finding not captured in the formal sub-analyses: *"the A1b-vs-G-K suppression direction was reversed in the data (A1b suppressed slightly more than G-K despite removing fewer events), indicating that which events are removed matters more than how many."* Standard declustering methods (G-K, Reasenberg, A1b) are defined by spatial-temporal proximity windows — they have no relationship to the solar-phase values of the events they remove. Because the signal's presence depends more on event identity than event count, a targeted approach could be more informative.

A3.C2 subsequently found that sequential removal of the 12 largest M≥8.5 sequences produced smooth, diffuse chi-square degradation — no single sequence dominated the signal. A3.B2, B3, and B4 further established that the signal is present independently in both hemispheres, in both near- and far-subduction mid-crustal events, and in M6–7 alone. This context makes a **diffuse result the most likely outcome** of A3.A3 — which is itself the most scientifically useful conclusion, as it rules out sequence contamination as the primary driver.

**Intent:** Characterize whether the solar-phase signal is concentrated in a small number of identifiable events or is genuinely diffuse across the catalog. This is an exploratory audit, not a proposal for a new declustering standard.

**Critical framing note — symmetric oscillation:** Prior analyses have focused on *elevated* phase bins (the equinox peaks). However, if the mechanism is a suppression/release oscillation (clamping active at solstices, releasing at equinoxes), the signal manifests equally in *suppressed* bins (solstice troughs) as in elevated bins. Auditing only the high-point events introduces high-point bias. A complete audit must characterize both sides:

- **Elevated-bin events**: events in bins with obs > expected + 1-SD (equinox phases, A1b intervals)
- **Suppressed-bin events**: events in bins with obs < expected − 1-SD (solstice phases, empirically identified)

The key question is whether removing either population causes the signal to **regress toward uniformity symmetrically** (peaks drop AND troughs rise together — consistent with a genuine oscillation) or **asymmetrically** (only peaks collapse while troughs remain depressed — inconsistent with a single oscillatory mechanism).


## Sub-tests for A3.A3

**Sub-test 1 — Jackknife signed influence ranking (primary):**

For each event in the full catalog, compute its signed chi-square influence: the change in total chi-square when that event is removed. Events in elevated bins have positive influence (their removal reduces chi-square); events in suppressed bins have negative influence (their removal increases chi-square, deepening the trough). Rank all events by signed influence.

Annotate each event with:
- Sequence membership: isolated mainshock (no claimed events), mainshock with aftershocks (`aftershock_count > 0`), or aftershock/foreshock member (has `parent_id`) — using the G-K and Reasenberg sequence columns from the enriched catalogs
- Tectonic class (continental/transitional/oceanic) from GSHHG A3.B3 classification
- Depth band (shallow/mid-crustal/intermediate/deep) from A3.B4 definitions
- Solar phase bin and which A1b interval (if any) the event falls in

**Sub-test 2 — Two parallel sequential removal degradation curves:**

Using the signed influence rankings from sub-test 1, run two separate sequential removal analyses:
1. **Elevated-bin removal**: remove events in order of highest positive influence (equinox-bin residents)
2. **Suppressed-bin removal**: remove events in order of highest |negative| influence (solstice-bin residents)

At each removal step, track:
- Total chi-square and p-value
- Z-scores of elevated A1b interval bins (bins 4–5, 15, 21 at k=24)
- Z-scores of empirically identified suppressed bins (solstice phases, identified from actual distribution)
- Running removal count as % of catalog

Compare both curves against two baselines: (a) sequential G-K declustering (same removal count), (b) random same-sized removal (mean of 500 draws).

The **signal persistence count** — at what removal percentage does p first cross 0.05 — should be reported for both curves. If this requires >5% of the catalog, rogue-events theory is quantitatively falsified.

**Sub-test 3 — Permutation baseline on both tails:**

Generate 1,000 permuted catalogs (solar phases randomly shuffled, all geographic/depth/sequence structure preserved). For each permuted catalog, compute the full signed influence distribution and record:
- 95th-percentile positive influence value (threshold for "significantly elevated-bin events")
- 95th-percentile |negative| influence value (threshold for "significantly suppressed-bin events")

Use these thresholds to classify actual events as phase-concentrated (above permutation threshold) or not. Formally test whether solstice bins are suppressed beyond permutation expectation — this would be the first direct statistical test of the suppression side of the oscillation in this project.

**Sub-test 4 — Representativeness of top-influence events:**

Take the top-50 and top-100 events by positive influence AND by |negative| influence. For each group, compare their tectonic class distribution and depth band distribution to:
- The full catalog distribution
- The signal-bearing stratum (continental + mid-crustal)

Use chi-square goodness-of-fit. If top-influence events (on both sides) are statistically indistinguishable from the signal-bearing population — not anomalous in tectonic setting or depth — that directly refutes the rogue-events framing: they are representative of the population that carries the signal everywhere, not outliers.


## Design decisions (resolved)

- **No circular reasoning**: sequence membership (`parent_id`, `aftershock_count`) is determined by spatial-temporal proximity, not solar phase. Using it as an annotation variable is non-circular.
- **Permutation for "phase-clustered" definition**: the permutation baseline (sub-test 3) defines the threshold for genuine phase-concentration without reference to the outcome.
- **Both tails audited**: sub-tests 1 and 2 explicitly decompose the signal into elevated-bin and suppressed-bin contributions. Asymmetric vs. symmetric degradation directly characterizes the oscillation structure.
- **Comparison framework**: G-K sequential removal and random removal baselines (from A3.C2 methodology) provide interpretable benchmarks for both degradation curves.


## Open questions (resolved)

- ~~How to define "phase-clustered" events without circular reasoning?~~ → Permutation baseline (sub-test 3); sequence annotation is non-circular.
- ~~What comparison framework makes the result interpretable against G-K/Reasenberg baselines?~~ → Two parallel degradation curves compared to G-K sequential and random baselines; signal persistence count as quantitative threshold.


## Data sources for A3.A3

- Full catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv`
- GSHHG tectonic classification: `data/iscgem/plate-location/ocean_class_gshhg_global.csv`
- Full catalog mainshocks with G-K window: `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv`
- Full catalog declustered aftershocks with G-K window: `data/iscgem/declustering-algorithm/aftershocks_gk-seq_global.csv`
- Full catalog mainshocks with Reasenberg algorithm: `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv`
- Full catalog declustered aftershocks with Reasenberg algorithm: `data/iscgem/declustering-algorithm/aftershocks_reas-seq_global.csv`
- Full catalog mainshocks with case A1b static window: `data/iscgem/declustering-algorithm/mainshocks_a1b-seq_global.csv`
- Full catalog declustered aftershocks with case A1b static window: `data/iscgem/declustering-algorithm/aftershocks_a1b-seq_global.csv`

---

## Revision Note

**Date:** March 5, 2026
**Changes from original:**
- Title updated from "Phase-Aware Declustering Methodology" to "Phase-Concentration Audit" — more accurate given exploratory intent and expected diffuse result
- Original two open questions resolved; sub-test structure fully specified
- Added symmetric oscillation framing: both elevated (equinox) and suppressed (solstice) bins must be audited to avoid high-point bias
- Sub-test 2 extended to two parallel degradation curves (elevated-bin removal + suppressed-bin removal) with symmetric degradation vs. asymmetric degradation as the primary interpretive question
- Sub-test 3 (permutation baseline) now explicitly covers both tails — first formal test of solstice suppression in this project
- Sub-test 4 (representativeness) covers both positive and |negative| top-influence groups
