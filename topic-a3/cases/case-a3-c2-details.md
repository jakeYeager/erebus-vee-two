# A3.C2: Targeted Major Sequence Phased Declustering Test

**Reference cases:** A2.B6, A2.A4, A3.B1


## Instructions

1. Read the reference case specifications file to understand how to structure a spec file.
2. Read the reference case results file to understand how to structure a `json` results file.
3. Read the expectations for this current case A3.C2, as well as the topic summary result for A3.B1.
4. Create a specifications file for this case at topic-a3/spec/case-a3-c2-spec.md :
   1. It should satisfy the expectations of this case
   2. Use the data sources for A3.C2 for analysis
   3. It should use the summary from A3.B1 to inform decisions on how to use sequence data from data sources
   4. It should incapsulate everthing required for an agent to execute the specifications outside of the main context thread


## Reference case A2.B6 files

**Specifications file:** @topic-a2/spec/case-b6-spec.md  
**Results file:** @topic-a2/output/case-b6-results.json  
**Topic Summary:** **Stop And Warn User:** *A3.B1 has not run yet! We need the results!*

## Expectations for A3.C2

**Gap or concern:**
A2.B6 identified that the most statistically significant rolling windows cluster in the 2003–2014 start-year range, contemporaneous with aftershock sequences from the 2004 Sumatra M9.1. A2.A4 showed that aftershock populations carry stronger solar-phase signals than mainshocks. Together, these findings raise a direct question the A2 framework did not test: is the global solar-phase signal diffuse across the catalog, or is it disproportionately driven by a small number of major event sequences? If removing the top 3–5 major sequences causes the global signal to collapse, the phenomenon is sequence-concentrated. If the signal survives, it is genuinely diffuse.

**Intent:** Sequentially remove the largest events in the ISC-GEM catalog and their associated aftershock sequences, recomputing the global solar-phase chi-square and interval structure after each removal:

1. Identify the 5 largest events in the catalog by magnitude (candidates: 1960 Chile M9.5, 1964 Alaska M9.2, 2004 Sumatra M9.1, 2011 Tohoku M9.0, 2010 Chile M8.8 or similar).
2. Remove each event and its aftershock sequence using a defined spatial-temporal window.
3. Recompute global chi-square, interval structure, and Cramér's V after each removal.
4. Assess whether the global signal degrades progressively, collapses at a specific removal, or remains robust.

This test should run on both the raw catalog and the declustered mainshock-only catalog, since A2.A4 showed the mainshock signal is also suppressed. If even the mainshock signal is concentrated in a few sequences, the declustering framing needs further revision.

**Report Requirement:** In the whitepaper "Results" block include a simple breakout to describe the major events removed, and basic metrics of their sequence train including: how many fore/after/complete, sequence window duration - mainshock to last aftershock, pre/post count of "last event" half-life (will inform subsequent "early" and "late" designation).

**Data source requirements:** Sequence-enriched mainshock and aftershock catalogs (G-K, Reasenberg) are **required**. `parent_id` on aftershock rows provides precise sequence membership for targeted removal without relying on a fixed spatial/temporal buffer. `foreshock_count`, `aftershock_count`, `window_secs`, and `window_km` on mainshock rows directly supply the sequence train metrics required by the Report Requirement block above, without additional derivation.


**Open questions:**
- Should aftershock windows for removal be defined using the existing G-K/Reasenberg declustering windows, or a fixed spatial/temporal buffer applied to each large event independently? *Decision: use G-K/Reasenberg declustering windows as typical in other cases*
- What magnitude threshold defines a "major sequence" for inclusion in the removal set? (Suggest M≥8.5 as a starting point, yielding approximately 10–15 events in the ISC-GEM catalog.) *Decision: Wait for the results of A3.B1 to inform, or use suggest M≥8.5 as fallback.*


## Data sources for A3.C2

- Full catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv`
- Full catalog declustered with G-K window: `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv`
- Declustered aftershocks of G-K window: `data/iscgem/declustering-algorithm/aftershocks_gk-seq_global.csv`
- Full catalog declustered with Reasenberg algorithm: `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv`
- Declustered aftershocks of  Reasenberg algorithm: `data/iscgem/declustering-algorithm/aftershocks_reas-seq_global.csv`
