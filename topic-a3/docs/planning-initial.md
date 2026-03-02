> **Status: Active**

# Topic A3 Initial Planning: Further Refinement of Topic A2 Cases

## Intent

All topic A2 cases are complete and reviewing them shows another liturature review base topic is not warrented as there are many areas to improve or expand topic A2 cases. Review the console report summaries and the inital review questions and create a `planning-refinement.md` document for us to iterate on, in building out topic A3 cases. Propose cases, their intent description, and data source requirements.

**Topic A2 console report summaries:** @topic-a3/docs/topic-a2-console-summaries.md

Only refer to the docuement if you need to know more details of the previous topics' case summaries not found in the above topic A2 console report summaries: `topic-a2/docs/topic-summary.md`

Do not look for specific whitepaper reports under A2 topic dir during this initial planning. If questions about analysis parameters, data source files, or data pipeline requirements arise list them under the proposed cases in the planning refinement document.

**Note:** Use canonical naming to refer to, and to create cases.

## Literature Informed Studies

### A2.A1: Schuster Spectrum and MFPA Periodicity Analysis

**Comparitive**
- A2.A4.4.2 sub-analysis "Post-Declustering Interval Structure" states "Interval 1 (~Mar 10–Apr 1) does not survive at k=24 under any method, though adjacent elevated bins appear near phase 0.25–0.33 in all three catalogs", but is shown in this A2.A1 study as the most robust under (~75.6 days). How are these two analysis compare/contrast and can their approach be strengthed by synthesis?
- Mainshocks tested only in this A2.A1 study, test aftershock population? A2.A4.4.3 show stronger signal in aftershock populations, how would they compare in the context of this A2.A1 framework?
  
**Structural**
- Better resolution with more bootstrap replicates? 
- Test threshold sensitivity with 2--3 other values? 3 day and 7 day?
- Bonferroni or FDR correction pos-hoc?

### A2.A2: b-Value Seasonal Variation

**Informative Refinement:**
- Any Vaulable at all until declustering methodology questions are better understood?

### A2.A3: Magnitude Stratification of the Solar Signal

**Informative Refinement:**
- Are mainshocks forcing the signal to propegate to lower magnitude bands? How to test?

**Comparitive**
- Would comparing declustered main/aftershock population be more informative, similar to A2.A4?

### A2.A4: Declustering Sensitivity Analysis
**Informative Refinement:**
- Seems to still be the central test case. Any refinements?

## Novel Studies

### A2.B1: Hemisphere Stratification — Phase Symmetry Test
**Informative Refinement:**
- Declustered population testing?

### A2.B2: Ocean vs. Continent Location — Hydrological Loading Discrimination

**Informative Refinement:**
- Leverage tectonic regime dataset from A2.B3 to find a better data-driven coastline/transitional metrics? or if tectonic regime dataset doens't have required attributes...
- Test sequentially tightened transitional metrics to determine where the signal appears?

**Structural**
- Enrich development environment with Python mapping libs? Other libs?

**Synthesis**
- After analysis refinements, create a synthesis test case with A2.B4 to attempt to isolate the signal within geological strata and region.

### A2.B3: Tectonic Regime Stratification 
**Informative Refinement:**
- Any value as too many unmatched events? Is there a way to match events in the pipeline?
- Keep as "also ran" test?

### A2.B4: Depth Stratification — Surface Loading Penetration Test

**Informative Refinement:**
- Isolate the Mohorovičić discontinuity (the Moho) as a specific band to review address p-wave and s-wave propogation noise/behaviors?
- Declustered population testing?

**Synthesis**
- After analysis refinements, create a synthesis test case with A2.B2 to attempt to isolate the signal within geological strata and region.

### A2.B5: Solar Declination Rate-of-Change vs. Position Test
**Informative Refinement:**
- Any value even with declustered population testing?

### A2.B6: Rolling Window Stationarity Test
**Informative Refinement:**
- Any value even with declustered population testing?
