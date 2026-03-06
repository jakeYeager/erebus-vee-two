### A3.B3: Ocean/Coast Sequential Threshold Sensitivity

**Source:** A2.B2  
**Reference:** A2.B3 (GCMT join refernce)


## Instructions

1. Read the source case specifications file to understand how the previous spec file was structured
2. Read the reference case specifications file to understand how the GCMT join was perfomed
3. Read the source case results file to understand how to structure a `json` results file
4. Read the expectations for this current case A3.B3
5. Create a specifications file for this case at topic-a3/spec/case-a3-b3-spec.md :
   1. It should satisfy the expectations of this case
   2. Use the data sources for A3.B3 for analysis
   3. It should incapsulate everthing required for an agent to execute the specifications outside of the main context thread

## Source case files

**Specifications file (primary):** @topic-a2/spec/case-b2-spec.md  
**Results file (primary):** @topic-a2/output/case-b2-results.json  
**Specifications file (GCMT join reference):** @topic-a2/spec/case-b3-spec.md  


## Expectations for A3.B3

**Gap or concern:**
A2.B2's result is highly sensitive to the coastline classification boundary. The oceanic subset (GSHHG primary) just misses significance at p=0.061, yet its Cramér's V (0.0276) is *larger* than the continental V (0.0245) — the significance difference is entirely a sample-size artifact, not an effect-size finding. When PB2002's broader oceanic class (absorbing back-arc and marginal-basin events) is used, the oceanic result flips to significant (p=5.24×10⁻³). The transitional zone (50–200 km offshore) is the most significant class under the primary GSHHG classification. The conclusion rests entirely on a boundary drawn at 50 km offshore, and this boundary has not been systematically tested.

**Intent:** Sequentially tighten the transitional zone boundary from 200 km offshore down to 0 km in defined incremental steps, recomputing oceanic/transitional/continental class sizes and chi-square statistics at each step. Identify the specific boundary at which the signal migrates from the transitional class into the oceanic class — or fails to — which directly quantifies the geographic extent of the mechanism. Additionally, assess whether GCMT-matched events can provide a data-driven tectonic setting attribute to supplement or replace the GSHHG proximity metric for the transitional zone definition.

**Data source requirements:** GSHHG ocean classification file, PB2002 classification file, focal mechanism join file, event lat/lon. Will require distance recalculation per event at each boundary increment. Python mapping libraries may be needed if not already present in the environment.

**Open questions:**
- What increment step for threshold tightening? (Suggest 25 km steps from 200 km → 0 km, yielding 8 test points.) *Decision: Use suggested 25km steps*
- Does the GCMT focal mechanism dataset (52.9% match rate) have sufficient geographic coverage to define subduction zone proximity as a tectonic setting metric? *Decision: Please review possible Python libs to support*

## Data sources for A3.B3

- GSHHG classification (primary): `data/iscgem/plate-location/ocean_class_gshhg_global.csv`
- PB2002 classification (coarse proxy): `data/iscgem/plate-location/ocean_class_pb2002_global.csv`
- ISC-GEM raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv`
- GCMT focal mechanism join: `data/iscgem/focal-mechanism/focal_join_global.csv`
