# A3.B1: Rolling-Window Chi-Square Repeat
**Source case:** A2.B6

## Instructions

1. Read the source case specifications file to understand how the previous case was constructed.
2. Read the source case results file to understand what results were gained from the analysis.
3. Read the expectations for this current case A3.B1 and use the knowledge from the source case to inform your approach to satisfy the expectations of this case.
4. Build a specification file for this case at topic-a3/spec/case-a3-b1-spec.md. It should use the same structure and pattern of deliverables as the source case specification file. Use the data sources for A3.B1 for analysis. The spec file should incapsulate everthing required for an agent to execute the specifications outside of the main context thread. 


## Source case A2.B6 files

**Specifications file:** @topic-a2/spec/case-b6-spec.md.

**Results file**: @topic-a2/output/case-b6-results.json.


## Expectations for A3.B1

**Gap or concern:**
A2.B6 used the Rayleigh test (unimodal) as its primary stationarity statistic. Chi-square (k=24) was run as secondary and reached significance in 71% of windows vs. Rayleigh's 38.7% — a 32-percentage-point divergence. The B6 console summary explicitly flagged this as "multi-modal within-window structure that the unimodal Rayleigh test misses, worth carrying forward into downstream cases." The stationarity conclusion in A2.B6 is based on Rayleigh; re-running with chi-square as primary may materially change the stationarity classification and its implications for all downstream mechanistic cases.

**Intent:** Repeat the A2.B6 rolling-window stationarity analysis with chi-square (k=24) as the primary statistic and Rayleigh as secondary. Add interval-level tracking within each window — which of the three A1b baseline intervals (Interval 1 ~0.19–0.25, Interval 2 ~0.625–0.656, Interval 3 ~0.875–0.917) are elevated — to determine whether the non-stationarity is global or interval-specific.

**Data source requirements:** Full ISC-GEM catalog with solar phase values. Sequence-enriched mainshock catalogs (G-K, Reasenberg, A1b) are **suggested** — `aftershock_count` per mainshock row enables sequence density to be computed per rolling window, providing a direct diagnostic check: windows with elevated chi-square can be tested for elevated sequence density, which would directly quantify the contribution of major aftershock sequences (e.g., 2004 Sumatra) to window-level significance.

**Open questions:**
- Should window size and stride remain unchanged from A2.B6 (5-year window, 1-year stride, 62 windows)? *Decision: apply what was used A2.B6 unless a change can materially benefit the chi-square tests*
- Should the stationarity threshold (70%) be retained or recalibrated for chi-square as primary? *Decision: apply what was used A2.B6 unless a change can materially benefit the chi-square tests*
- Should this run on declustered catalogs in addition to the raw catalog? *Decision: yes*


## Data sources for A3.B1

- Full catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv`
- Full catalog declustered with G-K window: `data/iscgem/declustering-algorithm/mainshocks_gk-seq_global.csv`
- Full catalog declustered with Reasenberg algorithm: `data/iscgem/declustering-algorithm/mainshocks_reas-seq_global.csv`
- Full catalog declustered with case A1b static window: `data/iscgem/declustering-algorithm/mainshocks_a1b-seq_global.csv`

Be aware that the mainshock data files are now enriched with cluster related columns that can be used in testing.