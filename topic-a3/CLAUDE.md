# Topic A3: Analysis Refinements of Topic A2 Seasonal Solar Signal Investigation

**STATUS: Active**

## Intent

Topic A2 focused on building analysis test cases to address a liturature review on the conventionally accepted view that hydrological/snow loading is the dominant mechanism in seasonal earthquake propagation. This mechanism is suspect from this studies perspective given the nature of the "solar signal" as seen in the historical population of strong earthquakes. Hydrological forcing peaks at summer or winter solstice depending on hemisphere and regional climate — not at the equinoxes where the solar signal is strongest. The equinox signal is more consistent with a direct solar-geometric forcing (e.g., solar declination rate of change, Earth-Sun distance, or tidal stress geometry) or with a globally symmetric effect that cancels hemisphere-specific hydrological loading and leaves a residual equinox peak. This distinction is the core scientific question. This topic A3 attempts to further refine the inital cases of A2 to better describe this solar signal.

**Asymetrical Topic Workflow**
Due to the iterative review process of the cases, as compared to other topics, the operational workflow for this specific topic will override any automated scaffolding that is discribed in any rules or agents. **Stop and report to user** if any topic scaffolding subagents are triggered. Allow any subagents to freely execute individual case-related work.


## Case Table

| Case | Status   | Title                                                        |
| ---- | -------- | ------------------------------------------------------------ |
| B1   | Complete | Rolling-Window Chi-Square Repeat                             |
| C2   | Complete | Targeted Major Sequence Phased Declustering Test             |
| B3   | Complete | Ocean/Coast Sequential Threshold Sensitivity                 |
| B4   | Complete | Depth × Magnitude Two-Way Stratification with Moho Isolation |
| B2   | Complete | Hemisphere Stratification Refinement                         |
| B5   | Complete | Corrected Null-Distribution Geometric Variable Test          |
| A3   | Complete | Phase-Concentration Audit                                    |
| A1   | Abandoned | Aftershock Phase-Preference Characterization                |
| C1   | Planning | Subduction Zone Subset Test                                  |
| A2   | Planning  | Stratified Schuster/MFPA Periodicity Audit                  |


### Execution Tiers

**Tier 1 — Foundational (run first)**
- **A3.B1** — Interval-level stationarity tracking; provides interpretive context for B2/B3/B4 and gates A3.C2

**Tier 2 — After Tier 1; independent cases may run alongside A3.C2**
- **A3.C2** — After A3.B1
- **A3.B3** — Standalone; needed for A3.C1
- **A3.B4** — Standalone; needed for A3.C1
- **A3.B2** — Standalone; A3.B1 results provide interpretive context but do not block it
- **A3.B5** — Standalone methodological correction
- **A3.A3** — Exploratory probe

**Tier 3 — After Tier 2 completes**
- **A3.C1** — After A3.B3 + A3.B4
- **A3.A2** — Standalone (A3.A1 dependency dissolved by reframe); informed by B2/B3/B4/A3 results

> **Note:** A3.A1 was abandoned — prior cases (A3.A3, A3.C2) provided sufficient coverage of the concentration-vs-diffuse and aftershock phase questions. Aftershock temporal decay analysis is deferred to a future Sequence Topic. A3.A2 was reframed as a stratified Schuster/MFPA periodicity audit; the aftershock focus and A3.A1 dependency are both dropped.


## Data Sources

### Data Pipeline Changes - Sequence Paradigm

To support the descriptive inquiry into the behavior of clusters and sequence event trains, new columns have been added to declustered datasets of mainshocks and aftershocks per the declustering algorithm. This applies to all declustered datasets available to this topic A3.

#### Aftershock output column expectations

The aftershock CSV retains all input columns plus four attribution columns appended at the end (no change from aftershock datasets common in topic A2):

| Column             | Description                                                                    |
| ------------------ | ------------------------------------------------------------------------------ |
| `parent_id`        | `usgs_id` of the mainshock whose window claimed this event                     |
| `parent_magnitude` | Magnitude of the parent mainshock                                              |
| `delta_t_sec`      | Signed elapsed seconds from the parent to this event (negative for foreshocks) |
| `delta_dist_km`    | Great-circle distance in km between this event and its parent                  |


#### Mainshock output columns

The mainshock CSV retains all input columns plus four summary columns appended at the end:

| Column             | Description                                                                           |
| ------------------ | ------------------------------------------------------------------------------------- |
| `foreshock_count`  | Number of claimed events with `delta_t_sec < 0` (occurred before the mainshock)       |
| `aftershock_count` | Number of claimed events with `delta_t_sec >= 0` (occurred at or after the mainshock) |
| `window_secs`      | Maximum `\|delta_t_sec\|` observed across all claimed events (seconds)                |
| `window_km`        | Maximum `delta_dist_km` observed across all claimed events (km)                       |

`window_secs` and `window_km` reflect the observed maximum temporal and spatial reach across all events claimed by that mainshock. For the G-K formula this equals the algorithm's theoretical window at the mainshock's magnitude. Mainshocks with no claimed events have `window_secs = 0` and `window_km = 0`.


**G-K window overlap behavior:** When two mainshock windows overlap and both could claim the same event, the parent is the mainshock with the smallest `|delta_t_sec|` (temporal proximity takes priority over spatial proximity).

**Reasenberg window overlap behavior:** For Reasenberg, each aftershock's parent is the highest-magnitude event in its cluster; no tie-breaking is required since each event belongs to exactly one cluster. `window_secs` and `window_km` on each mainshock row report the actual maximum temporal and spatial reach observed across its claimed events — Reasenberg's interaction radius and adaptive lookback window vary dynamically, so these observed maximums provide the most meaningful per-mainshock footprint.

