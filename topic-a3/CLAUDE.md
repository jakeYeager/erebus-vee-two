# Topic A3: Analysis Refinements of Topic A2 Seasonal Solar Signal Investigation

**STATUS: Planning**

## Intent

Topic A2 focused on building analysis test cases to address a liturature review on the conventionally accepted view that hydrological/snow loading is the dominant mechanism in seasonal earthquake propagation. This mechanism is suspect from this studies perspective given the nature of the "solar signal" as seen in the historical population of strong earthquakes. Hydrological forcing peaks at summer or winter solstice depending on hemisphere and regional climate — not at the equinoxes where the solar signal is strongest. The equinox signal is more consistent with a direct solar-geometric forcing (e.g., solar declination rate of change, Earth-Sun distance, or tidal stress geometry) or with a globally symmetric effect that cancels hemisphere-specific hydrological loading and leaves a residual equinox peak. This distinction is the core scientific question. This topic A3 attempts to further refine the inital cases of A2 to better describe this solar signal.

**Asymetrical Topic Workflow**
Due to the iterative review process of the cases, as compared to other topics, the operational workflow for this specific topic will override any automated scaffolding that is discribed in any rules or agents. **Stop and report to user** if any topic scaffolding subagents are triggered. Allow any subagents to freely execute individual case-related work.

## Case Table

| Case | Status   | Title                                                                       |
| ---- | -------- | --------------------------------------------------------------------------- |
| A1   | Planning | Aftershock Phase-Preference Characterization. Source: A2.A4 Sub-C           |
| B1   | Planning | Rolling-Window Chi-Square Repeat. Source: A2.B6                             |
| B3   | Planning | Ocean/Coast Sequential Threshold Sensitivity. Source: A2.B2                 |
| B4   | Planning | Depth × Magnitude Two-Way Stratification with Moho Isolation. Source: A2.B4 |
| B2   | Planning | Hemisphere Stratification Refinement. Source: A2.B1                         |
| B5   | Planning | Corrected Null-Distribution Geometric Variable Test. Source: A2.B5          |
| A2   | Planning | Aftershock Periodicity Analysis (Schuster/MFPA). Source: A2.A1, A2.A4       |
| A3   | Planning | Phase-Aware Declustering Methodology. Source: A2.A4 note                    |
| C1   | Planning | Subduction Zone Subset Test. Source: A2.B2, A2.B4                           |
| C2   | Planning | Major Sequence Removal Test. Source: A2.B6, A2.A4                           |

### Execution Tiers

**Tier 1 — Foundational (run first)**
- **A3.A1** — Most novel A2 finding; gates both A3.A2 and A3.C2
- **A3.B1** — Interval-level stationarity tracking; provides interpretive context for B2/B3/B4 and gates A3.C2

**Tier 2 — Independent (run in parallel after or alongside Tier 1)**
- **A3.B3** — Needed for A3.C1
- **A3.B4** — Needed for A3.C1
- **A3.B2** — Standalone; A3.B1 results provide interpretive context but do not block it
- **A3.B5** — Standalone methodological correction
- **A3.A3** — Exploratory probe; does not block any other case in the current list

**Tier 3 — After Tier 1 completes**
- **A3.A2** — After A3.A1
- **A3.C2** — After A3.A1 + A3.B1

**Tier 4 — After Tier 2 completes**
- **A3.C1** — After A3.B3 + A3.B4

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

