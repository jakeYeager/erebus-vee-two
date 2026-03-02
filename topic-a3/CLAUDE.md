# Topic A3: Analysis Refinements of Topic A2 Seasonal Solar Signal Investigation

**STATUS: Planning**

## Intent

Topic A2 focused on building analysis test cases to address a liturature review on the conventionally accepted view that hydrological/snow loading is the dominant mechanism in seasonal earthquake propagation. This mechanism is suspect from this studies perspective given the nature of the "solar signal" as seen in the historical population of strong earthquakes. Hydrological forcing peaks at summer or winter solstice depending on hemisphere and regional climate — not at the equinoxes where the solar signal is strongest. The equinox signal is more consistent with a direct solar-geometric forcing (e.g., solar declination rate of change, Earth-Sun distance, or tidal stress geometry) or with a globally symmetric effect that cancels hemisphere-specific hydrological loading and leaves a residual equinox peak. This distinction is the core scientific question. This topic A3 attempts to further refine the inital cases of A2 to better describe this solar signal.

## Data Sources
**TBD**

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

## Previous Case Analysis Reference Table

| Topic ID | Case ID | Title                                                              |
| -------- | ------- | ------------------------------------------------------------------ |
| Adhoc    | A0      | ComCat and ISC-GEM data file comparison.                           |
| Adhoc    | A0b     | Duplicate detection and cross-catalog event accounting.            |
| Adhoc    | A1      | Effects of binning increments on astronomical metrics.             |
| Adhoc    | A1b     | Elevated bin event characterization and declustering implications. |
| A2       | A1      | Schuster Spectrum and MFPA Periodicity Analysis                    |
| A2       | A2      | b-Value Seasonal Variation                                         |
| A2       | A3      | Magnitude Stratification of the Solar Signal                       |
| A2       | A4      | Declustering Sensitivity Analysis                                  |
| A2       | B1      | Hemisphere Stratification — Phase Symmetry Test                    |
| A2       | B2      | Ocean vs. Continent Location — Hydrological Loading Discrimination |
| A2       | B3      | Tectonic Regime Stratification                                     |
| A2       | B4      | Depth Stratification — Surface Loading Penetration Test            |
| A2       | B5      | Solar Declination Rate-of-Change vs. Position Test                 |
| A2       | B6      | Rolling Window Stationarity Test                                   |
