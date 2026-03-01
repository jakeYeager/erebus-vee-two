# Data Handling Rules

## Periodic Metric Binning

All periodic astronomical metrics (`solar_secs`, `lunar_secs`, `midnight_secs`) must use **phase-normalized binning** rather than absolute-seconds binning.

**Rationale:** Absolute-seconds binning against a fixed bin count creates artificial edge effects when the underlying period is variable (e.g., synodic month length varies ~29.18–29.93 days). Phase normalization maps each event's position to [0, 1) of its respective cycle before binning, eliminating period-length artifacts. Demonstrated in Approach Three Case 1B.3.1 for `lunar_secs`; applied as standard to all metrics from Case A1 forward.

**Implementation:** Divide each event's metric value by the actual cycle length for that event's period (e.g., actual synodic month length for `lunar_secs`, actual solar year length for `solar_secs`) to obtain a normalized phase in [0, 1), then multiply by the bin count and floor to assign bin index. Use Julian constant (31,557,600 s) for solar year normalization unless directed otherwise.

**Report disclosure:** The Methodology section of any whitepaper using binned astronomical metrics must state that phase-normalized binning was used and cite Case A1 as the reference for this standard.
