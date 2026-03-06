# Data Handling Rules

## Periodic Metric Binning

All periodic astronomical metrics (`solar_secs`, `lunar_secs`, `midnight_secs`) must use **phase-normalized binning** rather than absolute-seconds binning.

**Rationale:** Absolute-seconds binning against a fixed bin count creates artificial edge effects when the underlying period is variable (e.g., synodic month length varies ~29.18–29.93 days). Phase normalization maps each event's position to [0, 1) of its respective cycle before binning, eliminating period-length artifacts. Demonstrated in Approach Three Case 1B.3.1 for `lunar_secs`; applied as standard to all metrics from Case A1 forward.

**Implementation:** Divide each event's metric value by the actual cycle length for that event's period (e.g., actual synodic month length for `lunar_secs`, actual solar year length for `solar_secs`) to obtain a normalized phase in [0, 1), then multiply by the bin count and floor to assign bin index. Use Julian constant (31,557,600 s) for solar year normalization unless directed otherwise. Use the mean synodic month: 29.53059 days × 86,400 = 2,551,442.976 seconds for lunar month normalization unless directed otherwise.

**Report disclosure:** The Methodology section of any whitepaper using binned astronomical metrics must state that phase-normalized binning was used to prevent period-length artifacts.

## Data Pipeline Declustering Algorithms

### Gardner-Knopoff (1974)

The data pipeline uses a formulaic implementation of the Gardner-Knopoff (1974) empirical formulas for magnitude-dependent space-time windows:

- **Spatial window**: `d = 10^(0.1238 * M + 0.983)` km
- **Temporal window** (M < 6.5): `t = 10^(0.5409 * M - 0.547)` days
- **Temporal window** (M >= 6.5): `t = 10^(0.032 * M + 2.7389)` days

Distances are computed using the Haversine formula (spherical Earth, radius 6371 km). This introduces a minor approximation versus the WGS84 ellipsoid -- maximum error is ~0.3% (~0.5 km at the equator for a 150 km distance). At the spatial scales of the G-K windows (tens to hundreds of km), this is negligible relative to the uncertainty in the window parameters themselves.

**Report disclosure:** The Methodology section of any whitepaper using the G-K declustering method should state that the data treatment is a formulaic implementation and not the original published table. This implementation has a known maximum error is ~0.3% for spatial calculations due to the Haversine formula. If using sequence column `aftershock_count` the G-K window for two mainshock events can have temporal/spatial overlap, and both could claim the same event; the parent has been designated as is the mainshock with the smallest `|delta_t_sec|` (temporal proximity takes priority over spatial proximity).

## Reasenberg (1985)

The data pipeline uses a pure Python implementation of the Reasenberg (1985) formula with the following defaults:

| Default   | Description |
| --------- | ----------- | -------------------------------------- |
| `rfact`   | `10`        | Interaction radius scale factor        |
| `tau-min` | `1.0`       | Minimum cluster lookback window (days) |
| `tau-max` | `10.0`      | Maximum cluster lookback window (days) |
| `p-value` | `0.95`      | Omori decay probability threshold      |
| `xmeff`   | `1.5`       | Effective magnitude threshold          |

**Report disclosure:** The Methodology section of any whitepaper using the Reasenberg declustering method should state that the data treatment is preformed with an internally managed application with reasonable defaults. When using sequence columns from mainshock datasets `window_secs` and `window_km` on each mainshock row report the actual maximum temporal and spatial reach observed across its claimed events — Reasenberg's interaction radius and adaptive lookback window vary dynamically, so these observed maximums provide the most meaningful per-mainshock footprint.
