# Seismology Definitions

## Bath's Law

Bath's Law states that the largest aftershock of a mainshock is typically ~1.2 magnitude units smaller than the mainshock itself, regardless of mainshock magnitude. Like GR and Omori, this is an empirical seismic law without a complete physical explanation — it is an observed regularity.

## Gutenberg-Richter Law

**Conventional definition:** The GR Law states that earthquake magnitude and frequency follow a power law: `log₁₀(N) = a - b·M`, where N is the number of earthquakes at or above magnitude M, a is a regional productivity constant, and b is the slope (~1.0 globally). It is one of the most stable empirical relationships in seismology, observed across regions, time periods, and tectonic settings.

## Gardner-Knopoff (1974) Declustering Method

The Gardner-Knopoff declustering algorithm separates mainshocks from aftershocks using magnitude-dependent space-time windows:
- For each event, a spatial radius (km) and temporal window (days) are defined based on its magnitude
- Standard lookup tables: e.g., M6.0 ≈ 40 km / 22 days; M7.0 ≈ 150 km / 64 days
- Any smaller event falling within a larger event's window is classified as an aftershock
- Remaining events form the "declustered" background catalog
- Distance calculation uses Haversine formula for great-circle distances

## Omori's Law

**Conventional definition:** Omori's Law describes the temporal decay of aftershock rate following a mainshock: `n(t) = K / (t + c)^p`, where n(t) is the aftershock rate at time t, K scales productivity, c offsets the early period, and p governs the decay rate (typically ~1.0). It is a foundational empirical law of seismology with no complete physical explanation.

## Poissonian

**Conventional definition:** A Poisson process is a memoryless random process in which events occur independently at a constant average rate. In seismology, a catalog is considered Poissonian if earthquakes occur independently — meaning the occurrence of one earthquake provides no information about when the next will occur. Declustered catalogs (aftershocks removed) are expected to approximate a Poisson process; raw catalogs deviate due to aftershock clustering governed by Omori's Law.
