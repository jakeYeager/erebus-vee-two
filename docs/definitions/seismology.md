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

## Diurnal Tide

A gravitational tidal forcing with a period of approximately **24 hours**, driven primarily by the rotation of the Earth relative to the Moon and Sun. The diurnal tide produces one high-tide and one low-tide cycle per day at a given location. In the context of seismicity, diurnal tidal stress perturbations are on the order of 1–10 kPa and have been proposed as a triggering mechanism for earthquakes on near-critically-stressed faults. The Hao et al. (2018) analysis of Japanese seismicity identified a 24-hour periodicity that correlates with the solar (not lunar) diurnal cycle, suggesting thermoelastic or atmospheric loading contributions alongside gravitational forcing.

**Timescale:** ~24 hours
**Physical origin:** Earth rotation relative to tidal bulge
**Stress amplitude:** ~1–10 kPa at crustal depths
**Contrast with annual signal:** Faults responsive to diurnal tides have nucleation timescales of seconds to minutes; faults responsive to annual loading have nucleation timescales of days to months (Bettinelli et al. 2008).

## Semidiurnal Tide

A gravitational tidal forcing with a period of approximately **12 hours**, producing two high-tide and two low-tide cycles per day. It is the dominant tidal component in most ocean basins and results from the combined gravitational pull of the Moon and Sun. In seismology, the semidiurnal solid-earth tide deforms the crust twice daily, inducing cyclic normal and shear stress changes on fault planes. Tectonic earthquakes, volcanic earthquakes at mid-ocean ridges, and slow-slip tremor have all been shown to correlate with semidiurnal tidal phase (Yan et al. 2022). The semidiurnal tide is the most commonly detected periodic signal in tidal triggering studies.

**Timescale:** ~12 hours (12.42 hours for the principal lunar M2 constituent)
**Physical origin:** Lunar and solar gravitational pull; Earth rotation
**Stress amplitude:** ~1–10 kPa at crustal depths
**Key distinction from annual signal:** Semidiurnal tidal forcing operates 700× faster than annual loading; the two mechanisms are not competing — they act on different fault populations with different nucleation timescales. The absence of a tidal-period signal in Himalayan seismicity (Ader & Avouac 2013) alongside a strong annual signal is direct evidence that nucleation timescale governs which forcing is effective.

## Fortnightly (Bi-Weekly) Tide

A secondary tidal periodicity of approximately **14 days**, arising from the modulation of the semidiurnal and diurnal tides by the synodic lunar cycle (new moon to full moon). Stress perturbations are smaller than the principal semidiurnal/diurnal tides. The fortnightly tide has been detected in tectonic earthquake catalogs (Yan et al. 2022) and is sometimes used as an intermediate-frequency control when distinguishing diurnal/semidiurnal tidal effects from longer-period seasonal effects.

**Timescale:** ~14 days
**Physical origin:** Lunar synodic cycle modulating principal tidal constituents
**Role in this project:** Intermediate-frequency control; its absence or presence helps bracket the nucleation timescale of the events under study.
