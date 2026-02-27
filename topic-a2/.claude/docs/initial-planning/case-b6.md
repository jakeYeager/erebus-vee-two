# Case B6: Rolling Window Stationarity Test

> **Verdict: KEEP — elevated priority**

## Motivation

A 72-year signal could be driven by a concentrated subset of years with anomalous clustering (non-stationarity). The overall p-value would be misleading if the effect is not present across the full time span.

**Elevation reason:** Bradley & Hubbard (2024) found that tidal correlations detected in pre-2000 data consistently failed to replicate in post-2000 data — the canonical demonstration of non-stationarity as an artifact signature. Dutilleul et al. (2021) found that Parkfield's semiannual periodicity shifted phase after the 2004 M6.0 earthquake. Both findings support stationarity testing as a critical robustness check. A signal that holds stable phase and amplitude across independent decade-long windows is far more credible than one that concentrates in a single era.

## Test

Compute the Rayleigh statistic and mean phase angle of the solar-cycle distribution using a sliding 10-year window across the 72-year record. Plot the trajectory of p-value, mean vector length, and mean phase. Test whether the equinox phase is stable across decades or drifts.

## Novel Contribution

No published annual seismicity study has specifically tested temporal stationarity of the annual signal. Stationarity of the equinox signal across the 20th–21st century record would significantly strengthen the case for a persistent physical mechanism rather than a statistical artifact of a specific time period.

## Key References

- Bradley, K., & Hubbard, J. (2024). The great tidal earthquake hypothesis test. *Earthquake Insights* (multi-part series). https://earthquakeinsights.substack.com/p/the-great-tidal-earthquake-hypothesis
- Dutilleul, P., Johnson, C. W., & Burgmann, R. (2021). Periodicity analysis of earthquake occurrence and hypocenter depth near Parkfield, California. *Geophysical Research Letters*, 48, e2020GL089673. DOI: 10.1029/2020GL089673
- Ader, T., & Avouac, J. P. (2013). Detecting periodicities and declustering in earthquake catalogs using the Schuster spectrum. *Earth and Planetary Science Letters*, 377-378, 97-105.
