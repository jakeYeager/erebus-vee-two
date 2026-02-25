# Case A1: Schuster Spectrum and MFPA Periodicity Analysis

> **Verdict: REFRAME**

## Motivation

Ader & Avouac (2013) applied the Schuster spectrum to Himalayan seismicity and detected annual periodicity with 40% amplitude while detecting no tidal-period signal. Dutilleul et al. (2015) showed MFPA resolved sub-annual harmonics (6-month, 4-month) at Parkfield alongside the annual signal.

**Reframe reason:** Bradley & Hubbard (2024) specifically identified that the standard Schuster test is susceptible to false positives from two sources relevant here: (1) aftershock temporal clustering biases tidal phase distributions, and (2) a ~20% tidal asymmetry in the astronomical forcing itself creates spurious signal. While these critiques target tidal-period tests specifically, the same aftershock contamination risk applies to any periodic test on an unclustered catalog. The standard Schuster spectrum should not be used alone.

## Revised Test

Apply the cluster-robust Schuster Spectrum Test (Park et al. 2021) rather than the standard Schuster. Use MFPA (Dutilleul et al. 2015) as the primary periodicity method, scanning 6 hours to 18 months. Run both on the raw catalog and on GK-declustered catalog (after A4) for comparison.

## Expected Outcome

If literature holds: Annual period (365.25 days) reaches significance; tidal periods (12h, 24h, 14 days) do not. Sub-annual harmonics at 6 and 4 months would confirm the bimodal equinox structure as the fundamental signal.

## Key References

- Ader, T., & Avouac, J. P. (2013). Detecting periodicities and declustering in earthquake catalogs using the Schuster spectrum. *Earth and Planetary Science Letters*, 377-378, 97-105.
- Dutilleul, P., Johnson, C. W., Burgmann, R., Wan, Y., & Shen, Z. K. (2015). Multifrequential periodogram analysis of earthquake occurrence near Parkfield, California. *Journal of Geophysical Research: Solid Earth*, 120(12), 8494-8515. DOI: 10.1002/2015JB012467
- Park, J., Kiraly, E., & Bourne, S. (2021). Periodic seismicity detection without declustering. *arXiv:2101.11533*.
- Bradley, K., & Hubbard, J. (2024). The great tidal earthquake hypothesis test. *Earthquake Insights* (multi-part series). https://earthquakeinsights.substack.com/p/the-great-tidal-earthquake-hypothesis
