# Case A4: Declustering Sensitivity Analysis

> **Verdict: KEEP — must execute first**

## Motivation

Bradley & Hubbard (2024) specifically demonstrated that aftershock temporal clustering produces artificial tidal phase correlations, and that this artifact is large enough to explain several published positive results. The same mechanism applies directly to any periodic signal test. Park et al. (2021) showed that the standard Schuster test requires prior declustering and that failing to remove aftershocks produces false positives. This is the single most important methodological validation for Case 3A.

## Test

Apply two standard declustering algorithms (Gardner-Knopoff 1974 window method; Reasenberg 1985 spatial-temporal linking) to the catalog. Repeat the Case 3A chi-square and Rayleigh analysis on each declustered version. Compare p-values and bin distributions against the undeclustered result.

## Expected Outcome

If the equinox signal is genuine, it persists with comparable p-values across all three catalog versions (raw, GK-declustered, Reasenberg-declustered).

## Diagnostic Value

A signal that disappears after declustering would indicate the equinox structure is driven by temporal clustering of aftershock sequences — an artifact. Given the Bradley & Hubbard findings, this case is a prerequisite for trusting any result downstream.

## Key References

- Park, J., Kiraly, E., & Bourne, S. (2021). Periodic seismicity detection without declustering. *arXiv:2101.11533*.
- Bradley, K., & Hubbard, J. (2024). The great tidal earthquake hypothesis test. *Earthquake Insights* (multi-part series). https://earthquakeinsights.substack.com/p/the-great-tidal-earthquake-hypothesis
- Gardner, J. K., & Knopoff, L. (1974). Is the sequence of earthquakes in Southern California, with aftershocks removed, Poissonian? *Bulletin of the Seismological Society of America*, 64(5), 1363-1367.
- Ader, T., & Avouac, J. P. (2013). Detecting periodicities and declustering in earthquake catalogs using the Schuster spectrum. *Earth and Planetary Science Letters*, 377-378, 97-105.
