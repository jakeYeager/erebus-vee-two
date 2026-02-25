# Case B5: Solar Declination Rate-of-Change vs. Position Test

> **Verdict: KEEP**

## Motivation

Case 3A encoded seismic events by their position in the solar annual cycle (fraction of year from perihelion). The bimodal equinox structure could reflect two distinct solar variables that both peak at the equinoxes:

- **(a)** The rate of change of solar declination is maximum at the equinoxes.
- **(b)** The Earth-Sun distance is near mean value at equinoxes (perihelion is in January, aphelion in July).

These variables are correlated but not identical and have different physical implications.

## Test

Re-encode each event using three separate phase variables:

1. Solar declination angle (ranging –23.5° to +23.5°)
2. Rate of change of declination (degrees/day)
3. Earth-Sun distance (AU)

Compute bin distributions for each. Identify which variable produces the most coherent, highest-significance clustering.

## Novel Contribution

No published seismicity study has explicitly tested solar declination rate-of-change as a potential trigger variable, nor compared it to solar position and Earth-Sun distance within the same catalog.

## Key References

- Ader, T., & Avouac, J. P. (2013). Detecting periodicities and declustering in earthquake catalogs using the Schuster spectrum. *Earth and Planetary Science Letters*, 377-378, 97-105.
- Dutilleul, P., Johnson, C. W., Burgmann, R., Wan, Y., & Shen, Z. K. (2015). Multifrequential periodogram analysis of earthquake occurrence near Parkfield, California. *Journal of Geophysical Research: Solid Earth*, 120(12), 8494-8515. DOI: 10.1002/2015JB012467
