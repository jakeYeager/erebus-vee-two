# Case A2: b-Value Seasonal Variation

> **Verdict: KEEP**

## Motivation

Colledge et al. (2025) found that in Nepal, a 20-year dataset showed no robust seasonal rate variation but did show annual variation in the Gutenberg-Richter b-value (~0.1 per kPa), peaking during the monsoon loading phase. Additionally, Ide et al. (2016) found the b-value decreases as tidal shear stress amplitude increases for very large events — establishing that b-value responds to stress perturbations even when rate does not.

## Test

Divide the 72-year catalog into solar-cycle phase bins (matching Case 3A binning). Compute the Gutenberg-Richter b-value independently for each bin using maximum likelihood estimation. Test whether b-value varies significantly across the annual solar cycle. Compute confidence intervals via bootstrap resampling.

## Expected Outcome

If literature holds: b-value peaks near solstice positions (clamping phase), troughs near equinox positions — inverse of the rate signal.

## Diagnostic Value

If rate and b-value vary inversely in phase, the equinox excess is partly a magnitude-distribution effect (more events capable of reaching M 6.0 threshold during unclamping), not only a nucleation rate effect. This would be a new finding at global M >= 6 scale.

## Key References

- Colledge, M., Chanard, K., Duverger, C., Schubnel, A., Adhikari, L. B., & Bollinger, L. (2025). Annual variations in Nepalese seismicity: b-values and seismicity rates. *Geophysical Journal International*, 242(3), 259. DOI: 10.1093/gji/ggaf259
- Bettinelli, P., Avouac, J. P., et al. (2008). Seasonal variations of seismicity and geodetic strain in the Himalaya induced by surface hydrology. *Earth and Planetary Science Letters*. https://tectonics.caltech.edu/publications/pdf/bettinelli_avouacEPSLfeb08.pdf
- Ide, S., Yabe, S., & Tanaka, Y. (2016). Earthquake potential revealed by tidal influence on earthquake size-frequency statistics. *Nature Geoscience*, 9(11), 834-837. DOI: 10.1038/ngeo2796
