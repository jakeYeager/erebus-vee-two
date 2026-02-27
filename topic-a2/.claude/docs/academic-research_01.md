# Topic A2 Academic Research: Annual Seasonal Modulation of Earthquake Occurrence Rates

**Research Question:** What is the evidence for annual (solar declination / seasonal) modulation of earthquake occurrence rates, and how does this differ mechanistically and statistically from diurnal or semi-diurnal solid-earth tidal triggering?

**Focus areas:** (1) thermoelastic crustal loading, (2) hydrological/snow loading as confounds, (3) hemisphere asymmetry in annual seismicity signal, (4) spectral or circular-statistics methods for detecting annual periodicity, (5) key contrast papers distinguishing tidal (diurnal) vs. seasonal (annual) timescales.

---

## 1. Hydrological and Snow Loading as the Dominant Seasonal Forcing Mechanism

The strongest body of evidence for annual modulation of seismicity attributes the signal to seasonal redistribution of surface water mass (precipitation, snowpack, groundwater) rather than to solar declination or insolation directly.

**Johnson, Fu, and Burgmann (2017)** (*Science*, 356(6343), 1161-1164; DOI: 10.1126/science.aak9547) demonstrated that:

- Seasonal water storage in the Sierra Nevada and Coast Ranges depresses mountains by approximately 1 cm in winter, with elastic rebound in summer.
- Analysis of ~3,600 earthquakes (M >= 2.0) over 2006-2015 in northern and central California showed that earthquakes occur more frequently during stress conditions that favor rupture on specific fault orientations.
- The Central San Andreas Fault near Parkfield showed increased small earthquakes in late summer/early fall (hydrological unloading phase); eastern Sierra Nevada faults showed upticks in late spring/early summer.
- The modulation effect is modest -- seasonal stress perturbations increase earthquake probability by only a few percentage points.
- The authors explicitly cautioned that no increased risk of large-magnitude earthquakes could be inferred from these low-amplitude seasonal stresses.

**Johnson, Fu, and Burgmann (2017)** (*Journal of Geophysical Research: Solid Earth*, 122(12); DOI: 10.1002/2017JB014778) extended this work by modeling seven distinct annual loading sources on California faults:

- Annual hydrospheric loading produces the largest Coulomb stress changes at 0.5-2 kPa on regional fault systems.
- An ~8 kPa seasonal mean-normal-stress perturbation modestly modulates microseismicity timing.
- Among all loading sources examined (hydrospheric, atmospheric, thermal/thermoelastic, pole tide, solid-earth tide), hydrological loading dominated seasonal stress generation.
- Faults in California are near-critically stressed and respond to subsurface pressure changes at the kilopascal level.

**Heki (2003)** (*Earth and Planetary Science Letters*, 207(1-4), 159-164) established the snow-loading mechanism for Japan:

- Snow accumulation on the western flanks of the backbone mountain range of the Japanese Islands causes seasonal crustal deformation observable via GPS.
- Past large inland earthquakes (M >= 7) in snowy regions of Japan occurred more frequently in spring and summer than fall and winter, consistent with spring-thaw unclamping of fault planes.
- Statistical significance was acknowledged to be limited due to small sample sizes of large events.

**Ueda, Kato, Johnson, and Terakawa (2024)** (*Journal of Geophysical Research: Solid Earth*; DOI: 10.1029/2023JB028217) revisited the Heki hypothesis with modern data:

- Inland seismicity beneath northeastern Japan is modestly modulated by seasonal stress changes induced by the annual snow load.
- The seasonal modulation is weaker than observed in other regions (e.g., California, Nepal).
- The weak modulation may be due to the small perturbation of seasonal stress changes relative to the long-term-averaged stressing rate, and an inability of surface unloading to unclog fluids in crustal fractures.

---

## 2. The Himalayan Case: Monsoon Loading, Seasonal Rate Variation, and Recent Reappraisal

The Nepal Himalaya provides the most dramatic reported example of seasonal seismicity, and also the most instructive cautionary tale about confounds and reassessment.

**Bollinger, Perrier, Avouac, Sapkota, Gautam, and Tiwari (2007)** (*Geophysical Research Letters*; DOI: 10.1029/2006GL029192) reported:

- Over 1995-2000, the Nepal seismic network recorded 37 +/- 8% fewer earthquakes in summer than in winter.
- For ML > 2 to ML > 4, the suppression ranged from 31% to 63%, with probability of chance occurrence less than 1%.
- The proposed mechanism is direct mass loading from the summer monsoon (water accumulation in the Indo-Gangetic Plain), which clamps the Main Himalayan Thrust.
- The authors acknowledged that part of the summer deficit may reflect poorer detection capability due to enhanced seismic noise from monsoon-driven bedload transport and landsliding.

**Bettinelli, Avouac, et al. (2008)** (*Earth and Planetary Science Letters*) provided geodetic corroboration:

- GPS measurements showed seasonal strain and stress variations of ~2-4 kPa in the Nepal Himalaya, induced by monsoon water storage.
- The seismicity rate was approximately twice as high in winter as in summer, correlated with modeled stress rate variations.
- Because Earth tides exert little influence on Himalayan seismicity, the correlated seasonal variation constrains earthquake nucleation duration to the order of days to months -- a significant implication for friction laws.

**Colledge, Chanard, Duverger, Schubnel, Adhikari, and Bollinger (2025)** (*Geophysical Journal International*, 242(3), 259; DOI: 10.1093/gji/ggaf259) significantly refined the earlier Bollinger findings:

- Analyzing 20 years of pre-2015-Gorkha Nepalese seismicity, they found no statistically robust evidence of seasonal seismicity in rate for the eastern and central Main Himalayan Thrust.
- However, the Gutenberg-Richter b-value does exhibit annual variation: b-values peak during the monsoon when Coulomb stress rates are at their minimum.
- The b-value sensitivity is approximately 0.1 per kPa of stress variation, consistent with transient loading estimates.
- The interpretation is that periodic clamping reduces the likelihood of earthquakes growing to larger magnitudes rather than simply suppressing overall rates.

---

## 3. Nonlinear Elasticity and Hydro-Seismic Coupling

**Tarantino, Poli, D'Agostino, Vassallo, Festa, Ventafridda, and Zollo (2024)** (*Nature Communications*; DOI: 10.1038/s41467-024-54094-4) demonstrated a nonlinear elastic mechanism for seasonal triggering along the Irpinia Fault, Southern Italy:

- A 14-year natural "pump-probe" experiment (2008-2021) used coda wave interferometry, GPS strain, and spring discharge monitoring.
- Seismicity peaks during maximum hydrological forcing and minimum seismic velocity.
- Velocity variations of ~0.2% were observed near karst aquifers versus ~0.03% at 25 km distance.
- Three major seismicity peaks (2009, 2013, 2021) correlated with spring discharge increases > 2 m^3/s.
- A modified quantile-quantile analysis attributed ~40% of observed seismicity to hydrological forcing.
- The mechanism is cyclical crack softening leading to failure, providing an alternative to simple Coulomb stress models.

---

## 4. Spectral and Circular-Statistics Methods for Detecting Annual Periodicity

**Ader and Avouac (2013)** (*Earth and Planetary Science Letters*, 377-378, 97-105) developed the Schuster spectrum methodology:

- The Schuster test computes p-values for the probability that event timing relative to a periodic perturbation arises from a uniform random process.
- A low Schuster p-value is necessary but not sufficient: it does not guarantee periodicity at the tested period (false positives from aftershock clustering are possible).
- Applied to Nepalese seismicity, the Schuster spectrum detected annual variations in seismicity rate with amplitude up to 40%, while no variations at any tidal period were observed.
- The Schuster spectrum also serves as a diagnostic for whether a catalog has been properly declustered.

**Dutilleul, Johnson, Burgmann, Wan, and Shen (2015)** (*Journal of Geophysical Research: Solid Earth*, 120(12), 8494-8515; DOI: 10.1002/2015JB012467) proposed an improved alternative:

- Multifrequential periodogram analysis (MFPA) has an exact modified statistic at non-Fourier frequencies, unlike the Schuster spectrum.
- For the Central San Andreas Fault near Parkfield, MFPA resolved ~1-year periodicities with additional significant periods at ~6 and ~4 months.
- Peak annual earthquake occurrence at Parkfield falls during August-November with a secondary peak in April.
- The Sierra Nevada-Eastern California Shear Zone showed a ~14-month periodic component rather than standard annual cycles.

**Park, Kiraly, and Bourne (2021)** (*arXiv:2101.11533*) addressed the declustering confound:

- Traditional Schuster-type tests require prior removal of aftershock clusters, introducing uncertainty (removing independent events or failing to remove aftershocks produces false negatives or false positives).
- Their modified Schuster Spectrum Test accounts for clustering effects directly in the statistical framework without requiring identification and removal of clustered events.
- This eliminates a major source of error in periodic seismicity detection.

---

## 5. Diurnal/Semi-Diurnal Tidal Triggering vs. Annual Periodicity: Key Contrasts

**Hao, Zhang, and Yao (2018)** (*National Science Review*, 6(5), 1016; DOI: 10.1093/nsr/nwy117):

- Analysis of 410,642 earthquakes (M >= 1.0, depth < 30 km) in Japan from June 2002 to July 2018 revealed two dominant periods of 12 and 24 hours.
- Larger earthquakes demonstrated stronger periodicity (142% variation for M >= 4.5).
- The periodicity correlates with solar rather than lunar cycles, suggesting fundamentally different triggering mechanisms than traditional tidal forcing theory.

**Yan, Chen, Sun, Xu, and Zhou (2022)** (*Geodesy and Geodynamics*, 14(1), 35-42) reviewed tidal triggering globally:

- Tectonic earthquakes show high correlation with semidiurnal and diurnal tides and 14-day (fortnightly) tides.
- Volcanic earthquakes in near-shore and mid-ocean ridge settings correlate most strongly with semidiurnal and diurnal tidal periods.
- Slow earthquake tremor duration is highly correlated with semidiurnal and diurnal tides.
- The review emphasizes that short-period tidal forcing dominates tidal triggering evidence, while annual-period mechanisms operate through fundamentally different pathways (hydrological loading, thermoelastic strain).

The critical contrast established by Bettinelli et al. (2008) and Ader et al. (2013) for the Himalaya is particularly instructive: annual modulation is clearly detected (up to 40% amplitude) while no tidal-period modulation is observed at all. This implies that the nucleation timescale for Himalayan thrust earthquakes is days-to-months, making them responsive to seasonal but not tidal stress perturbations. This timescale argument provides the key physical distinction between the two forcing regimes.

---

## Summary Table

| Finding                                                                                                    | Evidence                                                              | Relevance to Project                                                                                                  |
| ---------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| Hydrological loading dominates annual forcing on California faults (0.5-2 kPa Coulomb stress)              | Johnson et al. (2017a, 2017b)                                         | Seasonal signal has a physical mechanism; test whether global signal matches hydrological loading geography           |
| Nepal seismicity shows up to 40% annual rate variation from monsoon loading, but no tidal-period signal    | Ader et al. (2013), Bollinger et al. (2007), Bettinelli et al. (2008) | Annual signal is real but mechanism is surface loading; nucleation timescale argument constrains interpretation       |
| Recent reappraisal finds no robust seasonal rate signal in Nepal; b-value varies instead (~0.1 per kPa)    | Colledge et al. (2025)                                                | Seasonal effect may manifest in magnitude distribution rather than rate -- project should examine b-value seasonality |
| Snow loading modulates seismicity in Japan (spring thaw unclamping) but effect is modest                   | Heki (2003), Ueda et al. (2024)                                       | Hemisphere-specific confound: strongest in mid-latitude Northern Hemisphere winter                                    |
| ~40% of Irpinia seismicity attributable to seasonal hydrological forcing via nonlinear elastic mechanism   | Tarantino et al. (2024)                                               | Alternative to Coulomb stress model; crack softening may explain why seasonal modulation varies by tectonic setting   |
| Schuster spectrum and MFPA detect annual periodicity at Parkfield (peak Aug-Nov) with sub-annual harmonics | Dutilleul et al. (2015), Ader et al. (2013)                           | Methods directly applicable to project's global catalog analysis                                                      |
| Diurnal (12h, 24h) seismicity periodicity in Japan correlates with solar cycles, not lunar                 | Hao et al. (2018)                                                     | Solar-driven diurnal signal exists alongside annual signal; both may have thermoelastic components                    |
| Declustering artifacts can produce false periodic signals; cluster-robust methods now available            | Park et al. (2021)                                                    | Critical methodological requirement: project must use declustering or cluster-robust tests                            |

---

## References

- Ader, T., & Avouac, J. P. (2013). Detecting periodicities and declustering in earthquake catalogs using the Schuster spectrum. *Earth and Planetary Science Letters*, 377-378, 97-105. https://www.sciencedirect.com/science/article/pii/S0012821X13003555

- Bettinelli, P., Avouac, J. P., et al. (2008). Seasonal variations of seismicity and geodetic strain in the Himalaya induced by surface hydrology. *Earth and Planetary Science Letters*. https://tectonics.caltech.edu/publications/pdf/bettinelli_avouacEPSLfeb08.pdf

- Bollinger, L., Perrier, F., Avouac, J. P., Sapkota, S., Gautam, U., & Tiwari, D. R. (2007). Seasonal modulation of seismicity in the Himalaya of Nepal. *Geophysical Research Letters*. DOI: 10.1029/2006GL029192

- Colledge, M., Chanard, K., Duverger, C., Schubnel, A., Adhikari, L. B., & Bollinger, L. (2025). Annual variations in Nepalese seismicity: b-values and seismicity rates. *Geophysical Journal International*, 242(3), 259. DOI: 10.1093/gji/ggaf259

- Dutilleul, P., Johnson, C. W., Burgmann, R., Wan, Y., & Shen, Z. K. (2015). Multifrequential periodogram analysis of earthquake occurrence near Parkfield, California, and in the Vicinity of the San Andreas Fault. *Journal of Geophysical Research: Solid Earth*, 120(12), 8494-8515. DOI: 10.1002/2015JB012467

- Hao, J., Zhang, J., & Yao, Z. (2018). Evidence for diurnal periodicity of earthquakes from midnight to daybreak. *National Science Review*, 6(5), 1016. DOI: 10.1093/nsr/nwy117. https://pmc.ncbi.nlm.nih.gov/articles/PMC8291617/

- Heki, K. (2003). Snow load and seasonal variation of earthquake occurrence in Japan. *Earth and Planetary Science Letters*, 207(1-4), 159-164.

- Johnson, C. W., Fu, Y., & Burgmann, R. (2017a). Seasonal water storage, stress modulation, and California seismicity. *Science*, 356(6343), 1161-1164. DOI: 10.1126/science.aak9547

- Johnson, C. W., Fu, Y., & Burgmann, R. (2017b). Stress Models of the Annual Hydrospheric, Atmospheric, Thermal, and Tidal Loading Cycles on California Faults. *Journal of Geophysical Research: Solid Earth*, 122(12). DOI: 10.1002/2017JB014778

- Park, J., Kiraly, E., & Bourne, S. (2021). Periodic seismicity detection without declustering. *arXiv:2101.11533*. https://arxiv.org/abs/2101.11533

- Tarantino, S., Poli, P., D'Agostino, N., Vassallo, M., Festa, G., Ventafridda, G., & Zollo, A. (2024). Non-linear elasticity, earthquake triggering and seasonal hydrological forcing along the Irpinia fault. *Nature Communications*. DOI: 10.1038/s41467-024-54094-4. https://pmc.ncbi.nlm.nih.gov/articles/PMC11561308/

- Ueda, T., Kato, A., Johnson, C. W., & Terakawa, T. (2024). Seasonal Modulation of Crustal Seismicity in Northeastern Japan Driven by Snow Load. *Journal of Geophysical Research: Solid Earth*. DOI: 10.1029/2023JB028217

- Yan, R., Chen, X., Sun, Q., Xu, J., & Zhou, Q. (2022). A review of tidal triggering of global earthquakes. *Geodesy and Geodynamics*, 14(1), 35-42. https://www.sciopen.com/article/10.1016/j.geog.2022.06.005

---

## Request Denied Log

The following sources were attempted during research but access was blocked. These may be retrievable via institutional access, open-access repositories (arXiv, PubMed Central, ESSOAr), or Unpaywall.

| URL                                                                     | Status        | Content Sought                                                              |
| ----------------------------------------------------------------------- | ------------- | --------------------------------------------------------------------------- |
| https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023JB028217   | 403 Forbidden | Ueda et al. (2024) full text — snow load seismicity, NE Japan               |
| https://www.science.org/doi/10.1126/science.aak9547                     | 403 Forbidden | Johnson et al. (2017a) full text — seasonal water storage, California       |
| https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012467   | 403 Forbidden | Dutilleul et al. (2015) full text — MFPA, Parkfield                         |
| https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB014778   | 403 Forbidden | Johnson et al. (2017b) full text — annual loading stress models, California |
| https://www.sciencedirect.com/science/article/abs/pii/S0012821X13004500 | 403 Forbidden | Roeloffs et al. — thermoelastic strain, Parkfield borehole dilatometers     |
| https://www.sciencedirect.com/science/article/abs/pii/S0012821X02011482 | 403 Forbidden | Heki (2003) full text — snow load, Japan                                    |
| https://academic.oup.com/gji/article/242/3/ggaf259/8191266              | 403 Forbidden | Colledge et al. (2025) full text — annual b-value variations, Nepal         |
| https://link.springer.com/article/10.3103/S0747923910030114             | 303 Redirect  | Diurnal periodicity and seasonal variations                                 |
| https://link.springer.com/article/10.1007/s00024-017-1592-0             | 303 Redirect  | Periodicity of strong seismicity in Italy via Schuster spectrum             |
