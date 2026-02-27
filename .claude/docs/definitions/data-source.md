# Data Source

Unless specifically stated otherwise data useed in this project was collected from the USGS / FDSN Event Web Service. Collection and enrichment of astronomical metrics using the JPL DE421 ephemeris was done with https://github.com/jakeYeager/nornir-urd-two. Catalog type will be documented in case reports.

## USGS ComCat Data Quality Notes

### Known ComCat Completeness Boundaries

Global ComCat datasets exhibit two well-documented record-count inflections caused by historical changes in global seismic network infrastructure — not by real changes in earthquake activity. Any dataset spanning these boundaries should be treated as non-uniform and potentially incomplete below the thresholds noted.

### ≥M6.0 — Inflection around ~1950

The cause is the transition from the pre-WWSSN (World-Wide Standardized Seismograph Network) era to the modern instrumental era. The ISC-GEM catalog's magnitude of completeness (Mc) for shallow global events was **Mc=6.8 for 1940–1954**, only dropping to **Mc=6.5 for 1955–1963**, and **Mc=6.0 for 1964–1975**. The catalog was therefore not reliably capturing M6.0–6.8 events before the mid-1950s. The underlying driver was a postwar surge in global seismic instrumentation investment, which accelerated in the late 1950s due to nuclear test ban treaty negotiations and culminated in the deployment of the WWSSN (~120 stations in 60 countries) in the early 1960s.

**Recommended cutoff:** Use ≥M6.0 global data only from **1964** onward (full WWSSN deployment). Treat 1955–1963 as a transitional period with incomplete M6.0–6.8 coverage.

### M4.0–5.9 — Inflection around ~1973

The cause is two compounding institutional and technical changes: (1) the NEIC (National Earthquake Information Center) was reorganized into the USGS in **1973**, significantly expanding its global monitoring mandate; and (2) global seismic networks began transitioning from analog photographic paper recordings to **digital recording through the 1970s**, which lowered detection thresholds enough to reliably capture M4–5 class events globally for the first time. The M4.0–4.5 sub-band in particular sits at the edge of global catalog completeness and is most affected.

**Recommended cutoff:** Use M4.0–5.9 global data only from **1976** onward for consistent completeness. The M4.0–4.4 sub-band should be treated with additional caution even post-1976 outside North America.

---

### Summary

| Magnitude Band    | Inflection Point | Root Cause                                           | Recommended Global Cutoff                             |
| ----------------- | ---------------- | ---------------------------------------------------- | ----------------------------------------------------- |
| ≥M6.0             | ~1950            | Pre-WWSSN era Mc=6.8; postwar network expansion      | **1964** (full WWSSN deployment)                      |
| M4.0–5.9          | ~1973            | NEIC→USGS integration + network digitization         | **1976**                                              |
| M4.0–4.4 sub-band | ~1973            | As above; sits at edge of global detection threshold | Use with caution even post-1976 outside North America |

---

*Sources: ISC-GEM Catalog Completeness Analysis (USGS, 2014); USGS NEIC History; WWSSN deployment records (USGS/ARPA, early 1960s); ANSS ComCat Documentation.*