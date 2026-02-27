# Topic Adhoc: Summary

## Completed Cases

### Case A0: ComCat and ISC-GEM Catalog Comparison Reference

**Status:** Complete

**Key Results:**
- ComCat contains 9,802 events; ISC-GEM contains 9,210 events (592 difference, +6.4%)
- ComCat is a hybrid catalog: 30.3% of records (2,973) carry ISC-GEM-prefixed IDs, concentrated pre-1980
- Magnitude precision differs substantially: ComCat 77.5% one-decimal vs ISC-GEM 82.9% two-decimal
- Both catalogs span 1950--2021 with identical 10-column schemas

### Case A0b: Duplicate Detection and Cross-Catalog Event Accounting

**Status:** Complete

**Key Results:**
- Within-ComCat duplication is negligible: only 1 duplicate candidate pair found (0.03% of iscgem-prefixed subset)
- Cross-catalog proximity matching: 8,413 matched, 1,389 ComCat-only, 797 ISC-GEM-only (primary tolerance)
- ComCat-only events are 92.5% concentrated in M 6.0--6.4, consistent with rounding artifacts
- The 592-event gap is entirely catalog divergence, not internal duplication
- ISC-GEM-only events spike in the 1970s (414 of 797), suggesting systematic catalog construction differences
- Pipeline recommendations R1--R6 documented; deduplication step (R4) not warranted

### Case A1: Binning Increment Sensitivity Analysis for Astronomical Metrics

**Status:** Complete

**Key Results:**
- `solar_secs` chi-square signal is robust: Bonferroni-significant at 3 of 4 bin counts (k=16, 24, 32); only k=8 falls short (p=0.015)
- `lunar_secs` shows no significant departure from uniformity at any bin count after phase normalization (all p > 0.18), confirming the legacy 16-bin result was a binning artifact
- `midnight_secs` remains non-significant across all bin counts (all p > 0.62), confirming expected behavior as a temporal control at global scale
- Cramer's V for `solar_secs` is small but consistent (0.016--0.018), indicating a modest effect size
- Phase normalization established as project standard; clamping required for 3 solar and 23 lunar boundary events

### Case A1b: Elevated Bin Event Characterization and Declustering Implications

**Status:** Complete

**Key Results:**
- Combined elevated phase set (solar phase intervals covered by top-3 bins in >=2 of 3 k values): 3 intervals at [0.1875, 0.25), [0.625, 0.656), [0.875, 0.917) — approximately early March, mid-August, and mid-November
- Elevated-bin population: 1,438 events (15.6% of catalog vs. 13.5% expected under null)
- Elevated-bin median IEI: 1.46 days — far below null 95% CI of 11.4–12.8 days, indicating strong temporal clustering consistent with aftershock sequence contamination
- Elevated-bin median NN distance: 35.1 km — below null 95% CI of 37.8–45.3 km, indicating modestly tighter spatial clustering than random
- Boundary proximity distribution nearly identical between elevated-bin set and full catalog; magnitude distribution also shows no skew toward lower magnitudes
- Proposed data-informed declustering reference window: 83.2 km spatial / 95.6 days temporal (p75 NN and p90 IEI of elevated-bin population); below G-K M6.0 spatial window (49 km at p50) but temporal clustering far tighter than any G-K temporal window
