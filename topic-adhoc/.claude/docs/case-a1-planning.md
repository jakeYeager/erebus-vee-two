# Topic Adhoc - Case A1 Planning

## Intent

Review the effects of binning on the three astronomical metrics: `solar_secs`, `lunar_secs` and `midnight_secs`.

**Data Source:** ISC-GEM catalog: `data/global-sets/iscgem_global_events.csv`

## Rationale

### Chi-square Results

In legacy Approach Three all three three astronomical metrics were given the statistical tests of Chi-square, Rayleigh and Cramér's V. Of the three test Chi-square provided the most informative results. For binning, each metric was binned into 16 bins

In Case 1.3.1 `solar_secs` had the following chi-square result:

> **Chi-square:** χ² = 46.07, p = 5.18 × 10⁻⁵ (**significant**)

In Case 1.3.2 `lunar_secs` had the following chi-square result:

> **Chi-square:** χ² = 28.64, p = 0.0179 (**significant**)

In Case 1.3.3 `midnight_secs` had the following chi-square result:

**Chi-square:** χ² = 16.35, p = 0.359 (**not significant**)

### lunar_secs Skew

The results were significant, however looking at the data viz they were suspicious. Subsquent tests in in legacy Approach Three Case 1B.3.1 Test A revealed that the chi-square results were skewed due to binning by max value and not accounting for variations in the synodic month.

> ### 3.1 Test A: Phase-Normalized Rebinning
>
> | Metric         | Original (abs. seconds) | Phase-Normalized    |
> | -------------- | ----------------------- | ------------------- |
> | Bin 16 count   | 517                     | 591                 |
> | Chi-square     | 28.64                   | 12.95               |
> | p-value        | 0.0179                  | 0.6061              |
> | Cramér's V     | 0.0140                  | 0.0094              |
> | Interpretation | Significant             | **Not significant** |
>
> Phase normalization increases the bin 16 count from 517 to 591 (+74 events, +14.3%), nearly eliminating the deficit. The chi-square p-value rises from 0.018 to 0.606 — well above the significance threshold. The entire distribution flattens, with the phase-normalized histogram showing no bin deviating meaningfully from the expected count of 612.6.
>
> **Conclusion: The bin 16 deficit is resolved by phase normalization.** Events that fell beyond the mean synodic month length were artificially excluded from bin 16 in the absolute-seconds binning scheme.

## Question

The `lunar_secs` solution points more to structural binning issues, however it raises the question: does binning suppress or enhance different patterns? Was 16 the lucky number for `solar_sec` or does it survive other divisions? What would the results of binning at sequential intervals look like for each metric: intervals of 8, 16, 24, 32. Perhaps these results could be informative to declustering requirements. The interval of 16 was selected to provide uniformity across metrics however beside it being base-2 divisible the increment was arbitrary.