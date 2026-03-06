"""
Case A3.A2 — MFPA (Modified Fourier Power Analysis) Module

Implements Dutilleul et al. (2015) MFPA over a range of test periods.
The bootstrap null uses uniform random phases in [0, 2π) — the correct null
hypothesis for testing periodicity (uniform seismicity in time at the test
period). This avoids the phase-count degeneracy issue encountered in A3.A3's
permutation baseline.
"""

import numpy as np


def mfpa_spectrum(
    event_time_days: np.ndarray,
    periods_days: np.ndarray,
    n_bootstrap: int = 10_000,
    rng_seed: int = 42,
) -> list[dict]:
    """
    Modified Fourier Power Analysis (Dutilleul et al. 2015) over a range of periods.

    Parameters
    ----------
    event_time_days : np.ndarray
        Event times in decimal days. Need not be sorted.
    periods_days : np.ndarray
        Test periods in days.
    n_bootstrap : int
        Number of bootstrap catalogs for significance estimation.
    rng_seed : int
        Random seed for bootstrap reproducibility.

    Returns
    -------
    List of dicts, one per period:
        period_days, power, p95_threshold, p99_threshold, p_mfpa,
        significant_95, significant_99
    """
    rng = np.random.default_rng(rng_seed)
    n = len(event_time_days)
    t = event_time_days

    # Pre-compute bootstrap null powers for all periods at once.
    # Each bootstrap: uniform random phases in [0, 2π).
    null_phases = rng.uniform(0, 2 * np.pi, size=(n_bootstrap, n))
    null_powers = (
        np.cos(null_phases).sum(axis=1) ** 2
        + np.sin(null_phases).sum(axis=1) ** 2
    ) / n

    # 95th and 99th percentiles of null powers (same for all periods under uniform null)
    p95 = float(np.percentile(null_powers, 95))
    p99 = float(np.percentile(null_powers, 99))

    results = []
    for T in periods_days:
        phi = 2 * np.pi * (t % T) / T
        C = np.cos(phi).sum()
        S = np.sin(phi).sum()
        power = float((C ** 2 + S ** 2) / n)
        p_mfpa = float((null_powers >= power).mean())
        results.append({
            "period_days": float(T),
            "power": power,
            "p95_threshold": p95,
            "p99_threshold": p99,
            "p_mfpa": p_mfpa,
            "significant_95": bool(power > p95),
            "significant_99": bool(power > p99),
        })

    return results


def apply_fdr_bh(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values : np.ndarray
        Array of p-values (one per test period).
    alpha : float
        FDR threshold.

    Returns
    -------
    np.ndarray of bool, same length as p_values. True = significant after FDR correction.
    """
    m = len(p_values)
    order = np.argsort(p_values)
    sorted_p = p_values[order]
    thresholds = alpha * np.arange(1, m + 1) / m
    # A period is significant if its sorted p-value is <= the BH threshold
    # AND all ranked-lower p-values are also <= their thresholds (step-up procedure)
    below = sorted_p <= thresholds
    # Find the largest k where the condition holds
    if below.any():
        max_k = np.where(below)[0].max()
        significant_sorted = np.zeros(m, dtype=bool)
        significant_sorted[:max_k + 1] = True
    else:
        significant_sorted = np.zeros(m, dtype=bool)
    # Map back to original order
    significant = np.zeros(m, dtype=bool)
    significant[order] = significant_sorted
    return significant
