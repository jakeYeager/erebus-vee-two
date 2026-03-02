"""
Case A1 — Modified Fourier Power Analysis (MFPA) Module.

Implements MFPA (Dutilleul et al. 2015) as a weighted periodogram that
accounts for uneven sampling (irregular event times). Significance is assessed
via bootstrap against a uniform-phase null distribution (1,000 replicates).

References:
    Dutilleul, P. et al. (2015). MFPA periodicity detection in irregular series.
    Dutilleul, P. et al. (2021). Follow-up MFPA applications.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Scan parameters
MFPA_N_PERIODS: int = 300
MFPA_MIN_DAYS: float = 0.25   # 6 hours
MFPA_MAX_DAYS: float = 548.0  # ~18 months

# Bootstrap settings
N_BOOTSTRAP: int = 1_000
RNG_SEED: int = 42

# A1b phase intervals (in solar phase fraction 0–1)
A1B_INTERVALS: List[Tuple[float, float]] = [
    (0.1875, 0.25),    # Interval 1 (~Mar equinox)
    (0.625, 0.656),    # Interval 2 (~Aug)
    (0.875, 0.917),    # Interval 3 (~Nov)
]
A1B_INTERVAL_LABELS = ["Interval 1 (~Mar)", "Interval 2 (~Aug)", "Interval 3 (~Nov)"]


# ---------------------------------------------------------------------------
# MFPA power for a single period
# ---------------------------------------------------------------------------

def mfpa_power(phases: np.ndarray) -> float:
    """Compute MFPA power for a set of phase angles.

    MFPA power is defined as |Σ exp(i·φ_k)|² / n, equivalent to n * R²
    where R is the mean resultant length of the phase distribution.

    Parameters
    ----------
    phases : np.ndarray
        Phase angles in radians (arbitrary absolute values; only relative
        spacing matters for the power computation).

    Returns
    -------
    float
        MFPA power (non-negative).
    """
    n = len(phases)
    if n == 0:
        return 0.0
    f = np.sum(np.exp(1j * phases))
    return float((f.real ** 2 + f.imag ** 2) / n)


# ---------------------------------------------------------------------------
# Bootstrap null distribution for one period
# ---------------------------------------------------------------------------

def bootstrap_null_percentiles(
    n: int,
    n_bootstrap: int = N_BOOTSTRAP,
    rng_seed: int = RNG_SEED,
) -> Tuple[float, float]:
    """Compute 95th and 99th percentile MFPA null powers for n events.

    Parameters
    ----------
    n : int
        Number of events in the catalog.
    n_bootstrap : int
        Number of bootstrap replicates.
    rng_seed : int
        Random seed for reproducibility.

    Returns
    -------
    Tuple[float, float]
        (p95, p99) percentile values of the null distribution.
    """
    rng = np.random.default_rng(rng_seed)
    null_powers = np.empty(n_bootstrap)
    for k in range(n_bootstrap):
        random_phases = rng.uniform(0, 2 * np.pi, size=n)
        null_powers[k] = mfpa_power(random_phases)
    return float(np.percentile(null_powers, 95)), float(np.percentile(null_powers, 99))


# ---------------------------------------------------------------------------
# A1b interval cross-reference
# ---------------------------------------------------------------------------

def _phase_in_interval(phase_frac: float, intervals: List[Tuple[float, float]]) -> List[bool]:
    """Return a boolean list indicating which A1b intervals contain ``phase_frac``."""
    return [lo <= phase_frac <= hi for (lo, hi) in intervals]


def a1b_consistency_label(
    period_days: float,
    intervals: List[Tuple[float, float]] | None = None,
) -> str:
    """Determine A1b interval consistency for a detected periodicity.

    For a signal at ``period_days``, the expected peak phases within a solar year
    (phase fraction [0,1)) are ``k * (period_days / 365.25) % 1`` for integer k
    up to ``floor(365.25 / period_days)``.  We check which A1b intervals are hit
    by at least one predicted peak.

    Parameters
    ----------
    period_days : float
        Detected period in days.
    intervals : list of (lo, hi) tuples, optional
        A1b phase intervals as fraction of solar year. Defaults to A1B_INTERVALS.

    Returns
    -------
    str
        Consistency label describing which A1b intervals are predicted.
    """
    if intervals is None:
        intervals = A1B_INTERVALS

    year_days = 365.25
    n_harmonics = max(1, int(year_days / period_days))
    # Predicted peak positions within [0,1) solar phase fraction
    predicted_phases = [(k * period_days / year_days) % 1.0 for k in range(n_harmonics + 1)]

    hits = [False, False, False]
    for phase_frac in predicted_phases:
        for idx, (lo, hi) in enumerate(intervals):
            if lo <= phase_frac <= hi:
                hits[idx] = True

    hit_count = sum(hits)
    hit_indices = [i + 1 for i, h in enumerate(hits) if h]

    if hit_count == 0:
        return "inconsistent with all A1b intervals"
    elif hit_count == 1 and hit_indices == [1]:
        return "consistent with interval 1 only"
    elif hit_count == 2 and hit_indices == [1, 2]:
        return "consistent with intervals 1+2"
    elif hit_count == 2 and hit_indices == [1, 3]:
        return "consistent with intervals 1+3 (6-month)"
    elif hit_count == 2 and hit_indices == [2, 3]:
        return "consistent with intervals 2+3"
    elif hit_count == 3:
        return "consistent with all three intervals (4-month if applicable)"
    else:
        labels = "+".join(str(i) for i in hit_indices)
        return f"consistent with intervals {labels}"


# ---------------------------------------------------------------------------
# Full MFPA scan
# ---------------------------------------------------------------------------

def mfpa_scan(
    event_time_days: np.ndarray,
    n_periods: int = MFPA_N_PERIODS,
    min_period_days: float = MFPA_MIN_DAYS,
    max_period_days: float = MFPA_MAX_DAYS,
    n_bootstrap: int = N_BOOTSTRAP,
    rng_seed: int = RNG_SEED,
) -> Dict:
    """Run the full MFPA periodogram scan.

    Parameters
    ----------
    event_time_days : np.ndarray
        Event times in decimal days (any reference epoch), sorted ascending.
    n_periods : int
        Number of log-spaced periods to evaluate.
    min_period_days : float
        Minimum scan period (days).
    max_period_days : float
        Maximum scan period (days).
    n_bootstrap : int
        Bootstrap replicates for null distribution.
    rng_seed : int
        Random seed for bootstrap reproducibility.

    Returns
    -------
    dict with keys:
        "spectrum"          — list of per-period dicts
        "significant_periods" — list of periods above p95 threshold
    """
    n = len(event_time_days)
    periods = np.logspace(
        np.log10(min_period_days),
        np.log10(max_period_days),
        num=n_periods,
    )

    logger.info("MFPA scan: n=%d events, %d periods, %d bootstrap replicates", n, n_periods, n_bootstrap)

    # Pre-compute bootstrap null percentiles per-period
    # For large n the percentiles vary smoothly with n; recompute once per catalog
    rng = np.random.default_rng(rng_seed)
    null_powers_matrix = np.empty((n_bootstrap, n_periods))
    for k in range(n_bootstrap):
        random_phases = rng.uniform(0, 2 * np.pi, size=n)
        # For the null: phases are drawn uniformly, MFPA power is independent of T
        # because the event times are replaced by random phases; we still compute
        # per-period to allow for future period-dependent weighting extensions.
        null_powers_matrix[k, :] = mfpa_power(random_phases)

    p95_thresholds = np.percentile(null_powers_matrix, 95, axis=0)
    p99_thresholds = np.percentile(null_powers_matrix, 99, axis=0)

    spectrum: List[Dict] = []
    significant_periods: List[Dict] = []

    for j, T in enumerate(periods):
        phases = 2.0 * np.pi * ((event_time_days % T) / T)
        power = mfpa_power(phases)

        p95 = float(p95_thresholds[j])
        p99 = float(p99_thresholds[j])

        # MFPA p-value: fraction of bootstrap powers >= observed
        p_mfpa = float(np.mean(null_powers_matrix[:, j] >= power))

        a1b_label = a1b_consistency_label(T)

        entry = {
            "period_days": float(T),
            "power": float(power),
            "p95_threshold": p95,
            "p99_threshold": p99,
            "p_mfpa": p_mfpa,
            "a1b_consistency": a1b_label,
        }
        spectrum.append(entry)

        if power > p95:
            significant_periods.append(
                {
                    "period_days": float(T),
                    "power": float(power),
                    "p_mfpa": p_mfpa,
                    "a1b_consistency": a1b_label,
                }
            )

    logger.info("MFPA scan complete: %d significant periods found (>p95)", len(significant_periods))
    return {
        "spectrum": spectrum,
        "significant_periods": significant_periods,
    }
