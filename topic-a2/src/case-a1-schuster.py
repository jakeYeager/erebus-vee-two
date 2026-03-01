"""
Case A1 — Cluster-robust Schuster Spectrum Module.

Implements the cluster-robust Schuster test (Park et al. 2021) to scan for
periodicity in earthquake catalogs while accounting for temporal clustering
(aftershock sequences). Also provides the standard Schuster test for comparison.

References:
    Park, S. et al. (2021). Cluster-robust periodicity detection.
    Bradley, D. & Hubbard, J. (2024). Failure modes of the standard Schuster test.
"""

from __future__ import annotations

import logging
from typing import Dict, List

import numpy as np

logger = logging.getLogger(__name__)

# Default inter-event time threshold for defining a temporal cluster (days)
DEFAULT_DT_CLUSTER_DAYS: float = 1.0

# Period scan range for full spectrum
SPECTRUM_N_PERIODS: int = 200
SPECTRUM_MIN_DAYS: float = 0.25   # 6 hours
SPECTRUM_MAX_DAYS: float = 548.0  # ~18 months

# Explicit named periods to test
EXPLICIT_PERIODS: Dict[str, float] = {
    "tidal_12h": 0.5,
    "tidal_24h": 1.0,
    "tidal_14d": 14.77,
    "tidal_27d": 27.32,
    "third_year_122": 121.75,
    "half_year_182": 182.625,
    "annual_365": 365.25,
}


# ---------------------------------------------------------------------------
# Core statistical helpers
# ---------------------------------------------------------------------------

def _compute_d2(cos_sum: float, sin_sum: float, n: int) -> float:
    """Compute the Schuster D² statistic from summed trigonometric components.

    Parameters
    ----------
    cos_sum : float
        Sum of cos(phi_i) over all events/clusters.
    sin_sum : float
        Sum of sin(phi_i) over all events/clusters.
    n : int
        Number of events or clusters used in the sum.

    Returns
    -------
    float
        D² value (non-negative).
    """
    if n == 0:
        return 0.0
    return (cos_sum ** 2 + sin_sum ** 2) / n


def _standard_p_from_d2(d2: float) -> float:
    """Convert D² to a standard Schuster p-value assuming independent events.

    Parameters
    ----------
    d2 : float
        Schuster D² statistic.

    Returns
    -------
    float
        p-value in [0, 1].
    """
    return float(np.exp(-d2))


# ---------------------------------------------------------------------------
# Cluster assignment
# ---------------------------------------------------------------------------

def assign_clusters(
    event_time_days: np.ndarray,
    dt_cluster: float = DEFAULT_DT_CLUSTER_DAYS,
) -> np.ndarray:
    """Assign temporal cluster IDs to sorted event times.

    A new cluster begins when the inter-event time exceeds ``dt_cluster``.
    Events must be sorted in ascending order before calling this function.

    Parameters
    ----------
    event_time_days : np.ndarray
        Sorted event times expressed as decimal days (any reference epoch).
    dt_cluster : float
        Minimum inter-event gap (days) required to start a new cluster.

    Returns
    -------
    np.ndarray of int
        Array of integer cluster IDs, same length as ``event_time_days``.
    """
    n = len(event_time_days)
    cluster_ids = np.zeros(n, dtype=int)
    if n == 0:
        return cluster_ids

    current_id = 0
    for i in range(1, n):
        iet = event_time_days[i] - event_time_days[i - 1]
        if iet >= dt_cluster:
            current_id += 1
        cluster_ids[i] = current_id

    return cluster_ids


# ---------------------------------------------------------------------------
# Phase computation
# ---------------------------------------------------------------------------

def compute_phases(
    event_time_days: np.ndarray,
    period_days: float,
) -> np.ndarray:
    """Compute circular phase angles (radians) for a given test period.

    Parameters
    ----------
    event_time_days : np.ndarray
        Event times in decimal days.
    period_days : float
        Test period in days.

    Returns
    -------
    np.ndarray
        Phase angles in [0, 2π).
    """
    fractional_phase = (event_time_days % period_days) / period_days
    return 2.0 * np.pi * fractional_phase


# ---------------------------------------------------------------------------
# Single-period Schuster test
# ---------------------------------------------------------------------------

def schuster_single_period(
    event_time_days: np.ndarray,
    period_days: float,
    dt_cluster: float = DEFAULT_DT_CLUSTER_DAYS,
) -> Dict[str, float | int]:
    """Run both standard and cluster-robust Schuster tests for one test period.

    Parameters
    ----------
    event_time_days : np.ndarray
        Event times in decimal days, sorted ascending.
    period_days : float
        Test period in days.
    dt_cluster : float
        Temporal clustering threshold (days).

    Returns
    -------
    dict with keys:
        n_events, n_clusters, D2, p_standard, p_cluster_robust
    """
    n_events = len(event_time_days)
    if n_events == 0:
        return {
            "n_events": 0,
            "n_clusters": 0,
            "D2": 0.0,
            "p_standard": 1.0,
            "p_cluster_robust": 1.0,
        }

    phases = compute_phases(event_time_days, period_days)

    # Standard test — use all events
    cos_sum = float(np.sum(np.cos(phases)))
    sin_sum = float(np.sum(np.sin(phases)))
    d2_standard = _compute_d2(cos_sum, sin_sum, n_events)
    p_standard = _standard_p_from_d2(d2_standard)

    # Cluster-robust modification (Park et al. 2021)
    cluster_ids = assign_clusters(event_time_days, dt_cluster)
    n_clusters = int(cluster_ids[-1] + 1) if n_events > 0 else 0

    # Compute representative phase per cluster as the mean resultant direction.
    # Each cluster contributes exactly one unit vector at the mean-angle direction
    # (via atan2 of the mean cos/sin components). This ensures that each temporal
    # cluster receives a single vote regardless of its size, giving the test its
    # cluster-robust property. Using unit vectors (not shrunk mean vectors) means
    # D²_standard >= D²_cluster * (n_clusters/n_events) when all cluster events
    # are perfectly phase-aligned — i.e., clustering inflates D² and deflates p.
    cluster_unit_cos = np.zeros(n_clusters)
    cluster_unit_sin = np.zeros(n_clusters)
    for cid in range(n_clusters):
        mask = cluster_ids == cid
        cluster_phases = phases[mask]
        mc = np.mean(np.cos(cluster_phases))
        ms = np.mean(np.sin(cluster_phases))
        angle = np.arctan2(ms, mc)
        cluster_unit_cos[cid] = np.cos(angle)
        cluster_unit_sin[cid] = np.sin(angle)

    cos_sum_c = float(np.sum(cluster_unit_cos))
    sin_sum_c = float(np.sum(cluster_unit_sin))
    d2_cluster = _compute_d2(cos_sum_c, sin_sum_c, n_clusters)
    p_cluster_robust = _standard_p_from_d2(d2_cluster)

    return {
        "n_events": int(n_events),
        "n_clusters": int(n_clusters),
        "D2": float(d2_standard),
        "p_standard": float(p_standard),
        "p_cluster_robust": float(p_cluster_robust),
    }


# ---------------------------------------------------------------------------
# Full spectrum scan
# ---------------------------------------------------------------------------

def schuster_spectrum(
    event_time_days: np.ndarray,
    dt_cluster: float = DEFAULT_DT_CLUSTER_DAYS,
    n_periods: int = SPECTRUM_N_PERIODS,
    min_period_days: float = SPECTRUM_MIN_DAYS,
    max_period_days: float = SPECTRUM_MAX_DAYS,
) -> List[Dict]:
    """Compute the Schuster power spectrum over a log-spaced period grid.

    Parameters
    ----------
    event_time_days : np.ndarray
        Sorted event times in decimal days.
    dt_cluster : float
        Temporal clustering threshold in days.
    n_periods : int
        Number of log-spaced periods to evaluate.
    min_period_days : float
        Minimum period (days).
    max_period_days : float
        Maximum period (days).

    Returns
    -------
    list of dicts, each with keys:
        period_days, D2, p_standard, p_cluster_robust
    """
    periods = np.logspace(
        np.log10(min_period_days),
        np.log10(max_period_days),
        num=n_periods,
    )

    spectrum: List[Dict] = []
    # Pre-compute cluster IDs once for efficiency (same for all periods)
    cluster_ids = assign_clusters(event_time_days, dt_cluster)
    n_events = len(event_time_days)
    n_clusters = int(cluster_ids[-1] + 1) if n_events > 0 else 0

    for T in periods:
        phases = compute_phases(event_time_days, T)

        # Standard
        cos_sum = float(np.sum(np.cos(phases)))
        sin_sum = float(np.sum(np.sin(phases)))
        d2_std = _compute_d2(cos_sum, sin_sum, n_events)
        p_std = _standard_p_from_d2(d2_std)

        # Cluster-robust (unit-vector per cluster at mean-resultant direction)
        cluster_unit_cos = np.zeros(n_clusters)
        cluster_unit_sin = np.zeros(n_clusters)
        for cid in range(n_clusters):
            mask = cluster_ids == cid
            cp = phases[mask]
            mc = np.mean(np.cos(cp))
            ms = np.mean(np.sin(cp))
            angle = np.arctan2(ms, mc)
            cluster_unit_cos[cid] = np.cos(angle)
            cluster_unit_sin[cid] = np.sin(angle)

        cos_sum_c = float(np.sum(cluster_unit_cos))
        sin_sum_c = float(np.sum(cluster_unit_sin))
        d2_cr = _compute_d2(cos_sum_c, sin_sum_c, n_clusters)
        p_cr = _standard_p_from_d2(d2_cr)

        spectrum.append(
            {
                "period_days": float(T),
                "D2": float(d2_std),
                "p_standard": float(p_std),
                "p_cluster_robust": float(p_cr),
            }
        )

    logger.debug("Schuster spectrum computed over %d periods", n_periods)
    return spectrum


# ---------------------------------------------------------------------------
# Explicit named period tests
# ---------------------------------------------------------------------------

def schuster_explicit_tests(
    event_time_days: np.ndarray,
    dt_cluster: float = DEFAULT_DT_CLUSTER_DAYS,
    periods: Dict[str, float] | None = None,
) -> Dict[str, Dict]:
    """Run standard + cluster-robust Schuster tests on named explicit periods.

    Parameters
    ----------
    event_time_days : np.ndarray
        Sorted event times in decimal days.
    dt_cluster : float
        Temporal clustering threshold in days.
    periods : dict, optional
        Mapping of label -> period in days. Defaults to EXPLICIT_PERIODS.

    Returns
    -------
    dict mapping label -> result dict (same structure as schuster_single_period).
    """
    if periods is None:
        periods = EXPLICIT_PERIODS

    results: Dict[str, Dict] = {}
    for label, T in periods.items():
        results[label] = schuster_single_period(event_time_days, T, dt_cluster)
    return results
