"""Sub-analysis C: Aftershock population phase preference.

Computes chi-square statistics for each aftershock catalog and classifies
whether aftershocks show phase preference, and if so, whether it matches
the mainshock pattern or the A1b baseline intervals.
"""

import logging
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats

logger = logging.getLogger(__name__)

SOLAR_YEAR_SECS = 31_557_600.0  # Julian constant: 365.25 * 86400

# A1b baseline intervals — from Adhoc Case A1b (same as sub-b)
A1B_BASELINE_INTERVALS = [
    {"id": 1, "phase_start": 0.1875, "phase_end": 0.25},
    {"id": 2, "phase_start": 0.625, "phase_end": 0.656},
    {"id": 3, "phase_start": 0.875, "phase_end": 0.917},
]


def compute_phase(solar_secs: pd.Series, year_secs: float = SOLAR_YEAR_SECS) -> pd.Series:
    """Compute solar phase in [0, 1).

    Args:
        solar_secs: Series of seconds elapsed since solar year start.
        year_secs: Length of solar year in seconds.

    Returns:
        Phase values in [0, 1).
    """
    return (solar_secs / year_secs) % 1.0


def overlap_fraction(a_start: float, a_end: float, b_start: float, b_end: float) -> float:
    """Compute what fraction of interval A overlaps with interval B.

    Args:
        a_start, a_end: Bounds of interval A.
        b_start, b_end: Bounds of interval B.

    Returns:
        Overlap length / width of A.
    """
    a_width = a_end - a_start
    if a_width <= 0:
        return 0.0
    overlap = max(0.0, min(a_end, b_end) - max(a_start, b_start))
    return overlap / a_width


def get_elevated_intervals(obs: np.ndarray, k: int, expected_per_bin: float) -> list[dict]:
    """Identify elevated intervals from bin counts.

    Args:
        obs: Observed bin counts.
        k: Number of bins.
        expected_per_bin: Expected events per bin (uniform null).

    Returns:
        List of interval dicts with phase_start and phase_end.
    """
    threshold = expected_per_bin + np.sqrt(expected_per_bin)
    elevated = [i for i in range(k) if obs[i] > threshold]

    if not elevated:
        return []

    # Merge adjacent bins
    merged = []
    start = elevated[0]
    prev = elevated[0]
    for idx in elevated[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            merged.append({"phase_start": float(start / k), "phase_end": float((prev + 1) / k)})
            start = idx
            prev = idx
    merged.append({"phase_start": float(start / k), "phase_end": float((prev + 1) / k)})
    return merged


def classify_aftershock_preference(
    aftershock_elevated: list[dict],
    mainshock_elevated: list[dict],
    p_values_all_k: list[float],
) -> str:
    """Classify aftershock phase preference.

    Args:
        aftershock_elevated: Elevated intervals from aftershock catalog (at k=24).
        mainshock_elevated: Elevated intervals from corresponding mainshock catalog (at k=24).
        p_values_all_k: p-values at k=16, 24, 32.

    Returns:
        Classification string.
    """
    any_significant = any(p < 0.05 for p in p_values_all_k)
    if not any_significant:
        return "no preference"

    # Check overlap with mainshock intervals
    overlap_with_mainshocks = False
    for ae in aftershock_elevated:
        for me in mainshock_elevated:
            frac = overlap_fraction(
                ae["phase_start"], ae["phase_end"],
                me["phase_start"], me["phase_end"],
            )
            if frac > 0.5:
                overlap_with_mainshocks = True
                break

    # Check overlap with A1b baseline
    overlap_with_a1b = False
    for ae in aftershock_elevated:
        for baseline in A1B_BASELINE_INTERVALS:
            frac = overlap_fraction(
                ae["phase_start"], ae["phase_end"],
                baseline["phase_start"], baseline["phase_end"],
            )
            if frac > 0.5:
                overlap_with_a1b = True
                break

    if overlap_with_a1b:
        return "significant, intervals match A1b baseline"
    if overlap_with_mainshocks:
        return "same intervals as mainshocks"
    return "different intervals from mainshocks"


def run_sub_c_single(
    aftershock_df: pd.DataFrame,
    aftershock_name: str,
    k: int,
) -> dict[str, Any]:
    """Run Sub-analysis C statistics for a single aftershock catalog at one bin count.

    Args:
        aftershock_df: DataFrame with solar_secs column.
        aftershock_name: Human-readable name for logging.
        k: Number of phase bins.

    Returns:
        Dict with chi2, p_chi2, cramer_v, rayleigh_R, p_rayleigh, elevated_intervals.
    """
    n = len(aftershock_df)
    if n == 0:
        raise ValueError(f"Aftershock catalog {aftershock_name} is empty.")

    phase = compute_phase(aftershock_df["solar_secs"])

    bin_idx = np.floor(phase.values * k).astype(int)
    bin_idx = np.clip(bin_idx, 0, k - 1)
    obs = np.bincount(bin_idx, minlength=k).astype(float)
    expected_per_bin = n / k
    expected = np.full(k, expected_per_bin)

    chi2_stat, p_chi2 = scipy.stats.chisquare(obs, expected)

    angles = 2.0 * np.pi * phase.values
    R = float(np.abs(np.mean(np.exp(1j * angles))))
    p_rayleigh = float(np.exp(-n * R ** 2))

    cramer_v = float(np.sqrt(chi2_stat / (n * (k - 1))))

    mean_sin = float(np.mean(np.sin(angles)))
    mean_cos = float(np.mean(np.cos(angles)))
    mean_angle_rad = float(np.arctan2(mean_sin, mean_cos))
    mean_phase_fraction = float((mean_angle_rad / (2.0 * np.pi)) % 1.0)

    elevated_intervals: list[dict] = []
    if p_chi2 < 0.05:
        elevated_intervals = get_elevated_intervals(obs, k, expected_per_bin)

    result = {
        "n": int(n),
        "k": int(k),
        "chi2": float(chi2_stat),
        "p_chi2": float(p_chi2),
        "cramer_v": cramer_v,
        "rayleigh_R": R,
        "p_rayleigh": p_rayleigh,
        "mean_phase_fraction": mean_phase_fraction,
        "bin_counts": obs.astype(int).tolist(),
        "elevated_intervals": elevated_intervals,
    }
    logger.info(
        "%s k=%d: n=%d chi2=%.4f p=%.4e",
        aftershock_name, k, n, chi2_stat, p_chi2,
    )
    return result


def run_sub_c(
    aftershock_catalogs: dict[str, pd.DataFrame],
    mainshock_sub_b_results: dict[str, Any],
) -> dict[str, Any]:
    """Run Sub-analysis C for all aftershock catalogs.

    Args:
        aftershock_catalogs: Dict with keys "gk_aftershocks", "reas_aftershocks", "a1b_aftershocks".
        mainshock_sub_b_results: Sub-B results for elevated interval comparison.

    Returns:
        Sub-analysis C result dict.
    """
    bin_counts = [16, 24, 32]

    catalog_map = {
        "gk_aftershocks": ("G-K Aftershocks", "gk_mainshocks"),
        "reas_aftershocks": ("Reasenberg Aftershocks", "reas_mainshocks"),
        "a1b_aftershocks": ("A1b Aftershocks", "a1b_mainshocks"),
    }

    results: dict[str, Any] = {}

    for key, (name, mainshock_key) in catalog_map.items():
        df = aftershock_catalogs[key]
        n = len(df)
        results[key] = {"n": int(n)}

        p_values_all_k = []
        for k in bin_counts:
            k_key = f"k{k}"
            results[key][k_key] = run_sub_c_single(df, name, k)
            p_values_all_k.append(results[key][k_key]["p_chi2"])

        # Elevated intervals from aftershock at k=24 for classification
        aftershock_elevated_k24 = results[key]["k24"]["elevated_intervals"]

        # Elevated intervals from mainshock at k=24 (from sub_b recovered_intervals)
        mainshock_recovered_k24 = mainshock_sub_b_results.get(mainshock_key, {}).get(
            "k24", {}
        ).get("recovered_intervals", [])
        mainshock_elevated_k24 = [
            {"phase_start": ri["phase_start"], "phase_end": ri["phase_end"]}
            for ri in mainshock_recovered_k24
        ]

        classification = classify_aftershock_preference(
            aftershock_elevated_k24,
            mainshock_elevated_k24,
            p_values_all_k,
        )
        results[key]["classification"] = classification
        logger.info("%s classification: %s", name, classification)

    return results
