"""Sub-analysis A: Scalar signal survival on mainshock catalogs.

Computes chi-square, Rayleigh, and Cramér's V statistics for each of
four catalogs (raw, G-K mainshocks, Reasenberg mainshocks, A1b mainshocks)
at bin counts k=16, 24, 32. Also computes chi-square suppression percentages
relative to the raw catalog.
"""

import logging
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats

logger = logging.getLogger(__name__)

SOLAR_YEAR_SECS = 31_557_600.0  # Julian constant: 365.25 * 86400


def compute_phase(solar_secs: pd.Series, year_secs: float = SOLAR_YEAR_SECS) -> pd.Series:
    """Compute solar phase in [0, 1) using the Julian year constant.

    Args:
        solar_secs: Series of seconds elapsed since solar year start.
        year_secs: Length of solar year in seconds (Julian constant by default).

    Returns:
        Phase values in [0, 1).
    """
    return (solar_secs / year_secs) % 1.0


def run_sub_a_single(
    catalog: pd.DataFrame,
    catalog_name: str,
    k: int,
) -> dict[str, Any]:
    """Run Sub-analysis A statistics for a single catalog at a single bin count.

    Args:
        catalog: DataFrame with solar_secs column.
        catalog_name: Human-readable name for logging.
        k: Number of phase bins.

    Returns:
        Dict with chi2, p_chi2, cramer_v, rayleigh_R, p_rayleigh,
        mean_phase_fraction, n, k, bin_counts.
    """
    n = len(catalog)
    if n == 0:
        raise ValueError(f"Catalog {catalog_name} is empty.")

    phase = compute_phase(catalog["solar_secs"])

    # Bin assignment
    bin_idx = np.floor(phase.values * k).astype(int)
    bin_idx = np.clip(bin_idx, 0, k - 1)

    # Observed and expected counts
    obs = np.bincount(bin_idx, minlength=k).astype(float)
    expected = np.full(k, n / k)

    # Chi-square
    chi2_stat, p_chi2 = scipy.stats.chisquare(obs, expected)

    # Rayleigh statistic (mean resultant length)
    angles = 2.0 * np.pi * phase.values
    R = float(np.abs(np.mean(np.exp(1j * angles))))
    p_rayleigh = float(np.exp(-n * R ** 2))

    # Cramér's V
    cramer_v = float(np.sqrt(chi2_stat / (n * (k - 1))))

    # Mean phase angle → fraction of year
    mean_sin = float(np.mean(np.sin(angles)))
    mean_cos = float(np.mean(np.cos(angles)))
    mean_angle_rad = float(np.arctan2(mean_sin, mean_cos))
    mean_phase_fraction = float((mean_angle_rad / (2.0 * np.pi)) % 1.0)

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
    }
    logger.info(
        "%s k=%d: n=%d chi2=%.4f p=%.4e V=%.4f R=%.4f",
        catalog_name, k, n, chi2_stat, p_chi2, cramer_v, R,
    )
    return result


def run_sub_a(
    catalogs: dict[str, pd.DataFrame],
) -> dict[str, Any]:
    """Run Sub-analysis A for all catalogs and bin counts.

    Args:
        catalogs: Dict mapping catalog key → DataFrame.
                  Expected keys: "raw", "gk_mainshocks", "reas_mainshocks", "a1b_mainshocks".

    Returns:
        Sub-analysis A result dict suitable for JSON serialisation.
    """
    bin_counts = [16, 24, 32]
    key_map = {
        "raw": "raw",
        "gk_mainshocks": "gk_mainshocks",
        "reas_mainshocks": "reas_mainshocks",
        "a1b_mainshocks": "a1b_mainshocks",
    }

    results: dict[str, Any] = {}
    for key, name in key_map.items():
        df = catalogs[key]
        results[key] = {}
        for k in bin_counts:
            k_key = f"k{k}"
            results[key][k_key] = run_sub_a_single(df, name, k)

    # Suppression summary
    suppression_summary = []
    for method_key, method_label in [
        ("gk_mainshocks", "gk"),
        ("reas_mainshocks", "reas"),
        ("a1b_mainshocks", "a1b"),
    ]:
        entry: dict[str, Any] = {"method": method_label}
        for k in bin_counts:
            k_key = f"k{k}"
            chi2_raw = results["raw"][k_key]["chi2"]
            chi2_ms = results[method_key][k_key]["chi2"]
            pct = float((chi2_raw - chi2_ms) / chi2_raw * 100.0) if chi2_raw > 0 else 0.0
            entry[k_key] = {
                "chi2_raw": chi2_raw,
                "chi2_mainshock": chi2_ms,
                "chi2_suppression_pct": pct,
            }
        suppression_summary.append(entry)

    results["suppression_summary"] = suppression_summary
    return results
