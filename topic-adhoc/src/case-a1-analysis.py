"""
Case A1: Binning Increment Sensitivity Analysis

Tests whether the solar_secs chi-square signal is robust across bin counts
(8, 16, 24, 32) using phase-normalized binning applied consistently to all
three astronomical metrics (solar_secs, lunar_secs, midnight_secs).
"""

import json
import logging
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import chisquare

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Path conventions
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent          # topic-adhoc/
PROJECT_ROOT = BASE_DIR.parent                             # erebus-vee-two/
DATA_PATH = PROJECT_ROOT / "data" / "global-sets" / "iscgem_global_events.csv"
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a1-results.json"

METRICS = ["solar_secs", "lunar_secs", "midnight_secs"]
BIN_COUNTS = [8, 16, 24, 32]
ALPHA = 0.05
N_TESTS = len(METRICS) * len(BIN_COUNTS)  # 12
ALPHA_BONFERRONI = ALPHA / N_TESTS
LUNAR_CYCLE_SECS = 29.53059 * 86400  # 2_551_442.976
DAY_SECS = 86400


def is_leap_year(year: int) -> bool:
    """Return True if *year* is a Gregorian leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def year_length_seconds(year: int) -> int:
    """Return the calendar-year length in seconds for *year*."""
    return 366 * DAY_SECS if is_leap_year(year) else 365 * DAY_SECS


def compute_phase(df: pd.DataFrame) -> pd.DataFrame:
    """Add phase columns for each metric using phase normalization.

    Solar secs can marginally exceed the calendar-year length for events
    near the solstice boundary (solaration year is astronomically defined
    and does not always align with 365/366-day calendar years).  Values
    within a small tolerance (0.002, ~17 h) are clamped to just below 1.0;
    values further out raise an error.
    """
    CLAMP_TOL = 0.01  # maximum overshoot to tolerate (~14 min for day, ~6 h for lunar)

    # solar_secs: divide by actual year length
    year_lengths = df["solaration_year"].apply(year_length_seconds).astype(float)
    df["solar_phase"] = df["solar_secs"] / year_lengths

    # lunar_secs: mean synodic month
    df["lunar_phase"] = df["lunar_secs"] / LUNAR_CYCLE_SECS

    # midnight_secs: fixed day
    df["midnight_phase"] = df["midnight_secs"] / float(DAY_SECS)

    # Validate bounds [0, 1) with clamping for marginal overshoot
    for col in ["solar_phase", "lunar_phase", "midnight_phase"]:
        below = df[df[col] < 0]
        if len(below) > 0:
            raise ValueError(
                f"Phase values < 0 for {col}. "
                f"Offending rows:\n{below.head(5)}"
            )
        over = df[df[col] >= 1.0]
        if len(over) > 0:
            max_over = df[col].max()
            if max_over >= 1.0 + CLAMP_TOL:
                raise ValueError(
                    f"Phase values significantly out of [0, 1) for {col} "
                    f"(max={max_over:.6f}). Offending rows:\n{over.head(5)}"
                )
            logger.warning(
                "%s: clamping %d events with phase in [1.0, %.6f) to 0.999999",
                col, len(over), max_over + 1e-9,
            )
            df.loc[df[col] >= 1.0, col] = 1.0 - 1e-9

    return df


def bin_phases(phases: np.ndarray, k: int) -> np.ndarray:
    """Return observed counts per bin for *k* equal-width bins on [0,1)."""
    bin_indices = np.floor(phases * k).astype(int)
    counts = np.bincount(bin_indices, minlength=k)
    return counts


def cramers_v(chi2_stat: float, n: int, k: int) -> float:
    """Compute Cramer's V for a one-dimensional chi-square test."""
    return math.sqrt(chi2_stat / (n * (k - 1)))


def run_analysis() -> dict[str, Any]:
    """Execute the full A1 analysis and return the results dict."""
    logger.info("Loading data from %s", DATA_PATH)
    df = pd.read_csv(DATA_PATH)
    n = len(df)
    logger.info("Loaded %d events", n)
    assert n == 9210, f"Expected 9210 events, got {n}"

    logger.info("Computing phase-normalized values")
    df = compute_phase(df)

    phase_columns = {
        "solar_secs": "solar_phase",
        "lunar_secs": "lunar_phase",
        "midnight_secs": "midnight_phase",
    }

    results: dict[str, Any] = {}

    for metric in METRICS:
        phase_col = phase_columns[metric]
        phases = df[phase_col].values
        metric_results: dict[str, Any] = {}

        sig_count = 0
        for k in BIN_COUNTS:
            observed = bin_phases(phases, k)
            assert observed.sum() == n, f"Bin sum mismatch for {metric} k={k}"

            chi2_stat, p_value = chisquare(observed)
            cv = cramers_v(chi2_stat, n, k)
            sig = bool(p_value < ALPHA_BONFERRONI)
            if sig:
                sig_count += 1

            metric_results[str(k)] = {
                "chi2": round(float(chi2_stat), 4),
                "p_value": float(p_value),
                "cramers_v": round(cv, 6),
                "significant_bonferroni": sig,
            }
            logger.info(
                "  %s k=%d: chi2=%.4f  p=%.6f  V=%.6f  sig=%s",
                metric, k, chi2_stat, p_value, cv, sig,
            )

        metric_results["n_significant"] = sig_count
        metric_results["robust"] = sig_count >= 3
        results[metric] = metric_results

    output: dict[str, Any] = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "data_source": "data/global-sets/iscgem_global_events.csv",
        "event_count": n,
        "methodology": {
            "phase_normalization": True,
            "normalization_factors": {
                "solar_secs": "actual_calendar_year_length_by_solaration_year",
                "lunar_secs": "mean_synodic_month_2551442.976s",
                "midnight_secs": "fixed_day_86400s",
            },
            "bin_counts": BIN_COUNTS,
            "n_tests": N_TESTS,
            "alpha_uncorrected": ALPHA,
            "alpha_bonferroni": round(ALPHA_BONFERRONI, 6),
            "robustness_criterion": "significant at >= 3 of 4 bin counts after Bonferroni correction",
        },
        "results": results,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(output, f, indent=2)
    logger.info("Results written to %s", RESULTS_PATH)

    return output


if __name__ == "__main__":
    run_analysis()
