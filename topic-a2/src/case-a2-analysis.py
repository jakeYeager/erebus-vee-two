"""
Case A2: b-Value Seasonal Variation — Main Analysis Script.

Divides the ISC-GEM catalog into solar-cycle phase bins and computes the
Gutenberg-Richter b-value independently for each bin using maximum likelihood
estimation (MLE), testing whether b-value varies significantly across the
annual solar cycle.

Usage:
    python topic-a2/src/case-a2-analysis.py

Outputs:
    topic-a2/output/case-a2-results.json
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats

# ---------------------------------------------------------------------------
# Project path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("case-a2-analysis")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
JULIAN_YEAR_SECS: float = 31_557_600.0   # Julian constant (confirmed pre-run)
MC: float = 6.0                           # Completeness magnitude
BIN_COUNTS: list[int] = [24, 32]          # Primary=24, secondary=32 per Adhoc A1
BOOTSTRAP_RESAMPLES: int = 1000
RANDOM_SEED: int = 42
LOW_N_THRESHOLD: int = 20

# Equinox and solstice reference phases (fraction of year)
EQUINOX_PHASES: list[float] = [0.19, 0.69]   # March, September
SOLSTICE_PHASES: list[float] = [0.44, 0.94]  # June, December
PHASE_CLASS_WINDOW: float = 0.1

# Data paths
RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
OUTPUT_PATH = BASE_DIR / "output" / "case-a2-results.json"

EXPECTED_N: int = 9210


# ---------------------------------------------------------------------------
# Phase normalization
# ---------------------------------------------------------------------------

def compute_solar_phase(solar_secs: pd.Series) -> pd.Series:
    """Compute normalized solar phase in [0, 1) using Julian year constant.

    Parameters
    ----------
    solar_secs : pd.Series
        Seconds elapsed since start of solar year for each event.

    Returns
    -------
    pd.Series
        Phase values in [0, 1).
    """
    return (solar_secs / JULIAN_YEAR_SECS) % 1.0


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_catalog(path: Path) -> pd.DataFrame:
    """Load ISC-GEM catalog, validate row count and magnitude column.

    Parameters
    ----------
    path : Path
        Filesystem path to the CSV file.

    Returns
    -------
    pd.DataFrame
        Loaded catalog.

    Raises
    ------
    AssertionError
        If row count != 9210 or usgs_mag contains NaN.
    """
    logger.info("Loading ISC-GEM catalog from %s", path)
    df = pd.read_csv(path)
    n = len(df)
    logger.info("  Loaded %d rows (expected %d)", n, EXPECTED_N)
    assert n == EXPECTED_N, f"Expected {EXPECTED_N} rows, got {n}"
    assert df["usgs_mag"].notna().all(), "usgs_mag contains NaN values"
    logger.info("  usgs_mag: no NaN values; range [%.2f, %.2f]",
                df["usgs_mag"].min(), df["usgs_mag"].max())
    return df


# ---------------------------------------------------------------------------
# b-Value MLE
# ---------------------------------------------------------------------------

def compute_b_mle(mags: np.ndarray, mc: float = MC) -> float:
    """Compute Gutenberg-Richter b-value using Aki (1965) MLE formula.

    Parameters
    ----------
    mags : np.ndarray
        Array of magnitudes for events in this bin.
    mc : float
        Completeness magnitude.

    Returns
    -------
    float
        MLE b-value estimate.
    """
    mean_mag = np.mean(mags)
    b = np.log10(np.e) / (mean_mag - mc)
    return float(b)


def compute_b_se(b_mle: float, n: int) -> float:
    """Compute Shi & Bolt (1982) standard error for b-value MLE.

    Parameters
    ----------
    b_mle : float
        Estimated b-value.
    n : int
        Number of events in sample.

    Returns
    -------
    float
        Standard error of b_mle.
    """
    return float(b_mle / np.sqrt(n))


def bootstrap_b_ci(
    mags: np.ndarray,
    mc: float = MC,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    rng: np.random.Generator | None = None,
) -> tuple[float, float]:
    """Compute 95% bootstrap confidence interval for b-value MLE.

    Parameters
    ----------
    mags : np.ndarray
        Array of magnitudes.
    mc : float
        Completeness magnitude.
    n_resamples : int
        Number of bootstrap resamples.
    rng : np.random.Generator or None
        Random number generator; if None, uses default.

    Returns
    -------
    tuple[float, float]
        (lower_2.5th_percentile, upper_97.5th_percentile)
    """
    if rng is None:
        rng = np.random.default_rng(RANDOM_SEED)
    n = len(mags)
    b_boots = np.empty(n_resamples)
    for i in range(n_resamples):
        sample = rng.choice(mags, size=n, replace=True)
        b_boots[i] = np.log10(np.e) / (np.mean(sample) - mc)
    lower = float(np.percentile(b_boots, 2.5))
    upper = float(np.percentile(b_boots, 97.5))
    return lower, upper


# ---------------------------------------------------------------------------
# Phase classification
# ---------------------------------------------------------------------------

def classify_phase(phase_center: float) -> str:
    """Classify a phase center as near_solstice, near_equinox, or other.

    Parameters
    ----------
    phase_center : float
        Phase fraction in [0, 1).

    Returns
    -------
    str
        One of "near_solstice", "near_equinox", "other".
    """
    for s in SOLSTICE_PHASES:
        if abs(phase_center - s) < PHASE_CLASS_WINDOW:
            return "near_solstice"
    for e in EQUINOX_PHASES:
        if abs(phase_center - e) < PHASE_CLASS_WINDOW:
            return "near_equinox"
    return "other"


# ---------------------------------------------------------------------------
# Per-bin b-value analysis
# ---------------------------------------------------------------------------

def analyze_bins(
    df: pd.DataFrame,
    k: int,
    bin_col: str,
) -> dict[str, Any]:
    """Compute b-value MLE, SE, bootstrap CI, and significance tests for all bins.

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with ``usgs_mag`` and a bin-index column ``bin_col``.
    k : int
        Number of phase bins.
    bin_col : str
        Column name containing bin index (0 to k-1).

    Returns
    -------
    dict
        Per-bin statistics and phase-variation test results.
    """
    logger.info("Analyzing b-values for k=%d bins", k)

    # Seed the bootstrap RNG once for this k
    rng = np.random.default_rng(RANDOM_SEED)

    bins_data: list[dict[str, Any]] = []
    groups_for_anova: list[np.ndarray] = []
    counts: list[int] = []
    b_values: list[float] = []

    for j in range(k):
        mask = df[bin_col] == j
        mags_j = df.loc[mask, "usgs_mag"].values
        n_j = len(mags_j)
        phase_center = (j + 0.5) / k
        low_n_flag = n_j < LOW_N_THRESHOLD

        if n_j == 0:
            logger.warning("Bin %d has 0 events; skipping b-value computation", j)
            bins_data.append({
                "bin_idx": j,
                "phase_center": float(phase_center),
                "n": 0,
                "b_mle": None,
                "b_se": None,
                "b_ci95_lower": None,
                "b_ci95_upper": None,
                "low_n_flag": True,
            })
            groups_for_anova.append(mags_j)
            counts.append(0)
            b_values.append(np.nan)
            continue

        b_j = compute_b_mle(mags_j)
        se_j = compute_b_se(b_j, n_j)
        ci_lower, ci_upper = bootstrap_b_ci(mags_j, rng=rng)

        logger.info(
            "  Bin %2d: n=%3d, phase=%.3f, b=%.4f ± %.4f, CI95=[%.4f, %.4f]%s",
            j, n_j, phase_center, b_j, se_j, ci_lower, ci_upper,
            " [LOW-N]" if low_n_flag else "",
        )

        bins_data.append({
            "bin_idx": j,
            "phase_center": float(phase_center),
            "n": int(n_j),
            "b_mle": float(b_j),
            "b_se": float(se_j),
            "b_ci95_lower": float(ci_lower),
            "b_ci95_upper": float(ci_upper),
            "low_n_flag": bool(low_n_flag),
        })
        groups_for_anova.append(mags_j)
        counts.append(int(n_j))
        b_values.append(float(b_j))

    # Aggregate statistics (excluding bins with no events)
    valid_b = [b for b in b_values if not np.isnan(b)]
    b_mean = float(np.mean(valid_b))
    b_std = float(np.std(valid_b))
    b_max_val = float(np.max(valid_b))
    b_min_val = float(np.min(valid_b))
    b_range = float(b_max_val - b_min_val)

    b_arr = np.array(b_values)
    b_max_bin = int(np.nanargmax(b_arr))
    b_min_bin = int(np.nanargmin(b_arr))

    b_max_phase = float((b_max_bin + 0.5) / k)
    b_min_phase = float((b_min_bin + 0.5) / k)
    b_max_phase_class = classify_phase(b_max_phase)
    b_min_phase_class = classify_phase(b_min_phase)

    logger.info(
        "  k=%d: b_mean=%.4f, b_std=%.4f, b_range=%.4f", k, b_mean, b_std, b_range
    )
    logger.info(
        "  k=%d: b_max at bin %d (phase=%.3f, class=%s)", k, b_max_bin, b_max_phase, b_max_phase_class
    )
    logger.info(
        "  k=%d: b_min at bin %d (phase=%.3f, class=%s)", k, b_min_bin, b_min_phase, b_min_phase_class
    )

    # One-way ANOVA across magnitude groups
    non_empty = [g for g in groups_for_anova if len(g) > 0]
    f_stat, p_anova = scipy.stats.f_oneway(*non_empty)
    logger.info("  k=%d: ANOVA F=%.4f, p=%.6f", k, f_stat, p_anova)

    # Pearson correlation between rate (count) and b-value
    counts_arr = np.array(counts, dtype=float)
    b_arr_valid_mask = ~np.isnan(b_arr)
    r_rate_b, p_rate_b = scipy.stats.pearsonr(
        counts_arr[b_arr_valid_mask], b_arr[b_arr_valid_mask]
    )
    inverse_phase_relationship = bool(r_rate_b < 0)
    logger.info(
        "  k=%d: Pearson r(rate,b)=%.4f, p=%.6f, inverse=%s",
        k, r_rate_b, p_rate_b, inverse_phase_relationship,
    )

    return {
        "bins": bins_data,
        "b_mean": b_mean,
        "b_std": b_std,
        "b_range": b_range,
        "b_max_bin": b_max_bin,
        "b_max_phase_center": b_max_phase,
        "b_max_phase_class": b_max_phase_class,
        "b_min_bin": b_min_bin,
        "b_min_phase_center": b_min_phase,
        "b_min_phase_class": b_min_phase_class,
        "f_stat": float(f_stat),
        "p_anova": float(p_anova),
        "r_rate_b": float(r_rate_b),
        "p_rate_b": float(p_rate_b),
        "inverse_phase_relationship": inverse_phase_relationship,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Execute Case A2 analysis: b-value seasonal variation across solar phase bins."""
    logger.info("=== Case A2: b-Value Seasonal Variation ===")

    # Load catalog
    df = load_catalog(RAW_PATH)

    # Compute solar phase
    df["solar_phase"] = compute_solar_phase(df["solar_secs"])
    logger.info(
        "Solar phase range: [%.6f, %.6f]",
        df["solar_phase"].min(), df["solar_phase"].max(),
    )
    assert (df["solar_phase"] >= 0.0).all() and (df["solar_phase"] < 1.0).all(), \
        "Solar phase values outside [0, 1)"

    # Assign bin indices for each k
    for k in BIN_COUNTS:
        col = f"bin_idx_k{k}"
        df[col] = np.floor(df["solar_phase"] * k).astype(int)
        # Clip to [0, k-1] in case of floating point edge cases
        df[col] = df[col].clip(0, k - 1)
        logger.info("k=%d: bin distribution min=%d, max=%d", k, df[col].min(), df[col].max())

    # Analyze each bin count
    phase_variation: dict[str, Any] = {}
    for k in BIN_COUNTS:
        col = f"bin_idx_k{k}"
        key = f"k{k}"
        result = analyze_bins(df, k, col)
        phase_variation[key] = result

    # Assemble results
    results: dict[str, Any] = {
        "case": "A2",
        "title": "b-Value Seasonal Variation",
        "julian_year_secs": JULIAN_YEAR_SECS,
        "mc": MC,
        "n_catalog": EXPECTED_N,
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
        "random_seed": RANDOM_SEED,
        "bin_counts": BIN_COUNTS,
        "phase_variation": phase_variation,
    }

    # Write output
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", OUTPUT_PATH)
    logger.info("=== Case A2 analysis complete ===")


if __name__ == "__main__":
    main()
