"""
Case B5: Solar Declination Rate-of-Change vs. Position Test

Tests whether solar declination, declination rate of change, or Earth-Sun distance
better explain the seismic clustering signal found in Case 3A's solar_secs variable.
Computes bin distributions, chi-square, Cramér's V, Rayleigh/KS tests, and variable
ranking. Performs A1b interval alignment analysis.
"""

import json
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from typing import Optional

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
SOLAR_GEO_PATH = BASE_DIR.parent / "data" / "iscgem" / "celestial-geometry" / "solar_geometry_global.csv"
OUTPUT_PATH = BASE_DIR / "output" / "case-b5-results.json"

# ── Logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────
JULIAN_YEAR_SECS: float = 31_557_600.0
BIN_COUNTS: list[int] = [16, 24, 32]

# A1b interval phase centers (from spec)
A1B_INTERVAL_PHASES: dict[str, float] = {
    "interval_1_phase_0.22": 0.22,
    "interval_2_phase_0.64": 0.64,
    "interval_3_phase_0.90": 0.90,
}

# Expected physical values at each A1b phase center (theoretical, from spec §3)
A1B_EXPECTED_VALUES: dict[str, dict] = {
    "interval_1_phase_0.22": {
        "solar_declination_expected": 0.0,
        "declination_rate_expected": 0.40,
        "earth_sun_distance_expected": 1.00,
    },
    "interval_2_phase_0.64": {
        "solar_declination_expected": 12.0,
        "declination_rate_expected": -0.25,
        "earth_sun_distance_expected": 1.012,
    },
    "interval_3_phase_0.90": {
        "solar_declination_expected": -21.0,
        "declination_rate_expected": -0.12,
        "earth_sun_distance_expected": 0.988,
    },
}


# ── Helper functions ───────────────────────────────────────────────────────────

def compute_cramer_v(chi2: float, n: int, k: int) -> float:
    """Compute Cramér's V from chi-square statistic.

    Parameters
    ----------
    chi2:
        Chi-square test statistic.
    n:
        Sample size.
    k:
        Number of bins (degrees of freedom = k - 1).

    Returns
    -------
    float
        Cramér's V in [0.0, 1.0].
    """
    return math.sqrt(chi2 / (n * (k - 1)))


def compute_rayleigh(phases: np.ndarray) -> tuple[float, float]:
    """Compute Rayleigh test R-statistic and p-value for cyclic data.

    Parameters
    ----------
    phases:
        Array of phase values in [0, 1).

    Returns
    -------
    tuple[float, float]
        (R, p_value) where R is the mean resultant length.
    """
    angles = 2.0 * np.pi * phases
    n = len(angles)
    c = np.cos(angles).sum()
    s = np.sin(angles).sum()
    r_bar = math.sqrt(c**2 + s**2) / n
    # Rayleigh p-value approximation (Mardia & Jupp 2000)
    z = n * r_bar**2
    p = math.exp(-z) * (1.0 + (2.0 * z - z**2) / (4.0 * n) - (24.0 * z - 132.0 * z**2 + 76.0 * z**3 - 9.0 * z**4) / (288.0 * n**2))
    p = max(0.0, min(1.0, p))
    return float(r_bar), float(p)


def compute_chi_square_uniform(
    bin_counts: np.ndarray, n: int, k: int
) -> tuple[float, float]:
    """Chi-square test against a uniform distribution.

    Parameters
    ----------
    bin_counts:
        Observed counts per bin.
    n:
        Total number of observations.
    k:
        Number of bins.

    Returns
    -------
    tuple[float, float]
        (chi2_stat, p_value)
    """
    expected = np.full(k, n / k)
    chi2_stat = float(np.sum((bin_counts - expected) ** 2 / expected))
    p_value = float(stats.chi2.sf(chi2_stat, df=k - 1))
    return chi2_stat, p_value


def compute_bins(
    values: np.ndarray, actual_min: float, actual_range: float, k: int
) -> np.ndarray:
    """Assign non-cyclic bin indices.

    Parameters
    ----------
    values:
        Raw variable values.
    actual_min:
        Minimum of the variable in the dataset.
    actual_range:
        Range (max - min) of the variable in the dataset.
    k:
        Number of bins.

    Returns
    -------
    np.ndarray
        Integer bin indices in [0, k-1].
    """
    raw = np.floor((values - actual_min) / actual_range * k).astype(int)
    return np.clip(raw, 0, k - 1)


def variable_stats_block(
    bin_counts: np.ndarray,
    n: int,
    k: int,
    is_cyclic: bool,
    phases: Optional[np.ndarray] = None,
    normalized_values: Optional[np.ndarray] = None,
) -> dict:
    """Compute all statistics for one variable at one k.

    Parameters
    ----------
    bin_counts:
        Observed bin counts (length k).
    n:
        Total observations.
    k:
        Number of bins.
    is_cyclic:
        True for solar_phase (Rayleigh test); False for others (KS test).
    phases:
        Phase values in [0, 1) — required when is_cyclic=True.
    normalized_values:
        Values normalized to [0, 1) — required when is_cyclic=False.

    Returns
    -------
    dict
        Statistics block ready for JSON serialization.
    """
    chi2, p_chi2 = compute_chi_square_uniform(bin_counts, n, k)
    cramer_v = compute_cramer_v(chi2, n, k)

    block: dict = {
        "chi2": chi2,
        "p_chi2": p_chi2,
        "cramer_v": cramer_v,
        "bin_counts": bin_counts.tolist(),
    }

    if is_cyclic:
        assert phases is not None, "phases required for cyclic variable"
        r_bar, p_rayleigh = compute_rayleigh(phases)
        block["rayleigh_R"] = r_bar
        block["p_rayleigh"] = p_rayleigh
    else:
        assert normalized_values is not None, "normalized_values required for non-cyclic"
        ks_stat, ks_p = stats.kstest(normalized_values, "uniform")
        block["ks_stat"] = float(ks_stat)
        block["ks_p"] = float(ks_p)

    return block


# ── Main analysis ──────────────────────────────────────────────────────────────

def main() -> None:
    """Execute Case B5 analysis."""

    # ── 1. Load and validate data ──────────────────────────────────────────────
    log.info("Loading solar geometry catalog: %s", SOLAR_GEO_PATH)
    df = pd.read_csv(SOLAR_GEO_PATH)
    n = len(df)
    log.info("Loaded %d rows", n)
    assert n == 9210, f"Expected 9210 rows, got {n}"

    required_cols = ["solar_secs", "solar_declination", "declination_rate", "earth_sun_distance"]
    for col in required_cols:
        assert col in df.columns, f"Missing required column: {col}"
    log.info("All required columns present")

    # Actual ranges
    dec_min, dec_max = float(df["solar_declination"].min()), float(df["solar_declination"].max())
    rate_min, rate_max = float(df["declination_rate"].min()), float(df["declination_rate"].max())
    dist_min, dist_max = float(df["earth_sun_distance"].min()), float(df["earth_sun_distance"].max())

    log.info("solar_declination actual range: %.6f to %.6f (spec: -23.5 to +23.5)", dec_min, dec_max)
    log.info("declination_rate actual range: %.6f to %.6f (spec: -0.40 to +0.40)", rate_min, rate_max)
    log.info("earth_sun_distance actual range: %.6f to %.6f (spec: 0.983 to 1.017)", dist_min, dist_max)

    # Warn if outside spec bounds
    if dec_min < -23.5 or dec_max > 23.5:
        log.warning("solar_declination outside expected bounds [-23.5, +23.5]")
    if rate_min < -0.40 or rate_max > 0.40:
        log.warning("declination_rate outside expected bounds [-0.40, +0.40]")
    if dist_min < 0.983 or dist_max > 1.017:
        log.warning("earth_sun_distance outside expected bounds [0.983, 1.017]")

    log.info("Using actual min/max for bin computation")

    dec_range = dec_max - dec_min
    rate_range = rate_max - rate_min
    dist_range = dist_max - dist_min

    # ── 2. Compute phases and bins ─────────────────────────────────────────────
    solar_phase = (df["solar_secs"].values / JULIAN_YEAR_SECS) % 1.0
    solar_declination = df["solar_declination"].values
    declination_rate = df["declination_rate"].values
    earth_sun_distance = df["earth_sun_distance"].values

    # Normalized values (0..1) for KS test
    norm_dec = (solar_declination - dec_min) / dec_range
    norm_rate = (declination_rate - rate_min) / rate_range
    norm_dist = (earth_sun_distance - dist_min) / dist_range

    log.info("solar_phase range: %.6f to %.6f", solar_phase.min(), solar_phase.max())

    # ── 3. Build stats for all variables at k=16, 24, 32 ──────────────────────
    variable_stats: dict = {
        "solar_phase": {},
        "solar_declination": {},
        "declination_rate": {},
        "earth_sun_distance": {},
    }

    for k in BIN_COUNTS:
        key = f"k{k}"

        # solar_phase (cyclic)
        bins_sp = (np.floor(solar_phase * k)).astype(int)
        bins_sp = np.clip(bins_sp, 0, k - 1)
        counts_sp = np.bincount(bins_sp, minlength=k)
        variable_stats["solar_phase"][key] = variable_stats_block(
            counts_sp, n, k, is_cyclic=True, phases=solar_phase
        )

        # solar_declination (non-cyclic)
        bins_dec = compute_bins(solar_declination, dec_min, dec_range, k)
        counts_dec = np.bincount(bins_dec, minlength=k)
        variable_stats["solar_declination"][key] = variable_stats_block(
            counts_dec, n, k, is_cyclic=False, normalized_values=norm_dec
        )

        # declination_rate (non-cyclic)
        bins_rate = compute_bins(declination_rate, rate_min, rate_range, k)
        counts_rate = np.bincount(bins_rate, minlength=k)
        variable_stats["declination_rate"][key] = variable_stats_block(
            counts_rate, n, k, is_cyclic=False, normalized_values=norm_rate
        )

        # earth_sun_distance (non-cyclic)
        bins_dist = compute_bins(earth_sun_distance, dist_min, dist_range, k)
        counts_dist = np.bincount(bins_dist, minlength=k)
        variable_stats["earth_sun_distance"][key] = variable_stats_block(
            counts_dist, n, k, is_cyclic=False, normalized_values=norm_dist
        )

        log.info(
            "k=%2d | solar_phase chi2=%.2f p=%.4f V=%.4f | "
            "declination chi2=%.2f p=%.4f | rate chi2=%.2f p=%.4f | dist chi2=%.2f p=%.4f",
            k,
            variable_stats["solar_phase"][key]["chi2"],
            variable_stats["solar_phase"][key]["p_chi2"],
            variable_stats["solar_phase"][key]["cramer_v"],
            variable_stats["solar_declination"][key]["chi2"],
            variable_stats["solar_declination"][key]["p_chi2"],
            variable_stats["declination_rate"][key]["chi2"],
            variable_stats["declination_rate"][key]["p_chi2"],
            variable_stats["earth_sun_distance"][key]["chi2"],
            variable_stats["earth_sun_distance"][key]["p_chi2"],
        )

    # ── 4. Variable ranking at k=24 ────────────────────────────────────────────
    ranking_data = [
        {
            "variable": var,
            "cramer_v_k24": variable_stats[var]["k24"]["cramer_v"],
            "p_chi2_k24": variable_stats[var]["k24"]["p_chi2"],
        }
        for var in variable_stats
    ]
    ranking_data.sort(key=lambda x: x["cramer_v_k24"], reverse=True)

    most_significant_variable = min(
        variable_stats.keys(),
        key=lambda v: variable_stats[v]["k24"]["p_chi2"]
    )

    log.info("Variable ranking (k=24, by Cramér's V):")
    for rank in ranking_data:
        log.info("  %s: V=%.4f, p=%.4e", rank["variable"], rank["cramer_v_k24"], rank["p_chi2_k24"])
    log.info("Most significant variable: %s", most_significant_variable)

    # ── 5. A1b interval alignment analysis ────────────────────────────────────
    # For each A1b interval, determine which physical variable has an elevated bin
    # at the expected physical value for that interval's phase center.
    # "Elevated" = bin count > mean + 1 SD

    def get_bin_stats_k24(counts: np.ndarray) -> tuple[float, float]:
        """Return (mean, std) of k=24 bin counts."""
        return float(counts.mean()), float(counts.std())

    # Recompute k=24 bins for alignment analysis
    bins_sp_24 = np.clip((np.floor(solar_phase * 24)).astype(int), 0, 23)
    counts_sp_24 = np.bincount(bins_sp_24, minlength=24)
    mean_sp, std_sp = get_bin_stats_k24(counts_sp_24)

    bins_dec_24 = compute_bins(solar_declination, dec_min, dec_range, 24)
    counts_dec_24 = np.bincount(bins_dec_24, minlength=24)
    mean_dec, std_dec = get_bin_stats_k24(counts_dec_24)

    bins_rate_24 = compute_bins(declination_rate, rate_min, rate_range, 24)
    counts_rate_24 = np.bincount(bins_rate_24, minlength=24)
    mean_rate, std_rate = get_bin_stats_k24(counts_rate_24)

    bins_dist_24 = compute_bins(earth_sun_distance, dist_min, dist_range, 24)
    counts_dist_24 = np.bincount(bins_dist_24, minlength=24)
    mean_dist, std_dist = get_bin_stats_k24(counts_dist_24)

    def is_elevated(value: float, val_min: float, val_range: float,
                    counts: np.ndarray, mean: float, std: float) -> bool:
        """Check whether the bin containing `value` has an elevated count."""
        bin_idx = int(np.clip(math.floor((value - val_min) / val_range * 24), 0, 23))
        return counts[bin_idx] > mean + std

    # For solar_phase, convert phase to bin directly
    def is_elevated_cyclic(phase: float, counts: np.ndarray, mean: float, std: float) -> bool:
        bin_idx = int(np.clip(math.floor(phase * 24), 0, 23))
        return counts[bin_idx] > mean + std

    a1b_alignment: dict = {}

    interval_configs = [
        ("interval_1_phase_0.22", 0.22),
        ("interval_2_phase_0.64", 0.64),
        ("interval_3_phase_0.90", 0.90),
    ]

    for interval_key, phase_center in interval_configs:
        expected = A1B_EXPECTED_VALUES[interval_key]
        elevated: list[str] = []

        # solar_phase (cyclic)
        if is_elevated_cyclic(phase_center, counts_sp_24, mean_sp, std_sp):
            elevated.append("solar_phase")

        # solar_declination
        exp_dec = expected["solar_declination_expected"]
        if is_elevated(exp_dec, dec_min, dec_range, counts_dec_24, mean_dec, std_dec):
            elevated.append("solar_declination")

        # declination_rate
        exp_rate = expected["declination_rate_expected"]
        if is_elevated(exp_rate, rate_min, rate_range, counts_rate_24, mean_rate, std_rate):
            elevated.append("declination_rate")

        # earth_sun_distance
        exp_dist = expected["earth_sun_distance_expected"]
        if is_elevated(exp_dist, dist_min, dist_range, counts_dist_24, mean_dist, std_dist):
            elevated.append("earth_sun_distance")

        a1b_alignment[interval_key] = {
            **expected,
            "variable_with_elevated_bin": elevated,
        }

        log.info(
            "A1b %s (phase=%.2f): elevated variables = %s",
            interval_key, phase_center, elevated
        )

    # ── 6. Collect data ranges for JSON ───────────────────────────────────────
    actual_ranges: dict = {
        "solar_declination": {"min": dec_min, "max": dec_max, "range": dec_range,
                               "spec_min": -23.5, "spec_max": 23.5, "used": "actual"},
        "declination_rate": {"min": rate_min, "max": rate_max, "range": rate_range,
                              "spec_min": -0.40, "spec_max": 0.40, "used": "actual"},
        "earth_sun_distance": {"min": dist_min, "max": dist_max, "range": dist_range,
                                "spec_min": 0.983, "spec_max": 1.017, "used": "actual"},
    }

    # ── 7. Assemble results JSON ───────────────────────────────────────────────
    results = {
        "case": "B5",
        "title": "Solar Declination Rate-of-Change vs. Position Test",
        "n": n,
        "julian_year_secs": JULIAN_YEAR_SECS,
        "actual_ranges": actual_ranges,
        "variable_stats": variable_stats,
        "variable_ranking": ranking_data,
        "most_significant_variable": most_significant_variable,
        "a1b_alignment": a1b_alignment,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as fh:
        json.dump(results, fh, indent=2)

    log.info("Results written to %s", OUTPUT_PATH)


if __name__ == "__main__":
    main()
