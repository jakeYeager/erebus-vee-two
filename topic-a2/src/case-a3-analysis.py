"""Case A3: Magnitude Stratification of the Solar Signal — Main Analysis Script.

Stratifies the ISC-GEM catalog into four magnitude bands and independently
computes solar-phase statistics (chi-square, Cramér's V, Rayleigh R) for each
band. Performs effect-size trend analysis across bands and evaluates competing
mechanistic predictions. Writes results to output/case-a3-results.json.
"""

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("case-a3-analysis")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SOLAR_YEAR_SECS = 31_557_600.0  # Julian year constant (confirmed: use uniformly)
EXPECTED_TOTAL = 9210
BIN_COUNTS = [16, 24, 32]
BOOTSTRAP_N = 1000
BOOTSTRAP_SEED = 42

BANDS = [
    {"label": "M6.0-6.4", "min": 6.0, "max": 6.5},
    {"label": "M6.5-6.9", "min": 6.5, "max": 7.0},
    {"label": "M7.0-7.4", "min": 7.0, "max": 7.5},
    {"label": "M7.5+",    "min": 7.5, "max": 99.0},
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_catalog(path: Path) -> pd.DataFrame:
    """Load the ISC-GEM catalog and assert row count.

    Args:
        path: Path to the raw CSV file.

    Returns:
        Loaded DataFrame with 9,210 rows.
    """
    logger.info("Loading catalog from %s", path)
    df = pd.read_csv(path)
    n = len(df)
    logger.info("Loaded %d rows", n)
    assert n == EXPECTED_TOTAL, f"Expected {EXPECTED_TOTAL} rows, got {n}"
    return df


def split_by_magnitude(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Split catalog into four magnitude bands.

    Args:
        df: Full catalog DataFrame with 'usgs_mag' column.

    Returns:
        Dict mapping band label to band DataFrame.
    """
    band_dfs: dict[str, pd.DataFrame] = {}
    for band in BANDS:
        mask = (df["usgs_mag"] >= band["min"]) & (df["usgs_mag"] < band["max"])
        df_band = df[mask].copy()
        n_band = len(df_band)
        logger.info("Band %s: %d events (mag in [%.1f, %.1f))",
                    band["label"], n_band, band["min"], band["max"])
        band_dfs[band["label"]] = df_band

    total = sum(len(v) for v in band_dfs.values())
    assert total == EXPECTED_TOTAL, (
        f"Band sizes sum to {total}, expected {EXPECTED_TOTAL}. "
        "Check for missing or overlapping band boundaries."
    )
    logger.info("Band partition check passed: total = %d", total)
    return band_dfs


# ---------------------------------------------------------------------------
# Phase computation
# ---------------------------------------------------------------------------
def compute_phase(solar_secs: pd.Series) -> np.ndarray:
    """Compute solar year phase using Julian year constant.

    Args:
        solar_secs: Series of solar seconds within the year.

    Returns:
        Phase array in [0, 1).
    """
    return ((solar_secs / SOLAR_YEAR_SECS) % 1.0).to_numpy()


# ---------------------------------------------------------------------------
# Rayleigh test
# ---------------------------------------------------------------------------
def compute_rayleigh(phases: np.ndarray) -> tuple[float, float, float]:
    """Compute Rayleigh R statistic, p-value, and mean phase angle fraction.

    Args:
        phases: Array of phase values in [0, 1).

    Returns:
        Tuple of (R, p_rayleigh, mean_phase_fraction).
    """
    n = len(phases)
    angles = 2.0 * np.pi * phases
    mean_cos = float(np.mean(np.cos(angles)))
    mean_sin = float(np.mean(np.sin(angles)))
    R = float(np.sqrt(mean_cos ** 2 + mean_sin ** 2))
    p_rayleigh = float(np.exp(-n * R ** 2))
    mean_angle_rad = float(np.arctan2(mean_sin, mean_cos))
    mean_phase = float((mean_angle_rad / (2.0 * np.pi)) % 1.0)
    return R, p_rayleigh, mean_phase


# ---------------------------------------------------------------------------
# Elevated-bin intervals
# ---------------------------------------------------------------------------
def find_elevated_intervals(
    bin_counts: np.ndarray, k: int, n: int
) -> list[dict[str, float]]:
    """Identify bins above the 1-SD threshold and merge into intervals.

    A bin is elevated if count > E + sqrt(E) where E = n/k.

    Args:
        bin_counts: Array of observed counts per bin.
        k: Number of bins.
        n: Total event count.

    Returns:
        List of dicts with phase_start, phase_end, mean_phase.
    """
    expected = n / k
    threshold = expected + np.sqrt(expected)
    elevated_mask = bin_counts > threshold

    bin_width = 1.0 / k
    intervals: list[dict[str, float]] = []

    i = 0
    while i < k:
        if elevated_mask[i]:
            start_bin = i
            while i < k and elevated_mask[i]:
                i += 1
            end_bin = i  # exclusive
            phase_start = start_bin * bin_width
            phase_end = end_bin * bin_width
            mean_phase = (phase_start + phase_end) / 2.0
            intervals.append({
                "phase_start": round(phase_start, 6),
                "phase_end": round(phase_end, 6),
                "mean_phase": round(mean_phase, 6),
            })
        else:
            i += 1

    return intervals


# ---------------------------------------------------------------------------
# Per-band bin statistics at one k
# ---------------------------------------------------------------------------
def compute_band_stats_at_k(phases: np.ndarray, k: int) -> dict[str, Any]:
    """Compute full bin statistics for one band at one bin count.

    Args:
        phases: Phase array in [0, 1).
        k: Number of bins.

    Returns:
        Dict with chi2, p_chi2, cramer_v, rayleigh_R, p_rayleigh,
        mean_phase, elevated_intervals, bin_counts.
    """
    n = len(phases)
    bin_indices = np.floor(phases * k).astype(int)
    bin_indices = np.clip(bin_indices, 0, k - 1)  # guard against phase=1.0
    O = np.bincount(bin_indices, minlength=k).astype(float)

    expected = n / k
    E = np.full(k, expected)

    chi2, p_chi2 = scipy.stats.chisquare(O, E)
    cramer_v = float(np.sqrt(chi2 / (n * (k - 1))))

    R, p_rayleigh, mean_phase = compute_rayleigh(phases)

    elevated_intervals = find_elevated_intervals(O, k, n)

    return {
        "chi2": float(chi2),
        "p_chi2": float(p_chi2),
        "cramer_v": float(cramer_v),
        "rayleigh_R": float(R),
        "p_rayleigh": float(p_rayleigh),
        "mean_phase": float(mean_phase),
        "bin_counts": O.astype(int).tolist(),
        "elevated_intervals": elevated_intervals,
    }


# ---------------------------------------------------------------------------
# Bootstrap Cramér's V CI
# ---------------------------------------------------------------------------
def compute_bootstrap_ci(
    phases: np.ndarray, k: int, rng: np.random.Generator
) -> tuple[float, float]:
    """Compute 95% bootstrap CI for Cramér's V at k bins.

    Uses 1000 bootstrap resamples drawn from the provided RNG.

    Args:
        phases: Phase array in [0, 1).
        k: Number of bins.
        rng: NumPy random generator (seeded externally for reproducibility).

    Returns:
        Tuple of (ci95_lower, ci95_upper).
    """
    n = len(phases)
    v_vals: list[float] = []

    for _ in range(BOOTSTRAP_N):
        resample_idx = rng.integers(0, n, size=n)
        resampled = phases[resample_idx]

        bin_indices = np.floor(resampled * k).astype(int)
        bin_indices = np.clip(bin_indices, 0, k - 1)
        O_boot = np.bincount(bin_indices, minlength=k).astype(float)
        E_boot = np.full(k, n / k)

        chi2_boot, _ = scipy.stats.chisquare(O_boot, E_boot)
        v_boot = float(np.sqrt(chi2_boot / (n * (k - 1))))
        v_vals.append(v_boot)

    lower = float(np.percentile(v_vals, 2.5))
    upper = float(np.percentile(v_vals, 97.5))
    return lower, upper


# ---------------------------------------------------------------------------
# Full per-band analysis
# ---------------------------------------------------------------------------
def compute_all_band_stats(
    band_dfs: dict[str, pd.DataFrame]
) -> dict[str, Any]:
    """Compute chi-square, Rayleigh, Cramér's V, and bootstrap CIs for all bands.

    A single RNG seeded with BOOTSTRAP_SEED is used for all bootstrap resamples
    across all bands, ensuring full reproducibility.

    Args:
        band_dfs: Dict mapping band label to band DataFrame.

    Returns:
        Dict structured as band_stats JSON output.
    """
    band_stats: dict[str, Any] = {}

    # Single RNG seeded once for all bootstrap resamples
    rng = np.random.default_rng(BOOTSTRAP_SEED)

    for band in BANDS:
        label = band["label"]
        df_band = band_dfs[label]
        n = len(df_band)
        phases = compute_phase(df_band["solar_secs"])

        logger.info("Band %s: n=%d, phase range [%.6f, %.6f]",
                    label, n, phases.min(), phases.max())

        entry: dict[str, Any] = {"n": n}

        for k in BIN_COUNTS:
            k_key = f"k{k}"
            stats = compute_band_stats_at_k(phases, k)

            # Bootstrap CI only at k=24
            if k == 24:
                ci_lower, ci_upper = compute_bootstrap_ci(phases, k, rng)
                stats["cramer_v_ci95_lower"] = ci_lower
                stats["cramer_v_ci95_upper"] = ci_upper
                logger.info(
                    "Band %s k=%d: chi2=%.3f p=%.4e V=%.4f CI95=[%.4f,%.4f]",
                    label, k,
                    stats["chi2"], stats["p_chi2"], stats["cramer_v"],
                    ci_lower, ci_upper,
                )
            else:
                logger.info(
                    "Band %s k=%d: chi2=%.3f p=%.4e V=%.4f",
                    label, k, stats["chi2"], stats["p_chi2"], stats["cramer_v"],
                )

            entry[k_key] = stats

        band_stats[label] = entry

    return band_stats


# ---------------------------------------------------------------------------
# Effect-size trend analysis
# ---------------------------------------------------------------------------
def compute_trend_analysis(band_stats: dict[str, Any]) -> dict[str, Any]:
    """Compute Spearman rank correlations and classify the effect-size trend.

    Uses k=24 as primary. Evaluates four mechanistic predictions.

    Args:
        band_stats: Per-band stats dict from compute_all_band_stats.

    Returns:
        Trend analysis dict matching the spec JSON schema.
    """
    band_labels = [b["label"] for b in BANDS]  # ascending magnitude order
    band_indices = [1, 2, 3, 4]

    cramer_v_by_band: dict[str, float] = {}
    rayleigh_R_by_band: dict[str, float] = {}

    for label in band_labels:
        cramer_v_by_band[label] = band_stats[label]["k24"]["cramer_v"]
        rayleigh_R_by_band[label] = band_stats[label]["k24"]["rayleigh_R"]

    v_vals = [cramer_v_by_band[l] for l in band_labels]
    r_vals = [rayleigh_R_by_band[l] for l in band_labels]

    rho_v, p_v = scipy.stats.spearmanr(band_indices, v_vals)
    rho_r, p_r = scipy.stats.spearmanr(band_indices, r_vals)

    # Trend classification based on Cramér's V Spearman rho
    if rho_v > 0.5:
        trend_classification = "increasing"
    elif rho_v < -0.5:
        trend_classification = "decreasing"
    else:
        trend_classification = "flat"

    logger.info(
        "Trend analysis: Cramér's V Spearman rho=%.4f p=%.4e → %s",
        rho_v, p_v, trend_classification,
    )
    logger.info(
        "Trend analysis: Rayleigh R Spearman rho=%.4f p=%.4e",
        rho_r, p_r,
    )

    # Significant bands (p_chi2 < 0.05 at k=24)
    significant_bands = [
        label for label in band_labels
        if band_stats[label]["k24"]["p_chi2"] < 0.05
    ]
    m75_significant = "M7.5+" in significant_bands

    logger.info("Significant bands at k=24: %s", significant_bands)
    logger.info("M7.5+ significant: %s", m75_significant)

    # Prediction matching
    # 1. Hydrological loading: weakens/disappears at M7.5+ → M7.5+ not significant AND trend decreasing
    if not m75_significant and trend_classification == "decreasing":
        hydro = "supported"
    elif not m75_significant or trend_classification == "decreasing":
        hydro = "partially"
    else:
        hydro = "not supported"

    # 2. Solar-geometric forcing: flat or increases with magnitude → trend flat or increasing AND M7.5+ significant
    if (trend_classification in ("flat", "increasing")) and m75_significant:
        solar_geo = "supported"
    elif trend_classification in ("flat", "increasing") or m75_significant:
        solar_geo = "partially"
    else:
        solar_geo = "not supported"

    # 3. Tidal literature pattern: M6.0-6.4 shows strongest signal; higher bands weaker or absent
    # Evaluate: is M6.0-6.4 the band with the highest Cramér's V?
    max_band = max(band_labels, key=lambda l: cramer_v_by_band[l])
    if max_band == "M6.0-6.4":
        tidal = "supported"
    elif v_vals[0] > v_vals[-1]:  # at least M6.0-6.4 > M7.5+
        tidal = "partially"
    else:
        tidal = "not supported"

    # 4. Magnitude-independent (A1b preliminary): max V / min V < 1.5
    v_max = max(v_vals)
    v_min = min(v_vals)
    ratio = v_max / v_min if v_min > 0 else float("inf")
    if ratio < 1.5:
        mag_indep = "supported"
    elif ratio < 2.0:
        mag_indep = "partially"
    else:
        mag_indep = "not supported"

    logger.info("Prediction support:")
    logger.info("  hydrological_loading: %s", hydro)
    logger.info("  solar_geometric: %s", solar_geo)
    logger.info("  tidal_literature_pattern: %s", tidal)
    logger.info("  magnitude_independent: %s", mag_indep)

    return {
        "cramer_v_by_band": cramer_v_by_band,
        "rayleigh_R_by_band": rayleigh_R_by_band,
        "spearman_rho_cramer_v": float(rho_v),
        "spearman_p_cramer_v": float(p_v),
        "spearman_rho_rayleigh": float(rho_r),
        "spearman_p_rayleigh": float(p_r),
        "trend_classification": trend_classification,
        "significant_bands": significant_bands,
        "m75_significant": bool(m75_significant),
        "cramer_v_ratio_max_min": float(ratio),
        "prediction_support": {
            "hydrological_loading": hydro,
            "solar_geometric": solar_geo,
            "tidal_literature_pattern": tidal,
            "magnitude_independent": mag_indep,
        },
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load catalog, compute per-band stats, run trend analysis, write JSON."""
    # --- Load and split ---
    df = load_catalog(RAW_PATH)
    band_dfs = split_by_magnitude(df)

    # --- Per-band statistics ---
    logger.info("Computing per-band statistics at k=16, 24, 32")
    band_stats = compute_all_band_stats(band_dfs)

    # --- Trend analysis ---
    logger.info("Computing effect-size trend analysis")
    trend_analysis = compute_trend_analysis(band_stats)

    # --- Assemble results ---
    results: dict[str, Any] = {
        "case": "A3",
        "title": "Magnitude Stratification of the Solar Signal",
        "solar_year_secs": SOLAR_YEAR_SECS,
        "band_stats": band_stats,
        "trend_analysis": trend_analysis,
    }

    # --- Write JSON ---
    output_path = OUTPUT_DIR / "case-a3-results.json"
    with open(output_path, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", output_path)


if __name__ == "__main__":
    main()
