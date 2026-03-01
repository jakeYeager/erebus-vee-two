"""Case B4: Depth Stratification — Surface Loading Penetration Test — Main Analysis Script.

Stratifies the ISC-GEM catalog by focal depth into four bands:
  - shallow       (0–20 km)
  - mid-crustal   (20–70 km)
  - intermediate  (70–300 km)
  - deep          (≥300 km, included if n ≥ 100)

Computes solar-phase bin statistics (chi-square, Cramér's V, Rayleigh R) for
each band at k=16, 24, 32. Performs Spearman monotonicity trend analysis.
Evaluates three mechanistic predictions (surface loading, geometric/deep forcing,
Zhan & Shearer 2015 pattern). Writes results to output/case-b4-results.json.
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
logger = logging.getLogger("case-b4-analysis")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SOLAR_YEAR_SECS = 31_557_600.0  # Julian year constant (confirmed: use uniformly)
EXPECTED_TOTAL = 9210
BIN_COUNTS = [16, 24, 32]
DEEP_BAND_MIN_N = 100  # minimum n for deep band to be included meaningfully

DEPTH_BANDS = [
    {"label": "shallow_0-20km",         "min": 0,   "max": 20},
    {"label": "midcrustal_20-70km",     "min": 20,  "max": 70},
    {"label": "intermediate_70-300km",  "min": 70,  "max": 300},
    {"label": "deep_300km+",            "min": 300, "max": 9999},
]

# A1b baseline elevated phase intervals (from Adhoc A1b analysis)
A1B_INTERVALS = [
    (0.1875, 0.25),   # Interval 1: ~March equinox
    (0.625, 0.656),   # Interval 2: ~mid-August
    (0.875, 0.917),   # Interval 3: ~mid-November
]

# April–September phase range (Zhan & Shearer 2015 comparison; Jan 1 = phase 0.0)
ZHAN_SHEARER_PHASE_MIN = 0.23  # ~late March
ZHAN_SHEARER_PHASE_MAX = 0.67  # ~early September


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


# ---------------------------------------------------------------------------
# Depth band partitioning
# ---------------------------------------------------------------------------
def split_by_depth(df: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], int]:
    """Split catalog into four depth bands, handling NaN depths.

    Events with NaN depth are excluded from band analysis. The sum of all
    band sizes plus n_depth_null must equal EXPECTED_TOTAL.

    Args:
        df: Full catalog DataFrame with 'depth' column.

    Returns:
        Tuple of (dict mapping band label to DataFrame, n_depth_null).
    """
    # Check for NaN depths
    n_depth_null = int(df["depth"].isna().sum())
    logger.info("NaN depth events: %d", n_depth_null)
    df_valid = df[df["depth"].notna()].copy()

    band_dfs: dict[str, pd.DataFrame] = {}
    for band in DEPTH_BANDS:
        label = band["label"]
        if band["max"] == 9999:
            mask = df_valid["depth"] >= band["min"]
        else:
            mask = (df_valid["depth"] >= band["min"]) & (df_valid["depth"] < band["max"])
        df_band = df_valid[mask].copy()
        n_band = len(df_band)
        logger.info(
            "Band %s: %d events (depth in [%.0f, %.0f))",
            label, n_band, band["min"], band["max"],
        )
        band_dfs[label] = df_band

    total_bands = sum(len(v) for v in band_dfs.values())
    expected_valid = EXPECTED_TOTAL - n_depth_null
    assert total_bands == expected_valid, (
        f"Band sizes sum to {total_bands}, expected {expected_valid} "
        "(total minus null depths). Check band boundaries."
    )

    # Confirm grand total
    grand_total = total_bands + n_depth_null
    assert grand_total == EXPECTED_TOTAL, (
        f"Grand total {grand_total} != {EXPECTED_TOTAL}"
    )
    logger.info("Depth band partition check passed: %d + %d null = %d",
                total_bands, n_depth_null, grand_total)

    # Warn if deep band has insufficient sample
    deep_n = len(band_dfs["deep_300km+"])
    if deep_n < DEEP_BAND_MIN_N:
        logger.warning(
            "Deep band (>300 km) has n=%d, below threshold %d; "
            "results flagged as insufficient",
            deep_n, DEEP_BAND_MIN_N,
        )
    else:
        logger.info("Deep band n=%d (≥ threshold %d); included in analysis",
                    deep_n, DEEP_BAND_MIN_N)

    return band_dfs, n_depth_null


# ---------------------------------------------------------------------------
# Phase computation
# ---------------------------------------------------------------------------
def compute_phase(solar_secs: pd.Series) -> np.ndarray:
    """Compute solar year phase using Julian year constant.

    Phase normalization: phase = (solar_secs / SOLAR_YEAR_SECS) % 1.0
    Consistent with data-handling.md standard and prior cases (A1, A3, B1, B2).

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
    """Identify bins above the 1-SD threshold and merge into contiguous intervals.

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
    """Compute full bin statistics for one depth band at one bin count.

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
# Full per-band analysis
# ---------------------------------------------------------------------------
def compute_all_band_stats(
    band_dfs: dict[str, pd.DataFrame], n_depth_null: int
) -> dict[str, Any]:
    """Compute chi-square, Rayleigh, and Cramér's V for all depth bands at k=16, 24, 32.

    Args:
        band_dfs: Dict mapping band label to band DataFrame.
        n_depth_null: Count of events excluded due to NaN depth.

    Returns:
        Dict structured as band_stats JSON output including n_depth_null key.
    """
    band_stats: dict[str, Any] = {}

    for band in DEPTH_BANDS:
        label = band["label"]
        df_band = band_dfs[label]
        n = len(df_band)
        sufficient_n = n >= DEEP_BAND_MIN_N if label == "deep_300km+" else True

        phases = compute_phase(df_band["solar_secs"])

        logger.info(
            "Band %s: n=%d, sufficient_n=%s, phase range [%.6f, %.6f]",
            label, n, sufficient_n, float(phases.min()), float(phases.max()),
        )

        entry: dict[str, Any] = {
            "n": n,
            "sufficient_n": sufficient_n,
        }

        for k in BIN_COUNTS:
            k_key = f"k{k}"
            stats = compute_band_stats_at_k(phases, k)
            logger.info(
                "  k=%d: chi2=%.3f p=%.4e V=%.4f R=%.4f mean_phase=%.4f",
                k, stats["chi2"], stats["p_chi2"], stats["cramer_v"],
                stats["rayleigh_R"], stats["mean_phase"],
            )
            entry[k_key] = stats

        band_stats[label] = entry

    return band_stats


# ---------------------------------------------------------------------------
# Depth trend analysis
# ---------------------------------------------------------------------------
def compute_trend_analysis(
    band_stats: dict[str, Any]
) -> dict[str, Any]:
    """Compute Spearman monotonicity analysis and evaluate mechanistic predictions.

    Uses k=24 as primary. Evaluates three competing hypotheses:
      1. Surface loading: signal strongest shallow, absent at >70 km
      2. Geometric/deep forcing: signal present at all depths
      3. Zhan & Shearer (2015): deep events show April–September seasonality

    Args:
        band_stats: Per-band stats dict from compute_all_band_stats.

    Returns:
        Trend analysis dict matching the spec JSON schema.
    """
    band_labels = [b["label"] for b in DEPTH_BANDS]
    band_indices = [1, 2, 3, 4]

    # Extract Cramér's V and Rayleigh R at k=24 for all bands
    cramer_v_by_band: dict[str, float] = {}
    for label in band_labels:
        cramer_v_by_band[label] = band_stats[label]["k24"]["cramer_v"]

    # Deep band mean phase (for Zhan & Shearer comparison)
    deep_label = "deep_300km+"
    deep_mean_phase = float(band_stats[deep_label]["k24"]["mean_phase"])
    deep_mean_phase_in_apr_sep = bool(
        ZHAN_SHEARER_PHASE_MIN <= deep_mean_phase <= ZHAN_SHEARER_PHASE_MAX
    )

    # Spearman rank correlation of Cramér's V vs depth band index
    v_vals = [cramer_v_by_band[label] for label in band_labels]
    rho, p_spearman = scipy.stats.spearmanr(band_indices, v_vals)
    rho = float(rho)
    p_spearman = float(p_spearman)

    # Monotonicity classification
    if rho < -0.5:
        monotonicity = "decreasing with depth"
    elif rho > 0.5:
        monotonicity = "increasing with depth"
    else:
        monotonicity = "non-monotonic"

    logger.info(
        "Spearman rho=%.4f p=%.4e → monotonicity: %s",
        rho, p_spearman, monotonicity,
    )

    # Significant bands (p_chi2 < 0.05 at k=24)
    significant_bands = [
        label for label in band_labels
        if band_stats[label]["k24"]["p_chi2"] < 0.05
    ]
    deep_band_significant = deep_label in significant_bands

    logger.info("Significant bands at k=24 (p<0.05): %s", significant_bands)
    logger.info("Deep band (>300 km) significant: %s", deep_band_significant)
    logger.info(
        "Deep mean phase=%.4f, in April–September range [%.2f, %.2f]: %s",
        deep_mean_phase, ZHAN_SHEARER_PHASE_MIN, ZHAN_SHEARER_PHASE_MAX,
        deep_mean_phase_in_apr_sep,
    )

    # ---------------------------------------------------------------------------
    # Prediction matching
    # ---------------------------------------------------------------------------

    # Prediction 1 — Surface loading hypothesis:
    # Signal strongest at 0–20 km and absent at >70 km
    # Evaluate: shallow significant AND intermediate not significant
    shallow_sig = "shallow_0-20km" in significant_bands
    intermediate_sig = "intermediate_70-300km" in significant_bands

    if shallow_sig and not intermediate_sig:
        surface_loading = "supported"
    elif shallow_sig and intermediate_sig:
        surface_loading = "partially"
    else:
        surface_loading = "not supported"

    # Prediction 2 — Geometric/deep forcing hypothesis:
    # Signal present at all depths including >300 km
    # Evaluate: all bands significant including deep
    all_significant = len(significant_bands) == 4
    if all_significant:
        geometric_deep = "supported"
    elif deep_band_significant:
        geometric_deep = "partially"
    else:
        geometric_deep = "not supported"

    # Prediction 3 — Zhan & Shearer pattern:
    # Deep events (>300 km) show seasonality with mean phase in April–September
    # Evaluate: deep band significant AND mean_phase in [0.23, 0.67]
    if deep_band_significant and deep_mean_phase_in_apr_sep:
        zhan_shearer = "supported"
    elif deep_band_significant or deep_mean_phase_in_apr_sep:
        zhan_shearer = "partially"
    else:
        zhan_shearer = "not supported"

    logger.info("Prediction support:")
    logger.info("  surface_loading: %s", surface_loading)
    logger.info("  geometric_deep_forcing: %s", geometric_deep)
    logger.info("  zhan_shearer_pattern: %s", zhan_shearer)

    return {
        "cramer_v_by_band": cramer_v_by_band,
        "spearman_rho": rho,
        "spearman_p": p_spearman,
        "monotonicity_classification": monotonicity,
        "significant_bands": significant_bands,
        "deep_band_significant": bool(deep_band_significant),
        "deep_mean_phase": deep_mean_phase,
        "deep_mean_phase_in_apr_sep": bool(deep_mean_phase_in_apr_sep),
        "prediction_support": {
            "surface_loading": surface_loading,
            "geometric_deep_forcing": geometric_deep,
            "zhan_shearer_pattern": zhan_shearer,
        },
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load catalog, split by depth, compute per-band stats, run trend analysis, write JSON."""
    # --- Load catalog ---
    df = load_catalog(RAW_PATH)

    # --- Split by depth ---
    logger.info("Partitioning catalog by depth band")
    band_dfs, n_depth_null = split_by_depth(df)

    # --- Per-band statistics ---
    logger.info("Computing per-band statistics at k=16, 24, 32")
    band_stats = compute_all_band_stats(band_dfs, n_depth_null)

    # --- Trend analysis ---
    logger.info("Computing depth trend and monotonicity analysis")
    trend_analysis = compute_trend_analysis(band_stats)

    # --- Assemble results ---
    results: dict[str, Any] = {
        "case": "B4",
        "title": "Depth Stratification — Surface Loading Penetration Test",
        "solar_year_secs": SOLAR_YEAR_SECS,
        "band_stats": band_stats,
        "n_depth_null": n_depth_null,
        "trend_analysis": trend_analysis,
    }

    # --- Write JSON ---
    output_path = OUTPUT_DIR / "case-b4-results.json"
    with open(output_path, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", output_path)


if __name__ == "__main__":
    main()
