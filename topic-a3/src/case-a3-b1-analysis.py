"""
Case A3.B1: Rolling-Window Chi-Square Repeat

Repeats the A2.B6 rolling-window stationarity analysis with chi-square (k=24)
promoted to primary statistic and Rayleigh demoted to secondary. Introduces
interval-level tracking within each window and sequence density analysis.

Window parameters: 10-year window, 1-year stride, 62 windows (1950–2011 starts).
"""

import json
import logging
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent

CATALOGS: dict[str, Path] = {
    "raw": BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv",
    "gk": BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv",
    "reasenberg": BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv",
    "a1b": BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_a1b-seq_global.csv",
}

EXPECTED_N: dict[str, int] = {
    "raw": 9210,
    "gk": 5883,
    "reasenberg": 8265,
    "a1b": 7137,
}

WINDOW_YEARS: int = 10
STEP_YEARS: int = 1
K_BINS: int = 24
JULIAN_YEAR_SECS: float = 31_557_600.0

# Window start years: 1950 through 2011 (last window = 2011–2020)
WINDOW_START_YEARS: list[int] = list(range(1950, 2021 - WINDOW_YEARS + 1))

# 1970s windows: start year in [1970, 1980)
YEARS_1970S = set(range(1970, 1980))

# A1b baseline interval bin ranges (0-based, k=24)
INTERVAL_BINS: dict[str, list[int]] = {
    "interval_1": [4, 5],   # phase [0.1667, 0.2500) — March equinox region
    "interval_2": [15],     # phase [0.6250, 0.6667) — mid-August region
    "interval_3": [21],     # phase [0.8750, 0.9167) — late-November region
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_catalog(key: str, path: Path) -> pd.DataFrame:
    """
    Load a catalog CSV, parse event_at as UTC datetime, compute event_year and phase.

    Parameters
    ----------
    key : str
        Catalog identifier used for logging and expected-n assertions.
    path : Path
        Absolute path to the catalog CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: all originals plus `event_year` (int) and `phase` (float in [0, 1)).
    """
    logger.info("Loading catalog '%s' from %s", key, path)
    df = pd.read_csv(path, parse_dates=["event_at"])
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)

    n_rows = len(df)
    logger.info("  Loaded %d rows (expected %d)", n_rows, EXPECTED_N[key])
    assert n_rows == EXPECTED_N[key], (
        f"Row count mismatch for '{key}': got {n_rows}, expected {EXPECTED_N[key]}"
    )

    df["event_year"] = df["event_at"].dt.year.astype(int)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    logger.info("  Phase range: [%.6f, %.6f]", df["phase"].min(), df["phase"].max())
    return df


# ---------------------------------------------------------------------------
# Per-window statistics
# ---------------------------------------------------------------------------

def compute_window_stats(
    key: str,
    subset: pd.DataFrame,
    window_start: int,
) -> dict[str, Any]:
    """
    Compute per-window statistics for one catalog window.

    Parameters
    ----------
    key : str
        Catalog identifier; used to determine whether sequence density is computed.
    subset : pd.DataFrame
        Rows of the catalog falling within the window.
    window_start : int
        Start year of the window (inclusive).

    Returns
    -------
    dict
        Per-window statistics record.
    """
    n_window = len(subset)
    window_end = window_start + WINDOW_YEARS

    # --- Phase-normalized bin counts ---
    bin_indices = np.floor(subset["phase"].values * K_BINS).astype(int) % K_BINS
    observed = np.bincount(bin_indices, minlength=K_BINS).astype(float)
    expected = np.full(K_BINS, n_window / K_BINS)

    # --- Chi-square (primary) ---
    chi2_stat, p_chi2 = scipy.stats.chisquare(observed, expected)

    # --- Rayleigh (secondary) ---
    angles = 2.0 * np.pi * subset["phase"].values
    mean_cos = float(np.mean(np.cos(angles)))
    mean_sin = float(np.mean(np.sin(angles)))
    R = float(np.sqrt(mean_cos ** 2 + mean_sin ** 2))
    rayleigh_z = float(n_window * R ** 2)
    p_rayleigh = float(np.exp(-rayleigh_z))
    mean_angle_rad = float(np.arctan2(mean_sin, mean_cos))
    mean_phase = float((mean_angle_rad / (2.0 * np.pi)) % 1.0)

    # --- Interval-level tracking ---
    interval_stats: dict[str, Any] = {}
    for iname, bins in INTERVAL_BINS.items():
        interval_count = int(observed[bins].sum())
        interval_expected = float(len(bins) * (n_window / K_BINS))
        interval_z = float(
            (interval_count - interval_expected) / math.sqrt(interval_expected)
        )
        interval_elevated = bool(interval_z > 1.96)
        interval_stats[f"{iname}_count"] = interval_count
        interval_stats[f"{iname}_z"] = interval_z
        interval_stats[f"{iname}_elevated"] = interval_elevated

    # --- Sequence density (declustered catalogs only) ---
    if key == "raw":
        mean_aftershock_count: float | None = None
        total_aftershock_count: int | None = None
        pct_isolated: float | None = None
    else:
        mean_aftershock_count = float(subset["aftershock_count"].mean())
        total_aftershock_count = int(subset["aftershock_count"].sum())
        pct_isolated = float((subset["aftershock_count"] == 0).mean())

    return {
        "window_start": window_start,
        "window_end": window_end,
        "n": n_window,
        "chi2_k24": float(chi2_stat),
        "p_chi2_k24": float(p_chi2),
        "rayleigh_R": R,
        "rayleigh_z": rayleigh_z,
        "p_rayleigh": p_rayleigh,
        "mean_phase": mean_phase,
        "interval_1_count": interval_stats["interval_1_count"],
        "interval_1_z": interval_stats["interval_1_z"],
        "interval_1_elevated": interval_stats["interval_1_elevated"],
        "interval_2_count": interval_stats["interval_2_count"],
        "interval_2_z": interval_stats["interval_2_z"],
        "interval_2_elevated": interval_stats["interval_2_elevated"],
        "interval_3_count": interval_stats["interval_3_count"],
        "interval_3_z": interval_stats["interval_3_z"],
        "interval_3_elevated": interval_stats["interval_3_elevated"],
        "mean_aftershock_count": mean_aftershock_count,
        "total_aftershock_count": total_aftershock_count,
        "pct_isolated": pct_isolated,
        "is_1970s_window": window_start in YEARS_1970S,
    }


def compute_rolling_windows(key: str, df: pd.DataFrame) -> list[dict[str, Any]]:
    """
    Compute statistics for all 62 rolling windows for one catalog.

    Parameters
    ----------
    key : str
        Catalog identifier.
    df : pd.DataFrame
        Full catalog DataFrame with `event_year` and `phase` columns.

    Returns
    -------
    list of dict
        One record per window.
    """
    logger.info("Computing rolling windows for catalog '%s'...", key)
    windows = []
    for y in WINDOW_START_YEARS:
        subset = df[(df["event_year"] >= y) & (df["event_year"] < y + WINDOW_YEARS)]
        if len(subset) < 100:
            logger.warning(
                "  Window %d–%d for '%s' has only %d events (< 100).",
                y, y + WINDOW_YEARS, key, len(subset),
            )
        record = compute_window_stats(key, subset, y)
        windows.append(record)
    logger.info("  %d windows computed for '%s'.", len(windows), key)
    return windows


# ---------------------------------------------------------------------------
# Stationarity classification
# ---------------------------------------------------------------------------

def classify_stationarity(
    key: str,
    windows: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Classify stationarity for one catalog from its rolling-window records.

    Uses chi-square as the primary statistic and circular SD as secondary.
    When criteria conflict, applies the more conservative classification.

    Parameters
    ----------
    key : str
        Catalog identifier.
    windows : list of dict
        Per-window statistics records.

    Returns
    -------
    dict
        Stationarity summary for this catalog.
    """
    n_windows = len(windows)

    # Chi-square significance counts
    n_significant_chi2_p05 = sum(1 for w in windows if w["p_chi2_k24"] < 0.05)
    bonferroni_threshold = 0.05 / n_windows
    n_bonferroni_significant_chi2 = sum(
        1 for w in windows if w["p_chi2_k24"] < bonferroni_threshold
    )
    n_significant_rayleigh_p05 = sum(1 for w in windows if w["p_rayleigh"] < 0.05)

    pct_chi2_sig = n_significant_chi2_p05 / n_windows
    pct_rayleigh_sig = n_significant_rayleigh_p05 / n_windows

    # Circular phase stability
    angles_w = [2.0 * np.pi * w["mean_phase"] for w in windows]
    mean_cos_w = float(np.mean(np.cos(angles_w)))
    mean_sin_w = float(np.mean(np.sin(angles_w)))
    R_of_means = float(np.sqrt(mean_cos_w ** 2 + mean_sin_w ** 2))
    circ_var = 1.0 - R_of_means
    if circ_var < 1.0:
        circ_std_deg = float(math.sqrt(-2.0 * math.log(max(1.0 - circ_var, 1e-300))) * 180.0 / math.pi)
    else:
        circ_std_deg = 180.0

    # Stationarity classification (chi-square primary)
    chi2_class: str
    if pct_chi2_sig >= 0.70:
        chi2_class = "stationary"
    elif pct_chi2_sig >= 0.30:
        chi2_class = "partially stationary"
    else:
        chi2_class = "non-stationary"

    circ_class: str
    if circ_std_deg < 20.0:
        circ_class = "stationary"
    elif circ_std_deg <= 40.0:
        circ_class = "partially stationary"
    else:
        circ_class = "non-stationary"

    # Rank ordering: stationary > partially stationary > non-stationary
    rank = {"stationary": 2, "partially stationary": 1, "non-stationary": 0}
    final_class = chi2_class if rank[chi2_class] <= rank[circ_class] else circ_class

    if chi2_class != circ_class:
        basis = (
            f"Conflict: chi2 criterion → '{chi2_class}' (pct={pct_chi2_sig:.3f}), "
            f"circ_std criterion → '{circ_class}' (circ_std_deg={circ_std_deg:.2f}°). "
            f"Applied more conservative classification: '{final_class}'."
        )
    elif chi2_class == "stationary":
        basis = (
            f"Both criteria agree: pct_chi2_sig={pct_chi2_sig:.3f} >= 0.70 "
            f"and circ_std_deg={circ_std_deg:.2f}° < 20°."
        )
    elif chi2_class == "partially stationary":
        basis = (
            f"Both criteria agree (partially stationary): pct_chi2_sig={pct_chi2_sig:.3f} "
            f"in [0.30, 0.70), circ_std_deg={circ_std_deg:.2f}°."
        )
    else:
        basis = (
            f"Both criteria agree (non-stationary): pct_chi2_sig={pct_chi2_sig:.3f} < 0.30 "
            f"or circ_std_deg={circ_std_deg:.2f}° > 40°."
        )

    # Interval stationarity summaries
    interval_stats: dict[str, Any] = {}
    for iname in ["interval_1", "interval_2", "interval_3"]:
        n_elev = sum(1 for w in windows if w[f"{iname}_elevated"])
        pct_elev = n_elev / n_windows
        if pct_elev >= 0.70:
            iclass = "globally elevated"
        elif pct_elev >= 0.30:
            iclass = "partially elevated"
        else:
            iclass = "absent"
        interval_stats[f"{iname}_n_elevated"] = n_elev
        interval_stats[f"{iname}_pct_elevated"] = pct_elev
        interval_stats[f"{iname}_classification"] = iclass

    # 1970s anomaly check (Rayleigh R for comparability with B6)
    r_1970s = [w["rayleigh_R"] for w in windows if w["is_1970s_window"]]
    r_non1970s = [w["rayleigh_R"] for w in windows if not w["is_1970s_window"]]
    mean_R_1970s = float(np.mean(r_1970s)) if r_1970s else 0.0
    mean_R_non1970s = float(np.mean(r_non1970s)) if r_non1970s else 0.0
    ratio_1970s = mean_R_1970s / mean_R_non1970s if mean_R_non1970s > 0 else 0.0
    flagged_1970s = bool(ratio_1970s > 1.5)

    logger.info(
        "  [%s] chi2_sig_pct=%.3f, circ_std=%.2f°, classification='%s'",
        key, pct_chi2_sig, circ_std_deg, final_class,
    )

    return {
        "n_windows": n_windows,
        "n_significant_chi2_p05": n_significant_chi2_p05,
        "pct_significant_chi2_p05": pct_chi2_sig,
        "n_bonferroni_significant_chi2": n_bonferroni_significant_chi2,
        "n_significant_rayleigh_p05": n_significant_rayleigh_p05,
        "pct_significant_rayleigh_p05": pct_rayleigh_sig,
        "circular_std_deg": circ_std_deg,
        "classification": final_class,
        "classification_basis": basis,
        "interval_1_pct_elevated": interval_stats["interval_1_pct_elevated"],
        "interval_2_pct_elevated": interval_stats["interval_2_pct_elevated"],
        "interval_3_pct_elevated": interval_stats["interval_3_pct_elevated"],
        "interval_1_classification": interval_stats["interval_1_classification"],
        "interval_2_classification": interval_stats["interval_2_classification"],
        "interval_3_classification": interval_stats["interval_3_classification"],
        "1970s_mean_R": mean_R_1970s,
        "non_1970s_mean_R": mean_R_non1970s,
        "1970s_anomaly_ratio": ratio_1970s,
        "1970s_anomaly_flagged": flagged_1970s,
    }


# ---------------------------------------------------------------------------
# Sequence density correlation
# ---------------------------------------------------------------------------

def compute_sequence_density(
    key: str,
    windows: list[dict[str, Any]],
) -> dict[str, Any] | None:
    """
    Compute Pearson correlation between chi-square p-value and mean aftershock count.

    Skips raw catalog (returns None). For declustered catalogs, correlates across
    all 62 windows, identifies the 10 most significant chi-square windows, and
    checks whether those windows have elevated aftershock density.

    Parameters
    ----------
    key : str
        Catalog identifier.
    windows : list of dict
        Per-window statistics records.

    Returns
    -------
    dict or None
        Sequence density statistics, or None for the raw catalog.
    """
    if key == "raw":
        return None

    p_values = np.array([w["p_chi2_k24"] for w in windows])
    aftershock_means = np.array([w["mean_aftershock_count"] for w in windows], dtype=float)

    r, p = scipy.stats.pearsonr(p_values, aftershock_means)

    # Top-10 most significant chi-square windows (lowest p-value)
    top10_idx = np.argsort(p_values)[:10]
    mean_aftershock_all = float(np.mean(aftershock_means))
    mean_aftershock_top10 = float(np.mean(aftershock_means[top10_idx]))
    elevated = bool(mean_aftershock_top10 > 1.5 * mean_aftershock_all)

    logger.info(
        "  [%s] Sequence density: r=%.4f, p=%.4f; top10 mean=%.3f vs overall=%.3f (elevated=%s)",
        key, r, p, mean_aftershock_top10, mean_aftershock_all, elevated,
    )

    return {
        "r_chi2_vs_aftershock": float(r),
        "p_chi2_vs_aftershock": float(p),
        "mean_aftershock_count_all_windows": mean_aftershock_all,
        "mean_aftershock_count_top10_chi2_windows": mean_aftershock_top10,
        "seq_density_elevated_in_high_chi2_windows": elevated,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the full A3.B1 rolling-window chi-square analysis and write results JSON."""
    logger.info("=== Case A3.B1: Rolling-Window Chi-Square Repeat ===")
    logger.info("Window parameters: %d-year window, %d-year stride, %d windows",
                WINDOW_YEARS, STEP_YEARS, len(WINDOW_START_YEARS))

    results: dict[str, Any] = {
        "meta": {
            "case": "A3.B1",
            "title": "Rolling-Window Chi-Square Repeat",
            "window_years": WINDOW_YEARS,
            "step_years": STEP_YEARS,
            "n_windows": len(WINDOW_START_YEARS),
            "k_bins": K_BINS,
            "julian_year_secs": JULIAN_YEAR_SECS,
            "bonferroni_threshold": 0.05 / len(WINDOW_START_YEARS),
            "window_start_range": [WINDOW_START_YEARS[0], WINDOW_START_YEARS[-1]],
            "interval_bins": INTERVAL_BINS,
        },
        "catalogs": {},
    }

    for cat_key, cat_path in CATALOGS.items():
        logger.info("--- Processing catalog: %s ---", cat_key)
        df = load_catalog(cat_key, cat_path)
        windows = compute_rolling_windows(cat_key, df)
        stationarity = classify_stationarity(cat_key, windows)
        seq_density = compute_sequence_density(cat_key, windows)

        cat_result: dict[str, Any] = {
            "stationarity": stationarity,
            "windows": windows,
        }
        if seq_density is not None:
            cat_result["sequence_density"] = seq_density

        results["catalogs"][cat_key] = cat_result

    # Write results JSON
    output_path = BASE_DIR / "output" / "case-a3-b1-results.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, allow_nan=False)
    logger.info("Results written to %s", output_path)
    logger.info("=== A3.B1 analysis complete ===")


if __name__ == "__main__":
    main()
