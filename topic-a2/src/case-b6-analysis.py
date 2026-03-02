"""
Case B6: Rolling Window Stationarity Test

Tests whether the solar-phase signal in the ISC-GEM catalog is stationary
across the 72-year record (1950-2021) using a sliding 10-year window.

Computes Rayleigh statistic, mean vector length (R), mean phase angle, and
chi-square supplementary test for each window position.
"""

import json
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"

JULIAN_YEAR_SECS = 31_557_600.0
WINDOW_YEARS = 10
STEP_YEARS = 1
K_BINS = 24


def load_catalog(path: Path) -> pd.DataFrame:
    """Load ISC-GEM catalog and compute event_year and solar phase.

    Args:
        path: Path to the raw CSV file.

    Returns:
        DataFrame with parsed event_at, event_year, and phase columns.
    """
    logger.info("Loading catalog from %s", path)
    df = pd.read_csv(path, parse_dates=["event_at"])
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    n = len(df)
    logger.info("Loaded %d rows", n)
    assert n == 9210, f"Expected 9210 rows, got {n}"

    nat_count = df["event_at"].isna().sum()
    if nat_count > 0:
        logger.warning("%d NaT values in event_at", nat_count)
    assert nat_count == 0, "NaT values found in event_at"

    df["event_year"] = df["event_at"].dt.year
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    logger.info("Phase range: [%.6f, %.6f]", df["phase"].min(), df["phase"].max())
    return df


def compute_rayleigh(phases: np.ndarray) -> dict:
    """Compute Rayleigh statistic for a set of circular phase values.

    Args:
        phases: Array of phase values in [0, 1).

    Returns:
        Dict with R, rayleigh_z, p_rayleigh, mean_angle_rad, mean_phase.
    """
    n = len(phases)
    angles = 2.0 * np.pi * phases
    mean_cos = np.mean(np.cos(angles))
    mean_sin = np.mean(np.sin(angles))
    R = np.sqrt(mean_cos**2 + mean_sin**2)
    rayleigh_z = n * R**2
    p_rayleigh = np.exp(-rayleigh_z)
    mean_angle_rad = np.arctan2(mean_sin, mean_cos)
    mean_phase = (mean_angle_rad / (2.0 * np.pi)) % 1.0
    return {
        "R": R,
        "rayleigh_z": rayleigh_z,
        "p_rayleigh": p_rayleigh,
        "mean_angle_rad": mean_angle_rad,
        "mean_phase": mean_phase,
    }


def compute_chi2_k24(phases: np.ndarray, k: int = K_BINS) -> tuple[float, float]:
    """Compute chi-square test with phase-normalized binning.

    Args:
        phases: Array of phase values in [0, 1).
        k: Number of bins.

    Returns:
        Tuple of (chi2_stat, p_chi2).
    """
    n = len(phases)
    bin_indices = (phases * k).astype(int)
    bin_indices = np.clip(bin_indices, 0, k - 1)
    observed, _ = np.histogram(bin_indices, bins=range(k + 1))
    expected = np.full(k, n / k)
    chi2_stat, p_chi2 = scipy.stats.chisquare(observed, expected)
    return float(chi2_stat), float(p_chi2)


def compute_circular_std(mean_phases: np.ndarray) -> tuple[float, float]:
    """Compute circular standard deviation of mean phase angles.

    Args:
        mean_phases: Array of mean phase values in [0, 1).

    Returns:
        Tuple of (circ_var, circ_std_deg).
    """
    angles = 2.0 * np.pi * mean_phases
    R_mean = np.sqrt(np.mean(np.cos(angles))**2 + np.mean(np.sin(angles))**2)
    circ_var = 1.0 - R_mean
    # Circular standard deviation in degrees
    circ_std_deg = math.sqrt(-2.0 * math.log(max(1.0 - circ_var, 1e-15))) * 180.0 / math.pi
    return float(circ_var), float(circ_std_deg)


def classify_stationarity(
    n_sig: int,
    n_windows: int,
    circ_std_deg: float,
) -> str:
    """Classify stationarity based on % significant windows and circular std.

    Args:
        n_sig: Number of significant windows (p < 0.05).
        n_windows: Total number of windows.
        circ_std_deg: Circular standard deviation of mean phases in degrees.

    Returns:
        Classification string: 'stationary', 'partially stationary', or 'non-stationary'.
    """
    pct_sig = n_sig / n_windows
    if pct_sig >= 0.70 and circ_std_deg < 20.0:
        return "stationary"
    elif pct_sig < 0.30 or circ_std_deg > 40.0:
        return "non-stationary"
    else:
        return "partially stationary"


def run_rolling_windows(df: pd.DataFrame) -> list[dict]:
    """Run sliding window analysis across the full catalog.

    Args:
        df: Catalog DataFrame with event_year and phase columns.

    Returns:
        List of per-window result dicts.
    """
    # range(1950, 2021 - WINDOW_YEARS + 1) = range(1950, 2012) → 62 windows per spec
    window_starts = list(range(1950, 2021 - WINDOW_YEARS + 1))
    logger.info("Running %d windows (start years %d–%d)", len(window_starts),
                window_starts[0], window_starts[-1])

    windows = []
    for y in window_starts:
        y_end = y + WINDOW_YEARS - 1
        subset = df[(df["event_year"] >= y) & (df["event_year"] < y + WINDOW_YEARS)]
        n_window = len(subset)

        phases = subset["phase"].values
        ray = compute_rayleigh(phases)
        chi2_stat, p_chi2 = compute_chi2_k24(phases)

        is_1970s = (1970 <= y <= 1979)

        windows.append({
            "window_start": int(y),
            "window_end": int(y_end),
            "n": int(n_window),
            "rayleigh_R": float(ray["R"]),
            "rayleigh_z": float(ray["rayleigh_z"]),
            "p_rayleigh": float(ray["p_rayleigh"]),
            "mean_phase": float(ray["mean_phase"]),
            "chi2_k24": float(chi2_stat),
            "p_chi2_k24": float(p_chi2),
            "is_1970s_window": bool(is_1970s),
        })
        logger.debug(
            "Window %d–%d: n=%d R=%.4f p=%.4f mean_phase=%.4f",
            y, y_end, n_window, ray["R"], ray["p_rayleigh"], ray["mean_phase"],
        )

    return windows


def compute_stationarity(windows: list[dict]) -> dict:
    """Compute stationarity classification and 1970s anomaly check.

    Args:
        windows: List of per-window result dicts.

    Returns:
        Stationarity summary dict.
    """
    n_windows = len(windows)
    bonferroni_threshold = 0.05 / n_windows

    n_sig_p05 = sum(1 for w in windows if w["p_rayleigh"] < 0.05)
    n_bonferroni = sum(1 for w in windows if w["p_rayleigh"] < bonferroni_threshold)

    mean_phases = np.array([w["mean_phase"] for w in windows])
    circ_var, circ_std_deg = compute_circular_std(mean_phases)

    classification = classify_stationarity(n_sig_p05, n_windows, circ_std_deg)

    # 1970s anomaly check
    r_1970s = [w["rayleigh_R"] for w in windows if w["is_1970s_window"]]
    r_non_1970s = [w["rayleigh_R"] for w in windows if not w["is_1970s_window"]]
    mean_r_1970s = float(np.mean(r_1970s)) if r_1970s else 0.0
    mean_r_non_1970s = float(np.mean(r_non_1970s)) if r_non_1970s else 0.0
    ratio = mean_r_1970s / mean_r_non_1970s if mean_r_non_1970s > 0 else 0.0
    anomaly_flagged = ratio > 1.5

    logger.info(
        "Stationarity: %s | n_sig_p05=%d/%d | circ_std=%.2f° | 1970s ratio=%.3f flagged=%s",
        classification, n_sig_p05, n_windows, circ_std_deg, ratio, anomaly_flagged,
    )

    return {
        "n_windows": int(n_windows),
        "n_significant_rayleigh_p05": int(n_sig_p05),
        "n_bonferroni_significant": int(n_bonferroni),
        "bonferroni_threshold": float(bonferroni_threshold),
        "circular_std_deg": float(circ_std_deg),
        "circular_variance": float(circ_var),
        "classification": classification,
        "1970s_mean_R": float(mean_r_1970s),
        "non_1970s_mean_R": float(mean_r_non_1970s),
        "1970s_anomaly_ratio": float(ratio),
        "1970s_anomaly_flagged": bool(anomaly_flagged),
    }


def main() -> None:
    """Main entry point for Case B6 analysis."""
    df = load_catalog(RAW_PATH)

    windows = run_rolling_windows(df)
    assert len(windows) == 62, f"Expected 62 windows, got {len(windows)}"

    stationarity = compute_stationarity(windows)

    results = {
        "case": "B6",
        "title": "Rolling Window Stationarity Test",
        "parameters": {
            "n_catalog": int(len(df)),
            "window_years": WINDOW_YEARS,
            "step_years": STEP_YEARS,
            "k_bins_chi2": K_BINS,
            "julian_year_secs": JULIAN_YEAR_SECS,
            "year_range": [1950, 2021],
            "window_start_range": [1950, 2012],
        },
        "stationarity": stationarity,
        "windows": windows,
    }

    out_path = BASE_DIR / "output" / "case-b6-results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", out_path)


if __name__ == "__main__":
    main()
