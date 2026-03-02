"""
Case A1: Schuster Spectrum and MFPA Periodicity Analysis — Main Analysis Script.

Applies cluster-robust Schuster Spectrum Test (Park et al. 2021) and Modified
Fourier Power Analysis (Dutilleul et al. 2015) to three ISC-GEM catalog versions
(raw, G-K declustered, A1b-informed declustered) to characterize periodic structure
in the solar seismic signal.

Usage:
    python topic-a2/src/case-a1-analysis.py

Outputs:
    topic-a2/output/case-a1-results.json
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import importlib.util
import numpy as np
import pandas as pd

# Project path setup
BASE_DIR = Path(__file__).resolve().parent.parent


def _load_module(name: str, filename: str):
    """Dynamically load a hyphen-named module from the src directory."""
    path = BASE_DIR / "src" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_schuster_mod = _load_module("case_a1_schuster", "case-a1-schuster.py")
schuster_spectrum = _schuster_mod.schuster_spectrum
schuster_explicit_tests = _schuster_mod.schuster_explicit_tests

_mfpa_mod = _load_module("case_a1_mfpa", "case-a1-mfpa.py")
mfpa_scan = _mfpa_mod.mfpa_scan

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("case-a1-analysis")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
JULIAN_YEAR_SECS: float = 31_557_600.0   # Julian constant (confirmed pre-run)
REFERENCE_EPOCH = pd.Timestamp("1950-01-01 00:00:00", tz="UTC")

# Data paths
RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GK_PATH = (
    BASE_DIR.parent
    / "data"
    / "iscgem"
    / "declustering-algorithm"
    / "mainshocks_G-K_global.csv"
)
A1B_PATH = (
    BASE_DIR.parent
    / "data"
    / "iscgem"
    / "declustering-algorithm"
    / "mainshocks_a1b_global.csv"
)

OUTPUT_PATH = BASE_DIR / "output" / "case-a1-results.json"

EXPECTED_COUNTS = {
    "raw": 9210,
    "gk_mainshocks": 5883,
    "a1b_mainshocks": 7137,
}

CATALOG_LABELS = {
    "raw": "Raw ISC-GEM catalog",
    "gk_mainshocks": "G-K mainshocks",
    "a1b_mainshocks": "A1b mainshocks",
}


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
# Time conversion
# ---------------------------------------------------------------------------

def to_event_time_days(event_at: pd.Series) -> pd.Series:
    """Convert UTC event timestamps to decimal days since 1950-01-01 00:00 UTC.

    Parameters
    ----------
    event_at : pd.Series
        Timezone-aware UTC datetime series.

    Returns
    -------
    pd.Series
        Decimal days since reference epoch.
    """
    delta = event_at - REFERENCE_EPOCH
    return delta.dt.total_seconds() / 86400.0


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_catalog(path: Path, label: str, expected_n: int) -> pd.DataFrame:
    """Load a catalog CSV, parse event_at as UTC, and validate row count.

    Parameters
    ----------
    path : Path
        Filesystem path to the CSV file.
    label : str
        Human-readable catalog label for log messages.
    expected_n : int
        Expected number of rows (assertion check).

    Returns
    -------
    pd.DataFrame
        Loaded catalog with ``event_at`` parsed as timezone-aware UTC datetime.

    Raises
    ------
    AssertionError
        If the row count does not match ``expected_n``.
    FileNotFoundError
        If the CSV file does not exist.
    """
    logger.info("Loading %s from %s", label, path)
    df = pd.read_csv(path, parse_dates=["event_at"])

    # Ensure timezone-aware UTC
    if df["event_at"].dt.tz is None:
        df["event_at"] = df["event_at"].dt.tz_localize("UTC")
    else:
        df["event_at"] = df["event_at"].dt.tz_convert("UTC")

    n = len(df)
    logger.info("  Loaded %d rows (expected %d)", n, expected_n)
    assert n == expected_n, (
        f"{label}: expected {expected_n} rows, got {n}"
    )
    return df


# ---------------------------------------------------------------------------
# Schuster analysis for one catalog
# ---------------------------------------------------------------------------

def run_schuster(
    event_time_days: np.ndarray,
    label: str,
) -> dict:
    """Run full Schuster spectrum and explicit period tests for one catalog.

    Parameters
    ----------
    event_time_days : np.ndarray
        Sorted decimal day times since reference epoch.
    label : str
        Catalog label for log messages.

    Returns
    -------
    dict with keys "spectrum" and "explicit_tests".
    """
    logger.info("Running Schuster spectrum for %s (%d events)", label, len(event_time_days))
    spectrum = schuster_spectrum(event_time_days)
    explicit = schuster_explicit_tests(event_time_days)
    # Enrich explicit tests with n_events / n_clusters at the catalog level
    n_events = len(event_time_days)
    n_clusters = explicit.get("annual_365", {}).get("n_clusters", None)
    logger.info(
        "  Schuster: n_events=%d, n_clusters=%s (annual period)",
        n_events,
        n_clusters,
    )
    return {
        "n_events": int(n_events),
        "n_clusters_at_annual": int(n_clusters) if n_clusters is not None else None,
        "spectrum": spectrum,
        "explicit_tests": explicit,
    }


# ---------------------------------------------------------------------------
# MFPA analysis for one catalog
# ---------------------------------------------------------------------------

def run_mfpa(
    event_time_days: np.ndarray,
    label: str,
) -> dict:
    """Run the full MFPA periodogram scan for one catalog.

    Parameters
    ----------
    event_time_days : np.ndarray
        Sorted decimal day times since reference epoch.
    label : str
        Catalog label for log messages.

    Returns
    -------
    dict with keys "significant_periods" and "spectrum".
    """
    logger.info("Running MFPA scan for %s (%d events)", label, len(event_time_days))
    result = mfpa_scan(event_time_days)
    logger.info(
        "  MFPA: %d significant periods (>p95)", len(result["significant_periods"])
    )
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Execute Case A1 analysis: Schuster spectrum and MFPA on all three catalogs."""
    logger.info("=== Case A1: Schuster Spectrum and MFPA Analysis ===")

    # Load catalogs
    df_raw = load_catalog(RAW_PATH, "Raw ISC-GEM", EXPECTED_COUNTS["raw"])
    df_gk = load_catalog(GK_PATH, "G-K mainshocks", EXPECTED_COUNTS["gk_mainshocks"])
    df_a1b = load_catalog(A1B_PATH, "A1b mainshocks", EXPECTED_COUNTS["a1b_mainshocks"])

    # Compute solar phases
    for df, name in [(df_raw, "raw"), (df_gk, "gk_mainshocks"), (df_a1b, "a1b_mainshocks")]:
        df["solar_phase"] = compute_solar_phase(df["solar_secs"])
        logger.info("%s: solar_phase range [%.4f, %.4f]", name, df["solar_phase"].min(), df["solar_phase"].max())

    # Compute event_time_days (sort each catalog by event_at first)
    for df in [df_raw, df_gk, df_a1b]:
        df.sort_values("event_at", inplace=True)
        df.reset_index(drop=True, inplace=True)

    df_raw["event_time_days"] = to_event_time_days(df_raw["event_at"])
    df_gk["event_time_days"] = to_event_time_days(df_gk["event_at"])
    df_a1b["event_time_days"] = to_event_time_days(df_a1b["event_at"])

    # Verify monotonicity
    for df, name in [(df_raw, "raw"), (df_gk, "gk"), (df_a1b, "a1b")]:
        diffs = df["event_time_days"].diff().dropna()
        assert (diffs >= 0).all(), f"{name}: event_time_days is not non-decreasing"
        logger.info("%s: event_time_days monotonic check passed", name)

    # Extract numpy arrays
    t_raw = df_raw["event_time_days"].to_numpy()
    t_gk = df_gk["event_time_days"].to_numpy()
    t_a1b = df_a1b["event_time_days"].to_numpy()

    # --- Schuster analysis ---
    schuster_raw = run_schuster(t_raw, "raw")
    schuster_gk = run_schuster(t_gk, "gk_mainshocks")
    schuster_a1b = run_schuster(t_a1b, "a1b_mainshocks")

    # --- MFPA analysis ---
    mfpa_raw = run_mfpa(t_raw, "raw")
    mfpa_gk = run_mfpa(t_gk, "gk_mainshocks")
    mfpa_a1b = run_mfpa(t_a1b, "a1b_mainshocks")

    # --- Assemble results ---
    results = {
        "case": "A1",
        "title": "Schuster Spectrum and MFPA Periodicity Analysis",
        "julian_year_secs": JULIAN_YEAR_SECS,
        "reference_epoch": str(REFERENCE_EPOCH),
        "catalog_counts": {
            "raw": int(len(df_raw)),
            "gk_mainshocks": int(len(df_gk)),
            "a1b_mainshocks": int(len(df_a1b)),
        },
        "schuster": {
            "raw": schuster_raw,
            "gk_mainshocks": schuster_gk,
            "a1b_mainshocks": schuster_a1b,
        },
        "mfpa": {
            "raw": mfpa_raw,
            "gk_mainshocks": mfpa_gk,
            "a1b_mainshocks": mfpa_a1b,
        },
    }

    # Write output
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", OUTPUT_PATH)
    logger.info("=== Case A1 analysis complete ===")


if __name__ == "__main__":
    main()
