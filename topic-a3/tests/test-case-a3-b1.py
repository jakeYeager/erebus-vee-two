"""
Test suite for Case A3.B1: Rolling-Window Chi-Square Repeat.

All tests load from the pre-computed results JSON where possible,
or re-run lightweight data operations for correctness verification.
"""

import json
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.stats

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent

CATALOGS_PATHS: dict[str, Path] = {
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

JULIAN_YEAR_SECS = 31_557_600.0
RESULTS_PATH = BASE_DIR / "output" / "case-a3-b1-results.json"

YEARS_1970S = set(range(1970, 1980))


@pytest.fixture(scope="session")
def results() -> dict:
    """Load the case results JSON once for all tests."""
    with open(RESULTS_PATH, encoding="utf-8") as fh:
        return json.load(fh)


def load_catalog(key: str) -> pd.DataFrame:
    """Load a catalog CSV and compute phase column."""
    df = pd.read_csv(CATALOGS_PATHS[key], parse_dates=["event_at"])
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


# ---------------------------------------------------------------------------
# Data loading tests
# ---------------------------------------------------------------------------

def test_catalog_load_raw() -> None:
    """Raw catalog must have exactly 9210 rows and no NaT in event_at."""
    df = load_catalog("raw")
    assert len(df) == 9210, f"Expected 9210 rows, got {len(df)}"
    assert df["event_at"].isna().sum() == 0, "Found NaT values in event_at"


def test_catalog_load_declustered() -> None:
    """Declustered catalogs must have expected row counts and aftershock_count column."""
    for key in ("gk", "reasenberg", "a1b"):
        df = load_catalog(key)
        assert len(df) == EXPECTED_N[key], (
            f"[{key}] Expected {EXPECTED_N[key]} rows, got {len(df)}"
        )
        assert "aftershock_count" in df.columns, (
            f"[{key}] Missing column 'aftershock_count'"
        )


# ---------------------------------------------------------------------------
# Phase range tests
# ---------------------------------------------------------------------------

def test_phase_range() -> None:
    """All computed phases must be in [0.0, 1.0) for all 4 catalogs."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        df = load_catalog(key)
        phases = df["phase"].values
        assert np.all(phases >= 0.0), f"[{key}] Negative phase values found"
        assert np.all(phases < 1.0), f"[{key}] Phase values >= 1.0 found"


# ---------------------------------------------------------------------------
# Window count tests
# ---------------------------------------------------------------------------

def test_window_count(results: dict) -> None:
    """Each catalog must have exactly 62 rolling windows."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        assert len(windows) == 62, f"[{key}] Expected 62 windows, got {len(windows)}"


# ---------------------------------------------------------------------------
# Window n > 0
# ---------------------------------------------------------------------------

def test_window_n_positive(results: dict) -> None:
    """All window n values must be > 0; log a warning if any window n < 100."""
    logger = logging.getLogger(__name__)
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        for w in windows:
            assert w["n"] > 0, f"[{key}] Window {w['window_start']} has n=0"
            if w["n"] < 100:
                logger.warning(
                    "[%s] Window %d–%d has only %d events (< 100).",
                    key, w["window_start"], w["window_end"], w["n"],
                )


# ---------------------------------------------------------------------------
# Chi-square value tests
# ---------------------------------------------------------------------------

def test_chi2_nonneg(results: dict) -> None:
    """All chi-square statistics must be non-negative."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        for w in windows:
            assert w["chi2_k24"] >= 0.0, (
                f"[{key}] Window {w['window_start']}: chi2_k24={w['chi2_k24']} < 0"
            )


def test_chi2_p_bounds(results: dict) -> None:
    """All chi-square p-values must be in [0.0, 1.0]."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        for w in windows:
            assert 0.0 <= w["p_chi2_k24"] <= 1.0, (
                f"[{key}] Window {w['window_start']}: p_chi2_k24={w['p_chi2_k24']} out of bounds"
            )


# ---------------------------------------------------------------------------
# Rayleigh R bounds
# ---------------------------------------------------------------------------

def test_rayleigh_R_bounds(results: dict) -> None:
    """All Rayleigh R values must be in [0.0, 1.0]."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        for w in windows:
            assert 0.0 <= w["rayleigh_R"] <= 1.0, (
                f"[{key}] Window {w['window_start']}: rayleigh_R={w['rayleigh_R']} out of [0,1]"
            )


# ---------------------------------------------------------------------------
# Mean phase bounds
# ---------------------------------------------------------------------------

def test_mean_phase_bounds(results: dict) -> None:
    """All mean_phase values must be in [0.0, 1.0)."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        for w in windows:
            assert 0.0 <= w["mean_phase"] < 1.0, (
                f"[{key}] Window {w['window_start']}: mean_phase={w['mean_phase']} out of [0,1)"
            )


# ---------------------------------------------------------------------------
# Interval z-score and elevated flag types
# ---------------------------------------------------------------------------

def test_interval_z_type(results: dict) -> None:
    """Interval z-scores must be floats in all window records for all catalogs."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        for w in windows:
            for iname in ("interval_1_z", "interval_2_z", "interval_3_z"):
                assert isinstance(w[iname], float), (
                    f"[{key}] Window {w['window_start']}: {iname} is not float"
                )


def test_interval_elevated_type(results: dict) -> None:
    """Interval elevated flags must be bools in all window records for all catalogs."""
    for key in ("raw", "gk", "reasenberg", "a1b"):
        windows = results["catalogs"][key]["windows"]
        for w in windows:
            for iname in ("interval_1_elevated", "interval_2_elevated", "interval_3_elevated"):
                assert isinstance(w[iname], bool), (
                    f"[{key}] Window {w['window_start']}: {iname} is not bool"
                )


# ---------------------------------------------------------------------------
# 1970s window flagging
# ---------------------------------------------------------------------------

def test_1970s_flagging(results: dict) -> None:
    """Windows with start year 1970–1979 must be flagged; 1980–2012 must not."""
    raw_windows = results["catalogs"]["raw"]["windows"]
    for w in raw_windows:
        ws = w["window_start"]
        if ws in YEARS_1970S:
            assert w["is_1970s_window"] is True, (
                f"Window start={ws} should be is_1970s_window=True"
            )
        elif 1980 <= ws <= 2012:
            assert w["is_1970s_window"] is False, (
                f"Window start={ws} should be is_1970s_window=False"
            )


# ---------------------------------------------------------------------------
# Stationarity classification values
# ---------------------------------------------------------------------------

def test_stationarity_classification_all_catalogs(results: dict) -> None:
    """Each catalog must have a valid classification string."""
    valid = {"stationary", "partially stationary", "non-stationary"}
    for key in ("raw", "gk", "reasenberg", "a1b"):
        cls = results["catalogs"][key]["stationarity"]["classification"]
        assert cls in valid, f"[{key}] classification='{cls}' not in {valid}"


# ---------------------------------------------------------------------------
# Sequence density: raw=null, declustered=float
# ---------------------------------------------------------------------------

def test_sequence_density_declustered_only(results: dict) -> None:
    """
    Raw catalog windows must have null mean_aftershock_count;
    declustered catalogs must have non-negative float values.
    """
    # Raw: all null
    for w in results["catalogs"]["raw"]["windows"]:
        assert w["mean_aftershock_count"] is None, (
            f"Raw window {w['window_start']}: expected null mean_aftershock_count"
        )

    # Declustered: non-negative float
    for key in ("gk", "reasenberg", "a1b"):
        for w in results["catalogs"][key]["windows"]:
            v = w["mean_aftershock_count"]
            assert isinstance(v, float), (
                f"[{key}] Window {w['window_start']}: mean_aftershock_count should be float"
            )
            assert v >= 0.0, (
                f"[{key}] Window {w['window_start']}: mean_aftershock_count={v} < 0"
            )


# ---------------------------------------------------------------------------
# Statistical sanity tests
# ---------------------------------------------------------------------------

def test_chi2_uniform() -> None:
    """
    Chi-square on truly uniform phases should yield p > 0.05 in >= 90% of bootstrap trials.

    Generates 100 bootstrap trials of 1000 uniform-random phases each and
    checks chi-square at k=24.
    """
    rng = np.random.default_rng(42)
    n_trials = 100
    k = 24
    n_per_trial = 1000
    n_significant = 0
    for _ in range(n_trials):
        phases = rng.uniform(0.0, 1.0, n_per_trial)
        bin_indices = np.floor(phases * k).astype(int) % k
        observed = np.bincount(bin_indices, minlength=k).astype(float)
        expected = np.full(k, n_per_trial / k)
        _, p = scipy.stats.chisquare(observed, expected)
        if p < 0.05:
            n_significant += 1
    pct_sig = n_significant / n_trials
    assert pct_sig <= 0.10, (
        f"Uniform chi-square: {pct_sig:.2%} significant (expected <= 10%)"
    )


def test_chi2_concentrated() -> None:
    """
    Chi-square on phases concentrated near 0.635 (bins 15–16) must yield p < 0.001.

    Uses 1000 samples drawn from a tight normal distribution centered at phase 0.635.
    """
    rng = np.random.default_rng(99)
    n = 1000
    k = 24
    # Tight concentration around phase 0.635, sigma=0.01
    raw_phases = rng.normal(loc=0.635, scale=0.01, size=n)
    phases = raw_phases % 1.0
    bin_indices = np.floor(phases * k).astype(int) % k
    observed = np.bincount(bin_indices, minlength=k).astype(float)
    expected = np.full(k, n / k)
    _, p = scipy.stats.chisquare(observed, expected)
    assert p < 0.001, f"Concentrated chi-square p={p:.6f} not < 0.001"


def test_b6_raw_chi2_pct_match(results: dict) -> None:
    """
    Raw catalog chi-square significance percentage must be between 65% and 80%
    (expected ~71% matching A2.B6's documented chi-square rate).
    """
    pct = results["catalogs"]["raw"]["stationarity"]["pct_significant_chi2_p05"]
    assert 0.65 <= pct <= 0.80, (
        f"Raw chi2 sig pct={pct:.4f} not in [0.65, 0.80] (expected ~0.71 per A2.B6)"
    )
