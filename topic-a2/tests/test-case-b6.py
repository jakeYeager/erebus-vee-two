"""
Test suite for Case B6: Rolling Window Stationarity Test

All tests must pass before the whitepaper may be written.
"""

import json
import math
from pathlib import Path

import numpy as np
import pytest

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_CSV = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
RESULTS_PATH = BASE_DIR / "output" / "case-b6-results.json"

JULIAN_YEAR_SECS = 31_557_600.0
K_BINS = 24


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def catalog():
    """Load ISC-GEM catalog once for all tests."""
    import pandas as pd
    df = pd.read_csv(RAW_CSV, parse_dates=["event_at"])
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def results():
    """Load case-b6-results.json once for all tests."""
    with open(RESULTS_PATH, "r", encoding="utf-8") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def windows(results):
    """Return per-window list from results."""
    return results["windows"]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_catalog_load(catalog):
    """Assert catalog has 9210 rows and no NaT values in event_at."""
    assert len(catalog) == 9210, f"Expected 9210 rows, got {len(catalog)}"
    nat_count = catalog["event_at"].isna().sum()
    assert nat_count == 0, f"Found {nat_count} NaT values in event_at"


def test_phase_range(catalog):
    """Assert all computed phases are in [0.0, 1.0)."""
    phases = catalog["phase"].values
    assert phases.min() >= 0.0, f"Phase below 0: {phases.min()}"
    assert phases.max() < 1.0, f"Phase >= 1: {phases.max()}"


def test_window_count(windows):
    """Assert exactly 62 windows are produced."""
    assert len(windows) == 62, f"Expected 62 windows, got {len(windows)}"


def test_window_n_positive(windows):
    """Assert all windows have n > 0 and n >= 50; warn if any n < 100."""
    for w in windows:
        assert w["n"] > 0, f"Window {w['window_start']} has n=0"
        assert w["n"] >= 50, f"Window {w['window_start']} has n={w['n']} < 50"
    low_n = [w for w in windows if w["n"] < 100]
    if low_n:
        import warnings
        warnings.warn(
            f"{len(low_n)} windows have n < 100: "
            + ", ".join(str(w["window_start"]) for w in low_n),
            UserWarning,
            stacklevel=2,
        )


def test_rayleigh_uniform():
    """Rayleigh p > 0.05 in at least 95% of 100 bootstrap trials for uniform data."""
    rng = np.random.default_rng(42)
    n_trials = 100
    n_samples = 500
    pass_count = 0
    for _ in range(n_trials):
        phases = rng.uniform(0.0, 1.0, size=n_samples)
        angles = 2.0 * np.pi * phases
        mc = np.mean(np.cos(angles))
        ms = np.mean(np.sin(angles))
        R = math.sqrt(mc**2 + ms**2)
        z = n_samples * R**2
        p = math.exp(-z)
        if p > 0.05:
            pass_count += 1
    pct = pass_count / n_trials
    assert pct >= 0.95, (
        f"Rayleigh uniform test: only {pass_count}/{n_trials} ({pct:.1%}) had p>0.05; expected >=95%"
    )


def test_rayleigh_concentrated():
    """Rayleigh p < 0.001 for phases concentrated at 0.2 with Gaussian noise σ=0.02."""
    rng = np.random.default_rng(99)
    phases = (rng.normal(loc=0.2, scale=0.02, size=500)) % 1.0
    angles = 2.0 * np.pi * phases
    mc = np.mean(np.cos(angles))
    ms = np.mean(np.sin(angles))
    R = math.sqrt(mc**2 + ms**2)
    z = len(phases) * R**2
    p = math.exp(-z)
    assert p < 0.001, f"Rayleigh concentrated test: p={p:.6f} not < 0.001"


def test_rayleigh_R_bounds(windows):
    """Assert all window Rayleigh R values are in [0.0, 1.0]."""
    for w in windows:
        r = w["rayleigh_R"]
        assert 0.0 <= r <= 1.0, (
            f"Window {w['window_start']}: rayleigh_R={r} out of [0, 1]"
        )


def test_mean_phase_bounds(windows):
    """Assert all window mean_phase values are in [0.0, 1.0)."""
    for w in windows:
        mp = w["mean_phase"]
        assert 0.0 <= mp < 1.0, (
            f"Window {w['window_start']}: mean_phase={mp} out of [0, 1)"
        )


def test_1970s_flagging(windows):
    """Assert 1970s windows (start 1970-1979) are flagged True; 1980-2012 flagged False."""
    for w in windows:
        start = w["window_start"]
        if 1970 <= start <= 1979:
            assert w["is_1970s_window"] is True, (
                f"Window start {start} should be flagged as 1970s window"
            )
        elif 1980 <= start <= 2012:
            assert w["is_1970s_window"] is False, (
                f"Window start {start} should NOT be flagged as 1970s window"
            )


def test_stationarity_classification_present(results):
    """Assert 'stationarity' key is present and classification is valid."""
    assert "stationarity" in results, "Key 'stationarity' missing from results JSON"
    stat = results["stationarity"]
    assert "classification" in stat, "Key 'classification' missing from stationarity block"
    valid_classes = {"stationary", "partially stationary", "non-stationary"}
    assert stat["classification"] in valid_classes, (
        f"Invalid classification: '{stat['classification']}'; expected one of {valid_classes}"
    )


def test_chi2_k24_all_windows(windows):
    """Assert all chi2_k24 values are non-negative floats."""
    for w in windows:
        val = w["chi2_k24"]
        assert isinstance(val, (int, float)), (
            f"Window {w['window_start']}: chi2_k24 is not numeric: {type(val)}"
        )
        assert val >= 0.0, (
            f"Window {w['window_start']}: chi2_k24={val} is negative"
        )
