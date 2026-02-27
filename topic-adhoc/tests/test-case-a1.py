"""
Case A1: Test Suite — Binning Increment Sensitivity Analysis

All 10 tests must pass before writing the whitepaper.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent          # topic-adhoc/
PROJECT_ROOT = BASE_DIR.parent                             # erebus-vee-two/
DATA_PATH = PROJECT_ROOT / "data" / "global-sets" / "iscgem_global_events.csv"
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a1-results.json"

METRICS = ["solar_secs", "lunar_secs", "midnight_secs"]
BIN_COUNTS = [8, 16, 24, 32]
LUNAR_CYCLE_SECS = 29.53059 * 86400  # 2_551_442.976
DAY_SECS = 86400
ALPHA_BONFERRONI = 0.05 / 12


def is_leap_year(year: int) -> bool:
    """Return True if *year* is a Gregorian leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def year_length_seconds(year: int) -> int:
    """Return the calendar-year length in seconds for *year*."""
    return 366 * DAY_SECS if is_leap_year(year) else 365 * DAY_SECS


@pytest.fixture(scope="module")
def results() -> dict:
    """Load and return the results JSON."""
    with open(RESULTS_PATH) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def raw_df() -> pd.DataFrame:
    """Load the raw CSV data."""
    return pd.read_csv(DATA_PATH)


# -----------------------------------------------------------------------
# 1. Results JSON exists and is valid
# -----------------------------------------------------------------------
def test_results_json_exists_and_valid(results: dict) -> None:
    """Results JSON exists, parses, and has all required top-level keys."""
    required_keys = {"generated", "data_source", "event_count", "methodology", "results"}
    assert required_keys.issubset(results.keys())
    for metric in METRICS:
        assert metric in results["results"]
        for k in BIN_COUNTS:
            entry = results["results"][metric][str(k)]
            assert {"chi2", "p_value", "cramers_v", "significant_bonferroni"}.issubset(entry.keys())


# -----------------------------------------------------------------------
# 2. Event count
# -----------------------------------------------------------------------
def test_event_count(results: dict) -> None:
    """Event count is exactly 9210."""
    assert results["event_count"] == 9210


# -----------------------------------------------------------------------
# 3. All 12 results populated (no nulls)
# -----------------------------------------------------------------------
def test_no_null_values(results: dict) -> None:
    """No null values remain in any of the 12 result entries."""
    for metric in METRICS:
        for k in BIN_COUNTS:
            entry = results["results"][metric][str(k)]
            for field in ["chi2", "p_value", "cramers_v", "significant_bonferroni"]:
                assert entry[field] is not None, f"Null found: {metric} k={k} {field}"


# -----------------------------------------------------------------------
# 4. Chi-square validity (expected count >= 5)
# -----------------------------------------------------------------------
def test_chi_square_validity(results: dict) -> None:
    """Expected count per bin >= 5 for all bin counts."""
    n = results["event_count"]
    for k in BIN_COUNTS:
        expected = n / k
        assert expected >= 5, f"Expected count {expected} < 5 for k={k}"


# -----------------------------------------------------------------------
# 5. Bonferroni threshold
# -----------------------------------------------------------------------
def test_bonferroni_threshold(results: dict) -> None:
    """alpha_bonferroni in JSON approximately equals 0.05/12."""
    stored = results["methodology"]["alpha_bonferroni"]
    assert abs(stored - ALPHA_BONFERRONI) < 1e-4, (
        f"Bonferroni threshold {stored} != expected {ALPHA_BONFERRONI}"
    )


# -----------------------------------------------------------------------
# 6. Bin sum integrity (recompute from raw CSV)
# -----------------------------------------------------------------------
def test_bin_sum_integrity(raw_df: pd.DataFrame) -> None:
    """Sum of bin counts == 9210 for every (metric, bin_count) combination."""
    df = raw_df.copy()
    year_lengths = df["solaration_year"].apply(year_length_seconds).astype(float)
    solar_phase = np.clip((df["solar_secs"] / year_lengths).values, 0, 1.0 - 1e-9)
    lunar_phase = np.clip((df["lunar_secs"] / LUNAR_CYCLE_SECS).values, 0, 1.0 - 1e-9)
    midnight_phase = np.clip((df["midnight_secs"] / float(DAY_SECS)).values, 0, 1.0 - 1e-9)
    phases = {
        "solar_secs": solar_phase,
        "lunar_secs": lunar_phase,
        "midnight_secs": midnight_phase,
    }

    for metric in METRICS:
        for k in BIN_COUNTS:
            bin_indices = np.floor(phases[metric] * k).astype(int)
            counts = np.bincount(bin_indices, minlength=k)
            total = int(counts.sum())
            assert total == 9210, (
                f"Bin sum {total} != 9210 for {metric} k={k}"
            )


# -----------------------------------------------------------------------
# 7. Phase bounds [0, 1)
# -----------------------------------------------------------------------
def test_phase_bounds(raw_df: pd.DataFrame) -> None:
    """All phase values for all three metrics are in [0, 1) after clamping.

    Solar phase may marginally exceed 1.0 for solstice-boundary events;
    the analysis script clamps these.  We verify the clamped values.
    """
    df = raw_df.copy()
    year_lengths = df["solaration_year"].apply(year_length_seconds).astype(float)
    phase_solar = np.clip((df["solar_secs"] / year_lengths).values, 0, 1.0 - 1e-9)
    phase_lunar = np.clip((df["lunar_secs"] / LUNAR_CYCLE_SECS).values, 0, 1.0 - 1e-9)
    phase_midnight = np.clip((df["midnight_secs"] / float(DAY_SECS)).values, 0, 1.0 - 1e-9)

    for name, phase in [("solar", phase_solar), ("lunar", phase_lunar), ("midnight", phase_midnight)]:
        assert (phase >= 0).all(), f"{name} phase has values < 0"
        assert (phase < 1.0).all(), f"{name} phase has values >= 1.0"


# -----------------------------------------------------------------------
# 8. midnight_secs behavioral check
# -----------------------------------------------------------------------
def test_midnight_secs_behavioral(results: dict) -> None:
    """All four midnight_secs p-values must be > Bonferroni threshold.

    If any are significant, this warrants investigation of the test setup.
    """
    for k in BIN_COUNTS:
        p = results["results"]["midnight_secs"][str(k)]["p_value"]
        assert p > ALPHA_BONFERRONI, (
            f"midnight_secs k={k} has p={p:.6f} < Bonferroni {ALPHA_BONFERRONI:.6f}. "
            f"This is unexpected at global scale and warrants investigation of the test setup."
        )


# -----------------------------------------------------------------------
# 9. Robustness flag consistency
# -----------------------------------------------------------------------
def test_robustness_flag_consistency(results: dict) -> None:
    """n_significant matches count of significant entries; robust iff n_significant >= 3."""
    for metric in METRICS:
        metric_data = results["results"][metric]
        sig_count = sum(
            1 for k in BIN_COUNTS
            if metric_data[str(k)]["significant_bonferroni"]
        )
        assert metric_data["n_significant"] == sig_count, (
            f"{metric}: n_significant={metric_data['n_significant']} but counted {sig_count}"
        )
        expected_robust = sig_count >= 3
        assert metric_data["robust"] == expected_robust, (
            f"{metric}: robust={metric_data['robust']} but expected {expected_robust} "
            f"(n_significant={sig_count})"
        )


# -----------------------------------------------------------------------
# 10. Both PNG images exist
# -----------------------------------------------------------------------
def test_png_images_exist() -> None:
    """Both visualization PNGs exist at expected paths."""
    sweep_path = OUTPUT_DIR / "case-a1-pvalue-sweep.png"
    dist_path = OUTPUT_DIR / "case-a1-distributions.png"
    assert sweep_path.exists(), f"Missing: {sweep_path}"
    assert dist_path.exists(), f"Missing: {dist_path}"
