"""
Test Suite: Case A1b — Elevated Bin Event Characterization and Declustering Implications

All 10 tests must pass.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent          # topic-adhoc/
PROJECT_ROOT = BASE_DIR.parent                             # erebus-vee-two/
DATA_PATH = PROJECT_ROOT / "data" / "global-sets" / "iscgem_global_events.csv"
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a1b-results.json"

# Phase normalization constants (matching analysis script)
DAY_SECS = 86400
BIN_COUNTS = [16, 24, 32]
TOP_N = 3


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def results() -> dict:
    """Load the A1b results JSON once for all tests."""
    with open(RESULTS_PATH) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def df() -> pd.DataFrame:
    """Load ISC-GEM catalog once for all tests."""
    return pd.read_csv(DATA_PATH)


@pytest.fixture(scope="module")
def elevated_df(df: pd.DataFrame, results: dict) -> pd.DataFrame:
    """Recompute the elevated-bin events from scratch for de-duplication tests."""
    phases = _compute_solar_phase(df)
    combined_intervals = results["phase_coherence"]["combined_elevated_intervals"]
    n = len(df)
    elevated_mask = np.zeros(n, dtype=bool)
    for lo, hi in combined_intervals:
        elevated_mask |= (phases >= lo) & (phases < hi)
    elev_df = df[elevated_mask].drop_duplicates(subset=["usgs_id"])
    return elev_df


# ---------------------------------------------------------------------------
# Utility: phase normalization (mirrors analysis script)
# ---------------------------------------------------------------------------

def _is_leap_year(year: int) -> bool:
    """Return True if year is a Gregorian leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def _year_length_seconds(year: int) -> int:
    """Return calendar-year length in seconds."""
    return 366 * DAY_SECS if _is_leap_year(year) else 365 * DAY_SECS


def _compute_solar_phase(df: pd.DataFrame) -> np.ndarray:
    """Compute solar phase in [0, 1) matching A1b analysis."""
    year_lengths = df["solaration_year"].apply(_year_length_seconds).astype(float)
    phase = df["solar_secs"] / year_lengths
    over_mask = phase >= 1.0
    if over_mask.any():
        phase = phase.copy()
        phase[over_mask] = 0.999999
    return phase.values


def _merge_intervals(intervals: list) -> list:
    """Merge overlapping or adjacent intervals."""
    if not intervals:
        return []
    sorted_ivs = sorted(intervals)
    merged = [list(sorted_ivs[0])]
    for lo, hi in sorted_ivs[1:]:
        if lo <= merged[-1][1] + 1e-12:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [tuple(m) for m in merged]


# ---------------------------------------------------------------------------
# Test 1: Results JSON exists and is valid; all top-level keys present; no nulls
# ---------------------------------------------------------------------------

EXPECTED_TOP_LEVEL_KEYS = {
    "generated", "data_source", "boundary_file", "event_count",
    "phase_coherence", "temporal", "spatial", "magnitude",
    "declustering_window_estimate",
}


def _check_no_nulls(obj, path: str = "") -> list:
    """Recursively find all null values in a nested dict/list.

    Returns list of paths to null values.
    """
    nulls = []
    if obj is None:
        nulls.append(path)
    elif isinstance(obj, dict):
        for k, v in obj.items():
            if k.startswith("_"):
                continue  # skip internal/private fields
            nulls.extend(_check_no_nulls(v, f"{path}.{k}" if path else str(k)))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            nulls.extend(_check_no_nulls(v, f"{path}[{i}]"))
    return nulls


def test_results_json_valid(results: dict) -> None:
    """Test 1: Results JSON exists, all top-level keys present, no null values."""
    assert RESULTS_PATH.exists(), f"Results JSON not found at {RESULTS_PATH}"

    for key in EXPECTED_TOP_LEVEL_KEYS:
        assert key in results, f"Missing top-level key: {key}"

    null_paths = _check_no_nulls(results)
    assert len(null_paths) == 0, f"Null values found at: {null_paths}"


# ---------------------------------------------------------------------------
# Test 2: Phase coherence — event count
# ---------------------------------------------------------------------------

def test_n_elevated_positive(results: dict) -> None:
    """Test 2: n_elevated is a positive integer < 9210."""
    n_elevated = results["phase_coherence"]["n_elevated"]
    assert isinstance(n_elevated, int), f"n_elevated is not int: {type(n_elevated)}"
    assert 0 < n_elevated < 9210, f"n_elevated out of expected range: {n_elevated}"


# ---------------------------------------------------------------------------
# Test 3: Phase coherence — de-duplication
# ---------------------------------------------------------------------------

def test_deduplication(elevated_df: pd.DataFrame, results: dict) -> None:
    """Test 3: Recomputed elevated count matches n_elevated; no duplicate usgs_ids."""
    n_elevated = results["phase_coherence"]["n_elevated"]
    assert len(elevated_df) == n_elevated, (
        f"Recomputed elevated count {len(elevated_df)} != stored n_elevated {n_elevated}"
    )
    assert elevated_df["usgs_id"].nunique() == len(elevated_df), (
        "Duplicate usgs_id values found in elevated-bin events"
    )


# ---------------------------------------------------------------------------
# Test 4: Phase coherence — combined interval coverage
# ---------------------------------------------------------------------------

def test_expected_pct_under_null(results: dict) -> None:
    """Test 4: expected_pct_under_null equals total fraction of [0,1) covered
    by the combined elevated intervals, within floating point tolerance."""
    combined_intervals = results["phase_coherence"]["combined_elevated_intervals"]
    stored_pct = results["phase_coherence"]["expected_pct_under_null"]

    total_coverage = sum(hi - lo for lo, hi in combined_intervals)
    expected_pct_recomputed = total_coverage * 100

    assert abs(stored_pct - expected_pct_recomputed) < 0.01, (
        f"expected_pct_under_null {stored_pct:.4f} != "
        f"recomputed {expected_pct_recomputed:.4f} (diff={abs(stored_pct - expected_pct_recomputed):.6f})"
    )


# ---------------------------------------------------------------------------
# Test 5: Temporal — IEI sum
# ---------------------------------------------------------------------------

def test_iei_count(results: dict) -> None:
    """Test 5: Number of IEI values equals n_elevated - 1."""
    n_elevated = results["phase_coherence"]["n_elevated"]
    # The analysis script stores _iei_count in the JSON
    iei_count = results.get("_iei_count")
    assert iei_count is not None, "Missing _iei_count in results JSON"
    assert iei_count == n_elevated - 1, (
        f"IEI count {iei_count} != n_elevated - 1 = {n_elevated - 1}"
    )


# ---------------------------------------------------------------------------
# Test 6: Spatial — NN count
# ---------------------------------------------------------------------------

def test_nn_count(results: dict) -> None:
    """Test 6: Number of nearest-neighbor distances equals n_elevated."""
    n_elevated = results["phase_coherence"]["n_elevated"]
    nn_count = results.get("_nn_count")
    assert nn_count is not None, "Missing _nn_count in results JSON"
    assert nn_count == n_elevated, (
        f"NN count {nn_count} != n_elevated {n_elevated}"
    )


# ---------------------------------------------------------------------------
# Test 7: Spatial — boundary proximity sums
# ---------------------------------------------------------------------------

def test_boundary_proximity_sums(results: dict) -> None:
    """Test 7: For both elevated and full catalog, proximity percentages sum to ~100."""
    bnd = results["spatial"]["boundary_proximity"]
    tolerance = 0.1

    for pop_name in ["elevated", "full_catalog"]:
        pop = bnd[pop_name]
        total = (pop["near_boundary_pct"] +
                 pop["transitional_pct"] +
                 pop["intraplate_pct"])
        assert abs(total - 100.0) < tolerance, (
            f"Boundary proximity sum for {pop_name}: {total:.4f} (expected ~100.0, "
            f"tolerance={tolerance})"
        )


# ---------------------------------------------------------------------------
# Test 8: Magnitude sums
# ---------------------------------------------------------------------------

def test_magnitude_sums(results: dict) -> None:
    """Test 8: Sum of magnitude band counts equals respective population size."""
    n_elevated = results["phase_coherence"]["n_elevated"]
    n_full = results["event_count"]

    mag = results["magnitude"]
    elev_sum = sum(mag["elevated"].values())
    full_sum = sum(mag["full_catalog"].values())

    assert elev_sum == n_elevated, (
        f"Elevated magnitude bin sum {elev_sum} != n_elevated {n_elevated}"
    )
    assert full_sum == n_full, (
        f"Full catalog magnitude bin sum {full_sum} != event_count {n_full}"
    )


# ---------------------------------------------------------------------------
# Test 9: Declustering proposal populated
# ---------------------------------------------------------------------------

def test_declustering_proposal_populated(results: dict) -> None:
    """Test 9: proposed_window spatial_km and temporal_days are positive numbers
    and basis is a non-empty string."""
    pw = results["declustering_window_estimate"]["proposed_window"]

    spatial_km = pw["spatial_km"]
    temporal_days = pw["temporal_days"]
    basis = pw["basis"]

    assert isinstance(spatial_km, (int, float)) and spatial_km > 0, (
        f"proposed_window.spatial_km must be positive number, got: {spatial_km}"
    )
    assert isinstance(temporal_days, (int, float)) and temporal_days > 0, (
        f"proposed_window.temporal_days must be positive number, got: {temporal_days}"
    )
    assert isinstance(basis, str) and len(basis.strip()) > 0, (
        "proposed_window.basis must be a non-empty string"
    )


# ---------------------------------------------------------------------------
# Test 10: All three PNG images exist
# ---------------------------------------------------------------------------

def test_png_images_exist() -> None:
    """Test 10: All three expected PNG files exist at output paths."""
    expected_images = [
        OUTPUT_DIR / "case-a1b-phase-coherence.png",
        OUTPUT_DIR / "case-a1b-temporal-magnitude.png",
        OUTPUT_DIR / "case-a1b-spatial.png",
    ]
    for img_path in expected_images:
        assert img_path.exists(), f"Expected PNG image not found: {img_path}"
