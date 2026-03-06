"""
Test suite for Case A3.C2: Targeted Major Sequence Phased Declustering Test.

All tests must pass before the case is marked Complete.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chisquare

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent

RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GK_MS_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv"
GK_AS_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_gk-seq_global.csv"
REAS_MS_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv"
REAS_AS_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_reas-seq_global.csv"
RESULTS_PATH = BASE_DIR / "output" / "case-a3-c2-results.json"

JULIAN_YEAR_SECS = 31_557_600.0
MAG_THRESHOLD = 8.5
K_BINS = 24

# All eight run keys: four magnitude-order + four chronological-order
ALL_RUN_KEYS = [
    "raw_gk", "raw_reas", "mainshock_gk", "mainshock_reas",
    "raw_gk_chron", "raw_reas_chron", "mainshock_gk_chron", "mainshock_reas_chron",
]
# Magnitude-order run keys only (used in some tests that pre-date chron runs)
MAG_RUN_KEYS = ["raw_gk", "raw_reas", "mainshock_gk", "mainshock_reas"]

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def raw_df() -> pd.DataFrame:
    """Load raw catalog with phase column."""
    df = pd.read_csv(RAW_PATH)
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def gk_ms_df() -> pd.DataFrame:
    """Load G-K mainshock catalog with phase column."""
    df = pd.read_csv(GK_MS_PATH)
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def gk_as_df() -> pd.DataFrame:
    """Load G-K aftershock catalog."""
    df = pd.read_csv(GK_AS_PATH)
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def reas_ms_df() -> pd.DataFrame:
    """Load Reasenberg mainshock catalog with phase column."""
    df = pd.read_csv(REAS_MS_PATH)
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def reas_as_df() -> pd.DataFrame:
    """Load Reasenberg aftershock catalog."""
    df = pd.read_csv(REAS_AS_PATH)
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def results() -> dict:
    """Load the results JSON output."""
    with open(RESULTS_PATH) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_catalog_loads(raw_df, gk_ms_df, gk_as_df, reas_ms_df, reas_as_df) -> None:
    """Assert expected row counts, timestamp parsing, and phase range for all catalogs."""
    assert len(raw_df) == 9210, f"raw: expected 9210, got {len(raw_df)}"
    assert len(gk_ms_df) == 5883, f"gk_ms: expected 5883, got {len(gk_ms_df)}"
    assert len(gk_as_df) == 3327, f"gk_as: expected 3327, got {len(gk_as_df)}"
    assert len(reas_ms_df) == 8265, f"reas_ms: expected 8265, got {len(reas_ms_df)}"
    assert len(reas_as_df) == 945, f"reas_as: expected 945, got {len(reas_as_df)}"

    for name, df in [
        ("raw", raw_df),
        ("gk_ms", gk_ms_df),
        ("gk_as", gk_as_df),
        ("reas_ms", reas_ms_df),
        ("reas_as", reas_as_df),
    ]:
        assert df["event_at"].isna().sum() == 0, f"{name}: NaT values found in event_at"
        assert (df["phase"] >= 0.0).all(), f"{name}: phase < 0.0 found"
        assert (df["phase"] < 1.0).all(), f"{name}: phase >= 1.0 found"


def test_major_event_count(results) -> None:
    """Assert exactly 12 M>=8.5 events identified after mainshock filter.

    The 1960-05-22 M8.6 Valdivia-2 event (iscgem879134) must be excluded because
    it is classified as aftershock in both G-K and Reasenberg. removal_order_magnitude
    must be sorted by usgs_mag descending with tie-breaking by event_at ascending.
    """
    major = results["major_events"]

    # Exactly 12 events after mainshock filter
    assert len(major) == 12, f"Expected exactly 12 major events, got {len(major)}"

    # All events must be above the magnitude threshold
    for ev in major:
        assert ev["usgs_mag"] >= MAG_THRESHOLD, (
            f"Event {ev['usgs_id']} has mag {ev['usgs_mag']} < {MAG_THRESHOLD}"
        )

    # Valdivia-2 (iscgem879134) must be excluded
    ids = [ev["usgs_id"] for ev in major]
    assert "iscgem879134" not in ids, (
        "iscgem879134 (Valdivia-2, aftershock in both algorithms) must be excluded"
    )

    # removal_order_magnitude must reflect magnitude-descending sort with tie-break by event_at
    ordered_by_mag_rank = sorted(major, key=lambda e: e["removal_order_magnitude"])
    mags = [ev["usgs_mag"] for ev in ordered_by_mag_rank]
    for i in range(len(mags) - 1):
        assert mags[i] >= mags[i + 1], (
            f"removal_order_magnitude not descending: position {i} mag={mags[i]} "
            f"< position {i+1} mag={mags[i+1]}"
        )
    # Check tie-breaking: events with same magnitude should be ordered by event_at ascending
    for i in range(len(ordered_by_mag_rank) - 1):
        ev_a = ordered_by_mag_rank[i]
        ev_b = ordered_by_mag_rank[i + 1]
        if ev_a["usgs_mag"] == ev_b["usgs_mag"]:
            assert ev_a["event_at"] <= ev_b["event_at"], (
                f"Tie in magnitude {ev_a['usgs_mag']}: event {ev_a['usgs_id']} at "
                f"{ev_a['event_at']} should precede {ev_b['usgs_id']} at {ev_b['event_at']}"
            )


def test_removal_monotonic(results) -> None:
    """Assert n_remaining is non-increasing across steps; step 0 equals full catalog size."""
    catalog_sizes = {
        "raw_gk": 9210,
        "raw_reas": 9210,
        "mainshock_gk": 5883,
        "mainshock_reas": 8265,
        "raw_gk_chron": 9210,
        "raw_reas_chron": 9210,
        "mainshock_gk_chron": 5883,
        "mainshock_reas_chron": 8265,
    }
    for run_key in ALL_RUN_KEYS:
        steps = results["runs"][run_key]["steps"]
        n_vals = [s["n_remaining"] for s in steps]

        # Step 0 must equal full catalog size
        assert n_vals[0] == catalog_sizes[run_key], (
            f"{run_key}: step 0 n_remaining={n_vals[0]}, expected {catalog_sizes[run_key]}"
        )

        # Non-increasing
        for i in range(len(n_vals) - 1):
            assert n_vals[i] >= n_vals[i + 1], (
                f"{run_key}: n_remaining increased from step {i} ({n_vals[i]}) "
                f"to step {i+1} ({n_vals[i+1]})"
            )


def test_step0_matches_full_catalog_stats(
    raw_df, gk_ms_df, reas_ms_df, results
) -> None:
    """Assert step 0 chi2 and Cramér's V match direct computation on unmodified base catalogs."""
    base_dfs = {
        "raw_gk": raw_df,
        "raw_reas": raw_df,
        "mainshock_gk": gk_ms_df,
        "mainshock_reas": reas_ms_df,
        "raw_gk_chron": raw_df,
        "raw_reas_chron": raw_df,
        "mainshock_gk_chron": gk_ms_df,
        "mainshock_reas_chron": reas_ms_df,
    }

    tol = 1e-6
    for run_key in ALL_RUN_KEYS:
        df = base_dfs[run_key]
        phases = df["phase"].values
        n = len(phases)
        bin_indices = np.floor(phases * K_BINS).astype(int) % K_BINS
        observed = np.bincount(bin_indices, minlength=K_BINS)
        expected = np.full(K_BINS, n / K_BINS)
        chi2_stat, _ = chisquare(observed, expected)
        cramers_v = float(np.sqrt(chi2_stat / (n * (K_BINS - 1))))

        step0 = results["runs"][run_key]["steps"][0]
        assert abs(step0["stats"]["chi2_k24"] - chi2_stat) < tol, (
            f"{run_key}: step0 chi2 mismatch: {step0['stats']['chi2_k24']} vs {chi2_stat}"
        )
        assert abs(step0["stats"]["cramers_v"] - cramers_v) < tol, (
            f"{run_key}: step0 Cramér's V mismatch: {step0['stats']['cramers_v']} vs {cramers_v}"
        )


def test_chi2_bounds(results) -> None:
    """Assert chi2 >= 0 and p-value in [0.0, 1.0] for all steps and runs."""
    for run_key in ALL_RUN_KEYS:
        for step in results["runs"][run_key]["steps"]:
            stats = step["stats"]
            chi2 = stats["chi2_k24"]
            p = stats["p_chi2_k24"]
            if chi2 is not None:
                assert chi2 >= 0, f"{run_key} step {step['step']}: chi2={chi2} < 0"
            if p is not None:
                assert 0.0 <= p <= 1.0, (
                    f"{run_key} step {step['step']}: p_chi2={p} out of [0,1]"
                )


def test_cramers_v_bounds(results) -> None:
    """Assert Cramér's V >= 0 for all steps and runs."""
    for run_key in ALL_RUN_KEYS:
        for step in results["runs"][run_key]["steps"]:
            v = step["stats"]["cramers_v"]
            if v is not None:
                assert v >= 0, (
                    f"{run_key} step {step['step']}: cramers_v={v} < 0"
                )


def test_interval_z_type(results) -> None:
    """Assert interval z-score fields are float-convertible for all steps."""
    for run_key in ALL_RUN_KEYS:
        for step in results["runs"][run_key]["steps"]:
            stats = step["stats"]
            for interval_key in ("interval_1_z", "interval_2_z", "interval_3_z"):
                val = stats[interval_key]
                assert isinstance(val, (int, float)), (
                    f"{run_key} step {step['step']} {interval_key}: type {type(val)} is not numeric"
                )


def test_parent_id_removal(raw_df, gk_as_df, results) -> None:
    """Assert that step 1 of raw_gk removes the major event and its G-K aftershocks."""
    step1 = results["runs"]["raw_gk"]["steps"][1]
    removed_event_id = step1["event_removed"]["usgs_id"]

    # Get G-K aftershocks for this event
    attributed_ids = set(gk_as_df[gk_as_df["parent_id"] == removed_event_id]["usgs_id"].tolist())
    removal_set = {removed_event_id} | attributed_ids

    remaining_df = raw_df[~raw_df["usgs_id"].isin(removal_set)]
    n_remaining = step1["n_remaining"]
    assert len(remaining_df) == n_remaining, (
        f"raw_gk step 1: expected {n_remaining} remaining, got {len(remaining_df)}"
    )

    # Assert removed event is absent
    assert removed_event_id not in remaining_df["usgs_id"].values, (
        f"Major event {removed_event_id} still present after step 1 removal"
    )

    # Assert attributed aftershocks are absent
    for as_id in attributed_ids:
        assert as_id not in remaining_df["usgs_id"].values, (
            f"G-K aftershock {as_id} of {removed_event_id} still present after step 1 removal"
        )


def test_mainshock_only_removal(gk_ms_df, results) -> None:
    """Assert that mainshock_gk step 1 removes exactly one event from the mainshock catalog.

    The first major event in magnitude order (M9.55, 1960 Valdivia) is present as a mainshock
    in G-K, so n_remaining for mainshock_gk step 1 should equal n_gk_mainshocks - 1.
    """
    step1 = results["runs"]["mainshock_gk"]["steps"][1]
    expected_n = len(gk_ms_df) - 1
    assert step1["n_remaining"] == expected_n, (
        f"mainshock_gk step 1: expected n_remaining={expected_n}, got {step1['n_remaining']}"
    )


def test_sequence_metrics_structure(results) -> None:
    """Assert sequence_metrics contains 'gk' and 'reasenberg' lists with correct structure."""
    seq_metrics = results["sequence_metrics"]
    assert "gk" in seq_metrics, "sequence_metrics missing 'gk' key"
    assert "reasenberg" in seq_metrics, "sequence_metrics missing 'reasenberg' key"

    n_major = len(results["major_events"])
    assert len(seq_metrics["gk"]) == n_major, (
        f"gk sequence_metrics length {len(seq_metrics['gk'])} != n_major {n_major}"
    )
    assert len(seq_metrics["reasenberg"]) == n_major, (
        f"reasenberg sequence_metrics length {len(seq_metrics['reasenberg'])} != n_major {n_major}"
    )

    required_keys = {
        "foreshock_count", "aftershock_count", "total_sequence",
        "window_days", "early_count", "late_count",
    }
    for label in ("gk", "reasenberg"):
        for entry in seq_metrics[label]:
            missing = required_keys - set(entry.keys())
            assert not missing, (
                f"sequence_metrics[{label}] entry {entry.get('usgs_id', '?')} "
                f"missing keys: {missing}"
            )


def test_chronological_order(results) -> None:
    """Assert removal_order_chronological is present and events are sorted by event_at ascending.

    Also asserts iscgem879134 is absent from major_events.
    """
    major = results["major_events"]

    # iscgem879134 must be absent (excluded as aftershock in both algorithms)
    ids = [ev["usgs_id"] for ev in major]
    assert "iscgem879134" not in ids, (
        "iscgem879134 must be absent from major_events (aftershock in both algorithms)"
    )

    # All events must have removal_order_chronological field
    for ev in major:
        assert "removal_order_chronological" in ev, (
            f"Event {ev['usgs_id']} missing removal_order_chronological"
        )

    # Events ordered by removal_order_chronological must be sorted by event_at ascending
    ordered_chron = sorted(major, key=lambda e: e["removal_order_chronological"])
    for i in range(len(ordered_chron) - 1):
        assert ordered_chron[i]["event_at"] <= ordered_chron[i + 1]["event_at"], (
            f"Chronological order violated at position {i}: "
            f"{ordered_chron[i]['event_at']} > {ordered_chron[i+1]['event_at']}"
        )


def test_results_json_completeness(results) -> None:
    """Assert all eight run keys are present and each has the correct number of steps."""
    n_major = len(results["major_events"])
    expected_steps = n_major + 1  # baseline + one per removal

    for run_key in ALL_RUN_KEYS:
        assert run_key in results["runs"], f"Run key '{run_key}' missing from results"
        actual_steps = len(results["runs"][run_key]["steps"])
        assert actual_steps == expected_steps, (
            f"{run_key}: expected {expected_steps} steps, got {actual_steps}"
        )


def test_raw_baseline_chi2(results) -> None:
    """Assert raw catalog step 0 chi2 approximately matches the A2.A4 global value (~69.37).

    Tolerance: ±2.0 (Julian year constant vs. per-year normalization differences).
    """
    step0_chi2 = results["runs"]["raw_gk"]["steps"][0]["stats"]["chi2_k24"]
    expected = 69.37
    tolerance = 2.0
    assert abs(step0_chi2 - expected) <= tolerance, (
        f"Raw baseline chi2={step0_chi2:.4f} deviates more than ±{tolerance} "
        f"from expected {expected}"
    )
