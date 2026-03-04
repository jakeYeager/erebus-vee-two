"""
Test suite for Case A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation

All 14 tests must pass. If test_midcrustal_m6_significant fails, it is a scientifically
significant null result and must be documented in the whitepaper, not silenced.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

logger = logging.getLogger(__name__)

BASE_DIR   = Path(__file__).resolve().parent.parent
LIB_DIR    = BASE_DIR.parent / "lib"
DATA_DIR   = BASE_DIR.parent / "data" / "iscgem"
OUTPUT_DIR = BASE_DIR / "output"

RESULTS_PATH  = OUTPUT_DIR / "case-a3-b4-results.json"
EVENTS_PKL    = OUTPUT_DIR / "case-a3-b4-events.pkl"
RAW_PATH      = DATA_DIR / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH    = DATA_DIR / "plate-location" / "ocean_class_gshhg_global.csv"
CRUST1_PATH   = LIB_DIR / "crust1.bnds"
STEPS_PATH    = LIB_DIR / "PB2002_steps.dat"

JULIAN_YEAR_SECS = 31_557_600.0

DEPTH_BANDS = [
    {"label": "shallow_0-20km",        "min": 0,   "max": 20},
    {"label": "midcrustal_20-70km",    "min": 20,  "max": 70},
    {"label": "intermediate_70-300km", "min": 70,  "max": 300},
    {"label": "deep_300km+",           "min": 300, "max": 9999},
]
MAG_BANDS = [
    {"label": "m6_6.9",  "min": 6.0, "max": 7.0},
    {"label": "m7_7.9",  "min": 7.0, "max": 8.0},
    {"label": "m8_plus", "min": 8.0, "max": 99.0},
]


@pytest.fixture(scope="module")
def results() -> dict:
    """Load the A3.B4 results JSON."""
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def events_df() -> pd.DataFrame:
    """Load the merged events dataframe."""
    return pd.read_pickle(EVENTS_PKL)


# ---------------------------------------------------------------------------
# Test 1: Catalog load
# ---------------------------------------------------------------------------

def test_catalog_load(events_df: pd.DataFrame) -> None:
    """Assert n=9210; assert phase in [0, 1); log any negative depths."""
    assert len(events_df) == 9210, f"Expected 9210 rows, got {len(events_df)}"

    phases = events_df["phase"].values
    assert np.all(phases >= 0.0), "Phase values below 0.0 found"
    assert np.all(phases < 1.0), "Phase values >= 1.0 found"

    n_neg_depth = int((events_df["depth"] < 0).sum())
    logger.info(f"Events with depth < 0: {n_neg_depth}")


# ---------------------------------------------------------------------------
# Test 2: Depth band partition
# ---------------------------------------------------------------------------

def test_depth_band_partition(results: dict, events_df: pd.DataFrame) -> None:
    """Assert sum of all band sizes + n_depth_null == 9210; assert no overlap."""
    n_depth_null = results["parameters"]["n_depth_null"]

    # Count events in each band from the dataframe
    df_valid = events_df[events_df["depth"].notna()].copy()
    total_in_bands = 0
    prev_sets: list[set] = []

    for db in DEPTH_BANDS:
        d_min, d_max = db["min"], db["max"]
        mask = (df_valid["depth"] >= d_min) & (df_valid["depth"] < d_max)
        band_ids = set(df_valid.loc[mask, "usgs_id"].tolist())

        # No overlap with previous bands
        for prev in prev_sets:
            assert len(band_ids & prev) == 0, (
                f"Overlap detected in depth band {db['label']}"
            )
        prev_sets.append(band_ids)
        total_in_bands += len(band_ids)

    assert total_in_bands + n_depth_null == 9210, (
        f"Band total {total_in_bands} + null {n_depth_null} != 9210"
    )


# ---------------------------------------------------------------------------
# Test 3: Magnitude band partition
# ---------------------------------------------------------------------------

def test_magnitude_band_partition(results: dict, events_df: pd.DataFrame) -> None:
    """For each depth band, assert sum across 3 mag bands equals depth band total."""
    sm = results["stratification_matrix"]
    db_totals = results["depth_band_totals"]

    for db in DEPTH_BANDS:
        d_label = db["label"]
        total_n = db_totals[d_label]["n"]

        mag_sum = sum(sm[d_label][mb["label"]]["n"] for mb in MAG_BANDS)
        assert mag_sum == total_n, (
            f"Depth band {d_label}: mag sum {mag_sum} != total {total_n}"
        )


# ---------------------------------------------------------------------------
# Test 4: Adaptive k assignment
# ---------------------------------------------------------------------------

def test_adaptive_k(results: dict) -> None:
    """For each cell, assert k matches adaptive rule: k=24 if n≥500, k=16 if 200≤n<500, k=12 otherwise."""
    sm = results["stratification_matrix"]

    for db in DEPTH_BANDS:
        d_label = db["label"]
        for mb in MAG_BANDS:
            m_label = mb["label"]
            cell = sm[d_label][m_label]
            n = cell["n"]
            k = cell["k"]

            if n >= 500:
                expected_k = 24
            elif n >= 200:
                expected_k = 16
            else:
                expected_k = 12

            assert k == expected_k, (
                f"Cell {d_label}×{m_label}: n={n}, expected k={expected_k}, got k={k}"
            )


# ---------------------------------------------------------------------------
# Test 5: Chi-square bounds
# ---------------------------------------------------------------------------

def test_chi2_bounds(results: dict) -> None:
    """For all cells in stratification_matrix, assert chi2 >= 0 and p_chi2 in [0, 1]."""
    sm = results["stratification_matrix"]

    for db in DEPTH_BANDS:
        d_label = db["label"]
        for mb in MAG_BANDS:
            m_label = mb["label"]
            cell = sm[d_label][m_label]
            assert cell["chi2"] >= 0.0, f"chi2 < 0 in {d_label}×{m_label}"
            assert 0.0 <= cell["p_chi2"] <= 1.0, f"p_chi2 out of [0,1] in {d_label}×{m_label}"


# ---------------------------------------------------------------------------
# Test 6: Cramér's V bounds
# ---------------------------------------------------------------------------

def test_cramer_v_bounds(results: dict) -> None:
    """For all cells, assert cramers_v >= 0."""
    sm = results["stratification_matrix"]

    for db in DEPTH_BANDS:
        d_label = db["label"]
        for mb in MAG_BANDS:
            m_label = mb["label"]
            cell = sm[d_label][m_label]
            assert cell["cramers_v"] >= 0.0, (
                f"cramers_v < 0 in {d_label}×{m_label}: {cell['cramers_v']}"
            )


# ---------------------------------------------------------------------------
# Test 7: CRUST1.0 load
# ---------------------------------------------------------------------------

def test_crust1_load() -> None:
    """Assert shape (180, 360, 9); assert Moho layer 8 all negative; assert range within [-90, 0]."""
    raw = np.loadtxt(CRUST1_PATH)
    assert raw.shape == (64800, 9), f"Unexpected raw shape: {raw.shape}"

    bnds = raw.reshape(180, 360, 9)
    assert bnds.shape == (180, 360, 9), f"Reshape failed: {bnds.shape}"

    moho_layer = bnds[:, :, 8]  # Should all be negative (depths below sea level)
    assert np.all(moho_layer <= 0.0), (
        f"Moho layer has positive values; min={moho_layer.min():.2f}, max={moho_layer.max():.2f}"
    )
    assert moho_layer.min() >= -90.0, f"Moho depth < -90 km: {moho_layer.min():.2f}"
    assert moho_layer.max() <= 0.0, f"Moho depth > 0 km: {moho_layer.max():.2f}"


# ---------------------------------------------------------------------------
# Test 8: Moho assignment range
# ---------------------------------------------------------------------------

def test_moho_assignment_range(events_df: pd.DataFrame) -> None:
    """Assert all moho_depth_km values in [3.0, 90.0] km."""
    moho_depths = events_df["moho_depth_km"].values
    assert np.all(moho_depths >= 3.0), f"Moho depth below 3.0 km found: min={moho_depths.min():.2f}"
    assert np.all(moho_depths <= 90.0), f"Moho depth above 90.0 km found: max={moho_depths.max():.2f}"


# ---------------------------------------------------------------------------
# Test 9: Moho proximal counts monotonicity
# ---------------------------------------------------------------------------

def test_moho_proximal_counts(results: dict) -> None:
    """Assert n_proximal at delta=10 > 0 and < 9210; assert counts increase with delta."""
    n_by_delta = results["moho_assignment"]["n_proximal_by_delta"]
    n5  = n_by_delta["5"]
    n10 = n_by_delta["10"]
    n15 = n_by_delta["15"]

    assert n10 > 0, f"n_proximal at delta=10 km is 0"
    assert n10 < 9210, f"n_proximal at delta=10 km equals total catalog ({n10})"
    assert n5 <= n10, f"n_proximal not monotone: delta=5 ({n5}) > delta=10 ({n10})"
    assert n10 <= n15, f"n_proximal not monotone: delta=10 ({n10}) > delta=15 ({n15})"


# ---------------------------------------------------------------------------
# Test 10: Sub boundary points count
# ---------------------------------------------------------------------------

def test_sub_boundary_points(results: dict) -> None:
    """Assert n_sub_boundary_points between 700 and 1100 (matching B3 test)."""
    n_pts = results["subduction_proximity"]["n_sub_boundary_points"]
    assert 700 <= n_pts <= 1100, (
        f"n_sub_boundary_points={n_pts} outside expected range [700, 1100]"
    )


# ---------------------------------------------------------------------------
# Test 11: Sub distances non-negative
# ---------------------------------------------------------------------------

def test_sub_distances_nonneg(events_df: pd.DataFrame) -> None:
    """Assert all dist_to_subduction_km >= 0."""
    dists = events_df["dist_to_subduction_km"].values
    assert np.all(dists >= 0.0), f"Negative subduction distances found: min={dists.min():.2f}"


# ---------------------------------------------------------------------------
# Test 12: A2.B4 regression
# ---------------------------------------------------------------------------

def test_a2b4_regression(results: dict) -> None:
    """
    Assert mid-crustal depth-band total chi2 is within ±2.0 of A2.B4 reference (85.48).
    Assert mid-crustal p_chi2 < 1e-6.
    """
    mid = results["depth_band_totals"]["midcrustal_20-70km"]
    chi2 = mid["chi2"]
    p    = mid["p_chi2"]

    assert abs(chi2 - 85.48) <= 2.0, (
        f"Mid-crustal chi2={chi2:.2f} deviates from A2.B4 reference 85.48 by more than ±2.0"
    )
    assert p < 1e-6, f"Mid-crustal p_chi2={p:.4e} not < 1e-6"


# ---------------------------------------------------------------------------
# Test 13: Mid-crustal M6 significance (scientific assertion)
# ---------------------------------------------------------------------------

def test_midcrustal_m6_significant(results: dict) -> None:
    """
    Assert mid-crustal × M6–7 cell p_chi2 < 0.05.

    This is a scientific assertion: if it fails, the mid-crustal signal does NOT
    persist in M6–7 alone, which is a scientifically significant null result.
    A failure here must be documented in the whitepaper, not silenced.
    """
    cell = results["stratification_matrix"]["midcrustal_20-70km"]["m6_6.9"]
    p = cell["p_chi2"]
    n = cell["n"]
    assert p < 0.05, (
        f"SCIENTIFIC NULL RESULT: Mid-crustal × M6–7 p_chi2={p:.4e} >= 0.05 (n={n}). "
        "The mid-crustal solar signal does NOT persist when magnitude-controlled to M6–7. "
        "This must be reported in the whitepaper as a null result."
    )


# ---------------------------------------------------------------------------
# Test 14: Results JSON completeness
# ---------------------------------------------------------------------------

def test_results_json_completeness(results: dict) -> None:
    """Assert all expected keys and sub-keys are present in results JSON."""
    sm = results["stratification_matrix"]
    expected_depth_labels = [
        "shallow_0-20km", "midcrustal_20-70km",
        "intermediate_70-300km", "deep_300km+",
    ]
    expected_mag_labels = ["m6_6.9", "m7_7.9", "m8_plus"]

    for d_label in expected_depth_labels:
        assert d_label in sm, f"Missing depth band in stratification_matrix: {d_label}"
        for m_label in expected_mag_labels:
            assert m_label in sm[d_label], (
                f"Missing mag band {m_label} in stratification_matrix[{d_label}]"
            )

    mi = results["moho_isolation"]
    for dk in ["delta_5km", "delta_10km", "delta_15km"]:
        assert dk in mi, f"Missing key in moho_isolation: {dk}"

    ct = results["subduction_crosstab_midcrustal"]
    assert "all_magnitudes" in ct, "Missing all_magnitudes in subduction_crosstab_midcrustal"
    for m_label in expected_mag_labels:
        assert m_label in ct, f"Missing {m_label} in subduction_crosstab_midcrustal"
    for key in ct:
        assert "near_sub" in ct[key], f"Missing near_sub in crosstab key {key}"
        assert "far_sub"  in ct[key], f"Missing far_sub in crosstab key {key}"
