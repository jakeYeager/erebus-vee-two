"""
Case A2: b-Value Seasonal Variation — Test Suite.

Tests:
    test_catalog_load               — Row count 9210; no NaN in usgs_mag
    test_phase_range                — All phases in [0.0, 1.0)
    test_bin_partition_k24          — Sum of per-bin counts at k=24 == 9210
    test_bin_partition_k32          — Sum of per-bin counts at k=32 == 9210
    test_b_mle_formula              — Synthetic b-value check ≈ 2.172
    test_b_mle_range                — All computed b_mle in [0.5, 3.0]
    test_bootstrap_ci_contains_mle  — CI lower <= b_mle <= CI upper for each bin
    test_low_n_flag_trigger         — low_n_flag matches n < 20 criterion
    test_anova_p_nonnegative        — p_anova in [0.0, 1.0] for both k values
    test_pearson_r_range            — r_rate_b in [-1.0, 1.0]
    test_phase_class_valid          — b_max and b_min phase classes are valid strings
    test_results_json_keys          — JSON has required top-level and sub-keys
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
RESULTS_PATH = BASE_DIR / "output" / "case-a2-results.json"

JULIAN_YEAR_SECS: float = 31_557_600.0
MC: float = 6.0
EXPECTED_N: int = 9210
VALID_PHASE_CLASSES = {"near_solstice", "near_equinox", "other"}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def df() -> pd.DataFrame:
    """Load and return the raw ISC-GEM catalog."""
    return pd.read_csv(RAW_PATH)


@pytest.fixture(scope="module")
def results() -> dict:
    """Load and return case-a2-results.json."""
    with open(RESULTS_PATH, "r") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_catalog_load(df: pd.DataFrame) -> None:
    """Assert catalog has 9210 rows and no NaN in usgs_mag."""
    assert len(df) == EXPECTED_N, f"Expected {EXPECTED_N} rows, got {len(df)}"
    assert df["usgs_mag"].notna().all(), "usgs_mag contains NaN values"


def test_phase_range(df: pd.DataFrame) -> None:
    """Assert all computed solar phases are in [0.0, 1.0)."""
    phases = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    assert (phases >= 0.0).all(), "Some phases are < 0.0"
    assert (phases < 1.0).all(), "Some phases are >= 1.0"


def test_bin_partition_k24(results: dict) -> None:
    """Assert sum of per-bin counts at k=24 equals 9210."""
    bins = results["phase_variation"]["k24"]["bins"]
    total = sum(b["n"] for b in bins)
    assert total == EXPECTED_N, f"k=24 bin sum {total} != {EXPECTED_N}"


def test_bin_partition_k32(results: dict) -> None:
    """Assert sum of per-bin counts at k=32 equals 9210."""
    bins = results["phase_variation"]["k32"]["bins"]
    total = sum(b["n"] for b in bins)
    assert total == EXPECTED_N, f"k=32 bin sum {total} != {EXPECTED_N}"


def test_b_mle_formula() -> None:
    """Assert b-value MLE formula gives expected result for synthetic data."""
    mags = np.array([6.0, 6.1, 6.2, 6.3, 6.4])
    mean_mag = np.mean(mags)
    b_computed = math.log10(math.e) / (mean_mag - MC)
    b_expected = math.log10(math.e) / 0.2   # mean - Mc = 0.2
    assert abs(b_computed - b_expected) < 0.01, (
        f"b_mle formula: got {b_computed:.4f}, expected ~{b_expected:.4f}"
    )
    # Also check the expected value itself is approximately 2.172
    assert abs(b_expected - 2.172) < 0.01, (
        f"Expected b_expected ≈ 2.172, got {b_expected:.4f}"
    )


def test_b_mle_range(results: dict) -> None:
    """Assert all computed b_mle values are in [0.5, 3.0]."""
    for k_key in ["k24", "k32"]:
        bins = results["phase_variation"][k_key]["bins"]
        for b_bin in bins:
            b_mle = b_bin["b_mle"]
            if b_mle is not None:
                assert 0.5 <= b_mle <= 3.0, (
                    f"{k_key} bin {b_bin['bin_idx']}: b_mle={b_mle:.4f} out of [0.5, 3.0]"
                )


def test_bootstrap_ci_contains_mle(results: dict) -> None:
    """Assert CI lower <= b_mle <= CI upper for each valid bin."""
    for k_key in ["k24", "k32"]:
        bins = results["phase_variation"][k_key]["bins"]
        for b_bin in bins:
            b_mle = b_bin["b_mle"]
            ci_lower = b_bin["b_ci95_lower"]
            ci_upper = b_bin["b_ci95_upper"]
            if b_mle is None or ci_lower is None or ci_upper is None:
                continue
            assert ci_lower <= b_mle, (
                f"{k_key} bin {b_bin['bin_idx']}: ci_lower={ci_lower:.4f} > b_mle={b_mle:.4f}"
            )
            assert b_mle <= ci_upper, (
                f"{k_key} bin {b_bin['bin_idx']}: b_mle={b_mle:.4f} > ci_upper={ci_upper:.4f}"
            )


def test_low_n_flag_trigger(results: dict) -> None:
    """Assert low_n_flag is True iff n < 20."""
    for k_key in ["k24", "k32"]:
        bins = results["phase_variation"][k_key]["bins"]
        for b_bin in bins:
            n = b_bin["n"]
            flag = b_bin["low_n_flag"]
            if n < 20:
                assert flag is True, (
                    f"{k_key} bin {b_bin['bin_idx']}: n={n} < 20 but low_n_flag=False"
                )
            else:
                assert flag is False, (
                    f"{k_key} bin {b_bin['bin_idx']}: n={n} >= 20 but low_n_flag=True"
                )


def test_anova_p_nonnegative(results: dict) -> None:
    """Assert p_anova is in [0.0, 1.0] for both k=24 and k=32."""
    for k_key in ["k24", "k32"]:
        p = results["phase_variation"][k_key]["p_anova"]
        assert 0.0 <= p <= 1.0, f"{k_key}: p_anova={p} out of [0.0, 1.0]"


def test_pearson_r_range(results: dict) -> None:
    """Assert r_rate_b is in [-1.0, 1.0]."""
    for k_key in ["k24", "k32"]:
        r = results["phase_variation"][k_key]["r_rate_b"]
        assert -1.0 <= r <= 1.0, f"{k_key}: r_rate_b={r} out of [-1.0, 1.0]"


def test_phase_class_valid(results: dict) -> None:
    """Assert b_max and b_min phase classes are valid strings."""
    for k_key in ["k24", "k32"]:
        pv = results["phase_variation"][k_key]
        max_class = pv["b_max_phase_class"]
        min_class = pv["b_min_phase_class"]
        assert max_class in VALID_PHASE_CLASSES, (
            f"{k_key}: b_max_phase_class '{max_class}' not in {VALID_PHASE_CLASSES}"
        )
        assert min_class in VALID_PHASE_CLASSES, (
            f"{k_key}: b_min_phase_class '{min_class}' not in {VALID_PHASE_CLASSES}"
        )


def test_results_json_keys(results: dict) -> None:
    """Assert results JSON has required top-level and sub-keys."""
    assert "phase_variation" in results, "Missing top-level key 'phase_variation'"
    pv = results["phase_variation"]
    assert "k24" in pv, "Missing sub-key 'k24' under 'phase_variation'"
    assert "k32" in pv, "Missing sub-key 'k32' under 'phase_variation'"
    # Verify per-k required fields
    required_k_fields = [
        "bins", "b_mean", "b_std", "b_range",
        "b_max_phase_center", "b_max_phase_class",
        "b_min_phase_center", "b_min_phase_class",
        "f_stat", "p_anova", "r_rate_b", "p_rate_b",
        "inverse_phase_relationship",
    ]
    for k_key in ["k24", "k32"]:
        for field in required_k_fields:
            assert field in pv[k_key], (
                f"Missing field '{field}' in phase_variation['{k_key}']"
            )
