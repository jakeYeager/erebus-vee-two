"""Test suite for Case B4: Depth Stratification — Surface Loading Penetration Test.

All tests load pre-computed results JSON where appropriate and validate
computed values, band statistics, and result correctness.
"""

import importlib.util
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
DATA_DIR = BASE_DIR.parent / "data" / "iscgem"
RESULTS_PATH = BASE_DIR / "output" / "case-b4-results.json"
SRC_DIR = BASE_DIR / "src"

SOLAR_YEAR_SECS = 31_557_600.0
EXPECTED_TOTAL = 9210
DEEP_BAND_MIN_N = 100
BAND_LABELS = [
    "shallow_0-20km",
    "midcrustal_20-70km",
    "intermediate_70-300km",
    "deep_300km+",
]

logger = logging.getLogger("test-case-b4")


# ---------------------------------------------------------------------------
# Helper: import hyphenated module
# ---------------------------------------------------------------------------
def _import_module(name: str, filename: str):
    """Import a module from a hyphenated filename.

    Args:
        name: Module alias.
        filename: .py filename (hyphenated allowed).

    Returns:
        Loaded module.
    """
    spec = importlib.util.spec_from_file_location(
        name, SRC_DIR / filename
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_analysis_mod = _import_module("case_b4_analysis", "case-b4-analysis.py")

DEPTH_BANDS = _analysis_mod.DEPTH_BANDS
compute_phase = _analysis_mod.compute_phase
split_by_depth = _analysis_mod.split_by_depth
ZHAN_SHEARER_PHASE_MIN = _analysis_mod.ZHAN_SHEARER_PHASE_MIN
ZHAN_SHEARER_PHASE_MAX = _analysis_mod.ZHAN_SHEARER_PHASE_MAX


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def raw_catalog() -> pd.DataFrame:
    """Load the raw ISC-GEM catalog.

    Returns:
        Full DataFrame.
    """
    path = DATA_DIR / "iscgem_global_6-9_1950-2021.csv"
    return pd.read_csv(path)


@pytest.fixture(scope="session")
def depth_bands(raw_catalog: pd.DataFrame) -> tuple[dict, int]:
    """Return split depth bands and n_depth_null.

    Returns:
        Tuple of (band_dfs dict, n_depth_null).
    """
    return split_by_depth(raw_catalog)


@pytest.fixture(scope="session")
def results() -> dict:
    """Load the case B4 results JSON.

    Returns:
        Parsed results dict.
    """
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def test_catalog_load(raw_catalog: pd.DataFrame) -> None:
    """Assert the raw catalog has exactly 9,210 rows."""
    assert len(raw_catalog) == EXPECTED_TOTAL, (
        f"Expected {EXPECTED_TOTAL} rows, got {len(raw_catalog)}"
    )


def test_depth_band_partition(raw_catalog: pd.DataFrame) -> None:
    """Assert sum of all band sizes + n_depth_null equals 9,210."""
    band_dfs, n_depth_null = split_by_depth(raw_catalog)
    total_bands = sum(len(v) for v in band_dfs.values())
    grand_total = total_bands + n_depth_null
    assert grand_total == EXPECTED_TOTAL, (
        f"Band sizes ({total_bands}) + null ({n_depth_null}) = {grand_total}, "
        f"expected {EXPECTED_TOTAL}"
    )


def test_depth_band_sizes_positive(raw_catalog: pd.DataFrame) -> None:
    """Assert all four depth bands have n > 0; warn if deep band < 100."""
    band_dfs, _ = split_by_depth(raw_catalog)
    for label, df_band in band_dfs.items():
        n = len(df_band)
        assert n > 0, f"Band {label} has zero events"
        if label == "deep_300km+" and n < DEEP_BAND_MIN_N:
            logger.warning("Deep band has only %d events (< %d)", n, DEEP_BAND_MIN_N)


def test_no_overlap_between_bands(raw_catalog: pd.DataFrame) -> None:
    """Assert no usgs_id appears in more than one depth band."""
    band_dfs, _ = split_by_depth(raw_catalog)
    all_ids: list[str] = []
    for df_band in band_dfs.values():
        all_ids.extend(df_band["usgs_id"].tolist())
    unique_ids = set(all_ids)
    assert len(all_ids) == len(unique_ids), (
        f"Duplicate usgs_id found across depth bands: "
        f"{len(all_ids)} total vs {len(unique_ids)} unique"
    )


def test_phase_range(raw_catalog: pd.DataFrame) -> None:
    """Assert all phases are in [0.0, 1.0)."""
    phases = compute_phase(raw_catalog["solar_secs"])
    assert float(phases.min()) >= 0.0, f"Phase below 0: {phases.min()}"
    assert float(phases.max()) < 1.0, f"Phase >= 1.0: {phases.max()}"


def test_chi_square_all_bands(results: dict) -> None:
    """Assert chi2 and p_chi2 are finite for all bands at k=16, 24, 32."""
    band_stats = results["band_stats"]
    for label in BAND_LABELS:
        for k_key in ["k16", "k24", "k32"]:
            entry = band_stats[label][k_key]
            assert np.isfinite(entry["chi2"]), (
                f"Band {label} {k_key} chi2 is not finite: {entry['chi2']}"
            )
            assert np.isfinite(entry["p_chi2"]), (
                f"Band {label} {k_key} p_chi2 is not finite: {entry['p_chi2']}"
            )


def test_cramer_v_range(results: dict) -> None:
    """Assert Cramér's V is in [0.0, 1.0] for all bands at all k values."""
    band_stats = results["band_stats"]
    for label in BAND_LABELS:
        for k_key in ["k16", "k24", "k32"]:
            v = band_stats[label][k_key]["cramer_v"]
            assert 0.0 <= v <= 1.0, (
                f"Band {label} {k_key} Cramér's V = {v} out of [0, 1]"
            )


def test_spearman_rho_range(results: dict) -> None:
    """Assert Spearman rho is in [-1.0, 1.0]."""
    rho = results["trend_analysis"]["spearman_rho"]
    assert -1.0 <= rho <= 1.0, f"Spearman rho out of range: {rho}"


def test_monotonicity_classification_valid(results: dict) -> None:
    """Assert monotonicity_classification is one of the three valid strings."""
    mono = results["trend_analysis"]["monotonicity_classification"]
    valid = {"decreasing with depth", "increasing with depth", "non-monotonic"}
    assert mono in valid, (
        f"monotonicity_classification '{mono}' not in {valid}"
    )


def test_zhan_shearer_phase_range(results: dict) -> None:
    """Assert deep_mean_phase_in_apr_sep is bool; check April-September range uses [0.23, 0.67]."""
    trend = results["trend_analysis"]
    flag = trend["deep_mean_phase_in_apr_sep"]
    assert isinstance(flag, bool), (
        f"deep_mean_phase_in_apr_sep should be bool, got {type(flag)}"
    )

    # Verify the April–September phase range constants
    assert ZHAN_SHEARER_PHASE_MIN == 0.23, (
        f"ZHAN_SHEARER_PHASE_MIN should be 0.23, got {ZHAN_SHEARER_PHASE_MIN}"
    )
    assert ZHAN_SHEARER_PHASE_MAX == 0.67, (
        f"ZHAN_SHEARER_PHASE_MAX should be 0.67, got {ZHAN_SHEARER_PHASE_MAX}"
    )

    # Verify the flag is consistent with the stored mean phase
    mean_phase = trend["deep_mean_phase"]
    expected_flag = bool(ZHAN_SHEARER_PHASE_MIN <= mean_phase <= ZHAN_SHEARER_PHASE_MAX)
    assert flag == expected_flag, (
        f"deep_mean_phase_in_apr_sep={flag} but mean_phase={mean_phase:.4f} "
        f"expected {expected_flag} for range [{ZHAN_SHEARER_PHASE_MIN}, {ZHAN_SHEARER_PHASE_MAX}]"
    )


def test_prediction_support_valid(results: dict) -> None:
    """Assert all prediction_support values are one of the three valid strings."""
    pred = results["trend_analysis"]["prediction_support"]
    valid = {"supported", "partially", "not supported"}
    for key, val in pred.items():
        assert val in valid, (
            f"prediction_support['{key}'] = '{val}' not in {valid}"
        )


def test_results_json_keys(results: dict) -> None:
    """Assert 'band_stats' and 'trend_analysis' keys present in results JSON."""
    assert "band_stats" in results, "Missing key 'band_stats' in results"
    assert "trend_analysis" in results, "Missing key 'trend_analysis' in results"
    band_stats = results["band_stats"]
    for label in BAND_LABELS:
        assert label in band_stats, (
            f"Missing band '{label}' in band_stats"
        )
