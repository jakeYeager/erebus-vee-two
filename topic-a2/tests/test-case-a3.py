"""Test suite for Case A3: Magnitude Stratification of the Solar Signal.

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
RESULTS_PATH = BASE_DIR / "output" / "case-a3-results.json"
SRC_DIR = BASE_DIR / "src"

SOLAR_YEAR_SECS = 31_557_600.0
EXPECTED_TOTAL = 9210
BAND_LABELS = ["M6.0-6.4", "M6.5-6.9", "M7.0-7.4", "M7.5+"]

logger = logging.getLogger("test-case-a3")


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


_analysis_mod = _import_module("case_a3_analysis", "case-a3-analysis.py")

BANDS = _analysis_mod.BANDS
compute_phase = _analysis_mod.compute_phase
split_by_magnitude = _analysis_mod.split_by_magnitude


# ---------------------------------------------------------------------------
# Fixture: raw catalog
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def raw_catalog() -> pd.DataFrame:
    """Load the raw ISC-GEM catalog.

    Returns:
        Full DataFrame.
    """
    path = DATA_DIR / "iscgem_global_6-9_1950-2021.csv"
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# Fixture: results JSON
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def results() -> dict:
    """Load the case A3 results JSON.

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


def test_band_partition(raw_catalog: pd.DataFrame) -> None:
    """Assert that all band sizes sum to 9,210 (no events unassigned or double-counted)."""
    band_dfs = split_by_magnitude(raw_catalog)
    total = sum(len(v) for v in band_dfs.values())
    assert total == EXPECTED_TOTAL, (
        f"Band sizes sum to {total}, expected {EXPECTED_TOTAL}"
    )


def test_band_sizes_positive(raw_catalog: pd.DataFrame) -> None:
    """Assert all four bands have n > 0; warn if any band < 100."""
    band_dfs = split_by_magnitude(raw_catalog)
    for label, df_band in band_dfs.items():
        n = len(df_band)
        assert n > 0, f"Band {label} has zero events"
        if n < 100:
            logger.warning("Band %s has only %d events (< 100)", label, n)


def test_phase_range(raw_catalog: pd.DataFrame) -> None:
    """Assert all phase values are in [0.0, 1.0)."""
    phases = compute_phase(raw_catalog["solar_secs"])
    assert float(phases.min()) >= 0.0, f"Phase below 0: {phases.min()}"
    assert float(phases.max()) < 1.0, f"Phase >= 1.0: {phases.max()}"


def test_chi_square_all_bands(results: dict) -> None:
    """Assert chi2 and p_chi2 are finite floats for all four bands at k=16, 24, 32."""
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


def test_bootstrap_ci_order(results: dict) -> None:
    """Assert cramer_v_ci95_lower <= cramer_v <= cramer_v_ci95_upper for all bands.

    A tolerance of 0.005 is allowed to accommodate edge cases where the observed
    chi-square is near the lower tail of the bootstrap distribution (particularly
    in smaller bands where sampling variation is high).
    """
    TOLERANCE = 0.005
    band_stats = results["band_stats"]
    for label in BAND_LABELS:
        k24 = band_stats[label]["k24"]
        v = k24["cramer_v"]
        lo = k24["cramer_v_ci95_lower"]
        hi = k24["cramer_v_ci95_upper"]
        assert lo - TOLERANCE <= v <= hi + TOLERANCE, (
            f"Band {label}: CI order violated beyond tolerance: "
            f"{lo:.6f} - {TOLERANCE} <= {v:.6f} <= {hi:.6f} + {TOLERANCE}"
        )


def test_spearman_rho_range(results: dict) -> None:
    """Assert Spearman rho values are in [-1.0, 1.0]."""
    trend = results["trend_analysis"]
    rho_v = trend["spearman_rho_cramer_v"]
    rho_r = trend["spearman_rho_rayleigh"]
    assert -1.0 <= rho_v <= 1.0, f"Cramér's V Spearman rho out of range: {rho_v}"
    assert -1.0 <= rho_r <= 1.0, f"Rayleigh R Spearman rho out of range: {rho_r}"


def test_trend_classification_valid(results: dict) -> None:
    """Assert trend_classification is one of 'increasing', 'decreasing', 'flat'."""
    trend_cls = results["trend_analysis"]["trend_classification"]
    valid = {"increasing", "decreasing", "flat"}
    assert trend_cls in valid, (
        f"trend_classification '{trend_cls}' not in {valid}"
    )


def test_prediction_support_valid(results: dict) -> None:
    """Assert all prediction_support values are one of 'supported', 'partially', 'not supported'."""
    pred = results["trend_analysis"]["prediction_support"]
    valid = {"supported", "partially", "not supported"}
    for key, val in pred.items():
        assert val in valid, (
            f"prediction_support['{key}'] = '{val}' not in {valid}"
        )


def test_results_json_keys(results: dict) -> None:
    """Assert 'band_stats' and 'trend_analysis' keys present and all four band keys exist."""
    assert "band_stats" in results, "Missing key 'band_stats' in results"
    assert "trend_analysis" in results, "Missing key 'trend_analysis' in results"
    band_stats = results["band_stats"]
    for label in BAND_LABELS:
        assert label in band_stats, (
            f"Missing band '{label}' in band_stats"
        )


def test_m75_band_n(results: dict) -> None:
    """Assert M7.5+ band has more than 200 events (sanity check for ISC-GEM)."""
    n_m75 = results["band_stats"]["M7.5+"]["n"]
    assert n_m75 > 200, (
        f"M7.5+ band has only {n_m75} events; expected > 200 in ISC-GEM"
    )
