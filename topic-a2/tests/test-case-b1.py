"""Test suite for Case B1: Hemisphere Stratification — Phase Symmetry Test.

All tests load pre-computed results JSON where appropriate and validate
computed values, symmetry test logic, and result correctness.
"""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
DATA_DIR = BASE_DIR.parent / "data" / "iscgem"
RESULTS_PATH = BASE_DIR / "output" / "case-b1-results.json"
SRC_DIR = BASE_DIR / "src"

SOLAR_YEAR_SECS = 31_557_600.0


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


_analysis_mod = _import_module("case_b1_analysis", "case-b1-analysis.py")

compute_phase = _analysis_mod.compute_phase
compute_half_cycle_offset = _analysis_mod.compute_half_cycle_offset


# ---------------------------------------------------------------------------
# Fixture: results JSON
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def results() -> dict:
    """Load the case B1 results JSON.

    Returns:
        Parsed results dictionary.
    """
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# test_hemisphere_split
# ---------------------------------------------------------------------------
def test_hemisphere_split(results: dict):
    """Assert NH + SH + equatorial = 9210; assert both hemispheres non-empty."""
    hemi = results["hemisphere_stats"]
    n_nh = hemi["n_nh"]
    n_sh = hemi["n_sh"]
    n_eq = hemi["n_equatorial"]

    assert n_nh + n_sh + n_eq == 9210, (
        f"NH ({n_nh}) + SH ({n_sh}) + equatorial ({n_eq}) != 9210"
    )
    assert n_nh > 0, "n_nh should be > 0"
    assert n_sh > 0, "n_sh should be > 0"


# ---------------------------------------------------------------------------
# test_phase_range
# ---------------------------------------------------------------------------
def test_phase_range():
    """Assert all computed phases are in [0.0, 1.0)."""
    df = pd.read_csv(DATA_DIR / "iscgem_global_6-9_1950-2021.csv")
    phases = compute_phase(df["solar_secs"])
    assert (phases >= 0.0).all(), "Some phases < 0.0"
    assert (phases < 1.0).all(), "Some phases >= 1.0"


# ---------------------------------------------------------------------------
# test_chi_square_both_hemispheres
# ---------------------------------------------------------------------------
def test_chi_square_both_hemispheres(results: dict):
    """Assert chi2 and p_chi2 are finite floats for both hemispheres at all k."""
    hemi = results["hemisphere_stats"]
    for hemi_key in ["nh", "sh"]:
        for k_key in ["k16", "k24", "k32"]:
            data = hemi[hemi_key][k_key]
            chi2 = data["chi2"]
            p_chi2 = data["p_chi2"]
            assert np.isfinite(chi2), (
                f"{hemi_key} {k_key} chi2 is not finite: {chi2}"
            )
            assert np.isfinite(p_chi2), (
                f"{hemi_key} {k_key} p_chi2 is not finite: {p_chi2}"
            )
            assert isinstance(chi2, float), (
                f"{hemi_key} {k_key} chi2 should be float, got {type(chi2)}"
            )
            assert isinstance(p_chi2, float), (
                f"{hemi_key} {k_key} p_chi2 should be float, got {type(p_chi2)}"
            )


# ---------------------------------------------------------------------------
# test_cramer_v_range
# ---------------------------------------------------------------------------
def test_cramer_v_range(results: dict):
    """Assert Cramér's V for both hemispheres at all k is in [0.0, 1.0]."""
    hemi = results["hemisphere_stats"]
    for hemi_key in ["nh", "sh"]:
        for k_key in ["k16", "k24", "k32"]:
            v = hemi[hemi_key][k_key]["cramer_v"]
            assert 0.0 <= v <= 1.0, (
                f"Cramér's V out of range for {hemi_key}/{k_key}: {v}"
            )


# ---------------------------------------------------------------------------
# test_interval_classification_keys
# ---------------------------------------------------------------------------
def test_interval_classification_keys(results: dict):
    """Assert symmetry_tests key is present with all four sub-keys."""
    assert "symmetry_tests" in results, "Missing 'symmetry_tests' key"
    st = results["symmetry_tests"]
    required_keys = [
        "test_1_global",
        "test_2_interval_1",
        "test_3_interval_23_specificity",
        "test_4_half_cycle_offset",
    ]
    for key in required_keys:
        assert key in st, f"Missing symmetry test key: '{key}'"


# ---------------------------------------------------------------------------
# test_prediction_support_valid_values
# ---------------------------------------------------------------------------
def test_prediction_support_valid_values(results: dict):
    """Assert prediction_support values are valid strings."""
    valid_support = {"supported", "partially supported", "not supported"}
    valid_conclusions = {"geometric", "hydrological", "mixed", "ambiguous"}

    assert "prediction_support" in results, "Missing 'prediction_support' key"
    ps = results["prediction_support"]

    for hypothesis in ["geometric", "hydrological", "mixed"]:
        assert ps[hypothesis] in valid_support, (
            f"prediction_support['{hypothesis}'] = '{ps[hypothesis]}' not in {valid_support}"
        )

    assert ps["primary_conclusion"] in valid_conclusions, (
        f"primary_conclusion = '{ps['primary_conclusion']}' not in {valid_conclusions}"
    )


# ---------------------------------------------------------------------------
# test_elevated_intervals_phase_range
# ---------------------------------------------------------------------------
def test_elevated_intervals_phase_range(results: dict):
    """Assert all elevated interval phase_start and phase_end values are in [0.0, 1.0]."""
    hemi = results["hemisphere_stats"]
    for hemi_key in ["nh", "sh"]:
        for k_key in ["k16", "k24", "k32"]:
            intervals = hemi[hemi_key][k_key]["elevated_intervals"]
            for interval in intervals:
                ps = interval["phase_start"]
                pe = interval["phase_end"]
                assert 0.0 <= ps <= 1.0, (
                    f"{hemi_key}/{k_key} interval phase_start={ps} out of [0,1]"
                )
                assert 0.0 <= pe <= 1.0, (
                    f"{hemi_key}/{k_key} interval phase_end={pe} out of [0,1]"
                )


# ---------------------------------------------------------------------------
# test_half_cycle_offset_logic
# ---------------------------------------------------------------------------
def test_half_cycle_offset_logic():
    """Synthetic NH interval at center 0.22: expected SH counterpart is 0.72."""
    # Create a synthetic NH elevated interval at mean_phase = 0.22
    elevated_nh = [{"phase_start": 0.20, "phase_end": 0.24, "mean_phase": 0.22}]
    # Create a synthetic SH elevated interval at the expected counterpart
    elevated_sh = [{"phase_start": 0.70, "phase_end": 0.74, "mean_phase": 0.72}]

    result = compute_half_cycle_offset(elevated_nh, elevated_sh)

    assert len(result["details"]) == 1, (
        f"Expected 1 detail entry, got {len(result['details'])}"
    )
    detail = result["details"][0]
    expected_counterpart = (0.22 + 0.5) % 1.0  # = 0.72
    assert abs(detail["expected_sh_counterpart"] - expected_counterpart) < 1e-6, (
        f"Expected counterpart {expected_counterpart}, got {detail['expected_sh_counterpart']}"
    )


# ---------------------------------------------------------------------------
# test_all_k_computed
# ---------------------------------------------------------------------------
def test_all_k_computed(results: dict):
    """Assert results JSON contains k16, k24, k32 for both hemispheres."""
    hemi = results["hemisphere_stats"]
    assert "nh" in hemi, "Missing 'nh' in hemisphere_stats"
    assert "sh" in hemi, "Missing 'sh' in hemisphere_stats"

    for hemi_key in ["nh", "sh"]:
        for k_key in ["k16", "k24", "k32"]:
            assert k_key in hemi[hemi_key], (
                f"Missing '{k_key}' in hemisphere_stats['{hemi_key}']"
            )
