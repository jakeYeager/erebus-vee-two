"""Test suite for Case B3: Tectonic Regime Stratification.

Tests load from the raw focal join file and/or the pre-computed results JSON at
topic-a2/output/case-b3-results.json. The analysis script must be run first.
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
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
DATA_DIR = BASE_DIR.parent / "data" / "iscgem"
FOCAL_DIR = DATA_DIR / "focal-mechanism"
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-b3-results.json"

SOLAR_YEAR_SECS: float = 31_557_600.0
EXPECTED_N: int = 9210
MECHANISM_CLASSES = {"thrust", "normal", "strike_slip"}
BIN_COUNTS = [16, 24, 32]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def focal_join() -> pd.DataFrame:
    """Load the GCMT focal join file."""
    return pd.read_csv(FOCAL_DIR / "focal_join_global.csv")


@pytest.fixture(scope="session")
def results() -> dict:
    """Load the case-b3-results.json output."""
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def focal_with_phase(focal_join: pd.DataFrame) -> pd.DataFrame:
    """Add solar_phase and tectonic_class to the focal join for testing."""
    df = focal_join.copy()
    df["solar_phase"] = (df["solar_secs"] / SOLAR_YEAR_SECS) % 1.0

    def classify_rake(rake: float) -> str:
        """Classify focal mechanism from rake angle."""
        while rake > 180:
            rake -= 360
        while rake <= -180:
            rake += 360
        if 45 <= rake <= 135:
            return "thrust"
        elif -135 <= rake <= -45:
            return "normal"
        elif (-45 < rake <= 45) or (135 < rake <= 180) or (-180 <= rake < -135):
            return "strike_slip"
        else:
            return "oblique"

    df["tectonic_class"] = "unmatched"
    matched_mask = df["match_confidence"] == "proximity"
    has_mech = matched_mask & df["mechanism"].notna()

    for mech in ["thrust", "normal", "strike_slip"]:
        df.loc[has_mech & (df["mechanism"] == mech), "tectonic_class"] = mech
    df.loc[has_mech & (df["mechanism"] == "oblique"), "tectonic_class"] = "oblique"

    needs_fallback = matched_mask & df["mechanism"].isna() & df["rake"].notna()
    for idx in df.index[needs_fallback]:
        df.at[idx, "tectonic_class"] = classify_rake(df.at[idx, "rake"])

    return df


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestFocalJoinLoad:
    """test_focal_join_load: Focal join file loads with n=9210."""

    def test_focal_join_load(self, focal_join: pd.DataFrame) -> None:
        """Assert focal join has exactly 9210 rows."""
        assert len(focal_join) == EXPECTED_N, (
            f"Expected {EXPECTED_N} rows, got {len(focal_join)}"
        )


class TestMechanismCoverage:
    """test_mechanism_coverage: Class count partition sums to 9210."""

    def test_mechanism_coverage(self, focal_with_phase: pd.DataFrame) -> None:
        """Assert n_thrust + n_normal + n_strike_slip + n_oblique + n_unmatched == 9210."""
        counts = focal_with_phase["tectonic_class"].value_counts()
        n_thrust = int(counts.get("thrust", 0))
        n_normal = int(counts.get("normal", 0))
        n_strike_slip = int(counts.get("strike_slip", 0))
        n_oblique = int(counts.get("oblique", 0))
        n_unmatched = int(counts.get("unmatched", 0))
        total = n_thrust + n_normal + n_strike_slip + n_oblique + n_unmatched
        assert total == EXPECTED_N, (
            f"mechanism counts sum to {total}, expected {EXPECTED_N}; "
            f"thrust={n_thrust} normal={n_normal} strike_slip={n_strike_slip} "
            f"oblique={n_oblique} unmatched={n_unmatched}"
        )


class TestMatchedCountRange:
    """test_matched_count_range: Total matched events in expected range."""

    def test_matched_count_range(self, results: dict) -> None:
        """Assert n_total_matched (thrust + normal + strike_slip + oblique) is 4800–5000."""
        cov = results["mechanism_stats"]["coverage"]
        n_total_matched = (
            cov["n_thrust"] + cov["n_normal"] + cov["n_strike_slip"] + cov["n_oblique"]
        )
        assert 4800 <= n_total_matched <= 5000, (
            f"n_total_matched={n_total_matched} is outside expected range [4800, 5000]"
        )


class TestPhaseRange:
    """test_phase_range: All solar phases are in [0.0, 1.0)."""

    def test_phase_range(self, focal_with_phase: pd.DataFrame) -> None:
        """Assert all solar_phase values are in [0, 1)."""
        phase = focal_with_phase["solar_phase"]
        assert (phase >= 0.0).all(), "Some solar_phase values are negative"
        assert (phase < 1.0).all(), "Some solar_phase values are >= 1.0"


class TestChiSquareAllMechanisms:
    """test_chi_square_all_mechanisms: chi2 and p_chi2 are finite for all classes."""

    def test_chi_square_all_mechanisms(self, results: dict) -> None:
        """Assert chi2 and p_chi2 are finite for all three mechanisms at k=16, 24, 32."""
        mech_stats = results["mechanism_stats"]
        for mech in MECHANISM_CLASSES:
            for k in BIN_COUNTS:
                k_key = f"k{k}"
                stats = mech_stats[mech][k_key]
                chi2 = stats["chi2"]
                p_chi2 = stats["p_chi2"]
                assert math.isfinite(chi2), (
                    f"{mech}/{k_key}: chi2={chi2} is not finite"
                )
                assert math.isfinite(p_chi2), (
                    f"{mech}/{k_key}: p_chi2={p_chi2} is not finite"
                )


class TestCramerVRange:
    """test_cramer_v_range: Cramér's V in [0.0, 1.0] for all mechanisms."""

    def test_cramer_v_range(self, results: dict) -> None:
        """Assert Cramér's V is in [0, 1] for all three mechanisms at all k."""
        mech_stats = results["mechanism_stats"]
        for mech in MECHANISM_CLASSES:
            for k in BIN_COUNTS:
                k_key = f"k{k}"
                v = mech_stats[mech][k_key]["cramer_v"]
                assert 0.0 <= v <= 1.0, (
                    f"{mech}/{k_key}: cramer_v={v} outside [0, 1]"
                )


class TestMechanismClassLabelsValid:
    """test_mechanism_class_labels_valid: cramer_v_rank_order contains valid labels."""

    def test_mechanism_class_labels_valid(self, results: dict) -> None:
        """Assert cramer_v_rank_order contains exactly thrust, normal, strike_slip."""
        rank_order = results["prediction_evaluation"]["cramer_v_rank_order"]
        assert len(rank_order) == 3, (
            f"cramer_v_rank_order has {len(rank_order)} elements, expected 3"
        )
        assert set(rank_order) == MECHANISM_CLASSES, (
            f"cramer_v_rank_order contains invalid labels: {set(rank_order)}"
        )


class TestRakeClassificationLogic:
    """test_rake_classification_logic: rake quadrant boundary logic is correct."""

    def test_rake_classification_logic(self) -> None:
        """Assert specific rake values classify correctly per spec quadrant boundaries."""
        def classify_rake(rake: float) -> str:
            while rake > 180:
                rake -= 360
            while rake <= -180:
                rake += 360
            if 45 <= rake <= 135:
                return "thrust"
            elif -135 <= rake <= -45:
                return "normal"
            elif (-45 < rake <= 45) or (135 < rake <= 180) or (-180 <= rake < -135):
                return "strike_slip"
            else:
                return "oblique"

        assert classify_rake(90) == "thrust", "rake=90 should be thrust"
        assert classify_rake(-90) == "normal", "rake=-90 should be normal"
        assert classify_rake(0) == "strike_slip", "rake=0 should be strike_slip"
        assert classify_rake(170) == "strike_slip", "rake=170 should be strike_slip"


class TestPhaseOffsetRange:
    """test_phase_offset_range: phase_offset_thrust_normal in [0.0, 0.5]."""

    def test_phase_offset_range(self, results: dict) -> None:
        """Assert thrust-normal phase offset is in [0.0, 0.5] (wrapped)."""
        offset = results["prediction_evaluation"]["phase_offset_thrust_normal"]
        assert 0.0 <= offset <= 0.5, (
            f"phase_offset_thrust_normal={offset} outside [0.0, 0.5]"
        )


class TestUnmatchedCheckPresent:
    """test_unmatched_check_present: unmatched_check key present with n_unmatched."""

    def test_unmatched_check_present(self, results: dict) -> None:
        """Assert unmatched_check key exists and contains n_unmatched field."""
        assert "unmatched_check" in results, (
            "Key 'unmatched_check' missing from results JSON"
        )
        uc = results["unmatched_check"]
        assert "n_unmatched" in uc, (
            "Key 'n_unmatched' missing from unmatched_check"
        )
        assert isinstance(uc["n_unmatched"], int), (
            f"n_unmatched should be int, got {type(uc['n_unmatched'])}"
        )


class TestPrimaryConclusionValid:
    """test_primary_conclusion_valid: primary_conclusion is one of the valid values."""

    def test_primary_conclusion_valid(self, results: dict) -> None:
        """Assert primary_conclusion is loading, geometric, metivier_consistent, or ambiguous."""
        valid = {"loading", "geometric", "metivier_consistent", "ambiguous"}
        conclusion = results["prediction_evaluation"]["primary_conclusion"]
        assert conclusion in valid, (
            f"primary_conclusion='{conclusion}' is not one of {valid}"
        )
