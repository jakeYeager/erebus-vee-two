"""Test suite for Case B2: Ocean vs. Continent Location — Hydrological Loading Discrimination.

All tests load from the raw data files and/or the pre-computed results JSON at
topic-a2/output/case-b2-results.json. The analysis script must be run first.
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
PLATE_LOC_DIR = DATA_DIR / "plate-location"
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-b2-results.json"

SOLAR_YEAR_SECS: float = 31_557_600.0
EXPECTED_N: int = 9210
LOCATION_CLASSES = {"oceanic", "transitional", "continental"}
METHODS = ["gshhg", "ne", "pb2002"]
BIN_COUNTS = [16, 24, 32]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def raw_catalog() -> pd.DataFrame:
    """Load the raw ISC-GEM catalog."""
    return pd.read_csv(DATA_DIR / "iscgem_global_6-9_1950-2021.csv")


@pytest.fixture(scope="session")
def class_gshhg() -> pd.DataFrame:
    """Load GSHHG classification file."""
    return pd.read_csv(PLATE_LOC_DIR / "ocean_class_gshhg_global.csv")


@pytest.fixture(scope="session")
def class_ne() -> pd.DataFrame:
    """Load Natural Earth classification file."""
    return pd.read_csv(PLATE_LOC_DIR / "ocean_class_ne_global.csv")


@pytest.fixture(scope="session")
def class_pb2002() -> pd.DataFrame:
    """Load PB2002 classification file."""
    return pd.read_csv(PLATE_LOC_DIR / "ocean_class_pb2002_global.csv")


@pytest.fixture(scope="session")
def merged_df(raw_catalog, class_gshhg, class_ne, class_pb2002) -> pd.DataFrame:
    """Merge all classification files onto raw catalog and compute solar phase."""
    df = raw_catalog.copy()
    for key, cls_df in [
        ("gshhg", class_gshhg),
        ("ne",    class_ne),
        ("pb2002", class_pb2002),
    ]:
        renamed = cls_df.rename(columns={
            "ocean_class":      f"ocean_class_{key}",
            "dist_to_coast_km": f"dist_km_{key}",
        })
        df = df.merge(
            renamed[["usgs_id", f"ocean_class_{key}", f"dist_km_{key}"]],
            on="usgs_id", how="left",
        )
    df["solar_phase"] = (df["solar_secs"] / SOLAR_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="session")
def results() -> dict:
    """Load the case-b2-results.json output."""
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestCatalogLoad:
    """Test 1: Raw catalog loads with expected row count."""

    def test_catalog_load(self, raw_catalog: pd.DataFrame) -> None:
        """Assert raw catalog has exactly 9210 rows."""
        assert len(raw_catalog) == EXPECTED_N, (
            f"Expected {EXPECTED_N} rows, got {len(raw_catalog)}"
        )


class TestClassificationJoins:
    """Test 2: All three classification joins produce 9210 rows with no NaN."""

    def test_classification_joins(self, merged_df: pd.DataFrame) -> None:
        """Assert merged dataframe has 9210 rows and no NaN in class columns."""
        assert len(merged_df) == EXPECTED_N, (
            f"Merged dataframe: expected {EXPECTED_N} rows, got {len(merged_df)}"
        )
        for method in METHODS:
            col = f"ocean_class_{method}"
            nan_count = merged_df[col].isna().sum()
            assert nan_count == 0, (
                f"Column {col} has {nan_count} NaN values after join"
            )


class TestClassCountsGshhg:
    """Test 3: GSHHG class counts match spec (±1)."""

    def test_class_counts_gshhg(self, merged_df: pd.DataFrame) -> None:
        """Assert GSHHG oceanic=1952, transitional=3459, continental=3799 (±1)."""
        expected = {
            "oceanic":      1952,
            "transitional": 3459,
            "continental":  3799,
        }
        counts = merged_df["ocean_class_gshhg"].value_counts().to_dict()
        for class_label, exp_n in expected.items():
            actual_n = counts.get(class_label, 0)
            assert abs(actual_n - exp_n) <= 1, (
                f"GSHHG {class_label}: expected {exp_n} ± 1, got {actual_n}"
            )


class TestClassCountsPb2002:
    """Test 4: PB2002 class counts match spec (±1)."""

    def test_class_counts_pb2002(self, merged_df: pd.DataFrame) -> None:
        """Assert PB2002 oceanic=2682, transitional=3677, continental=2851 (±1)."""
        expected = {
            "oceanic":      2682,
            "transitional": 3677,
            "continental":  2851,
        }
        counts = merged_df["ocean_class_pb2002"].value_counts().to_dict()
        for class_label, exp_n in expected.items():
            actual_n = counts.get(class_label, 0)
            assert abs(actual_n - exp_n) <= 1, (
                f"PB2002 {class_label}: expected {exp_n} ± 1, got {actual_n}"
            )


class TestPhaseRange:
    """Test 5: All solar_phase values are in [0.0, 1.0)."""

    def test_phase_range(self, merged_df: pd.DataFrame) -> None:
        """Assert all phases are in [0, 1)."""
        phase = merged_df["solar_phase"]
        assert (phase >= 0.0).all(), "Some phases are negative"
        assert (phase < 1.0).all(), "Some phases are >= 1.0"


class TestClassPartition:
    """Test 6: For each method, the three classes sum to 9210."""

    def test_class_partition(self, merged_df: pd.DataFrame) -> None:
        """Assert oceanic + transitional + continental == 9210 for each method."""
        for method in METHODS:
            col = f"ocean_class_{method}"
            total = 0
            for class_label in LOCATION_CLASSES:
                total += (merged_df[col] == class_label).sum()
            assert total == EXPECTED_N, (
                f"Method {method}: class partition sum={total}, expected {EXPECTED_N}"
            )


class TestChiSquareAllClasses:
    """Test 7: Chi-square and p-value are finite for all method × class × k combos."""

    def test_chi_square_all_classes(self, results: dict) -> None:
        """Assert chi2 and p_chi2 are finite for all 3 methods × 3 classes × 3 bin counts."""
        class_stats = results["class_stats"]
        for method in METHODS:
            for class_label in LOCATION_CLASSES:
                for k in BIN_COUNTS:
                    k_key = f"k{k}"
                    stats = class_stats[method][class_label][k_key]
                    chi2 = stats["chi2"]
                    p_chi2 = stats["p_chi2"]
                    assert math.isfinite(chi2), (
                        f"{method}/{class_label}/{k_key}: chi2={chi2} is not finite"
                    )
                    assert math.isfinite(p_chi2), (
                        f"{method}/{class_label}/{k_key}: p_chi2={p_chi2} is not finite"
                    )


class TestCramerVRange:
    """Test 8: All Cramér's V values are in [0.0, 1.0]."""

    def test_cramer_v_range(self, results: dict) -> None:
        """Assert Cramér's V is in [0, 1] for all method × class × k combos."""
        class_stats = results["class_stats"]
        for method in METHODS:
            for class_label in LOCATION_CLASSES:
                for k in BIN_COUNTS:
                    k_key = f"k{k}"
                    v = class_stats[method][class_label][k_key]["cramer_v"]
                    assert 0.0 <= v <= 1.0, (
                        f"{method}/{class_label}/{k_key}: cramer_v={v} outside [0,1]"
                    )


class TestMethodSensitivityStructure:
    """Test 9: method_sensitivity key exists with correct structure."""

    def test_method_sensitivity_structure(self, results: dict) -> None:
        """Assert method_sensitivity present with per-class entries and agreement field."""
        assert "method_sensitivity" in results, "Key 'method_sensitivity' missing from results"
        sensitivity = results["method_sensitivity"]
        assert isinstance(sensitivity, list), "method_sensitivity should be a list"
        class_labels_found = {entry["class"] for entry in sensitivity}
        for class_label in LOCATION_CLASSES:
            assert class_label in class_labels_found, (
                f"Class '{class_label}' missing from method_sensitivity"
            )
            entry = next(e for e in sensitivity if e["class"] == class_label)
            for method in METHODS:
                assert f"{method}_p" in entry, (
                    f"method_sensitivity entry for {class_label} missing '{method}_p'"
                )
                assert f"{method}_cramer_v" in entry, (
                    f"method_sensitivity entry for {class_label} missing '{method}_cramer_v'"
                )
            assert "agreement" in entry, (
                f"method_sensitivity entry for {class_label} missing 'agreement'"
            )
            assert isinstance(entry["agreement"], bool), (
                f"agreement field for {class_label} is not bool"
            )


class TestPredictionEvaluationKeys:
    """Test 10: prediction_evaluation key present with required fields."""

    def test_prediction_evaluation_keys(self, results: dict) -> None:
        """Assert prediction_evaluation present with valid primary_conclusion."""
        assert "prediction_evaluation" in results, (
            "Key 'prediction_evaluation' missing from results"
        )
        pred = results["prediction_evaluation"]
        required_keys = [
            "primary_method",
            "oceanic_significant",
            "continental_significant",
            "continental_stronger",
            "oceanic_matches_a1b_intervals",
            "geometric_supported",
            "hydrological_supported",
            "primary_conclusion",
        ]
        for key in required_keys:
            assert key in pred, (
                f"prediction_evaluation missing key: '{key}'"
            )
        valid_conclusions = {"geometric", "hydrological", "ambiguous"}
        conclusion = pred["primary_conclusion"]
        assert conclusion in valid_conclusions, (
            f"primary_conclusion='{conclusion}' is not one of {valid_conclusions}"
        )


class TestClassificationLabelsValid:
    """Test 11: All ocean_class values are valid labels."""

    def test_classification_labels_valid(self, merged_df: pd.DataFrame) -> None:
        """Assert all ocean_class_* values are oceanic, transitional, or continental."""
        for method in METHODS:
            col = f"ocean_class_{method}"
            unique_labels = set(merged_df[col].dropna().unique())
            invalid = unique_labels - LOCATION_CLASSES
            assert not invalid, (
                f"Column {col} contains invalid labels: {invalid}"
            )
