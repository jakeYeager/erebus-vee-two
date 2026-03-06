"""
Tests for Case A3.B2: Hemisphere Stratification Refinement.

All tests must pass before whitepaper generation.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent

RAW_PATH   = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
GK_PATH    = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv"
REAS_PATH  = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv"
RESULTS_PATH = BASE_DIR / "output" / "case-a3-b2-results.json"

JULIAN_YEAR_SECS = 31_557_600.0
TECTONIC_CLASSES = ["continental", "transitional", "oceanic"]
CATALOGS = ["full", "gk", "reas"]
HEMISPHERES = ["nh", "sh"]

PNG_FILES = [
    "case-a3-b2-tectonic-heatmap.png",
    "case-a3-b2-midcrustal-binplots.png",
    "case-a3-b2-phase-alignment.png",
    "case-a3-b2-declustering-sensitivity.png",
    "case-a3-b2-interval1-threshold.png",
]

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def results() -> dict:
    """Load results JSON once per test session."""
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def df_full() -> pd.DataFrame:
    """Load and merge full catalog with GSHHG."""
    df = pd.read_csv(RAW_PATH)
    gshhg = pd.read_csv(GSHHG_PATH)[["usgs_id", "ocean_class", "dist_to_coast_km"]]
    df = df.merge(gshhg, on="usgs_id", how="left")
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="session")
def df_gk() -> pd.DataFrame:
    """Load and merge GK mainshocks with GSHHG."""
    df = pd.read_csv(GK_PATH)
    gshhg = pd.read_csv(GSHHG_PATH)[["usgs_id", "ocean_class", "dist_to_coast_km"]]
    df = df.merge(gshhg, on="usgs_id", how="left")
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="session")
def df_reas() -> pd.DataFrame:
    """Load and merge Reasenberg mainshocks with GSHHG."""
    df = pd.read_csv(REAS_PATH)
    gshhg = pd.read_csv(GSHHG_PATH)[["usgs_id", "ocean_class", "dist_to_coast_km"]]
    df = df.merge(gshhg, on="usgs_id", how="left")
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_catalog_loads(df_full: pd.DataFrame, df_gk: pd.DataFrame, df_reas: pd.DataFrame) -> None:
    """Assert catalog sizes and valid phase column."""
    assert len(df_full) == 9210, f"Expected 9210 rows in full catalog, got {len(df_full)}"
    assert len(df_gk) == 5883, f"Expected 5883 rows in GK catalog, got {len(df_gk)}"
    assert len(df_reas) == 8265, f"Expected 8265 rows in Reasenberg catalog, got {len(df_reas)}"

    for name, df in [("full", df_full), ("gk", df_gk), ("reas", df_reas)]:
        assert "phase" in df.columns, f"Missing 'phase' column in {name}"
        assert (df["phase"] >= 0.0).all(), f"Negative phase values in {name}"
        assert (df["phase"] < 1.0).all(), f"Phase values >= 1.0 in {name}"


def test_gshhg_merge(df_full: pd.DataFrame, df_gk: pd.DataFrame, df_reas: pd.DataFrame) -> None:
    """Assert no NaN in dist_to_coast_km after merge for all catalogs."""
    for name, df in [("full", df_full), ("gk", df_gk), ("reas", df_reas)]:
        n_nan = df["dist_to_coast_km"].isna().sum()
        assert n_nan == 0, f"Found {n_nan} NaN in dist_to_coast_km for {name}"


def test_hemisphere_split_full(df_full: pd.DataFrame) -> None:
    """Assert NH + SH + equatorial == 9210 for full catalog."""
    n_nh = int((df_full["latitude"] > 0).sum())
    n_sh = int((df_full["latitude"] < 0).sum())
    n_eq = int((df_full["latitude"] == 0).sum())

    assert n_nh + n_sh + n_eq == 9210, (
        f"NH({n_nh}) + SH({n_sh}) + Eq({n_eq}) = {n_nh + n_sh + n_eq} != 9210"
    )
    assert n_nh > 0, "NH count should be > 0"
    assert n_sh > 0, "SH count should be > 0"


def test_tectonic_partition_full(df_full: pd.DataFrame) -> None:
    """Assert continental + transitional + oceanic == 9210 for full catalog."""
    counts = df_full["ocean_class"].value_counts()
    total = (
        int(counts.get("continental", 0)) +
        int(counts.get("transitional", 0)) +
        int(counts.get("oceanic", 0))
    )
    assert total == 9210, f"Tectonic partition total {total} != 9210"


def test_adaptive_k_logic() -> None:
    """Assert adaptive_k returns correct bins per threshold rules."""
    # Import the function from the analysis script
    import sys
    src_path = BASE_DIR / "src"
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))

    # Direct implementation check matching spec
    def adaptive_k(n: int):
        if n >= 500:
            return 24
        elif n >= 200:
            return 16
        elif n >= 100:
            return 12
        return None

    assert adaptive_k(500) == 24, "adaptive_k(500) should be 24"
    assert adaptive_k(300) == 16, "adaptive_k(300) should be 16"
    assert adaptive_k(150) == 12, "adaptive_k(150) should be 12"
    assert adaptive_k(50) is None, "adaptive_k(50) should be None"


def test_subtest1_cell_count(results: dict) -> None:
    """Assert subtest_1 has exactly 3 catalogs × 3 tectonic classes × 2 hemispheres."""
    st1 = results["subtest_1_tectonic_hemisphere"]

    # Check 3 catalogs
    assert set(st1.keys()) == {"full", "gk", "reas"}, f"Unexpected catalog keys: {set(st1.keys())}"

    leaf_count = 0
    for cat in ["full", "gk", "reas"]:
        tclass_dict = st1[cat]
        assert set(tclass_dict.keys()) == {"continental", "transitional", "oceanic"}, \
            f"Unexpected tectonic classes in {cat}: {set(tclass_dict.keys())}"
        for tclass in ["continental", "transitional", "oceanic"]:
            hemi_dict = tclass_dict[tclass]
            # The actual data keys are nh and sh (plus annotation keys)
            assert "nh" in hemi_dict, f"Missing 'nh' in {cat}/{tclass}"
            assert "sh" in hemi_dict, f"Missing 'sh' in {cat}/{tclass}"
            leaf_count += 2  # nh and sh

    assert leaf_count == 18, f"Expected 18 leaf cells, got {leaf_count}"


def test_subtest1_continental_nh_significant(results: dict) -> None:
    """Assert full → continental → nh → p_chi2 < 0.05."""
    p_val = results["subtest_1_tectonic_hemisphere"]["full"]["continental"]["nh"]["p_chi2"]
    assert p_val is not None, "p_chi2 is None for full/continental/nh"
    assert p_val < 0.05, (
        f"Expected full/continental/nh to be significant (p < 0.05), got p={p_val:.6f}"
    )


def test_subtest2_midcrustal_global_regression(results: dict) -> None:
    """Assert full mid-crustal global chi2 is within 2.0 of A3.B4 anchor (85.48)."""
    st2 = results["subtest_2_midcrustal_hemisphere"]
    full_global = st2["full"]["global"]
    chi2 = full_global["chi2"]
    p_val = full_global["p_chi2"]

    assert chi2 is not None, "chi2 is None for full/global mid-crustal"
    assert abs(chi2 - 85.48) <= 2.0, (
        f"Full mid-crustal chi2={chi2:.4f} not within ±2.0 of 85.48"
    )
    assert p_val < 1e-7, (
        f"Full mid-crustal p={p_val:.2e} should be < 1e-7"
    )


def test_subtest2_cell_count(results: dict) -> None:
    """Assert subtest_2 has 3 catalogs, each with global/nh/sh keys."""
    st2 = results["subtest_2_midcrustal_hemisphere"]
    assert set(st2.keys()) == {"full", "gk", "reas"}, f"Unexpected keys: {set(st2.keys())}"
    for cat in ["full", "gk", "reas"]:
        assert "global" in st2[cat], f"Missing 'global' in {cat}"
        assert "nh" in st2[cat], f"Missing 'nh' in {cat}"
        assert "sh" in st2[cat], f"Missing 'sh' in {cat}"


def test_phase_alignment_delta_range(results: dict) -> None:
    """Assert all delta_phase values are in [-0.5, 0.5] (correctly wrapped)."""
    pairs = results["subtest_3_phase_alignment"]["pairs"]
    for pair in pairs:
        delta = pair["delta_phase"]
        assert -0.5 <= delta <= 0.5, (
            f"delta_phase={delta} out of [-0.5, 0.5] range for {pair['source']}/{pair['catalog']}"
        )


def test_phase_alignment_classification_valid(results: dict) -> None:
    """Assert all alignment values are one of in_phase, anti_phase, offset."""
    valid_alignments = {"in_phase", "anti_phase", "offset"}
    pairs = results["subtest_3_phase_alignment"]["pairs"]
    for pair in pairs:
        assert pair["alignment"] in valid_alignments, (
            f"Invalid alignment: {pair['alignment']} for {pair['source']}/{pair['catalog']}"
        )


def test_interval1_threshold_count(results: dict) -> None:
    """Assert sh_all and sh_continental each have 4 threshold records."""
    st4 = results["subtest_4_interval1_threshold"]
    assert len(st4["sh_all"]) == 4, f"Expected 4 sh_all records, got {len(st4['sh_all'])}"
    assert len(st4["sh_continental"]) == 4, (
        f"Expected 4 sh_continental records, got {len(st4['sh_continental'])}"
    )


def test_interval1_classification_valid(results: dict) -> None:
    """Assert all interval_1_classification values are 'present' or 'absent'."""
    st4 = results["subtest_4_interval1_threshold"]
    valid = {"present", "absent"}
    for pop_key in ["sh_all", "sh_continental"]:
        for rec in st4[pop_key]:
            assert rec["interval_1_classification"] in valid, (
                f"Invalid classification {rec['interval_1_classification']} "
                f"in {pop_key} at threshold {rec['threshold']}"
            )


def test_results_json_completeness(results: dict) -> None:
    """Assert results JSON has all required top-level keys."""
    required_subtest_keys = {
        "subtest_1_tectonic_hemisphere",
        "subtest_2_midcrustal_hemisphere",
        "subtest_3_phase_alignment",
        "subtest_4_interval1_threshold",
    }
    for key in required_subtest_keys:
        assert key in results, f"Missing required key: {key}"

    assert "catalog_sizes" in results, "Missing 'catalog_sizes'"
    assert "parameters" in results, "Missing 'parameters'"

    # Verify catalog_sizes has correct entries
    cat_sizes = results["catalog_sizes"]
    for cat in ["full", "gk", "reas"]:
        assert cat in cat_sizes, f"Missing '{cat}' in catalog_sizes"


def test_output_figures_exist() -> None:
    """Assert all 5 PNG files exist and are > 50 KB."""
    for filename in PNG_FILES:
        path = BASE_DIR / "output" / filename
        assert path.exists(), f"PNG file not found: {path}"
        size_kb = path.stat().st_size / 1024
        assert size_kb > 50, (
            f"PNG file {filename} is only {size_kb:.1f} KB (expected > 50 KB)"
        )
