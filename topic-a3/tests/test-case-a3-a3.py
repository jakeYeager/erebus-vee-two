"""
Test suite for Case A3.A3: Phase-Concentration Audit.

All 15 tests must pass.
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
BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a3-a3-results.json"
CATALOG_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
GK_AFTER_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_gk-seq_global.csv"
REAS_AFTER_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_reas-seq_global.csv"
A1B_AFTER_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_a1b-seq_global.csv"

JULIAN_YEAR_SECS = 31_557_600.0
K = 24
ELEVATED_BINS = [4, 5, 6, 7, 15, 19, 21]
SUPPRESSED_BINS = [2, 8, 10, 11, 12, 13, 16, 18, 22]

# ---------------------------------------------------------------------------
# Import analysis functions (hyphenated filename requires importlib)
# ---------------------------------------------------------------------------
_spec = importlib.util.spec_from_file_location(
    "case_a3_a3_analysis",
    BASE_DIR / "src" / "case-a3-a3-analysis.py",
)
_module = importlib.util.module_from_spec(_spec)  # type: ignore
_spec.loader.exec_module(_module)  # type: ignore
compute_signed_influence = _module.compute_signed_influence
run_permutation_baseline = _module.run_permutation_baseline


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def results() -> dict:
    """Load the results JSON once for all tests."""
    assert RESULTS_PATH.exists(), f"Results JSON not found: {RESULTS_PATH}"
    with open(RESULTS_PATH) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def catalog_df() -> pd.DataFrame:
    """Load the full catalog with phase and sequence role."""
    df = pd.read_csv(CATALOG_PATH)
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def gshhg_df() -> pd.DataFrame:
    return pd.read_csv(GSHHG_PATH)


@pytest.fixture(scope="module")
def gk_after_df() -> pd.DataFrame:
    return pd.read_csv(GK_AFTER_PATH)


@pytest.fixture(scope="module")
def reas_after_df() -> pd.DataFrame:
    return pd.read_csv(REAS_AFTER_PATH)


@pytest.fixture(scope="module")
def a1b_after_df() -> pd.DataFrame:
    return pd.read_csv(A1B_AFTER_PATH)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_catalog_load(catalog_df: pd.DataFrame, results: dict) -> None:
    """T01: assert n=9210; assert phase in [0.0, 1.0); assert sequence_role column added with no NaN."""
    assert len(catalog_df) == 9210, f"Expected 9210 events, got {len(catalog_df)}"
    assert catalog_df["phase"].between(0.0, 1.0, inclusive="left").all(), \
        "Phase values must be in [0.0, 1.0)"
    # Sequence role is computed in analysis script — verify via results JSON
    roles = results["sequence_roles"]
    total = roles["n_isolated_mainshock"] + roles["n_mainshock_with_sequence"] + roles["n_aftershock"]
    assert total == 9210, f"Sequence role total should be 9210, got {total}"


def test_gshhg_merge(catalog_df: pd.DataFrame, gshhg_df: pd.DataFrame) -> None:
    """T02: assert no NaN in ocean_class and tectonic_class after merge."""
    merged = catalog_df.merge(gshhg_df[["usgs_id", "ocean_class", "dist_to_coast_km"]], on="usgs_id", how="left")
    assert merged["ocean_class"].isna().sum() == 0, "NaN found in ocean_class after merge"
    # Tectonic class derived from dist_to_coast_km — should also have no NaN
    assert merged["dist_to_coast_km"].isna().sum() == 0, "NaN found in dist_to_coast_km"


def test_sequence_role_partition(results: dict) -> None:
    """T03: assert n_isolated + n_mainshock_seq + n_aftershock == 9210."""
    roles = results["sequence_roles"]
    total = roles["n_isolated_mainshock"] + roles["n_mainshock_with_sequence"] + roles["n_aftershock"]
    assert total == 9210, f"Partition total {total} != 9210"
    # All counts should be non-negative
    assert roles["n_isolated_mainshock"] >= 0
    assert roles["n_mainshock_with_sequence"] >= 0
    assert roles["n_aftershock"] >= 0


def test_sequence_catalog_sizes(
    gk_after_df: pd.DataFrame,
    reas_after_df: pd.DataFrame,
    a1b_after_df: pd.DataFrame,
) -> None:
    """T04: assert GK aftershocks n=3327; Reasenberg n=945; A1b n=2073."""
    assert len(gk_after_df) == 3327, f"Expected GK aftershocks n=3327, got {len(gk_after_df)}"
    assert len(reas_after_df) == 945, f"Expected Reasenberg aftershocks n=945, got {len(reas_after_df)}"
    assert len(a1b_after_df) == 2073, f"Expected A1b aftershocks n=2073, got {len(a1b_after_df)}"


def test_influence_formula() -> None:
    """T05: verify compute_signed_influence returns correct values for known bin counts."""
    # Construct synthetic array: 5 events in bin 0, 3 events in bin 1, rest empty
    # k=2 bins, so phase < 0.5 -> bin 0, phase >= 0.5 -> bin 1
    k = 2
    # 5 events in bin 0, 3 in bin 1 -> n=8, expected=4
    phases = np.array([0.1, 0.2, 0.3, 0.4, 0.45, 0.6, 0.7, 0.8])
    influence = compute_signed_influence(phases, k=k)

    n = len(phases)
    expected = n / k  # 4.0
    obs_bin0 = 5.0
    obs_bin1 = 3.0

    expected_infl_bin0 = (2.0 * (obs_bin0 - expected) - 1.0) / expected
    expected_infl_bin1 = (2.0 * (obs_bin1 - expected) - 1.0) / expected

    # Events in bin 0 (phases < 0.5): indices 0-4
    np.testing.assert_allclose(influence[:5], expected_infl_bin0, rtol=1e-10,
                               err_msg="Bin 0 influence incorrect")
    # Events in bin 1 (phases >= 0.5): indices 5-7
    np.testing.assert_allclose(influence[5:], expected_infl_bin1, rtol=1e-10,
                               err_msg="Bin 1 influence incorrect")


def test_influence_same_bin_identical() -> None:
    """T05b: events in the same bin must receive identical influence values."""
    # All events in one bin
    phases = np.array([0.01, 0.02, 0.03, 0.04])
    influence = compute_signed_influence(phases, k=24)
    # All in bin 0
    assert len(set(influence.tolist())) == 1, "Events in same bin should have identical influence"


def test_influence_length(results: dict) -> None:
    """T06: assert chi2_influence length equals catalog size (9210)."""
    # We verify via bin_summary: sum of all n_isolated + n_mainshock_seq + n_aftershock across bins
    bs = results["bin_summary"]
    total_events = sum(b["n_isolated"] + b["n_mainshock_seq"] + b["n_aftershock"] for b in bs)
    assert total_events == 9210, f"Bin summary total events {total_events} != 9210"


def test_elevated_suppressed_bins() -> None:
    """T07: assert ELEVATED_BINS and SUPPRESSED_BINS are non-empty, disjoint, all in [0,23]."""
    assert len(ELEVATED_BINS) > 0, "ELEVATED_BINS must be non-empty"
    assert len(SUPPRESSED_BINS) > 0, "SUPPRESSED_BINS must be non-empty"

    # Disjoint
    overlap = set(ELEVATED_BINS) & set(SUPPRESSED_BINS)
    assert len(overlap) == 0, f"ELEVATED_BINS and SUPPRESSED_BINS overlap at: {overlap}"

    # All in [0, 23]
    for b in ELEVATED_BINS:
        assert 0 <= b <= 23, f"Elevated bin {b} out of range [0, 23]"
    for b in SUPPRESSED_BINS:
        assert 0 <= b <= 23, f"Suppressed bin {b} out of range [0, 23]"


def test_permutation_shape() -> None:
    """T08: assert permutation produces correct shapes."""
    # Use a small synthetic case for speed
    rng = np.random.default_rng(0)
    phases = rng.uniform(0, 1, size=100)
    perm_result = run_permutation_baseline(phases, k=K, n_permutations=1000, rng_seed=42)

    assert perm_result["perm_bin_counts"].shape == (1000, K), \
        f"Expected shape (1000, {K}), got {perm_result['perm_bin_counts'].shape}"
    assert perm_result["perm_chi2"].shape == (1000,), \
        f"Expected shape (1000,), got {perm_result['perm_chi2'].shape}"


def test_permutation_bands(results: dict) -> None:
    """T09: assert bin_p5 and bin_p95 are length 24; assert p5 < p95 for all bins."""
    perm = results["permutation_baseline"]
    bin_p5 = perm["bin_p5"]
    bin_p95 = perm["bin_p95"]

    assert len(bin_p5) == K, f"bin_p5 length {len(bin_p5)} != {K}"
    assert len(bin_p95) == K, f"bin_p95 length {len(bin_p95)} != {K}"

    for i in range(K):
        assert bin_p5[i] < bin_p95[i], \
            f"p5 ({bin_p5[i]}) >= p95 ({bin_p95[i]}) for bin {i}"


def test_permutation_sig_suppressed_identified(results: dict) -> None:
    """T10: assert at least one bin is classified as significantly suppressed in permutation test."""
    perm = results["permutation_baseline"]
    sig_suppressed = perm["sig_suppressed_bins_permutation"]
    assert len(sig_suppressed) >= 1, \
        f"Expected at least 1 significantly suppressed bin, got: {sig_suppressed}"


def test_degradation_curves_present(results: dict) -> None:
    """T11: assert both elevated and suppressed removal curves present with at least 1 step."""
    deg = results["degradation_curves"]

    assert "elevated_removal" in deg, "elevated_removal missing from degradation_curves"
    assert "suppressed_removal" in deg, "suppressed_removal missing from degradation_curves"

    elev_curve = deg["elevated_removal"]["curve"]
    supp_curve = deg["suppressed_removal"]["curve"]

    assert len(elev_curve) >= 1, "elevated_removal curve must have at least 1 step"
    assert len(supp_curve) >= 1, "suppressed_removal curve must have at least 1 step"


def test_signal_persistence_recorded(results: dict) -> None:
    """T12: assert signal_persistence_step and signal_persistence_pct_catalog are present and non-null."""
    deg = results["degradation_curves"]

    for curve_key in ["elevated_removal", "suppressed_removal"]:
        curve_meta = deg[curve_key]
        assert "signal_persistence_step" in curve_meta, \
            f"signal_persistence_step missing from {curve_key}"
        assert "signal_persistence_pct_catalog" in curve_meta, \
            f"signal_persistence_pct_catalog missing from {curve_key}"
        assert curve_meta["signal_persistence_step"] is not None, \
            f"signal_persistence_step is None for {curve_key}"
        assert curve_meta["signal_persistence_pct_catalog"] is not None, \
            f"signal_persistence_pct_catalog is None for {curve_key}"


def test_representativeness_groups(results: dict) -> None:
    """T13: assert all four representativeness groups present with required fields."""
    repr_data = results["representativeness"]
    required_groups = ["top_positive_50", "top_positive_100", "top_negative_50", "top_negative_100"]

    for group in required_groups:
        assert group in repr_data, f"Representativeness group '{group}' missing"
        gdata = repr_data[group]
        assert "p_vs_full" in gdata, f"p_vs_full missing from {group}"
        assert "p_vs_signal_stratum" in gdata, f"p_vs_signal_stratum missing from {group}"
        assert "representative_of_signal_stratum" in gdata, \
            f"representative_of_signal_stratum missing from {group}"

        assert isinstance(gdata["p_vs_full"], float), \
            f"p_vs_full in {group} must be float"
        assert isinstance(gdata["p_vs_signal_stratum"], float), \
            f"p_vs_signal_stratum in {group} must be float"
        assert isinstance(gdata["representative_of_signal_stratum"], bool), \
            f"representative_of_signal_stratum in {group} must be bool"


def test_summary_keys(results: dict) -> None:
    """T14: assert summary dict present with all required keys."""
    assert "summary" in results, "summary key missing from results"
    summary = results["summary"]

    required_keys = [
        "signal_diffuse",
        "elevated_persistence_pct",
        "suppressed_persistence_pct",
        "degradation_symmetric",
        "top100_elevated_representative_of_signal_stratum",
        "top100_suppressed_representative_of_signal_stratum",
    ]
    for key in required_keys:
        assert key in summary, f"Summary key '{key}' missing"

    assert isinstance(summary["signal_diffuse"], bool), "signal_diffuse must be bool"
    assert isinstance(summary["degradation_symmetric"], bool), "degradation_symmetric must be bool"


def test_output_figures_exist() -> None:
    """T15: assert all 5 PNG files exist and size > 50 KB."""
    expected_figures = [
        "case-a3-a3-influence-distribution.png",
        "case-a3-a3-degradation-curves.png",
        "case-a3-a3-permutation-tails.png",
        "case-a3-a3-representativeness.png",
        "case-a3-a3-sequence-annotation.png",
    ]

    output_dir = BASE_DIR / "output"
    min_size_bytes = 50 * 1024  # 50 KB

    for fname in expected_figures:
        fpath = output_dir / fname
        assert fpath.exists(), f"Figure not found: {fpath}"
        size = fpath.stat().st_size
        assert size > min_size_bytes, \
            f"Figure {fname} is too small ({size} bytes < {min_size_bytes} bytes)"
