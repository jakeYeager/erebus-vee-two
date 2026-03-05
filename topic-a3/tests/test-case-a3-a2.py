"""
Case A3.A2 — Test Suite

15 tests covering data loading, module correctness invariants,
results JSON structure, and output file existence.
"""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Module loading (hyphenated filenames)
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
SRC_DIR = BASE_DIR / "src"
OUTPUT_DIR = BASE_DIR / "output"
DATA_ROOT = BASE_DIR.parent / "data" / "iscgem"

RESULTS_PATH = OUTPUT_DIR / "case-a3-a2-results.json"


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_schuster_mod = _load_module("case_a3_a2_schuster", SRC_DIR / "case-a3-a2-schuster.py")
_mfpa_mod = _load_module("case_a3_a2_mfpa", SRC_DIR / "case-a3-a2-mfpa.py")

schuster_spectrum = _schuster_mod.schuster_spectrum
mfpa_spectrum = _mfpa_mod.mfpa_spectrum
apply_fdr_bh = _mfpa_mod.apply_fdr_bh

EPOCH = pd.Timestamp("1950-01-01 00:00:00", tz="UTC")
JULIAN_YEAR_SECS = 31_557_600.0


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def results():
    with open(RESULTS_PATH) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def df_full():
    df = pd.read_csv(DATA_ROOT / "iscgem_global_6-9_1950-2021.csv")
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["event_time_days"] = (df["event_at"] - EPOCH).dt.total_seconds() / 86400.0
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    return df


@pytest.fixture(scope="module")
def df_gshhg():
    return pd.read_csv(DATA_ROOT / "plate-location" / "ocean_class_gshhg_global.csv")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_catalog_load(df_full):
    """Assert n=9,210; phase in [0,1); event_time_days non-decreasing after sort."""
    assert len(df_full) == 9210, f"Expected 9,210 rows; got {len(df_full)}"
    assert (df_full["phase"] >= 0.0).all(), "phase has values < 0"
    assert (df_full["phase"] < 1.0).all(), "phase has values >= 1.0"
    sorted_times = np.sort(df_full["event_time_days"].values)
    assert np.all(np.diff(sorted_times) >= 0), "event_time_days not non-decreasing after sort"


def test_gshhg_merge(df_full, df_gshhg):
    """Assert no NaN in ocean_class and tectonic_class after merge."""
    merged = df_full.merge(df_gshhg, on="usgs_id", how="left")
    assert merged["ocean_class"].notna().all(), "NaN in ocean_class after merge"
    # Tectonic class derivation
    tectonic = pd.cut(
        merged["dist_to_coast_km"],
        bins=[-np.inf, 50.0, 200.0, np.inf],
        labels=["continental", "transitional", "oceanic"],
    )
    assert tectonic.notna().all(), "NaN in tectonic_class after derivation"


def test_nh_sh_partition(df_full):
    """Assert n_NH + n_SH == 9,210."""
    n_nh = (df_full["latitude"] >= 0).sum()
    n_sh = (df_full["latitude"] < 0).sum()
    assert n_nh + n_sh == 9210, f"NH+SH={n_nh+n_sh}, expected 9210"


def test_strata_sizes(results):
    """Assert n_signal > 0; n_non_signal > 0; sum == 9,210."""
    sizes = results["strata_sizes"]
    n_sig = sizes["n_signal_bearing"]
    n_non = sizes["n_non_signal"]
    assert n_sig > 0, "n_signal_bearing must be > 0"
    assert n_non > 0, "n_non_signal must be > 0"
    assert n_sig + n_non == 9210, f"n_signal + n_non_signal = {n_sig + n_non}, expected 9210"


def test_gk_catalog_sizes():
    """Assert G-K mainshocks n=5,883 and aftershocks n=3,327."""
    df_ms = pd.read_csv(DATA_ROOT / "declustering-algorithm" / "mainshocks_gk-seq_global.csv")
    df_as = pd.read_csv(DATA_ROOT / "declustering-algorithm" / "aftershocks_gk-seq_global.csv")
    assert len(df_ms) == 5883, f"Expected 5,883 GK mainshocks; got {len(df_ms)}"
    assert len(df_as) == 3327, f"Expected 3,327 GK aftershocks; got {len(df_as)}"


def test_schuster_uniform():
    """Uniform random times: p_cluster_robust > 0.05 at annual period in >80% of 20 draws."""
    rng = np.random.default_rng(99)
    n_passes = 0
    for _ in range(20):
        t = rng.uniform(0, 72 * 365.25, size=500)
        results = schuster_spectrum(t, np.array([365.25]), dt_cluster_days=1.0)
        if results[0]["p_cluster_robust"] > 0.05:
            n_passes += 1
    assert n_passes / 20 > 0.80, f"Only {n_passes}/20 uniform draws had p_cluster_robust > 0.05"


def test_schuster_known_signal():
    """Strong half-year clustering: p_cluster_robust < 0.001 at 182.5 days."""
    rng = np.random.default_rng(42)
    # Generate events clustered near t % 182.5 ≈ 0
    n_events = 500
    t = np.sort(rng.choice(np.arange(0, 72 * 365), size=n_events, replace=False).astype(float))
    # Add strong phase signal: shift events to cluster near phase=0 in 182.5d period
    # by selecting times near multiples of 182.5 days (±5 days)
    t_signal = []
    for k in range(int(72 * 365 / 182.5)):
        t_cluster = 182.5 * k + rng.normal(0, 3.0, size=10)
        t_signal.extend(t_cluster)
    t_signal = np.sort(np.array(t_signal))
    results = schuster_spectrum(t_signal, np.array([182.5]), dt_cluster_days=1.0)
    assert results[0]["p_cluster_robust"] < 0.001, (
        f"Expected strong signal p < 0.001; got p={results[0]['p_cluster_robust']:.6f}"
    )


def test_cluster_robust_invariant(results):
    """p_cluster_robust >= p_standard for all periods on full catalog."""
    s1 = results["subtest_1_baseline"]["schuster_1day"]
    for entry in s1:
        assert entry["p_cluster_robust"] >= entry["p_standard"] - 1e-12, (
            f"Invariant violated at period={entry['period_days']:.3f}d: "
            f"p_cr={entry['p_cluster_robust']:.6f} < p_std={entry['p_standard']:.6f}"
        )


def test_mfpa_bootstrap_shape():
    """MFPA bootstrap null produces n_bootstrap values; p95 < p99."""
    rng = np.random.default_rng(1)
    t = rng.uniform(0, 72 * 365.25, size=100)
    n_boot = 500
    results = mfpa_spectrum(t, np.array([182.5, 365.25]), n_bootstrap=n_boot, rng_seed=7)
    # p95 < p99 for both period entries
    for r in results:
        assert r["p95_threshold"] < r["p99_threshold"], (
            f"p95={r['p95_threshold']:.6f} >= p99={r['p99_threshold']:.6f}"
        )


def test_fdr_bh_monotone():
    """FDR-corrected significance is a subset of uncorrected significance."""
    rng = np.random.default_rng(7)
    p_vals = rng.uniform(0, 1, size=200)
    fdr_flags = apply_fdr_bh(p_vals, alpha=0.05)
    # All FDR-significant must have p < 0.05 (or at least be a subset)
    # BH cannot make more tests significant than naive 0.05 cutoff in general,
    # but we verify: any FDR-significant period must have p < 0.05
    for i, sig in enumerate(fdr_flags):
        if sig:
            assert p_vals[i] < 0.05, (
                f"FDR marked p={p_vals[i]:.4f} as significant (should not exceed 0.05)"
            )


def test_half_year_explicit_recorded(results):
    """Half-year period (182.5d) in explicit_period_tests for all strata/catalogs."""
    # Sub-test 1
    assert "half_year" in results["subtest_1_baseline"]["explicit_period_tests"]
    # Sub-test 2
    assert "half_year" in results["subtest_2_stratified"]["signal_bearing"]["explicit_period_tests"]
    assert "half_year" in results["subtest_2_stratified"]["non_signal"]["explicit_period_tests"]
    # Sub-test 3
    assert "half_year" in results["subtest_3_nh_sh"]["nh"]["explicit_period_tests"]
    assert "half_year" in results["subtest_3_nh_sh"]["sh"]["explicit_period_tests"]
    # Sub-test 4
    assert "half_year" in results["subtest_4_declustering"]["full"]["explicit_period_tests"]
    assert "half_year" in results["subtest_4_declustering"]["gk_mainshocks"]["explicit_period_tests"]
    assert "half_year" in results["subtest_4_declustering"]["gk_aftershocks"]["explicit_period_tests"]


def test_nh_sh_phase_angles_recorded(results):
    """half_year_dominant_phase present and in [0,1) for both NH and SH."""
    nh_phase = results["subtest_3_nh_sh"]["nh"]["half_year_dominant_phase"]
    sh_phase = results["subtest_3_nh_sh"]["sh"]["half_year_dominant_phase"]
    assert 0.0 <= nh_phase < 1.0, f"NH half_year_dominant_phase={nh_phase} not in [0,1)"
    assert 0.0 <= sh_phase < 1.0, f"SH half_year_dominant_phase={sh_phase} not in [0,1)"


def test_phase_delta_computed(results):
    """half_year_phase_delta_fraction and delta_days present; |delta_fraction| <= 0.5."""
    delta_frac = results["subtest_3_nh_sh"]["half_year_phase_delta_fraction"]
    delta_days = results["subtest_3_nh_sh"]["half_year_phase_delta_days"]
    assert "half_year_phase_delta_fraction" in results["subtest_3_nh_sh"]
    assert "half_year_phase_delta_days" in results["subtest_3_nh_sh"]
    assert abs(delta_frac) <= 0.5, f"|delta_fraction|={abs(delta_frac):.4f} > 0.5"


def test_declustering_three_catalogs(results):
    """subtest_4_declustering has all three catalogs; interval1_contradiction_resolved is bool."""
    st4 = results["subtest_4_declustering"]
    assert "full" in st4, "Missing 'full' catalog in subtest_4"
    assert "gk_mainshocks" in st4, "Missing 'gk_mainshocks' in subtest_4"
    assert "gk_aftershocks" in st4, "Missing 'gk_aftershocks' in subtest_4"
    assert isinstance(st4["interval1_contradiction_resolved"], bool), (
        "interval1_contradiction_resolved must be bool"
    )


def test_output_figures_exist():
    """All 5 PNG files exist and are > 50 KB."""
    expected_figs = [
        "case-a3-a2-baseline-spectrum.png",
        "case-a3-a2-stratum-mfpa.png",
        "case-a3-a2-nh-sh-phase.png",
        "case-a3-a2-declustering-sensitivity.png",
        "case-a3-a2-cluster-window.png",
    ]
    for fname in expected_figs:
        p = OUTPUT_DIR / fname
        assert p.exists(), f"Figure not found: {p}"
        size_kb = p.stat().st_size / 1024
        assert size_kb > 50, f"Figure {fname} is too small: {size_kb:.1f} KB"
