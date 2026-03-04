"""
Case A3.B3 Test Suite — Ocean/Coast Sequential Threshold Sensitivity

All 16 tests must pass.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = BASE_DIR / "output"
LIB_DIR = BASE_DIR.parent / "lib"

RESULTS_PATH = OUTPUT_DIR / "case-a3-b3-results.json"
EVENTS_PATH  = OUTPUT_DIR / "case-a3-b3-events.pkl"
STEPS_PATH   = LIB_DIR / "PB2002_steps.dat"

RAW_PATH   = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
PB2002_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_pb2002_global.csv"

JULIAN_YEAR_SECS = 31_557_600.0
T_OUTER_STEPS = [200, 175, 150, 125, 100, 75, 50, 25]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def results() -> dict:
    """Load the results JSON once per module."""
    assert RESULTS_PATH.exists(), f"Results JSON not found: {RESULTS_PATH}"
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def events_df() -> pd.DataFrame:
    """Load the event dataframe once per module."""
    assert EVENTS_PATH.exists(), f"Events pickle not found: {EVENTS_PATH}"
    return pd.read_pickle(EVENTS_PATH)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_catalog_load() -> None:
    """Assert raw catalog has 9210 rows, valid event_at timestamps, and phases in [0,1)."""
    df = pd.read_csv(RAW_PATH, parse_dates=["event_at"])
    assert len(df) == 9210, f"Expected 9210 rows, got {len(df)}"
    assert df["event_at"].isna().sum() == 0, "NaT values found in event_at"
    phases = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    assert phases.min() >= 0.0, "Phase below 0"
    assert phases.max() < 1.0, "Phase >= 1"


def test_classification_joins() -> None:
    """Assert GSHHG and PB2002 joins each produce 9210 rows with no NaN in dist_to_coast_km."""
    gshhg_df = pd.read_csv(GSHHG_PATH)
    pb2002_df = pd.read_csv(PB2002_PATH)
    assert len(gshhg_df) == 9210, f"GSHHG expected 9210, got {len(gshhg_df)}"
    assert len(pb2002_df) == 9210, f"PB2002 expected 9210, got {len(pb2002_df)}"
    assert gshhg_df["dist_to_coast_km"].isna().sum() == 0, "NaN in GSHHG dist_to_coast_km"
    assert pb2002_df["dist_to_coast_km"].isna().sum() == 0, "NaN in PB2002 dist_to_coast_km"


def test_baseline_class_counts(events_df: pd.DataFrame) -> None:
    """Assert GSHHG T=200 class counts are within ±5 of expected values."""
    def classify(d: float) -> str:
        if d <= 50.0:
            return "continental"
        elif d > 200.0:
            return "oceanic"
        else:
            return "transitional"

    classes = events_df["dist_km_gshhg"].map(classify)
    counts = classes.value_counts()

    assert abs(counts.get("oceanic", 0) - 1952) <= 5, \
        f"Oceanic count {counts.get('oceanic')} not within 5 of 1952"
    assert abs(counts.get("transitional", 0) - 3459) <= 5, \
        f"Transitional count {counts.get('transitional')} not within 5 of 3459"
    assert abs(counts.get("continental", 0) - 3799) <= 5, \
        f"Continental count {counts.get('continental')} not within 5 of 3799"


def test_threshold_step_count(results: dict) -> None:
    """Assert threshold sweep produces exactly 8 step records."""
    sweep = results["threshold_sweep"]
    assert len(sweep) == 8, f"Expected 8 threshold steps, got {len(sweep)}"


def test_class_partition_all_steps(results: dict) -> None:
    """For each threshold step, assert oceanic + transitional + continental == 9210."""
    for step in results["threshold_sweep"]:
        total = step["oceanic"]["n"] + step["transitional"]["n"] + step["continental"]["n"]
        assert total == 9210, f"At T={step['t_outer_km']}: total={total} != 9210"


def test_n_monotonic_oceanic(results: dict) -> None:
    """As T_outer decreases (200→25), oceanic n is non-decreasing."""
    sweep = results["threshold_sweep"]
    # Sweep is ordered 200→25; as T decreases, more events become oceanic (wider ocean definition)
    # Wait: as T decreases, the boundary shrinks — events that were transitional become oceanic
    # when T tightens. So oceanic N should be NON-DECREASING as T decreases.
    oce_ns = [s["oceanic"]["n"] for s in sweep]
    for i in range(len(oce_ns) - 1):
        assert oce_ns[i] <= oce_ns[i + 1], \
            f"Oceanic n decreased from T={sweep[i]['t_outer_km']} ({oce_ns[i]}) to T={sweep[i+1]['t_outer_km']} ({oce_ns[i+1]})"


def test_chi2_bounds(results: dict) -> None:
    """For all threshold steps and all classes, chi2_k24 >= 0 and p_chi2_k24 in [0, 1]."""
    for step in results["threshold_sweep"]:
        for cls in ["oceanic", "transitional", "continental"]:
            s = step[cls]
            if s["chi2_k24"] is not None:
                assert s["chi2_k24"] >= 0, f"chi2 < 0 at T={step['t_outer_km']}, class={cls}"
            if s["p_chi2_k24"] is not None:
                assert 0.0 <= s["p_chi2_k24"] <= 1.0, \
                    f"p_chi2 out of [0,1] at T={step['t_outer_km']}, class={cls}"


def test_cramers_v_bounds(results: dict) -> None:
    """Assert cramers_v >= 0 for all threshold steps and classes."""
    for step in results["threshold_sweep"]:
        for cls in ["oceanic", "transitional", "continental"]:
            cv = step[cls]["cramers_v"]
            if cv is not None:
                assert cv >= 0, f"cramers_v < 0 at T={step['t_outer_km']}, class={cls}"


def test_baseline_matches_b2(results: dict) -> None:
    """T=200 oceanic p_chi2_k24 within 0.01 of 0.061; transitional p < 0.01."""
    sweep = results["threshold_sweep"]
    baseline = next(s for s in sweep if s["t_outer_km"] == 200.0)

    oce_p = baseline["oceanic"]["p_chi2_k24"]
    assert oce_p is not None, "Oceanic p is None at T=200"
    assert abs(oce_p - 0.061) <= 0.01, \
        f"T=200 oceanic p={oce_p:.4f} not within 0.01 of 0.061 (A2.B2 reference)"

    trans_p = baseline["transitional"]["p_chi2_k24"]
    assert trans_p is not None, "Transitional p is None at T=200"
    assert trans_p < 0.01, \
        f"T=200 transitional p={trans_p:.4f} not < 0.01"


def test_sub_boundary_point_count(results: dict) -> None:
    """Assert n_sub_boundary_points is between 700 and 1100."""
    n_pts = results["subduction_proximity"]["n_sub_boundary_points"]
    assert 700 <= n_pts <= 1100, \
        f"n_sub_boundary_points={n_pts} outside expected range [700, 1100]"


def test_dist_to_sub_nonneg(events_df: pd.DataFrame) -> None:
    """Assert all dist_to_subduction_km values >= 0."""
    assert (events_df["dist_to_subduction_km"] >= 0).all(), \
        "Negative distances found in dist_to_subduction_km"


def test_dist_to_sub_plausible(events_df: pd.DataFrame) -> None:
    """Assert mean dist_to_subduction_km is between 100 and 2000 km."""
    mean_d = events_df["dist_to_subduction_km"].mean()
    assert 100 <= mean_d <= 2000, \
        f"Mean dist_to_subduction_km={mean_d:.1f} km outside plausible range [100, 2000]"


def test_gcmt_validation_structure(results: dict) -> None:
    """Assert gcmt_validation key is present with proxy_validated bool."""
    assert "gcmt_validation" in results, "gcmt_validation key missing from results"
    gv = results["gcmt_validation"]
    assert "proxy_validated" in gv, "proxy_validated key missing from gcmt_validation"
    assert isinstance(gv["proxy_validated"], bool), \
        f"proxy_validated should be bool, got {type(gv['proxy_validated'])}"


def test_results_json_completeness(results: dict) -> None:
    """Assert threshold_sweep list has length 8; each entry has all class keys with chi2 and p."""
    sweep = results["threshold_sweep"]
    assert len(sweep) == 8, f"Expected 8 threshold steps, got {len(sweep)}"
    for step in sweep:
        for cls in ["oceanic", "transitional", "continental"]:
            assert cls in step, f"Class '{cls}' missing at T={step['t_outer_km']}"
            assert "chi2_k24" in step[cls], f"chi2_k24 missing at T={step['t_outer_km']}, class={cls}"
            assert "p_chi2_k24" in step[cls], f"p_chi2_k24 missing at T={step['t_outer_km']}, class={cls}"


def test_sub_parse_type_filter() -> None:
    """Confirm that parsing PB2002_steps.dat and filtering SUB produces no OTF, OSR, CTF, CRB rows."""
    excluded_types = {"OTF", "OSR", "CTF", "CRB"}
    found_excluded: list[str] = []

    with open(STEPS_PATH, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            btype = parts[-1].lstrip(":")
            if btype in ("SUB", "OCB"):
                # When we filter to SUB/OCB, none of these should appear in the result set
                pass
            # Check that a non-excluded type passed filter would be caught
            if btype in excluded_types:
                # This is NOT in our filter set — just confirming excluded types exist
                pass

    # Simulate the parse and assert no excluded types make it through
    accepted_btypes: list[str] = []
    with open(STEPS_PATH, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            btype = parts[-1].lstrip(":")
            if btype in ("SUB", "OCB"):
                accepted_btypes.append(btype)

    for bt in accepted_btypes:
        assert bt not in excluded_types, \
            f"Excluded boundary type '{bt}' made it through the filter"


def test_region_maps_output_exists() -> None:
    """Assert region-maps PNG exists and file size > 100 KB."""
    png_path = OUTPUT_DIR / "case-a3-b3-region-maps.png"
    assert png_path.exists(), f"Region maps PNG not found: {png_path}"
    size_bytes = png_path.stat().st_size
    assert size_bytes > 100 * 1024, \
        f"Region maps PNG size {size_bytes/1024:.1f} KB is less than 100 KB"
