"""
Test suite for Case B5: Solar Declination Rate-of-Change vs. Position Test.

Tests cover data loading, column ranges, bin assignments, statistical
outputs, and JSON structure as specified in the case spec section 5.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
SOLAR_GEO_PATH = (
    BASE_DIR.parent / "data" / "iscgem" / "celestial-geometry" / "solar_geometry_global.csv"
)
RESULTS_PATH = BASE_DIR / "output" / "case-b5-results.json"

JULIAN_YEAR_SECS = 31_557_600.0
K24 = 24


# ── Fixtures ───────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def df() -> pd.DataFrame:
    """Load the solar geometry catalog once for the test module."""
    return pd.read_csv(SOLAR_GEO_PATH)


@pytest.fixture(scope="module")
def results() -> dict:
    """Load the case B5 results JSON once for the test module."""
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


# ── Tests ──────────────────────────────────────────────────────────────────────

class TestDataLoading:
    """Tests for solar geometry catalog loading and schema."""

    def test_solar_geo_load(self, df: pd.DataFrame) -> None:
        """Assert n=9210 and all four required columns present."""
        assert len(df) == 9210, f"Expected 9210 rows, got {len(df)}"
        required_cols = [
            "solar_secs",
            "solar_declination",
            "declination_rate",
            "earth_sun_distance",
        ]
        for col in required_cols:
            assert col in df.columns, f"Missing required column: {col}"

    def test_column_ranges_approx(self, df: pd.DataFrame) -> None:
        """Assert variable ranges within tolerance of REQ-3 spec ranges."""
        # solar_declination within ±2° of [-23.5, +23.5]
        assert df["solar_declination"].min() >= -23.5 - 2.0, (
            f"solar_declination min {df['solar_declination'].min()} too low"
        )
        assert df["solar_declination"].max() <= 23.5 + 2.0, (
            f"solar_declination max {df['solar_declination'].max()} too high"
        )

        # declination_rate within ±0.05 deg/day of [-0.40, +0.40]
        assert df["declination_rate"].min() >= -0.40 - 0.05, (
            f"declination_rate min {df['declination_rate'].min()} too low"
        )
        assert df["declination_rate"].max() <= 0.40 + 0.05, (
            f"declination_rate max {df['declination_rate'].max()} too high"
        )

        # earth_sun_distance within ±0.002 AU of [0.983, 1.017]
        assert df["earth_sun_distance"].min() >= 0.983 - 0.002, (
            f"earth_sun_distance min {df['earth_sun_distance'].min()} too low"
        )
        assert df["earth_sun_distance"].max() <= 1.017 + 0.002, (
            f"earth_sun_distance max {df['earth_sun_distance'].max()} too high"
        )


class TestPhaseAndBins:
    """Tests for phase computation and bin assignments."""

    def test_phase_range(self, df: pd.DataFrame) -> None:
        """Assert solar_phase is in [0.0, 1.0)."""
        solar_phase = (df["solar_secs"].values / JULIAN_YEAR_SECS) % 1.0
        assert float(solar_phase.min()) >= 0.0, f"solar_phase min {solar_phase.min()} < 0"
        assert float(solar_phase.max()) < 1.0, f"solar_phase max {solar_phase.max()} >= 1.0"

    def test_non_cyclic_bins_in_range(self, df: pd.DataFrame) -> None:
        """Assert all non-cyclic bin indices are in [0, k-1] for k=24."""
        k = K24

        # Declination
        dec_vals = df["solar_declination"].values
        dec_min, dec_max = dec_vals.min(), dec_vals.max()
        dec_range = dec_max - dec_min
        bins_dec = np.clip(
            np.floor((dec_vals - dec_min) / dec_range * k).astype(int), 0, k - 1
        )
        assert int(bins_dec.min()) >= 0, "bin_declination below 0"
        assert int(bins_dec.max()) <= k - 1, f"bin_declination above {k - 1}"

        # Rate
        rate_vals = df["declination_rate"].values
        rate_min, rate_max = rate_vals.min(), rate_vals.max()
        rate_range = rate_max - rate_min
        bins_rate = np.clip(
            np.floor((rate_vals - rate_min) / rate_range * k).astype(int), 0, k - 1
        )
        assert int(bins_rate.min()) >= 0, "bin_rate below 0"
        assert int(bins_rate.max()) <= k - 1, f"bin_rate above {k - 1}"

        # Distance
        dist_vals = df["earth_sun_distance"].values
        dist_min, dist_max = dist_vals.min(), dist_vals.max()
        dist_range = dist_max - dist_min
        bins_dist = np.clip(
            np.floor((dist_vals - dist_min) / dist_range * k).astype(int), 0, k - 1
        )
        assert int(bins_dist.min()) >= 0, "bin_distance below 0"
        assert int(bins_dist.max()) <= k - 1, f"bin_distance above {k - 1}"


class TestStatisticalOutputs:
    """Tests for statistical measures in the results JSON."""

    ALL_VARIABLES = ["solar_phase", "solar_declination", "declination_rate", "earth_sun_distance"]
    ALL_K_KEYS = ["k16", "k24", "k32"]

    def test_chi_square_all_variables(self, results: dict) -> None:
        """Assert chi2 and p_chi2 are finite for all four variables at k=16, 24, 32."""
        var_stats = results["variable_stats"]
        for var in self.ALL_VARIABLES:
            for k_key in self.ALL_K_KEYS:
                block = var_stats[var][k_key]
                assert math.isfinite(block["chi2"]), (
                    f"{var} {k_key} chi2 is not finite: {block['chi2']}"
                )
                assert math.isfinite(block["p_chi2"]) or block["p_chi2"] == 0.0, (
                    f"{var} {k_key} p_chi2 is not finite: {block['p_chi2']}"
                )
                assert block["chi2"] >= 0.0, f"{var} {k_key} chi2 < 0"
                assert 0.0 <= block["p_chi2"] <= 1.0, (
                    f"{var} {k_key} p_chi2 {block['p_chi2']} not in [0, 1]"
                )

    def test_cramer_v_range(self, results: dict) -> None:
        """Assert Cramér's V is in [0.0, 1.0] for all variables and k values."""
        var_stats = results["variable_stats"]
        for var in self.ALL_VARIABLES:
            for k_key in self.ALL_K_KEYS:
                v = var_stats[var][k_key]["cramer_v"]
                assert 0.0 <= v <= 1.0, (
                    f"{var} {k_key} cramer_v {v} not in [0, 1]"
                )

    def test_rayleigh_solar_phase_only(self, results: dict) -> None:
        """Assert rayleigh_R is present for solar_phase but not for other variables."""
        var_stats = results["variable_stats"]

        # solar_phase must have rayleigh_R
        for k_key in self.ALL_K_KEYS:
            assert "rayleigh_R" in var_stats["solar_phase"][k_key], (
                f"solar_phase {k_key} missing rayleigh_R"
            )

        # Other variables must NOT have rayleigh_R
        for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
            for k_key in self.ALL_K_KEYS:
                assert "rayleigh_R" not in var_stats[var][k_key], (
                    f"{var} {k_key} should not have rayleigh_R"
                )

    def test_ks_non_cyclic_only(self, results: dict) -> None:
        """Assert ks_stat is present for non-cyclic variables but not solar_phase."""
        var_stats = results["variable_stats"]

        # Non-cyclic variables must have ks_stat
        for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
            for k_key in self.ALL_K_KEYS:
                assert "ks_stat" in var_stats[var][k_key], (
                    f"{var} {k_key} missing ks_stat"
                )

        # solar_phase must NOT have ks_stat
        for k_key in self.ALL_K_KEYS:
            assert "ks_stat" not in var_stats["solar_phase"][k_key], (
                f"solar_phase {k_key} should not have ks_stat"
            )


class TestRankingAndAlignment:
    """Tests for variable ranking and A1b alignment structure."""

    ALL_VARIABLES = ["solar_phase", "solar_declination", "declination_rate", "earth_sun_distance"]

    def test_variable_ranking_complete(self, results: dict) -> None:
        """Assert variable_ranking has exactly 4 entries covering all four variable names."""
        ranking = results["variable_ranking"]
        assert len(ranking) == 4, f"Expected 4 ranking entries, got {len(ranking)}"
        ranked_vars = {r["variable"] for r in ranking}
        expected_vars = set(self.ALL_VARIABLES)
        assert ranked_vars == expected_vars, (
            f"Ranked variables {ranked_vars} != expected {expected_vars}"
        )

    def test_most_significant_valid(self, results: dict) -> None:
        """Assert most_significant_variable is one of the four variable names."""
        msv = results["most_significant_variable"]
        assert msv in self.ALL_VARIABLES, (
            f"most_significant_variable '{msv}' not in {self.ALL_VARIABLES}"
        )

    def test_a1b_alignment_keys(self, results: dict) -> None:
        """Assert a1b_alignment has keys for all three intervals."""
        alignment = results["a1b_alignment"]
        expected_keys = {
            "interval_1_phase_0.22",
            "interval_2_phase_0.64",
            "interval_3_phase_0.90",
        }
        assert set(alignment.keys()) == expected_keys, (
            f"a1b_alignment keys {set(alignment.keys())} != expected {expected_keys}"
        )
        # Each interval entry must have a variable_with_elevated_bin list
        for key in expected_keys:
            assert "variable_with_elevated_bin" in alignment[key], (
                f"a1b_alignment[{key}] missing 'variable_with_elevated_bin'"
            )
            assert isinstance(alignment[key]["variable_with_elevated_bin"], list)
