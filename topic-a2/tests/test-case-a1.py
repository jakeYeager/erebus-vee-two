"""
Case A1 — Test Suite.

Tests catalog loading, time conversion, Schuster spectrum, MFPA, and JSON output
structure as specified in the Case A1 spec.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
SRC_DIR = BASE_DIR / "src"
DATA_DIR = BASE_DIR.parent / "data"

RAW_PATH = DATA_DIR / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GK_PATH = DATA_DIR / "iscgem" / "declustering-algorithm" / "mainshocks_G-K_global.csv"
A1B_PATH = DATA_DIR / "iscgem" / "declustering-algorithm" / "mainshocks_a1b_global.csv"
RESULTS_PATH = BASE_DIR / "output" / "case-a1-results.json"

REFERENCE_EPOCH = pd.Timestamp("1950-01-01 00:00:00", tz="UTC")
JULIAN_YEAR_SECS = 31_557_600.0


# ---------------------------------------------------------------------------
# Module loaders
# ---------------------------------------------------------------------------

def _load_module(name: str, filename: str):
    """Dynamically load a hyphen-named module from src directory."""
    path = SRC_DIR / filename
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def schuster_mod():
    return _load_module("case_a1_schuster", "case-a1-schuster.py")


@pytest.fixture(scope="session")
def mfpa_mod():
    return _load_module("case_a1_mfpa", "case-a1-mfpa.py")


# ---------------------------------------------------------------------------
# Shared data fixture
# ---------------------------------------------------------------------------

def _load_catalog(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["event_at"])
    if df["event_at"].dt.tz is None:
        df["event_at"] = df["event_at"].dt.tz_localize("UTC")
    else:
        df["event_at"] = df["event_at"].dt.tz_convert("UTC")
    df.sort_values("event_at", inplace=True)
    df.reset_index(drop=True, inplace=True)
    delta = df["event_at"] - REFERENCE_EPOCH
    df["event_time_days"] = delta.dt.total_seconds() / 86400.0
    return df


@pytest.fixture(scope="session")
def catalogs():
    return {
        "raw": _load_catalog(RAW_PATH),
        "gk": _load_catalog(GK_PATH),
        "a1b": _load_catalog(A1B_PATH),
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestCatalogLoads:
    def test_catalog_loads(self, catalogs):
        """All three catalogs load with correct row counts."""
        assert len(catalogs["raw"]) == 9210, f"Raw: expected 9210, got {len(catalogs['raw'])}"
        assert len(catalogs["gk"]) == 5883, f"G-K: expected 5883, got {len(catalogs['gk'])}"
        assert len(catalogs["a1b"]) == 7137, f"A1b: expected 7137, got {len(catalogs['a1b'])}"


class TestEventTimeDays:
    def test_event_time_days_monotonic(self, catalogs):
        """event_time_days is non-decreasing after sort by event_at."""
        for name, df in catalogs.items():
            diffs = df["event_time_days"].diff().dropna()
            assert (diffs >= 0).all(), (
                f"{name}: event_time_days is not non-decreasing; "
                f"min diff = {diffs.min():.6f}"
            )


class TestSchusterUniform:
    def test_schuster_uniform(self, schuster_mod):
        """Standard Schuster p > 0.05 at annual period for uniform random data in >90% of trials."""
        rng = np.random.default_rng(42)
        year_days = 365.25
        span_days = 72 * year_days  # 72 years

        n_trials = 100
        n_events = 1000
        above_threshold = 0

        for _ in range(n_trials):
            t = np.sort(rng.uniform(0, span_days, size=n_events))
            result = schuster_mod.schuster_single_period(t, 365.25)
            if result["p_standard"] > 0.05:
                above_threshold += 1

        frac = above_threshold / n_trials
        assert frac > 0.90, (
            f"Standard Schuster p > 0.05 in only {frac*100:.1f}% of uniform trials "
            f"(expected >90%)"
        )


class TestSchusterAnnualSignal:
    def test_schuster_annual_signal(self, schuster_mod):
        """Strongly clustered annual signal yields p < 0.001 at 365.25 days."""
        rng = np.random.default_rng(123)
        year_days = 365.25
        span_days = 72 * year_days
        n_events = 1000

        # 90% of events near phase 0.2 (within ±0.02 of year), 10% uniform
        n_clustered = int(0.9 * n_events)
        n_uniform = n_events - n_clustered

        t_clustered = []
        for _ in range(n_clustered):
            year_start = rng.integers(0, 72) * year_days
            offset = rng.normal(0.2 * year_days, 0.02 * year_days)
            offset = np.clip(offset, 0, year_days)
            t_clustered.append(year_start + offset)

        t_uniform = rng.uniform(0, span_days, size=n_uniform)
        t_all = np.sort(np.concatenate([t_clustered, t_uniform]))

        result = schuster_mod.schuster_single_period(t_all, 365.25)
        assert result["p_standard"] < 0.001, (
            f"Expected p_standard < 0.001 for strongly annual signal, "
            f"got {result['p_standard']:.6f}"
        )


class TestClusterRobust:
    def test_cluster_robust_n_clusters(self, schuster_mod, catalogs):
        """n_clusters <= n_events for all catalogs and all spectrum periods."""
        # Check at a subset of periods for speed
        test_periods = [0.5, 1.0, 14.77, 182.625, 365.25]
        for cat_name, df in catalogs.items():
            t = df["event_time_days"].to_numpy()
            for T in test_periods:
                result = schuster_mod.schuster_single_period(t, T)
                assert result["n_clusters"] <= result["n_events"], (
                    f"{cat_name}, T={T}d: n_clusters={result['n_clusters']} > "
                    f"n_events={result['n_events']}"
                )

    def test_cluster_robust_p_ge_standard(self, schuster_mod):
        """cluster-robust p >= standard p for synthetic data with strong aftershock clustering.

        After-shock clustering inflates the standard Schuster test statistic by
        giving temporally clustered events multiple votes. The cluster-robust test
        collapses each temporal cluster (inter-event gap < 1 day) into a single
        representative, so for a synthetic catalog where clustering is the dominant
        source of phase coherence, the standard p-value should be lower (more
        significant) than the cluster-robust p-value.

        This test uses a synthetic catalog with explicit mainshock-aftershock
        structure: 200 mainshocks randomly distributed in time, each followed
        by a tight cluster of 5 aftershocks within hours. The cluster-robust test
        should return p >= standard p at the annual period.
        """
        rng = np.random.default_rng(42)
        year_days = 365.25
        n_mainshocks = 200
        n_aftershocks_per = 5

        # Mainshocks uniformly distributed over 72 years
        t_main = np.sort(rng.uniform(0, 72 * year_days, n_mainshocks))

        # Aftershocks: each within 0–12 hours (0–0.5 days) of their mainshock
        t_after = []
        for tm in t_main:
            offsets = rng.uniform(0.01, 0.5, n_aftershocks_per)
            t_after.extend(tm + offsets)

        t_all = np.sort(np.concatenate([t_main, t_after]))

        # At the annual period the aftershock clusters inflate the standard test
        result = schuster_mod.schuster_single_period(t_all, year_days)
        p_std = result["p_standard"]
        p_cr = result["p_cluster_robust"]
        # n_clusters should be much less than n_events
        assert result["n_clusters"] < result["n_events"], (
            "Synthetic catalog should have n_clusters < n_events with aftershock structure"
        )
        # cluster-robust p >= standard p: the standard test is inflated by aftershock clusters
        assert p_cr >= p_std - 1e-12, (
            f"cluster-robust p ({p_cr:.6f}) < standard p ({p_std:.6f}) for "
            f"synthetic catalog with {result['n_events']} events, "
            f"{result['n_clusters']} clusters"
        )


class TestMFPASpectrum:
    def test_mfpa_spectrum_length(self, mfpa_mod):
        """MFPA spectrum has exactly 300 entries."""
        rng = np.random.default_rng(7)
        t = np.sort(rng.uniform(0, 72 * 365.25, size=500))
        result = mfpa_mod.mfpa_scan(t)
        assert len(result["spectrum"]) == 300, (
            f"Expected 300 MFPA spectrum entries, got {len(result['spectrum'])}"
        )

    def test_mfpa_period_range(self, mfpa_mod):
        """Min MFPA period >= 0.25 days and max <= 548 days."""
        rng = np.random.default_rng(8)
        t = np.sort(rng.uniform(0, 72 * 365.25, size=500))
        result = mfpa_mod.mfpa_scan(t)
        periods = [e["period_days"] for e in result["spectrum"]]
        assert min(periods) >= 0.25, f"Min period {min(periods):.4f} < 0.25"
        assert max(periods) <= 548.0, f"Max period {max(periods):.4f} > 548"

    def test_mfpa_power_positive(self, mfpa_mod):
        """All MFPA power values >= 0."""
        rng = np.random.default_rng(9)
        t = np.sort(rng.uniform(0, 72 * 365.25, size=500))
        result = mfpa_mod.mfpa_scan(t)
        for entry in result["spectrum"]:
            assert entry["power"] >= 0, (
                f"Negative power {entry['power']} at period {entry['period_days']:.2f}d"
            )

    def test_a1b_crossref_format(self, mfpa_mod):
        """All MFPA spectrum entries have an 'a1b_consistency' field."""
        rng = np.random.default_rng(11)
        t = np.sort(rng.uniform(0, 72 * 365.25, size=500))
        result = mfpa_mod.mfpa_scan(t)
        for entry in result["spectrum"]:
            assert "a1b_consistency" in entry, (
                f"Entry for period {entry.get('period_days', '?')} missing 'a1b_consistency'"
            )
            assert isinstance(entry["a1b_consistency"], str), (
                "a1b_consistency should be a string"
            )


class TestResultsJSON:
    def test_results_json_structure(self):
        """Results JSON has correct top-level structure."""
        assert RESULTS_PATH.exists(), f"Results JSON not found at {RESULTS_PATH}"
        with open(RESULTS_PATH) as fh:
            data = json.load(fh)

        assert "schuster" in data, "Key 'schuster' missing from results JSON"
        assert "mfpa" in data, "Key 'mfpa' missing from results JSON"

        for section in ["schuster", "mfpa"]:
            for sub in ["raw", "gk_mainshocks", "a1b_mainshocks"]:
                assert sub in data[section], (
                    f"Key '{sub}' missing from results JSON section '{section}'"
                )
