"""Test suite for Case A0: Catalog Comparison Reference Report.

All tests run against the computed results JSON.
"""

import csv
import json
import os

import pytest

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_ROOT = os.path.dirname(BASE_DIR)
RESULTS_PATH = os.path.join(BASE_DIR, "output", "case-a0-results.json")
COMCAT_PATH = os.path.join(PROJECT_ROOT, "data", "global-sets", "comcat_global_6-9_1949-2021.csv")
ISCGEM_PATH = os.path.join(PROJECT_ROOT, "data", "global-sets", "iscgem_global_events.csv")

EXPECTED_COLUMNS = {
    "usgs_id", "usgs_mag", "event_at", "solaration_year", "solar_secs",
    "lunar_secs", "midnight_secs", "latitude", "longitude", "depth",
}

REQUIRED_TOP_KEYS = [
    "generated", "catalogs", "comcat_prefix_breakdown",
    "comcat_iscgem_prefix_temporal", "magnitude_precision", "magnitude_bins",
]


@pytest.fixture(scope="module")
def results() -> dict:
    """Load the results JSON once for all tests."""
    with open(RESULTS_PATH, "r") as f:
        return json.load(f)


def test_schema_validation_comcat() -> None:
    """Both CSV files contain exactly the 10 expected columns."""
    with open(COMCAT_PATH, "r", newline="") as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames or [])
    assert columns == EXPECTED_COLUMNS, f"ComCat columns mismatch: {columns}"


def test_schema_validation_iscgem() -> None:
    """Both CSV files contain exactly the 10 expected columns."""
    with open(ISCGEM_PATH, "r", newline="") as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames or [])
    assert columns == EXPECTED_COLUMNS, f"ISC-GEM columns mismatch: {columns}"


def test_comcat_event_count(results: dict) -> None:
    """ComCat event count equals 9,802."""
    assert results["catalogs"]["comcat"]["event_count"] == 9802


def test_iscgem_event_count(results: dict) -> None:
    """ISC-GEM event count equals 9,210."""
    assert results["catalogs"]["iscgem"]["event_count"] == 9210


def test_comcat_prefix_counts_sum(results: dict) -> None:
    """ComCat prefix counts sum to total event count."""
    breakdown = results["comcat_prefix_breakdown"]
    total = sum(breakdown[g]["count"] for g in ["us_native", "iscgem", "other"])
    assert total == 9802, f"Prefix counts sum to {total}, expected 9802"


def test_comcat_iscgem_prefix_count(results: dict) -> None:
    """ComCat iscgem prefix count equals 2,973 (+/-1 tolerance)."""
    count = results["comcat_prefix_breakdown"]["iscgem"]["count"]
    assert abs(count - 2973) <= 1, f"iscgem prefix count {count}, expected ~2973"


def test_magnitude_precision_sums(results: dict) -> None:
    """For each catalog, one_decimal + two_decimal == event_count."""
    for catalog in ["comcat", "iscgem"]:
        prec = results["magnitude_precision"][catalog]
        total = prec["one_decimal"]["count"] + prec["two_decimal"]["count"]
        expected = results["catalogs"][catalog]["event_count"]
        assert total == expected, (
            f"{catalog}: precision sum {total} != event_count {expected}"
        )


def test_magnitude_bin_coverage(results: dict) -> None:
    """All events accounted for in bins (sum of bin counts == event count)."""
    bins_data = results["magnitude_bins"]
    for catalog, key in [("comcat", "comcat_counts"), ("iscgem", "iscgem_counts")]:
        bin_sum = sum(bins_data[key])
        expected = results["catalogs"][catalog]["event_count"]
        assert bin_sum == expected, (
            f"{catalog}: bin sum {bin_sum} != event_count {expected}"
        )


def test_results_json_exists_and_valid(results: dict) -> None:
    """Results JSON exists, parses, and contains all required top-level keys."""
    assert os.path.isfile(RESULTS_PATH), f"Results file not found: {RESULTS_PATH}"
    for key in REQUIRED_TOP_KEYS:
        assert key in results, f"Missing top-level key: {key}"
