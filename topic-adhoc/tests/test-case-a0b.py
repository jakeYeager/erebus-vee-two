"""
Test suite for Case A0b: Duplicate Detection and Cross-Catalog Event Accounting

Tests:
  1. Results JSON exists and is valid
  2. ID accounting: id_matched + id_unmatched == 2973
  3. Within-ComCat duplicate accounting: 0 <= duplicate_candidate_count <= 2973
  4. Net effective count: net_effective_comcat_count == 9802 - duplicate_candidate_count
  5. Cross-catalog primary accounting (ComCat): matched + comcat_only == 9802
  6. Cross-catalog primary accounting (ISC-GEM): matched + iscgem_only == 9210
  7. Strict tolerance <= primary: strict.matched <= primary.matched
  8. All three PNG images exist
"""

import json
from pathlib import Path

import pytest

BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a0b-results.json"


@pytest.fixture(scope="module")
def results() -> dict:
    """Load results JSON once for all tests."""
    assert RESULTS_PATH.exists(), f"Results file not found: {RESULTS_PATH}"
    with open(RESULTS_PATH) as f:
        data = json.load(f)
    return data


def test_results_json_exists_and_valid(results: dict) -> None:
    """Results JSON exists, parses as JSON, all top-level keys present."""
    expected_keys = {
        "generated",
        "prerequisites",
        "id_cross_reference",
        "within_comcat_duplicates",
        "cross_catalog_matching",
    }
    assert expected_keys.issubset(set(results.keys())), (
        f"Missing top-level keys: {expected_keys - set(results.keys())}"
    )


def test_id_accounting(results: dict) -> None:
    """id_matched + id_unmatched == 2973 (ComCat iscgem-prefix total)."""
    id_ref = results["id_cross_reference"]
    matched = id_ref["id_matched_in_iscgem"]["count"]
    unmatched = id_ref["id_unmatched_in_iscgem"]["count"]
    total = id_ref["comcat_iscgem_prefix_count"]
    assert matched + unmatched == total == 2973, (
        f"ID accounting failed: {matched} + {unmatched} = {matched + unmatched}, expected {total}"
    )


def test_within_comcat_duplicate_accounting(results: dict) -> None:
    """duplicate_candidate_count is a non-negative integer <= 2973."""
    dup_count = results["within_comcat_duplicates"]["duplicate_candidate_count"]
    assert isinstance(dup_count, int), f"duplicate_candidate_count is not int: {type(dup_count)}"
    assert 0 <= dup_count <= 2973, f"duplicate_candidate_count out of range: {dup_count}"


def test_net_effective_count(results: dict) -> None:
    """net_effective_comcat_count == 9802 - duplicate_candidate_count."""
    dups = results["within_comcat_duplicates"]
    expected = 9802 - dups["duplicate_candidate_count"]
    assert dups["net_effective_comcat_count"] == expected, (
        f"Net effective count: {dups['net_effective_comcat_count']} != {expected}"
    )


def test_cross_catalog_primary_comcat(results: dict) -> None:
    """Primary: matched + comcat_only == 9802."""
    primary = results["cross_catalog_matching"]["primary"]
    total = primary["matched"] + primary["comcat_only"]
    assert total == 9802, f"ComCat primary accounting: {total} != 9802"


def test_cross_catalog_primary_iscgem(results: dict) -> None:
    """Primary: matched + iscgem_only == 9210."""
    primary = results["cross_catalog_matching"]["primary"]
    total = primary["matched"] + primary["iscgem_only"]
    assert total == 9210, f"ISC-GEM primary accounting: {total} != 9210"


def test_strict_lte_primary(results: dict) -> None:
    """Strict tolerance matched count <= primary matched count."""
    primary_matched = results["cross_catalog_matching"]["primary"]["matched"]
    strict_matched = results["cross_catalog_matching"]["strict"]["matched"]
    assert strict_matched <= primary_matched, (
        f"Strict matched ({strict_matched}) > primary matched ({primary_matched})"
    )


def test_all_png_images_exist() -> None:
    """All three PNG images exist at expected paths."""
    expected_images = [
        "case-a0b-id-overlap.png",
        "case-a0b-event-accounting.png",
        "case-a0b-unmatched-temporal.png",
    ]
    for img in expected_images:
        path = OUTPUT_DIR / img
        assert path.exists(), f"Image not found: {path}"
