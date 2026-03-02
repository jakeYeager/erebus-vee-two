"""Case A4: Declustering Sensitivity Analysis — Main Analysis Script.

Loads all seven ISC-GEM catalogs (raw + three pairs of mainshock/aftershock
files), runs Sub-analyses A, B, and C, and writes the combined results to
output/case-a4-results.json.
"""

import json
import logging
import sys
from pathlib import Path
from typing import Any

import importlib.util
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
DATA_DIR = BASE_DIR.parent / "data" / "iscgem"
DECLUSTER_DIR = DATA_DIR / "declustering-algorithm"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SRC_DIR = BASE_DIR / "src"


def _import_hyphenated(name: str, filepath: Path):
    """Import a module from a hyphenated filename.

    Args:
        name: Module alias to register.
        filepath: Absolute path to the .py file.

    Returns:
        Loaded module object.
    """
    spec = importlib.util.spec_from_file_location(name, filepath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_sub_a_mod = _import_hyphenated("case_a4_sub_a", SRC_DIR / "case-a4-sub-a.py")
_sub_b_mod = _import_hyphenated("case_a4_sub_b", SRC_DIR / "case-a4-sub-b.py")
_sub_c_mod = _import_hyphenated("case_a4_sub_c", SRC_DIR / "case-a4-sub-c.py")

run_sub_a = _sub_a_mod.run_sub_a
run_sub_b = _sub_b_mod.run_sub_b
run_sub_c = _sub_c_mod.run_sub_c

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("case-a4-analysis")

# ---------------------------------------------------------------------------
# Julian constant (confirmed: use uniformly for all solar year normalization)
# ---------------------------------------------------------------------------
SOLAR_YEAR_SECS = 31_557_600.0  # 365.25 * 86400 seconds

# ---------------------------------------------------------------------------
# Expected row counts from spec
# ---------------------------------------------------------------------------
EXPECTED_COUNTS: dict[str, int] = {
    "raw": 9210,
    "gk_mainshocks": 5883,
    "gk_aftershocks": 3327,
    "reas_mainshocks": 8265,
    "reas_aftershocks": 945,
    "a1b_mainshocks": 7137,
    "a1b_aftershocks": 2073,
}


def load_catalog(path: Path, label: str) -> pd.DataFrame:
    """Load a CSV catalog and validate required columns.

    Args:
        path: Path to CSV file.
        label: Human-readable label for logging.

    Returns:
        Loaded DataFrame.
    """
    logger.info("Loading %s from %s", label, path)
    df = pd.read_csv(path)
    logger.info("%s: %d rows loaded", label, len(df))
    return df


def validate_counts(catalogs: dict[str, pd.DataFrame]) -> None:
    """Assert row counts match expected values and log discrepancies.

    Args:
        catalogs: Dict mapping catalog key → DataFrame.
    """
    for key, expected in EXPECTED_COUNTS.items():
        actual = len(catalogs[key])
        if actual != expected:
            logger.warning(
                "Count mismatch for %s: expected %d, got %d", key, expected, actual
            )
        else:
            logger.info("Count OK: %s = %d", key, actual)


def validate_partition_integrity(catalogs: dict[str, pd.DataFrame]) -> None:
    """Check that mainshocks + aftershocks = raw for each method, with no usgs_id overlap.

    Args:
        catalogs: Full catalog dict.
    """
    pairs = [
        ("gk_mainshocks", "gk_aftershocks"),
        ("reas_mainshocks", "reas_aftershocks"),
        ("a1b_mainshocks", "a1b_aftershocks"),
    ]
    total_raw = len(catalogs["raw"])
    for ms_key, as_key in pairs:
        n_ms = len(catalogs[ms_key])
        n_as = len(catalogs[as_key])
        total = n_ms + n_as
        if total != total_raw:
            logger.warning(
                "Partition mismatch for %s+%s: %d + %d = %d (expected %d)",
                ms_key, as_key, n_ms, n_as, total, total_raw,
            )
        else:
            logger.info("Partition OK: %s + %s = %d", ms_key, as_key, total)

        # Check for usgs_id overlap
        ms_ids = set(catalogs[ms_key]["usgs_id"].tolist())
        as_ids = set(catalogs[as_key]["usgs_id"].tolist())
        overlap = ms_ids & as_ids
        if overlap:
            logger.warning(
                "usgs_id overlap between %s and %s: %d events", ms_key, as_key, len(overlap)
            )
        else:
            logger.info("No usgs_id overlap between %s and %s", ms_key, as_key)


def main() -> None:
    """Load all catalogs, run all sub-analyses, write results JSON."""
    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    RAW_PATH = DATA_DIR / "iscgem_global_6-9_1950-2021.csv"

    catalogs: dict[str, pd.DataFrame] = {
        "raw": load_catalog(RAW_PATH, "raw"),
        "gk_mainshocks": load_catalog(DECLUSTER_DIR / "mainshocks_G-K_global.csv", "G-K mainshocks"),
        "gk_aftershocks": load_catalog(DECLUSTER_DIR / "aftershocks_G-K_global.csv", "G-K aftershocks"),
        "reas_mainshocks": load_catalog(DECLUSTER_DIR / "mainshocks_reas_global.csv", "Reasenberg mainshocks"),
        "reas_aftershocks": load_catalog(DECLUSTER_DIR / "aftershocks_reas_global.csv", "Reasenberg aftershocks"),
        "a1b_mainshocks": load_catalog(DECLUSTER_DIR / "mainshocks_a1b_global.csv", "A1b mainshocks"),
        "a1b_aftershocks": load_catalog(DECLUSTER_DIR / "aftershocks_a1b_global.csv", "A1b aftershocks"),
    }

    validate_counts(catalogs)
    validate_partition_integrity(catalogs)

    # ------------------------------------------------------------------
    # Sub-analysis A
    # ------------------------------------------------------------------
    logger.info("--- Sub-analysis A: Scalar signal survival ---")
    mainshock_catalogs_a = {
        "raw": catalogs["raw"],
        "gk_mainshocks": catalogs["gk_mainshocks"],
        "reas_mainshocks": catalogs["reas_mainshocks"],
        "a1b_mainshocks": catalogs["a1b_mainshocks"],
    }
    sub_a_results = run_sub_a(mainshock_catalogs_a)

    # ------------------------------------------------------------------
    # Sub-analysis B
    # ------------------------------------------------------------------
    logger.info("--- Sub-analysis B: Interval structure ---")
    mainshock_catalogs_b = {
        "gk_mainshocks": catalogs["gk_mainshocks"],
        "reas_mainshocks": catalogs["reas_mainshocks"],
        "a1b_mainshocks": catalogs["a1b_mainshocks"],
    }
    sub_b_results = run_sub_b(mainshock_catalogs_b)

    # ------------------------------------------------------------------
    # Sub-analysis C
    # ------------------------------------------------------------------
    logger.info("--- Sub-analysis C: Aftershock phase preference ---")
    aftershock_catalogs = {
        "gk_aftershocks": catalogs["gk_aftershocks"],
        "reas_aftershocks": catalogs["reas_aftershocks"],
        "a1b_aftershocks": catalogs["a1b_aftershocks"],
    }
    sub_c_results = run_sub_c(aftershock_catalogs, sub_b_results)

    # ------------------------------------------------------------------
    # Assemble and write results JSON
    # ------------------------------------------------------------------
    results: dict[str, Any] = {
        "case": "A4",
        "title": "Declustering Sensitivity Analysis",
        "catalog_counts": {k: len(v) for k, v in catalogs.items()},
        "solar_year_secs": SOLAR_YEAR_SECS,
        "sub_a": sub_a_results,
        "sub_b": sub_b_results,
        "sub_c": sub_c_results,
    }

    output_path = OUTPUT_DIR / "case-a4-results.json"
    with open(output_path, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", output_path)


if __name__ == "__main__":
    main()
