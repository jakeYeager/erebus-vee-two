"""Case A0 Analysis: Catalog Comparison Reference Report

Loads ComCat and ISC-GEM CSVs, computes descriptive statistics,
and writes case-a0-results.json. No visualizations.
"""

import csv
import json
import logging
import os
from datetime import datetime, timezone
from typing import Any

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_ROOT = os.path.dirname(BASE_DIR)

COMCAT_PATH = os.path.join(PROJECT_ROOT, "data", "global-sets", "comcat_global_6-9_1949-2021.csv")
ISCGEM_PATH = os.path.join(PROJECT_ROOT, "data", "global-sets", "iscgem_global_events.csv")
OUTPUT_PATH = os.path.join(BASE_DIR, "output", "case-a0-results.json")

EXPECTED_COLUMNS = {
    "usgs_id", "usgs_mag", "event_at", "solaration_year", "solar_secs",
    "lunar_secs", "midnight_secs", "latitude", "longitude", "depth",
}


def load_csv(path: str) -> list[dict[str, str]]:
    """Load a CSV file and return list of row dicts."""
    logger.info("Loading %s", path)
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    logger.info("Loaded %d rows from %s", len(rows), os.path.basename(path))
    return rows


def population_summary(rows: list[dict[str, str]], label: str) -> dict[str, Any]:
    """Compute event count, year range, and magnitude range."""
    years = [int(float(r["solaration_year"])) for r in rows]
    mags = [float(r["usgs_mag"]) for r in rows]
    summary = {
        "event_count": len(rows),
        "year_range": {"min": min(years), "max": max(years)},
        "mag_range": {"min": round(min(mags), 2), "max": round(max(mags), 2)},
    }
    logger.info("%s population: %d events, years %d-%d, mag %.2f-%.2f",
                label, summary["event_count"],
                summary["year_range"]["min"], summary["year_range"]["max"],
                summary["mag_range"]["min"], summary["mag_range"]["max"])
    return summary


def classify_prefix(usgs_id: str) -> str:
    """Classify a ComCat usgs_id into us_native, iscgem, or other."""
    if usgs_id.startswith("iscgem"):
        return "iscgem"
    elif usgs_id.startswith("us"):
        return "us_native"
    else:
        return "other"


def comcat_prefix_breakdown(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    """Count ComCat rows by ID prefix group."""
    counts: dict[str, int] = {"us_native": 0, "iscgem": 0, "other": 0}
    total = len(rows)
    for r in rows:
        group = classify_prefix(r["usgs_id"])
        counts[group] += 1

    result = {}
    for group in ["us_native", "iscgem", "other"]:
        result[group] = {
            "count": counts[group],
            "pct": round(counts[group] / total * 100, 1),
        }
    logger.info("ComCat prefix breakdown: %s", result)
    return result


def iscgem_prefix_temporal(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    """Temporal distribution of iscgem-prefixed ComCat records by decade."""
    iscgem_rows = [r for r in rows if r["usgs_id"].startswith("iscgem")]
    total = len(iscgem_rows)
    logger.info("ISC-GEM prefixed ComCat records: %d", total)

    decades = ["1940s", "1950s", "1960s", "1970s", "1980s", "1990s", "2000s", "2010s", "2020s"]
    decade_starts = [1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020]

    decade_counts: dict[str, int] = {d: 0 for d in decades}
    for r in iscgem_rows:
        year = int(float(r["solaration_year"]))
        for i in range(len(decade_starts) - 1, -1, -1):
            if year >= decade_starts[i]:
                decade_counts[decades[i]] += 1
                break

    result = {}
    for d in decades:
        result[d] = {
            "count": decade_counts[d],
            "pct": round(decade_counts[d] / total * 100, 1) if total > 0 else 0.0,
        }
    return result


def classify_precision(mag: float) -> str:
    """Classify magnitude as one_decimal or two_decimal.

    Multiply by 100, round to integer, check if divisible by 10.
    """
    mag100 = round(mag * 100)
    if mag100 % 10 == 0:
        return "one_decimal"
    else:
        return "two_decimal"


def magnitude_precision(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, Any]]:
    """Classify magnitude precision for a catalog."""
    counts = {"one_decimal": 0, "two_decimal": 0}
    for r in rows:
        mag = float(r["usgs_mag"])
        counts[classify_precision(mag)] += 1

    total = len(rows)
    result = {}
    for cls in ["one_decimal", "two_decimal"]:
        result[cls] = {
            "count": counts[cls],
            "pct": round(counts[cls] / total * 100, 1),
        }
    logger.info("%s magnitude precision: %s", label, result)
    return result


def magnitude_bins(rows: list[dict[str, str]]) -> list[int]:
    """Compute event counts per 0.1-mag bin from 6.0 to 9.6."""
    # Bins: [6.0, 6.1), [6.1, 6.2), ..., [9.5, 9.6]
    bin_edges = [round(6.0 + i * 0.1, 1) for i in range(37)]  # 6.0 to 9.6
    counts = [0] * 36
    for r in rows:
        mag = float(r["usgs_mag"])
        idx = int(round((mag - 6.0) * 10))
        if idx < 0:
            idx = 0
        if idx >= 36:
            idx = 35
        counts[idx] += 1
    return counts


def main() -> None:
    """Run all analyses and write results JSON."""
    comcat_rows = load_csv(COMCAT_PATH)
    iscgem_rows = load_csv(ISCGEM_PATH)

    comcat_summary = population_summary(comcat_rows, "ComCat")
    iscgem_summary = population_summary(iscgem_rows, "ISC-GEM")

    prefix_breakdown = comcat_prefix_breakdown(comcat_rows)
    temporal = iscgem_prefix_temporal(comcat_rows)

    comcat_precision = magnitude_precision(comcat_rows, "ComCat")
    iscgem_precision = magnitude_precision(iscgem_rows, "ISC-GEM")

    comcat_bins = magnitude_bins(comcat_rows)
    iscgem_bins = magnitude_bins(iscgem_rows)

    bin_edges = [round(6.0 + i * 0.1, 1) for i in range(36)]

    results = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "catalogs": {
            "comcat": {
                "file": "data/global-sets/comcat_global_6-9_1949-2021.csv",
                **comcat_summary,
            },
            "iscgem": {
                "file": "data/global-sets/iscgem_global_events.csv",
                **iscgem_summary,
            },
        },
        "comcat_prefix_breakdown": prefix_breakdown,
        "comcat_iscgem_prefix_temporal": temporal,
        "magnitude_precision": {
            "comcat": comcat_precision,
            "iscgem": iscgem_precision,
        },
        "magnitude_bins": {
            "bin_width": 0.1,
            "bins": bin_edges,
            "comcat_counts": comcat_bins,
            "iscgem_counts": iscgem_bins,
        },
    }

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Results written to %s", OUTPUT_PATH)


if __name__ == "__main__":
    main()
