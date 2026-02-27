"""
Case A0b Analysis: Duplicate Detection and Cross-Catalog Event Accounting

Extends Case A0 by performing event-level matching between ComCat and ISC-GEM
catalogs to quantify within-ComCat duplication and cross-catalog divergence.

Three analyses:
  1. Direct ID cross-reference (iscgem-prefixed ComCat IDs vs ISC-GEM IDs)
  2. Within-ComCat duplicate detection (iscgem vs us_native subsets)
  3. Full cross-catalog proximity matching at two tolerance levels
"""

import json
import logging
import math
import os
from bisect import bisect_left, bisect_right
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ---------- paths ----------
BASE_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = BASE_DIR.parent
OUTPUT_DIR = BASE_DIR / "output"

COMCAT_PATH = PROJECT_ROOT / "data" / "global-sets" / "comcat_global_6-9_1949-2021.csv"
ISCGEM_PATH = PROJECT_ROOT / "data" / "global-sets" / "iscgem_global_events.csv"
PREREQ_PATH = OUTPUT_DIR / "case-a0-results.json"
RESULTS_PATH = OUTPUT_DIR / "case-a0b-results.json"


# ---------- helpers ----------

def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Compute great-circle distance in km using the Haversine formula."""
    R = 6371.0  # Earth radius in km
    rlat1, rlon1, rlat2, rlon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = rlat2 - rlat1
    dlon = rlon2 - rlon1
    a = math.sin(dlat / 2) ** 2 + math.cos(rlat1) * math.cos(rlat2) * math.sin(dlon / 2) ** 2
    return R * 2 * math.asin(math.sqrt(a))


def parse_event_at(s: str) -> float:
    """Parse ISO-8601 datetime string to Unix timestamp (float seconds)."""
    dt = datetime.fromisoformat(s.replace("Z", "+00:00"))
    return dt.timestamp()


def decade_label(year: int) -> str:
    """Return decade label like '1950s', '1960s', etc."""
    dec = (year // 10) * 10
    return f"{dec}s"


def mag_band_label(mag: float) -> str:
    """Return magnitude band: '6.0-6.4', '6.5-6.9', '7.0-7.4', '7.5+'."""
    if mag < 6.5:
        return "6.0-6.4"
    elif mag < 7.0:
        return "6.5-6.9"
    elif mag < 7.5:
        return "7.0-7.4"
    else:
        return "7.5+"


def load_catalog(path: Path) -> pd.DataFrame:
    """Load a catalog CSV and add derived columns."""
    df = pd.read_csv(path)
    df["timestamp"] = df["event_at"].apply(parse_event_at)
    return df


def get_id_prefix(usgs_id: str) -> str:
    """Classify ComCat ID prefix."""
    if usgs_id.startswith("iscgem"):
        return "iscgem"
    elif usgs_id.startswith("us"):
        return "us_native"
    else:
        return "other"


def find_proximity_matches(
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    time_sec: float,
    dist_km: float,
    mag_tol: float,
    one_to_one: bool = False,
) -> List[Tuple[int, int]]:
    """
    Find proximity matches between source and target DataFrames.

    For each source record, finds the best (closest in time) target record
    within the specified tolerances.

    Args:
        source_df: DataFrame sorted by timestamp with columns timestamp, latitude, longitude, usgs_mag
        target_df: DataFrame sorted by timestamp with columns timestamp, latitude, longitude, usgs_mag
        time_sec: Maximum absolute time difference in seconds
        dist_km: Maximum Haversine distance in km
        mag_tol: Maximum absolute magnitude difference
        one_to_one: If True, each target can match at most one source (greedy)

    Returns:
        List of (source_index, target_index) pairs
    """
    source_sorted = source_df.sort_values("timestamp").reset_index(drop=True)
    target_sorted = target_df.sort_values("timestamp").reset_index(drop=True)

    target_ts = target_sorted["timestamp"].values
    target_lats = target_sorted["latitude"].values
    target_lons = target_sorted["longitude"].values
    target_mags = target_sorted["usgs_mag"].values

    matches: List[Tuple[int, int, float]] = []  # (src_idx, tgt_idx, time_diff)

    for src_idx in range(len(source_sorted)):
        src_ts = source_sorted.iloc[src_idx]["timestamp"]
        src_lat = source_sorted.iloc[src_idx]["latitude"]
        src_lon = source_sorted.iloc[src_idx]["longitude"]
        src_mag = source_sorted.iloc[src_idx]["usgs_mag"]

        # Binary search for time window
        lo = bisect_left(target_ts, src_ts - time_sec)
        hi = bisect_right(target_ts, src_ts + time_sec)

        best_tgt_idx: Optional[int] = None
        best_time_diff = float("inf")

        for tgt_idx in range(lo, hi):
            # Magnitude filter
            if abs(src_mag - target_mags[tgt_idx]) > mag_tol:
                continue
            # Distance filter
            d = haversine_km(src_lat, src_lon, target_lats[tgt_idx], target_lons[tgt_idx])
            if d > dist_km:
                continue
            # Best temporal match
            td = abs(src_ts - target_ts[tgt_idx])
            if td < best_time_diff:
                best_time_diff = td
                best_tgt_idx = tgt_idx

        if best_tgt_idx is not None:
            matches.append((src_idx, best_tgt_idx, best_time_diff))

    if one_to_one:
        # Greedy: sort by time_diff, assign each target at most once
        matches.sort(key=lambda x: x[2])
        used_targets: set = set()
        used_sources: set = set()
        final: List[Tuple[int, int]] = []
        for src_idx, tgt_idx, _ in matches:
            if tgt_idx not in used_targets and src_idx not in used_sources:
                final.append((src_idx, tgt_idx))
                used_targets.add(tgt_idx)
                used_sources.add(src_idx)
        return final
    else:
        return [(s, t) for s, t, _ in matches]


# ---------- main ----------

def main() -> None:
    logger.info("Loading prerequisite: %s", PREREQ_PATH)
    assert PREREQ_PATH.exists(), f"Prerequisite not found: {PREREQ_PATH}"

    logger.info("Loading catalogs")
    comcat = load_catalog(COMCAT_PATH)
    iscgem = load_catalog(ISCGEM_PATH)

    comcat["id_prefix"] = comcat["usgs_id"].apply(get_id_prefix)

    comcat_iscgem = comcat[comcat["id_prefix"] == "iscgem"].copy()
    comcat_us_native = comcat[comcat["id_prefix"] == "us_native"].copy()

    logger.info(
        "ComCat: %d total, %d iscgem-prefix, %d us_native",
        len(comcat), len(comcat_iscgem), len(comcat_us_native),
    )

    # ======== Analysis 1: Direct ID Cross-Reference ========
    logger.info("Analysis 1: Direct ID cross-reference")

    comcat_iscgem_ids = set(comcat_iscgem["usgs_id"].values)
    iscgem_ids = set(iscgem["usgs_id"].values)

    id_matched = comcat_iscgem_ids & iscgem_ids
    id_unmatched = comcat_iscgem_ids - iscgem_ids
    iscgem_not_in_comcat = iscgem_ids - set(comcat["usgs_id"].values)

    id_cross_ref = {
        "comcat_iscgem_prefix_count": len(comcat_iscgem),
        "id_matched_in_iscgem": {
            "count": len(id_matched),
            "pct_of_comcat_iscgem": round(100 * len(id_matched) / len(comcat_iscgem), 1),
        },
        "id_unmatched_in_iscgem": {
            "count": len(id_unmatched),
            "pct_of_comcat_iscgem": round(100 * len(id_unmatched) / len(comcat_iscgem), 1),
        },
        "iscgem_ids_not_in_comcat": {
            "count": len(iscgem_not_in_comcat),
            "pct_of_iscgem": round(100 * len(iscgem_not_in_comcat) / len(iscgem), 1),
        },
    }

    logger.info("ID matched: %d, unmatched: %d, ISC-GEM not in ComCat: %d",
                len(id_matched), len(id_unmatched), len(iscgem_not_in_comcat))

    # ======== Analysis 2: Within-ComCat Duplicate Detection ========
    logger.info("Analysis 2: Within-ComCat duplicate detection")

    dup_matches = find_proximity_matches(
        comcat_iscgem, comcat_us_native,
        time_sec=60, dist_km=50, mag_tol=0.3,
        one_to_one=False,
    )

    dup_count = len(dup_matches)
    logger.info("Duplicate candidates found: %d", dup_count)

    # Temporal distribution of duplicates by decade
    comcat_iscgem_sorted = comcat_iscgem.sort_values("timestamp").reset_index(drop=True)
    dup_decades: Dict[str, int] = {"1950s": 0, "1960s": 0, "1970s": 0, "1980s_plus": 0}
    for src_idx, _ in dup_matches:
        year = int(comcat_iscgem_sorted.iloc[src_idx]["solaration_year"])
        if year < 1960:
            dup_decades["1950s"] += 1
        elif year < 1970:
            dup_decades["1960s"] += 1
        elif year < 1980:
            dup_decades["1970s"] += 1
        else:
            dup_decades["1980s_plus"] += 1

    within_comcat_dups = {
        "tolerances": {"time_seconds": 60, "distance_km": 50, "magnitude": 0.3},
        "duplicate_candidate_count": dup_count,
        "pct_of_iscgem_prefix_subset": round(100 * dup_count / len(comcat_iscgem), 1),
        "net_effective_comcat_count": 9802 - dup_count,
        "temporal_distribution": dup_decades,
    }

    # ======== Analysis 3: Full Cross-Catalog Proximity Matching ========
    logger.info("Analysis 3: Full cross-catalog proximity matching")

    def run_cross_catalog(
        time_sec: float, dist_km: float, mag_tol: float, label: str
    ) -> Dict[str, Any]:
        logger.info("Running %s tolerance: %ds, %dkm, %.1f mag",
                     label, time_sec, dist_km, mag_tol)
        matches = find_proximity_matches(
            comcat, iscgem,
            time_sec=time_sec, dist_km=dist_km, mag_tol=mag_tol,
            one_to_one=True,
        )
        matched_count = len(matches)
        comcat_only_count = len(comcat) - matched_count
        iscgem_only_count = len(iscgem) - matched_count

        logger.info("%s: matched=%d, comcat_only=%d, iscgem_only=%d",
                     label, matched_count, comcat_only_count, iscgem_only_count)

        result: Dict[str, Any] = {
            "tolerances": {
                "time_seconds": time_sec,
                "distance_km": dist_km,
                "magnitude": mag_tol,
            },
            "matched": matched_count,
            "comcat_only": comcat_only_count,
            "iscgem_only": iscgem_only_count,
        }

        if label == "primary":
            # Get the matched source indices to identify unmatched
            comcat_sorted = comcat.sort_values("timestamp").reset_index(drop=True)
            iscgem_sorted = iscgem.sort_values("timestamp").reset_index(drop=True)

            matched_comcat_idx = set(s for s, _ in matches)
            matched_iscgem_idx = set(t for _, t in matches)

            comcat_unmatched = comcat_sorted[~comcat_sorted.index.isin(matched_comcat_idx)]
            iscgem_unmatched = iscgem_sorted[~iscgem_sorted.index.isin(matched_iscgem_idx)]

            # Temporal distribution by decade
            comcat_only_temporal: Dict[str, int] = {}
            for _, row in comcat_unmatched.iterrows():
                dec = decade_label(int(row["solaration_year"]))
                comcat_only_temporal[dec] = comcat_only_temporal.get(dec, 0) + 1

            iscgem_only_temporal: Dict[str, int] = {}
            for _, row in iscgem_unmatched.iterrows():
                dec = decade_label(int(row["solaration_year"]))
                iscgem_only_temporal[dec] = iscgem_only_temporal.get(dec, 0) + 1

            # Magnitude distribution by band
            comcat_only_mag: Dict[str, int] = {}
            for _, row in comcat_unmatched.iterrows():
                band = mag_band_label(row["usgs_mag"])
                comcat_only_mag[band] = comcat_only_mag.get(band, 0) + 1

            iscgem_only_mag: Dict[str, int] = {}
            for _, row in iscgem_unmatched.iterrows():
                band = mag_band_label(row["usgs_mag"])
                iscgem_only_mag[band] = iscgem_only_mag.get(band, 0) + 1

            result["comcat_only_temporal"] = dict(sorted(comcat_only_temporal.items()))
            result["iscgem_only_temporal"] = dict(sorted(iscgem_only_temporal.items()))
            result["comcat_only_by_mag_band"] = comcat_only_mag
            result["iscgem_only_by_mag_band"] = iscgem_only_mag

        return result

    primary_result = run_cross_catalog(60, 50, 0.3, "primary")
    strict_result = run_cross_catalog(30, 25, 0.2, "strict")

    # ======== Write results ========
    results = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "prerequisites": {
            "case_a0_results": "topic-adhoc/output/case-a0-results.json",
        },
        "id_cross_reference": id_cross_ref,
        "within_comcat_duplicates": within_comcat_dups,
        "cross_catalog_matching": {
            "primary": primary_result,
            "strict": strict_result,
        },
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)

    logger.info("Results written to %s", RESULTS_PATH)


if __name__ == "__main__":
    main()
