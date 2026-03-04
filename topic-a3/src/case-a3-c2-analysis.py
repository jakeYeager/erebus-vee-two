"""
Case A3.C2: Targeted Major Sequence Phased Declustering Test

Sequential removal of M>=8.5 event sequences from the ISC-GEM catalog to test
whether the global solar-phase signal is diffuse across events or concentrated
in a small number of major aftershock sequences.

Run order: largest magnitude first (descending by usgs_mag).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import chisquare

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent

RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GK_MS_PATH = (
    BASE_DIR.parent
    / "data"
    / "iscgem"
    / "declustering-algorithm"
    / "mainshocks_gk-seq_global.csv"
)
GK_AS_PATH = (
    BASE_DIR.parent
    / "data"
    / "iscgem"
    / "declustering-algorithm"
    / "aftershocks_gk-seq_global.csv"
)
REAS_MS_PATH = (
    BASE_DIR.parent
    / "data"
    / "iscgem"
    / "declustering-algorithm"
    / "mainshocks_reas-seq_global.csv"
)
REAS_AS_PATH = (
    BASE_DIR.parent
    / "data"
    / "iscgem"
    / "declustering-algorithm"
    / "aftershocks_reas-seq_global.csv"
)

OUTPUT_PATH = BASE_DIR / "output" / "case-a3-c2-results.json"

K_BINS: int = 24
MAG_THRESHOLD: float = 8.5
JULIAN_YEAR_SECS: float = 31_557_600.0

INTERVAL_BINS: dict[str, list[int]] = {
    "interval_1": [4, 5],
    "interval_2": [15],
    "interval_3": [21],
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Section 1: Data loading
# ---------------------------------------------------------------------------


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load all five ISC-GEM catalog files and compute solar phase.

    Returns:
        Tuple of (raw, gk_ms, gk_as, reas_ms, reas_as) DataFrames, each with
        a ``phase`` column added.

    Raises:
        AssertionError: if row counts do not match expected values.
    """
    log.info("Loading raw catalog from %s", RAW_PATH)
    raw = pd.read_csv(RAW_PATH)
    raw["event_at"] = pd.to_datetime(raw["event_at"], utc=True)
    log.info("raw: %d rows", len(raw))
    assert len(raw) == 9210, f"raw row count mismatch: {len(raw)}"

    log.info("Loading G-K mainshocks from %s", GK_MS_PATH)
    gk_ms = pd.read_csv(GK_MS_PATH)
    gk_ms["event_at"] = pd.to_datetime(gk_ms["event_at"], utc=True)
    log.info("gk_ms: %d rows", len(gk_ms))
    assert len(gk_ms) == 5883, f"gk_ms row count mismatch: {len(gk_ms)}"

    log.info("Loading G-K aftershocks from %s", GK_AS_PATH)
    gk_as = pd.read_csv(GK_AS_PATH)
    gk_as["event_at"] = pd.to_datetime(gk_as["event_at"], utc=True)
    log.info("gk_as: %d rows", len(gk_as))
    assert len(gk_as) == 3327, f"gk_as row count mismatch: {len(gk_as)}"

    log.info("Loading Reasenberg mainshocks from %s", REAS_MS_PATH)
    reas_ms = pd.read_csv(REAS_MS_PATH)
    reas_ms["event_at"] = pd.to_datetime(reas_ms["event_at"], utc=True)
    log.info("reas_ms: %d rows", len(reas_ms))
    assert len(reas_ms) == 8265, f"reas_ms row count mismatch: {len(reas_ms)}"

    log.info("Loading Reasenberg aftershocks from %s", REAS_AS_PATH)
    reas_as = pd.read_csv(REAS_AS_PATH)
    reas_as["event_at"] = pd.to_datetime(reas_as["event_at"], utc=True)
    log.info("reas_as: %d rows", len(reas_as))
    assert len(reas_as) == 945, f"reas_as row count mismatch: {len(reas_as)}"

    for df in (raw, gk_ms, gk_as, reas_ms, reas_as):
        df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0

    return raw, gk_ms, gk_as, reas_ms, reas_as


# ---------------------------------------------------------------------------
# Section 2: Major event identification and classification
# ---------------------------------------------------------------------------


def identify_major_events(raw_df: pd.DataFrame, mag_threshold: float) -> pd.DataFrame:
    """Identify events at or above the magnitude threshold, sorted largest first.

    Args:
        raw_df: Full raw catalog DataFrame.
        mag_threshold: Minimum magnitude for inclusion (e.g. 8.5).

    Returns:
        DataFrame containing major events sorted by usgs_mag descending.
        Columns: usgs_id, usgs_mag, event_at, latitude, longitude, depth.
    """
    major = (
        raw_df[raw_df["usgs_mag"] >= mag_threshold]
        .sort_values("usgs_mag", ascending=False)
        .reset_index(drop=True)
    )
    major = major[["usgs_id", "usgs_mag", "event_at", "latitude", "longitude", "depth"]].copy()
    log.info(
        "Identified %d major events at M>=%.1f (expected ~10–15)",
        len(major),
        mag_threshold,
    )
    return major


def classify_event_in_catalog(
    usgs_id: str,
    mainshock_df: pd.DataFrame,
    aftershock_df: pd.DataFrame,
) -> dict:
    """Classify a single event as mainshock, aftershock, or absent in a declustered catalog.

    Args:
        usgs_id: USGS identifier string for the event.
        mainshock_df: Mainshock catalog DataFrame.
        aftershock_df: Aftershock catalog DataFrame.

    Returns:
        Dict with keys ``is_mainshock``, ``is_aftershock``, ``parent_id``.
    """
    is_mainshock = usgs_id in mainshock_df["usgs_id"].values
    is_aftershock = usgs_id in aftershock_df["usgs_id"].values

    parent_id: Optional[str] = None
    if is_aftershock:
        row = aftershock_df[aftershock_df["usgs_id"] == usgs_id].iloc[0]
        parent_id = str(row["parent_id"])

    if not is_mainshock and not is_aftershock:
        log.warning(
            "Event %s is neither mainshock nor aftershock in this catalog — edge case",
            usgs_id,
        )

    return {
        "is_mainshock": bool(is_mainshock),
        "is_aftershock": bool(is_aftershock),
        "parent_id": parent_id,
    }


# ---------------------------------------------------------------------------
# Section 3: Sequence metrics computation
# ---------------------------------------------------------------------------


def compute_sequence_metrics(
    major_event: pd.Series,
    mainshock_df: pd.DataFrame,
    aftershock_df: pd.DataFrame,
    label: str,
) -> dict:
    """Compute per-event sequence attribution metrics from one declustered catalog.

    Args:
        major_event: A row from the major events DataFrame.
        mainshock_df: Mainshock catalog for this declustering.
        aftershock_df: Aftershock catalog for this declustering.
        label: Identifier string for logging (e.g. "gk" or "reasenberg").

    Returns:
        Dict with sequence metric fields as specified in the case spec.
    """
    usgs_id = major_event["usgs_id"]
    classification_result = classify_event_in_catalog(usgs_id, mainshock_df, aftershock_df)
    is_mainshock = classification_result["is_mainshock"]
    is_aftershock = classification_result["is_aftershock"]

    event_at_str = (
        major_event["event_at"].isoformat()
        if hasattr(major_event["event_at"], "isoformat")
        else str(major_event["event_at"])
    )

    if is_mainshock:
        ms_row = mainshock_df[mainshock_df["usgs_id"] == usgs_id].iloc[0]
        foreshock_count = int(ms_row["foreshock_count"])
        aftershock_count = int(ms_row["aftershock_count"])
        total_sequence = foreshock_count + aftershock_count + 1
        window_secs = float(ms_row["window_secs"])
        window_days = window_secs / 86400.0
        window_km = float(ms_row["window_km"])

        as_rows = aftershock_df[aftershock_df["parent_id"] == usgs_id]
        early_count = int((as_rows["delta_t_sec"] <= window_secs / 2).sum())
        late_count = int((as_rows["delta_t_sec"] > window_secs / 2).sum())
        early_pct = float(early_count / max(aftershock_count, 1))
        classification = "mainshock"

    elif is_aftershock:
        log.info(
            "Event %s classified as aftershock in %s — no sequence attributed", usgs_id, label
        )
        foreshock_count = 0
        aftershock_count = 0
        total_sequence = 1
        window_secs = 0.0
        window_days = 0.0
        window_km = 0.0
        early_count = 0
        late_count = 0
        early_pct = None
        classification = "aftershock"

    else:
        log.info(
            "Event %s absent from %s — may have been below M>=6.0 threshold", usgs_id, label
        )
        foreshock_count = 0
        aftershock_count = 0
        total_sequence = 1
        window_secs = 0.0
        window_days = 0.0
        window_km = 0.0
        early_count = 0
        late_count = 0
        early_pct = None
        classification = "absent"

    return {
        "usgs_id": usgs_id,
        "usgs_mag": float(major_event["usgs_mag"]),
        "event_at": event_at_str,
        "latitude": float(major_event["latitude"]),
        "longitude": float(major_event["longitude"]),
        "declustering": label,
        "classification": classification,
        "foreshock_count": foreshock_count,
        "aftershock_count": aftershock_count,
        "total_sequence": total_sequence,
        "window_secs": window_secs,
        "window_days": window_days,
        "window_km": window_km,
        "early_count": early_count,
        "late_count": late_count,
        "early_pct": early_pct,
    }


# ---------------------------------------------------------------------------
# Section 4: Chi-square and interval statistics
# ---------------------------------------------------------------------------


def compute_chi2_stats(phases: np.ndarray, k: int = 24) -> dict:
    """Compute chi-square uniformity statistics and A1b interval z-scores.

    Args:
        phases: Array of solar phase values in [0, 1).
        k: Number of bins for chi-square test.

    Returns:
        Dict with n, chi2_k24, p_chi2_k24, cramers_v, bin_counts,
        interval_1_z, interval_2_z, interval_3_z.
    """
    n = len(phases)
    if n < k:
        log.warning("n=%d < k=%d bins — degenerate case, returning null chi2", n, k)
        return {
            "n": n,
            "chi2_k24": None,
            "p_chi2_k24": None,
            "cramers_v": None,
            "bin_counts": [],
            "interval_1_z": 0.0,
            "interval_2_z": 0.0,
            "interval_3_z": 0.0,
        }

    bin_indices = np.floor(phases * k).astype(int) % k
    observed = np.bincount(bin_indices, minlength=k)
    expected = np.full(k, n / k)

    chi2_stat, p_chi2 = chisquare(observed, expected)
    cramers_v = float(np.sqrt(chi2_stat / (n * (k - 1)))) if n > 0 else 0.0

    interval_z: dict[str, float] = {}
    for interval_name, bins in INTERVAL_BINS.items():
        interval_count = float(observed[bins].sum())
        interval_expected = float(len(bins) * (n / k))
        if interval_expected > 0:
            z = (interval_count - interval_expected) / np.sqrt(interval_expected)
        else:
            z = 0.0
        interval_z[interval_name] = float(z)

    return {
        "n": n,
        "chi2_k24": float(chi2_stat),
        "p_chi2_k24": float(p_chi2),
        "cramers_v": cramers_v,
        "bin_counts": observed.tolist(),
        "interval_1_z": interval_z["interval_1"],
        "interval_2_z": interval_z["interval_2"],
        "interval_3_z": interval_z["interval_3"],
    }


# ---------------------------------------------------------------------------
# Section 5: Phased removal analysis
# ---------------------------------------------------------------------------


def run_phased_removal(
    base_df: pd.DataFrame,
    major_events: pd.DataFrame,
    aftershock_df: Optional[pd.DataFrame],
    run_key: str,
) -> list[dict]:
    """Run sequential phased removal of major events and their attributed sequences.

    For raw_gk and raw_reas run keys, the attributed aftershocks (via parent_id)
    are also removed at each step. For mainshock_gk and mainshock_reas, only the
    major event row itself is removed (mainshock-only catalogs).

    Args:
        base_df: The starting catalog DataFrame (with ``phase`` column).
        major_events: DataFrame of major events sorted by magnitude descending.
        aftershock_df: Aftershock catalog for `parent_id` lookup; None for
            mainshock-only runs.
        run_key: One of 'raw_gk', 'raw_reas', 'mainshock_gk', 'mainshock_reas'.

    Returns:
        List of step dicts with removal metadata and chi-square statistics.
    """
    use_aftershocks = run_key in ("raw_gk", "raw_reas") and aftershock_df is not None

    steps: list[dict] = []
    accumulated_ids: set[str] = set()

    # Step 0: baseline — no removals
    baseline_stats = compute_chi2_stats(base_df["phase"].values)
    steps.append(
        {
            "step": 0,
            "event_removed": None,
            "n_removed_cumulative": 0,
            "ids_removed_cumulative_count": 0,
            "n_remaining": int(len(base_df)),
            "stats": {
                "chi2_k24": baseline_stats["chi2_k24"],
                "p_chi2_k24": baseline_stats["p_chi2_k24"],
                "cramers_v": baseline_stats["cramers_v"],
                "interval_1_z": baseline_stats["interval_1_z"],
                "interval_2_z": baseline_stats["interval_2_z"],
                "interval_3_z": baseline_stats["interval_3_z"],
            },
        }
    )
    log.info(
        "[%s] Step 0 baseline: n=%d, chi2=%.4f, p=%.6f, V=%.6f",
        run_key,
        baseline_stats["n"],
        baseline_stats["chi2_k24"] or 0.0,
        baseline_stats["p_chi2_k24"] or 1.0,
        baseline_stats["cramers_v"] or 0.0,
    )

    for i, (_, major_event) in enumerate(major_events.iterrows(), start=1):
        ev_id = major_event["usgs_id"]

        # Always remove the major event itself
        ids_this_step: set[str] = {ev_id}

        # For raw runs, also remove attributed aftershocks
        n_as_removed = 0
        if use_aftershocks:
            attributed_as_ids = aftershock_df[
                aftershock_df["parent_id"] == ev_id
            ]["usgs_id"].tolist()
            ids_this_step.update(attributed_as_ids)
            n_as_removed = len(attributed_as_ids)

        n_removed_this_step = len(ids_this_step - accumulated_ids)
        accumulated_ids.update(ids_this_step)

        remaining_df = base_df[~base_df["usgs_id"].isin(accumulated_ids)]
        stats = compute_chi2_stats(remaining_df["phase"].values)

        event_at_val = major_event["event_at"]
        event_at_str = (
            event_at_val.isoformat()
            if hasattr(event_at_val, "isoformat")
            else str(event_at_val)
        )

        step_record: dict = {
            "step": i,
            "event_removed": {
                "usgs_id": ev_id,
                "usgs_mag": float(major_event["usgs_mag"]),
                "event_at": event_at_str,
                "latitude": float(major_event["latitude"]),
                "longitude": float(major_event["longitude"]),
                "n_removed_this_step": n_removed_this_step,
            },
            "n_removed_cumulative": len(accumulated_ids),
            "ids_removed_cumulative_count": len(accumulated_ids),
            "n_remaining": int(len(remaining_df)),
            "stats": {
                "chi2_k24": stats["chi2_k24"],
                "p_chi2_k24": stats["p_chi2_k24"],
                "cramers_v": stats["cramers_v"],
                "interval_1_z": stats["interval_1_z"],
                "interval_2_z": stats["interval_2_z"],
                "interval_3_z": stats["interval_3_z"],
            },
        }
        steps.append(step_record)

        log.info(
            "[%s] Step %d (M%.2f %s): removed %d (this step), %d (cumulative) → n=%d, "
            "chi2=%.4f, p=%.6f, V=%.6f",
            run_key,
            i,
            major_event["usgs_mag"],
            event_at_str[:10],
            n_removed_this_step,
            len(accumulated_ids),
            len(remaining_df),
            stats["chi2_k24"] or 0.0,
            stats["p_chi2_k24"] or 1.0,
            stats["cramers_v"] or 0.0,
        )

    return steps


# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------


def main() -> None:
    """Execute full A3.C2 analysis and write results JSON."""
    log.info("=== Case A3.C2: Targeted Major Sequence Phased Declustering Test ===")

    # 1. Load data
    raw, gk_ms, gk_as, reas_ms, reas_as = load_data()

    # 2. Identify major events
    major_events = identify_major_events(raw, MAG_THRESHOLD)
    log.info("Major events (removal order, largest first):")
    for idx, row in major_events.iterrows():
        log.info("  %d. M%.2f  %s  id=%s", idx + 1, row["usgs_mag"], str(row["event_at"])[:10], row["usgs_id"])

    # Build major_events list for JSON
    major_events_list = []
    for order_idx, (_, row) in enumerate(major_events.iterrows(), start=1):
        event_at_str = (
            row["event_at"].isoformat()
            if hasattr(row["event_at"], "isoformat")
            else str(row["event_at"])
        )
        major_events_list.append(
            {
                "usgs_id": row["usgs_id"],
                "usgs_mag": float(row["usgs_mag"]),
                "event_at": event_at_str,
                "latitude": float(row["latitude"]),
                "longitude": float(row["longitude"]),
                "removal_order": order_idx,
            }
        )

    # 3. Sequence metrics
    log.info("Computing sequence metrics for G-K declustering...")
    gk_seq_metrics = [
        compute_sequence_metrics(row, gk_ms, gk_as, "gk")
        for _, row in major_events.iterrows()
    ]

    log.info("Computing sequence metrics for Reasenberg declustering...")
    reas_seq_metrics = [
        compute_sequence_metrics(row, reas_ms, reas_as, "reasenberg")
        for _, row in major_events.iterrows()
    ]

    # 4. Phased removal — four runs
    log.info("Running phased removal: raw_gk ...")
    steps_raw_gk = run_phased_removal(raw, major_events, gk_as, "raw_gk")

    log.info("Running phased removal: raw_reas ...")
    steps_raw_reas = run_phased_removal(raw, major_events, reas_as, "raw_reas")

    log.info("Running phased removal: mainshock_gk ...")
    steps_ms_gk = run_phased_removal(gk_ms, major_events, None, "mainshock_gk")

    log.info("Running phased removal: mainshock_reas ...")
    steps_ms_reas = run_phased_removal(reas_ms, major_events, None, "mainshock_reas")

    # 5. Assemble results JSON
    results = {
        "case": "A3.C2",
        "title": "Targeted Major Sequence Phased Declustering Test",
        "parameters": {
            "n_raw": 9210,
            "n_gk_mainshocks": 5883,
            "n_gk_aftershocks": 3327,
            "n_reas_mainshocks": 8265,
            "n_reas_aftershocks": 945,
            "mag_threshold": MAG_THRESHOLD,
            "k_bins": K_BINS,
            "julian_year_secs": JULIAN_YEAR_SECS,
        },
        "major_events": major_events_list,
        "sequence_metrics": {
            "gk": gk_seq_metrics,
            "reasenberg": reas_seq_metrics,
        },
        "runs": {
            "raw_gk": {"steps": steps_raw_gk},
            "raw_reas": {"steps": steps_raw_reas},
            "mainshock_gk": {"steps": steps_ms_gk},
            "mainshock_reas": {"steps": steps_ms_reas},
        },
    }

    # Write output
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(results, f, indent=2, default=str)

    log.info("Results written to %s", OUTPUT_PATH)


if __name__ == "__main__":
    main()
