"""Sub-analysis B: Post-declustering interval structure comparison.

Identifies elevated phase bins in each mainshock catalog and compares the
recovered intervals against the A1b baseline intervals established in
Adhoc Case A1b.
"""

import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SOLAR_YEAR_SECS = 31_557_600.0  # Julian constant: 365.25 * 86400

# A1b baseline intervals — from Adhoc Case A1b
A1B_BASELINE_INTERVALS = [
    {"id": 1, "phase_start": 0.1875, "phase_end": 0.25, "calendar": "~Mar 10 - Apr 1"},
    {"id": 2, "phase_start": 0.625, "phase_end": 0.656, "calendar": "~Aug 16 - Aug 28"},
    {"id": 3, "phase_start": 0.875, "phase_end": 0.917, "calendar": "~Nov 16 - Dec 1"},
]


def compute_phase(solar_secs: pd.Series, year_secs: float = SOLAR_YEAR_SECS) -> pd.Series:
    """Compute solar phase in [0, 1).

    Args:
        solar_secs: Series of seconds elapsed since solar year start.
        year_secs: Length of solar year in seconds.

    Returns:
        Phase values in [0, 1).
    """
    return (solar_secs / year_secs) % 1.0


def identify_elevated_bins(obs: np.ndarray, expected_per_bin: float) -> list[int]:
    """Return list of bin indices exceeding E + sqrt(E) threshold.

    Args:
        obs: Observed bin counts array.
        expected_per_bin: Uniform expected count per bin.

    Returns:
        Sorted list of elevated bin indices.
    """
    threshold = expected_per_bin + np.sqrt(expected_per_bin)
    return [int(i) for i in range(len(obs)) if obs[i] > threshold]


def merge_adjacent_bins(elevated: list[int], k: int) -> list[tuple[int, int]]:
    """Merge consecutive bin indices into contiguous ranges.

    Args:
        elevated: Sorted list of elevated bin indices.
        k: Total number of bins (for phase conversion).

    Returns:
        List of (start_bin, end_bin) tuples (inclusive).
    """
    if not elevated:
        return []
    merged = []
    start = elevated[0]
    prev = elevated[0]
    for idx in elevated[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            merged.append((start, prev))
            start = idx
            prev = idx
    merged.append((start, prev))
    return merged


def bins_to_phase_interval(start_bin: int, end_bin: int, k: int) -> tuple[float, float]:
    """Convert bin range to phase interval.

    Args:
        start_bin: First bin index (inclusive).
        end_bin: Last bin index (inclusive).
        k: Total bin count.

    Returns:
        (phase_start, phase_end) as floats.
    """
    return float(start_bin / k), float((end_bin + 1) / k)


def overlap_fraction(a_start: float, a_end: float, b_start: float, b_end: float) -> float:
    """Compute what fraction of interval A overlaps with interval B.

    Args:
        a_start, a_end: Bounds of interval A.
        b_start, b_end: Bounds of interval B.

    Returns:
        Overlap length / width of A.
    """
    a_width = a_end - a_start
    if a_width <= 0:
        return 0.0
    overlap = max(0.0, min(a_end, b_end) - max(a_start, b_start))
    return overlap / a_width


def classify_interval(phase_start: float, phase_end: float) -> str:
    """Classify a recovered interval against A1b baseline intervals.

    Args:
        phase_start: Start of recovered interval.
        phase_end: End of recovered interval.

    Returns:
        Classification string.
    """
    for baseline in A1B_BASELINE_INTERVALS:
        frac = overlap_fraction(
            phase_start, phase_end,
            baseline["phase_start"], baseline["phase_end"],
        )
        if frac > 0.5:
            return f"matches interval {baseline['id']}"
    return "new interval"


def compute_interval_coherence(
    phase: np.ndarray,
    phase_start: float,
    phase_end: float,
) -> float:
    """Compute mean resultant length R for phases within a given interval.

    Args:
        phase: All event phases in [0, 1).
        phase_start: Start of interval.
        phase_end: End of interval.

    Returns:
        Rayleigh R for the subset, or 0.0 if no events fall in interval.
    """
    mask = (phase >= phase_start) & (phase < phase_end)
    sub = phase[mask]
    if len(sub) == 0:
        return 0.0
    angles = 2.0 * np.pi * sub
    return float(np.abs(np.mean(np.exp(1j * angles))))


def run_sub_b_single(
    catalog: pd.DataFrame,
    catalog_name: str,
    k: int,
) -> dict[str, Any]:
    """Run Sub-analysis B for a single catalog at a single bin count.

    Args:
        catalog: DataFrame with solar_secs column.
        catalog_name: Human-readable name for logging.
        k: Number of phase bins.

    Returns:
        Dict with recovered_intervals and baseline_survival.
    """
    n = len(catalog)
    phase = compute_phase(catalog["solar_secs"]).values

    bin_idx = np.floor(phase * k).astype(int)
    bin_idx = np.clip(bin_idx, 0, k - 1)
    obs = np.bincount(bin_idx, minlength=k).astype(float)
    expected_per_bin = n / k

    elevated = identify_elevated_bins(obs, expected_per_bin)
    merged_bins = merge_adjacent_bins(elevated, k)

    recovered_intervals = []
    for (sb, eb) in merged_bins:
        ps, pe = bins_to_phase_interval(sb, eb, k)
        classification = classify_interval(ps, pe)
        r_coherence = compute_interval_coherence(phase, ps, pe)
        n_in = int(np.sum((phase >= ps) & (phase < pe)))
        recovered_intervals.append({
            "phase_start": ps,
            "phase_end": pe,
            "classification": classification,
            "n_events": n_in,
            "R_coherence": float(r_coherence),
        })

    # Baseline survival: does each baseline interval appear in recovered intervals?
    baseline_survival: dict[str, Any] = {}
    for baseline in A1B_BASELINE_INTERVALS:
        matched = [
            ri for ri in recovered_intervals
            if ri["classification"] == f"matches interval {baseline['id']}"
        ]
        baseline_survival[f"interval_{baseline['id']}"] = {
            "survives": len(matched) > 0,
            "matching_recovered": matched,
        }

    logger.info(
        "%s k=%d: %d elevated bins, %d merged intervals",
        catalog_name, k, len(elevated), len(merged_bins),
    )
    return {
        "recovered_intervals": recovered_intervals,
        "baseline_survival": baseline_survival,
    }


def run_sub_b(
    mainshock_catalogs: dict[str, pd.DataFrame],
) -> dict[str, Any]:
    """Run Sub-analysis B for all mainshock catalogs and bin counts.

    Args:
        mainshock_catalogs: Dict with keys "gk_mainshocks", "reas_mainshocks", "a1b_mainshocks".

    Returns:
        Sub-analysis B result dict.
    """
    bin_counts = [16, 24, 32]

    results: dict[str, Any] = {
        "a1b_baseline_intervals": A1B_BASELINE_INTERVALS,
    }

    catalog_map = {
        "gk_mainshocks": "G-K Mainshocks",
        "reas_mainshocks": "Reasenberg Mainshocks",
        "a1b_mainshocks": "A1b Mainshocks",
    }

    for key, name in catalog_map.items():
        results[key] = {}
        df = mainshock_catalogs[key]
        for k in bin_counts:
            k_key = f"k{k}"
            results[key][k_key] = run_sub_b_single(df, name, k)

    # Interval survival summary across methods
    survival_summary: dict[str, Any] = {}
    for bsl in A1B_BASELINE_INTERVALS:
        interval_id = f"interval_{bsl['id']}"
        entry: dict[str, str] = {}
        for method_key, method_label in [
            ("gk_mainshocks", "gk"),
            ("reas_mainshocks", "reas"),
            ("a1b_mainshocks", "a1b"),
        ]:
            # Use k=24 as the canonical bin count for the summary
            survival_at_k24 = results[method_key]["k24"]["baseline_survival"][interval_id]["survives"]
            entry[method_label] = "survives" if survival_at_k24 else "absent"
        survival_summary[interval_id] = entry
    results["interval_survival_summary"] = survival_summary

    return results
