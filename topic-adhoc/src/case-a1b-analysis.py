"""
Case A1b: Elevated Bin Event Characterization and Declustering Implications

Drills into the events driving the solar_secs elevated bins identified in
Case A1 (k=16, 24, 32) to determine whether they show temporal/spatial
clustering (sequence residuals) or are dispersed (genuine phase-preferring
events).
"""

import json
import logging
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Path conventions
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent          # topic-adhoc/
PROJECT_ROOT = BASE_DIR.parent                             # erebus-vee-two/
DATA_PATH = PROJECT_ROOT / "data" / "global-sets" / "iscgem_global_events.csv"
BOUNDARY_PATH = PROJECT_ROOT / "lib" / "pb2002_boundaries.dig"
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a1b-results.json"

# Phase normalization constants (identical to A1)
LUNAR_CYCLE_SECS = 29.53059 * 86400  # 2_551_442.976
DAY_SECS = 86400

# Analysis parameters
BIN_COUNTS = [16, 24, 32]
TOP_N = 3           # top bins to select per k
N_SAMPLES = 1000    # null distribution samples
RANDOM_SEED = 42

# Boundary proximity thresholds (km)
NEAR_BOUNDARY_KM = 100.0
TRANSITIONAL_KM = 300.0

# G-K reference windows
GK_REFERENCE = {
    "m6": {"spatial_km": 49, "temporal_days": 295},
    "m7": {"spatial_km": 156, "temporal_days": 790},
    "m8": {"spatial_km": 493, "temporal_days": 2117},
}


# ---------------------------------------------------------------------------
# Utility: phase normalization (reuse from A1)
# ---------------------------------------------------------------------------

def is_leap_year(year: int) -> bool:
    """Return True if year is a Gregorian leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def year_length_seconds(year: int) -> int:
    """Return the calendar-year length in seconds for year."""
    return 366 * DAY_SECS if is_leap_year(year) else 365 * DAY_SECS


def compute_solar_phase(df: pd.DataFrame) -> pd.Series:
    """Compute solar phase in [0, 1) using actual calendar year length.

    Clamps any marginal overshoot (within 1% of cycle) to 0.999999.
    """
    year_lengths = df["solaration_year"].apply(year_length_seconds).astype(float)
    phase = df["solar_secs"] / year_lengths
    # Clamp boundary values
    over_mask = phase >= 1.0
    if over_mask.any():
        logger.warning("Clamping %d solar_phase values >= 1.0 to 0.999999", over_mask.sum())
        phase = phase.copy()
        phase[over_mask] = 0.999999
    return phase


# ---------------------------------------------------------------------------
# Utility: Haversine distance
# ---------------------------------------------------------------------------

EARTH_RADIUS_KM = 6371.0


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Compute great-circle distance in km between two points (degrees)."""
    r = EARTH_RADIUS_KM
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlon / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def haversine_array(lat: np.ndarray, lon: np.ndarray,
                    lat_ref: np.ndarray, lon_ref: np.ndarray) -> np.ndarray:
    """Vectorized Haversine: distance from each (lat[i], lon[i]) to each row of ref array.

    Returns a (n, m) matrix where n = len(lat), m = len(lat_ref).
    """
    r = EARTH_RADIUS_KM
    lat_r = np.radians(lat)[:, np.newaxis]
    lon_r = np.radians(lon)[:, np.newaxis]
    lat_ref_r = np.radians(lat_ref)[np.newaxis, :]
    lon_ref_r = np.radians(lon_ref)[np.newaxis, :]

    dlat = lat_ref_r - lat_r
    dlon = lon_ref_r - lon_r

    a = np.sin(dlat / 2) ** 2 + np.cos(lat_r) * np.cos(lat_ref_r) * np.sin(dlon / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


# ---------------------------------------------------------------------------
# Analysis 1: Phase coherence and combined elevated phase set
# ---------------------------------------------------------------------------

def phase_coherence_analysis(df: pd.DataFrame, phases: np.ndarray) -> dict[str, Any]:
    """Identify top-3 bins per k and build the combined elevated phase set.

    The combined elevated phase set is built by finding regions of the phase
    axis [0, 1) that are covered by top-3 bins from >=2 of the 3 bin counts.
    Ranges from different bin counts partially overlap due to differing bin
    widths; this is handled by union-ing all top-3 intervals, then for each
    resulting sub-interval, counting how many k values have a top-3 bin that
    contains that sub-interval.

    Parameters
    ----------
    df : DataFrame with usgs_id and solar_phase
    phases : np.ndarray of solar phase values in [0, 1)

    Returns
    -------
    dict with phase_coherence sub-results and the elevated event DataFrame
    """
    n = len(phases)
    per_k: dict[str, Any] = {}
    # top3_intervals_by_k: list of (lo, hi) for each k
    top3_intervals_by_k: dict[int, list[tuple]] = {}

    for k in BIN_COUNTS:
        bin_indices = np.floor(phases * k).astype(int)
        counts = np.bincount(bin_indices, minlength=k)
        expected = n / k
        deviations = counts - expected

        top3_idx = np.argsort(deviations)[::-1][:TOP_N]
        top3_idx_sorted = sorted(top3_idx.tolist())
        top3_ranges = [(int(i) / k, (int(i) + 1) / k) for i in top3_idx_sorted]
        top3_devs = [float(deviations[i]) for i in top3_idx_sorted]

        per_k[f"k{k}"] = {
            "top3_bins": top3_idx_sorted,
            "top3_phase_ranges": top3_ranges,
            "top3_deviations": top3_devs,
        }
        top3_intervals_by_k[k] = top3_ranges
        logger.info("k=%d top-3 bins: %s  deviations: %s", k, top3_idx_sorted,
                    [round(d, 2) for d in top3_devs])

    # Build the set of phase-axis breakpoints from all top-3 interval boundaries
    breakpoints = set([0.0, 1.0])
    for intervals in top3_intervals_by_k.values():
        for lo, hi in intervals:
            breakpoints.add(lo)
            breakpoints.add(hi)
    breakpoints_sorted = sorted(breakpoints)

    # For each sub-interval between consecutive breakpoints, count how many k
    # values have a top-3 bin covering that sub-interval's midpoint
    def k_covers_point(x: float, k: int) -> bool:
        """Return True if any top-3 bin for this k covers point x."""
        for lo, hi in top3_intervals_by_k[k]:
            if lo <= x < hi:
                return True
        return False

    # Select sub-intervals covered by >=2 of 3 k values
    candidate_intervals = []
    for i in range(len(breakpoints_sorted) - 1):
        sub_lo = breakpoints_sorted[i]
        sub_hi = breakpoints_sorted[i + 1]
        midpoint = (sub_lo + sub_hi) / 2.0
        vote_count = sum(1 for k in BIN_COUNTS if k_covers_point(midpoint, k))
        if vote_count >= 2:
            candidate_intervals.append((sub_lo, sub_hi))

    logger.info("Sub-intervals covered by >=2 of 3 k: %s", candidate_intervals)

    # Merge contiguous candidate sub-intervals
    merged_intervals = _merge_intervals(candidate_intervals)
    logger.info("Merged combined elevated intervals: %s", merged_intervals)

    # Extract events in combined intervals, de-duplicate by usgs_id
    elevated_mask = np.zeros(n, dtype=bool)
    for lo, hi in merged_intervals:
        elevated_mask |= (phases >= lo) & (phases < hi)

    elevated_df = df[elevated_mask].copy()
    # De-duplicate by usgs_id (keep first occurrence)
    elevated_df = elevated_df.drop_duplicates(subset=["usgs_id"])
    n_elevated = len(elevated_df)

    # Expected pct under null: fraction of [0,1) covered by combined intervals
    total_interval_coverage = sum(hi - lo for lo, hi in merged_intervals)
    elevated_pct = n_elevated / n * 100
    expected_pct = total_interval_coverage * 100

    # Full catalog count in combined intervals (for baseline comparison)
    full_catalog_in_intervals = int(elevated_mask.sum())

    logger.info(
        "Combined elevated set: n=%d (%.1f%% of catalog); "
        "expected under null: %.1f%%",
        n_elevated, elevated_pct, expected_pct,
    )

    result: dict[str, Any] = {**per_k}
    result["combined_elevated_intervals"] = merged_intervals
    result["n_elevated"] = n_elevated
    result["elevated_pct_of_catalog"] = round(elevated_pct, 4)
    result["expected_pct_under_null"] = round(expected_pct, 4)
    result["full_catalog_in_intervals"] = full_catalog_in_intervals

    return result, elevated_df


def _merge_intervals(intervals: list[tuple]) -> list[tuple]:
    """Merge overlapping or adjacent intervals.

    Parameters
    ----------
    intervals : sorted list of (lo, hi) tuples

    Returns
    -------
    list of merged (lo, hi) tuples
    """
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for lo, hi in intervals[1:]:
        if lo <= merged[-1][1] + 1e-12:  # overlapping or adjacent
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [tuple(m) for m in merged]


# ---------------------------------------------------------------------------
# Analysis 2: Temporal distribution
# ---------------------------------------------------------------------------

def temporal_analysis(elevated_df: pd.DataFrame, full_df: pd.DataFrame,
                      rng: np.random.Generator) -> dict[str, Any]:
    """Compute IEI for elevated-bin events and compare to null distribution.

    Parameters
    ----------
    elevated_df : elevated-bin events
    full_df : full catalog
    rng : numpy random generator (fixed seed)

    Returns
    -------
    dict with temporal analysis results
    """
    n_elevated = len(elevated_df)

    # Sort by event_at
    elevated_sorted = elevated_df.copy()
    elevated_sorted["event_dt"] = pd.to_datetime(elevated_sorted["event_at"], utc=True)
    elevated_sorted = elevated_sorted.sort_values("event_dt")

    # Compute IEI in days
    diffs = elevated_sorted["event_dt"].diff().dropna()
    iei_days = diffs.dt.total_seconds() / 86400.0
    iei_days = iei_days.values

    elevated_iei = {
        "median": round(float(np.median(iei_days)), 4),
        "mean": round(float(np.mean(iei_days)), 4),
        "p10": round(float(np.percentile(iei_days, 10)), 4),
        "p90": round(float(np.percentile(iei_days, 90)), 4),
    }
    logger.info("Elevated IEI: median=%.2f  mean=%.2f  p10=%.2f  p90=%.2f",
                elevated_iei["median"], elevated_iei["mean"],
                elevated_iei["p10"], elevated_iei["p90"])

    # Null distribution: 1,000 random samples of size n_elevated
    full_df_sorted = full_df.copy()
    full_df_sorted["event_dt"] = pd.to_datetime(full_df_sorted["event_at"], utc=True)
    full_df_sorted = full_df_sorted.sort_values("event_dt")
    event_times = full_df_sorted["event_dt"].values  # sorted array of timestamps

    null_medians = []
    for _ in range(N_SAMPLES):
        idx = rng.choice(len(event_times), size=n_elevated, replace=False)
        sample_times = np.sort(event_times[idx])
        sample_diffs = np.diff(sample_times.astype(np.int64)) / 1e9 / 86400.0
        null_medians.append(float(np.median(sample_diffs)))

    null_medians_arr = np.array(null_medians)
    null_iei_ci = {
        "p2_5": round(float(np.percentile(null_medians_arr, 2.5)), 4),
        "p50": round(float(np.percentile(null_medians_arr, 50)), 4),
        "p97_5": round(float(np.percentile(null_medians_arr, 97.5)), 4),
    }
    logger.info("Null IEI 95%% CI: [%.2f, %.2f]",
                null_iei_ci["p2_5"], null_iei_ci["p97_5"])

    elevated_median_outside = bool(
        elevated_iei["median"] < null_iei_ci["p2_5"] or
        elevated_iei["median"] > null_iei_ci["p97_5"]
    )
    logger.info("Elevated IEI median outside null CI: %s", elevated_median_outside)

    return {
        "elevated_iei_days": elevated_iei,
        "null_iei_median_ci": null_iei_ci,
        "elevated_median_outside_null_ci": elevated_median_outside,
        "iei_values": iei_days.tolist(),  # stored for test validation
    }


# ---------------------------------------------------------------------------
# Analysis 3: Spatial distribution
# ---------------------------------------------------------------------------

def parse_plate_boundaries(boundary_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse PB2002 .dig file and return (lat_array, lon_array) of all boundary vertices.

    The .dig format has coordinates as (longitude, latitude) per line.
    Header lines and separator lines are skipped.
    """
    lons, lats = [], []
    with open(boundary_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "***" in line:
                continue
            # Try to parse as two numbers (lon, lat)
            parts = line.split(",")
            if len(parts) == 2:
                try:
                    lon_val = float(parts[0].strip())
                    lat_val = float(parts[1].strip())
                    lons.append(lon_val)
                    lats.append(lat_val)
                except ValueError:
                    pass  # header or unrecognized line
    lons_arr = np.array(lons)
    lats_arr = np.array(lats)
    logger.info("Parsed %d boundary vertices from %s", len(lons_arr), boundary_path)
    return lats_arr, lons_arr


def nearest_neighbor_analysis(elevated_df: pd.DataFrame, full_df: pd.DataFrame,
                               rng: np.random.Generator) -> dict[str, Any]:
    """Compute nearest-neighbor distances for elevated-bin events and null comparison.

    Parameters
    ----------
    elevated_df : elevated-bin events
    full_df : full catalog
    rng : numpy random generator

    Returns
    -------
    dict with NN statistics
    """
    n_elevated = len(elevated_df)
    elev_lat = elevated_df["latitude"].values
    elev_lon = elevated_df["longitude"].values

    # Compute pairwise NN distances within elevated set
    logger.info("Computing NN distances for %d elevated events...", n_elevated)
    nn_km = _compute_nn_distances(elev_lat, elev_lon)

    elev_nn_stats = {
        "p25": round(float(np.percentile(nn_km, 25)), 4),
        "p50": round(float(np.percentile(nn_km, 50)), 4),
        "p75": round(float(np.percentile(nn_km, 75)), 4),
    }
    logger.info("Elevated NN: p25=%.1f  p50=%.1f  p75=%.1f km",
                elev_nn_stats["p25"], elev_nn_stats["p50"], elev_nn_stats["p75"])

    # Null distribution
    full_lat = full_df["latitude"].values
    full_lon = full_df["longitude"].values
    null_nn_medians = []
    logger.info("Computing null NN distribution (%d samples)...", N_SAMPLES)
    for i in range(N_SAMPLES):
        if i % 100 == 0:
            logger.info("  NN null sample %d/%d", i, N_SAMPLES)
        idx = rng.choice(len(full_lat), size=n_elevated, replace=False)
        samp_lat = full_lat[idx]
        samp_lon = full_lon[idx]
        nn = _compute_nn_distances(samp_lat, samp_lon)
        null_nn_medians.append(float(np.median(nn)))

    null_nn_arr = np.array(null_nn_medians)
    null_nn_ci = {
        "p2_5": round(float(np.percentile(null_nn_arr, 2.5)), 4),
        "p50": round(float(np.percentile(null_nn_arr, 50)), 4),
        "p97_5": round(float(np.percentile(null_nn_arr, 97.5)), 4),
    }
    logger.info("Null NN 95%% CI: [%.1f, %.1f] km",
                null_nn_ci["p2_5"], null_nn_ci["p97_5"])

    elevated_nn_outside = bool(
        elev_nn_stats["p50"] < null_nn_ci["p2_5"] or
        elev_nn_stats["p50"] > null_nn_ci["p97_5"]
    )
    logger.info("Elevated NN median outside null CI: %s", elevated_nn_outside)

    return {
        "elevated_nn_km": elev_nn_stats,
        "null_nn_median_ci": null_nn_ci,
        "elevated_nn_outside_null_ci": elevated_nn_outside,
        "nn_values": nn_km.tolist(),  # stored for test validation
    }


def _compute_nn_distances(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Return nearest-neighbor distances (km) for a set of points.

    Uses chunked vectorized Haversine to manage memory.
    """
    n = len(lat)
    nn_dist = np.full(n, np.inf)
    chunk_size = 500

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        lat_chunk = lat[start:end]
        lon_chunk = lon[start:end]
        # Distance from chunk rows to all other points
        dist_mat = haversine_array(lat_chunk, lon_chunk, lat, lon)
        # Zero out self-distances
        for ci, gi in enumerate(range(start, end)):
            dist_mat[ci, gi] = np.inf
        nn_dist[start:end] = dist_mat.min(axis=1)

    return nn_dist


def boundary_proximity_analysis(elevated_df: pd.DataFrame, full_df: pd.DataFrame,
                                 boundary_path: Path) -> dict[str, Any]:
    """Classify events by distance to nearest plate boundary vertex.

    Parameters
    ----------
    elevated_df : elevated-bin events
    full_df : full catalog
    boundary_path : path to pb2002_boundaries.dig

    Returns
    -------
    dict with proximity classification percentages
    """
    bnd_lat, bnd_lon = parse_plate_boundaries(boundary_path)

    def classify_proximity(event_df: pd.DataFrame) -> dict[str, float]:
        """Classify events in df by boundary proximity."""
        e_lat = event_df["latitude"].values
        e_lon = event_df["longitude"].values
        n = len(e_lat)
        logger.info("Computing boundary proximity for %d events...", n)

        # Process in chunks to manage memory
        chunk_size = 200
        min_dists = np.full(n, np.inf)
        for start in range(0, n, chunk_size):
            end = min(start + chunk_size, n)
            lat_chunk = e_lat[start:end]
            lon_chunk = e_lon[start:end]
            dist_mat = haversine_array(lat_chunk, lon_chunk, bnd_lat, bnd_lon)
            min_dists[start:end] = dist_mat.min(axis=1)

        near = (min_dists <= NEAR_BOUNDARY_KM).sum()
        trans = ((min_dists > NEAR_BOUNDARY_KM) & (min_dists <= TRANSITIONAL_KM)).sum()
        intra = (min_dists > TRANSITIONAL_KM).sum()
        total = len(min_dists)

        return {
            "near_boundary_pct": round(float(near / total * 100), 4),
            "transitional_pct": round(float(trans / total * 100), 4),
            "intraplate_pct": round(float(intra / total * 100), 4),
        }

    logger.info("Classifying elevated-bin events by boundary proximity...")
    elevated_prox = classify_proximity(elevated_df)
    logger.info("Classifying full catalog by boundary proximity...")
    full_prox = classify_proximity(full_df)

    logger.info("Elevated proximity: %s", elevated_prox)
    logger.info("Full catalog proximity: %s", full_prox)

    return {
        "elevated": elevated_prox,
        "full_catalog": full_prox,
    }


def latitude_banding(elevated_df: pd.DataFrame, full_df: pd.DataFrame) -> dict[str, Any]:
    """Count events per 30-degree latitude band for both populations.

    Bands: 90S-60S, 60S-30S, 30S-0, 0-30N, 30N-60N, 60N-90N

    Parameters
    ----------
    elevated_df : elevated-bin events
    full_df : full catalog

    Returns
    -------
    dict with latitude band percentages for each population
    """
    bands = [
        ("90S-60S", -90, -60),
        ("60S-30S", -60, -30),
        ("30S-0", -30, 0),
        ("0-30N", 0, 30),
        ("30N-60N", 30, 60),
        ("60N-90N", 60, 90),
    ]

    def band_counts(df: pd.DataFrame) -> dict[str, float]:
        lat = df["latitude"].values
        total = len(lat)
        result = {}
        for name, lo, hi in bands:
            count = int(((lat >= lo) & (lat < hi)).sum())
            result[name] = round(count / total * 100, 4)
        return result

    return {
        "elevated": band_counts(elevated_df),
        "full_catalog": band_counts(full_df),
    }


# ---------------------------------------------------------------------------
# Analysis 4: Magnitude distribution
# ---------------------------------------------------------------------------

def magnitude_analysis(elevated_df: pd.DataFrame, full_df: pd.DataFrame) -> dict[str, Any]:
    """Compare magnitude distributions in 0.5-mag bins.

    Bins: 6.0-6.4, 6.5-6.9, 7.0-7.4, 7.5+

    Parameters
    ----------
    elevated_df : elevated-bin events
    full_df : full catalog

    Returns
    -------
    dict with magnitude bin counts for each population
    """
    def mag_bins(df: pd.DataFrame) -> dict[str, int]:
        mag = df["usgs_mag"].values
        return {
            "6.0-6.4": int(((mag >= 6.0) & (mag < 6.5)).sum()),
            "6.5-6.9": int(((mag >= 6.5) & (mag < 7.0)).sum()),
            "7.0-7.4": int(((mag >= 7.0) & (mag < 7.5)).sum()),
            "7.5+": int((mag >= 7.5).sum()),
        }

    elevated_mag = mag_bins(elevated_df)
    full_mag = mag_bins(full_df)

    logger.info("Elevated magnitude bins: %s", elevated_mag)
    logger.info("Full catalog magnitude bins: %s", full_mag)

    return {
        "elevated": elevated_mag,
        "full_catalog": full_mag,
    }


# ---------------------------------------------------------------------------
# Analysis 5: Declustering window estimation
# ---------------------------------------------------------------------------

def declustering_window_estimate(nn_km: dict[str, Any], iei_days: dict[str, Any]) -> dict[str, Any]:
    """Summarize observed clustering footprint vs G-K reference windows.

    Parameters
    ----------
    nn_km : nearest-neighbor stats (p25, p50, p75)
    iei_days : IEI stats (median, mean, p10, p90)

    Returns
    -------
    dict with G-K comparison and proposed window
    """
    p50_spatial = nn_km["p50"]
    p75_spatial = nn_km["p75"]
    p50_temporal = iei_days["median"]
    p75_temporal = iei_days["p90"]

    # Proposed window: use p75 of observed spatial (captures most clustering)
    # and p75 of IEI (p90) as a conservative temporal estimate
    proposed_spatial = round(p75_spatial, 1)
    proposed_temporal = round(p75_temporal, 1)

    basis = (
        f"Spatial: 75th percentile of elevated-bin nearest-neighbor distance "
        f"({proposed_spatial} km); "
        f"Temporal: 90th percentile of elevated-bin IEI ({proposed_temporal} days). "
        f"These values capture the observed clustering footprint of the elevated-bin "
        f"population; they are a data-informed reference point, not a validated "
        f"declustering algorithm."
    )

    return {
        "observed_spatial_p50_km": round(p50_spatial, 4),
        "observed_spatial_p75_km": round(p75_spatial, 4),
        "observed_temporal_p50_days": round(p50_temporal, 4),
        "observed_temporal_p75_days": round(p75_temporal, 4),
        "gk_reference_m6": GK_REFERENCE["m6"],
        "gk_reference_m7": GK_REFERENCE["m7"],
        "gk_reference_m8": GK_REFERENCE["m8"],
        "proposed_window": {
            "spatial_km": proposed_spatial,
            "temporal_days": proposed_temporal,
            "basis": basis,
        },
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_analysis() -> dict[str, Any]:
    """Execute the full A1b analysis and return the results dict."""
    logger.info("Loading data from %s", DATA_PATH)
    df = pd.read_csv(DATA_PATH)
    n = len(df)
    logger.info("Loaded %d events", n)
    assert n == 9210, f"Expected 9210 events, got {n}"

    # Compute solar phase
    logger.info("Computing solar phase normalization")
    df["solar_phase"] = compute_solar_phase(df)
    phases = df["solar_phase"].values

    rng = np.random.default_rng(RANDOM_SEED)

    # Analysis 1: Phase coherence
    logger.info("=== Analysis 1: Phase coherence ===")
    phase_coh, elevated_df = phase_coherence_analysis(df, phases)

    # Analysis 2: Temporal
    logger.info("=== Analysis 2: Temporal distribution ===")
    temporal = temporal_analysis(elevated_df, df, rng)

    # Analysis 3: Spatial
    logger.info("=== Analysis 3: Spatial distribution ===")
    rng2 = np.random.default_rng(RANDOM_SEED + 1)
    nn_results = nearest_neighbor_analysis(elevated_df, df, rng2)
    bnd_prox = boundary_proximity_analysis(elevated_df, df, BOUNDARY_PATH)
    lat_bands = latitude_banding(elevated_df, df)

    # Analysis 4: Magnitude
    logger.info("=== Analysis 4: Magnitude distribution ===")
    mag_results = magnitude_analysis(elevated_df, df)

    # Analysis 5: Declustering window estimate
    logger.info("=== Analysis 5: Declustering window estimate ===")
    declust = declustering_window_estimate(nn_results["elevated_nn_km"],
                                           temporal["elevated_iei_days"])

    # Assemble output
    output: dict[str, Any] = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "data_source": "data/global-sets/iscgem_global_events.csv",
        "boundary_file": "lib/pb2002_boundaries.dig",
        "event_count": n,
        "phase_coherence": phase_coh,
        "temporal": {
            "elevated_iei_days": temporal["elevated_iei_days"],
            "null_iei_median_ci": temporal["null_iei_median_ci"],
            "elevated_median_outside_null_ci": temporal["elevated_median_outside_null_ci"],
        },
        "spatial": {
            "elevated_nn_km": nn_results["elevated_nn_km"],
            "null_nn_median_ci": nn_results["null_nn_median_ci"],
            "elevated_nn_outside_null_ci": nn_results["elevated_nn_outside_null_ci"],
            "boundary_proximity": bnd_prox,
            "latitude_bands": lat_bands,
        },
        "magnitude": mag_results,
        "declustering_window_estimate": declust,
        # Store raw arrays for tests (not in the spec JSON schema but useful)
        "_iei_count": len(temporal["iei_values"]),
        "_nn_count": len(nn_results["nn_values"]),
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(output, f, indent=2)
    logger.info("Results written to %s", RESULTS_PATH)

    return output


if __name__ == "__main__":
    run_analysis()
