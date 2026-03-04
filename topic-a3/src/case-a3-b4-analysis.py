"""
Case A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation

Three-component analysis to disentangle confounds in A2.B4's mid-crustal solar signal:
1. Two-way depth × magnitude stratification
2. CRUST1.0 Moho isolation sub-test
3. Subduction proximity cross-tabulation within mid-crustal band

CRUST1.0 data from https://igppweb.ucsd.edu/~gabi/crust1.html (Laske et al. 2013)
"""

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats
from scipy.spatial import cKDTree
from scipy.interpolate import RectBivariateSpline
import pyproj

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
LIB_DIR = BASE_DIR.parent / "lib"

RAW_PATH    = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH  = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
FOCAL_PATH  = BASE_DIR.parent / "data" / "iscgem" / "focal-mechanism" / "focal_join_global.csv"
STEPS_PATH  = LIB_DIR / "PB2002_steps.dat"
CRUST1_PATH = LIB_DIR / "crust1.bnds"

OUTPUT_PATH = BASE_DIR / "output" / "case-a3-b4-results.json"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
K_BINS: int = 24
JULIAN_YEAR_SECS: float = 31_557_600.0
SUB_PROXIMITY_THRESHOLD_KM: float = 200.0
MOHO_DELTAS_KM: list[float] = [5.0, 10.0, 15.0]

MAG_BANDS: list[dict[str, Any]] = [
    {"label": "m6_6.9",  "min": 6.0, "max": 7.0},
    {"label": "m7_7.9",  "min": 7.0, "max": 8.0},
    {"label": "m8_plus", "min": 8.0, "max": 99.0},
]
DEPTH_BANDS: list[dict[str, Any]] = [
    {"label": "shallow_0-20km",        "min": 0,   "max": 20},
    {"label": "midcrustal_20-70km",    "min": 20,  "max": 70},
    {"label": "intermediate_70-300km", "min": 70,  "max": 300},
    {"label": "deep_300km+",           "min": 300, "max": 9999},
]
INTERVAL_BINS: dict[str, list[int]] = {
    "interval_1": [4, 5],
    "interval_2": [15],
    "interval_3": [21],
}


# ---------------------------------------------------------------------------
# Section 1 helper: adaptive k selection
# ---------------------------------------------------------------------------

def select_k(n: int) -> int | None:
    """
    Select adaptive bin count based on sample size.

    Returns:
        24 if n >= 500; 16 if 200 <= n < 500; 12 if 100 <= n < 200;
        None if n < 100 (low-n flag; caller should use k=12 for computation).
    """
    if n >= 500:
        return 24
    if n >= 200:
        return 16
    if n >= 100:
        return 12
    return None  # low-n; flag but still compute at k=12 for record


# ---------------------------------------------------------------------------
# Section 2: CRUST1.0 Moho depth loading and assignment
# CRUST1.0 data from https://igppweb.ucsd.edu/~gabi/crust1.html (Laske et al. 2013)
# ---------------------------------------------------------------------------

def load_crust1_moho(
    crust1_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, RectBivariateSpline]:
    """
    Load CRUST1.0 boundary file and build Moho interpolation surface.

    The file has 64,800 rows × 9 columns. Grid order: lat 89.5→-89.5°N (step -1),
    lon -179.5→179.5°E (step +1). Column index 8 (0-based) is Moho depth in km,
    stored as negative values (depth below sea level). We negate to get positive km.

    Parameters
    ----------
    crust1_path : Path to crust1.bnds file

    Returns
    -------
    tuple of (moho_grid shape (180,360), lats (180,), lons (360,), moho_interp)
    """
    logger.info(f"Loading CRUST1.0 from {crust1_path}...")
    raw = np.loadtxt(crust1_path)  # shape (64800, 9)
    assert raw.shape == (64800, 9), f"Unexpected CRUST1.0 shape: {raw.shape}"

    bnds = raw.reshape(180, 360, 9)
    # Layer 8 (0-based) = Moho; values are negative (depth below sea level)
    # Negate to get positive depth in km
    moho_grid = -bnds[:, :, 8]  # shape (180, 360), positive km

    # Grid cell centers
    lats = np.arange(89.5, -90.5, -1.0)   # 89.5 → -89.5 (descending)
    lons = np.arange(-179.5, 180.5, 1.0)  # -179.5 → 179.5

    logger.info(
        f"CRUST1.0 Moho grid: min={moho_grid.min():.2f} km, "
        f"max={moho_grid.max():.2f} km, mean={moho_grid.mean():.2f} km"
    )

    # RectBivariateSpline requires ascending lat axis; flip
    moho_interp = RectBivariateSpline(
        lats[::-1], lons, moho_grid[::-1], kx=1, ky=1
    )

    return moho_grid, lats, lons, moho_interp


def assign_moho_depth(
    event_lats: np.ndarray,
    event_lons: np.ndarray,
    moho_interp: RectBivariateSpline,
) -> np.ndarray:
    """
    Assign per-event Moho depth (km) via bilinear interpolation of CRUST1.0 grid.

    Clips result to [3.0, 90.0] km to handle edge cells.

    Parameters
    ----------
    event_lats : array of event latitudes (degrees)
    event_lons : array of event longitudes (degrees)
    moho_interp : RectBivariateSpline interpolator

    Returns
    -------
    np.ndarray of shape (n_events,) with per-event Moho depth in km
    """
    # RectBivariateSpline.ev() evaluates at paired (xi, yi) points vectorized
    moho_depths = moho_interp.ev(event_lats, event_lons)

    # Clip to plausible range
    moho_depths = np.clip(moho_depths, 3.0, 90.0)

    logger.info(
        f"Moho assignment: mean={moho_depths.mean():.2f} km, "
        f"median={np.median(moho_depths):.2f} km, "
        f"std={moho_depths.std():.2f} km"
    )
    return moho_depths


# ---------------------------------------------------------------------------
# Section 3: Subduction zone proximity (reuse B3 pipeline exactly)
# ---------------------------------------------------------------------------

def parse_sub_boundaries(steps_path: Path) -> tuple[np.ndarray, int]:
    """
    Parse PB2002_steps.dat; return sample points for SUB and OCB segments.

    Generates 3 points per segment: endpoint1, endpoint2, midpoint.

    Returns
    -------
    tuple of (np.ndarray shape (N_points, 2) [lon, lat], n_sub_segments int)
    """
    points: list[tuple[float, float]] = []
    n_sub_segments: int = 0
    n_ocb_segments: int = 0

    with open(steps_path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue

            btype = parts[-1].lstrip(":")
            if btype not in ("SUB", "OCB"):
                continue

            try:
                lon1 = float(parts[2].lstrip(":"))
                lat1 = float(parts[3].lstrip(":"))
                lon2 = float(parts[4].lstrip(":"))
                lat2 = float(parts[5].lstrip(":"))
            except (ValueError, IndexError):
                continue

            mid_lon = (lon1 + lon2) / 2.0
            mid_lat = (lat1 + lat2) / 2.0

            points.append((lon1, lat1))
            points.append((lon2, lat2))
            points.append((mid_lon, mid_lat))

            if btype == "SUB":
                n_sub_segments += 1
            else:
                n_ocb_segments += 1

    result = np.array(points, dtype=float)
    logger.info(
        f"Parsed {n_sub_segments} SUB segments + {n_ocb_segments} OCB segments "
        f"({len(result)} total sample points) from PB2002_steps.dat"
    )
    return result, n_sub_segments


def latlon_to_unit_sphere(lon_deg: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    """
    Convert lon/lat arrays (degrees) to unit-sphere 3D Cartesian.

    Returns np.ndarray shape (N, 3): [x, y, z].
    """
    lon_r = np.radians(lon_deg)
    lat_r = np.radians(lat_deg)
    x = np.cos(lat_r) * np.cos(lon_r)
    y = np.cos(lat_r) * np.sin(lon_r)
    z = np.sin(lat_r)
    return np.column_stack([x, y, z])


def compute_subduction_distances(
    event_lons: np.ndarray,
    event_lats: np.ndarray,
    sub_boundary_lonlat: np.ndarray,
) -> np.ndarray:
    """
    Compute geodesic distance (km) from each event to nearest SUB/OCB point.

    Uses unit-sphere cKDTree (k=5 neighbors) then pyproj WGS84 refinement.

    Returns np.ndarray of shape (n_events,).
    """
    bnd_lons = sub_boundary_lonlat[:, 0]
    bnd_lats = sub_boundary_lonlat[:, 1]
    bnd_xyz = latlon_to_unit_sphere(bnd_lons, bnd_lats)
    evt_xyz = latlon_to_unit_sphere(event_lons, event_lats)

    tree = cKDTree(bnd_xyz)
    k = min(5, len(bnd_xyz))
    _, indices = tree.query(evt_xyz, k=k)

    geod = pyproj.Geod(ellps="WGS84")
    n_events = len(event_lons)
    dist_to_sub = np.empty(n_events, dtype=float)

    for i in range(n_events):
        min_dist_km = np.inf
        idx_arr = indices[i] if k > 1 else [indices[i]]
        for idx in idx_arr:
            cand_lon = bnd_lons[idx]
            cand_lat = bnd_lats[idx]
            _, _, dist_m = geod.inv(event_lons[i], event_lats[i], cand_lon, cand_lat)
            dist_km = abs(dist_m) / 1000.0
            if dist_km < min_dist_km:
                min_dist_km = dist_km
        dist_to_sub[i] = min_dist_km

    logger.info(
        f"Subduction distances: mean={dist_to_sub.mean():.1f} km, "
        f"median={np.median(dist_to_sub):.1f} km, "
        f"frac within {SUB_PROXIMITY_THRESHOLD_KM} km="
        f"{(dist_to_sub <= SUB_PROXIMITY_THRESHOLD_KM).mean():.3f}"
    )
    return dist_to_sub


# ---------------------------------------------------------------------------
# Section 4: Chi-square statistics
# ---------------------------------------------------------------------------

def compute_chi2_stats(phases: np.ndarray, k: int) -> dict[str, Any]:
    """
    Compute chi-square uniformity test, Cramér's V, and interval z-scores.

    Parameters
    ----------
    phases : array of phase values in [0, 1)
    k : number of bins

    Returns
    -------
    dict with keys: n, k, low_n, chi2, p_chi2, cramers_v, bin_counts,
                    interval_1_z, interval_2_z, interval_3_z
    """
    n = len(phases)
    low_n = n < 100

    if n == 0:
        return {
            "n": 0, "k": k, "low_n": True,
            "chi2": 0.0, "p_chi2": 1.0, "cramers_v": 0.0,
            "bin_counts": [],
            "interval_1_z": None, "interval_2_z": None, "interval_3_z": None,
        }

    observed = np.bincount(
        (np.floor(phases * k).astype(int)) % k,
        minlength=k,
    )
    expected_per_bin = n / k
    chi2_stat, p_chi2 = scipy.stats.chisquare(observed, np.full(k, expected_per_bin))
    cramers_v = np.sqrt(chi2_stat / (n * (k - 1)))

    # Per-interval z-scores (only meaningful at k=24)
    if k == 24:
        interval_z: dict[str, float | None] = {}
        for iname, bins in INTERVAL_BINS.items():
            count = int(observed[bins].sum())
            expected = len(bins) * expected_per_bin
            z = (count - expected) / np.sqrt(len(bins) * expected_per_bin)
            interval_z[iname] = float(z)
    else:
        interval_z = {"interval_1": None, "interval_2": None, "interval_3": None}

    return {
        "n": int(n),
        "k": int(k),
        "low_n": bool(low_n),
        "chi2": float(chi2_stat),
        "p_chi2": float(p_chi2),
        "cramers_v": float(cramers_v),
        "bin_counts": observed.tolist(),
        "interval_1_z": interval_z["interval_1"],
        "interval_2_z": interval_z["interval_2"],
        "interval_3_z": interval_z["interval_3"],
    }


# ---------------------------------------------------------------------------
# Section 5: Two-way stratification analysis
# ---------------------------------------------------------------------------

def run_stratification_matrix(df: pd.DataFrame) -> tuple[dict, dict]:
    """
    Compute chi-square statistics for each depth × magnitude cell.

    Also computes depth-band totals (all magnitudes) as regression anchors.

    Returns
    -------
    tuple of (stratification_matrix dict, depth_band_totals dict)
    """
    stratification_matrix: dict[str, Any] = {}
    depth_band_totals: dict[str, Any] = {}

    # Filter to events with valid depth
    df_valid = df[df["depth"].notna()].copy()

    for db in DEPTH_BANDS:
        d_label = db["label"]
        d_min, d_max = db["min"], db["max"]

        depth_subset = df_valid[
            (df_valid["depth"] >= d_min) & (df_valid["depth"] < d_max)
        ]

        # Depth-band total (all magnitudes combined, k=24)
        total_k = 24
        total_stats = compute_chi2_stats(depth_subset["phase"].values, total_k)
        depth_band_totals[d_label] = {
            "n": total_stats["n"],
            "k": total_k,
            "chi2": total_stats["chi2"],
            "p_chi2": total_stats["p_chi2"],
            "cramers_v": total_stats["cramers_v"],
        }
        logger.info(
            f"Depth total {d_label}: n={total_stats['n']}, "
            f"chi2={total_stats['chi2']:.2f}, p={total_stats['p_chi2']:.4e}, "
            f"V={total_stats['cramers_v']:.4f}"
        )

        stratification_matrix[d_label] = {}

        for mb in MAG_BANDS:
            m_label = mb["label"]
            m_min, m_max = mb["min"], mb["max"]

            cell_df = depth_subset[
                (depth_subset["usgs_mag"] >= m_min) & (depth_subset["usgs_mag"] < m_max)
            ]
            n_cell = len(cell_df)

            k_sel = select_k(n_cell)
            low_n = k_sel is None
            k_use = k_sel if k_sel is not None else 12

            stats = compute_chi2_stats(cell_df["phase"].values, k_use)
            # Ensure low_n reflects spec rule (n < 100)
            stats["low_n"] = low_n

            stratification_matrix[d_label][m_label] = stats
            logger.info(
                f"Cell {d_label} × {m_label}: n={n_cell}, k={k_use}, "
                f"chi2={stats['chi2']:.2f}, p={stats['p_chi2']:.4e}, "
                f"V={stats['cramers_v']:.4f}"
                + (" [LOW-N]" if low_n else "")
            )

    return stratification_matrix, depth_band_totals


# ---------------------------------------------------------------------------
# Section 6: Moho isolation analysis
# ---------------------------------------------------------------------------

def run_moho_isolation(df: pd.DataFrame) -> dict[str, Any]:
    """
    Compute chi-square statistics for Moho-proximal events at each delta value.

    Splits results by ocean_class (continental, oceanic, transitional).

    Returns
    -------
    dict keyed by "delta_5km", "delta_10km", "delta_15km"
    """
    moho_isolation: dict[str, Any] = {}
    df_valid = df[df["depth"].notna()].copy()

    for delta in MOHO_DELTAS_KM:
        col = f"moho_proximal_{int(delta)}km"
        moho_df = df_valid[df_valid[col]].copy()
        n_total = len(moho_df)
        logger.info(f"Moho isolation delta={delta} km: n={n_total}")

        k_sel = select_k(n_total)
        k_use = k_sel if k_sel is not None else 12

        overall_stats = compute_chi2_stats(moho_df["phase"].values, k_use)

        # Split by ocean_class
        sub_results: dict[str, Any] = {}
        for class_label in ["continental", "oceanic", "transitional"]:
            class_df = moho_df[moho_df["ocean_class"] == class_label]
            n_class = len(class_df)
            k_class_sel = select_k(n_class)
            k_class = k_class_sel if k_class_sel is not None else 12
            class_stats = compute_chi2_stats(class_df["phase"].values, k_class)
            sub_results[class_label] = class_stats
            logger.info(
                f"  {class_label}: n={n_class}, k={k_class}, "
                f"chi2={class_stats['chi2']:.2f}, p={class_stats['p_chi2']:.4e}"
            )

        delta_key = f"delta_{int(delta)}km"
        moho_isolation[delta_key] = {
            "n_total": n_total,
            "k": k_use,
            "chi2": overall_stats["chi2"],
            "p_chi2": overall_stats["p_chi2"],
            "cramers_v": overall_stats["cramers_v"],
            "bin_counts": overall_stats["bin_counts"],
            "interval_1_z": overall_stats["interval_1_z"],
            "interval_2_z": overall_stats["interval_2_z"],
            "interval_3_z": overall_stats["interval_3_z"],
            "continental": sub_results["continental"],
            "oceanic": sub_results["oceanic"],
            "transitional": sub_results["transitional"],
        }

    return moho_isolation


# ---------------------------------------------------------------------------
# Section 7: Subduction proximity cross-tabulation (mid-crustal only)
# ---------------------------------------------------------------------------

def run_subduction_crosstab(df: pd.DataFrame) -> dict[str, Any]:
    """
    Cross-tabulate solar-phase signal by near/far subduction within mid-crustal band.

    Tests all magnitudes combined and each magnitude band separately.

    Returns
    -------
    dict with keys: all_magnitudes, m6_6.9, m7_7.9, m8_plus
    """
    df_valid = df[df["depth"].notna()].copy()

    # Mid-crustal band
    mid_df = df_valid[
        (df_valid["depth"] >= 20) & (df_valid["depth"] < 70)
    ].copy()

    crosstab: dict[str, Any] = {}

    # All magnitudes
    near_all = mid_df[mid_df["near_subduction"]]
    far_all  = mid_df[~mid_df["near_subduction"]]
    k_near = select_k(len(near_all)) or 12
    k_far  = select_k(len(far_all))  or 12
    crosstab["all_magnitudes"] = {
        "near_sub": compute_chi2_stats(near_all["phase"].values, k_near),
        "far_sub":  compute_chi2_stats(far_all["phase"].values, k_far),
    }
    logger.info(
        f"Subduction crosstab all-mag: near n={len(near_all)} p={crosstab['all_magnitudes']['near_sub']['p_chi2']:.4e}, "
        f"far n={len(far_all)} p={crosstab['all_magnitudes']['far_sub']['p_chi2']:.4e}"
    )

    # Per magnitude band
    for mb in MAG_BANDS:
        m_label = mb["label"]
        m_min, m_max = mb["min"], mb["max"]
        mag_df = mid_df[(mid_df["usgs_mag"] >= m_min) & (mid_df["usgs_mag"] < m_max)]

        near_mag = mag_df[mag_df["near_subduction"]]
        far_mag  = mag_df[~mag_df["near_subduction"]]
        k_near_m = select_k(len(near_mag)) or 12
        k_far_m  = select_k(len(far_mag))  or 12

        crosstab[m_label] = {
            "near_sub": compute_chi2_stats(near_mag["phase"].values, k_near_m),
            "far_sub":  compute_chi2_stats(far_mag["phase"].values, k_far_m),
        }
        logger.info(
            f"Subduction crosstab {m_label}: near n={len(near_mag)}, "
            f"far n={len(far_mag)}"
        )

    return crosstab


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation."""

    # -----------------------------------------------------------------------
    # Section 1: Load and merge data
    # -----------------------------------------------------------------------
    logger.info("Loading raw catalog...")
    raw_df = pd.read_csv(RAW_PATH, parse_dates=["event_at"])
    assert len(raw_df) == 9210, f"Expected 9210 rows, got {len(raw_df)}"
    logger.info(f"Raw catalog: {len(raw_df)} rows")

    logger.info("Loading GSHHG classification...")
    gshhg_df = pd.read_csv(GSHHG_PATH)

    logger.info("Loading focal join...")
    focal_df = pd.read_csv(FOCAL_PATH)

    # Merge GSHHG
    df = raw_df.merge(
        gshhg_df[["usgs_id", "ocean_class", "dist_to_coast_km"]],
        on="usgs_id", how="left"
    )
    # Merge focal (mechanism, match_confidence)
    df = df.merge(
        focal_df[["usgs_id", "mechanism", "rake", "match_confidence"]],
        on="usgs_id", how="left"
    )

    assert df["dist_to_coast_km"].isna().sum() == 0, "NaN in dist_to_coast_km after merge"
    logger.info("All merges completed successfully")

    # Compute solar phase
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0

    # Check depth nulls
    n_depth_null = int(df["depth"].isna().sum())
    logger.info(f"Depth null count: {n_depth_null}")

    # -----------------------------------------------------------------------
    # Section 2: CRUST1.0 Moho depth assignment
    # -----------------------------------------------------------------------
    moho_grid, lats, lons, moho_interp = load_crust1_moho(CRUST1_PATH)

    # Assign to all events (use lat/lon from full df)
    logger.info("Assigning Moho depth to all events...")
    moho_depths = assign_moho_depth(
        df["latitude"].values,
        df["longitude"].values,
        moho_interp,
    )
    df["moho_depth_km"] = moho_depths

    n_events_assigned = int(df["moho_depth_km"].notna().sum())
    mean_moho = float(df["moho_depth_km"].mean())
    median_moho = float(df["moho_depth_km"].median())
    std_moho = float(df["moho_depth_km"].std())

    # Moho proximity flags
    n_proximal_by_delta: dict[str, int] = {}
    for delta in MOHO_DELTAS_KM:
        col = f"moho_proximal_{int(delta)}km"
        # Only flag events with valid depth
        df[col] = (
            df["depth"].notna() &
            (np.abs(df["depth"] - df["moho_depth_km"]) <= delta)
        )
        n_prox = int(df[col].sum())
        n_proximal_by_delta[str(int(delta))] = n_prox
        logger.info(f"Moho proximal at delta={delta} km: n={n_prox}")

    moho_assignment_info = {
        "n_events_assigned": n_events_assigned,
        "mean_moho_depth_km": mean_moho,
        "median_moho_depth_km": median_moho,
        "std_moho_depth_km": std_moho,
        "n_proximal_by_delta": n_proximal_by_delta,
    }

    # -----------------------------------------------------------------------
    # Section 3: Subduction zone proximity
    # -----------------------------------------------------------------------
    logger.info("Parsing SUB+OCB boundary points from PB2002_steps.dat...")
    sub_points, n_sub_segments = parse_sub_boundaries(STEPS_PATH)

    logger.info("Computing per-event distance to nearest subduction boundary...")
    dist_to_sub = compute_subduction_distances(
        df["longitude"].values,
        df["latitude"].values,
        sub_points,
    )
    df["dist_to_subduction_km"] = dist_to_sub
    df["near_subduction"] = df["dist_to_subduction_km"] <= SUB_PROXIMITY_THRESHOLD_KM

    n_near_sub = int(df["near_subduction"].sum())
    pct_near_sub = float(n_near_sub / len(df))
    mean_dist_sub = float(dist_to_sub.mean())

    subduction_proximity_info = {
        "n_sub_boundary_points": n_sub_segments,
        "n_near_subduction": n_near_sub,
        "pct_near_subduction": pct_near_sub,
        "mean_dist_to_sub_km": mean_dist_sub,
    }

    # -----------------------------------------------------------------------
    # Section 5: Two-way stratification matrix
    # -----------------------------------------------------------------------
    logger.info("Running two-way stratification matrix...")
    stratification_matrix, depth_band_totals = run_stratification_matrix(df)

    # -----------------------------------------------------------------------
    # Section 6: Moho isolation
    # -----------------------------------------------------------------------
    logger.info("Running Moho isolation analysis...")
    moho_isolation = run_moho_isolation(df)

    # -----------------------------------------------------------------------
    # Section 7: Subduction proximity cross-tabulation
    # -----------------------------------------------------------------------
    logger.info("Running subduction proximity cross-tabulation (mid-crustal)...")
    subduction_crosstab = run_subduction_crosstab(df)

    # -----------------------------------------------------------------------
    # Section 8: Write results JSON
    # -----------------------------------------------------------------------
    results: dict[str, Any] = {
        "case": "A3.B4",
        "title": "Depth × Magnitude Two-Way Stratification with Moho Isolation",
        "parameters": {
            "n_catalog": 9210,
            "n_depth_null": n_depth_null,
            "k_bins_primary": 24,
            "adaptive_k_thresholds": {"k24": 500, "k16": 200, "k12": 100},
            "moho_deltas_km": MOHO_DELTAS_KM,
            "sub_proximity_threshold_km": SUB_PROXIMITY_THRESHOLD_KM,
            "julian_year_secs": JULIAN_YEAR_SECS,
            "crust1_source": "https://igppweb.ucsd.edu/~gabi/crust1.html",
        },
        "moho_assignment": moho_assignment_info,
        "subduction_proximity": subduction_proximity_info,
        "depth_band_totals": depth_band_totals,
        "stratification_matrix": stratification_matrix,
        "moho_isolation": moho_isolation,
        "subduction_crosstab_midcrustal": subduction_crosstab,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info(f"Results written to {OUTPUT_PATH}")

    # Save merged dataframe for visualization use
    df_out_path = BASE_DIR / "output" / "case-a3-b4-events.pkl"
    df.to_pickle(df_out_path)
    logger.info(f"Event dataframe written to {df_out_path}")

    # Save moho grid for visualization
    np.save(str(BASE_DIR / "output" / "case-a3-b4-moho-grid.npy"), moho_grid)
    np.save(str(BASE_DIR / "output" / "case-a3-b4-sub-points.npy"), sub_points)
    logger.info("Moho grid and sub boundary points saved.")


if __name__ == "__main__":
    main()
