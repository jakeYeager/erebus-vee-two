"""
Case A3.B5: Corrected Null-Distribution Geometric Variable Test

Tests solar_declination, declination_rate, earth_sun_distance against solar_phase
using time-weighted analytic null distributions rather than a uniform null.
Corrects the methodological error in A2.B5 where non-cyclic variables were tested
against a uniform null, which is incorrect for variables with non-uniform temporal
distributions.

Strata: full, continental, mid-crustal 20-70km, continental × mid-crustal.
Adaptive k: 24 if n>=500, 16 if 200<=n<500, 12 if 100<=n<200, skip if n<100.
"""

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent

SOLAR_PATH = BASE_DIR.parent / "data" / "iscgem" / "celestial-geometry" / "solar_geometry_global.csv"
GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
OUTPUT_PATH = BASE_DIR / "output" / "case-a3-b5-results.json"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
JULIAN_YEAR_SECS: float = 31_557_600.0
MIDCRUSTAL_MIN_KM: float = 20.0
MIDCRUSTAL_MAX_KM: float = 70.0
VARIABLES: list[str] = ["solar_phase", "solar_declination", "declination_rate", "earth_sun_distance"]

STRATA: dict[str, str | None] = {
    "full": None,
    "continental": "ocean_class == 'continental'",
    "midcrustal": "20 <= depth < 70",
    "continental_midcrustal": "ocean_class == 'continental' and 20 <= depth < 70",
}


# ---------------------------------------------------------------------------
# Section 2: Analytic null distribution generation
# ---------------------------------------------------------------------------

def generate_analytic_null(
    year_start: int,
    year_end: int,
    k: int,
    var_ranges: dict[str, tuple[float, float]] | None = None,
) -> dict[str, np.ndarray]:
    """
    Generate time-weighted expected bin fractions for solar_declination,
    declination_rate, and earth_sun_distance using a dense daily analytic model.

    Parameters
    ----------
    year_start : int
        First year of catalog coverage (inclusive).
    year_end : int
        Last year of catalog coverage (inclusive).
    k : int
        Number of bins.
    var_ranges : dict, optional
        Maps variable name to (min, max) tuple from actual data.
        If None, uses synthetic series range (not recommended for chi-square).

    Returns
    -------
    dict with keys "solar_declination", "declination_rate", "earth_sun_distance",
    each mapping to a np.ndarray of shape (k,) representing the expected fraction
    of time in each bin (sums to 1.0).
    """
    # Generate daily time steps spanning the catalog coverage period
    n_days = int((year_end - year_start + 1) * 365.25)
    days = np.arange(0, n_days, dtype=float)

    # J2000.0 offset: days since Jan 1.5 2000 (J2000.0 epoch)
    D = days + (year_start - 2000) * 365.25 - 0.5

    # Solar mean longitude (degrees) → radians
    L = np.radians((280.46 + 0.9856474 * D) % 360)

    # Mean anomaly (degrees) → radians
    g = np.radians((357.528 + 0.9856003 * D) % 360)

    # Ecliptic longitude with first-order aberration correction
    lam = L + np.radians(1.915 * np.sin(g) + 0.020 * np.sin(2 * g))

    # Solar declination (degrees)
    dec = np.degrees(np.arcsin(np.sin(np.radians(23.439)) * np.sin(lam)))

    # Earth-Sun distance (AU) — Kepler approximation
    dist = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)

    # Declination rate (degrees/day) via finite difference
    rate = np.gradient(dec, days)

    logger.info(
        f"Synthetic series: n={len(days)} points, "
        f"dec range=[{dec.min():.4f}, {dec.max():.4f}] deg, "
        f"dist range=[{dist.min():.6f}, {dist.max():.6f}] AU, "
        f"rate range=[{rate.min():.4f}, {rate.max():.4f}] deg/day"
    )

    result: dict[str, np.ndarray] = {}

    for var_name, synthetic_values in [
        ("solar_declination", dec),
        ("declination_rate", rate),
        ("earth_sun_distance", dist),
    ]:
        if var_ranges is not None and var_name in var_ranges:
            vmin, vmax = var_ranges[var_name]
        else:
            vmin = float(synthetic_values.min())
            vmax = float(synthetic_values.max())

        edges = np.linspace(vmin, vmax, k + 1)
        null_counts, _ = np.histogram(synthetic_values, bins=edges)
        # Ensure no division by zero
        total = null_counts.sum()
        if total == 0:
            null_fractions = np.full(k, 1.0 / k)
        else:
            null_fractions = null_counts / total
        result[var_name] = null_fractions

    return result


def get_null_qc_stats(year_start: int = 1950, year_end: int = 2021) -> dict[str, Any]:
    """
    Return QC statistics for the synthetic series.

    Parameters
    ----------
    year_start : int
    year_end : int

    Returns
    -------
    dict with n_synthetic_points, dec_synthetic_range, dist_synthetic_range
    """
    n_days = int((year_end - year_start + 1) * 365.25)
    days = np.arange(0, n_days, dtype=float)
    D = days + (year_start - 2000) * 365.25 - 0.5
    L = np.radians((280.46 + 0.9856474 * D) % 360)
    g = np.radians((357.528 + 0.9856003 * D) % 360)
    lam = L + np.radians(1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    dec = np.degrees(np.arcsin(np.sin(np.radians(23.439)) * np.sin(lam)))
    dist = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    return {
        "n_synthetic_points": int(len(days)),
        "dec_synthetic_range": [float(dec.min()), float(dec.max())],
        "dist_synthetic_range": [float(dist.min()), float(dist.max())],
    }


# ---------------------------------------------------------------------------
# Section 3: Per-event binning
# ---------------------------------------------------------------------------

def assign_bins(
    df: pd.DataFrame,
    k: int,
    var_ranges: dict[str, tuple[float, float]],
) -> pd.DataFrame:
    """
    Assign bin indices for all four variables using equal-width bins over observed ranges.

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with solar_phase, solar_declination, declination_rate, earth_sun_distance.
    k : int
        Number of bins.
    var_ranges : dict
        Maps variable name to (min, max) tuple from actual data.

    Returns
    -------
    DataFrame with added columns bin_solar_phase, bin_solar_declination,
    bin_declination_rate, bin_earth_sun_distance.
    """
    df = df.copy()

    # solar_phase: cyclic [0, 1)
    df["bin_solar_phase"] = np.floor(df["solar_phase"] * k).astype(int).clip(0, k - 1)

    # Non-cyclic variables: equal-width bins over observed range
    for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
        vmin, vmax = var_ranges[var]
        bin_col = f"bin_{var}"
        df[bin_col] = (
            np.floor((df[var] - vmin) / (vmax - vmin) * k).astype(int).clip(0, k - 1)
        )

    return df


# ---------------------------------------------------------------------------
# Section 4: Corrected chi-square statistics
# ---------------------------------------------------------------------------

def compute_corrected_chi2(
    bin_col: np.ndarray,
    k: int,
    null_fractions: np.ndarray,
    n: int,
) -> dict[str, Any]:
    """
    Chi-square test using time-weighted expected counts.

    Parameters
    ----------
    bin_col : np.ndarray
        Integer bin assignments for each event.
    k : int
        Number of bins.
    null_fractions : np.ndarray
        Expected fraction of time in each bin (shape k, sums to 1.0).
    n : int
        Total events in the subset.

    Returns
    -------
    dict with: n, k, chi2_corrected, p_corrected, cramers_v_corrected,
               chi2_uniform, p_uniform, cramers_v_uniform,
               bin_counts, expected_corrected, expected_uniform.
    """
    if n == 0:
        return {
            "n": 0, "k": k, "low_n": True,
            "chi2_corrected": 0.0, "p_corrected": 1.0, "cramers_v_corrected": 0.0,
            "chi2_uniform": 0.0, "p_uniform": 1.0, "cramers_v_uniform": 0.0,
            "bin_counts": [],
            "expected_corrected": [],
            "expected_uniform": [],
            "peak_bin": 0, "peak_phase_fraction": 0.0,
        }

    observed = np.bincount(bin_col, minlength=k)
    expected_corrected = null_fractions * n
    expected_uniform = np.full(k, n / k)

    # Avoid zero expected counts (can cause issues with chisquare)
    # Clip expected to very small value to prevent divide-by-zero
    expected_corrected_safe = np.where(expected_corrected > 0, expected_corrected, 1e-10)

    chi2_c, p_c = scipy.stats.chisquare(observed, expected_corrected_safe)
    chi2_u, p_u = scipy.stats.chisquare(observed, expected_uniform)

    cramers_v_c = float(np.sqrt(chi2_c / (n * (k - 1))))
    cramers_v_u = float(np.sqrt(chi2_u / (n * (k - 1))))

    peak_bin = int(np.argmax(observed))
    peak_phase_fraction = float(peak_bin / k)

    return {
        "n": int(n),
        "k": int(k),
        "low_n": bool(n < 100),
        "chi2_corrected": float(chi2_c),
        "p_corrected": float(p_c),
        "cramers_v_corrected": float(cramers_v_c),
        "chi2_uniform": float(chi2_u),
        "p_uniform": float(p_u),
        "cramers_v_uniform": float(cramers_v_u),
        "bin_counts": observed.tolist(),
        "expected_corrected": expected_corrected.tolist(),
        "expected_uniform": expected_uniform.tolist(),
        "peak_bin": peak_bin,
        "peak_phase_fraction": peak_phase_fraction,
    }


# ---------------------------------------------------------------------------
# Adaptive k selection
# ---------------------------------------------------------------------------

def select_k(n: int) -> int | None:
    """
    Select adaptive bin count based on sample size.

    Returns
    -------
    24 if n >= 500; 16 if 200 <= n < 500; 12 if 100 <= n < 200;
    None if n < 100 (low-n flag).
    """
    if n >= 500:
        return 24
    if n >= 200:
        return 16
    if n >= 100:
        return 12
    return None


# ---------------------------------------------------------------------------
# Section 5: Full analysis loop
# ---------------------------------------------------------------------------

def run_stratum(
    df_stratum: pd.DataFrame,
    k: int,
    var_ranges: dict[str, tuple[float, float]],
    null_fracs: dict[str, np.ndarray],
) -> dict[str, Any]:
    """
    Run corrected chi-square for all four variables in a given stratum.

    Parameters
    ----------
    df_stratum : pd.DataFrame
        Filtered subset for this stratum.
    k : int
        Bin count (adaptive).
    var_ranges : dict
        Observed (min, max) for non-cyclic variables.
    null_fracs : dict
        Corrected null fractions for each non-cyclic variable at this k.

    Returns
    -------
    dict with keys: n, k, solar_phase, solar_declination, declination_rate,
                    earth_sun_distance
    """
    n = len(df_stratum)

    # Assign bins
    df_binned = assign_bins(df_stratum, k, var_ranges)

    # solar_phase: uniform null (cyclic)
    uniform_null = np.full(k, 1.0 / k)
    sp_result = compute_corrected_chi2(
        df_binned["bin_solar_phase"].values, k, uniform_null, n
    )
    # For solar_phase corrected == uniform; ensure identical
    sp_result["chi2_corrected"] = sp_result["chi2_uniform"]
    sp_result["p_corrected"] = sp_result["p_uniform"]
    sp_result["cramers_v_corrected"] = sp_result["cramers_v_uniform"]

    # Non-cyclic variables
    results: dict[str, Any] = {
        "n": int(n),
        "k": int(k),
        "solar_phase": sp_result,
    }

    for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
        null_f = null_fracs[var]
        bin_col = df_binned[f"bin_{var}"].values
        var_result = compute_corrected_chi2(bin_col, k, null_f, n)
        results[var] = var_result

    return results


# ---------------------------------------------------------------------------
# Section 6: Variable ranking
# ---------------------------------------------------------------------------

def compute_variable_ranking(
    full_k24: dict[str, Any],
    all_strata: dict[str, Any],
) -> dict[str, Any]:
    """
    Rank variables by corrected Cramér's V at k=24 full catalog.

    Parameters
    ----------
    full_k24 : dict
        k=24 full stratum results (four variable dicts).
    all_strata : dict
        All strata results.

    Returns
    -------
    dict with variable ranking information.
    """
    ranking = []
    for var in VARIABLES:
        stats = full_k24[var]
        ranking.append({
            "variable": var,
            "cramers_v_corrected": stats["cramers_v_corrected"],
            "p_corrected": stats["p_corrected"],
            "cramers_v_uniform": stats["cramers_v_uniform"],
            "p_uniform": stats["p_uniform"],
        })

    ranking.sort(key=lambda x: x["cramers_v_corrected"], reverse=True)

    most_sig_corrected = min(VARIABLES, key=lambda v: full_k24[v]["p_corrected"])
    most_sig_uniform = min(VARIABLES, key=lambda v: full_k24[v]["p_uniform"])

    per_stratum_top: dict[str, str] = {}
    for stratum_name, stratum_data in all_strata.items():
        if stratum_data is None:
            per_stratum_top[stratum_name] = "skipped"
            continue
        top_var = max(VARIABLES, key=lambda v: stratum_data[v]["cramers_v_corrected"])
        per_stratum_top[stratum_name] = top_var

    return {
        "by_cramers_v_corrected": ranking,
        "most_significant_variable_corrected": most_sig_corrected,
        "most_significant_variable_uniform": most_sig_uniform,
        "per_stratum_top_variable": per_stratum_top,
    }


# ---------------------------------------------------------------------------
# Section 6 continued: A1b interval alignment
# ---------------------------------------------------------------------------

def compute_a1b_alignment(
    full_k24: dict[str, Any],
    var_ranges: dict[str, tuple[float, float]],
    null_fracs_k24: dict[str, np.ndarray],
) -> dict[str, Any]:
    """
    Compute expected physical variable values at A1b interval phase centers.
    Check whether those bins are elevated in the corrected distribution.

    Parameters
    ----------
    full_k24 : dict
        k=24 full stratum results.
    var_ranges : dict
        Observed (min, max) for each non-cyclic variable.
    null_fracs_k24 : dict
        Corrected null fractions at k=24.

    Returns
    -------
    dict with interval_1, interval_2, interval_3 alignment details.
    """
    k = 24

    # A1b interval phase centers and approximate DOYs
    intervals = [
        {"name": "interval_1", "phase_center": 0.22, "doy_approx": 80},
        {"name": "interval_2", "phase_center": 0.64, "doy_approx": 234},
        {"name": "interval_3", "phase_center": 0.90, "doy_approx": 329},
    ]

    # Compute expected variable values at given DOY using the analytic model
    # Use a reference year (2000) for computing values at given DOY
    def compute_solar_vars_at_doy(doy: float) -> dict[str, float]:
        """Compute solar geometric variables at given day of year (1-based)."""
        # D = days since J2000.0 for Jan 1, 2000 + doy
        # J2000.0 = Jan 1.5, 2000 → Jan 1, 2000 is D = -0.5
        D = float(doy) - 1.0 - 0.5  # doy=1 → D = -0.5
        L = np.radians((280.46 + 0.9856474 * D) % 360)
        g = np.radians((357.528 + 0.9856003 * D) % 360)
        lam = L + np.radians(1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
        dec = float(np.degrees(np.arcsin(np.sin(np.radians(23.439)) * np.sin(lam))))
        dist = float(1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g))
        # Approximate rate via finite difference at this DOY ± 1 day
        D_p1 = D + 1.0
        D_m1 = D - 1.0
        L_p1 = np.radians((280.46 + 0.9856474 * D_p1) % 360)
        g_p1 = np.radians((357.528 + 0.9856003 * D_p1) % 360)
        lam_p1 = L_p1 + np.radians(1.915 * np.sin(g_p1) + 0.020 * np.sin(2 * g_p1))
        dec_p1 = float(np.degrees(np.arcsin(np.sin(np.radians(23.439)) * np.sin(lam_p1))))
        L_m1 = np.radians((280.46 + 0.9856474 * D_m1) % 360)
        g_m1 = np.radians((357.528 + 0.9856003 * D_m1) % 360)
        lam_m1 = L_m1 + np.radians(1.915 * np.sin(g_m1) + 0.020 * np.sin(2 * g_m1))
        dec_m1 = float(np.degrees(np.arcsin(np.sin(np.radians(23.439)) * np.sin(lam_m1))))
        rate = (dec_p1 - dec_m1) / 2.0
        return {
            "solar_declination": dec,
            "declination_rate": rate,
            "earth_sun_distance": dist,
        }

    alignment: dict[str, Any] = {}

    for ivl in intervals:
        name = ivl["name"]
        phase_center = ivl["phase_center"]
        doy = ivl["doy_approx"]

        expected_vals = compute_solar_vars_at_doy(doy)

        elevated_variables: list[str] = []

        var_details: dict[str, Any] = {
            "phase_center": phase_center,
            "doy_approx": doy,
        }

        for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
            val = expected_vals[var]
            var_details[f"{var}_expected"] = float(val)

            # Find which bin this expected value falls in
            vmin, vmax = var_ranges[var]
            bin_idx = int(np.floor((val - vmin) / (vmax - vmin) * k))
            bin_idx = np.clip(bin_idx, 0, k - 1)

            # Check if that bin is elevated in the corrected distribution
            bin_counts = np.array(full_k24[var]["bin_counts"])
            exp_corrected = np.array(full_k24[var]["expected_corrected"])

            obs = bin_counts[bin_idx]
            exp_c = exp_corrected[bin_idx]
            threshold = exp_c + np.sqrt(max(exp_c, 1.0))
            is_elevated = bool(obs > threshold)

            var_details[f"{var}_bin_idx"] = int(bin_idx)
            var_details[f"{var}_elevated"] = is_elevated

            if is_elevated:
                elevated_variables.append(var)

        var_details["elevated_variables"] = elevated_variables
        alignment[name] = var_details

    return alignment


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run A3.B5: Corrected Null-Distribution Geometric Variable Test."""

    # -----------------------------------------------------------------------
    # Section 1: Load and validate data
    # -----------------------------------------------------------------------
    logger.info("Loading solar geometry catalog...")
    solar_df = pd.read_csv(SOLAR_PATH)
    assert len(solar_df) == 9210, f"Expected 9210 rows, got {len(solar_df)}"
    logger.info(f"Solar catalog columns: {list(solar_df.columns)}")

    required_cols = ["solar_secs", "solar_declination", "declination_rate",
                     "earth_sun_distance", "depth", "latitude", "longitude", "usgs_id"]
    for col in required_cols:
        assert col in solar_df.columns, f"Missing required column: {col}"

    # Log and validate variable ranges
    confirmed_ranges = {
        "solar_declination": (-24.0, 24.0),
        "declination_rate": (-0.45, 0.45),
        "earth_sun_distance": (0.975, 1.025),
    }
    for var, (expected_min, expected_max) in confirmed_ranges.items():
        actual_min = float(solar_df[var].min())
        actual_max = float(solar_df[var].max())
        logger.info(f"{var}: actual=[{actual_min:.6f}, {actual_max:.6f}] "
                    f"confirmed=[{expected_min}, {expected_max}]")
        if actual_min < expected_min or actual_max > expected_max:
            logger.warning(f"  {var} outside confirmed bounds!")

    logger.info("Loading GSHHG classification...")
    gshhg_df = pd.read_csv(GSHHG_PATH)
    assert len(gshhg_df) == 9210, f"Expected 9210 GSHHG rows, got {len(gshhg_df)}"

    # Merge
    df = solar_df.merge(
        gshhg_df[["usgs_id", "ocean_class", "dist_to_coast_km"]],
        on="usgs_id", how="left"
    )
    assert df["dist_to_coast_km"].isna().sum() == 0, "NaN in dist_to_coast_km after merge"
    logger.info("Merge complete, no NaN in dist_to_coast_km")

    # Compute solar phase
    df["solar_phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0

    # Record observed variable ranges for bin-edge alignment
    var_ranges: dict[str, tuple[float, float]] = {}
    for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
        vmin = float(df[var].min())
        vmax = float(df[var].max())
        var_ranges[var] = (vmin, vmax)
        logger.info(f"Observed range {var}: [{vmin:.6f}, {vmax:.6f}]")

    # -----------------------------------------------------------------------
    # Section 2: Generate analytic null distributions (at k=24)
    # -----------------------------------------------------------------------
    logger.info("Generating analytic null distributions...")
    null_qc = get_null_qc_stats(1950, 2021)
    logger.info(f"Null QC: {null_qc}")

    # -----------------------------------------------------------------------
    # Section 5: Full analysis loop over strata
    # -----------------------------------------------------------------------
    strata_results: dict[str, Any] = {}

    for stratum_name, filter_expr in STRATA.items():
        if filter_expr is None:
            df_stratum = df.copy()
        else:
            # Parse compound filter expression
            if " and " in filter_expr:
                parts = filter_expr.split(" and ")
                mask = pd.Series([True] * len(df), index=df.index)
                for part in parts:
                    part = part.strip()
                    if part == "ocean_class == 'continental'":
                        mask &= df["ocean_class"] == "continental"
                    elif part == "20 <= depth < 70":
                        mask &= (df["depth"] >= 20) & (df["depth"] < 70)
                df_stratum = df[mask].copy()
            elif filter_expr == "ocean_class == 'continental'":
                df_stratum = df[df["ocean_class"] == "continental"].copy()
            elif filter_expr == "20 <= depth < 70":
                df_stratum = df[(df["depth"] >= 20) & (df["depth"] < 70)].copy()
            else:
                df_stratum = df.query(filter_expr).copy()

        n = len(df_stratum)
        k_sel = select_k(n)

        if k_sel is None:
            logger.warning(f"Stratum '{stratum_name}': n={n} < 100, skipping")
            strata_results[stratum_name] = None
            continue

        logger.info(f"Stratum '{stratum_name}': n={n}, k={k_sel}")

        # Generate null at this k using observed var ranges
        null_fracs = generate_analytic_null(1950, 2021, k_sel, var_ranges)

        stratum_data = run_stratum(df_stratum, k_sel, var_ranges, null_fracs)
        strata_results[stratum_name] = stratum_data

        logger.info(
            f"  solar_phase: V_corr={stratum_data['solar_phase']['cramers_v_corrected']:.4f}, "
            f"p_corr={stratum_data['solar_phase']['p_corrected']:.4e}"
        )
        for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
            vs = stratum_data[var]
            logger.info(
                f"  {var}: V_corr={vs['cramers_v_corrected']:.4f}, "
                f"p_corr={vs['p_corrected']:.4e}, "
                f"V_unif={vs['cramers_v_uniform']:.4f}, p_unif={vs['p_uniform']:.4e}"
            )

    # -----------------------------------------------------------------------
    # Multi-k comparison for full catalog (k=16, 24, 32)
    # -----------------------------------------------------------------------
    full_multibin: dict[str, Any] = {}

    for k_val, k_label in [(16, "k16"), (24, "k24"), (32, "k32")]:
        null_fracs_k = generate_analytic_null(1950, 2021, k_val, var_ranges)
        df_full = df.copy()
        n_full = len(df_full)

        multibin_entry: dict[str, Any] = {}

        # solar_phase uniform
        uniform_null_k = np.full(k_val, 1.0 / k_val)
        df_binned_k = assign_bins(df_full, k_val, var_ranges)
        sp_res = compute_corrected_chi2(
            df_binned_k["bin_solar_phase"].values, k_val, uniform_null_k, n_full
        )
        sp_res["chi2_corrected"] = sp_res["chi2_uniform"]
        sp_res["p_corrected"] = sp_res["p_uniform"]
        sp_res["cramers_v_corrected"] = sp_res["cramers_v_uniform"]
        multibin_entry["solar_phase"] = sp_res

        for var in ["solar_declination", "declination_rate", "earth_sun_distance"]:
            null_f = null_fracs_k[var]
            bin_col = df_binned_k[f"bin_{var}"].values
            multibin_entry[var] = compute_corrected_chi2(bin_col, k_val, null_f, n_full)

        full_multibin[k_label] = multibin_entry

    # -----------------------------------------------------------------------
    # Section 6: Variable ranking
    # -----------------------------------------------------------------------
    # Use full k=24 stratum for ranking
    full_k24_strata = strata_results.get("full") or full_multibin["k24"]

    # Build a k=24 full result to use for ranking (use full_multibin["k24"])
    full_k24_for_ranking = full_multibin["k24"]

    variable_ranking = compute_variable_ranking(full_k24_for_ranking, strata_results)

    # -----------------------------------------------------------------------
    # A1b interval alignment
    # -----------------------------------------------------------------------
    null_fracs_k24 = generate_analytic_null(1950, 2021, 24, var_ranges)
    a1b_alignment = compute_a1b_alignment(full_k24_for_ranking, var_ranges, null_fracs_k24)

    # -----------------------------------------------------------------------
    # Section 7: Build results JSON
    # -----------------------------------------------------------------------
    results: dict[str, Any] = {
        "case": "A3.B5",
        "title": "Corrected Null-Distribution Geometric Variable Test",
        "parameters": {
            "n_catalog": 9210,
            "julian_year_secs": JULIAN_YEAR_SECS,
            "midcrustal_depth_range_km": [MIDCRUSTAL_MIN_KM, MIDCRUSTAL_MAX_KM],
            "tectonic_class_boundaries_km": {
                "continental_max": 50,
                "transitional_max": 200,
            },
            "analytic_null_year_range": [1950, 2021],
            "adaptive_k_thresholds": {"k24": 500, "k16": 200, "k12": 100},
            "var_ranges": {
                var: [float(var_ranges[var][0]), float(var_ranges[var][1])]
                for var in ["solar_declination", "declination_rate", "earth_sun_distance"]
            },
        },
        "null_generation": {
            "n_synthetic_points": null_qc["n_synthetic_points"],
            "dec_synthetic_range": null_qc["dec_synthetic_range"],
            "dist_synthetic_range": null_qc["dist_synthetic_range"],
            "note": (
                "Analytic solar model: mean longitude + first-order aberration + "
                "Kepler distance; daily resolution"
            ),
        },
        "strata": strata_results,
        "full_multibin": full_multibin,
        "variable_ranking": variable_ranking,
        "a1b_alignment": a1b_alignment,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info(f"Results written to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
