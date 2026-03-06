"""
Case A3.B2: Hemisphere Stratification Refinement

Four sub-tests to determine whether the NH/SH solar-phase signal asymmetry
observed in A2.B1 is a tectonic-composition artifact or a genuine hemispheric
phase difference:
1. Tectonic-matched hemisphere comparison (18 cells: 3 catalogs × 3 classes × 2 hemispheres)
2. Mid-crustal hemisphere split (locked to 20-70 km depth band from A3.B4)
3. Phase alignment comparison (wrapped NH/SH peak-phase offset)
4. Revised Interval 1 SH threshold sensitivity (all-SH vs. continental-SH)
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

RAW_PATH   = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
GK_PATH    = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv"
REAS_PATH  = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv"

OUTPUT_PATH = BASE_DIR / "output" / "case-a3-b2-results.json"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
JULIAN_YEAR_SECS: float = 31_557_600.0
ADAPTIVE_K_THRESHOLDS: dict[str, int] = {"k24": 500, "k16": 200, "k12": 100}
TECTONIC_CLASSES: list[str] = ["continental", "transitional", "oceanic"]
MIDCRUSTAL_DEPTH_MIN: float = 20.0
MIDCRUSTAL_DEPTH_MAX: float = 70.0
INTERVAL_BINS: dict[str, list[int]] = {
    "interval_1": [4, 5],
    "interval_2": [15],
    "interval_3": [21],
}
INTERVAL1_THRESHOLDS: list[float] = [0.33, 0.40, 0.45, 0.50]


# ---------------------------------------------------------------------------
# Section 2: Core statistics helper
# ---------------------------------------------------------------------------

def adaptive_k(n: int) -> int | None:
    """
    Return bin count per adaptive-k rule; None if n < 100.

    Parameters
    ----------
    n : int
        Sample size.

    Returns
    -------
    int | None
        24 if n >= 500; 16 if 200 <= n < 500; 12 if 100 <= n < 200; None if n < 100.
    """
    if n >= 500:
        return 24
    elif n >= 200:
        return 16
    elif n >= 100:
        return 12
    return None


def compute_chi2_stats(phases: np.ndarray) -> dict[str, Any]:
    """
    Compute chi-square uniformity test on solar phase distribution.

    Parameters
    ----------
    phases : np.ndarray
        Array of phase values in [0, 1).

    Returns
    -------
    dict with keys: n, k, low_n, chi2, p_chi2, cramers_v, bin_counts,
                    peak_bin, peak_phase, interval_1_z, interval_2_z, interval_3_z.
    """
    n = len(phases)
    k = adaptive_k(n)
    if k is None:
        return {
            "n": n,
            "k": None,
            "low_n": True,
            "chi2": None,
            "p_chi2": None,
            "cramers_v": None,
            "bin_counts": None,
            "peak_bin": None,
            "peak_phase": None,
            "interval_1_z": None,
            "interval_2_z": None,
            "interval_3_z": None,
        }

    bin_counts = np.bincount(
        np.floor(phases * k).astype(int) % k, minlength=k
    )
    expected = n / k
    chi2_stat, p_chi2 = scipy.stats.chisquare(bin_counts, np.full(k, expected))
    cramers_v = np.sqrt(chi2_stat / (n * (k - 1)))
    peak_bin = int(np.argmax(bin_counts))
    peak_phase = float((peak_bin + 0.5) / k)

    interval_z: dict[str, float | None] = {}
    for iname, bins in INTERVAL_BINS.items():
        valid_bins = [b for b in bins if b < k]
        if not valid_bins or expected == 0:
            interval_z[iname] = None
        else:
            obs_i = sum(int(bin_counts[b]) for b in valid_bins)
            exp_i = len(valid_bins) * expected
            interval_z[iname] = float((obs_i - exp_i) / np.sqrt(exp_i))

    return {
        "n": int(n),
        "k": int(k),
        "low_n": False,
        "chi2": float(chi2_stat),
        "p_chi2": float(p_chi2),
        "cramers_v": float(cramers_v),
        "bin_counts": bin_counts.tolist(),
        "peak_bin": int(peak_bin),
        "peak_phase": peak_phase,
        "interval_1_z": interval_z["interval_1"],
        "interval_2_z": interval_z["interval_2"],
        "interval_3_z": interval_z["interval_3"],
    }


# ---------------------------------------------------------------------------
# Section 1: Hemisphere split helper
# ---------------------------------------------------------------------------

def split_hemisphere(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """
    Return (df_nh, df_sh, n_equatorial).

    Parameters
    ----------
    df : pd.DataFrame
        Event catalog with 'latitude' column.

    Returns
    -------
    tuple of (Northern Hemisphere df, Southern Hemisphere df, equatorial count)
    """
    n_equatorial = int((df["latitude"] == 0).sum())
    return (
        df[df["latitude"] > 0].copy(),
        df[df["latitude"] < 0].copy(),
        n_equatorial,
    )


# ---------------------------------------------------------------------------
# Section 3: Sub-test 1 — Tectonic-matched hemisphere comparison
# ---------------------------------------------------------------------------

def run_tectonic_hemisphere_comparison(
    catalogs: dict[str, pd.DataFrame],
    tectonic_classes: list[str],
) -> dict[str, Any]:
    """
    Compute chi-square stats for each catalog × tectonic class × hemisphere cell.

    Parameters
    ----------
    catalogs : dict mapping catalog name to DataFrame (with 'ocean_class', 'latitude', 'phase')
    tectonic_classes : list of tectonic class labels

    Returns
    -------
    Nested dict: catalog → tectonic_class → hemisphere → stats dict
    """
    results: dict[str, Any] = {}

    for cat_name, df in catalogs.items():
        results[cat_name] = {}
        for tclass in tectonic_classes:
            tclass_df = df[df["ocean_class"] == tclass]
            nh_df = tclass_df[tclass_df["latitude"] > 0]
            sh_df = tclass_df[tclass_df["latitude"] < 0]

            nh_stats = compute_chi2_stats(nh_df["phase"].values)
            sh_stats = compute_chi2_stats(sh_df["phase"].values)

            nh_sig = (nh_stats["p_chi2"] is not None and nh_stats["p_chi2"] < 0.05)
            sh_sig = (sh_stats["p_chi2"] is not None and sh_stats["p_chi2"] < 0.05)

            logger.info(
                f"Sub-test 1 | {cat_name} | {tclass} | NH: n={nh_stats['n']}, "
                f"p={nh_stats['p_chi2']}, sig={nh_sig} | "
                f"SH: n={sh_stats['n']}, p={sh_stats['p_chi2']}, sig={sh_sig}"
            )

            results[cat_name][tclass] = {
                "nh": nh_stats,
                "sh": sh_stats,
                "nh_significant": nh_sig,
                "sh_significant": sh_sig,
                "both_significant": nh_sig and sh_sig,
                "neither_significant": not nh_sig and not sh_sig,
            }

    return results


# ---------------------------------------------------------------------------
# Section 4: Sub-test 2 — Mid-crustal hemisphere split
# ---------------------------------------------------------------------------

def run_midcrustal_hemisphere_split(catalogs: dict[str, pd.DataFrame]) -> dict[str, Any]:
    """
    Compute chi-square stats for mid-crustal band (20-70 km) by hemisphere.

    Parameters
    ----------
    catalogs : dict mapping catalog name to DataFrame (with 'depth', 'latitude', 'phase')

    Returns
    -------
    Nested dict: catalog → global/nh/sh → stats dict
    """
    results: dict[str, Any] = {}

    for cat_name, df in catalogs.items():
        mid_df = df[
            (df["depth"] >= MIDCRUSTAL_DEPTH_MIN) & (df["depth"] < MIDCRUSTAL_DEPTH_MAX)
        ].copy()

        global_stats = compute_chi2_stats(mid_df["phase"].values)

        nh_df = mid_df[mid_df["latitude"] > 0]
        sh_df = mid_df[mid_df["latitude"] < 0]

        nh_stats = compute_chi2_stats(nh_df["phase"].values)
        sh_stats = compute_chi2_stats(sh_df["phase"].values)

        logger.info(
            f"Sub-test 2 | {cat_name} | global: n={global_stats['n']}, "
            f"chi2={global_stats['chi2']}, p={global_stats['p_chi2']}"
        )
        logger.info(
            f"Sub-test 2 | {cat_name} | NH: n={nh_stats['n']}, "
            f"p={nh_stats['p_chi2']} | SH: n={sh_stats['n']}, p={sh_stats['p_chi2']}"
        )

        results[cat_name] = {
            "global": global_stats,
            "nh": nh_stats,
            "sh": sh_stats,
        }

    # Regression anchor check for full catalog
    full_global_chi2 = results.get("full", {}).get("global", {}).get("chi2")
    if full_global_chi2 is not None:
        anchor_target = 85.48
        anchor_tol = 2.0
        deviation = abs(full_global_chi2 - anchor_target)
        logger.info(
            f"Regression anchor: full mid-crustal chi2={full_global_chi2:.4f}, "
            f"target={anchor_target} ± {anchor_tol}, deviation={deviation:.4f}, "
            f"within_tolerance={deviation <= anchor_tol}"
        )

    return results


# ---------------------------------------------------------------------------
# Section 5: Sub-test 3 — Phase alignment comparison
# ---------------------------------------------------------------------------

def run_phase_alignment(
    subtest1: dict[str, Any],
    subtest2: dict[str, Any],
) -> dict[str, Any]:
    """
    Collect significant NH/SH pairs and classify phase alignment.

    Parameters
    ----------
    subtest1 : sub-test 1 results (catalog → tectonic_class → hemisphere → stats)
    subtest2 : sub-test 2 results (catalog → global/nh/sh → stats)

    Returns
    -------
    Summary dict with pair list and alignment counts.
    """
    pairs: list[dict[str, Any]] = []

    # Pairs from sub-test 1: catalog × tectonic_class
    for cat_name, tclass_dict in subtest1.items():
        for tclass, hemi_dict in tclass_dict.items():
            nh_stats = hemi_dict["nh"]
            sh_stats = hemi_dict["sh"]
            nh_sig = hemi_dict["nh_significant"]
            sh_sig = hemi_dict["sh_significant"]

            # Include pair if at least one is significant and has valid peak_phase
            if (nh_sig or sh_sig) and nh_stats.get("peak_phase") is not None and sh_stats.get("peak_phase") is not None:
                nh_peak = nh_stats["peak_phase"]
                sh_peak = sh_stats["peak_phase"]
                delta = (nh_peak - sh_peak + 0.5) % 1.0 - 0.5

                abs_delta = abs(delta)
                if abs_delta < 0.083:
                    alignment = "in_phase"
                elif abs_delta > 0.417:
                    alignment = "anti_phase"
                else:
                    alignment = "offset"

                pairs.append({
                    "source": f"subtest1_{tclass}",
                    "catalog": cat_name,
                    "nh_peak_phase": float(nh_peak),
                    "sh_peak_phase": float(sh_peak),
                    "delta_phase": float(delta),
                    "alignment": alignment,
                    "nh_significant": nh_sig,
                    "sh_significant": sh_sig,
                })

    # Pairs from sub-test 2: catalog × mid-crustal
    for cat_name, hemi_dict in subtest2.items():
        nh_stats = hemi_dict.get("nh", {})
        sh_stats = hemi_dict.get("sh", {})
        nh_sig = (nh_stats.get("p_chi2") is not None and nh_stats["p_chi2"] < 0.05)
        sh_sig = (sh_stats.get("p_chi2") is not None and sh_stats["p_chi2"] < 0.05)

        if (nh_sig or sh_sig) and nh_stats.get("peak_phase") is not None and sh_stats.get("peak_phase") is not None:
            nh_peak = nh_stats["peak_phase"]
            sh_peak = sh_stats["peak_phase"]
            delta = (nh_peak - sh_peak + 0.5) % 1.0 - 0.5

            abs_delta = abs(delta)
            if abs_delta < 0.083:
                alignment = "in_phase"
            elif abs_delta > 0.417:
                alignment = "anti_phase"
            else:
                alignment = "offset"

            pairs.append({
                "source": "subtest2_midcrustal",
                "catalog": cat_name,
                "nh_peak_phase": float(nh_peak),
                "sh_peak_phase": float(sh_peak),
                "delta_phase": float(delta),
                "alignment": alignment,
                "nh_significant": nh_sig,
                "sh_significant": sh_sig,
            })

    n_in_phase = sum(1 for p in pairs if p["alignment"] == "in_phase")
    n_anti_phase = sum(1 for p in pairs if p["alignment"] == "anti_phase")
    n_offset = sum(1 for p in pairs if p["alignment"] == "offset")

    # Determine dominant alignment
    counts = {"in_phase": n_in_phase, "anti_phase": n_anti_phase, "offset": n_offset}
    max_count = max(counts.values()) if counts else 0
    dominant_candidates = [k for k, v in counts.items() if v == max_count]
    dominant_alignment = dominant_candidates[0] if len(dominant_candidates) == 1 else "mixed"

    logger.info(
        f"Phase alignment: n_pairs={len(pairs)}, in_phase={n_in_phase}, "
        f"anti_phase={n_anti_phase}, offset={n_offset}, dominant={dominant_alignment}"
    )

    return {
        "n_pairs_evaluated": len(pairs),
        "n_in_phase": n_in_phase,
        "n_anti_phase": n_anti_phase,
        "n_offset": n_offset,
        "dominant_alignment": dominant_alignment,
        "pairs": pairs,
    }


# ---------------------------------------------------------------------------
# Section 6: Sub-test 4 — Revised Interval 1 SH threshold sensitivity
# ---------------------------------------------------------------------------

def run_interval1_threshold_sensitivity(df_full: pd.DataFrame) -> dict[str, Any]:
    """
    Test whether A2.B1 Interval 1 SH absence persists under tectonic composition control.

    Parameters
    ----------
    df_full : pd.DataFrame
        Full catalog with 'latitude', 'phase', 'ocean_class' columns.

    Returns
    -------
    dict with sh_all and sh_continental threshold records and flip point info.
    """
    populations: dict[str, pd.DataFrame] = {
        "sh_all": df_full[df_full["latitude"] < 0].copy(),
        "sh_continental": df_full[
            (df_full["latitude"] < 0) & (df_full["ocean_class"] == "continental")
        ].copy(),
    }

    results: dict[str, Any] = {}
    flip_thresholds: dict[str, float | None] = {}

    for pop_name, pop_df in populations.items():
        n_sh = len(pop_df)
        phases = pop_df["phase"].values

        # Compute k=24 bin counts (use k=24 regardless of n for consistency)
        k = 24
        if n_sh == 0:
            bin_counts = np.zeros(k, dtype=int)
        else:
            bin_counts = np.bincount(
                np.floor(phases * k).astype(int) % k, minlength=k
            )

        expected = n_sh / k if n_sh > 0 else 0.0

        # Interval 1 bins at k=24: bins 4 and 5
        interval1_bins = [4, 5]
        elevated_bins = [
            b for b in interval1_bins
            if n_sh > 0 and bin_counts[b] > expected + np.sqrt(expected)
        ]
        overlap_fraction = len(elevated_bins) / 2.0

        threshold_records: list[dict[str, Any]] = []
        prev_classification: str | None = None
        flip_threshold: float | None = None

        for t in INTERVAL1_THRESHOLDS:
            classification = "present" if overlap_fraction >= t else "absent"
            threshold_records.append({
                "threshold": float(t),
                "n_sh": int(n_sh),
                "overlap_fraction": float(overlap_fraction),
                "interval_1_classification": classification,
                "bin_4_obs": int(bin_counts[4]) if n_sh > 0 else 0,
                "bin_5_obs": int(bin_counts[5]) if n_sh > 0 else 0,
                "expected": float(expected),
            })

            if prev_classification is not None and classification != prev_classification:
                if flip_threshold is None:
                    flip_threshold = float(t)

            prev_classification = classification

        results[pop_name] = threshold_records
        flip_thresholds[pop_name] = flip_threshold

        logger.info(
            f"Sub-test 4 | {pop_name}: n={n_sh}, overlap_fraction={overlap_fraction:.3f}, "
            f"elevated_bins={elevated_bins}, flip_threshold={flip_threshold}"
        )

    return {
        "sh_all": results["sh_all"],
        "sh_continental": results["sh_continental"],
        "sh_all_flip_threshold": flip_thresholds["sh_all"],
        "sh_continental_flip_threshold": flip_thresholds["sh_continental"],
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run A3.B2: Hemisphere Stratification Refinement."""

    # -----------------------------------------------------------------------
    # Section 1: Load and merge data
    # -----------------------------------------------------------------------
    logger.info("Loading raw catalog...")
    df_full = pd.read_csv(RAW_PATH, parse_dates=["event_at"])
    assert len(df_full) == 9210, f"Expected 9210 rows, got {len(df_full)}"
    logger.info(f"Raw catalog: {len(df_full)} rows")

    logger.info("Loading GSHHG classification...")
    df_gshhg = pd.read_csv(GSHHG_PATH)
    assert len(df_gshhg) == 9210, f"Expected 9210 GSHHG rows, got {len(df_gshhg)}"
    assert df_gshhg["dist_to_coast_km"].isna().sum() == 0, "NaN in dist_to_coast_km"

    logger.info("Loading G-K mainshocks...")
    df_gk = pd.read_csv(GK_PATH, parse_dates=["event_at"])
    assert len(df_gk) == 5883, f"Expected 5883 GK rows, got {len(df_gk)}"
    logger.info(f"G-K mainshocks: {len(df_gk)} rows")

    logger.info("Loading Reasenberg mainshocks...")
    df_reas = pd.read_csv(REAS_PATH, parse_dates=["event_at"])
    assert len(df_reas) == 8265, f"Expected 8265 Reasenberg rows, got {len(df_reas)}"
    logger.info(f"Reasenberg mainshocks: {len(df_reas)} rows")

    # Merge GSHHG classification onto all catalogs
    gshhg_cols = df_gshhg[["usgs_id", "ocean_class", "dist_to_coast_km"]]

    for name, df in [("full", df_full), ("gk", df_gk), ("reas", df_reas)]:
        merged = df.merge(gshhg_cols, on="usgs_id", how="left")
        n_nan = merged["dist_to_coast_km"].isna().sum()
        assert n_nan == 0, f"NaN in dist_to_coast_km after merge for {name}: {n_nan}"
        if name == "full":
            df_full = merged
        elif name == "gk":
            df_gk = merged
        else:
            df_reas = merged

    # Compute phase for all catalogs
    for df in [df_full, df_gk, df_reas]:
        df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0

    # Verify phase ranges
    for name, df in [("full", df_full), ("gk", df_gk), ("reas", df_reas)]:
        assert df["phase"].between(0.0, 1.0, inclusive="left").all(), \
            f"Phase out of [0, 1) range in {name}"
        logger.info(
            f"Phase computed for {name}: min={df['phase'].min():.4f}, "
            f"max={df['phase'].max():.4f}"
        )

    # -----------------------------------------------------------------------
    # Hemisphere split for all catalogs
    # -----------------------------------------------------------------------
    nh_full, sh_full, eq_full = split_hemisphere(df_full)
    nh_gk, sh_gk, eq_gk = split_hemisphere(df_gk)
    nh_reas, sh_reas, eq_reas = split_hemisphere(df_reas)

    catalog_sizes: dict[str, dict[str, int]] = {
        "full": {"n_nh": len(nh_full), "n_sh": len(sh_full), "n_equatorial": eq_full},
        "gk":   {"n_nh": len(nh_gk),   "n_sh": len(sh_gk),   "n_equatorial": eq_gk},
        "reas": {"n_nh": len(nh_reas),  "n_sh": len(sh_reas),  "n_equatorial": eq_reas},
    }
    logger.info(f"Catalog sizes: {catalog_sizes}")

    # Verify tectonic partition sums
    for name, df in [("full", df_full), ("gk", df_gk), ("reas", df_reas)]:
        tectonic_counts = df["ocean_class"].value_counts()
        total_tectonic = int(tectonic_counts.get("continental", 0) +
                             tectonic_counts.get("transitional", 0) +
                             tectonic_counts.get("oceanic", 0))
        logger.info(f"Tectonic partition {name}: {dict(tectonic_counts)}, total={total_tectonic}")

    # -----------------------------------------------------------------------
    # Sub-test 1: Tectonic-matched hemisphere comparison
    # -----------------------------------------------------------------------
    logger.info("Running Sub-test 1: Tectonic-matched hemisphere comparison...")
    CATALOGS: dict[str, pd.DataFrame] = {"full": df_full, "gk": df_gk, "reas": df_reas}

    subtest1 = run_tectonic_hemisphere_comparison(CATALOGS, TECTONIC_CLASSES)

    # -----------------------------------------------------------------------
    # Sub-test 2: Mid-crustal hemisphere split
    # -----------------------------------------------------------------------
    logger.info("Running Sub-test 2: Mid-crustal hemisphere split...")
    subtest2 = run_midcrustal_hemisphere_split(CATALOGS)

    # -----------------------------------------------------------------------
    # Sub-test 3: Phase alignment comparison
    # -----------------------------------------------------------------------
    logger.info("Running Sub-test 3: Phase alignment comparison...")
    subtest3 = run_phase_alignment(subtest1, subtest2)

    # -----------------------------------------------------------------------
    # Sub-test 4: Interval 1 SH threshold sensitivity
    # -----------------------------------------------------------------------
    logger.info("Running Sub-test 4: Interval 1 SH threshold sensitivity...")
    subtest4 = run_interval1_threshold_sensitivity(df_full)

    # -----------------------------------------------------------------------
    # Section 7: Write results JSON
    # -----------------------------------------------------------------------
    results: dict[str, Any] = {
        "case": "A3.B2",
        "title": "Hemisphere Stratification Refinement",
        "parameters": {
            "n_catalog": 9210,
            "n_gk": 5883,
            "n_reas": 8265,
            "julian_year_secs": JULIAN_YEAR_SECS,
            "adaptive_k_thresholds": ADAPTIVE_K_THRESHOLDS,
            "tectonic_class_boundaries_km": {
                "continental_max": 50,
                "transitional_max": 200,
            },
            "midcrustal_depth_range_km": [MIDCRUSTAL_DEPTH_MIN, MIDCRUSTAL_DEPTH_MAX],
            "interval1_threshold_sweep": INTERVAL1_THRESHOLDS,
            "phase_alignment_bin_tolerance": 2,
        },
        "catalog_sizes": catalog_sizes,
        "subtest_1_tectonic_hemisphere": subtest1,
        "subtest_2_midcrustal_hemisphere": subtest2,
        "subtest_3_phase_alignment": subtest3,
        "subtest_4_interval1_threshold": subtest4,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info(f"Results written to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
