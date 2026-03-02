"""Case B2: Ocean vs. Continent Location — Hydrological Loading Discrimination.

Loads ISC-GEM catalog and three ocean/continent classification files (GSHHG,
Natural Earth, PB2002), computes per-class solar-phase bin statistics, runs
classification method sensitivity comparison, performs prediction evaluation,
and writes results to output/case-b2-results.json.
"""

import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
DATA_DIR = BASE_DIR.parent / "data" / "iscgem"
PLATE_LOC_DIR = DATA_DIR / "plate-location"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("case-b2-analysis")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SOLAR_YEAR_SECS: float = 31_557_600.0  # Julian constant: 365.25 * 86400
EXPECTED_CATALOG_N: int = 9210
LOCATION_CLASSES: list[str] = ["oceanic", "transitional", "continental"]
BIN_COUNTS: list[int] = [16, 24, 32]

# A1b baseline elevated phase intervals (from case-a4-results.json sub_b)
A1B_INTERVALS: list[tuple[float, float]] = [
    (0.1875, 0.25),   # Interval 1: ~Mar 10 – Apr 1
    (0.625, 0.656),   # Interval 2: ~Aug 16 – Aug 28
    (0.875, 0.917),   # Interval 3: ~Nov 16 – Dec 1
]

# Classification methods: (key, column name)
CLASSIFICATION_METHODS: list[tuple[str, str]] = [
    ("gshhg",  "ocean_class_gshhg"),
    ("ne",     "ocean_class_ne"),
    ("pb2002", "ocean_class_pb2002"),
]


# ---------------------------------------------------------------------------
# Section 1: Data loading
# ---------------------------------------------------------------------------
def load_catalog() -> pd.DataFrame:
    """Load the ISC-GEM raw catalog and assert expected row count.

    Returns:
        DataFrame with all catalog columns.

    Raises:
        AssertionError: If row count does not equal EXPECTED_CATALOG_N.
    """
    path = DATA_DIR / "iscgem_global_6-9_1950-2021.csv"
    logger.info("Loading raw catalog from %s", path)
    df = pd.read_csv(path)
    n = len(df)
    assert n == EXPECTED_CATALOG_N, (
        f"Raw catalog: expected {EXPECTED_CATALOG_N} rows, got {n}"
    )
    logger.info("Raw catalog loaded: n=%d", n)
    return df


def load_classification_file(name: str, filename: str) -> pd.DataFrame:
    """Load a classification file and assert expected row count.

    Args:
        name: Human-readable label for logging.
        filename: CSV filename inside PLATE_LOC_DIR.

    Returns:
        DataFrame with usgs_id, ocean_class, dist_to_coast_km.

    Raises:
        AssertionError: If row count does not equal EXPECTED_CATALOG_N.
    """
    path = PLATE_LOC_DIR / filename
    logger.info("Loading %s classification from %s", name, path)
    df = pd.read_csv(path)
    n = len(df)
    assert n == EXPECTED_CATALOG_N, (
        f"{name} classification: expected {EXPECTED_CATALOG_N} rows, got {n}"
    )
    logger.info("%s classification loaded: n=%d", name, n)
    return df


def build_merged_dataframe(
    raw: pd.DataFrame,
    class_gshhg: pd.DataFrame,
    class_ne: pd.DataFrame,
    class_pb2002: pd.DataFrame,
) -> pd.DataFrame:
    """Merge all three classification files onto the raw catalog by usgs_id.

    Each classification's ocean_class is renamed with a method suffix to avoid
    column name conflicts.

    Args:
        raw: Raw ISC-GEM catalog DataFrame.
        class_gshhg: GSHHG classification DataFrame.
        class_ne: Natural Earth classification DataFrame.
        class_pb2002: PB2002 classification DataFrame.

    Returns:
        Merged DataFrame with columns ocean_class_gshhg, ocean_class_ne,
        ocean_class_pb2002 (and corresponding dist_km_* columns).

    Raises:
        AssertionError: If merge row count differs or NaN values appear.
    """
    df = raw.copy()

    for suffix, cls_df in [
        ("gshhg", class_gshhg),
        ("ne", class_ne),
        ("pb2002", class_pb2002),
    ]:
        renamed = cls_df.rename(columns={
            "ocean_class":      f"ocean_class_{suffix}",
            "dist_to_coast_km": f"dist_km_{suffix}",
        })
        df = df.merge(
            renamed[["usgs_id", f"ocean_class_{suffix}", f"dist_km_{suffix}"]],
            on="usgs_id",
            how="left",
        )
        n_after = len(df)
        assert n_after == EXPECTED_CATALOG_N, (
            f"After {suffix} merge: expected {EXPECTED_CATALOG_N} rows, got {n_after}"
        )
        nan_count = df[f"ocean_class_{suffix}"].isna().sum()
        assert nan_count == 0, (
            f"ocean_class_{suffix} has {nan_count} NaN values after merge"
        )
        logger.info(
            "Merged %s classification: n=%d, NaN in class=%d",
            suffix, n_after, nan_count,
        )

    return df


def compute_phase(df: pd.DataFrame) -> pd.DataFrame:
    """Add solar phase column to merged dataframe using Julian year constant.

    Phase is computed as (solar_secs / SOLAR_YEAR_SECS) % 1.0, mapping each
    event to [0, 1) of the annual cycle. Consistent with prior cases A4, B1.

    Args:
        df: DataFrame with solar_secs column.

    Returns:
        DataFrame with added solar_phase column.
    """
    df = df.copy()
    df["solar_phase"] = (df["solar_secs"] / SOLAR_YEAR_SECS) % 1.0
    logger.info(
        "Solar phase computed: min=%.6f max=%.6f",
        df["solar_phase"].min(), df["solar_phase"].max(),
    )
    return df


# ---------------------------------------------------------------------------
# Section 2: Per-class bin statistics
# ---------------------------------------------------------------------------
def compute_bin_stats(
    subset: pd.DataFrame,
    class_label: str,
    method: str,
    k: int,
) -> dict[str, Any]:
    """Compute chi-square, Rayleigh, Cramér's V, and elevated-bin stats.

    Args:
        subset: DataFrame filtered to one location class.
        class_label: e.g. "oceanic", "transitional", "continental".
        method: Classification method name for logging.
        k: Number of phase bins.

    Returns:
        Dict with n, k, chi2, p_chi2, cramer_v, rayleigh_R, p_rayleigh,
        mean_phase_fraction, bin_counts, elevated_bins, elevated_intervals.
    """
    n = len(subset)
    if n == 0:
        raise ValueError(
            f"Empty subset for method={method}, class={class_label}, k={k}"
        )
    if n < 200:
        logger.warning(
            "Small sample: method=%s class=%s k=%d n=%d (< 200)",
            method, class_label, k, n,
        )

    phase = subset["solar_phase"].values

    # Bin assignment using phase-normalized approach
    bin_idx = np.floor(phase * k).astype(int)
    bin_idx = np.clip(bin_idx, 0, k - 1)

    obs = np.bincount(bin_idx, minlength=k).astype(float)
    expected = np.full(k, n / k)

    # Chi-square
    chi2_stat, p_chi2 = scipy.stats.chisquare(obs, expected)

    # Rayleigh statistic
    angles = 2.0 * np.pi * phase
    R = float(np.abs(np.mean(np.exp(1j * angles))))
    p_rayleigh = float(np.exp(-n * R ** 2))

    # Cramér's V
    cramer_v = float(np.sqrt(chi2_stat / (n * (k - 1))))

    # Mean phase angle → fraction of year
    mean_sin = float(np.mean(np.sin(angles)))
    mean_cos = float(np.mean(np.cos(angles)))
    mean_angle_rad = float(np.arctan2(mean_sin, mean_cos))
    mean_phase_fraction = float((mean_angle_rad / (2.0 * np.pi)) % 1.0)

    # Elevated bins at 1-SD threshold
    threshold = expected[0] + np.sqrt(expected[0])
    elevated_bin_indices = [int(i) for i in np.where(obs > threshold)[0]]

    # Merge adjacent elevated bins into intervals
    elevated_intervals = _bins_to_intervals(elevated_bin_indices, k)

    result: dict[str, Any] = {
        "n": int(n),
        "k": int(k),
        "chi2": float(chi2_stat),
        "p_chi2": float(p_chi2),
        "cramer_v": cramer_v,
        "rayleigh_R": R,
        "p_rayleigh": p_rayleigh,
        "mean_phase_fraction": mean_phase_fraction,
        "bin_counts": obs.astype(int).tolist(),
        "expected_count": float(expected[0]),
        "threshold_1sd": float(threshold),
        "elevated_bin_indices": elevated_bin_indices,
        "elevated_intervals": elevated_intervals,
    }
    logger.info(
        "%s/%s k=%d: n=%d chi2=%.4f p=%.4e V=%.4f R=%.4f",
        method, class_label, k, n, chi2_stat, p_chi2, cramer_v, R,
    )
    return result


def _bins_to_intervals(
    elevated_indices: list[int],
    k: int,
) -> list[dict[str, Any]]:
    """Merge adjacent elevated bin indices into phase intervals.

    Args:
        elevated_indices: Sorted list of elevated bin indices.
        k: Total number of bins.

    Returns:
        List of dicts with phase_start, phase_end for each contiguous run.
    """
    if not elevated_indices:
        return []

    intervals: list[dict[str, Any]] = []
    run_start = elevated_indices[0]
    run_end = elevated_indices[0]

    for idx in elevated_indices[1:]:
        if idx == run_end + 1:
            run_end = idx
        else:
            intervals.append({
                "phase_start": round(run_start / k, 6),
                "phase_end": round((run_end + 1) / k, 6),
            })
            run_start = idx
            run_end = idx

    intervals.append({
        "phase_start": round(run_start / k, 6),
        "phase_end": round((run_end + 1) / k, 6),
    })
    return intervals


def compute_all_class_stats(df: pd.DataFrame) -> dict[str, Any]:
    """Compute per-class bin statistics for all methods and bin counts.

    Args:
        df: Merged DataFrame with solar_phase and all ocean_class_* columns.

    Returns:
        Nested dict: class_stats[method][class_label][k_key] = stat dict.
    """
    class_stats: dict[str, Any] = {}

    for method_key, col_name in CLASSIFICATION_METHODS:
        class_stats[method_key] = {}

        # Validate labels
        unique_labels = set(df[col_name].unique())
        invalid = unique_labels - set(LOCATION_CLASSES)
        assert not invalid, (
            f"Unexpected labels in {col_name}: {invalid}"
        )

        for class_label in LOCATION_CLASSES:
            subset = df[df[col_name] == class_label]
            n_class = len(subset)
            logger.info(
                "Method=%s class=%s n=%d", method_key, class_label, n_class
            )

            class_stats[method_key][class_label] = {"n": int(n_class)}
            for k in BIN_COUNTS:
                k_key = f"k{k}"
                class_stats[method_key][class_label][k_key] = compute_bin_stats(
                    subset, class_label, method_key, k
                )

    return class_stats


# ---------------------------------------------------------------------------
# Section 3: Classification method sensitivity comparison
# ---------------------------------------------------------------------------
def compute_method_sensitivity(
    class_stats: dict[str, Any],
) -> list[dict[str, Any]]:
    """Compare classification method results at k=24 across all classes.

    For each class, tabulates chi2, p_chi2, Cramér's V for each method and
    flags whether the three methods agree on statistical significance.

    Args:
        class_stats: Nested dict from compute_all_class_stats.

    Returns:
        List of dicts: one per class with per-method stats and agreement flag.
    """
    sensitivity: list[dict[str, Any]] = []

    for class_label in LOCATION_CLASSES:
        entry: dict[str, Any] = {"class": class_label}
        p_values: list[float] = []
        cramer_vs: list[float] = []

        for method_key, _ in CLASSIFICATION_METHODS:
            k24 = class_stats[method_key][class_label]["k24"]
            entry[f"{method_key}_chi2"] = k24["chi2"]
            entry[f"{method_key}_p"] = k24["p_chi2"]
            entry[f"{method_key}_cramer_v"] = k24["cramer_v"]
            p_values.append(k24["p_chi2"])
            cramer_vs.append(k24["cramer_v"])

        # Agreement: all methods agree on significance or all agree on non-significance
        sig_flags = [p < 0.05 for p in p_values]
        agreement = (all(sig_flags) or not any(sig_flags))
        entry["agreement"] = bool(agreement)
        entry["avg_cramer_v"] = float(np.mean(cramer_vs))

        logger.info(
            "Sensitivity class=%s p=[%.4e, %.4e, %.4e] agreement=%s",
            class_label, p_values[0], p_values[1], p_values[2], agreement,
        )
        sensitivity.append(entry)

    return sensitivity


# ---------------------------------------------------------------------------
# Section 4: Prediction evaluation
# ---------------------------------------------------------------------------
def _intervals_overlap(
    obs_intervals: list[dict[str, Any]],
    ref_intervals: list[tuple[float, float]],
    overlap_threshold: float = 0.5,
) -> bool:
    """Check if any observed elevated interval overlaps a reference interval.

    Args:
        obs_intervals: Elevated intervals from compute_bin_stats output.
        ref_intervals: Reference A1b intervals as (phase_start, phase_end) tuples.
        overlap_threshold: Fraction of reference interval width that must overlap.

    Returns:
        True if at least one reference interval is substantially overlapped.
    """
    for ref_start, ref_end in ref_intervals:
        ref_width = ref_end - ref_start
        if ref_width <= 0:
            continue
        for obs in obs_intervals:
            overlap = max(
                0.0,
                min(ref_end, obs["phase_end"]) - max(ref_start, obs["phase_start"])
            )
            if overlap / ref_width > overlap_threshold:
                return True
    return False


def evaluate_predictions(
    class_stats: dict[str, Any],
) -> dict[str, Any]:
    """Evaluate geometric vs. hydrological predictions using GSHHG at k=24.

    Args:
        class_stats: Nested dict from compute_all_class_stats.

    Returns:
        Dict with prediction evaluation flags and primary_conclusion.
    """
    gshhg_oceanic = class_stats["gshhg"]["oceanic"]["k24"]
    gshhg_continental = class_stats["gshhg"]["continental"]["k24"]

    oceanic_significant = bool(gshhg_oceanic["p_chi2"] < 0.05)
    continental_significant = bool(gshhg_continental["p_chi2"] < 0.05)

    # Hydrological: continental Cramér's V > 2× oceanic, OR oceanic not sig
    continental_v = gshhg_continental["cramer_v"]
    oceanic_v = gshhg_oceanic["cramer_v"]
    continental_stronger = bool(
        continental_v > 2.0 * oceanic_v or not oceanic_significant
    )

    # Geometric: signal appears in both oceanic and continental
    geometric_supported = bool(oceanic_significant and continental_significant)
    hydrological_supported = bool(continental_stronger)

    # Check if oceanic elevated intervals overlap A1b baseline intervals
    oceanic_intervals = gshhg_oceanic.get("elevated_intervals", [])
    oceanic_matches_a1b = _intervals_overlap(oceanic_intervals, A1B_INTERVALS)

    # Primary conclusion
    if geometric_supported and not hydrological_supported:
        primary_conclusion = "geometric"
    elif hydrological_supported and not geometric_supported:
        primary_conclusion = "hydrological"
    elif not oceanic_significant and not continental_significant:
        primary_conclusion = "ambiguous"
    else:
        primary_conclusion = "ambiguous"

    result: dict[str, Any] = {
        "primary_method": "gshhg",
        "oceanic_p_chi2": float(gshhg_oceanic["p_chi2"]),
        "oceanic_cramer_v": float(oceanic_v),
        "continental_p_chi2": float(gshhg_continental["p_chi2"]),
        "continental_cramer_v": float(continental_v),
        "oceanic_significant": oceanic_significant,
        "continental_significant": continental_significant,
        "continental_cramer_v_ratio": float(continental_v / oceanic_v) if oceanic_v > 0 else None,
        "continental_stronger": continental_stronger,
        "oceanic_matches_a1b_intervals": bool(oceanic_matches_a1b),
        "geometric_supported": geometric_supported,
        "hydrological_supported": hydrological_supported,
        "primary_conclusion": primary_conclusion,
    }

    logger.info(
        "Prediction evaluation: oceanic_sig=%s continental_sig=%s "
        "continental_stronger=%s primary_conclusion=%s",
        oceanic_significant, continental_significant,
        continental_stronger, primary_conclusion,
    )
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Run the full Case B2 analysis and write results JSON."""
    logger.info("=== Case B2: Ocean vs. Continent Analysis ===")
    logger.info("Julian year constant: %.1f s", SOLAR_YEAR_SECS)

    # --- Section 1: Load data ---
    raw = load_catalog()
    class_gshhg = load_classification_file(
        "GSHHG", "ocean_class_gshhg_global.csv"
    )
    class_ne = load_classification_file(
        "Natural Earth", "ocean_class_ne_global.csv"
    )
    class_pb2002 = load_classification_file(
        "PB2002", "ocean_class_pb2002_global.csv"
    )

    # Merge and compute phase
    df = build_merged_dataframe(raw, class_gshhg, class_ne, class_pb2002)
    df = compute_phase(df)

    assert df["solar_phase"].between(0.0, 1.0, inclusive="left").all(), (
        "Some solar_phase values outside [0, 1)"
    )

    # --- Section 2: Per-class bin statistics ---
    logger.info("Computing per-class bin statistics for all methods...")
    class_stats = compute_all_class_stats(df)

    # --- Section 3: Method sensitivity comparison ---
    logger.info("Computing classification method sensitivity comparison...")
    method_sensitivity = compute_method_sensitivity(class_stats)

    # --- Section 4: Prediction evaluation ---
    logger.info("Evaluating geometric vs. hydrological predictions...")
    prediction_evaluation = evaluate_predictions(class_stats)

    # --- Write results JSON ---
    results: dict[str, Any] = {
        "case": "B2",
        "title": "Ocean vs. Continent Location — Hydrological Loading Discrimination",
        "solar_year_secs": SOLAR_YEAR_SECS,
        "catalog_n": EXPECTED_CATALOG_N,
        "location_classes": LOCATION_CLASSES,
        "classification_methods": [m for m, _ in CLASSIFICATION_METHODS],
        "a1b_intervals": [
            {"phase_start": s, "phase_end": e} for s, e in A1B_INTERVALS
        ],
        "class_stats": class_stats,
        "method_sensitivity": method_sensitivity,
        "prediction_evaluation": prediction_evaluation,
    }

    out_path = OUTPUT_DIR / "case-b2-results.json"
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", out_path)

    # Summary log
    pred = prediction_evaluation
    logger.info(
        "SUMMARY: oceanic sig=%s (p=%.4e V=%.4f) | "
        "continental sig=%s (p=%.4e V=%.4f) | conclusion=%s",
        pred["oceanic_significant"],
        pred["oceanic_p_chi2"],
        pred["oceanic_cramer_v"],
        pred["continental_significant"],
        pred["continental_p_chi2"],
        pred["continental_cramer_v"],
        pred["primary_conclusion"],
    )


if __name__ == "__main__":
    main()
