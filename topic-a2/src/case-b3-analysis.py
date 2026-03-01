"""Case B3: Tectonic Regime Stratification.

Loads ISC-GEM catalog and GCMT focal mechanism join file, classifies events by
tectonic mechanism (thrust, normal, strike-slip), computes per-class solar-phase
bin statistics, evaluates Métivier (2009) tidal-pattern reference ranking, tests
the loading and geometric hypotheses, performs an unmatched-event sensitivity
check, and writes all results to output/case-b3-results.json.
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
FOCAL_DIR = DATA_DIR / "focal-mechanism"
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
logger = logging.getLogger("case-b3-analysis")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SOLAR_YEAR_SECS: float = 31_557_600.0  # Julian constant: 365.25 × 86400
EXPECTED_CATALOG_N: int = 9210
BIN_COUNTS: list[int] = [16, 24, 32]
MECHANISM_CLASSES: list[str] = ["thrust", "normal", "strike_slip"]

# A1b baseline elevated phase intervals (from case-a4-results.json sub_b)
A1B_INTERVALS: list[tuple[float, float]] = [
    (0.1875, 0.25),   # Interval 1: ~Mar 10 – Apr 1  (March equinox)
    (0.625,  0.656),  # Interval 2: ~Aug 16 – Aug 28
    (0.875,  0.917),  # Interval 3: ~Nov 16 – Dec 1
]


# ---------------------------------------------------------------------------
# Section 1: Data loading and classification
# ---------------------------------------------------------------------------
def classify_rake(rake: float) -> str:
    """Classify focal mechanism from rake angle using standard quadrant boundaries.

    Args:
        rake: Rake angle in degrees.

    Returns:
        One of "thrust", "normal", "strike_slip", or "oblique".
    """
    # Normalise to (-180, 180]
    rake = float(rake)
    while rake > 180:
        rake -= 360
    while rake <= -180:
        rake += 360

    if 45 <= rake <= 135:
        return "thrust"
    elif -135 <= rake <= -45:
        return "normal"
    elif (-45 < rake <= 45) or (135 < rake <= 180) or (-180 <= rake < -135):
        return "strike_slip"
    else:
        return "oblique"


def load_focal_join() -> pd.DataFrame:
    """Load the GCMT focal mechanism join file (enriched catalog, n=9210).

    The file includes all base catalog columns (usgs_id, solar_secs, …) plus
    GCMT columns. solar_secs is used directly; no separate join required.

    Returns:
        DataFrame with 9210 rows.

    Raises:
        AssertionError: If row count is not EXPECTED_CATALOG_N.
    """
    path = FOCAL_DIR / "focal_join_global.csv"
    logger.info("Loading focal join from %s", path)
    df = pd.read_csv(path)
    n = len(df)
    assert n == EXPECTED_CATALOG_N, (
        f"Focal join: expected {EXPECTED_CATALOG_N} rows, got {n}"
    )
    logger.info("Focal join loaded: n=%d", n)
    logger.info("Columns: %s", list(df.columns))
    return df


def classify_mechanisms(df: pd.DataFrame) -> pd.DataFrame:
    """Classify each event by tectonic mechanism and compute solar phase.

    Classification logic (per spec §1):
    - Primary: use existing `mechanism` column if not null and not "oblique"
    - Fallback: re-classify from `rake` column if `mechanism` is null but
      `rake` is not null (log count)
    - Unmatched (`match_confidence` null): set class = "unmatched"
    - Oblique (from either column): set class = "oblique"

    Args:
        df: Raw focal join DataFrame.

    Returns:
        DataFrame with added `tectonic_class` and `solar_phase` columns.
    """
    df = df.copy()

    # Solar phase: (solar_secs / Julian_year) mod 1
    df["solar_phase"] = (df["solar_secs"] / SOLAR_YEAR_SECS) % 1.0

    # Initialise tectonic_class
    df["tectonic_class"] = "unmatched"

    # Matched events (match_confidence == "proximity")
    matched_mask = df["match_confidence"] == "proximity"
    n_matched = matched_mask.sum()
    logger.info("Events with proximity match: n=%d", n_matched)

    # Step 1: use existing mechanism column where available and not oblique
    has_mech = matched_mask & df["mechanism"].notna()
    for mech in ["thrust", "normal", "strike_slip"]:
        mech_mask = has_mech & (df["mechanism"] == mech)
        df.loc[mech_mask, "tectonic_class"] = mech
    # Mark oblique from mechanism column
    oblique_from_mech = has_mech & (df["mechanism"] == "oblique")
    df.loc[oblique_from_mech, "tectonic_class"] = "oblique"

    # Step 2: fallback re-classify from rake for matched events where mechanism is null
    needs_fallback = matched_mask & df["mechanism"].isna() & df["rake"].notna()
    n_reclassified = needs_fallback.sum()
    if n_reclassified > 0:
        logger.info(
            "Re-classifying %d matched events from rake (mechanism column was null)",
            n_reclassified,
        )
        for idx in df.index[needs_fallback]:
            rake_val = df.at[idx, "rake"]
            df.at[idx, "tectonic_class"] = classify_rake(rake_val)
    else:
        logger.info("No fallback rake re-classification needed")

    # Summary
    counts = df["tectonic_class"].value_counts()
    n_thrust = int(counts.get("thrust", 0))
    n_normal = int(counts.get("normal", 0))
    n_strike_slip = int(counts.get("strike_slip", 0))
    n_oblique = int(counts.get("oblique", 0))
    n_unmatched = int(counts.get("unmatched", 0))
    n_total_matched = n_thrust + n_normal + n_strike_slip + n_oblique

    logger.info(
        "Classification summary: thrust=%d  normal=%d  strike_slip=%d  "
        "oblique=%d  unmatched=%d  total_matched=%d",
        n_thrust, n_normal, n_strike_slip, n_oblique, n_unmatched, n_total_matched,
    )

    assert n_total_matched + n_unmatched == EXPECTED_CATALOG_N, (
        f"n_total_matched({n_total_matched}) + n_unmatched({n_unmatched}) "
        f"!= {EXPECTED_CATALOG_N}"
    )

    return df


# ---------------------------------------------------------------------------
# Section 2: Bin statistics helpers
# ---------------------------------------------------------------------------
def _bins_to_intervals(elevated_indices: list[int], k: int) -> list[dict[str, Any]]:
    """Merge adjacent elevated bin indices into contiguous phase intervals.

    Args:
        elevated_indices: Sorted list of elevated bin indices.
        k: Total number of bins.

    Returns:
        List of dicts with phase_start and phase_end for each contiguous run.
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
                "phase_end":   round((run_end + 1) / k, 6),
            })
            run_start = idx
            run_end = idx

    intervals.append({
        "phase_start": round(run_start / k, 6),
        "phase_end":   round((run_end + 1) / k, 6),
    })
    return intervals


def compute_bin_stats(phase_values: np.ndarray, k: int, label: str) -> dict[str, Any]:
    """Compute chi-square, Rayleigh, Cramér's V, and elevated-bin statistics.

    Uses the same procedure as Cases A4, B1, B2.

    Args:
        phase_values: Array of solar phase values in [0, 1).
        k: Number of phase bins.
        label: Human-readable label for logging.

    Returns:
        Dict with n, k, chi2, p_chi2, cramer_v, rayleigh_R, p_rayleigh,
        mean_phase, bin_counts, elevated_intervals.
    """
    n = len(phase_values)
    if n == 0:
        raise ValueError(f"Empty subset for label={label}, k={k}")
    if n < 200:
        logger.warning("Small sample: label=%s k=%d n=%d (< 200)", label, k, n)

    # Bin assignment
    bin_idx = np.floor(phase_values * k).astype(int)
    bin_idx = np.clip(bin_idx, 0, k - 1)

    obs = np.bincount(bin_idx, minlength=k).astype(float)
    expected = np.full(k, n / k)

    # Chi-square
    chi2_stat, p_chi2 = scipy.stats.chisquare(obs, expected)

    # Rayleigh statistic
    angles = 2.0 * np.pi * phase_values
    R = float(np.abs(np.mean(np.exp(1j * angles))))
    p_rayleigh = float(np.exp(-n * R ** 2))

    # Cramér's V
    cramer_v = float(np.sqrt(chi2_stat / (n * (k - 1))))

    # Mean phase angle → fraction of year
    mean_sin = float(np.mean(np.sin(angles)))
    mean_cos = float(np.mean(np.cos(angles)))
    mean_angle_rad = float(np.arctan2(mean_sin, mean_cos))
    mean_phase = float((mean_angle_rad / (2.0 * np.pi)) % 1.0)

    # Elevated bins at 1-SD threshold
    threshold = expected[0] + np.sqrt(expected[0])
    elevated_bin_indices = [int(i) for i in np.where(obs > threshold)[0]]
    elevated_intervals = _bins_to_intervals(elevated_bin_indices, k)

    # Overlap check with A1b intervals
    a1b_overlaps = []
    for ivl in elevated_intervals:
        for ref_start, ref_end in A1B_INTERVALS:
            overlap = max(
                0.0,
                min(ref_end, ivl["phase_end"]) - max(ref_start, ivl["phase_start"]),
            )
            ref_width = ref_end - ref_start
            if ref_width > 0 and overlap / ref_width > 0.5:
                a1b_overlaps.append({
                    "elevated": ivl,
                    "a1b_interval": {"phase_start": ref_start, "phase_end": ref_end},
                })

    result: dict[str, Any] = {
        "n": int(n),
        "k": int(k),
        "chi2": float(chi2_stat),
        "p_chi2": float(p_chi2),
        "cramer_v": cramer_v,
        "rayleigh_R": R,
        "p_rayleigh": p_rayleigh,
        "mean_phase": mean_phase,
        "bin_counts": obs.astype(int).tolist(),
        "expected_count": float(expected[0]),
        "threshold_1sd": float(threshold),
        "elevated_bin_indices": elevated_bin_indices,
        "elevated_intervals": elevated_intervals,
        "a1b_overlapping_intervals": a1b_overlaps,
    }

    logger.info(
        "%s k=%d: n=%d chi2=%.4f p=%.4e V=%.4f R=%.4f",
        label, k, n, chi2_stat, p_chi2, cramer_v, R,
    )
    return result


# ---------------------------------------------------------------------------
# Section 2: Per-mechanism bin statistics
# ---------------------------------------------------------------------------
def compute_mechanism_stats(df: pd.DataFrame) -> dict[str, Any]:
    """Compute per-mechanism solar phase bin statistics for k=16, 24, 32.

    Args:
        df: DataFrame with tectonic_class and solar_phase columns.

    Returns:
        Nested dict: mechanism_stats[class][k_key] = stat dict.
        Includes coverage counts.
    """
    counts = df["tectonic_class"].value_counts()

    coverage: dict[str, Any] = {
        "n_thrust":      int(counts.get("thrust", 0)),
        "n_normal":      int(counts.get("normal", 0)),
        "n_strike_slip": int(counts.get("strike_slip", 0)),
        "n_oblique":     int(counts.get("oblique", 0)),
        "n_unmatched":   int(counts.get("unmatched", 0)),
        "n_total":       EXPECTED_CATALOG_N,
    }

    mechanism_stats: dict[str, Any] = {"coverage": coverage}

    for mech in MECHANISM_CLASSES:
        subset = df[df["tectonic_class"] == mech]
        n_mech = len(subset)
        logger.info("Computing stats for mechanism=%s n=%d", mech, n_mech)
        phase_vals = subset["solar_phase"].values

        mech_entry: dict[str, Any] = {"n": int(n_mech)}
        for k in BIN_COUNTS:
            k_key = f"k{k}"
            mech_entry[k_key] = compute_bin_stats(phase_vals, k, f"{mech}")

        mechanism_stats[mech] = mech_entry

    return mechanism_stats


# ---------------------------------------------------------------------------
# Section 3: Métivier pattern comparison and prediction evaluation
# ---------------------------------------------------------------------------
def _intervals_overlap_a1b(elevated_intervals: list[dict]) -> bool:
    """Check whether any elevated interval overlaps Interval 1 (March equinox).

    Args:
        elevated_intervals: List of dicts with phase_start, phase_end.

    Returns:
        True if any elevated interval overlaps the March equinox A1b interval.
    """
    ref_start, ref_end = A1B_INTERVALS[0]  # 0.1875–0.25
    ref_width = ref_end - ref_start
    for ivl in elevated_intervals:
        overlap = max(
            0.0,
            min(ref_end, ivl["phase_end"]) - max(ref_start, ivl["phase_start"]),
        )
        if ref_width > 0 and overlap / ref_width > 0.5:
            return True
    return False


def evaluate_predictions(mechanism_stats: dict[str, Any]) -> dict[str, Any]:
    """Evaluate loading, geometric, and Métivier predictions using k=24 results.

    Args:
        mechanism_stats: Dict from compute_mechanism_stats (includes coverage).

    Returns:
        Dict with prediction_evaluation fields as specified in §3.
    """
    # k=24 stats for each mechanism
    thrust_k24 = mechanism_stats["thrust"]["k24"]
    normal_k24 = mechanism_stats["normal"]["k24"]
    ss_k24 = mechanism_stats["strike_slip"]["k24"]

    # 1. Cramér's V rank order
    cramer_vs = {
        "thrust":      thrust_k24["cramer_v"],
        "normal":      normal_k24["cramer_v"],
        "strike_slip": ss_k24["cramer_v"],
    }
    rank_order = sorted(cramer_vs, key=lambda m: cramer_vs[m], reverse=True)

    # 2. Métivier pattern: normal >= strike_slip >= thrust
    matches_metivier = (
        cramer_vs["normal"] >= cramer_vs["strike_slip"] >= cramer_vs["thrust"]
    )

    # 3. Anti-phase test for loading hypothesis
    mean_phase_thrust = thrust_k24["mean_phase"]
    mean_phase_normal = normal_k24["mean_phase"]
    raw_offset = abs(mean_phase_thrust - mean_phase_normal)
    # Wrap to [0, 0.5]
    phase_offset = raw_offset if raw_offset <= 0.5 else 1.0 - raw_offset

    if 0.4 <= phase_offset <= 0.6:
        thrust_normal_relationship = "anti_phased"
    elif phase_offset < 0.2:
        thrust_normal_relationship = "in_phase"
    else:
        thrust_normal_relationship = "other"

    # 4. All-classes equinox excess (overlap with A1b Interval 1)
    thrust_equinox = _intervals_overlap_a1b(thrust_k24["elevated_intervals"])
    normal_equinox = _intervals_overlap_a1b(normal_k24["elevated_intervals"])
    ss_equinox = _intervals_overlap_a1b(ss_k24["elevated_intervals"])
    geometric_equinox_in_all = bool(thrust_equinox and normal_equinox and ss_equinox)

    # Prediction support flags
    loading_supported = bool(
        thrust_normal_relationship == "anti_phased"
        and thrust_k24["p_chi2"] < 0.05  # some signal present in thrust to detect deficit
    )
    geometric_supported = bool(geometric_equinox_in_all)
    metivier_pattern_supported = bool(matches_metivier)

    # Primary conclusion
    n_supported = sum([loading_supported, geometric_supported, metivier_pattern_supported])
    if loading_supported and not geometric_supported and not metivier_pattern_supported:
        primary_conclusion = "loading"
    elif geometric_supported and not loading_supported and not metivier_pattern_supported:
        primary_conclusion = "geometric"
    elif metivier_pattern_supported and not loading_supported:
        primary_conclusion = "metivier_consistent"
    elif n_supported == 0:
        primary_conclusion = "ambiguous"
    else:
        primary_conclusion = "ambiguous"

    result: dict[str, Any] = {
        "cramer_v_rank_order": rank_order,
        "cramer_v_by_mechanism": cramer_vs,
        "matches_metivier_pattern": bool(matches_metivier),
        "phase_offset_thrust_normal": float(phase_offset),
        "mean_phase_thrust": float(mean_phase_thrust),
        "mean_phase_normal": float(mean_phase_normal),
        "thrust_normal_phase_relationship": thrust_normal_relationship,
        "equinox_in_thrust": bool(thrust_equinox),
        "equinox_in_normal": bool(normal_equinox),
        "equinox_in_strike_slip": bool(ss_equinox),
        "geometric_equinox_in_all": bool(geometric_equinox_in_all),
        "loading_supported": bool(loading_supported),
        "geometric_supported": bool(geometric_supported),
        "metivier_pattern_supported": bool(metivier_pattern_supported),
        "primary_conclusion": primary_conclusion,
    }

    logger.info(
        "Prediction evaluation: V_rank=%s metivier=%s phase_offset=%.4f "
        "relationship=%s equinox_all=%s conclusion=%s",
        rank_order, matches_metivier, phase_offset,
        thrust_normal_relationship, geometric_equinox_in_all, primary_conclusion,
    )
    return result


# ---------------------------------------------------------------------------
# Section 4: Unmatched event sensitivity check
# ---------------------------------------------------------------------------
def compute_unmatched_check(df: pd.DataFrame, full_catalog_cramer_v: float) -> dict[str, Any]:
    """Compute chi-square and Cramér's V for unmatched events at k=24.

    Assesses whether the GCMT matched subsample is representative of the full
    catalog. If unmatched events show similar signal strength, selection bias
    from 47% null rate is minimal.

    Args:
        df: Full DataFrame with tectonic_class and solar_phase.
        full_catalog_cramer_v: Cramér's V for all 9210 events at k=24
            (computed from full solar_phase column) for comparison.

    Returns:
        Dict with unmatched_check fields.
    """
    unmatched = df[df["tectonic_class"] == "unmatched"]
    n_unmatched = len(unmatched)
    logger.info("Unmatched event check: n=%d", n_unmatched)

    k = 24
    phase_vals = unmatched["solar_phase"].values
    stats = compute_bin_stats(phase_vals, k, "unmatched")

    # Similar signal: Cramér's V within 20% of full catalog V
    similar = bool(
        abs(stats["cramer_v"] - full_catalog_cramer_v) < 0.2 * full_catalog_cramer_v
        or stats["cramer_v"] < 0.05  # essentially no signal
    )

    result: dict[str, Any] = {
        "n_unmatched": int(n_unmatched),
        "full_catalog_cramer_v_k24": float(full_catalog_cramer_v),
        "k24_chi2":    float(stats["chi2"]),
        "k24_p":       float(stats["p_chi2"]),
        "k24_cramer_v": float(stats["cramer_v"]),
        "k24_rayleigh_R": float(stats["rayleigh_R"]),
        "k24_p_rayleigh": float(stats["p_rayleigh"]),
        "similar_to_full_catalog": similar,
    }
    logger.info(
        "Unmatched check: n=%d V=%.4f similar=%s",
        n_unmatched, stats["cramer_v"], similar,
    )
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Run the full Case B3 analysis and write results JSON."""
    logger.info("=== Case B3: Tectonic Regime Stratification ===")
    logger.info("Julian year constant: %.1f s", SOLAR_YEAR_SECS)

    # --- Section 1: Load and classify ---
    df = load_focal_join()
    df = classify_mechanisms(df)

    # Validate phase range
    assert df["solar_phase"].between(0.0, 1.0, inclusive="left").all(), (
        "Some solar_phase values outside [0, 1)"
    )

    # Full-catalog Cramér's V at k=24 for unmatched comparison
    full_phase = df["solar_phase"].values
    full_k24 = compute_bin_stats(full_phase, 24, "full_catalog")
    full_cramer_v = full_k24["cramer_v"]
    logger.info("Full catalog k=24 Cramér's V=%.4f", full_cramer_v)

    # --- Section 2: Per-mechanism bin statistics ---
    logger.info("Computing per-mechanism bin statistics...")
    mechanism_stats = compute_mechanism_stats(df)

    # --- Section 3: Prediction evaluation ---
    logger.info("Evaluating predictions...")
    prediction_evaluation = evaluate_predictions(mechanism_stats)

    # --- Section 4: Unmatched event check ---
    logger.info("Computing unmatched event sensitivity check...")
    unmatched_check = compute_unmatched_check(df, full_cramer_v)

    # --- Write results JSON ---
    coverage = mechanism_stats["coverage"]
    results: dict[str, Any] = {
        "case": "B3",
        "title": "Tectonic Regime Stratification",
        "solar_year_secs": SOLAR_YEAR_SECS,
        "catalog_n": EXPECTED_CATALOG_N,
        "a1b_intervals": [
            {"phase_start": s, "phase_end": e} for s, e in A1B_INTERVALS
        ],
        "full_catalog_k24": full_k24,
        "mechanism_stats": {
            "coverage": coverage,
            "thrust":      mechanism_stats["thrust"],
            "normal":      mechanism_stats["normal"],
            "strike_slip": mechanism_stats["strike_slip"],
        },
        "prediction_evaluation": prediction_evaluation,
        "unmatched_check": unmatched_check,
    }

    out_path = OUTPUT_DIR / "case-b3-results.json"
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", out_path)

    # Summary
    cov = coverage
    pred = prediction_evaluation
    logger.info(
        "SUMMARY: thrust=%d  normal=%d  strike_slip=%d  "
        "oblique=%d  unmatched=%d",
        cov["n_thrust"], cov["n_normal"], cov["n_strike_slip"],
        cov["n_oblique"], cov["n_unmatched"],
    )
    logger.info(
        "V rank=%s  métivier=%s  phase_offset=%.4f  conclusion=%s",
        pred["cramer_v_rank_order"],
        pred["matches_metivier_pattern"],
        pred["phase_offset_thrust_normal"],
        pred["primary_conclusion"],
    )


if __name__ == "__main__":
    main()
