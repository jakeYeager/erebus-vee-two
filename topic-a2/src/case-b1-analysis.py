"""Case B1: Hemisphere Stratification — Phase Symmetry Test.

Splits the ISC-GEM catalog by hemisphere and independently computes solar-phase
bin distributions for Northern Hemisphere (NH, latitude > 0°) and Southern
Hemisphere (SH, latitude < 0°) events. Performs four symmetry tests based on
the three-interval structure identified in Adhoc A1b.

Four symmetry tests:
  1. Global symmetry — do all three A1b intervals appear in both hemispheres?
  2. Interval 1 symmetry — is the March-equinox interval (0.1875–0.25) present
     in both hemispheres at the same phase position?
  3. Intervals 2 and 3 hemisphere-specificity — are they NH-only, SH-only,
     both, or neither?
  4. Half-cycle offset test — are any NH intervals offset by ~0.5 phase units
     in SH (hydrological hypothesis)?
"""

import json
import logging
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import scipy.stats

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
RAW_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
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
logger = logging.getLogger("case-b1-analysis")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SOLAR_YEAR_SECS = 31_557_600.0  # Julian year constant (confirmed: use uniformly)
EXPECTED_TOTAL = 9210
BIN_COUNTS = [16, 24, 32]

# A1b baseline elevated phase intervals (from Adhoc A1b)
A1B_INTERVALS = [
    (0.1875, 0.25),   # Interval 1: March equinox
    (0.625, 0.656),   # Interval 2: ~mid-August
    (0.875, 0.917),   # Interval 3: ~mid-November
]

# Bin width at k=24 (used for phase offset tolerance in symmetry tests)
BIN_WIDTH_K24 = 1.0 / 24  # ≈ 0.04167


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_catalog(path: Path) -> pd.DataFrame:
    """Load the ISC-GEM catalog and assert row count.

    Args:
        path: Path to the raw CSV file.

    Returns:
        Loaded DataFrame.
    """
    logger.info("Loading catalog from %s", path)
    df = pd.read_csv(path)
    n = len(df)
    logger.info("Loaded %d rows", n)
    assert n == EXPECTED_TOTAL, f"Expected {EXPECTED_TOTAL} rows, got {n}"
    return df


def split_hemispheres(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Split catalog into NH and SH subsets by latitude.

    Args:
        df: Full catalog DataFrame with 'latitude' column.

    Returns:
        Tuple of (df_nh, df_sh, n_equatorial).
    """
    df_nh = df[df["latitude"] > 0].copy()
    df_sh = df[df["latitude"] < 0].copy()
    n_equatorial = int((df["latitude"] == 0).sum())

    n_nh = len(df_nh)
    n_sh = len(df_sh)

    logger.info("NH events (latitude > 0): %d", n_nh)
    logger.info("SH events (latitude < 0): %d", n_sh)
    logger.info("Equatorial events (latitude == 0): %d", n_equatorial)

    assert n_nh + n_sh + n_equatorial == EXPECTED_TOTAL, (
        f"Partition mismatch: {n_nh} + {n_sh} + {n_equatorial} != {EXPECTED_TOTAL}"
    )
    return df_nh, df_sh, n_equatorial


def compute_phase(solar_secs: pd.Series) -> pd.Series:
    """Compute solar year phase using Julian year constant.

    Args:
        solar_secs: Series of solar seconds within the year.

    Returns:
        Phase values in [0, 1).
    """
    return (solar_secs / SOLAR_YEAR_SECS) % 1.0


# ---------------------------------------------------------------------------
# Per-hemisphere bin statistics
# ---------------------------------------------------------------------------
def compute_rayleigh(phases: np.ndarray) -> tuple[float, float, float]:
    """Compute Rayleigh R statistic, p-value, and mean phase.

    Args:
        phases: Array of phase values in [0, 1).

    Returns:
        Tuple of (R, p_rayleigh, mean_phase_fraction).
    """
    n = len(phases)
    angles = 2.0 * np.pi * phases
    mean_cos = float(np.mean(np.cos(angles)))
    mean_sin = float(np.mean(np.sin(angles)))
    R = float(np.sqrt(mean_cos**2 + mean_sin**2))
    p_rayleigh = float(np.exp(-n * R**2))
    mean_angle_rad = float(np.arctan2(mean_sin, mean_cos))
    mean_phase = float((mean_angle_rad / (2.0 * np.pi)) % 1.0)
    return R, p_rayleigh, mean_phase


def find_elevated_intervals(
    bin_counts: np.ndarray, k: int, n: int
) -> list[dict[str, float]]:
    """Identify contiguous elevated bins and merge them into intervals.

    A bin is elevated if count > E + sqrt(E) where E = n/k.

    Args:
        bin_counts: Array of observed counts per bin.
        k: Number of bins.
        n: Total event count.

    Returns:
        List of dicts with phase_start, phase_end, mean_phase.
    """
    expected = n / k
    threshold = expected + np.sqrt(expected)
    elevated_mask = bin_counts > threshold

    bin_width = 1.0 / k
    intervals: list[dict[str, float]] = []

    i = 0
    while i < k:
        if elevated_mask[i]:
            start_bin = i
            # Extend while adjacent bins are also elevated
            while i < k and elevated_mask[i]:
                i += 1
            end_bin = i  # exclusive
            phase_start = start_bin * bin_width
            phase_end = end_bin * bin_width
            mean_phase = (phase_start + phase_end) / 2.0
            intervals.append({
                "phase_start": round(phase_start, 6),
                "phase_end": round(phase_end, 6),
                "mean_phase": round(mean_phase, 6),
            })
        else:
            i += 1

    return intervals


def compute_hemisphere_stats(
    phases: np.ndarray, k: int
) -> dict[str, Any]:
    """Compute full bin statistics for one hemisphere at one bin count.

    Args:
        phases: Phase array in [0, 1).
        k: Number of bins.

    Returns:
        Dict with chi2, p_chi2, cramer_v, rayleigh_R, p_rayleigh,
        mean_phase, elevated_intervals, bin_counts.
    """
    n = len(phases)
    # Build bin counts using floor(phase * k) for phase-normalized binning
    bin_indices = np.floor(phases * k).astype(int)
    bin_indices = np.clip(bin_indices, 0, k - 1)  # guard against phase=1.0
    O = np.bincount(bin_indices, minlength=k).astype(float)

    expected = n / k
    E = np.full(k, expected)

    chi2, p_chi2 = scipy.stats.chisquare(O, E)
    cramer_v = float(np.sqrt(chi2 / (n * (k - 1))))

    R, p_rayleigh, mean_phase = compute_rayleigh(phases)

    elevated_intervals = find_elevated_intervals(O, k, n)

    return {
        "n": n,
        "k": k,
        "chi2": float(chi2),
        "p_chi2": float(p_chi2),
        "cramer_v": float(cramer_v),
        "rayleigh_R": float(R),
        "p_rayleigh": float(p_rayleigh),
        "mean_phase": float(mean_phase),
        "elevated_intervals": elevated_intervals,
        "bin_counts": O.astype(int).tolist(),
    }


def compute_all_hemisphere_stats(
    phases_nh: np.ndarray, phases_sh: np.ndarray
) -> dict[str, Any]:
    """Compute per-hemisphere statistics at k=16, 24, 32.

    Args:
        phases_nh: NH phase array.
        phases_sh: SH phase array.

    Returns:
        Dict keyed by hemisphere label with k16/k24/k32 results.
    """
    result: dict[str, Any] = {}
    for label, phases in [("nh", phases_nh), ("sh", phases_sh)]:
        result[label] = {}
        for k in BIN_COUNTS:
            k_key = f"k{k}"
            result[label][k_key] = compute_hemisphere_stats(phases, k)
            logger.info(
                "%s k=%d: chi2=%.3f p=%.4e V=%.4f intervals=%d",
                label.upper(), k,
                result[label][k_key]["chi2"],
                result[label][k_key]["p_chi2"],
                result[label][k_key]["cramer_v"],
                len(result[label][k_key]["elevated_intervals"]),
            )
    return result


# ---------------------------------------------------------------------------
# Symmetry tests
# ---------------------------------------------------------------------------
def intervals_overlap(
    s1: float, e1: float, s2: float, e2: float
) -> float:
    """Compute overlap fraction of [s1,e1] covered by [s2,e2].

    Args:
        s1: Start of first interval.
        e1: End of first interval.
        s2: Start of second interval.
        e2: End of second interval.

    Returns:
        Overlap as fraction of [s1,e1] width.
    """
    width = e1 - s1
    if width <= 0:
        return 0.0
    overlap = max(0.0, min(e1, e2) - max(s1, s2))
    return overlap / width


def interval_in_elevated(
    baseline_start: float,
    baseline_end: float,
    elevated_intervals: list[dict[str, float]],
    overlap_threshold: float = 0.5,
) -> bool:
    """Check whether a baseline interval is covered by any elevated interval.

    Args:
        baseline_start: Baseline interval start phase.
        baseline_end: Baseline interval end phase.
        elevated_intervals: List of elevated interval dicts.
        overlap_threshold: Minimum overlap fraction to count as matching.

    Returns:
        True if any elevated interval overlaps baseline by > threshold.
    """
    for ei in elevated_intervals:
        frac = intervals_overlap(
            baseline_start, baseline_end,
            ei["phase_start"], ei["phase_end"],
        )
        if frac > overlap_threshold:
            return True
    return False


def test_1_global_symmetry(
    elevated_nh: list[dict], elevated_sh: list[dict]
) -> dict[str, Any]:
    """Test 1: Do all three A1b intervals appear in both hemispheres?

    Args:
        elevated_nh: NH elevated intervals at k=24.
        elevated_sh: SH elevated intervals at k=24.

    Returns:
        Per-interval presence dict and overall classification.
    """
    interval_results = {}
    for idx, (ps, pe) in enumerate(A1B_INTERVALS, start=1):
        in_nh = interval_in_elevated(ps, pe, elevated_nh)
        in_sh = interval_in_elevated(ps, pe, elevated_sh)
        symmetric = in_nh and in_sh
        interval_results[f"interval_{idx}"] = {
            "in_nh": bool(in_nh),
            "in_sh": bool(in_sh),
            "symmetric": bool(symmetric),
        }

    all_symmetric = all(v["symmetric"] for v in interval_results.values())
    any_symmetric = any(v["symmetric"] for v in interval_results.values())

    if all_symmetric:
        classification = "fully symmetric"
    elif any_symmetric:
        classification = "partially symmetric"
    else:
        classification = "asymmetric"

    return {
        "interval_1": interval_results["interval_1"],
        "interval_2": interval_results["interval_2"],
        "interval_3": interval_results["interval_3"],
        "classification": classification,
    }


def test_2_interval_1_symmetry(
    elevated_nh: list[dict], elevated_sh: list[dict]
) -> dict[str, Any]:
    """Test 2: Is interval 1 (March equinox) symmetric in both hemispheres?

    Args:
        elevated_nh: NH elevated intervals at k=24.
        elevated_sh: SH elevated intervals at k=24.

    Returns:
        Dict with presence flags and phase offset.
    """
    ps, pe = A1B_INTERVALS[0]  # Interval 1: 0.1875–0.25
    tolerance = BIN_WIDTH_K24  # ±1 bin width at k=24

    in_nh = interval_in_elevated(ps, pe, elevated_nh)
    in_sh = interval_in_elevated(ps, pe, elevated_sh)

    # Find mean phase of the matching elevated interval in each hemisphere
    def get_interval_mean_phase(elevated: list[dict]) -> Optional[float]:
        best_overlap = 0.0
        best_mean: float | None = None
        for ei in elevated:
            frac = intervals_overlap(ps, pe, ei["phase_start"], ei["phase_end"])
            if frac > 0.5 and frac > best_overlap:
                best_overlap = frac
                best_mean = ei["mean_phase"]
        return best_mean

    nh_mean = get_interval_mean_phase(elevated_nh) if in_nh else None
    sh_mean = get_interval_mean_phase(elevated_sh) if in_sh else None

    if nh_mean is not None and sh_mean is not None:
        phase_offset = float(nh_mean - sh_mean)
    else:
        phase_offset = None

    return {
        "interval_1_nh": bool(in_nh),
        "interval_1_sh": bool(in_sh),
        "phase_offset": phase_offset,
    }


def test_3_interval_23_specificity(
    elevated_nh: list[dict], elevated_sh: list[dict]
) -> dict[str, Any]:
    """Test 3: Are intervals 2 and 3 hemisphere-specific?

    Args:
        elevated_nh: NH elevated intervals at k=24.
        elevated_sh: SH elevated intervals at k=24.

    Returns:
        Dict with hemisphere classification for intervals 2 and 3.
    """
    def classify(ps: float, pe: float) -> str:
        in_nh = interval_in_elevated(ps, pe, elevated_nh)
        in_sh = interval_in_elevated(ps, pe, elevated_sh)
        if in_nh and in_sh:
            return "both"
        elif in_nh:
            return "nh_only"
        elif in_sh:
            return "sh_only"
        else:
            return "neither"

    ps2, pe2 = A1B_INTERVALS[1]  # Interval 2: 0.625–0.656
    ps3, pe3 = A1B_INTERVALS[2]  # Interval 3: 0.875–0.917

    return {
        "interval_2": {"hemisphere": classify(ps2, pe2)},
        "interval_3": {"hemisphere": classify(ps3, pe3)},
    }


def compute_half_cycle_offset(
    elevated_nh: list[dict], elevated_sh: list[dict]
) -> dict[str, Any]:
    """Test 4: Are any NH intervals offset by ~0.5 cycles in SH?

    Args:
        elevated_nh: NH elevated intervals at k=24.
        elevated_sh: SH elevated intervals at k=24.

    Returns:
        Dict with any_half_cycle_offset_found and details list.
    """
    tolerance = BIN_WIDTH_K24  # 1 bin width at k=24

    details = []
    any_found = False

    for ei_nh in elevated_nh:
        nh_center = ei_nh["mean_phase"]
        counterpart = (nh_center + 0.5) % 1.0

        # Check whether any SH interval is within tolerance of counterpart
        for ei_sh in elevated_sh:
            sh_center = ei_sh["mean_phase"]
            offset = abs(sh_center - counterpart)
            # Wrap-around distance on [0,1)
            offset = min(offset, 1.0 - offset)
            within_tolerance = bool(offset <= tolerance)
            if within_tolerance:
                any_found = True

            details.append({
                "nh_interval_center": round(float(nh_center), 6),
                "expected_sh_counterpart": round(float(counterpart), 6),
                "sh_interval_center": round(float(sh_center), 6),
                "offset": round(float(offset), 6),
                "within_tolerance": within_tolerance,
            })

    return {
        "any_half_cycle_offset_found": bool(any_found),
        "details": details,
    }


def run_symmetry_tests(
    hemi_stats: dict[str, Any]
) -> dict[str, Any]:
    """Run all four symmetry tests using k=24 elevated intervals.

    Args:
        hemi_stats: Dict with 'nh' and 'sh' stats keyed by k16/k24/k32.

    Returns:
        Dict with all four test results.
    """
    elevated_nh_k24 = hemi_stats["nh"]["k24"]["elevated_intervals"]
    elevated_sh_k24 = hemi_stats["sh"]["k24"]["elevated_intervals"]

    logger.info("Running symmetry tests with k=24 elevated intervals")
    logger.info("NH k=24 elevated intervals: %d", len(elevated_nh_k24))
    logger.info("SH k=24 elevated intervals: %d", len(elevated_sh_k24))

    t1 = test_1_global_symmetry(elevated_nh_k24, elevated_sh_k24)
    t2 = test_2_interval_1_symmetry(elevated_nh_k24, elevated_sh_k24)
    t3 = test_3_interval_23_specificity(elevated_nh_k24, elevated_sh_k24)
    t4 = compute_half_cycle_offset(elevated_nh_k24, elevated_sh_k24)

    logger.info("Test 1 (global symmetry): %s", t1["classification"])
    logger.info(
        "Test 2 (interval 1): NH=%s SH=%s offset=%s",
        t2["interval_1_nh"], t2["interval_1_sh"], t2["phase_offset"]
    )
    logger.info("Test 3 interval_2: %s", t3["interval_2"]["hemisphere"])
    logger.info("Test 3 interval_3: %s", t3["interval_3"]["hemisphere"])
    logger.info("Test 4 half-cycle offset found: %s", t4["any_half_cycle_offset_found"])

    return {
        "test_1_global": t1,
        "test_2_interval_1": t2,
        "test_3_interval_23_specificity": t3,
        "test_4_half_cycle_offset": t4,
    }


# ---------------------------------------------------------------------------
# Prediction matching
# ---------------------------------------------------------------------------
def evaluate_predictions(symmetry_tests: dict[str, Any]) -> dict[str, Any]:
    """Evaluate which mechanistic hypothesis is best supported.

    Geometric: all three intervals in both hemispheres at the same phase.
    Hydrological: NH/SH peaks offset by ~0.5 cycle; intervals 2+3 hemisphere-specific.
    Mixed: Interval 1 symmetric; intervals 2 and/or 3 hemisphere-specific.

    Args:
        symmetry_tests: Full symmetry test results dict.

    Returns:
        Dict with support classification for each hypothesis and primary conclusion.
    """
    t1 = symmetry_tests["test_1_global"]
    t2 = symmetry_tests["test_2_interval_1"]
    t3 = symmetry_tests["test_3_interval_23_specificity"]
    t4 = symmetry_tests["test_4_half_cycle_offset"]

    # Geometric: all three intervals symmetric AND no half-cycle offset
    all_symmetric = (
        t1["interval_1"]["symmetric"]
        and t1["interval_2"]["symmetric"]
        and t1["interval_3"]["symmetric"]
    )
    no_half_cycle = not t4["any_half_cycle_offset_found"]

    if all_symmetric and no_half_cycle:
        geometric = "supported"
    elif all_symmetric or (
        t1["interval_1"]["symmetric"] and not t4["any_half_cycle_offset_found"]
    ):
        geometric = "partially supported"
    else:
        geometric = "not supported"

    # Hydrological: half-cycle offset found AND interval 2 or 3 is hemisphere-specific
    i2_specific = t3["interval_2"]["hemisphere"] in ("nh_only", "sh_only")
    i3_specific = t3["interval_3"]["hemisphere"] in ("nh_only", "sh_only")
    half_cycle_found = t4["any_half_cycle_offset_found"]

    if half_cycle_found and (i2_specific or i3_specific):
        hydrological = "supported"
    elif half_cycle_found or (i2_specific and i3_specific):
        hydrological = "partially supported"
    else:
        hydrological = "not supported"

    # Mixed: interval 1 symmetric AND intervals 2 and/or 3 hemisphere-specific
    i1_symmetric = t1["interval_1"]["symmetric"] or (
        t2["interval_1_nh"] and t2["interval_1_sh"]
    )

    if i1_symmetric and (i2_specific or i3_specific):
        mixed = "supported"
    elif i1_symmetric or (i2_specific or i3_specific):
        mixed = "partially supported"
    else:
        mixed = "not supported"

    # Primary conclusion
    support_scores = {
        "geometric": {"supported": 3, "partially supported": 2, "not supported": 1}[geometric],
        "hydrological": {"supported": 3, "partially supported": 2, "not supported": 1}[hydrological],
        "mixed": {"supported": 3, "partially supported": 2, "not supported": 1}[mixed],
    }
    top_score = max(support_scores.values())
    top_hypotheses = [h for h, s in support_scores.items() if s == top_score]

    if len(top_hypotheses) == 1:
        primary = top_hypotheses[0]
    elif len(top_hypotheses) > 1:
        primary = "ambiguous"
    else:
        primary = "ambiguous"

    logger.info(
        "Prediction support: geometric=%s hydrological=%s mixed=%s primary=%s",
        geometric, hydrological, mixed, primary,
    )

    return {
        "geometric": geometric,
        "hydrological": hydrological,
        "mixed": mixed,
        "primary_conclusion": primary,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load catalog, split hemispheres, compute stats, run tests, write JSON."""
    # --- Load and split ---
    df = load_catalog(RAW_PATH)
    df_nh, df_sh, n_equatorial = split_hemispheres(df)

    # --- Compute phases ---
    phases_nh = compute_phase(df_nh["solar_secs"]).to_numpy()
    phases_sh = compute_phase(df_sh["solar_secs"]).to_numpy()

    logger.info("NH phase range: [%.6f, %.6f]", phases_nh.min(), phases_nh.max())
    logger.info("SH phase range: [%.6f, %.6f]", phases_sh.min(), phases_sh.max())

    # --- Per-hemisphere bin statistics ---
    logger.info("Computing per-hemisphere bin statistics at k=16, 24, 32")
    hemi_stats = compute_all_hemisphere_stats(phases_nh, phases_sh)

    # --- Symmetry tests ---
    symmetry_tests = run_symmetry_tests(hemi_stats)

    # --- Prediction matching ---
    prediction_support = evaluate_predictions(symmetry_tests)

    # --- Assemble results ---
    results: dict[str, Any] = {
        "case": "B1",
        "title": "Hemisphere Stratification — Phase Symmetry Test",
        "solar_year_secs": SOLAR_YEAR_SECS,
        "hemisphere_stats": {
            "n_nh": int(len(df_nh)),
            "n_sh": int(len(df_sh)),
            "n_equatorial": int(n_equatorial),
            "nh": hemi_stats["nh"],
            "sh": hemi_stats["sh"],
        },
        "symmetry_tests": symmetry_tests,
        "prediction_support": prediction_support,
    }

    # --- Write JSON ---
    output_path = OUTPUT_DIR / "case-b1-results.json"
    with open(output_path, "w") as fh:
        json.dump(results, fh, indent=2)
    logger.info("Results written to %s", output_path)


if __name__ == "__main__":
    main()
