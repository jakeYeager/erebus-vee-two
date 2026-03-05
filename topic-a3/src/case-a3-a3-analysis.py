"""
Case A3.A3: Phase-Concentration Audit

Computes signed chi-square influence for each event in the full catalog,
runs a permutation baseline, sequential removal curves (elevated and suppressed
bins), and representativeness tests to characterize whether the solar-phase
signal is concentrated in a small number of identifiable events or genuinely
diffuse across the catalog.
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

CATALOG_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"
GK_MAINSHOCK_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_gk-seq_global.csv"
GK_AFTERSHOCK_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_gk-seq_global.csv"
REAS_MAINSHOCK_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_reas-seq_global.csv"
REAS_AFTERSHOCK_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_reas-seq_global.csv"
A1B_MAINSHOCK_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "mainshocks_a1b-seq_global.csv"
A1B_AFTERSHOCK_PATH = BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_a1b-seq_global.csv"

OUTPUT_PATH = BASE_DIR / "output" / "case-a3-a3-results.json"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
K: int = 24
JULIAN_YEAR_SECS: float = 31_557_600.0
N_PERMUTATIONS: int = 1000
N_RANDOM_BASELINES: int = 500
ELEVATED_BINS: list[int] = [4, 5, 6, 7, 15, 19, 21]
SUPPRESSED_BINS: list[int] = [2, 8, 10, 11, 12, 13, 16, 18, 22]
A1B_INTERVALS: dict[str, list[int]] = {
    "interval_1": [4, 5],
    "interval_2": [15],
    "interval_3": [21],
}
TOP_N_REPRESENTATIVENESS: int = 100


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def compute_signed_influence(phases: np.ndarray, k: int = 24) -> np.ndarray:
    """
    Compute each event's signed chi-square influence: the change in total chi-square
    when that event is removed.

    Parameters
    ----------
    phases : np.ndarray
        Phase values in [0, 1) for all events (length n).
    k : int
        Number of bins.

    Returns
    -------
    np.ndarray of length n. Positive = event is in an elevated bin (removal
    reduces chi-square). Negative = event is in a suppressed bin (removal
    increases chi-square, deepening the trough).

    Notes
    -----
    All events in the same bin receive identical influence values.
    Influence for an event in bin j:
        influence_j = (2 * (obs_j - expected) - 1) / expected
    where expected = n / k. This is the exact first-order change in chi-square.
    """
    n = len(phases)
    bin_indices = np.floor(phases * k).astype(int) % k
    obs = np.bincount(bin_indices, minlength=k).astype(float)
    expected = n / k
    # Per-bin influence value
    bin_influence = (2.0 * (obs - expected) - 1.0) / expected
    # Map each event to its bin's influence
    event_influence = bin_influence[bin_indices]
    return event_influence


def run_permutation_baseline(
    phases: np.ndarray,
    k: int = 24,
    n_permutations: int = 1000,
    rng_seed: int = 42,
) -> dict:
    """
    Permute solar phases and compute per-bin count distributions.

    Parameters
    ----------
    phases : np.ndarray
        Observed phase values (length n).
    k : int
        Number of bins.
    n_permutations : int
        Number of random phase shuffles.
    rng_seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict with:
        "perm_bin_counts": np.ndarray shape (n_permutations, k)
        "perm_chi2": np.ndarray shape (n_permutations,)
        "bin_p5": np.ndarray shape (k,) — 5th percentile per bin
        "bin_p95": np.ndarray shape (k,) — 95th percentile per bin
        "sig_elevated_bins": list of bin indices where actual obs > bin_p95
        "sig_suppressed_bins": list of bin indices where actual obs < bin_p5
    """
    rng = np.random.default_rng(rng_seed)
    n = len(phases)
    expected = n / k

    # Actual bin counts
    actual_bin_indices = np.floor(phases * k).astype(int) % k
    actual_obs = np.bincount(actual_bin_indices, minlength=k).astype(float)

    perm_bin_counts = np.zeros((n_permutations, k), dtype=float)
    perm_chi2 = np.zeros(n_permutations, dtype=float)

    for i in range(n_permutations):
        # Draw random phases from Uniform[0, 1) — proper null hypothesis of uniform distribution.
        # Shuffling observed phases would preserve bin counts identically (since each observed
        # phase maps deterministically to a single bin), so we use random draws instead.
        rand_phases = rng.uniform(0.0, 1.0, size=n)
        bin_idx = np.floor(rand_phases * k).astype(int) % k
        counts = np.bincount(bin_idx, minlength=k).astype(float)
        perm_bin_counts[i] = counts
        perm_chi2[i] = np.sum((counts - expected) ** 2 / expected)

    bin_p5 = np.percentile(perm_bin_counts, 5, axis=0)
    bin_p95 = np.percentile(perm_bin_counts, 95, axis=0)

    sig_elevated_bins = [j for j in range(k) if actual_obs[j] > bin_p95[j]]
    sig_suppressed_bins = [j for j in range(k) if actual_obs[j] < bin_p5[j]]

    return {
        "perm_bin_counts": perm_bin_counts,
        "perm_chi2": perm_chi2,
        "bin_p5": bin_p5,
        "bin_p95": bin_p95,
        "sig_elevated_bins": sig_elevated_bins,
        "sig_suppressed_bins": sig_suppressed_bins,
    }


def sequential_removal_curve(
    df: pd.DataFrame,
    removal_order: list[str],
    k: int = 24,
    track_bins: list[int] | None = None,
) -> list[dict]:
    """
    Sequentially remove events and track chi-square and bin z-scores at each step.

    Parameters
    ----------
    df : pd.DataFrame
        Full catalog with 'phase' and 'usgs_id' columns.
    removal_order : list[str]
        usgs_id values in the order they should be removed.
    k : int
        Number of bins.
    track_bins : list[int] | None
        Specific bins to track z-scores for. If None, tracks ELEVATED_BINS and SUPPRESSED_BINS.

    Returns
    -------
    List of dicts, one per removal step, each containing:
        step, n_removed, n_remaining, chi2, p_chi2,
        elevated_bin_zscores (dict bin->z), suppressed_bin_zscores (dict bin->z),
        pct_catalog_removed.
    """
    if track_bins is None:
        track_bins = ELEVATED_BINS + SUPPRESSED_BINS

    n_total = len(df)
    remaining_ids: set[str] = set(df["usgs_id"].tolist())

    # Vectorized approach: pre-compute bin assignment for all events
    id_to_bin: dict[str, int] = {
        uid: int(np.floor(phase_val * k)) % k
        for uid, phase_val in zip(df["usgs_id"].tolist(), df["phase"].tolist())
    }

    # Start with all bins populated
    bin_counts = np.bincount(
        [id_to_bin[uid] for uid in df["usgs_id"].tolist()],
        minlength=k
    ).astype(float)

    def make_step_dict(step: int, n_removed: int, bc: np.ndarray) -> dict:
        n_rem = n_total - n_removed
        if n_rem <= 0:
            return {
                "step": step, "n_removed": n_removed, "n_remaining": 0,
                "chi2": 0.0, "p_chi2": 1.0,
                "elevated_bin_zscores": {b: 0.0 for b in ELEVATED_BINS},
                "suppressed_bin_zscores": {b: 0.0 for b in SUPPRESSED_BINS},
                "pct_catalog_removed": 100.0,
            }
        expected_rem = n_rem / k
        chi2_val = float(np.sum((bc - expected_rem) ** 2 / expected_rem))
        dof = k - 1
        p_val = float(scipy.stats.chi2.sf(chi2_val, dof))
        z_arr = (bc - expected_rem) / np.sqrt(expected_rem)
        return {
            "step": step, "n_removed": n_removed, "n_remaining": n_rem,
            "chi2": chi2_val, "p_chi2": p_val,
            "elevated_bin_zscores": {b2: float(z_arr[b2]) for b2 in ELEVATED_BINS},
            "suppressed_bin_zscores": {b2: float(z_arr[b2]) for b2 in SUPPRESSED_BINS},
            "pct_catalog_removed": (n_removed / n_total) * 100.0,
        }

    results = [make_step_dict(0, 0, bin_counts.copy())]

    # Filter removal_order to only those in the dataframe
    removal_order_filtered = [uid for uid in removal_order if uid in remaining_ids]

    for i, uid in enumerate(removal_order_filtered):
        b = id_to_bin.get(uid)
        if b is not None:
            bin_counts[b] -= 1
        results.append(make_step_dict(i + 1, i + 1, bin_counts.copy()))

    return results


def representativeness_test(
    df: pd.DataFrame,
    top_n: int = 100,
) -> dict:
    """
    Test whether top-influence events are anomalous relative to the full catalog
    and the signal-bearing stratum.

    Parameters
    ----------
    df : pd.DataFrame
        Full catalog with chi2_influence, tectonic_class, depth_band columns.
    top_n : int
        Number of top positive and top |negative| influence events to test.

    Returns
    -------
    dict with representativeness results for four groups:
        top_positive_50, top_positive_100, top_negative_50, top_negative_100.
    """
    n_full = len(df)

    # Signal-bearing stratum: continental + mid-crustal
    signal_stratum = df[
        (df["tectonic_class"] == "continental") | (df["depth_band"] == "mid_crustal")
    ]
    n_stratum = len(signal_stratum)

    # Sort by influence
    sorted_pos = df.sort_values("chi2_influence", ascending=False)
    sorted_neg = df.sort_values("chi2_influence", ascending=True)

    tectonic_cats = ["continental", "transitional", "oceanic"]
    depth_cats = ["shallow", "mid_crustal", "intermediate", "deep"]
    role_cats = ["isolated_mainshock", "mainshock_with_sequence", "aftershock"]

    def full_cat_tectonic() -> np.ndarray:
        return np.array([
            (df["tectonic_class"] == c).sum() for c in tectonic_cats
        ], dtype=float)

    def full_cat_depth() -> np.ndarray:
        return np.array([
            (df["depth_band"] == c).sum() for c in depth_cats
        ], dtype=float)

    def stratum_tectonic() -> np.ndarray:
        return np.array([
            (signal_stratum["tectonic_class"] == c).sum() for c in tectonic_cats
        ], dtype=float)

    def stratum_depth() -> np.ndarray:
        return np.array([
            (signal_stratum["depth_band"] == c).sum() for c in depth_cats
        ], dtype=float)

    def test_group(group_df: pd.DataFrame, n: int) -> dict:
        tect_dist = {c: int((group_df["tectonic_class"] == c).sum()) for c in tectonic_cats}
        depth_dist = {c: int((group_df["depth_band"] == c).sum()) for c in depth_cats}
        role_dist = {c: int((group_df["sequence_role"] == c).sum()) for c in role_cats}

        group_tect = np.array([tect_dist[c] for c in tectonic_cats], dtype=float)
        group_depth = np.array([depth_dist[c] for c in depth_cats], dtype=float)

        # Chi-square vs full catalog (tectonic)
        full_tect_expected = full_cat_tectonic() * (n / n_full)
        # Remove zero-expected categories
        mask_t = full_tect_expected > 0
        if mask_t.sum() > 1:
            chi2_t_full, p_t_full = scipy.stats.chisquare(
                group_tect[mask_t], full_tect_expected[mask_t]
            )
        else:
            chi2_t_full, p_t_full = 0.0, 1.0

        # Chi-square vs full catalog (depth)
        full_depth_expected = full_cat_depth() * (n / n_full)
        mask_d = full_depth_expected > 0
        if mask_d.sum() > 1:
            chi2_d_full, p_d_full = scipy.stats.chisquare(
                group_depth[mask_d], full_depth_expected[mask_d]
            )
        else:
            chi2_d_full, p_d_full = 0.0, 1.0

        # Combined p vs full: use tectonic (primary)
        p_vs_full = float(p_t_full)

        # Chi-square vs signal stratum (tectonic)
        strat_tect_expected = stratum_tectonic() * (n / n_stratum)
        mask_st = strat_tect_expected > 0
        if mask_st.sum() > 1:
            chi2_t_strat, p_t_strat = scipy.stats.chisquare(
                group_tect[mask_st], strat_tect_expected[mask_st]
            )
        else:
            chi2_t_strat, p_t_strat = 0.0, 1.0

        # Chi-square vs signal stratum (depth)
        strat_depth_expected = stratum_depth() * (n / n_stratum)
        mask_sd = strat_depth_expected > 0
        if mask_sd.sum() > 1:
            chi2_d_strat, p_d_strat = scipy.stats.chisquare(
                group_depth[mask_sd], strat_depth_expected[mask_sd]
            )
        else:
            chi2_d_strat, p_d_strat = 0.0, 1.0

        p_vs_signal_stratum = float(p_t_strat)
        representative_of_signal_stratum = p_vs_signal_stratum > 0.05

        return {
            "n": n,
            "tectonic_dist": tect_dist,
            "depth_dist": depth_dist,
            "sequence_role_dist": role_dist,
            "p_vs_full_tectonic": float(p_t_full),
            "p_vs_full_depth": float(p_d_full),
            "p_vs_full": p_vs_full,
            "p_vs_signal_stratum_tectonic": float(p_t_strat),
            "p_vs_signal_stratum_depth": float(p_d_strat),
            "p_vs_signal_stratum": p_vs_signal_stratum,
            "representative_of_signal_stratum": representative_of_signal_stratum,
        }

    results = {
        "top_positive_50": test_group(sorted_pos.head(50), 50),
        "top_positive_100": test_group(sorted_pos.head(100), 100),
        "top_negative_50": test_group(sorted_neg.head(50), 50),
        "top_negative_100": test_group(sorted_neg.head(100), 100),
    }
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    logger.info("=== Case A3.A3: Phase-Concentration Audit ===")

    # ------------------------------------------------------------------
    # Section 1: Load data
    # ------------------------------------------------------------------
    logger.info("Loading full catalog...")
    df = pd.read_csv(CATALOG_PATH)
    assert len(df) == 9210, f"Expected 9210 events, got {len(df)}"
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    assert df["phase"].between(0.0, 1.0, inclusive="left").all(), "Phase values out of [0, 1)"
    logger.info(f"Full catalog loaded: n={len(df)}")

    logger.info("Loading GSHHG classification...")
    gshhg = pd.read_csv(GSHHG_PATH)
    assert len(gshhg) == 9210, f"Expected 9210 GSHHG rows, got {len(gshhg)}"
    df = df.merge(gshhg[["usgs_id", "ocean_class", "dist_to_coast_km"]], on="usgs_id", how="left")
    assert df["ocean_class"].isna().sum() == 0, "NaN values found in ocean_class after merge"
    logger.info("GSHHG merged successfully")

    # Tectonic class
    def classify_tectonic(dist: float) -> str:
        if dist <= 50.0:
            return "continental"
        elif dist <= 200.0:
            return "transitional"
        else:
            return "oceanic"

    df["tectonic_class"] = df["dist_to_coast_km"].apply(classify_tectonic)

    # Depth band
    def classify_depth(d: float) -> str:
        if d < 20.0:
            return "shallow"
        elif d < 70.0:
            return "mid_crustal"
        elif d < 300.0:
            return "intermediate"
        else:
            return "deep"

    df["depth_band"] = df["depth"].apply(classify_depth)

    # ------------------------------------------------------------------
    # Load sequence catalogs
    # ------------------------------------------------------------------
    logger.info("Loading sequence catalogs...")
    gk_main = pd.read_csv(GK_MAINSHOCK_PATH)
    gk_after = pd.read_csv(GK_AFTERSHOCK_PATH)
    reas_main = pd.read_csv(REAS_MAINSHOCK_PATH)
    reas_after = pd.read_csv(REAS_AFTERSHOCK_PATH)
    a1b_main = pd.read_csv(A1B_MAINSHOCK_PATH)
    a1b_after = pd.read_csv(A1B_AFTERSHOCK_PATH)

    logger.info(f"GK aftershocks: n={len(gk_after)}")
    logger.info(f"Reasenberg aftershocks: n={len(reas_after)}")
    logger.info(f"A1b aftershocks: n={len(a1b_after)}")

    # Build sequence membership sets
    # "mainshock_with_aftershocks": parent_id values in any aftershock file
    mainshock_with_aftershocks: set[str] = set()
    aftershock_member: set[str] = set()

    for adf in [gk_after, reas_after, a1b_after]:
        mainshock_with_aftershocks.update(adf["parent_id"].dropna().astype(str).tolist())
        aftershock_member.update(adf["usgs_id"].dropna().astype(str).tolist())

    def classify_role(uid: str) -> str:
        if uid in aftershock_member:
            return "aftershock"
        elif uid in mainshock_with_aftershocks:
            return "mainshock_with_sequence"
        else:
            return "isolated_mainshock"

    df["sequence_role"] = df["usgs_id"].astype(str).apply(classify_role)
    assert df["sequence_role"].isna().sum() == 0, "NaN in sequence_role"

    n_isolated = int((df["sequence_role"] == "isolated_mainshock").sum())
    n_mainshock_seq = int((df["sequence_role"] == "mainshock_with_sequence").sum())
    n_aftershock = int((df["sequence_role"] == "aftershock").sum())
    assert n_isolated + n_mainshock_seq + n_aftershock == 9210, "Partition does not sum to 9210"

    logger.info(f"Sequence roles: isolated={n_isolated}, mainshock_seq={n_mainshock_seq}, aftershock={n_aftershock}")

    # ------------------------------------------------------------------
    # Section 2: Signed chi-square influence
    # ------------------------------------------------------------------
    logger.info("Computing signed chi-square influence...")
    phases = df["phase"].values
    df["chi2_influence"] = compute_signed_influence(phases, k=K)

    # Bin-level summary
    n = len(df)
    expected = n / K
    bin_indices = np.floor(phases * K).astype(int) % K
    obs = np.bincount(bin_indices, minlength=K).astype(float)
    z_scores = (obs - expected) / np.sqrt(expected)
    bin_influence_vals = (2.0 * (obs - expected) - 1.0) / expected

    bin_summary = []
    for b in range(K):
        bin_type = "neutral"
        if b in ELEVATED_BINS:
            bin_type = "elevated"
        elif b in SUPPRESSED_BINS:
            bin_type = "suppressed"

        bin_mask = bin_indices == b
        n_iso = int((df.loc[bin_mask, "sequence_role"] == "isolated_mainshock").sum())
        n_ms = int((df.loc[bin_mask, "sequence_role"] == "mainshock_with_sequence").sum())
        n_as = int((df.loc[bin_mask, "sequence_role"] == "aftershock").sum())

        bin_summary.append({
            "bin": b,
            "obs": int(obs[b]),
            "expected": float(expected),
            "z": float(z_scores[b]),
            "bin_influence": float(bin_influence_vals[b]),
            "bin_type": bin_type,
            "n_isolated": n_iso,
            "n_mainshock_seq": n_ms,
            "n_aftershock": n_as,
        })

    logger.info("Bin summary computed.")
    for entry in bin_summary:
        if entry["bin_type"] != "neutral":
            logger.info(f"  Bin {entry['bin']}: obs={entry['obs']}, z={entry['z']:.3f}, type={entry['bin_type']}")

    # ------------------------------------------------------------------
    # Section 3: Permutation baseline
    # ------------------------------------------------------------------
    logger.info(f"Running permutation baseline (n={N_PERMUTATIONS})...")
    perm_result = run_permutation_baseline(phases, k=K, n_permutations=N_PERMUTATIONS, rng_seed=42)
    logger.info(f"Sig elevated bins (permutation): {perm_result['sig_elevated_bins']}")
    logger.info(f"Sig suppressed bins (permutation): {perm_result['sig_suppressed_bins']}")

    permutation_baseline_out = {
        "n_permutations": N_PERMUTATIONS,
        "bin_p5": perm_result["bin_p5"].tolist(),
        "bin_p95": perm_result["bin_p95"].tolist(),
        "sig_elevated_bins_permutation": perm_result["sig_elevated_bins"],
        "sig_suppressed_bins_permutation": perm_result["sig_suppressed_bins"],
        "perm_chi2_p95": float(np.percentile(perm_result["perm_chi2"], 95)),
    }

    # ------------------------------------------------------------------
    # Section 4: Sequential removal curves
    # ------------------------------------------------------------------
    logger.info("Building removal orders...")

    # Elevated-bin removal order: sort bins by influence descending, then usgs_id for determinism
    elevated_bin_order_sorted = sorted(ELEVATED_BINS, key=lambda b: -bin_influence_vals[b])
    logger.info(f"Elevated bins sorted by influence desc: {elevated_bin_order_sorted}")

    elevated_removal_order: list[str] = []
    for b in elevated_bin_order_sorted:
        bin_events = df[bin_indices == b][["usgs_id"]].sort_values("usgs_id")
        elevated_removal_order.extend(bin_events["usgs_id"].tolist())

    # Suppressed-bin removal order: sort bins by |negative influence| descending
    suppressed_bin_order_sorted = sorted(SUPPRESSED_BINS, key=lambda b: bin_influence_vals[b])  # most negative first
    logger.info(f"Suppressed bins sorted by |influence| desc: {suppressed_bin_order_sorted}")

    suppressed_removal_order: list[str] = []
    for b in suppressed_bin_order_sorted:
        bin_events = df[bin_indices == b][["usgs_id"]].sort_values("usgs_id")
        suppressed_removal_order.extend(bin_events["usgs_id"].tolist())

    logger.info(f"Elevated removal group size: {len(elevated_removal_order)}")
    logger.info(f"Suppressed removal group size: {len(suppressed_removal_order)}")

    # Run elevated-bin curve
    logger.info("Running elevated-bin sequential removal curve...")
    elevated_curve = sequential_removal_curve(df, elevated_removal_order, k=K)

    # Find signal persistence step for elevated
    elevated_persistence_step = None
    elevated_persistence_pct = None
    for step_data in elevated_curve:
        if step_data["p_chi2"] >= 0.05:
            elevated_persistence_step = step_data["step"]
            elevated_persistence_pct = step_data["pct_catalog_removed"]
            break

    if elevated_persistence_step is None:
        # Signal persists throughout
        elevated_persistence_step = elevated_curve[-1]["step"]
        elevated_persistence_pct = elevated_curve[-1]["pct_catalog_removed"]

    logger.info(f"Elevated-bin signal persistence: step={elevated_persistence_step}, pct={elevated_persistence_pct:.2f}%")

    # Run suppressed-bin curve
    logger.info("Running suppressed-bin sequential removal curve...")
    suppressed_curve = sequential_removal_curve(df, suppressed_removal_order, k=K)

    # Find signal persistence step for suppressed
    suppressed_persistence_step = None
    suppressed_persistence_pct = None
    for step_data in suppressed_curve:
        if step_data["p_chi2"] >= 0.05:
            suppressed_persistence_step = step_data["step"]
            suppressed_persistence_pct = step_data["pct_catalog_removed"]
            break

    if suppressed_persistence_step is None:
        suppressed_persistence_step = suppressed_curve[-1]["step"]
        suppressed_persistence_pct = suppressed_curve[-1]["pct_catalog_removed"]

    logger.info(f"Suppressed-bin signal persistence: step={suppressed_persistence_step}, pct={suppressed_persistence_pct:.2f}%")

    # G-K sequential baseline: GK aftershocks ordered by |delta_t_sec| ascending
    logger.info("Building G-K sequential baseline...")
    gk_after_sorted = gk_after.copy()
    gk_after_sorted["abs_delta_t"] = gk_after_sorted["delta_t_sec"].abs()
    gk_after_sorted = gk_after_sorted.sort_values("abs_delta_t")
    n_gk_baseline_elevated = min(len(gk_after_sorted), len(elevated_removal_order))
    n_gk_baseline_suppressed = min(len(gk_after_sorted), len(suppressed_removal_order))
    gk_removal_order_elevated = gk_after_sorted["usgs_id"].tolist()[:n_gk_baseline_elevated]
    gk_removal_order_suppressed = gk_after_sorted["usgs_id"].tolist()[:n_gk_baseline_suppressed]

    logger.info(f"Running G-K baseline curve (n={n_gk_baseline_elevated} for elevated, {n_gk_baseline_suppressed} for suppressed)...")
    gk_curve_elevated = sequential_removal_curve(df, gk_removal_order_elevated, k=K)
    gk_curve_suppressed = sequential_removal_curve(df, gk_removal_order_suppressed, k=K)

    # Random removal baseline — vectorized fast computation
    # Rather than running full sequential_removal_curve for each draw,
    # we compute chi2 at each step using vectorized bin-count updates.
    def fast_random_baseline(
        bin_idx_arr: np.ndarray,
        n_remove: int,
        n_draws: int,
        rng_seed: int,
        k: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute random removal chi2 curves efficiently.

        Returns mean, p10, p90 arrays of length (n_remove + 1).
        """
        n_total_local = len(bin_idx_arr)
        rng_local = np.random.default_rng(rng_seed)
        chi2_matrix = np.zeros((n_draws, n_remove + 1), dtype=float)

        # Initial chi2 (step 0) — same for all draws
        init_counts = np.bincount(bin_idx_arr, minlength=k).astype(float)
        init_expected = n_total_local / k
        init_chi2 = float(np.sum((init_counts - init_expected) ** 2 / init_expected))
        chi2_matrix[:, 0] = init_chi2

        for draw in range(n_draws):
            # Draw random removal indices (positions in bin_idx_arr)
            rand_positions = rng_local.choice(n_total_local, size=n_remove, replace=False)
            # Get bin assignments for those positions
            rand_bins = bin_idx_arr[rand_positions]

            # Running bin counts
            bc = init_counts.copy()
            for step_i, b_rm in enumerate(rand_bins):
                bc[b_rm] -= 1
                n_rem = n_total_local - (step_i + 1)
                if n_rem > 0:
                    exp_rem = n_rem / k
                    chi2_matrix[draw, step_i + 1] = float(np.sum((bc - exp_rem) ** 2 / exp_rem))
                else:
                    chi2_matrix[draw, step_i + 1] = 0.0

        return (
            chi2_matrix.mean(axis=0),
            np.percentile(chi2_matrix, 10, axis=0),
            np.percentile(chi2_matrix, 90, axis=0),
        )

    logger.info(f"Building random removal baseline (n_draws={N_RANDOM_BASELINES}) for elevated...")
    n_elevate_remove = len(elevated_removal_order)
    rand_mean_elev_arr, rand_p10_elev_arr, rand_p90_elev_arr = fast_random_baseline(
        bin_idx_arr=bin_indices,
        n_remove=n_elevate_remove,
        n_draws=N_RANDOM_BASELINES,
        rng_seed=42,
        k=K,
    )
    rand_mean_elev = rand_mean_elev_arr.tolist()
    rand_p10_elev = rand_p10_elev_arr.tolist()
    rand_p90_elev = rand_p90_elev_arr.tolist()

    logger.info(f"Building random removal baseline (n_draws={N_RANDOM_BASELINES}) for suppressed...")
    n_suppress_remove = len(suppressed_removal_order)
    rand_mean_supp_arr, rand_p10_supp_arr, rand_p90_supp_arr = fast_random_baseline(
        bin_idx_arr=bin_indices,
        n_remove=n_suppress_remove,
        n_draws=N_RANDOM_BASELINES,
        rng_seed=43,
        k=K,
    )
    rand_mean_supp = rand_mean_supp_arr.tolist()
    rand_p10_supp = rand_p10_supp_arr.tolist()
    rand_p90_supp = rand_p90_supp_arr.tolist()

    # Trim curves to remove step 0 duplication if needed and convert to JSON-serializable
    def curve_to_json(curve: list[dict]) -> list[dict]:
        out = []
        for step_data in curve:
            d = dict(step_data)
            d["elevated_bin_zscores"] = {str(k): v for k, v in d["elevated_bin_zscores"].items()}
            d["suppressed_bin_zscores"] = {str(k): v for k, v in d["suppressed_bin_zscores"].items()}
            out.append(d)
        return out

    degradation_curves = {
        "elevated_removal": {
            "n_events_in_group": len(elevated_removal_order),
            "signal_persistence_step": elevated_persistence_step,
            "signal_persistence_pct_catalog": float(elevated_persistence_pct) if elevated_persistence_pct is not None else None,
            "curve": curve_to_json(elevated_curve),
        },
        "suppressed_removal": {
            "n_events_in_group": len(suppressed_removal_order),
            "signal_persistence_step": suppressed_persistence_step,
            "signal_persistence_pct_catalog": float(suppressed_persistence_pct) if suppressed_persistence_pct is not None else None,
            "curve": curve_to_json(suppressed_curve),
        },
        "gk_sequential_baseline_elevated": {
            "n_events": n_gk_baseline_elevated,
            "curve": curve_to_json(gk_curve_elevated),
        },
        "gk_sequential_baseline_suppressed": {
            "n_events": n_gk_baseline_suppressed,
            "curve": curve_to_json(gk_curve_suppressed),
        },
        "random_baseline_elevated": {
            "n_draws": N_RANDOM_BASELINES,
            "curve_mean_chi2": rand_mean_elev,
            "curve_p10_chi2": rand_p10_elev,
            "curve_p90_chi2": rand_p90_elev,
        },
        "random_baseline_suppressed": {
            "n_draws": N_RANDOM_BASELINES,
            "curve_mean_chi2": rand_mean_supp,
            "curve_p10_chi2": rand_p10_supp,
            "curve_p90_chi2": rand_p90_supp,
        },
    }

    # ------------------------------------------------------------------
    # Section 5: Representativeness test
    # ------------------------------------------------------------------
    logger.info("Running representativeness tests...")
    repr_results = representativeness_test(df, top_n=TOP_N_REPRESENTATIVENESS)
    logger.info(f"Top-100 positive representative of signal stratum: {repr_results['top_positive_100']['representative_of_signal_stratum']}")
    logger.info(f"Top-100 negative representative of signal stratum: {repr_results['top_negative_100']['representative_of_signal_stratum']}")

    # ------------------------------------------------------------------
    # Section 6: Summary flags
    # ------------------------------------------------------------------
    # Signal diffuse: persistence requires removing >5% of catalog
    signal_diffuse = (elevated_persistence_pct is not None and elevated_persistence_pct > 5.0)

    # Symmetric degradation: check if elevated and suppressed z-scores converge at similar rates
    # Compare persistence percentages — if within 5% of each other, call it symmetric
    elev_pct = elevated_persistence_pct if elevated_persistence_pct is not None else 100.0
    supp_pct = suppressed_persistence_pct if suppressed_persistence_pct is not None else 100.0
    degradation_symmetric = abs(elev_pct - supp_pct) < 5.0

    summary = {
        "signal_diffuse": signal_diffuse,
        "elevated_persistence_pct": float(elev_pct),
        "suppressed_persistence_pct": float(supp_pct),
        "degradation_symmetric": degradation_symmetric,
        "top100_elevated_representative_of_signal_stratum": repr_results["top_positive_100"]["representative_of_signal_stratum"],
        "top100_suppressed_representative_of_signal_stratum": repr_results["top_negative_100"]["representative_of_signal_stratum"],
    }

    # ------------------------------------------------------------------
    # Serialize results JSON
    # ------------------------------------------------------------------
    results = {
        "case": "A3.A3",
        "title": "Phase-Concentration Audit",
        "parameters": {
            "n_catalog": 9210,
            "k": K,
            "julian_year_secs": JULIAN_YEAR_SECS,
            "elevated_bins": ELEVATED_BINS,
            "suppressed_bins": SUPPRESSED_BINS,
            "n_permutations": N_PERMUTATIONS,
            "n_random_baselines": N_RANDOM_BASELINES,
            "top_n_representativeness": TOP_N_REPRESENTATIVENESS,
            "rng_seed": 42,
        },
        "sequence_roles": {
            "n_isolated_mainshock": n_isolated,
            "n_mainshock_with_sequence": n_mainshock_seq,
            "n_aftershock": n_aftershock,
        },
        "bin_summary": bin_summary,
        "permutation_baseline": permutation_baseline_out,
        "degradation_curves": degradation_curves,
        "representativeness": repr_results,
        "summary": summary,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results written to {OUTPUT_PATH}")
    logger.info("=== A3.A3 analysis complete ===")


if __name__ == "__main__":
    main()
