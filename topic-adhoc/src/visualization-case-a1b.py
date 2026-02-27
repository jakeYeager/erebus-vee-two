"""
Case A1b: Visualization script

Reads case-a1b-results.json and produces three PNG images:
  1. case-a1b-phase-coherence.png  — elevated bin positions across k=16, 24, 32
  2. case-a1b-temporal-magnitude.png — IEI distribution + magnitude distribution
  3. case-a1b-spatial.png — geographic scatter, NN distances, boundary proximity
"""

import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Path conventions
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent          # topic-adhoc/
PROJECT_ROOT = BASE_DIR.parent                             # erebus-vee-two/
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a1b-results.json"
DATA_PATH = PROJECT_ROOT / "data" / "global-sets" / "iscgem_global_events.csv"

# Phase normalization constants (matching A1 and A1b analysis)
DAY_SECS = 86400


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def load_results() -> dict:
    """Load the A1b results JSON."""
    with open(RESULTS_PATH) as f:
        return json.load(f)


def is_leap_year(year: int) -> bool:
    """Return True if year is a Gregorian leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def year_length_seconds(year: int) -> int:
    """Return calendar-year length in seconds."""
    return 366 * DAY_SECS if is_leap_year(year) else 365 * DAY_SECS


def compute_solar_phase(df: pd.DataFrame) -> pd.Series:
    """Compute solar phase in [0, 1) matching A1b analysis."""
    year_lengths = df["solaration_year"].apply(year_length_seconds).astype(float)
    phase = df["solar_secs"] / year_lengths
    over_mask = phase >= 1.0
    if over_mask.any():
        phase = phase.copy()
        phase[over_mask] = 0.999999
    return phase


# ---------------------------------------------------------------------------
# Image 1: Phase coherence
# ---------------------------------------------------------------------------

def plot_phase_coherence(results: dict) -> None:
    """Three-row panel showing elevated bin positions at k=16, 24, 32,
    plus a fourth row for the combined elevated phase intervals."""
    pc = results["phase_coherence"]
    bin_counts = [16, 24, 32]
    combined_intervals = pc["combined_elevated_intervals"]

    fig, axes = plt.subplots(4, 1, figsize=(12, 7))
    fig.suptitle("Solar Phase Elevated Bins by Bin Count", fontsize=14, fontweight="bold")

    # Month annotations (approximate positions)
    month_positions = {
        "Jan": 0.0, "Apr": 0.25, "Jul": 0.5, "Oct": 0.75
    }

    for row_idx, k in enumerate(bin_counts):
        ax = axes[row_idx]
        key = f"k{k}"
        top3_bins = pc[key]["top3_bins"]
        bin_width = 1.0 / k

        for i in range(k):
            lo = i * bin_width
            hi = (i + 1) * bin_width
            color = "steelblue" if i in top3_bins else "lightgray"
            rect = mpatches.FancyBboxPatch(
                (lo, 0.1), bin_width - 0.002, 0.8,
                boxstyle="square,pad=0",
                facecolor=color, edgecolor="white", linewidth=0.5
            )
            ax.add_patch(rect)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(f"k={k}", fontsize=10, rotation=0, ha="right", va="center")

        if row_idx < len(bin_counts) - 1:
            ax.set_xticks([])
        else:
            ax.set_xticks([v for v in month_positions.values()])
            ax.set_xticklabels([k_label for k_label in month_positions.keys()], fontsize=9)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    # Row 4: combined elevated intervals
    ax = axes[3]
    for lo, hi in combined_intervals:
        rect = mpatches.FancyBboxPatch(
            (lo, 0.1), hi - lo - 0.001, 0.8,
            boxstyle="square,pad=0",
            facecolor="red", edgecolor="darkred", linewidth=0.8, alpha=0.8
        )
        ax.add_patch(rect)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel("Combined", fontsize=10, rotation=0, ha="right", va="center")
    ax.set_xticks(list(month_positions.values()))
    ax.set_xticklabels(list(month_positions.keys()), fontsize=9)
    ax.set_xlabel("Solar Phase (approximate calendar month)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Legend
    elevated_patch = mpatches.Patch(color="steelblue", label="Top-3 elevated bin")
    other_patch = mpatches.Patch(color="lightgray", label="Other bins")
    combined_patch = mpatches.Patch(color="red", alpha=0.8, label="Combined elevated intervals (>=2 of 3 k)")
    fig.legend(handles=[elevated_patch, other_patch, combined_patch],
               loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout(rect=[0, 0.04, 1, 0.97])

    out_path = OUTPUT_DIR / "case-a1b-phase-coherence.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Saved %s", out_path)


# ---------------------------------------------------------------------------
# Image 2: Temporal and magnitude
# ---------------------------------------------------------------------------

def plot_temporal_magnitude(results: dict) -> None:
    """Two-panel: IEI distribution (left) and magnitude distribution (right)."""
    df = pd.read_csv(DATA_PATH)
    df["solar_phase"] = compute_solar_phase(df)
    phases = df["solar_phase"].values
    n = len(df)

    # Reconstruct elevated-bin events
    combined_intervals = results["phase_coherence"]["combined_elevated_intervals"]
    elevated_mask = np.zeros(n, dtype=bool)
    for lo, hi in combined_intervals:
        elevated_mask |= (phases >= lo) & (phases < hi)
    elevated_df = df[elevated_mask].drop_duplicates(subset=["usgs_id"])

    # Parse event timestamps for IEI
    elevated_df = elevated_df.copy()
    elevated_df["event_dt"] = pd.to_datetime(elevated_df["event_at"], utc=True)
    elevated_df_sorted = elevated_df.sort_values("event_dt")
    diffs = elevated_df_sorted["event_dt"].diff().dropna()
    iei_days_elev = diffs.dt.total_seconds() / 86400.0

    # Generate a single representative random sample of equal size (fixed seed)
    rng = np.random.default_rng(42)
    n_elevated = len(elevated_df)
    df_sorted = df.copy()
    df_sorted["event_dt"] = pd.to_datetime(df_sorted["event_at"], utc=True)
    df_sorted = df_sorted.sort_values("event_dt")
    idx = rng.choice(len(df_sorted), size=n_elevated, replace=False)
    sample_times = np.sort(df_sorted["event_dt"].values[idx])
    iei_sample = np.diff(sample_times.astype(np.int64)) / 1e9 / 86400.0

    # Null CI from results
    null_ci = results["temporal"]["null_iei_median_ci"]
    elev_stats = results["temporal"]["elevated_iei_days"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Temporal and Magnitude Analysis: Elevated Bins vs Full Catalog",
                 fontsize=13, fontweight="bold")

    # Left panel: IEI distribution
    # Clip extreme values for display
    iei_clip = 1000.0
    iei_elev_clipped = np.clip(iei_days_elev, 0.01, iei_clip)
    iei_samp_clipped = np.clip(iei_sample, 0.01, iei_clip)

    bins = np.logspace(np.log10(0.01), np.log10(iei_clip), 50)
    ax1.hist(iei_samp_clipped, bins=bins, histtype="step", color="gray", linewidth=1.5,
             label=f"Random sample (n={n_elevated})", alpha=0.8)
    ax1.hist(iei_elev_clipped, bins=bins, histtype="step", color="steelblue", linewidth=1.5,
             label=f"Elevated-bin events (n={n_elevated})")
    ax1.set_xscale("log")
    ax1.set_xlabel("Inter-Event Interval (days, log scale)", fontsize=10)
    ax1.set_ylabel("Count", fontsize=10)
    ax1.set_title("Inter-Event Intervals:\nElevated Bins vs Random Baseline", fontsize=11)

    elev_median = elev_stats["median"]
    samp_median = float(np.median(iei_sample))
    ax1.axvline(elev_median, color="steelblue", linestyle="--", linewidth=1.5,
                label=f"Elevated median: {elev_median:.2f} d")
    ax1.axvline(samp_median, color="gray", linestyle="--", linewidth=1.5,
                label=f"Sample median: {samp_median:.2f} d")

    # Null CI band
    ax1.axvspan(null_ci["p2_5"], null_ci["p97_5"], alpha=0.15, color="gray",
                label=f"Null 95% CI [{null_ci['p2_5']:.1f}–{null_ci['p97_5']:.1f} d]")
    ax1.legend(fontsize=8, loc="upper right")

    # Right panel: Magnitude distribution
    mag = results["magnitude"]
    mag_bins = ["6.0-6.4", "6.5-6.9", "7.0-7.4", "7.5+"]
    n_elev = sum(mag["elevated"].values())
    n_full = sum(mag["full_catalog"].values())

    elev_pcts = [mag["elevated"][b] / n_elev * 100 for b in mag_bins]
    full_pcts = [mag["full_catalog"][b] / n_full * 100 for b in mag_bins]

    x = np.arange(len(mag_bins))
    width = 0.35
    ax2.bar(x - width / 2, full_pcts, width, color="lightgray", label=f"Full catalog (n={n_full})")
    ax2.bar(x + width / 2, elev_pcts, width, color="steelblue", label=f"Elevated-bin events (n={n_elev})")
    ax2.set_xticks(x)
    ax2.set_xticklabels(mag_bins, fontsize=9)
    ax2.set_xlabel("Magnitude Bin", fontsize=10)
    ax2.set_ylabel("% of Population", fontsize=10)
    ax2.set_title("Magnitude Distribution:\nElevated Bins vs Full Catalog", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.yaxis.grid(True, linestyle="--", alpha=0.5)
    ax2.set_axisbelow(True)

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a1b-temporal-magnitude.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Saved %s", out_path)


# ---------------------------------------------------------------------------
# Image 3: Spatial
# ---------------------------------------------------------------------------

def plot_spatial(results: dict) -> None:
    """Three-panel: geographic scatter, NN distances, boundary proximity."""
    df = pd.read_csv(DATA_PATH)
    df["solar_phase"] = compute_solar_phase(df)
    phases = df["solar_phase"].values
    n = len(df)

    # Reconstruct elevated-bin events
    combined_intervals = results["phase_coherence"]["combined_elevated_intervals"]
    elevated_mask = np.zeros(n, dtype=bool)
    for lo, hi in combined_intervals:
        elevated_mask |= (phases >= lo) & (phases < hi)
    elevated_df = df[elevated_mask].drop_duplicates(subset=["usgs_id"])

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Spatial Analysis: Elevated Bin Events", fontsize=13, fontweight="bold")

    # Left panel: geographic scatter
    ax1.scatter(df["longitude"], df["latitude"],
                c="lightgray", s=1, alpha=0.4, rasterized=True, label="Full catalog")
    ax1.scatter(elevated_df["longitude"], elevated_df["latitude"],
                c="steelblue", s=3, alpha=0.6, rasterized=True,
                label=f"Elevated-bin events (n={len(elevated_df)})")
    ax1.set_xlabel("Longitude", fontsize=10)
    ax1.set_ylabel("Latitude", fontsize=10)
    ax1.set_title("Elevated Bin Events:\nGeographic Distribution", fontsize=11)
    ax1.set_xlim(-180, 180)
    ax1.set_ylim(-90, 90)
    ax1.legend(fontsize=8, loc="lower left", markerscale=3)
    ax1.set_aspect("equal")

    # Center panel: NN distance distribution
    nn_stats = results["spatial"]["elevated_nn_km"]
    null_nn_ci = results["spatial"]["null_nn_median_ci"]

    # We need the actual NN distances — recompute from elevated_df
    # (stored in the analysis script's internal list but not serialized to JSON spec)
    # Use the analysis module's NN logic inline here
    from typing import TYPE_CHECKING
    elev_lat = elevated_df["latitude"].values
    elev_lon = elevated_df["longitude"].values
    nn_km_vals = _compute_nn_for_plot(elev_lat, elev_lon)

    bins_nn = np.linspace(0, 1000, 80)
    ax2.hist(nn_km_vals, bins=bins_nn, histtype="step", color="steelblue",
             linewidth=1.5, label=f"Elevated-bin NN (n={len(nn_km_vals)})")
    ax2.axvline(nn_stats["p50"], color="steelblue", linestyle="--", linewidth=1.5,
                label=f"Elevated p50: {nn_stats['p50']:.1f} km")
    ax2.axvspan(null_nn_ci["p2_5"], null_nn_ci["p97_5"], alpha=0.2, color="gray",
                label=f"Null 95% CI [{null_nn_ci['p2_5']:.0f}–{null_nn_ci['p97_5']:.0f} km]")
    ax2.axvline(null_nn_ci["p50"], color="gray", linestyle="--", linewidth=1.2,
                label=f"Null median: {null_nn_ci['p50']:.1f} km")
    ax2.axvline(49, color="red", linestyle="--", linewidth=1.5, alpha=0.8,
                label="G-K M6.0 spatial window (49 km)")
    ax2.set_xlabel("Nearest-Neighbor Distance (km)", fontsize=10)
    ax2.set_ylabel("Count", fontsize=10)
    ax2.set_title("Nearest-Neighbor Distance:\nElevated Bins vs Null", fontsize=11)
    ax2.legend(fontsize=8)
    ax2.set_xlim(0, 800)

    # Right panel: boundary proximity stacked horizontal bars
    bnd = results["spatial"]["boundary_proximity"]
    categories = ["Near Boundary\n(<=100 km)", "Transitional\n(100-300 km)", "Intraplate\n(>300 km)"]
    elev_vals = [
        bnd["elevated"]["near_boundary_pct"],
        bnd["elevated"]["transitional_pct"],
        bnd["elevated"]["intraplate_pct"],
    ]
    full_vals = [
        bnd["full_catalog"]["near_boundary_pct"],
        bnd["full_catalog"]["transitional_pct"],
        bnd["full_catalog"]["intraplate_pct"],
    ]

    # Steelblue shades for elevated, gray shades for full catalog
    elev_colors = ["#2378ae", "#5ba3cf", "#9ecae1"]  # dark to light steelblue
    full_colors = ["#525252", "#969696", "#d4d4d4"]   # dark to light gray

    y_positions = [0.7, 0.3]
    bar_height = 0.25

    # Full catalog bar (bottom row)
    left = 0
    for i, (val, color) in enumerate(zip(full_vals, full_colors)):
        ax3.barh(y_positions[1], val, height=bar_height, left=left,
                 color=color, edgecolor="white", linewidth=0.5)
        if val > 5:
            ax3.text(left + val / 2, y_positions[1], f"{val:.1f}%",
                     ha="center", va="center", fontsize=8, color="white" if i < 2 else "black")
        left += val

    # Elevated bar (top row)
    left = 0
    for i, (val, color) in enumerate(zip(elev_vals, elev_colors)):
        ax3.barh(y_positions[0], val, height=bar_height, left=left,
                 color=color, edgecolor="white", linewidth=0.5)
        if val > 5:
            ax3.text(left + val / 2, y_positions[0], f"{val:.1f}%",
                     ha="center", va="center", fontsize=8, color="white" if i < 2 else "black")
        left += val

    ax3.set_yticks(y_positions)
    ax3.set_yticklabels(["Elevated-bin", "Full catalog"], fontsize=10)
    ax3.set_xlabel("% of Population", fontsize=10)
    ax3.set_title("Boundary Proximity:\nElevated Bins vs Full Catalog", fontsize=11)
    ax3.set_xlim(0, 100)

    # Legend for proximity classes
    legend_handles = [
        mpatches.Patch(color=elev_colors[0], label="Near boundary (<=100 km)"),
        mpatches.Patch(color=elev_colors[1], label="Transitional (100-300 km)"),
        mpatches.Patch(color=elev_colors[2], label="Intraplate (>300 km)"),
    ]
    ax3.legend(handles=legend_handles, fontsize=8, loc="lower right")

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a1b-spatial.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Saved %s", out_path)


def _compute_nn_for_plot(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Compute nearest-neighbor distances (km) using chunked vectorized Haversine.

    Duplicate of analysis script logic for visualization use.
    """
    EARTH_RADIUS_KM = 6371.0
    n = len(lat)
    nn_dist = np.full(n, np.inf)
    chunk_size = 500

    lat_r = np.radians(lat)
    lon_r = np.radians(lon)

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        lat_chunk = lat_r[start:end, np.newaxis]
        lon_chunk = lon_r[start:end, np.newaxis]
        lat_all = lat_r[np.newaxis, :]
        lon_all = lon_r[np.newaxis, :]

        dlat = lat_all - lat_chunk
        dlon = lon_all - lon_chunk

        a = np.sin(dlat / 2) ** 2 + np.cos(lat_chunk) * np.cos(lat_all) * np.sin(dlon / 2) ** 2
        dist_mat = 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

        for ci, gi in enumerate(range(start, end)):
            dist_mat[ci, gi] = np.inf
        nn_dist[start:end] = dist_mat.min(axis=1)

    return nn_dist


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Generate all three visualization images for Case A1b."""
    logger.info("Loading results from %s", RESULTS_PATH)
    results = load_results()

    logger.info("Generating phase coherence plot...")
    plot_phase_coherence(results)

    logger.info("Generating temporal/magnitude plot...")
    plot_temporal_magnitude(results)

    logger.info("Generating spatial plot...")
    plot_spatial(results)

    logger.info("All visualization images written to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
