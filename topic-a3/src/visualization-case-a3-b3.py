"""
Case A3.B3: Visualization script

Generates 5 publication-quality figures:
1. case-a3-b3-threshold-sweep.png — chi2 p + Cramér's V + class n across 8 thresholds
2. case-a3-b3-global-map.png — Cartopy Robinson global map with SUB overlay
3. case-a3-b3-subduction-crosstab.png — cross-tab of ocean class × subduction proximity
4. case-a3-b3-binplots.png — bin distributions at baseline and key threshold
5. case-a3-b3-region-maps.png — 3-panel regional maps (Japan, Philippines, New Zealand)
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = BASE_DIR / "output"
LIB_DIR = BASE_DIR.parent / "lib"
STEPS_PATH = LIB_DIR / "PB2002_steps.dat"

RESULTS_PATH = OUTPUT_DIR / "case-a3-b3-results.json"
EVENTS_PATH  = OUTPUT_DIR / "case-a3-b3-events.pkl"
SUB_PTS_PATH = OUTPUT_DIR / "case-a3-b3-sub-points.npy"

CLASS_COLORS = {
    "oceanic":      "steelblue",
    "transitional": "orange",
    "continental":  "red",
}

INTERVAL_BINS = {
    "interval_1": [4, 5],
    "interval_2": [15],
    "interval_3": [21],
}
K_BINS = 24
JULIAN_YEAR_SECS = 31_557_600.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def classify_at_threshold(
    dist_series: pd.Series, t_outer: float, t_inner: float = 50.0
) -> pd.Series:
    """Classify events into oceanic/transitional/continental at given threshold."""
    def _label(d: float) -> str:
        if d <= t_inner:
            return "continental"
        elif d > t_outer:
            return "oceanic"
        else:
            return "transitional"
    return dist_series.map(_label)


def parse_sub_segments_for_plot(steps_path: Path) -> list[tuple[float, float, float, float]]:
    """
    Parse PB2002_steps.dat returning SUB/OCB segments as (lon1, lat1, lon2, lat2) tuples.
    """
    segments: list[tuple[float, float, float, float]] = []
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
                segments.append((lon1, lat1, lon2, lat2))
            except (ValueError, IndexError):
                continue
    return segments


def draw_sub_segments(ax: plt.Axes, segments: list, linewidth: float = 1.0, color: str = "black") -> None:
    """Draw SUB boundary segments on a cartopy axes."""
    for lon1, lat1, lon2, lat2 in segments:
        ax.plot(
            [lon1, lon2], [lat1, lat2],
            color=color, linewidth=linewidth,
            transform=ccrs.PlateCarree(), zorder=3,
        )


# ---------------------------------------------------------------------------
# Figure 1: Threshold sweep trajectory
# ---------------------------------------------------------------------------

def figure1_threshold_sweep(results: dict) -> None:
    """3-row stacked plot: chi2 p (log), Cramér's V, class n vs. T_outer."""
    sweep = results["threshold_sweep"]
    thresholds = [s["t_outer_km"] for s in sweep]
    classes = ["oceanic", "transitional", "continental"]
    colors = [CLASS_COLORS[c] for c in classes]
    labels = ["Oceanic", "Transitional", "Continental"]

    p_vals   = {c: [s[c]["p_chi2_k24"] for s in sweep] for c in classes}
    cv_vals  = {c: [s[c]["cramers_v"]   for s in sweep] for c in classes}
    n_vals   = {c: [s[c]["n"]           for s in sweep] for c in classes}

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.suptitle("Threshold Sensitivity — GSHHG Classification", fontsize=14, fontweight="bold")

    # Row 1: chi2 p-value (log scale)
    ax0 = axes[0]
    for c, col, lbl in zip(classes, colors, labels):
        ax0.plot(thresholds, p_vals[c], "o-", color=col, label=lbl, linewidth=1.8, markersize=5)
    ax0.axhline(0.05, color="gray", linestyle="--", linewidth=1.2, label="p=0.05")
    ax0.set_yscale("log")
    ax0.set_ylabel("Chi-square p-value (log)")
    ax0.legend(fontsize=9, loc="upper right")
    ax0.grid(True, which="both", alpha=0.3)

    # Vertical reference lines
    ax0.axvline(200, color="gray", linestyle="--", linewidth=1.0)

    # Mark oceanic significance crossing (if any)
    first_sig = results["summary"]["first_oceanic_significant_threshold_km"]
    if first_sig is not None:
        ax0.axvline(first_sig, color="steelblue", linestyle=":", linewidth=1.5, label=f"Oceanic sig. T={int(first_sig)} km")
        ax0.legend(fontsize=9, loc="upper right")

    # Row 2: Cramér's V
    ax1 = axes[1]
    for c, col, lbl in zip(classes, colors, labels):
        ax1.plot(thresholds, cv_vals[c], "o-", color=col, label=lbl, linewidth=1.8, markersize=5)
    ax1.set_ylabel("Cramér's V")
    ax1.grid(True, alpha=0.3)
    ax1.axvline(200, color="gray", linestyle="--", linewidth=1.0)
    if first_sig is not None:
        ax1.axvline(first_sig, color="steelblue", linestyle=":", linewidth=1.5)

    # Row 3: class n
    ax2 = axes[2]
    for c, col, lbl in zip(classes, colors, labels):
        ax2.plot(thresholds, n_vals[c], "o-", color=col, label=lbl, linewidth=1.8, markersize=5)
    ax2.set_ylabel("Class N")
    ax2.set_xlabel("Outer threshold (km from coast)")
    ax2.grid(True, alpha=0.3)
    ax2.axvline(200, color="gray", linestyle="--", linewidth=1.0, label="A2.B2 baseline")
    ax2.legend(fontsize=9, loc="upper left")
    if first_sig is not None:
        ax2.axvline(first_sig, color="steelblue", linestyle=":", linewidth=1.5)

    # Reverse x-axis (200 on left, 25 on right)
    ax2.set_xlim(max(thresholds) + 10, min(thresholds) - 10)
    ax2.set_xticks(thresholds)
    ax2.set_xticklabels([str(int(t)) for t in thresholds])

    # Annotate baseline
    ax0.text(200 + 2, ax0.get_ylim()[0] * 1.5, "A2.B2\nbaseline", fontsize=7, color="gray", va="bottom")

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b3-threshold-sweep.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Global map
# ---------------------------------------------------------------------------

def figure2_global_map(df: pd.DataFrame, sub_segs: list, results: dict) -> None:
    """Cartopy Robinson global map with GSHHG baseline classification and SUB overlay."""
    # Baseline classification at T=200
    classes = classify_at_threshold(df["dist_km_gshhg"], 200.0)
    n_counts = classes.value_counts()

    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()

    # Coastlines
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="gray", zorder=2)

    # Plot events by class
    for cls in ["continental", "transitional", "oceanic"]:
        mask = classes == cls
        sub_df = df[mask]
        sizes = ((sub_df["usgs_mag"] - 5.5) ** 2 * 1.5).clip(lower=0.5)
        ax.scatter(
            sub_df["longitude"].values,
            sub_df["latitude"].values,
            c=CLASS_COLORS[cls],
            s=sizes,
            alpha=0.4,
            linewidths=0,
            transform=ccrs.PlateCarree(),
            zorder=4,
            label=f"{cls.capitalize()} (n={n_counts.get(cls, 0):,})",
        )

    # SUB boundary overlay
    draw_sub_segments(ax, sub_segs, linewidth=1.0, color="black")

    # Legend
    handles, labels_leg = ax.get_legend_handles_labels()
    bnd_patch = mpatches.Patch(color="black", label="PB2002 subduction boundaries")
    handles.append(bnd_patch)
    labels_leg.append("Black lines = PB2002 subduction boundaries")
    ax.legend(handles=handles, labels=labels_leg, loc="lower left", fontsize=8, framealpha=0.8)

    ax.set_title("ISC-GEM Events — GSHHG Classification at T=200 km (A2.B2 baseline)", fontsize=12, fontweight="bold")

    out_path = OUTPUT_DIR / "case-a3-b3-global-map.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Subduction proximity cross-tabulation
# ---------------------------------------------------------------------------

def figure3_subduction_crosstab(results: dict) -> None:
    """2-row stacked: stacked bar (class × near/far sub) + oceanic pct_near_sub line."""
    sweep = results["threshold_sweep"]
    thresholds = [s["t_outer_km"] for s in sweep]
    classes = ["oceanic", "transitional", "continental"]

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle("Subduction Zone Proximity vs. Threshold Reclassification", fontsize=13, fontweight="bold")

    # Row 1: Stacked bars — for each class, near vs. far at each threshold
    x = np.arange(len(thresholds))
    width = 0.25
    offsets = [-0.25, 0.0, 0.25]

    for i, cls in enumerate(classes):
        n_near = [s[cls]["n_near_sub"] for s in sweep]
        n_total = [s[cls]["n"] for s in sweep]
        n_far  = [n_total[j] - n_near[j] for j in range(len(sweep))]

        base_color = CLASS_COLORS[cls]
        dark_color = mcolors.to_rgba(base_color, alpha=0.85)
        light_color = mcolors.to_rgba(base_color, alpha=0.35)

        ax0.bar(x + offsets[i], n_near, width, color=dark_color, label=f"{cls.capitalize()} near sub")
        ax0.bar(x + offsets[i], n_far,  width, bottom=n_near, color=light_color, label=f"{cls.capitalize()} far sub")

    ax0.set_ylabel("N events")
    ax0.legend(fontsize=7, ncol=3, loc="upper right")
    ax0.grid(True, axis="y", alpha=0.3)
    ax0.set_xticks(x)
    ax0.set_xticklabels([str(int(t)) for t in thresholds])

    # Row 2: Oceanic pct_near_sub vs. threshold
    oce_pct = [s["oceanic"]["pct_near_sub"] for s in sweep]
    catalog_frac = results["subduction_proximity"]["pct_near_subduction"]

    ax1.plot(x, oce_pct, "o-", color="steelblue", linewidth=2, markersize=6, label="Oceanic pct near subduction")
    ax1.axhline(catalog_frac, color="gray", linestyle="--", linewidth=1.2, label=f"Full catalog ({catalog_frac:.2%})")
    ax1.set_ylabel("Fraction near subduction")
    ax1.set_xlabel("Outer threshold (km from coast, reversed)")
    ax1.set_ylim(0, 1)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(x)
    ax1.set_xticklabels([str(int(t)) for t in thresholds])

    # Reverse x-axis
    ax0.set_xlim(x[-1] + 0.5, x[0] - 0.5)

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b3-subduction-crosstab.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 4: Bin distributions at baseline and key threshold
# ---------------------------------------------------------------------------

def figure4_binplots(results: dict) -> None:
    """2×3 grid of bin distribution plots at baseline and key threshold."""
    sweep = results["threshold_sweep"]
    classes = ["oceanic", "transitional", "continental"]
    col_labels = ["Oceanic", "Transitional", "Continental"]

    # Determine key threshold
    first_sig = results["summary"]["first_oceanic_significant_threshold_km"]
    key_t = first_sig if first_sig is not None else 100.0

    # Get the two steps
    baseline_step = next(s for s in sweep if s["t_outer_km"] == 200.0)
    key_step = next((s for s in sweep if s["t_outer_km"] == key_t), sweep[-1])

    rows = [baseline_step, key_step]
    row_labels = [
        f"T=200 km (A2.B2 baseline)",
        f"T={int(key_t)} km",
    ]

    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    fig.suptitle("Solar Phase Bin Distributions — Baseline vs. Key Threshold", fontsize=13, fontweight="bold")

    bin_positions = np.arange(K_BINS)

    for row_idx, (step, row_lbl) in enumerate(zip(rows, row_labels)):
        for col_idx, (cls, col_lbl) in enumerate(zip(classes, col_labels)):
            ax = axes[row_idx, col_idx]
            stats = step[cls]
            n = stats["n"]
            bin_counts = stats.get("bin_counts", [])

            if not bin_counts or n < K_BINS:
                ax.text(0.5, 0.5, f"n={n}\n(insufficient)", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{col_lbl}")
                continue

            expected = n / K_BINS
            sd_expected = np.sqrt(expected)

            # Gray band for A1b intervals
            for iname, bins in INTERVAL_BINS.items():
                for b in bins:
                    ax.axvspan(b - 0.5, b + 0.5, color="lightgray", alpha=0.5, zorder=0)

            # Expected line and SD threshold
            ax.axhline(expected, color="black", linestyle="--", linewidth=1.0, label="Expected")
            ax.axhline(expected + sd_expected, color="darkgray", linestyle=":", linewidth=0.8, label="±1 SD")
            ax.axhline(expected - sd_expected, color="darkgray", linestyle=":", linewidth=0.8)

            # Horizontal bar chart
            ax.barh(bin_positions, bin_counts, color="steelblue", alpha=0.8, height=0.85)

            # Annotations
            p_val = stats.get("p_chi2_k24")
            chi2_val = stats.get("chi2_k24")
            cv = stats.get("cramers_v")
            ann_txt = f"n={n:,}\nχ²={chi2_val:.2f}\np={p_val:.4f}\nV={cv:.4f}"
            ax.text(0.97, 0.97, ann_txt, transform=ax.transAxes,
                    ha="right", va="top", fontsize=7.5,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))

            if row_idx == 0:
                ax.set_title(f"{col_lbl}", fontsize=11, fontweight="bold")
            ax.set_xlabel("Count" if row_idx == 1 else "")
            ax.set_ylabel(f"Bin (k=24)\n{row_lbl}" if col_idx == 0 else "")
            ax.set_ylim(-0.5, K_BINS - 0.5)
            ax.set_yticks([0, 6, 12, 18, 23])
            ax.grid(True, axis="x", alpha=0.3)

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b3-binplots.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 5: Regional detail maps
# ---------------------------------------------------------------------------

def figure5_region_maps(df: pd.DataFrame, sub_segs: list) -> None:
    """1×3 Cartopy regional inset maps: Japan, Philippines, New Zealand."""
    # Baseline classification at T=200
    classes_series = classify_at_threshold(df["dist_km_gshhg"], 200.0)

    regions = [
        {
            "name": "Japan",
            "extent": [128, 148, 30, 46],
        },
        {
            "name": "Philippines",
            "extent": [114, 130, 4, 22],
        },
        {
            "name": "New Zealand",
            "extent": [164, 180, -48, -32],
        },
    ]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(18, 7),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    fig.suptitle(
        "Regional Coastline Resolution — GSHHG Classification at T=200 km (A2.B2 baseline)",
        fontsize=13, fontweight="bold",
    )

    # 50m resolution features
    land_50m = cfeature.NaturalEarthFeature("physical", "land", "50m", facecolor="whitesmoke")
    coast_50m = cfeature.NaturalEarthFeature("physical", "coastline", "50m", edgecolor="dimgray", facecolor="none")

    for ax, region in zip(axes, regions):
        ext = region["extent"]
        lon_min, lon_max, lat_min, lat_max = ext

        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

        # Background features
        ax.add_feature(land_50m, zorder=0)
        ax.add_feature(coast_50m, linewidth=0.6, zorder=2)

        # Filter events in extent
        mask_ext = (
            (df["longitude"] >= lon_min) & (df["longitude"] <= lon_max) &
            (df["latitude"] >= lat_min) & (df["latitude"] <= lat_max)
        )
        df_ext = df[mask_ext].copy()
        cls_ext = classes_series[mask_ext]

        n_total_ext = len(df_ext)
        n_trans = (cls_ext == "transitional").sum()
        n_oce   = (cls_ext == "oceanic").sum()

        # Plot events
        for cls in ["continental", "transitional", "oceanic"]:
            m = cls_ext == cls
            sub_evt = df_ext[m]
            sizes = ((sub_evt["usgs_mag"] - 5.5) ** 2 * 1.5).clip(lower=0.5)
            ax.scatter(
                sub_evt["longitude"].values,
                sub_evt["latitude"].values,
                c=CLASS_COLORS[cls],
                s=sizes,
                alpha=0.5,
                linewidths=0,
                transform=ccrs.PlateCarree(),
                zorder=4,
            )

        # SUB overlay — filter to region extent (with small buffer)
        buf = 2.0
        for lon1, lat1, lon2, lat2 in sub_segs:
            if (
                max(lon1, lon2) >= lon_min - buf and min(lon1, lon2) <= lon_max + buf and
                max(lat1, lat2) >= lat_min - buf and min(lat1, lat2) <= lat_max + buf
            ):
                ax.plot(
                    [lon1, lon2], [lat1, lat2],
                    color="black", linewidth=1.2,
                    transform=ccrs.PlateCarree(), zorder=3,
                )

        # Gridlines
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle=":")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": 7}
        gl.ylabel_style = {"size": 7}

        # Annotation
        ann = f"{region['name']}\nn={n_total_ext}\ntransitional={n_trans}, oceanic={n_oce}"
        ax.text(0.02, 0.98, ann, transform=ax.transAxes,
                fontsize=8, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Shared legend below panels
    legend_patches = [
        mpatches.Patch(color=CLASS_COLORS["oceanic"], label="Oceanic"),
        mpatches.Patch(color=CLASS_COLORS["transitional"], label="Transitional"),
        mpatches.Patch(color=CLASS_COLORS["continental"], label="Continental"),
        mpatches.Patch(color="black", label="PB2002 subduction boundaries"),
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center",
        ncol=4,
        fontsize=9,
        framealpha=0.9,
        bbox_to_anchor=(0.5, -0.02),
    )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    out_path = OUTPUT_DIR / "case-a3-b3-region-maps.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Generate all A3.B3 figures."""
    logger.info("Loading results JSON...")
    with open(RESULTS_PATH) as fh:
        results = json.load(fh)

    logger.info("Loading event dataframe...")
    df = pd.read_pickle(EVENTS_PATH)

    logger.info("Loading subduction boundary points...")
    sub_pts = np.load(str(SUB_PTS_PATH))

    logger.info("Parsing subduction segments for plotting...")
    sub_segs = parse_sub_segments_for_plot(STEPS_PATH)
    logger.info(f"  {len(sub_segs)} SUB+OCB segments loaded for plotting")

    logger.info("Generating Figure 1: threshold sweep trajectory...")
    figure1_threshold_sweep(results)

    logger.info("Generating Figure 2: global map...")
    figure2_global_map(df, sub_segs, results)

    logger.info("Generating Figure 3: subduction crosstab...")
    figure3_subduction_crosstab(results)

    logger.info("Generating Figure 4: bin plots...")
    figure4_binplots(results)

    logger.info("Generating Figure 5: regional maps...")
    figure5_region_maps(df, sub_segs)

    logger.info("All figures generated successfully.")


if __name__ == "__main__":
    main()
