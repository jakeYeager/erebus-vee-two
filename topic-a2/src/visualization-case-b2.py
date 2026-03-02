"""Case B2: Ocean vs. Continent Location — Visualization Script.

Generates seven PNG figures:
  1. case-b2-global-map-gshhg.png  — global scatter map, GSHHG classification
  2. case-b2-global-map-ne.png     — global scatter map, Natural Earth classification
  3. case-b2-global-map-pb2002.png — global scatter map, PB2002 classification
  4. case-b2-regional-maps.png     — 4-row × 3-col regional zoom grid
  5. case-b2-binplots-gshhg.png    — 3-panel bin distributions (GSHHG, k=24)
  6. case-b2-binplots-ne.png       — 3-panel bin distributions (NE, k=24)
  7. case-b2-binplots-pb2002.png   — 3-panel bin distributions (PB2002, k=24)

cartopy/basemap/geopandas are unavailable in this environment. Global and
regional maps use lat/lon scatter with a simple continental outline drawn from
a built-in world-land boundary path derived from matplotlib's own test data.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
DATA_DIR = BASE_DIR.parent / "data" / "iscgem"
PLATE_LOC_DIR = DATA_DIR / "plate-location"
OUTPUT_DIR = BASE_DIR / "output"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("visualization-case-b2")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DPI = 300
SOLAR_YEAR_SECS: float = 31_557_600.0
K = 24

# Color scheme: oceanic=blue, transitional=green, continental=red
CLASS_COLORS: dict[str, str] = {
    "oceanic":      "#1f77b4",   # blue
    "transitional": "#2ca02c",   # green
    "continental":  "#d62728",   # red
}
CLASS_ORDER = ["oceanic", "transitional", "continental"]

# A1b baseline elevated phase intervals
A1B_INTERVALS: list[tuple[float, float]] = [
    (0.1875, 0.25),
    (0.625,  0.656),
    (0.875,  0.917),
]

# Solar calendar reference positions (fraction of year)
EQUINOX_SPRING  = 0.19
SOLSTICE_SUMMER = 0.44
EQUINOX_AUTUMN  = 0.69
SOLSTICE_WINTER = 0.94

# Regional map extents: (name, lat_min, lat_max, lon_min, lon_max)
REGIONS: list[tuple[str, float, float, float, float]] = [
    ("Philippines", 5,   22,  115, 130),
    ("Japan",       28,  47,  128, 148),
    ("Chile",       -58, -15, -80, -60),
    ("Java",        -12, -4,  100, 115),
]

# Classification methods
METHODS: list[tuple[str, str, str]] = [
    ("gshhg",  "ocean_class_gshhg",  "GSHHG"),
    ("ne",     "ocean_class_ne",     "Natural Earth"),
    ("pb2002", "ocean_class_pb2002", "PB2002"),
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_data() -> tuple[pd.DataFrame, dict]:
    """Load merged catalog and results JSON.

    Returns:
        Tuple of (merged DataFrame with all ocean_class_* columns and solar_phase,
                  results dict from JSON).
    """
    raw_path = DATA_DIR / "iscgem_global_6-9_1950-2021.csv"
    raw = pd.read_csv(raw_path)

    for key, filename in [
        ("gshhg",  "ocean_class_gshhg_global.csv"),
        ("ne",     "ocean_class_ne_global.csv"),
        ("pb2002", "ocean_class_pb2002_global.csv"),
    ]:
        cls_df = pd.read_csv(PLATE_LOC_DIR / filename).rename(columns={
            "ocean_class":      f"ocean_class_{key}",
            "dist_to_coast_km": f"dist_km_{key}",
        })
        raw = raw.merge(
            cls_df[["usgs_id", f"ocean_class_{key}", f"dist_km_{key}"]],
            on="usgs_id", how="left",
        )

    raw["solar_phase"] = (raw["solar_secs"] / SOLAR_YEAR_SECS) % 1.0

    results_path = OUTPUT_DIR / "case-b2-results.json"
    with open(results_path) as fh:
        results = json.load(fh)

    logger.info("Data loaded: n=%d", len(raw))
    return raw, results


# ---------------------------------------------------------------------------
# Simple world coastline (approximate continental outline via scatter fallback)
# ---------------------------------------------------------------------------
def draw_world_background(ax: plt.Axes) -> None:
    """Draw a light gray world background rectangle for a global map.

    Since cartopy/geopandas/basemap are unavailable, draws a simple
    ocean-colored background rectangle and sets axis limits.

    Args:
        ax: Matplotlib axes to configure.
    """
    # Light ocean background
    ax.set_facecolor("#d6eaf8")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude", fontsize=8)
    ax.set_ylabel("Latitude", fontsize=8)
    ax.tick_params(labelsize=7)
    # Grid lines
    ax.grid(color="white", linewidth=0.4, alpha=0.6)


def draw_region_background(
    ax: plt.Axes,
    lat_min: float, lat_max: float,
    lon_min: float, lon_max: float,
) -> None:
    """Set axis limits and background for a regional zoom map.

    Args:
        ax: Matplotlib axes.
        lat_min: Southern latitude bound.
        lat_max: Northern latitude bound.
        lon_min: Western longitude bound.
        lon_max: Eastern longitude bound.
    """
    ax.set_facecolor("#d6eaf8")
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.grid(color="white", linewidth=0.5, alpha=0.7)
    ax.tick_params(labelsize=6)


# ---------------------------------------------------------------------------
# Figures 1–3: Global maps
# ---------------------------------------------------------------------------
def plot_global_map(
    df: pd.DataFrame,
    method_key: str,
    col_name: str,
    method_label: str,
    filename: str,
) -> None:
    """Generate a global scatter map for one classification method.

    Events are plotted as scatter points colored by class. Point size is
    proportional to usgs_mag: size = (usgs_mag - 5.5) ** 2 * 2.

    Args:
        df: Merged DataFrame.
        method_key: e.g. "gshhg".
        col_name: Column name for ocean_class, e.g. "ocean_class_gshhg".
        method_label: Display label for the map title.
        filename: Output PNG filename.
    """
    fig, ax = plt.subplots(figsize=(14, 7))
    draw_world_background(ax)

    ax.set_title(
        f"ISC-GEM Events — {method_label} Classification\n"
        f"(M≥6.0, 1950–2021, n={len(df):,})",
        fontsize=12, fontweight="bold",
    )

    legend_handles = []
    for class_label in CLASS_ORDER:
        subset = df[df[col_name] == class_label]
        n_class = len(subset)
        color = CLASS_COLORS[class_label]
        sizes = (subset["usgs_mag"] - 5.5) ** 2 * 2
        ax.scatter(
            subset["longitude"], subset["latitude"],
            s=sizes, c=color, alpha=0.4, linewidths=0,
            rasterized=True,
            zorder=2,
        )
        legend_handles.append(mpatches.Patch(
            color=color,
            label=f"{class_label.capitalize()} (n={n_class:,})",
        ))

    ax.legend(
        handles=legend_handles,
        loc="lower left", fontsize=8,
        framealpha=0.9,
    )

    out_path = OUTPUT_DIR / filename
    plt.tight_layout()
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Global map saved: %s", out_path)


def plot_all_global_maps(df: pd.DataFrame) -> None:
    """Generate Figures 1–3: global maps for GSHHG, NE, and PB2002.

    Args:
        df: Merged DataFrame with all classification columns.
    """
    for method_key, col_name, label in METHODS:
        plot_global_map(
            df, method_key, col_name, label,
            f"case-b2-global-map-{method_key}.png",
        )


# ---------------------------------------------------------------------------
# Figure 4: Regional maps (4 rows × 3 cols)
# ---------------------------------------------------------------------------
def plot_regional_maps(df: pd.DataFrame) -> None:
    """Generate Figure 4: 4-row × 3-column regional zoom grid.

    Rows: Philippines, Japan, Chile, Java.
    Columns: GSHHG, Natural Earth, PB2002.

    Args:
        df: Merged DataFrame with all classification columns.
    """
    n_rows = len(REGIONS)
    n_cols = len(METHODS)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(15, 20),
    )
    fig.suptitle(
        "Case B2: Regional Event Maps by Classification Method",
        fontsize=13, fontweight="bold",
    )

    for col_idx, (method_key, col_name, method_label) in enumerate(METHODS):
        # Column header
        axes[0][col_idx].set_title(method_label, fontsize=10, fontweight="bold", pad=8)

        for row_idx, (region_name, lat_min, lat_max, lon_min, lon_max) in enumerate(REGIONS):
            ax = axes[row_idx][col_idx]
            draw_region_background(ax, lat_min, lat_max, lon_min, lon_max)

            # Filter to region
            mask = (
                (df["latitude"] >= lat_min) & (df["latitude"] <= lat_max) &
                (df["longitude"] >= lon_min) & (df["longitude"] <= lon_max)
            )
            df_region = df[mask]

            for class_label in CLASS_ORDER:
                subset = df_region[df_region[col_name] == class_label]
                if len(subset) == 0:
                    continue
                color = CLASS_COLORS[class_label]
                sizes = (subset["usgs_mag"] - 5.5) ** 2 * 2
                ax.scatter(
                    subset["longitude"], subset["latitude"],
                    s=sizes, c=color, alpha=0.55, linewidths=0.2,
                    edgecolors="white",
                    zorder=2,
                )

            # Row label on leftmost column
            if col_idx == 0:
                ax.set_ylabel(
                    f"{region_name}\n[lat {lat_min}–{lat_max}]",
                    fontsize=8, fontweight="bold",
                )
            else:
                ax.set_ylabel("")

            # Count annotation
            n_region = len(df_region)
            ax.text(
                0.02, 0.97, f"n={n_region}",
                transform=ax.transAxes,
                fontsize=7, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
            )

    # Shared legend at bottom
    legend_handles = [
        mpatches.Patch(color=CLASS_COLORS[cl], label=cl.capitalize())
        for cl in CLASS_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center", ncol=3, fontsize=9,
        bbox_to_anchor=(0.5, 0.0),
    )

    plt.tight_layout(rect=[0, 0.025, 1, 0.98])
    out_path = OUTPUT_DIR / "case-b2-regional-maps.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Regional maps saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figures 5–7: Bin distribution plots (3-panel per method)
# ---------------------------------------------------------------------------
def _compute_threshold(n: int, k: int) -> float:
    """Compute 1-SD elevated-bin threshold (expected + sqrt(expected)).

    Args:
        n: Total event count.
        k: Number of bins.

    Returns:
        Float threshold value.
    """
    expected = n / k
    return expected + np.sqrt(expected)


def plot_single_class_binplot(
    ax: plt.Axes,
    bin_counts: list[int],
    n: int,
    k: int,
    class_label: str,
    chi2: float,
    p_chi2: float,
    cramer_v: float,
) -> None:
    """Draw a horizontal bar chart of phase bin distribution.

    Same format as A4/B1 bin plots: steelblue bars, dashed expected-count line,
    1-SD threshold, A1b baseline interval gray bands.

    Args:
        ax: Matplotlib axes.
        bin_counts: Observed bin counts list (length k).
        n: Total events in this class.
        k: Number of bins.
        class_label: e.g. "Oceanic".
        chi2: Chi-square statistic.
        p_chi2: Chi-square p-value.
        cramer_v: Cramér's V effect size.
    """
    counts = np.array(bin_counts, dtype=float)
    expected = n / k
    threshold = _compute_threshold(n, k)

    bin_edges = np.linspace(0, 1, k + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = 1.0 / k

    colors = ["orange" if c > threshold else "steelblue" for c in counts]

    ax.barh(
        bin_centers, counts,
        height=bin_width * 0.85,
        color=colors, edgecolor="white", linewidth=0.4,
    )

    # Expected count line
    ax.axvline(
        expected, color="black", linestyle="--", linewidth=1.2,
        label=f"Expected ({expected:.1f})",
    )
    # 1-SD threshold line
    ax.axvline(
        threshold, color="gray", linestyle=":", linewidth=1.0,
        label=f"1-SD threshold ({threshold:.1f})",
    )

    # A1b baseline interval shaded bands
    for ps, pe in A1B_INTERVALS:
        ax.axhspan(ps, pe, color="lightgray", alpha=0.55, zorder=0)

    # Calendar reference lines
    for pos, lbl in [
        (EQUINOX_SPRING,  "Spr Eq"),
        (SOLSTICE_SUMMER, "Sum Sol"),
        (EQUINOX_AUTUMN,  "Aut Eq"),
        (SOLSTICE_WINTER, "Win Sol"),
    ]:
        ax.axhline(pos, color="dimgray", linestyle=":", linewidth=0.7, alpha=0.7)

    # Axes formatting
    ax.set_title(f"{class_label.capitalize()} (n={n:,})", fontsize=10, fontweight="bold")
    ax.set_xlabel("Event Count", fontsize=9)
    ax.set_ylabel("Solar Phase (fraction of year)", fontsize=9)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0, 1.05, 0.25))
    ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=8)

    # Statistical annotation
    p_str = f"{p_chi2:.2e}" if p_chi2 < 0.001 else f"{p_chi2:.4f}"
    annot = (
        f"n={n:,}\n"
        f"\u03c7\u00b2={chi2:.2f}\n"
        f"p={p_str}\n"
        f"V={cramer_v:.4f}"
    )
    ax.text(
        0.97, 0.97, annot,
        transform=ax.transAxes,
        ha="right", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85),
    )

    # Legend
    legend_handles = [
        mpatches.Patch(color="steelblue", label="Normal bin"),
        mpatches.Patch(color="orange", label=">1 SD elevated"),
        mpatches.Patch(color="lightgray", alpha=0.55, label="A1b baseline interval"),
        mlines.Line2D([0], [0], color="black", linestyle="--", linewidth=1.2,
                      label="Expected count"),
        mlines.Line2D([0], [0], color="gray", linestyle=":", linewidth=1.0,
                      label="1-SD threshold"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="lower right")


def plot_binplot_figure(
    results: dict,
    method_key: str,
    method_label: str,
    filename: str,
) -> None:
    """Generate a 3-panel bin distribution figure for one classification method.

    Panels are ordered: oceanic (top), transitional (middle), continental (bottom).

    Args:
        results: Loaded results JSON dict.
        method_key: e.g. "gshhg".
        method_label: Display label for figure title.
        filename: Output PNG filename.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 9))
    fig.suptitle(
        f"Case B2: Solar Phase Bin Distributions by Location Class — {method_label} (k={K})",
        fontsize=13, fontweight="bold",
    )

    for col_idx, class_label in enumerate(CLASS_ORDER):
        ax = axes[col_idx]
        stats = results["class_stats"][method_key][class_label]
        k_stats = stats[f"k{K}"]

        plot_single_class_binplot(
            ax=ax,
            bin_counts=k_stats["bin_counts"],
            n=k_stats["n"],
            k=K,
            class_label=class_label,
            chi2=k_stats["chi2"],
            p_chi2=k_stats["p_chi2"],
            cramer_v=k_stats["cramer_v"],
        )

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = OUTPUT_DIR / filename
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Binplot figure saved: %s", out_path)


def plot_all_binplots(results: dict) -> None:
    """Generate Figures 5–7: bin distribution plots for all three methods.

    Args:
        results: Loaded results JSON dict.
    """
    for method_key, _, method_label in METHODS:
        plot_binplot_figure(
            results,
            method_key=method_key,
            method_label=method_label,
            filename=f"case-b2-binplots-{method_key}.png",
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load data and generate all seven PNG figures."""
    logger.info("=== Case B2 Visualization ===")

    df, results = load_data()

    # Figures 1–3: Global maps
    logger.info("Generating global maps...")
    plot_all_global_maps(df)

    # Figure 4: Regional maps
    logger.info("Generating regional maps...")
    plot_regional_maps(df)

    # Figures 5–7: Bin distribution plots
    logger.info("Generating bin distribution plots...")
    plot_all_binplots(results)

    logger.info("All figures written to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
