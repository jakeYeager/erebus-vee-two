"""
Case A3.B1: Visualization Script

Generates all four publication-quality PNG figures from the A3.B1 results JSON:
  1. case-a3-b1-trajectory.png      — Chi-square + Rayleigh + mean phase (raw catalog)
  2. case-a3-b1-catalog-comparison.png — Chi-square trajectory overlay, all 4 catalogs
  3. case-a3-b1-interval-heatmap.png   — Per-window interval z-score heatmap, all catalogs
  4. case-a3-b1-phase-stability.png    — Circular mean-phase stability plot (raw catalog)
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a3-b1-results.json"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DPI = 300

# A1b interval centers (phase fraction)
INTERVAL_CENTERS = {
    "Interval 1": 0.208,
    "Interval 2": 0.646,
    "Interval 3": 0.896,
}

CATALOG_STYLES: dict[str, dict] = {
    "raw":        {"color": "steelblue",  "linestyle": "-",    "label": "Raw"},
    "gk":         {"color": "red",        "linestyle": "--",   "label": "G-K"},
    "reasenberg": {"color": "green",      "linestyle": ":",    "label": "Reasenberg"},
    "a1b":        {"color": "purple",     "linestyle": "-.",   "label": "A1b"},
}

DECADE_COLORS: dict[str, str] = {
    "1950s": "blue",
    "1960s": "steelblue",
    "1970s": "orange",
    "1980s": "green",
    "1990s": "purple",
    "2000s": "red",
    "2010s": "black",
}


def decade_label(window_start: int) -> str:
    """Return decade label string for a given window start year."""
    decade = (window_start // 10) * 10
    return f"{decade}s"


def mark_1970s_band(ax: plt.Axes, center_years: list[float]) -> None:
    """
    Shade 1970s windows (start years 1970–1979 → center years 1975–1984) on ax.

    Parameters
    ----------
    ax : matplotlib Axes
        Axes to annotate.
    center_years : list of float
        Window center years for all 62 windows.
    """
    in_band = [cy for cy in center_years if 1975 <= cy <= 1984]
    if in_band:
        ax.axvspan(min(in_band) - 0.5, max(in_band) + 0.5, color="lightyellow", zorder=0, alpha=0.9)


# ---------------------------------------------------------------------------
# Figure 1 — Chi-square + Rayleigh trajectory, raw catalog
# ---------------------------------------------------------------------------

def plot_trajectory(results: dict) -> None:
    """
    3-row stacked subplot: chi-square p (primary), Rayleigh p (secondary),
    mean phase — all for the raw catalog.

    Parameters
    ----------
    results : dict
        Full results dictionary loaded from the JSON output.
    """
    raw_windows = results["catalogs"]["raw"]["windows"]
    center_years = [w["window_start"] + 5 for w in raw_windows]
    p_chi2 = [w["p_chi2_k24"] for w in raw_windows]
    p_ray = [w["p_rayleigh"] for w in raw_windows]
    mean_ph = [w["mean_phase"] for w in raw_windows]

    bonferroni = results["meta"]["bonferroni_threshold"]

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    fig.suptitle("Raw Catalog — Chi-Square Primary (k=24)", fontsize=14, fontweight="bold", y=0.98)

    # Row 1: Chi-square p-value
    ax0 = axes[0]
    mark_1970s_band(ax0, center_years)
    ax0.semilogy(center_years, p_chi2, color="steelblue", linewidth=1.5)
    ax0.axhline(0.05, color="gray", linestyle="--", linewidth=1.0, label="p = 0.05")
    ax0.axhline(bonferroni, color="orange", linestyle="--", linewidth=1.0,
                label=f"Bonferroni p ≈ {bonferroni:.5f}")
    ax0.set_ylim(1e-12, 1.0)
    ax0.set_ylabel("p-value (log scale)", fontsize=9)
    ax0.set_title("Chi-square p (k=24) — primary", fontsize=10, loc="left")
    ax0.legend(fontsize=8, loc="upper right")
    ax0.grid(True, which="both", linestyle=":", alpha=0.4)

    # Row 2: Rayleigh p-value
    ax1 = axes[1]
    mark_1970s_band(ax1, center_years)
    ax1.semilogy(center_years, p_ray, color="red", linewidth=1.5)
    ax1.axhline(0.05, color="gray", linestyle="--", linewidth=1.0, label="p = 0.05")
    ax1.set_ylim(0.001, 1.0)
    ax1.set_ylabel("p-value (log scale)", fontsize=9)
    ax1.set_title("Rayleigh p — secondary", fontsize=10, loc="left")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(True, which="both", linestyle=":", alpha=0.4)

    # Row 3: Mean phase fraction
    ax2 = axes[2]
    mark_1970s_band(ax2, center_years)
    ax2.plot(center_years, mean_ph, color="black", linewidth=1.2, marker="o", markersize=3)
    for label, center in INTERVAL_CENTERS.items():
        ax2.axhline(center, color="gray", linestyle="--", linewidth=0.8, alpha=0.6, label=label)
    ax2.set_ylim(0.0, 1.0)
    ax2.set_ylabel("Mean phase (0–1)", fontsize=9)
    ax2.set_title("Mean phase fraction", fontsize=10, loc="left")
    ax2.legend(fontsize=7, loc="upper right", ncol=3)
    ax2.grid(True, linestyle=":", alpha=0.4)

    ax2.set_xlabel("Window center year", fontsize=10)
    ax2.set_xticks(range(1955, 2025, 5))
    ax2.tick_params(axis="x", labelsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out_path = OUTPUT_DIR / "case-a3-b1-trajectory.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 2 — Multi-catalog chi-square comparison
# ---------------------------------------------------------------------------

def plot_catalog_comparison(results: dict) -> None:
    """
    Single panel with chi-square p-value trajectories for all 4 catalogs overlaid.

    Parameters
    ----------
    results : dict
        Full results dictionary loaded from the JSON output.
    """
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_title("Chi-Square p-Value Trajectory by Catalog (k=24)", fontsize=13, fontweight="bold")

    # Determine 1970s band from raw windows
    raw_windows = results["catalogs"]["raw"]["windows"]
    center_years_raw = [w["window_start"] + 5 for w in raw_windows]
    mark_1970s_band(ax, center_years_raw)

    for cat_key, style in CATALOG_STYLES.items():
        windows = results["catalogs"][cat_key]["windows"]
        cx = [w["window_start"] + 5 for w in windows]
        py = [w["p_chi2_k24"] for w in windows]
        ax.semilogy(cx, py, color=style["color"], linestyle=style["linestyle"],
                    linewidth=1.5, label=style["label"])

    ax.axhline(0.05, color="gray", linestyle="--", linewidth=1.0, label="p = 0.05")
    ax.set_ylim(1e-12, 1.0)
    ax.set_ylabel("Chi-square p-value (log scale)", fontsize=10)
    ax.set_xlabel("Window center year", fontsize=10)
    ax.set_xticks(range(1955, 2025, 5))
    ax.tick_params(axis="x", labelsize=8)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, which="both", linestyle=":", alpha=0.4)

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b1-catalog-comparison.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 3 — Interval-level elevation heatmap
# ---------------------------------------------------------------------------

def plot_interval_heatmap(results: dict) -> None:
    """
    4-panel stacked heatmap (one per catalog). Each panel shows a 3-row heatmap
    (one row per A1b interval) of per-window interval z-scores.

    Parameters
    ----------
    results : dict
        Full results dictionary loaded from the JSON output.
    """
    cat_keys = list(CATALOG_STYLES.keys())
    n_cats = len(cat_keys)
    n_windows = 62
    interval_keys = ["interval_1", "interval_2", "interval_3"]
    interval_labels = ["Interval 1\n(Mar equinox)", "Interval 2\n(mid-Aug)", "Interval 3\n(late-Nov)"]

    # Gather center years from raw catalog
    raw_windows = results["catalogs"]["raw"]["windows"]
    center_years = [w["window_start"] + 5 for w in raw_windows]

    fig, axes = plt.subplots(n_cats, 1, figsize=(16, 10), sharex=True)
    fig.suptitle("Interval-Level Z-Score Heatmap by Catalog", fontsize=13, fontweight="bold", y=0.99)

    # Determine global z-score max for consistent colormap scaling
    all_z: list[float] = []
    for cat_key in cat_keys:
        for w in results["catalogs"][cat_key]["windows"]:
            for ik in interval_keys:
                all_z.append(w[f"{ik}_z"])
    vmax = max(max(all_z), 1.0)

    # White-to-red colormap (white=0/min, red=max)
    cmap = mcolors.LinearSegmentedColormap.from_list("white_red", ["white", "red"])
    norm = mcolors.Normalize(vmin=0.0, vmax=vmax)

    for ax_idx, cat_key in enumerate(cat_keys):
        ax = axes[ax_idx]
        windows = results["catalogs"][cat_key]["windows"]
        stationarity = results["catalogs"][cat_key]["stationarity"]

        # Build 3 x n_windows z-score matrix (clip negative z to 0 for colormap)
        z_matrix = np.zeros((3, n_windows))
        for wi, w in enumerate(windows):
            for ii, ik in enumerate(interval_keys):
                z_matrix[ii, wi] = max(0.0, w[f"{ik}_z"])

        im = ax.imshow(
            z_matrix,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            extent=[center_years[0] - 0.5, center_years[-1] + 0.5, -0.5, 2.5],
            origin="lower",
        )

        # Contour lines at z=1.96 and z=3.0
        x_fine = np.array(center_years, dtype=float)
        for ii in range(3):
            for z_thresh, lw in [(1.96, 1.0), (3.0, 1.5)]:
                # Draw horizontal band boundaries as line overlays
                # Use a filled step-like indicator above threshold
                passes = z_matrix[ii, :] >= z_thresh
                for wi_start in range(n_windows):
                    if passes[wi_start]:
                        # Draw a thin border at the cell boundaries
                        x0 = center_years[wi_start] - 0.5
                        x1 = center_years[wi_start] + 0.5
                        y0 = ii - 0.5
                        y1 = ii + 0.5
                        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0],
                                color="black" if z_thresh < 3.0 else "navy",
                                linewidth=lw * 0.4, alpha=0.5)

        # 1970s band
        in_band = [cy for cy in center_years if 1975 <= cy <= 1984]
        if in_band:
            ax.axvspan(min(in_band) - 0.5, max(in_band) + 0.5, color="lightyellow", zorder=0, alpha=0.7)

        # y-axis labels
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(interval_labels, fontsize=7)

        # Title with interval classifications
        i1c = stationarity["interval_1_classification"]
        i2c = stationarity["interval_2_classification"]
        i3c = stationarity["interval_3_classification"]
        style_label = CATALOG_STYLES[cat_key]["label"]
        ax.set_title(
            f"{style_label} — Int1: {i1c}, Int2: {i2c}, Int3: {i3c}",
            fontsize=9, loc="left",
        )
        ax.set_xlim(center_years[0] - 0.5, center_years[-1] + 0.5)

        # Colorbar on the right
        cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
        cbar.set_label("z-score", fontsize=7)
        cbar.ax.axhline(1.96, color="black", linewidth=0.8)
        cbar.ax.axhline(3.0, color="navy", linewidth=0.8)
        cbar.ax.tick_params(labelsize=7)

    axes[-1].set_xlabel("Window center year", fontsize=10)
    axes[-1].set_xticks(range(1955, 2025, 5))
    axes[-1].tick_params(axis="x", labelsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    out_path = OUTPUT_DIR / "case-a3-b1-interval-heatmap.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 4 — Phase stability circular plot, raw catalog
# ---------------------------------------------------------------------------

def plot_phase_stability(results: dict) -> None:
    """
    Polar plot of mean phase angle per 10-year window for the raw catalog.
    Points colored by decade; radius proportional to Rayleigh R.
    Gray wedges mark the three A1b baseline intervals.

    Parameters
    ----------
    results : dict
        Full results dictionary loaded from the JSON output.
    """
    raw_windows = results["catalogs"]["raw"]["windows"]

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw={"projection": "polar"})
    ax.set_title("Mean Phase Angle per 10-Year Window — Raw Catalog",
                 fontsize=12, fontweight="bold", pad=20)

    # Draw A1b interval wedges on the outer rim
    interval_ranges = [
        (2 * np.pi * 0.1667, 2 * np.pi * 0.2500),
        (2 * np.pi * 0.6250, 2 * np.pi * 0.6667),
        (2 * np.pi * 0.8750, 2 * np.pi * 0.9167),
    ]
    for theta_start, theta_end in interval_ranges:
        theta_arc = np.linspace(theta_start, theta_end, 60)
        r_inner = 0.90
        r_outer = 1.00
        ax.fill_between(theta_arc, r_inner, r_outer, color="lightgray", alpha=0.8, zorder=1)
        ax.plot(theta_arc, [r_outer] * len(theta_arc), color="darkgray", linewidth=0.8)

    # Plot window points
    decade_handles: dict[str, mpatches.Patch] = {}
    for w in raw_windows:
        theta = 2.0 * np.pi * w["mean_phase"]
        r = w["rayleigh_R"]
        dlabel = decade_label(w["window_start"])
        color = DECADE_COLORS.get(dlabel, "black")
        ax.scatter(theta, r, color=color, s=30, zorder=3, alpha=0.85)
        if dlabel not in decade_handles:
            decade_handles[dlabel] = mpatches.Patch(color=color, label=dlabel)

    # Sorted legend by decade
    sorted_handles = [decade_handles[d] for d in sorted(decade_handles.keys())]
    ax.legend(handles=sorted_handles, loc="upper left", bbox_to_anchor=(-0.18, 1.12),
              fontsize=8, title="Decade", title_fontsize=8)

    ax.set_rmax(1.0)
    ax.set_rticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.set_rlabel_position(90)
    ax.tick_params(labelsize=8)

    # Angle labels as month-of-year fractions
    ax.set_thetagrids(
        [0, 45, 90, 135, 180, 225, 270, 315],
        labels=["0.0", "0.125", "0.25", "0.375", "0.5", "0.625", "0.75", "0.875"],
        fontsize=7,
    )
    ax.set_xlabel("Rayleigh R (radius) · Mean phase (angle)", labelpad=15, fontsize=9)

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b1-phase-stability.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Load results JSON and render all four figures."""
    print(f"Loading results from {RESULTS_PATH}")
    with open(RESULTS_PATH, encoding="utf-8") as fh:
        results = json.load(fh)

    plot_trajectory(results)
    plot_catalog_comparison(results)
    plot_interval_heatmap(results)
    plot_phase_stability(results)
    print("All figures saved.")


if __name__ == "__main__":
    main()
