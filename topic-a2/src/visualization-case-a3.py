"""Case A3: Magnitude Stratification — Visualization Script.

Generates two PNG figures:
  1. case-a3-binplots.png      — 2×2 panel, one per magnitude band at k=24
  2. case-a3-effect-trend.png  — Dual-axis Cramér's V and Rayleigh R vs band
"""

import json
import logging
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a3-results.json"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("visualization-case-a3")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DPI = 300

# A1b baseline elevated phase intervals (from Adhoc A1b)
A1B_INTERVALS = [
    (0.1875, 0.25),   # Interval 1: ~March equinox
    (0.625, 0.656),   # Interval 2: ~mid-August
    (0.875, 0.917),   # Interval 3: ~mid-November
]

# Solar calendar reference positions (fraction of Julian year, starting Jan 1)
EQUINOX_SPRING  = 0.19   # ~March 11
SOLSTICE_SUMMER = 0.44   # ~June 11
EQUINOX_AUTUMN  = 0.69   # ~September 10
SOLSTICE_WINTER = 0.94   # ~December 10

CALENDAR_REFS = [
    (EQUINOX_SPRING,  "Spr Eq"),
    (SOLSTICE_SUMMER, "Sum Sol"),
    (EQUINOX_AUTUMN,  "Aut Eq"),
    (SOLSTICE_WINTER, "Win Sol"),
]

BAND_LABELS = ["M6.0-6.4", "M6.5-6.9", "M7.0-7.4", "M7.5+"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_results() -> dict:
    """Load the case A3 results JSON.

    Returns:
        Parsed results dictionary.
    """
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


def compute_threshold(n: int, k: int) -> float:
    """Compute 1-SD elevated-bin threshold (E + sqrt(E)).

    Args:
        n: Total events.
        k: Number of bins.

    Returns:
        Threshold value.
    """
    expected = n / k
    return expected + np.sqrt(expected)


# ---------------------------------------------------------------------------
# Figure 1 — Four-panel bin distributions (2×2, k=24)
# ---------------------------------------------------------------------------
def plot_band_panel(
    ax: plt.Axes,
    band_label: str,
    k24_data: dict,
    n: int,
) -> None:
    """Draw a single magnitude-band horizontal bar chart panel.

    Args:
        ax: Matplotlib axes.
        band_label: Band label string.
        k24_data: k=24 statistics dict for this band.
        n: Band event count.
    """
    k = 24
    counts = np.array(k24_data["bin_counts"], dtype=float)
    expected = n / k
    threshold = compute_threshold(n, k)

    bin_edges = np.linspace(0, 1, k + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = 1.0 / k

    colors = ["orange" if c > threshold else "steelblue" for c in counts]

    ax.barh(
        bin_centers, counts,
        height=bin_width * 0.85,
        color=colors, edgecolor="white", linewidth=0.4,
    )

    # Expected count dashed line
    ax.axvline(
        expected, color="black", linestyle="--", linewidth=1.0,
        label=f"Expected ({expected:.1f})",
    )
    # 1-SD threshold dotted line
    ax.axvline(
        threshold, color="gray", linestyle=":", linewidth=0.9,
        label=f"1-SD threshold ({threshold:.1f})",
    )

    # A1b baseline interval shaded bands
    for ps, pe in A1B_INTERVALS:
        ax.axhspan(ps, pe, color="lightgray", alpha=0.55, zorder=0)

    # Calendar reference lines
    for pos, lbl in CALENDAR_REFS:
        ax.axhline(pos, color="dimgray", linestyle=":", linewidth=0.6, alpha=0.7)
        ax.text(
            0.01, pos + 0.005, lbl,
            transform=ax.get_yaxis_transform(),
            fontsize=5.5, color="dimgray", va="bottom",
        )

    # Statistical annotation
    p_val = k24_data["p_chi2"]
    p_str = f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.4f}"
    annot = (
        f"n={n:,}\n"
        f"\u03c7\u00b2={k24_data['chi2']:.2f}\n"
        f"p={p_str}\n"
        f"V={k24_data['cramer_v']:.4f}"
    )
    ax.text(
        0.97, 0.97, annot,
        transform=ax.transAxes,
        ha="right", va="top", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
    )

    ax.set_title(band_label, fontsize=10, fontweight="bold")
    ax.set_xlabel("Event Count", fontsize=8)
    ax.set_ylabel("Solar Phase (0–1)", fontsize=8)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0, 1.05, 0.25))
    ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=7)
    ax.tick_params(axis="x", labelsize=7)


def plot_binplots(results: dict) -> None:
    """Generate Figure 1: 2×2 panel bin distributions at k=24.

    Args:
        results: Full results dict.
    """
    band_stats = results["band_stats"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 14))
    fig.suptitle(
        "Case A3: Solar Phase Distribution by Magnitude Band (k=24)",
        fontsize=13, fontweight="bold",
    )

    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for (row, col), band_label in zip(positions, BAND_LABELS):
        ax = axes[row][col]
        band_data = band_stats[band_label]
        n = band_data["n"]
        k24_data = band_data["k24"]
        plot_band_panel(ax, band_label, k24_data, n)

    # Shared legend
    legend_patches = [
        mpatches.Patch(color="steelblue", label="Normal bin"),
        mpatches.Patch(color="orange", label=">1 SD elevated"),
        mpatches.Patch(color="lightgray", alpha=0.55, label="A1b baseline interval"),
        plt.Line2D([0], [0], color="black", linestyle="--", linewidth=1.0,
                   label="Expected count"),
        plt.Line2D([0], [0], color="gray", linestyle=":", linewidth=0.9,
                   label="1-SD threshold"),
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center", ncol=5, fontsize=7.5,
        bbox_to_anchor=(0.5, 0.0),
    )

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    out_path = OUTPUT_DIR / "case-a3-binplots.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 1 (binplots) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 2 — Effect-size trend plot (dual axis)
# ---------------------------------------------------------------------------
def plot_effect_trend(results: dict) -> None:
    """Generate Figure 2: Dual-axis Cramér's V and Rayleigh R vs magnitude band.

    Args:
        results: Full results dict.
    """
    trend = results["trend_analysis"]
    band_stats = results["band_stats"]

    x_positions = [1, 2, 3, 4]
    x_labels = BAND_LABELS

    v_vals = [trend["cramer_v_by_band"][lbl] for lbl in BAND_LABELS]
    r_vals = [trend["rayleigh_R_by_band"][lbl] for lbl in BAND_LABELS]

    # Bootstrap CIs for Cramér's V
    ci_lowers = [band_stats[lbl]["k24"]["cramer_v_ci95_lower"] for lbl in BAND_LABELS]
    ci_uppers = [band_stats[lbl]["k24"]["cramer_v_ci95_upper"] for lbl in BAND_LABELS]

    # Significance at k=24
    sig_bands = set(trend["significant_bands"])
    v_sig_markers = ["o" if lbl in sig_bands else "^" for lbl in BAND_LABELS]
    r_sig_markers = ["o" if lbl in sig_bands else "^" for lbl in BAND_LABELS]

    fig, ax1 = plt.subplots(figsize=(8, 5))

    # Cramér's V on left axis (steelblue)
    ax1.set_xlabel("Magnitude Band", fontsize=10)
    ax1.set_ylabel("Cramér's V", color="steelblue", fontsize=10)
    ax1.tick_params(axis="y", labelcolor="steelblue")

    # Draw line
    ax1.plot(x_positions, v_vals, color="steelblue", linewidth=1.5, zorder=3)

    # Error bars: clamp to zero to handle edge case where point estimate
    # falls near the boundary of the bootstrap CI distribution
    v_err_low = [max(0.0, v - lo) for v, lo in zip(v_vals, ci_lowers)]
    v_err_high = [max(0.0, hi - v) for v, hi in zip(v_vals, ci_uppers)]
    ax1.errorbar(
        x_positions, v_vals,
        yerr=[v_err_low, v_err_high],
        fmt="none", color="steelblue", capsize=4, linewidth=1.0, zorder=4,
    )

    # Filled vs open circles for significance
    for xi, vi, marker, lbl in zip(x_positions, v_vals, v_sig_markers, BAND_LABELS):
        if lbl in sig_bands:
            ax1.scatter(xi, vi, color="steelblue", s=50, zorder=5,
                        marker="o", label="p<0.05 (V)")
        else:
            ax1.scatter(xi, vi, color="steelblue", s=50, zorder=5,
                        marker="o", facecolors="none", edgecolors="steelblue",
                        linewidths=1.5)

    # Rayleigh R on right axis (orange)
    ax2 = ax1.twinx()
    ax2.set_ylabel("Rayleigh R", color="darkorange", fontsize=10)
    ax2.tick_params(axis="y", labelcolor="darkorange")

    ax2.plot(x_positions, r_vals, color="darkorange", linewidth=1.5, zorder=3)
    for xi, ri, lbl in zip(x_positions, r_vals, BAND_LABELS):
        if lbl in sig_bands:
            ax2.scatter(xi, ri, color="darkorange", s=50, zorder=5,
                        marker="o")
        else:
            ax2.scatter(xi, ri, color="darkorange", s=50, zorder=5,
                        marker="o", facecolors="none", edgecolors="darkorange",
                        linewidths=1.5)

    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(x_labels, fontsize=9)
    ax1.set_xlim(0.5, 4.5)

    # Spearman annotation
    rho_v = trend["spearman_rho_cramer_v"]
    p_v = trend["spearman_p_cramer_v"]
    p_str = f"{p_v:.2e}" if p_v < 0.001 else f"{p_v:.4f}"
    trend_cls = trend["trend_classification"]
    annot = (
        f"Cramér's V Spearman \u03c1={rho_v:.3f}, p={p_str}\n"
        f"Trend: {trend_cls}"
    )
    ax1.text(
        0.03, 0.97, annot,
        transform=ax1.transAxes,
        ha="left", va="top", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
    )

    # Legend
    legend_elements = [
        plt.Line2D([0], [0], color="steelblue", linewidth=1.5,
                   label="Cramér's V (left axis)"),
        plt.Line2D([0], [0], color="darkorange", linewidth=1.5,
                   label="Rayleigh R (right axis)"),
        plt.scatter([], [], color="dimgray", marker="o", s=40,
                    label="p < 0.05 (filled)"),
        plt.scatter([], [], color="dimgray", marker="o", s=40,
                    facecolors="none", edgecolors="dimgray", linewidths=1.5,
                    label="p ≥ 0.05 (open)"),
    ]
    ax1.legend(handles=legend_elements, fontsize=8, loc="upper right")

    fig.suptitle(
        "Case A3: Effect Size vs Magnitude Band (k=24)",
        fontsize=12, fontweight="bold",
    )

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-effect-trend.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 2 (effect trend) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load results and generate all figures."""
    logger.info("Loading results from %s", RESULTS_PATH)
    results = load_results()

    plot_binplots(results)
    plot_effect_trend(results)

    logger.info("All figures written to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
