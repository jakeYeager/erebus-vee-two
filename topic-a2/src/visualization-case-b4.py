"""Case B4: Depth Stratification — Visualization Script.

Generates three PNG figures:
  1. case-b4-binplots.png          — 2×2 panel, one per depth band at k=24
  2. case-b4-depth-trend.png       — Dual-axis Cramér's V and Rayleigh R vs depth band
  3. case-b4-deep-events-phase.png — Detailed phase distribution for deep (>300 km) events
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
RESULTS_PATH = OUTPUT_DIR / "case-b4-results.json"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("visualization-case-b4")

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

BAND_LABELS = [
    "shallow_0-20km",
    "midcrustal_20-70km",
    "intermediate_70-300km",
    "deep_300km+",
]

BAND_DISPLAY = {
    "shallow_0-20km":        "Shallow (0–20 km)",
    "midcrustal_20-70km":    "Mid-Crustal (20–70 km)",
    "intermediate_70-300km": "Intermediate (70–300 km)",
    "deep_300km+":           "Deep (>300 km)",
}

# Zhan & Shearer (2015) April–September phase range
ZHAN_PHASE_MIN = 0.23
ZHAN_PHASE_MAX = 0.67


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_results() -> dict:
    """Load the case B4 results JSON.

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
# Figure 1 — Multi-panel bin distributions (2×2, k=24)
# ---------------------------------------------------------------------------
def plot_band_panel(
    ax: plt.Axes,
    band_label: str,
    k24_data: dict,
    n: int,
    sufficient_n: bool,
) -> None:
    """Draw a single depth-band horizontal bar chart panel.

    Args:
        ax: Matplotlib axes.
        band_label: Internal band label string.
        k24_data: k=24 statistics dict for this band.
        n: Band event count.
        sufficient_n: Whether band meets minimum sample threshold.
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
    low_n_note = "\n(low n)" if not sufficient_n else ""
    annot = (
        f"n={n:,}\n"
        f"\u03c7\u00b2={k24_data['chi2']:.2f}\n"
        f"p={p_str}\n"
        f"V={k24_data['cramer_v']:.4f}"
        f"{low_n_note}"
    )
    ax.text(
        0.97, 0.97, annot,
        transform=ax.transAxes,
        ha="right", va="top", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
    )

    title = BAND_DISPLAY[band_label]
    if not sufficient_n:
        title += " (low n)"
    ax.set_title(title, fontsize=10, fontweight="bold")
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
        "Case B4: Solar Phase Distribution by Depth Band (k=24)",
        fontsize=13, fontweight="bold",
    )

    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for (row, col), band_label in zip(positions, BAND_LABELS):
        ax = axes[row][col]
        band_data = band_stats[band_label]
        n = band_data["n"]
        sufficient_n = band_data["sufficient_n"]
        k24_data = band_data["k24"]
        plot_band_panel(ax, band_label, k24_data, n, sufficient_n)

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
    out_path = OUTPUT_DIR / "case-b4-binplots.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 1 (binplots) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 2 — Depth trend plot (dual axis)
# ---------------------------------------------------------------------------
def plot_depth_trend(results: dict) -> None:
    """Generate Figure 2: Dual-axis Cramér's V and Rayleigh R vs depth band.

    Args:
        results: Full results dict.
    """
    trend = results["trend_analysis"]
    band_stats = results["band_stats"]

    x_positions = [1, 2, 3, 4]
    x_labels = [BAND_DISPLAY[lbl] for lbl in BAND_LABELS]

    v_vals = [trend["cramer_v_by_band"][lbl] for lbl in BAND_LABELS]
    r_vals = [band_stats[lbl]["k24"]["rayleigh_R"] for lbl in BAND_LABELS]

    # Significance at k=24
    sig_bands = set(trend["significant_bands"])
    sufficient_bands = {lbl for lbl in BAND_LABELS if band_stats[lbl]["sufficient_n"]}

    fig, ax1 = plt.subplots(figsize=(9, 5.5))

    # Cramér's V on left axis (steelblue)
    ax1.set_xlabel("Depth Band", fontsize=10)
    ax1.set_ylabel("Cramér's V", color="steelblue", fontsize=10)
    ax1.tick_params(axis="y", labelcolor="steelblue")

    ax1.plot(x_positions, v_vals, color="steelblue", linewidth=1.5, zorder=3)

    for xi, vi, lbl in zip(x_positions, v_vals, BAND_LABELS):
        if lbl not in sufficient_bands:
            # Insufficient n: 'x' marker
            ax1.scatter(xi, vi, color="steelblue", s=80, zorder=5,
                        marker="x", linewidths=2.0)
        elif lbl in sig_bands:
            # Significant: filled circle
            ax1.scatter(xi, vi, color="steelblue", s=60, zorder=5, marker="o")
        else:
            # Not significant: open circle
            ax1.scatter(xi, vi, color="steelblue", s=60, zorder=5,
                        marker="o", facecolors="none", edgecolors="steelblue",
                        linewidths=1.5)

    # Rayleigh R on right axis (orange)
    ax2 = ax1.twinx()
    ax2.set_ylabel("Rayleigh R", color="darkorange", fontsize=10)
    ax2.tick_params(axis="y", labelcolor="darkorange")

    ax2.plot(x_positions, r_vals, color="darkorange", linewidth=1.5, zorder=3)
    for xi, ri, lbl in zip(x_positions, r_vals, BAND_LABELS):
        if lbl not in sufficient_bands:
            ax2.scatter(xi, ri, color="darkorange", s=80, zorder=5,
                        marker="x", linewidths=2.0)
        elif lbl in sig_bands:
            ax2.scatter(xi, ri, color="darkorange", s=60, zorder=5, marker="o")
        else:
            ax2.scatter(xi, ri, color="darkorange", s=60, zorder=5,
                        marker="o", facecolors="none", edgecolors="darkorange",
                        linewidths=1.5)

    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(x_labels, fontsize=8, rotation=10, ha="right")
    ax1.set_xlim(0.5, 4.5)

    # Spearman annotation
    rho = trend["spearman_rho"]
    p_val = trend["spearman_p"]
    p_str = f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.4f}"
    mono_cls = trend["monotonicity_classification"]
    annot = (
        f"Spearman \u03c1={rho:.3f}, p={p_str}\n"
        f"Monotonicity: {mono_cls}"
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
                    label="p \u2265 0.05 (open)"),
        plt.scatter([], [], color="dimgray", marker="x", s=60, linewidths=2.0,
                    label="low n (< 100)"),
    ]
    ax1.legend(handles=legend_elements, fontsize=8, loc="upper right")

    fig.suptitle(
        "Case B4: Effect Size vs Depth Band (k=24)",
        fontsize=12, fontweight="bold",
    )

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-b4-depth-trend.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 2 (depth trend) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 3 — Deep events phase detail with Zhan & Shearer highlight
# ---------------------------------------------------------------------------
def plot_deep_events_phase(results: dict) -> None:
    """Generate Figure 3: Detailed phase distribution for deep (>300 km) events.

    Shades the April–September range (phase 0.23–0.67) in light orange
    and annotates the Zhan & Shearer (2015) finding.

    Args:
        results: Full results dict.
    """
    band_stats = results["band_stats"]
    trend = results["trend_analysis"]
    deep_data = band_stats["deep_300km+"]
    n = deep_data["n"]
    sufficient_n = deep_data["sufficient_n"]
    k24_data = deep_data["k24"]

    k = 24
    counts = np.array(k24_data["bin_counts"], dtype=float)
    expected = n / k
    threshold = compute_threshold(n, k)

    bin_edges = np.linspace(0, 1, k + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = 1.0 / k

    fig, ax = plt.subplots(figsize=(9, 6))

    # April–September shade (Zhan & Shearer range) — light orange background
    ax.axhspan(ZHAN_PHASE_MIN, ZHAN_PHASE_MAX,
               color="moccasin", alpha=0.55, zorder=0,
               label=f"April–September (phase {ZHAN_PHASE_MIN}–{ZHAN_PHASE_MAX})")

    # A1b baseline intervals
    for ps, pe in A1B_INTERVALS:
        ax.axhspan(ps, pe, color="lightgray", alpha=0.45, zorder=1)

    colors = ["orange" if c > threshold else "steelblue" for c in counts]
    ax.barh(
        bin_centers, counts,
        height=bin_width * 0.85,
        color=colors, edgecolor="white", linewidth=0.4, zorder=2,
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

    # Calendar reference lines
    for pos, lbl in CALENDAR_REFS:
        ax.axhline(pos, color="dimgray", linestyle=":", linewidth=0.6, alpha=0.7)
        ax.text(
            0.01, pos + 0.005, lbl,
            transform=ax.get_yaxis_transform(),
            fontsize=6.5, color="dimgray", va="bottom",
        )

    # Statistical annotation
    p_val = k24_data["p_chi2"]
    p_str = f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.4f}"
    mean_phase = k24_data["mean_phase"]
    in_apr_sep = trend["deep_mean_phase_in_apr_sep"]
    low_n_note = "\n(low n)" if not sufficient_n else ""

    annot = (
        f"n={n:,}\n"
        f"\u03c7\u00b2={k24_data['chi2']:.2f}\n"
        f"p={p_str}\n"
        f"V={k24_data['cramer_v']:.4f}\n"
        f"Mean phase={mean_phase:.4f}\n"
        f"In Apr–Sep: {in_apr_sep}"
        f"{low_n_note}"
    )
    ax.text(
        0.97, 0.97, annot,
        transform=ax.transAxes,
        ha="right", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
    )

    # Zhan & Shearer citation note
    ax.text(
        0.03, 0.03,
        "Zhan & Shearer (2015): 70% of deep M>7 events in April\u2013September",
        transform=ax.transAxes,
        ha="left", va="bottom", fontsize=7.5, style="italic",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="lightyellow", alpha=0.85),
    )

    ax.set_title(
        "Case B4: Solar Phase Distribution — Deep Events (>300 km, k=24)",
        fontsize=11, fontweight="bold",
    )
    ax.set_xlabel("Event Count", fontsize=9)
    ax.set_ylabel("Solar Phase (0–1)", fontsize=9)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0, 1.05, 0.25))
    ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=8)
    ax.tick_params(axis="x", labelsize=8)

    # Legend
    legend_patches = [
        mpatches.Patch(color="steelblue", label="Normal bin"),
        mpatches.Patch(color="orange", label=">1 SD elevated"),
        mpatches.Patch(color="moccasin", alpha=0.55,
                       label="April–September (Zhan & Shearer 2015)"),
        mpatches.Patch(color="lightgray", alpha=0.45, label="A1b baseline interval"),
        plt.Line2D([0], [0], color="black", linestyle="--", linewidth=1.0,
                   label="Expected count"),
        plt.Line2D([0], [0], color="gray", linestyle=":", linewidth=0.9,
                   label="1-SD threshold"),
    ]
    ax.legend(
        handles=legend_patches,
        loc="upper left", fontsize=7.5,
        bbox_to_anchor=(0.0, 0.95),
    )

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-b4-deep-events-phase.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 3 (deep events phase) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load results and generate all figures."""
    logger.info("Loading results from %s", RESULTS_PATH)
    results = load_results()

    plot_binplots(results)
    plot_depth_trend(results)
    plot_deep_events_phase(results)

    logger.info("All figures written to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
