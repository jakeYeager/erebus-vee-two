"""Case B1: Hemisphere Stratification — Visualization Script.

Generates three PNG figures:
  1. case-b1-binplot-nh.png     — NH solar-phase bin distribution at k=24
  2. case-b1-binplot-sh.png     — SH solar-phase bin distribution at k=24
  3. case-b1-interval-comparison.png — Side-by-side NH vs SH interval maps
                                       across k=16, 24, 32
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
RESULTS_PATH = OUTPUT_DIR / "case-b1-results.json"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("visualization-case-b1")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DPI = 300

# A1b baseline elevated phase intervals
A1B_INTERVALS = [
    (0.1875, 0.25),   # Interval 1: March equinox
    (0.625, 0.656),   # Interval 2: ~mid-August
    (0.875, 0.917),   # Interval 3: ~mid-November
]

# Solar calendar reference positions (fraction of year)
EQUINOX_SPRING = 0.19
SOLSTICE_SUMMER = 0.44
EQUINOX_AUTUMN = 0.69
SOLSTICE_WINTER = 0.94


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_results() -> dict:
    """Load the case B1 results JSON.

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
# Figure helper: single hemisphere bin distribution
# ---------------------------------------------------------------------------
def plot_hemisphere_binplot(
    ax: plt.Axes,
    bin_counts: list[int],
    n: int,
    k: int,
    title: str,
) -> None:
    """Draw a horizontal bar chart of phase bin distribution on the given axes.

    Args:
        ax: Matplotlib axes to draw on.
        bin_counts: Observed bin counts.
        n: Total events.
        k: Number of bins.
        title: Axes title.
    """
    counts = np.array(bin_counts, dtype=float)
    expected = n / k
    threshold = compute_threshold(n, k)

    bin_edges = np.linspace(0, 1, k + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = 1.0 / k

    # Color elevated bars orange
    colors = ["orange" if c > threshold else "steelblue" for c in counts]

    ax.barh(
        bin_centers, counts,
        height=bin_width * 0.85,
        color=colors, edgecolor="white", linewidth=0.4,
    )

    # Expected dashed line
    ax.axvline(
        expected, color="black", linestyle="--", linewidth=1.2,
        label=f"Expected ({expected:.1f})",
    )
    # 1-SD threshold dashed line
    ax.axvline(
        threshold, color="gray", linestyle=":", linewidth=1.0,
        label=f"1-SD threshold ({threshold:.1f})",
    )

    # A1b baseline interval shaded bands
    for i, (ps, pe) in enumerate(A1B_INTERVALS):
        ax.axhspan(ps, pe, color="lightgray", alpha=0.55, zorder=0)

    # Calendar reference lines
    for pos, lbl in [
        (EQUINOX_SPRING, "Spr Eq"),
        (SOLSTICE_SUMMER, "Sum Sol"),
        (EQUINOX_AUTUMN, "Aut Eq"),
        (SOLSTICE_WINTER, "Win Sol"),
    ]:
        ax.axhline(pos, color="dimgray", linestyle=":", linewidth=0.7, alpha=0.7)

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("Event Count", fontsize=9)
    ax.set_ylabel("Solar Phase (fraction of year)", fontsize=9)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0, 1.05, 0.25))
    ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=8)


def add_binplot_legend_and_annotations(
    ax: plt.Axes,
    n: int,
    chi2: float,
    p_chi2: float,
    cramer_v: float,
) -> None:
    """Add statistical annotation text box and legend to a binplot axes.

    Args:
        ax: Matplotlib axes.
        n: Event count.
        chi2: Chi-square statistic.
        p_chi2: Chi-square p-value.
        cramer_v: Cramér's V effect size.
    """
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

    legend_patches = [
        mpatches.Patch(color="steelblue", label="Normal bin"),
        mpatches.Patch(color="orange", label=">1 SD elevated"),
        mpatches.Patch(color="lightgray", alpha=0.55, label="A1b baseline interval"),
        plt.Line2D([0], [0], color="black", linestyle="--", linewidth=1.2,
                   label="Expected count"),
        plt.Line2D([0], [0], color="gray", linestyle=":", linewidth=1.0,
                   label="1-SD threshold"),
    ]
    ax.legend(handles=legend_patches, fontsize=7, loc="lower right")


# ---------------------------------------------------------------------------
# Figure 1 — NH bin distribution at k=24
# ---------------------------------------------------------------------------
def plot_nh_binplot(results: dict) -> None:
    """Generate Figure 1: NH solar-phase bin distribution at k=24.

    Args:
        results: Full results dict.
    """
    nh_stats = results["hemisphere_stats"]["nh"]
    n_nh = results["hemisphere_stats"]["n_nh"]
    k24 = nh_stats["k24"]

    fig, ax = plt.subplots(figsize=(8, 9))
    fig.suptitle(
        "Case B1: Northern Hemisphere Solar Phase Distribution",
        fontsize=13, fontweight="bold",
    )

    plot_hemisphere_binplot(
        ax,
        k24["bin_counts"],
        n=k24["n"],
        k=24,
        title=f"Northern Hemisphere (n={n_nh:,})",
    )
    add_binplot_legend_and_annotations(
        ax,
        n=k24["n"],
        chi2=k24["chi2"],
        p_chi2=k24["p_chi2"],
        cramer_v=k24["cramer_v"],
    )

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-b1-binplot-nh.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 1 (NH binplot) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 2 — SH bin distribution at k=24
# ---------------------------------------------------------------------------
def plot_sh_binplot(results: dict) -> None:
    """Generate Figure 2: SH solar-phase bin distribution at k=24.

    Args:
        results: Full results dict.
    """
    sh_stats = results["hemisphere_stats"]["sh"]
    n_sh = results["hemisphere_stats"]["n_sh"]
    k24 = sh_stats["k24"]

    fig, ax = plt.subplots(figsize=(8, 9))
    fig.suptitle(
        "Case B1: Southern Hemisphere Solar Phase Distribution",
        fontsize=13, fontweight="bold",
    )

    plot_hemisphere_binplot(
        ax,
        k24["bin_counts"],
        n=k24["n"],
        k=24,
        title=f"Southern Hemisphere (n={n_sh:,})",
    )
    add_binplot_legend_and_annotations(
        ax,
        n=k24["n"],
        chi2=k24["chi2"],
        p_chi2=k24["p_chi2"],
        cramer_v=k24["cramer_v"],
    )

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-b1-binplot-sh.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 2 (SH binplot) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 3 — Side-by-side interval comparison (2-col × 3-row)
# ---------------------------------------------------------------------------
def _interval_recovered_labels(elevated: list[dict]) -> str:
    """Summarize which A1b intervals were recovered.

    Args:
        elevated: List of elevated interval dicts.

    Returns:
        Multi-line string labeling each A1b interval as recovered or not.
    """
    def overlaps(ps: float, pe: float) -> bool:
        for ei in elevated:
            width = pe - ps
            if width <= 0:
                continue
            overlap = max(0.0, min(pe, ei["phase_end"]) - max(ps, ei["phase_start"]))
            if overlap / width > 0.5:
                return True
        return False

    parts = []
    for i, (ps, pe) in enumerate(A1B_INTERVALS, start=1):
        symbol = "Y" if overlaps(ps, pe) else "N"
        parts.append(f"Int{i}: {symbol}")
    return "\n".join(parts)


def plot_interval_comparison(results: dict) -> None:
    """Generate Figure 3: 2-col x 3-row side-by-side interval maps.

    Columns: NH (left), SH (right).
    Rows: k=16, k=24, k=32.

    Args:
        results: Full results dict.
    """
    bin_counts = [16, 24, 32]
    hemispheres = [
        ("nh", "Northern Hemisphere"),
        ("sh", "Southern Hemisphere"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(12, 14))
    fig.suptitle(
        "Case B1: Elevated Phase Intervals vs A1b Baseline — NH vs SH",
        fontsize=13, fontweight="bold",
    )

    for row_idx, k in enumerate(bin_counts):
        k_key = f"k{k}"
        bin_width = 1.0 / k

        for col_idx, (hemi_key, hemi_label) in enumerate(hemispheres):
            ax = axes[row_idx][col_idx]
            hemi_data = results["hemisphere_stats"][hemi_key][k_key]
            n = hemi_data["n"]
            bin_counts_arr = np.array(hemi_data["bin_counts"], dtype=float)
            elevated_intervals = hemi_data["elevated_intervals"]
            expected = n / k

            # Phase axis (y-axis): 0 to 1
            bin_edges = np.linspace(0, 1, k + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            # Draw A1b baseline intervals as gray shaded bands
            for ps, pe in A1B_INTERVALS:
                ax.axhspan(ps, pe, color="lightgray", alpha=0.55, zorder=0,
                           label="A1b baseline")

            # Draw elevated bins as steelblue vertical bars (horizontal bars in
            # phase orientation: use axhspan for bin spans)
            threshold = expected + np.sqrt(expected)
            for bin_i, count in enumerate(bin_counts_arr):
                if count > threshold:
                    bp_start = bin_i * bin_width
                    bp_end = (bin_i + 1) * bin_width
                    ax.axhspan(bp_start, bp_end, color="steelblue", alpha=0.7,
                               zorder=1)

            # Calendar reference lines
            for pos in [EQUINOX_SPRING, SOLSTICE_SUMMER, EQUINOX_AUTUMN,
                        SOLSTICE_WINTER]:
                ax.axhline(pos, color="navy", linestyle=":", linewidth=0.7,
                           alpha=0.5)

            ax.set_ylim(0, 1)
            ax.set_xlim(0, 1)
            ax.set_title(
                f"{hemi_label} — k={k}", fontsize=9, fontweight="bold"
            )
            ax.set_ylabel("Solar Phase (0–1)", fontsize=8)
            ax.tick_params(labelbottom=False, bottom=False)
            ax.set_yticks(np.arange(0, 1.05, 0.25))
            ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=7)

            # Annotate with A1b interval recovery status
            recovery_label = _interval_recovered_labels(elevated_intervals)
            ax.text(
                0.03, 0.97,
                recovery_label,
                transform=ax.transAxes,
                va="top", ha="left", fontsize=7.5,
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          alpha=0.9),
            )

    # Shared legend
    legend_patches = [
        mpatches.Patch(color="lightgray", alpha=0.55,
                       label="A1b baseline interval"),
        mpatches.Patch(color="steelblue", alpha=0.7,
                       label="Elevated bin (>1 SD)"),
        plt.Line2D([0], [0], color="navy", linestyle=":", linewidth=0.7,
                   label="Equinox / Solstice"),
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center", ncol=3, fontsize=8,
        bbox_to_anchor=(0.5, -0.01),
    )

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    out_path = OUTPUT_DIR / "case-b1-interval-comparison.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 3 (interval comparison) saved: %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load results and generate all three figures."""
    logger.info("Loading results from %s", RESULTS_PATH)
    results = load_results()

    plot_nh_binplot(results)
    plot_sh_binplot(results)
    plot_interval_comparison(results)

    logger.info("All figures written to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
