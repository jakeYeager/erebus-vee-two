"""Case A4: Visualization script.

Generates three PNG figures:
  1. case-a4-sub-a-binplot.png  — Sub-A phase bin distributions (4 panels, k=24)
  2. case-a4-sub-b-intervals.png — Sub-B interval maps (3x3 grid)
  3. case-a4-sub-c-aftershock.png — Sub-C aftershock distributions (3 panels, k=24)
"""

import json
import logging
import sys
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
RESULTS_PATH = OUTPUT_DIR / "case-a4-results.json"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("visualization-case-a4")

# Solar calendar reference positions (fraction of year)
EQUINOX_SPRING = 0.19   # ~March 10–20
SOLSTICE_SUMMER = 0.44  # ~June 20
EQUINOX_AUTUMN = 0.69   # ~September 22
SOLSTICE_WINTER = 0.94  # ~December 21

# A1b baseline intervals
A1B_INTERVALS = [
    (0.1875, 0.25),
    (0.625, 0.656),
    (0.875, 0.917),
]

DPI = 300


def load_results() -> dict:
    """Load the case A4 results JSON.

    Returns:
        Parsed results dictionary.
    """
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


def compute_threshold(n: int, k: int) -> float:
    """Compute 1-SD elevated-bin threshold.

    Args:
        n: Total events.
        k: Number of bins.

    Returns:
        Threshold = E + sqrt(E).
    """
    expected = n / k
    return expected + np.sqrt(expected)


# ---------------------------------------------------------------------------
# Figure 1 — Sub-analysis A bin distributions (k=24)
# ---------------------------------------------------------------------------
def plot_sub_a(results: dict) -> None:
    """Generate Figure 1: 4-panel phase bin distributions.

    Args:
        results: Full results dict.
    """
    sub_a = results["sub_a"]

    catalog_info = [
        ("raw", "Raw Catalog"),
        ("gk_mainshocks", "G-K Mainshocks"),
        ("reas_mainshocks", "Reasenberg Mainshocks"),
        ("a1b_mainshocks", "A1b Mainshocks"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(20, 6), sharey=False)
    fig.suptitle(
        "Case A4 Sub-A: Solar Phase Bin Distributions at k=24",
        fontsize=14, fontweight="bold", y=1.02,
    )

    for ax, (key, label) in zip(axes, catalog_info):
        data = sub_a[key]["k24"]
        n = data["n"]
        k = data["k"]
        bin_counts = np.array(data["bin_counts"])
        expected = n / k
        threshold = compute_threshold(n, k)

        # Build phase bin centers
        bin_edges = np.linspace(0, 1, k + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_width = 1.0 / k

        # Color bars
        colors = ["orange" if bc > threshold else "steelblue" for bc in bin_counts]

        ax.barh(
            bin_centers, bin_counts,
            height=bin_width * 0.85,
            color=colors, edgecolor="white", linewidth=0.4,
        )

        # Expected dashed line
        ax.axvline(expected, color="black", linestyle="--", linewidth=1.0, label=f"E={expected:.1f}")

        # Calendar reference lines
        for pos, lbl in [
            (EQUINOX_SPRING, "Spr Eq"),
            (SOLSTICE_SUMMER, "Sum Sol"),
            (EQUINOX_AUTUMN, "Aut Eq"),
            (SOLSTICE_WINTER, "Win Sol"),
        ]:
            ax.axhline(pos, color="gray", linestyle=":", linewidth=0.8)

        # Annotations
        chi2 = data["chi2"]
        p_chi2 = data["p_chi2"]
        cramer_v = data["cramer_v"]
        p_str = f"{p_chi2:.2e}" if p_chi2 < 0.001 else f"{p_chi2:.4f}"
        annot = f"n={n:,}\n\u03c7\u00b2={chi2:.2f}\np={p_str}\nV={cramer_v:.3f}"
        ax.text(
            0.97, 0.97, annot,
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.set_xlabel("Event Count", fontsize=9)
        ax.set_ylabel("Solar Phase (fraction of year)", fontsize=9)
        ax.set_ylim(0, 1)
        ax.invert_yaxis()

        # Y-axis tick labels as fractions
        ax.set_yticks(np.arange(0, 1.1, 0.25))
        ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=7)

        # Legend
        legend_patches = [
            mpatches.Patch(color="steelblue", label="Normal"),
            mpatches.Patch(color="orange", label=">1 SD elevated"),
            plt.Line2D([0], [0], color="black", linestyle="--", label=f"Expected"),
            plt.Line2D([0], [0], color="gray", linestyle=":", label="Equinox/Solstice"),
        ]
        ax.legend(handles=legend_patches, fontsize=6, loc="lower right")

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a4-sub-a-binplot.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 1 saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 2 — Sub-analysis B interval maps (3x3 grid)
# ---------------------------------------------------------------------------
def plot_sub_b(results: dict) -> None:
    """Generate Figure 2: 3x3 grid of interval maps.

    Args:
        results: Full results dict.
    """
    sub_b = results["sub_b"]

    catalog_keys = ["gk_mainshocks", "reas_mainshocks", "a1b_mainshocks"]
    catalog_labels = ["G-K Mainshocks", "Reasenberg Mainshocks", "A1b Mainshocks"]
    bin_counts = [16, 24, 32]

    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    fig.suptitle(
        "Case A4 Sub-B: Elevated Phase Interval Maps vs A1b Baseline",
        fontsize=14, fontweight="bold",
    )

    for row_idx, (cat_key, cat_label) in enumerate(zip(catalog_keys, catalog_labels)):
        for col_idx, k in enumerate(bin_counts):
            ax = axes[row_idx][col_idx]
            k_key = f"k{k}"
            k_data = sub_b[cat_key][k_key]

            # Draw A1b baseline intervals as gray bands
            for (ps, pe) in A1B_INTERVALS:
                ax.axhspan(ps, pe, color="lightgray", alpha=0.6, label="A1b baseline")

            # Draw recovered intervals
            recovered = k_data.get("recovered_intervals", [])
            for ri in recovered:
                ps = ri["phase_start"]
                pe = ri["phase_end"]
                classification = ri["classification"]
                color = "green" if classification.startswith("matches") else "red"
                ax.axhspan(ps, pe, color=color, alpha=0.4)

            # Calendar reference lines
            for pos in [EQUINOX_SPRING, SOLSTICE_SUMMER, EQUINOX_AUTUMN, SOLSTICE_WINTER]:
                ax.axhline(pos, color="navy", linestyle=":", linewidth=0.7, alpha=0.6)

            ax.set_title(f"{cat_label} — k={k}", fontsize=9, fontweight="bold")
            ax.set_ylim(0, 1)
            ax.set_xlim(0, 1)
            ax.set_xlabel("Arbitrary x", fontsize=7)
            ax.set_ylabel("Solar Phase (0–1)", fontsize=7)
            ax.tick_params(labelbottom=False, bottom=False)
            ax.set_yticks(np.arange(0, 1.1, 0.25))
            ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=7)

            # Annotate recovered count
            n_rec = len(recovered)
            n_match = sum(1 for ri in recovered if ri["classification"].startswith("matches"))
            ax.text(
                0.03, 0.97,
                f"{n_rec} intervals\n{n_match} match baseline",
                transform=ax.transAxes,
                va="top", ha="left", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85),
            )

    # Add legend to last axes
    legend_patches = [
        mpatches.Patch(color="lightgray", alpha=0.6, label="A1b baseline interval"),
        mpatches.Patch(color="green", alpha=0.4, label="Recovered — matches baseline"),
        mpatches.Patch(color="red", alpha=0.4, label="Recovered — new interval"),
        plt.Line2D([0], [0], color="navy", linestyle=":", label="Equinox/Solstice"),
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center", ncol=4, fontsize=8,
        bbox_to_anchor=(0.5, -0.02),
    )

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    out_path = OUTPUT_DIR / "case-a4-sub-b-intervals.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 2 saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 3 — Sub-analysis C aftershock distributions (k=24)
# ---------------------------------------------------------------------------
def plot_sub_c(results: dict) -> None:
    """Generate Figure 3: 3-panel aftershock phase distributions.

    Args:
        results: Full results dict.
    """
    sub_c = results["sub_c"]

    catalog_info = [
        ("gk_aftershocks", "G-K Aftershocks"),
        ("reas_aftershocks", "Reasenberg Aftershocks"),
        ("a1b_aftershocks", "A1b Aftershocks"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=False)
    fig.suptitle(
        "Case A4 Sub-C: Aftershock Solar Phase Distributions at k=24",
        fontsize=14, fontweight="bold", y=1.02,
    )

    for ax, (key, label) in zip(axes, catalog_info):
        data = sub_c[key]["k24"]
        classification = sub_c[key]["classification"]
        n = data["n"]
        k = data["k"]
        bin_counts_arr = np.array(data["bin_counts"])
        expected = n / k
        threshold = compute_threshold(n, k)

        bin_edges = np.linspace(0, 1, k + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_width = 1.0 / k

        colors = ["orange" if bc > threshold else "steelblue" for bc in bin_counts_arr]

        ax.barh(
            bin_centers, bin_counts_arr,
            height=bin_width * 0.85,
            color=colors, edgecolor="white", linewidth=0.4,
        )
        ax.axvline(expected, color="black", linestyle="--", linewidth=1.0)

        # Calendar reference lines
        for pos in [EQUINOX_SPRING, SOLSTICE_SUMMER, EQUINOX_AUTUMN, SOLSTICE_WINTER]:
            ax.axhline(pos, color="gray", linestyle=":", linewidth=0.8)

        chi2 = data["chi2"]
        p_chi2 = data["p_chi2"]
        p_str = f"{p_chi2:.2e}" if p_chi2 < 0.001 else f"{p_chi2:.4f}"
        annot = (
            f"n={n:,}\n\u03c7\u00b2={chi2:.2f}\np={p_str}\n\n"
            f"Classification:\n{classification}"
        )
        ax.text(
            0.97, 0.97, annot,
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.set_xlabel("Event Count", fontsize=9)
        ax.set_ylabel("Solar Phase (fraction of year)", fontsize=9)
        ax.set_ylim(0, 1)
        ax.invert_yaxis()
        ax.set_yticks(np.arange(0, 1.1, 0.25))
        ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=7)

        legend_patches = [
            mpatches.Patch(color="steelblue", label="Normal"),
            mpatches.Patch(color="orange", label=">1 SD elevated"),
            plt.Line2D([0], [0], color="black", linestyle="--", label="Expected"),
            plt.Line2D([0], [0], color="gray", linestyle=":", label="Equinox/Solstice"),
        ]
        ax.legend(handles=legend_patches, fontsize=6, loc="lower right")

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a4-sub-c-aftershock.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Figure 3 saved: %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load results and generate all three figures."""
    logger.info("Loading results from %s", RESULTS_PATH)
    results = load_results()

    plot_sub_a(results)
    plot_sub_b(results)
    plot_sub_c(results)

    logger.info("All figures written to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
