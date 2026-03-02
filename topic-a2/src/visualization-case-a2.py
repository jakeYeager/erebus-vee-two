"""
Case A2: b-Value Seasonal Variation — Visualization Script.

Generates three figures:
  1. b-value vs solar phase at k=24 (with CI bars)
  2. b-value vs solar phase at k=32 (with CI bars)
  3. Rate and b-value dual-axis overlay at k=24

Usage:
    python topic-a2/src/visualization-case-a2.py

Outputs:
    topic-a2/output/case-a2-bvalue-phase.png
    topic-a2/output/case-a2-bvalue-k32.png
    topic-a2/output/case-a2-rate-bvalue-overlay.png
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------------------
# Project path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("visualization-case-a2")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RESULTS_PATH = BASE_DIR / "output" / "case-a2-results.json"
OUT_BVALUE_K24 = BASE_DIR / "output" / "case-a2-bvalue-phase.png"
OUT_BVALUE_K32 = BASE_DIR / "output" / "case-a2-bvalue-k32.png"
OUT_OVERLAY = BASE_DIR / "output" / "case-a2-rate-bvalue-overlay.png"

# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------
EQUINOX_PHASES = [0.19, 0.69]
SOLSTICE_PHASES = [0.44, 0.94]

# A1b baseline highlight intervals (phase ranges)
A1B_BANDS = [(0.1875, 0.25), (0.625, 0.656), (0.875, 0.917)]

DPI = 300
STEELBLUE = "steelblue"
ORANGE = "orange"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_results() -> dict:
    """Load case-a2-results.json."""
    with open(RESULTS_PATH, "r") as fh:
        return json.load(fh)


def extract_bin_arrays(bins_data: list[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extract arrays from per-bin dicts.

    Returns
    -------
    phases, b_values, b_lower, b_upper, counts, low_n_flags
    """
    phases = np.array([b["phase_center"] for b in bins_data])
    b_vals = np.array([b["b_mle"] if b["b_mle"] is not None else np.nan for b in bins_data])
    b_lower = np.array([b["b_ci95_lower"] if b["b_ci95_lower"] is not None else np.nan for b in bins_data])
    b_upper = np.array([b["b_ci95_upper"] if b["b_ci95_upper"] is not None else np.nan for b in bins_data])
    counts = np.array([b["n"] for b in bins_data])
    low_n_flags = np.array([b["low_n_flag"] for b in bins_data])
    return phases, b_vals, b_lower, b_upper, counts, low_n_flags


# ---------------------------------------------------------------------------
# Figure 1 & 2: b-value vs phase (generic)
# ---------------------------------------------------------------------------

def plot_bvalue_phase(
    results: dict,
    k: int,
    out_path: Path,
) -> None:
    """Plot b-value vs solar phase bin for a given bin count k.

    Parameters
    ----------
    results : dict
        Loaded results JSON.
    k : int
        Bin count (24 or 32).
    out_path : Path
        Output PNG path.
    """
    key = f"k{k}"
    pv = results["phase_variation"][key]
    bins_data = pv["bins"]
    b_mean = pv["b_mean"]
    f_stat = pv["f_stat"]
    p_anova = pv["p_anova"]
    r_rate_b = pv["r_rate_b"]
    p_rate_b = pv["p_rate_b"]
    b_max_bin = pv["b_max_bin"]
    b_min_bin = pv["b_min_bin"]
    b_max_phase_class = pv["b_max_phase_class"]
    b_min_phase_class = pv["b_min_phase_class"]

    phases, b_vals, b_lower, b_upper, counts, low_n_flags = extract_bin_arrays(bins_data)

    # Error bar half-widths
    err_lower = b_vals - b_lower
    err_upper = b_upper - b_vals

    fig, ax = plt.subplots(figsize=(12, 6))

    # A1b baseline shaded bands
    for (lo, hi) in A1B_BANDS:
        ax.axvspan(lo, hi, alpha=0.12, color="gray", label=None)

    # Equinox and solstice reference lines
    for ph in EQUINOX_PHASES:
        ax.axvline(ph, color="darkorange", linestyle=":", linewidth=1.0, alpha=0.8)
    for ph in SOLSTICE_PHASES:
        ax.axvline(ph, color="darkgreen", linestyle=":", linewidth=1.0, alpha=0.8)

    # Mean b-value line
    ax.axhline(b_mean, color="gray", linestyle="--", linewidth=1.2, label=f"b_mean = {b_mean:.4f}")

    # Split into normal and low-n bins
    normal_mask = ~low_n_flags & ~np.isnan(b_vals)
    lown_mask = low_n_flags & ~np.isnan(b_vals)

    # Main line (all valid bins connected)
    valid_mask = ~np.isnan(b_vals)
    ax.plot(
        phases[valid_mask], b_vals[valid_mask],
        color=STEELBLUE, linewidth=1.5, zorder=3,
    )

    # Error bars — normal bins (filled circles)
    if normal_mask.any():
        ax.errorbar(
            phases[normal_mask], b_vals[normal_mask],
            yerr=[err_lower[normal_mask], err_upper[normal_mask]],
            fmt="o", color=STEELBLUE, markersize=5,
            ecolor=STEELBLUE, elinewidth=1.0, capsize=3,
            label="b-value (95% CI)", zorder=4,
        )

    # Error bars — low-n bins (open circles)
    if lown_mask.any():
        ax.errorbar(
            phases[lown_mask], b_vals[lown_mask],
            yerr=[err_lower[lown_mask], err_upper[lown_mask]],
            fmt="o", color="white", markeredgecolor=STEELBLUE,
            markeredgewidth=1.5, markersize=5,
            ecolor=STEELBLUE, elinewidth=1.0, capsize=3,
            label="b-value, low n < 20 (95% CI)", zorder=4,
        )

    # Annotate b_max bin
    b_max_phase_c = phases[b_max_bin]
    b_max_val = b_vals[b_max_bin]
    ax.annotate(
        f"b_max\n(bin {b_max_bin}, {b_max_phase_class})",
        xy=(b_max_phase_c, b_max_val),
        xytext=(b_max_phase_c + 0.05, b_max_val + 0.02),
        fontsize=7,
        arrowprops=dict(arrowstyle="->", color="navy", lw=0.8),
        color="navy",
    )

    # Annotate b_min bin
    b_min_phase_c = phases[b_min_bin]
    b_min_val = b_vals[b_min_bin]
    ax.annotate(
        f"b_min\n(bin {b_min_bin}, {b_min_phase_class})",
        xy=(b_min_phase_c, b_min_val),
        xytext=(b_min_phase_c + 0.05, b_min_val - 0.04),
        fontsize=7,
        arrowprops=dict(arrowstyle="->", color="darkred", lw=0.8),
        color="darkred",
    )

    # Legend patches for shaded bands and reference lines
    band_patch = mpatches.Patch(color="gray", alpha=0.3, label="A1b baseline intervals")
    equinox_line = plt.Line2D([0], [0], color="darkorange", linestyle=":", label="Equinox (0.19, 0.69)")
    solstice_line = plt.Line2D([0], [0], color="darkgreen", linestyle=":", label="Solstice (0.44, 0.94)")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles + [band_patch, equinox_line, solstice_line],
        labels + [band_patch.get_label(), equinox_line.get_label(), solstice_line.get_label()],
        fontsize=7,
        loc="upper right",
    )

    # Annotations for statistics
    stats_text = (
        f"ANOVA: F={f_stat:.3f}, p={p_anova:.4f}\n"
        f"Pearson r(rate, b)={r_rate_b:.4f}, p={p_rate_b:.4f}"
    )
    ax.text(
        0.02, 0.98, stats_text,
        transform=ax.transAxes, fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    ax.set_xlabel("Solar Phase (fraction of year)", fontsize=10)
    ax.set_ylabel("Gutenberg-Richter b-value (MLE)", fontsize=10)
    ax.set_title(
        f"Case A2: b-Value vs Solar Phase (k={k} bins)\n"
        f"ISC-GEM Catalog, M ≥ 6.0, 1950–2021 (n=9,210)",
        fontsize=11,
    )
    ax.set_xlim(-0.02, 1.02)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 3: Rate and b-value overlay
# ---------------------------------------------------------------------------

def plot_rate_bvalue_overlay(results: dict, out_path: Path) -> None:
    """Plot dual-axis overlay of event rate (count) and b-value at k=24.

    Parameters
    ----------
    results : dict
        Loaded results JSON.
    out_path : Path
        Output PNG path.
    """
    key = "k24"
    pv = results["phase_variation"][key]
    bins_data = pv["bins"]
    r_rate_b = pv["r_rate_b"]
    inverse = pv["inverse_phase_relationship"]

    phases, b_vals, b_lower, b_upper, counts, low_n_flags = extract_bin_arrays(bins_data)

    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Left axis: event count (steelblue bars)
    bin_width = 1.0 / len(phases)
    ax1.bar(
        phases, counts,
        width=bin_width * 0.85,
        color=STEELBLUE, alpha=0.7,
        label="Event count per bin",
        zorder=2,
    )
    ax1.set_xlabel("Solar Phase (fraction of year)", fontsize=10)
    ax1.set_ylabel("Event Count per Bin", color=STEELBLUE, fontsize=10)
    ax1.tick_params(axis="y", labelcolor=STEELBLUE)

    # Right axis: b-value (orange line with points)
    ax2 = ax1.twinx()
    valid_mask = ~np.isnan(b_vals)
    ax2.plot(
        phases[valid_mask], b_vals[valid_mask],
        color=ORANGE, linewidth=2.0, zorder=5,
        label="b-value (MLE)",
    )
    ax2.scatter(
        phases[valid_mask], b_vals[valid_mask],
        color=ORANGE, s=40, zorder=6,
    )
    ax2.set_ylabel("Gutenberg-Richter b-value (MLE)", color=ORANGE, fontsize=10)
    ax2.tick_params(axis="y", labelcolor=ORANGE)

    # Equinox and solstice reference lines
    for ph in EQUINOX_PHASES:
        ax1.axvline(ph, color="darkorange", linestyle=":", linewidth=1.0, alpha=0.8)
    for ph in SOLSTICE_PHASES:
        ax1.axvline(ph, color="darkgreen", linestyle=":", linewidth=1.0, alpha=0.8)

    # Direction annotation
    direction_str = "negative (inverse)" if inverse else "positive (in-phase)"
    stats_text = (
        f"Pearson r(rate, b) = {r_rate_b:.4f}\n"
        f"Direction: {direction_str}"
    )
    ax1.text(
        0.02, 0.98, stats_text,
        transform=ax1.transAxes, fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    # Combined legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    equinox_line = plt.Line2D([0], [0], color="darkorange", linestyle=":", label="Equinox (0.19, 0.69)")
    solstice_line = plt.Line2D([0], [0], color="darkgreen", linestyle=":", label="Solstice (0.44, 0.94)")
    ax1.legend(
        handles1 + handles2 + [equinox_line, solstice_line],
        labels1 + labels2 + [equinox_line.get_label(), solstice_line.get_label()],
        fontsize=8,
        loc="upper right",
    )

    ax1.set_xlim(-0.02, 1.02)
    ax1.set_title(
        "Case A2: Event Rate and b-Value Overlay (k=24 bins)\n"
        "ISC-GEM Catalog, M ≥ 6.0, 1950–2021 (n=9,210)",
        fontsize=11,
    )

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved: %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Generate all Case A2 visualizations."""
    logger.info("=== Case A2 Visualization ===")

    results = load_results()

    # Figure 1: k=24
    plot_bvalue_phase(results, k=24, out_path=OUT_BVALUE_K24)

    # Figure 2: k=32
    plot_bvalue_phase(results, k=32, out_path=OUT_BVALUE_K32)

    # Figure 3: Rate-b overlay at k=24
    plot_rate_bvalue_overlay(results, out_path=OUT_OVERLAY)

    logger.info("=== Case A2 visualization complete ===")


if __name__ == "__main__":
    main()
