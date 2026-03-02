"""
Case B6: Visualization — Rolling Window Stationarity Test

Generates:
  - case-b6-trajectory.png: Time trajectory of Rayleigh R, p-value, mean phase
  - case-b6-phase-stability.png: Circular polar plot of mean phase by window, colored by decade
"""

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-b6-results.json"
OUTPUT_DIR = BASE_DIR / "output"

# A1b baseline phase centers (fraction)
A1B_PHASES = [0.22, 0.64, 0.90]

# Decade color map for circular plot
DECADE_COLORS = {
    1950: "blue",
    1960: "steelblue",
    1970: "orange",
    1980: "green",
    1990: "purple",
    2000: "red",
    2010: "black",
}


def get_decade(start_year: int) -> int:
    """Return the decade bucket for a window start year."""
    return (start_year // 10) * 10


def load_results() -> dict:
    """Load case-b6-results.json."""
    with open(RESULTS_PATH, "r", encoding="utf-8") as fh:
        return json.load(fh)


def plot_trajectory(windows: list[dict]) -> None:
    """Figure 1: 3-row stacked time trajectory plot.

    Args:
        windows: List of per-window result dicts.
    """
    # Window center year = start + WINDOW_YEARS/2 - 0.5
    centers = [w["window_start"] + 4.5 for w in windows]
    r_values = [w["rayleigh_R"] for w in windows]
    p_values = [w["p_rayleigh"] for w in windows]
    phases = [w["mean_phase"] for w in windows]
    is_1970s = [w["is_1970s_window"] for w in windows]

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    # --- 1970s shading helper ---
    def shade_1970s(ax: plt.Axes) -> None:
        """Draw light-yellow vertical band over 1970s window centers."""
        shaded = False
        for i, (c, flag) in enumerate(zip(centers, is_1970s)):
            if flag and not shaded:
                x_start = c
                shaded = True
            if shaded and (not flag or i == len(centers) - 1):
                x_end = c if not flag else centers[-1]
                ax.axvspan(x_start - 0.5, x_end + 0.5, color="lightyellow",
                           alpha=0.9, zorder=0, label="1970s windows" if ax == axes[0] else None)
                shaded = False

    # For the 1970s shading: shade between center 1974.5 (start 1970) and 1983.5 (start 1979)
    # Find range of 1970s-flagged centers
    flag_centers = [c for c, f in zip(centers, is_1970s) if f]
    if flag_centers:
        x_lo = min(flag_centers) - 0.5
        x_hi = max(flag_centers) + 0.5

    # Row 1: Rayleigh R
    ax1 = axes[0]
    ax1.plot(centers, r_values, color="steelblue", linewidth=1.5, zorder=2)
    if flag_centers:
        ax1.axvspan(x_lo, x_hi, color="lightyellow", alpha=0.9, zorder=0, label="1970s windows")
    ax1.set_ylim(0.0, 0.1)
    ax1.set_ylabel("Rayleigh R", fontsize=11)
    ax1.set_title("Case B6: Rolling Window Stationarity — Solar Phase Signal (1950–2021)", fontsize=13, pad=10)
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(axis="y", alpha=0.3, linestyle="--")

    # Row 2: p-value (log scale)
    ax2 = axes[1]
    ax2.plot(centers, p_values, color="red", linewidth=1.5, zorder=2)
    if flag_centers:
        ax2.axvspan(x_lo, x_hi, color="lightyellow", alpha=0.9, zorder=0)
    ax2.axhline(0.05, color="gray", linestyle="--", linewidth=1.2, label="p = 0.05")
    ax2.axhline(0.001, color="darkgray", linestyle=":", linewidth=1.2, label="p = 0.001")
    ax2.set_yscale("log")
    ax2.set_ylim(0.001, 1.0)
    ax2.set_ylabel("Rayleigh p-value (log)", fontsize=11)
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(axis="y", alpha=0.3, linestyle="--")

    # Row 3: mean phase fraction
    ax3 = axes[2]
    ax3.plot(centers, phases, color="black", linewidth=1.2, zorder=2)
    ax3.scatter(centers, phases, color="black", s=18, zorder=3)
    if flag_centers:
        ax3.axvspan(x_lo, x_hi, color="lightyellow", alpha=0.9, zorder=0)
    for ph, lbl in zip(A1B_PHASES, ["A1b 0.22", "A1b 0.64", "A1b 0.90"]):
        ax3.axhline(ph, color="dimgray", linestyle="--", linewidth=1.0, label=lbl, alpha=0.7)
    ax3.set_ylim(-0.05, 1.05)
    ax3.set_ylabel("Mean Phase Fraction", fontsize=11)
    ax3.set_xlabel("Window center year", fontsize=11)
    ax3.legend(loc="upper right", fontsize=8, ncol=3)
    ax3.grid(axis="y", alpha=0.3, linestyle="--")

    # x-axis ticks every 5 years
    tick_centers = [y + 4.5 for y in range(1950, 2020, 5)]
    tick_labels = [str(y) for y in range(1950, 2020, 5)]
    ax3.set_xticks(tick_centers)
    ax3.set_xticklabels(tick_labels, fontsize=9)

    out_path = OUTPUT_DIR / "case-b6-trajectory.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_phase_stability(windows: list[dict]) -> None:
    """Figure 2: Circular polar plot of mean phase angle per window, colored by decade.

    Args:
        windows: List of per-window result dicts.
    """
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="polar")

    # Scale Rayleigh R to visible radius range [0.4, 1.0]
    r_vals = np.array([w["rayleigh_R"] for w in windows])
    r_min, r_max = r_vals.min(), r_vals.max()
    r_range = r_max - r_min if r_max > r_min else 1.0
    radii = 0.4 + 0.6 * (r_vals - r_min) / r_range

    # Plot A1b baseline arcs as gray wedges (±0.04 width around center)
    arc_width = 0.08  # phase fraction
    for ph in A1B_PHASES:
        theta_center = 2.0 * math.pi * ph
        theta_lo = 2.0 * math.pi * (ph - arc_width / 2.0)
        theta_hi = 2.0 * math.pi * (ph + arc_width / 2.0)
        theta_arc = np.linspace(theta_lo, theta_hi, 50)
        ax.fill_between(theta_arc, 0.85, 1.0, color="lightgray", alpha=0.5, zorder=0)
        ax.text(theta_center, 1.08, f"A1b\n{ph:.2f}",
                ha="center", va="center", fontsize=7.5, color="gray")

    # Plot windows colored by decade
    plotted_decades = set()
    for w, radius in zip(windows, radii):
        decade = get_decade(w["window_start"])
        color = DECADE_COLORS.get(decade, "gray")
        theta = 2.0 * math.pi * w["mean_phase"]
        ax.scatter(theta, radius, color=color, s=55, zorder=3,
                   alpha=0.85, label=f"{decade}s" if decade not in plotted_decades else None)
        plotted_decades.add(decade)

    # Legend for decades
    legend_handles = []
    for decade in sorted(DECADE_COLORS.keys()):
        if decade in plotted_decades:
            patch = mpatches.Patch(color=DECADE_COLORS[decade], label=f"{decade}s")
            legend_handles.append(patch)
    ax.legend(handles=legend_handles, loc="lower right",
              bbox_to_anchor=(1.25, -0.05), fontsize=9, title="Decade")

    # Polar plot aesthetics
    ax.set_theta_zero_location("N")  # phase 0 at top
    ax.set_theta_direction(-1)       # clockwise (calendar convention)
    ax.set_rlim(0.0, 1.2)
    ax.set_rticks([0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["min R", "", "", "max R"], fontsize=7)
    ax.set_xticks([2.0 * math.pi * p for p in [0.0, 0.25, 0.5, 0.75]])
    ax.set_xticklabels(["Jan 1\n(phase 0)", "Apr\n(0.25)", "Jul\n(0.50)", "Oct\n(0.75)"], fontsize=8)
    ax.set_title("Mean Phase Angle per 10-Year Window\n(radius ∝ Rayleigh R; shaded arcs = A1b baseline intervals)",
                 fontsize=11, pad=22)

    ax.grid(alpha=0.3)

    out_path = OUTPUT_DIR / "case-b6-phase-stability.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def main() -> None:
    """Main entry point for Case B6 visualization."""
    results = load_results()
    windows = results["windows"]

    plot_trajectory(windows)
    plot_phase_stability(windows)
    print("All Case B6 figures saved.")


if __name__ == "__main__":
    main()
