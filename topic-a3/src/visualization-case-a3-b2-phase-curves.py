"""
Phase curve visualization for A3.B2 — multi-line smoothed solar-phase distributions.

Produces case-a3-b2-phase-curves.png: one smooth curve per population from Table 3
(phase alignment pairs), normalized to z-score, with A1b interval shading.
NH populations in warm tones, SH in cool tones.
"""

import json
import logging
import pathlib

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.ndimage import gaussian_filter1d

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

BASE_DIR = pathlib.Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a3-b2-results.json"
OUT_PATH = BASE_DIR / "output" / "case-a3-b2-phase-curves.png"

# A1b baseline intervals (bin edges as phase fractions, k=24)
INTERVALS = {
    "Interval 1\n(Mar equinox)":  (4 / 24, 6 / 24),
    "Interval 2\n(mid-Aug)":      (15 / 24, 16 / 24),
    "Interval 3\n(late Nov)":     (21 / 24, 22 / 24),
}
INTERVAL_COLOR = "#dddddd"

# Approximate perihelion-relative solar phase → calendar month
MONTH_PHASES = [0.0, 0.085, 0.162, 0.247, 0.329, 0.413,
                0.496, 0.578, 0.661, 0.745, 0.827, 0.912]
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# Populations: (label, key_path tuple, color, linestyle, linewidth)
POPULATIONS = [
    ("Continental NH",  ("subtest_1_tectonic_hemisphere", "full", "continental", "nh"), "#c0392b", "-",  2.0),
    ("Continental SH",  ("subtest_1_tectonic_hemisphere", "full", "continental", "sh"), "#2980b9", "-",  2.0),
    ("Transitional NH", ("subtest_1_tectonic_hemisphere", "full", "transitional", "nh"), "#e67e22", "--", 1.6),
    ("Transitional SH", ("subtest_1_tectonic_hemisphere", "full", "transitional", "sh"), "#27ae60", "--", 1.6),
    ("Oceanic NH",      ("subtest_1_tectonic_hemisphere", "full", "oceanic", "nh"), "#f39c12", ":",  1.4),
    ("Oceanic SH",      ("subtest_1_tectonic_hemisphere", "full", "oceanic", "sh"), "#8e44ad", ":",  1.4),
    ("Mid-crustal NH",  ("subtest_2_midcrustal_hemisphere", "full", "nh"), "#922b21", "-",  2.6),
    ("Mid-crustal SH",  ("subtest_2_midcrustal_hemisphere", "full", "sh"), "#154360", "-",  2.6),
]

SIGMA = 1.3  # gaussian smoothing kernel (bins)


def get_cell(data: dict, key_path: tuple) -> dict:
    """Traverse nested dict by key_path tuple."""
    node = data
    for key in key_path:
        node = node[key]
    return node


def zscores(bin_counts: list) -> np.ndarray:
    """Convert bin counts to z-scores relative to uniform expected."""
    counts = np.array(bin_counts, dtype=float)
    expected = counts.sum() / len(counts)
    return (counts - expected) / np.sqrt(expected)


def smooth_wrap(z: np.ndarray, sigma: float) -> tuple:
    """
    Gaussian-smooth z-scores with circular wrapping. Returns (phases, smoothed).

    Parameters
    ----------
    z : np.ndarray
        Z-score array of length k.
    sigma : float
        Smoothing width in bins.
    """
    k = len(z)
    tiled = np.tile(z, 3)
    smoothed = gaussian_filter1d(tiled, sigma=sigma)[k:2 * k]
    phases = (np.arange(k) + 0.5) / k
    # Wrap last → first for closed curve
    return np.append(phases, 1.0), np.append(smoothed, smoothed[0])


def main() -> None:
    log.info("Loading results JSON")
    with open(RESULTS_PATH) as f:
        data = json.load(f)

    # --- Pre-compute all smoothed curves and y range ---
    curves = []
    for label, key_path, color, ls, lw in POPULATIONS:
        cell = get_cell(data, key_path)
        if not cell.get("bin_counts") or cell.get("low_n"):
            log.warning("Skipping %s", label)
            continue
        z = zscores(cell["bin_counts"])
        phases, sz = smooth_wrap(z, SIGMA)
        peak_bin = int(np.argmax(cell["bin_counts"]))
        peak_phase = (peak_bin + 0.5) / 24
        curves.append((label, color, ls, lw, phases, sz, z, peak_phase, peak_bin))

    all_z = np.concatenate([c[6] for c in curves])
    y_min = all_z.min() - 0.4
    y_max = all_z.max() + 0.8  # headroom for interval labels

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.set_xlim(0, 1)
    ax.set_ylim(y_min, y_max)

    # --- Interval shading (drawn first, behind everything) ---
    label_y = y_max - 0.15
    for iname, (start, end) in INTERVALS.items():
        ax.axvspan(start, end, color=INTERVAL_COLOR, alpha=0.7, zorder=0)
        ax.text(
            (start + end) / 2, label_y,
            iname,
            ha="center", va="top", fontsize=8, color="#777777",
            linespacing=1.3,
        )

    # --- Reference lines ---
    ax.axhline(0, color="#333333", linewidth=0.8, linestyle="-", zorder=1, alpha=0.5)
    ax.axhline(1.0, color="#888888", linewidth=0.7, linestyle="--", zorder=1, alpha=0.5)
    ax.axhline(-1.0, color="#888888", linewidth=0.7, linestyle="--", zorder=1, alpha=0.5)
    ax.text(0.995, 1.05, "1-SD", ha="right", va="bottom", fontsize=7.5,
            color="#888888", transform=ax.get_xaxis_transform())

    # --- Draw curves ---
    for label, color, ls, lw, phases, sz, z, peak_phase, peak_bin in curves:
        ax.plot(phases, sz, color=color, linestyle=ls, linewidth=lw,
                label=label, zorder=3, alpha=0.88)
        # Peak dot on smoothed curve at peak bin phase
        peak_smoothed = sz[np.argmin(np.abs(phases - peak_phase))]
        ax.plot(peak_phase, peak_smoothed, "o", color=color, markersize=6,
                markeredgewidth=0.9, markeredgecolor="white", zorder=5)

    # --- Axes formatting ---
    ax.set_xlabel("Solar phase (0 = perihelion ≈ Jan 3)", fontsize=11)
    ax.set_ylabel("Z-score (deviation from uniform)", fontsize=11)
    ax.set_title(
        "A3.B2 — Solar Phase Distribution by Population and Hemisphere (Full Catalog)\n"
        "Warm tones = Northern Hemisphere · Cool tones = Southern Hemisphere · "
        "Dots = peak bin · Shading = A1b elevated intervals",
        fontsize=10.5,
    )
    ax.set_xticks(MONTH_PHASES)
    ax.set_xticklabels(MONTH_LABELS, fontsize=9.5)
    ax.set_xticks([i / 24 for i in range(25)], minor=True)
    ax.tick_params(axis="x", which="minor", length=2.5, color="#aaaaaa")
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6, zorder=0)

    # --- Legend ---
    nh_handles = [(h, l) for h, l in zip(*ax.get_legend_handles_labels()) if "NH" in l]
    sh_handles = [(h, l) for h, l in zip(*ax.get_legend_handles_labels()) if "SH" in l]
    interval_patch = mpatches.Patch(color=INTERVAL_COLOR, alpha=0.7,
                                    label="A1b intervals")
    all_handles = [h for h, _ in nh_handles] + [h for h, _ in sh_handles] + [interval_patch]
    all_labels  = [l for _, l in nh_handles] + [l for _, l in sh_handles] + ["A1b intervals"]
    ax.legend(all_handles, all_labels, loc="lower right", fontsize=8.5,
              ncol=3, framealpha=0.9, edgecolor="#cccccc")

    # --- Stats annotation ---
    stats = (
        "NH peak phases  —  Continental: ≈0.27  Transitional: ≈0.23  Mid-crustal: ≈0.23\n"
        "SH peak phases  —  Continental: ≈0.65  Transitional: ≈0.65  Mid-crustal: ≈0.65\n"
        "NH→SH offset: ~0.4 phase units (~5 months)"
    )
    ax.text(0.01, 0.03, stats, transform=ax.transAxes,
            fontsize=8, verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#cccccc", alpha=0.88))

    fig.tight_layout(rect=[0, 0, 1, 1])
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight")
    log.info("Saved: %s (%.0f KB)", OUT_PATH, OUT_PATH.stat().st_size / 1024)
    plt.close()


if __name__ == "__main__":
    main()
