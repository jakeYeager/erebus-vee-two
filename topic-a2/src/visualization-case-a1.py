"""
Case A1 — Visualization Script.

Generates three figures from the Case A1 results JSON:
  1. Schuster power spectrum (3-panel, one per catalog)
  2. MFPA periodogram scan (3-panel, one per catalog)
  3. Harmonic-interval overlay (single panel)

Outputs:
    topic-a2/output/case-a1-schuster-spectrum.png
    topic-a2/output/case-a1-mfpa-scan.png
    topic-a2/output/case-a1-harmonic-intervals.png
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("visualization-case-a1")

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a1-results.json"
OUTPUT_DIR = BASE_DIR / "output"

# ---------------------------------------------------------------------------
# Named period annotations (days)
# ---------------------------------------------------------------------------
ANNOTATION_PERIODS = {
    "12h": 0.5,
    "24h": 1.0,
    "14d": 14.77,
    "27d": 27.32,
    "91d": 91.3,
    "182d": 182.625,
    "365d": 365.25,
    "548d": 548.0,
}

ANNUAL_PERIODS = {
    "annual_365": (365.25, "red", "Annual (365d)"),
    "half_year_182": (182.625, "darkorange", "Half-year (182d)"),
    "third_year_122": (121.75, "green", "Third-year (122d)"),
}

TIDAL_PERIODS = [0.5, 1.0, 14.77, 27.32]

CATALOG_KEYS = ["raw", "gk_mainshocks", "a1b_mainshocks"]
CATALOG_TITLES = {
    "raw": "Raw ISC-GEM (n=9,210)",
    "gk_mainshocks": "G-K Mainshocks",
    "a1b_mainshocks": "A1b Mainshocks",
}

# A1b intervals (phase fraction)
A1B_INTERVALS: List[Tuple[float, float, str]] = [
    (0.1875, 0.25, "Interval 1 (~Mar)"),
    (0.625, 0.656, "Interval 2 (~Aug)"),
    (0.875, 0.917, "Interval 3 (~Nov)"),
]


# ---------------------------------------------------------------------------
# Figure 1 — Schuster power spectrum
# ---------------------------------------------------------------------------

def plot_schuster_spectrum(results: dict, out_path: Path) -> None:
    """3-panel stacked Schuster spectrum plot.

    Parameters
    ----------
    results : dict
        Full results JSON loaded as dict.
    out_path : Path
        Output PNG path.
    """
    schuster = results["schuster"]
    fig, axes = plt.subplots(3, 1, figsize=(12, 13), dpi=300)
    fig.suptitle("Cluster-Robust Schuster Power Spectrum", fontsize=14, fontweight="bold", y=1.01)

    for ax, key in zip(axes, CATALOG_KEYS):
        cat_data = schuster[key]
        spectrum = cat_data["spectrum"]
        n_events = cat_data["n_events"]
        n_clusters = cat_data.get("n_clusters_at_annual", "?")

        periods = np.array([e["period_days"] for e in spectrum])
        p_std = np.array([e["p_standard"] for e in spectrum])
        p_cr = np.array([e["p_cluster_robust"] for e in spectrum])

        # Clamp to avoid log(0)
        p_std_plot = np.clip(p_std, 1e-20, 1.0)
        p_cr_plot = np.clip(p_cr, 1e-20, 1.0)

        # Standard p as thin gray line
        ax.plot(periods, p_std_plot, color="gray", linewidth=0.8, alpha=0.6, label="Standard p-value")
        # Cluster-robust as thick steelblue
        ax.plot(periods, p_cr_plot, color="steelblue", linewidth=1.8, label="Cluster-robust p-value")

        # Significance thresholds
        ax.axhline(0.05, color="black", linestyle="--", linewidth=0.9, alpha=0.7, label="p=0.05")
        ax.axhline(0.001, color="black", linestyle="--", linewidth=0.6, alpha=0.5, label="p=0.001")

        # Annual / sub-annual period markers
        for akey, (T, color, lbl) in ANNUAL_PERIODS.items():
            ax.axvline(T, color=color, linestyle=":", linewidth=1.2, alpha=0.85, label=lbl)

        # Tidal period markers
        for T in TIDAL_PERIODS:
            ax.axvline(T, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.25, 548)
        ax.set_ylim(1e-20, 2.0)
        ax.set_ylabel("p-value (log scale)", fontsize=9)

        # x-tick labels at named periods
        xtick_vals = list(ANNOTATION_PERIODS.values())
        xtick_labs = list(ANNOTATION_PERIODS.keys())
        ax.set_xticks(xtick_vals)
        ax.set_xticklabels(xtick_labs, fontsize=8)

        title_n = f"{CATALOG_TITLES[key]} — n={n_events:,}, n_clusters={n_clusters:,}" if isinstance(n_clusters, int) else f"{CATALOG_TITLES[key]} — n={n_events:,}"
        ax.set_title(title_n, fontsize=10)

        if ax is axes[0]:
            handles, labels_list = ax.get_legend_handles_labels()
            # Deduplicate
            seen = {}
            for h, l in zip(handles, labels_list):
                if l not in seen:
                    seen[l] = h
            ax.legend(seen.values(), seen.keys(), fontsize=7, loc="upper left", ncol=3)

    axes[-1].set_xlabel("Period (days)", fontsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Saved Schuster spectrum to %s", out_path)


# ---------------------------------------------------------------------------
# Figure 2 — MFPA periodogram
# ---------------------------------------------------------------------------

def plot_mfpa_scan(results: dict, out_path: Path) -> None:
    """3-panel stacked MFPA periodogram plot.

    Parameters
    ----------
    results : dict
        Full results JSON loaded as dict.
    out_path : Path
        Output PNG path.
    """
    mfpa = results["mfpa"]
    fig, axes = plt.subplots(3, 1, figsize=(12, 13), dpi=300)
    fig.suptitle("MFPA Periodogram (6 hours – 18 months)", fontsize=14, fontweight="bold", y=1.01)

    for ax, key in zip(axes, CATALOG_KEYS):
        cat_data = mfpa[key]
        spectrum = cat_data["spectrum"]
        sig_periods = cat_data["significant_periods"]

        periods = np.array([e["period_days"] for e in spectrum])
        power = np.array([e["power"] for e in spectrum])
        p95 = np.array([e["p95_threshold"] for e in spectrum])
        p99 = np.array([e["p99_threshold"] for e in spectrum])

        # Power curve
        ax.plot(periods, power, color="steelblue", linewidth=1.2, label="MFPA power")
        # Significance thresholds as envelope curves
        ax.plot(periods, p95, color="darkorange", linestyle="--", linewidth=0.9, label="95th percentile threshold")
        ax.plot(periods, p99, color="red", linestyle="--", linewidth=0.9, label="99th percentile threshold")

        # Mark significant peaks with filled orange triangles
        for sp in sig_periods:
            T = sp["period_days"]
            P = sp["power"]
            ax.plot(T, P, marker="^", color="darkorange", markersize=8, zorder=5)
            ax.annotate(
                f"{T:.1f}d",
                xy=(T, P),
                xytext=(0, 8),
                textcoords="offset points",
                fontsize=7,
                ha="center",
                color="darkorange",
            )

        # Annual / sub-annual period markers
        for akey, (T, color, lbl) in ANNUAL_PERIODS.items():
            ax.axvline(T, color=color, linestyle=":", linewidth=1.0, alpha=0.7)

        ax.set_xscale("log")
        ax.set_xlim(0.25, 548)
        ax.set_ylabel("MFPA Power", fontsize=9)

        xtick_vals = list(ANNOTATION_PERIODS.values())
        xtick_labs = list(ANNOTATION_PERIODS.keys())
        ax.set_xticks(xtick_vals)
        ax.set_xticklabels(xtick_labs, fontsize=8)

        n_sig = len(sig_periods)
        ax.set_title(f"{CATALOG_TITLES[key]} — {n_sig} significant period(s) >p95", fontsize=10)

        if ax is axes[0]:
            ax.legend(fontsize=8, loc="upper left")

    axes[-1].set_xlabel("Period (days)", fontsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Saved MFPA scan to %s", out_path)


# ---------------------------------------------------------------------------
# Figure 3 — Harmonic-interval overlay
# ---------------------------------------------------------------------------

def _predicted_phases(period_days: float) -> List[float]:
    """Return predicted peak solar-phase fractions for a given period."""
    year_days = 365.25
    n_harmonics = max(1, int(year_days / period_days))
    return [(k * period_days / year_days) % 1.0 for k in range(n_harmonics + 1)]


def plot_harmonic_intervals(results: dict, out_path: Path) -> None:
    """Single-panel harmonic-interval overlay figure.

    Parameters
    ----------
    results : dict
        Full results JSON loaded as dict.
    out_path : Path
        Output PNG path.
    """
    # Collect all significant MFPA periods across all catalogs
    all_sig_periods_raw = results["mfpa"]["raw"]["significant_periods"]
    all_sig_periods_gk = results["mfpa"]["gk_mainshocks"]["significant_periods"]
    all_sig_periods_a1b = results["mfpa"]["a1b_mainshocks"]["significant_periods"]

    # Combine and deduplicate by period (round to 2 decimal places)
    seen_periods = {}
    for sp in (all_sig_periods_raw + all_sig_periods_gk + all_sig_periods_a1b):
        key = round(sp["period_days"], 1)
        if key not in seen_periods:
            seen_periods[key] = sp["period_days"]
    unique_periods = sorted(seen_periods.values())

    # Color palette for detected periods
    period_colors = plt.cm.tab10(np.linspace(0, 0.9, max(len(unique_periods), 1)))

    fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Solar Phase Fraction (0 = Jan 1)", fontsize=11)
    ax.set_title("A1b Phase Intervals vs. Detected MFPA Harmonic Phases", fontsize=12, fontweight="bold")

    # Draw A1b interval bands
    for lo, hi, lbl in A1B_INTERVALS:
        ax.axvspan(lo, hi, alpha=0.18, color="gray", zorder=1)
        ax.text(
            (lo + hi) / 2,
            0.85,
            lbl,
            ha="center",
            va="center",
            fontsize=8.5,
            color="dimgray",
            fontstyle="italic",
        )

    # Draw predicted peak phases for each significant period
    legend_handles = []
    if len(unique_periods) == 0:
        # No significant periods — add informational text
        ax.text(
            0.5,
            0.5,
            "No MFPA significant periods detected above p95 threshold",
            ha="center",
            va="center",
            fontsize=10,
            color="gray",
            transform=ax.transAxes,
        )
    else:
        for i, T in enumerate(unique_periods):
            color = period_colors[i]
            phases = _predicted_phases(T)
            for phase in phases:
                ax.axvline(phase, color=color, linewidth=1.5, alpha=0.75, zorder=3)
            label = f"{T:.1f}d harmonic"
            legend_handles.append(mpatches.Patch(color=color, label=label))

    # A1b interval legend entry
    legend_handles.append(mpatches.Patch(color="gray", alpha=0.35, label="A1b baseline interval"))

    ax.legend(handles=legend_handles, fontsize=8, loc="upper right", ncol=2)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=9)

    # Add equinox / solstice markers
    equinox_phases = {
        "Mar equinox (~0.19)": 0.19,
        "Jun solstice (~0.45)": 0.45,
        "Sep equinox (~0.72)": 0.72,
        "Dec solstice (~0.95)": 0.95,
    }
    for lbl, ph in equinox_phases.items():
        ax.axvline(ph, color="black", linewidth=0.5, linestyle="--", alpha=0.3, zorder=2)
        ax.text(ph, 0.05, lbl.split("(")[0].strip(), ha="center", fontsize=6, color="black", alpha=0.5, rotation=90, va="bottom")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Saved harmonic-interval overlay to %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Generate all Case A1 figures from results JSON."""
    logger.info("=== Case A1 Visualization ===")

    with open(RESULTS_PATH) as fh:
        results = json.load(fh)

    plot_schuster_spectrum(results, OUTPUT_DIR / "case-a1-schuster-spectrum.png")
    plot_mfpa_scan(results, OUTPUT_DIR / "case-a1-mfpa-scan.png")
    plot_harmonic_intervals(results, OUTPUT_DIR / "case-a1-harmonic-intervals.png")

    logger.info("=== All figures written ===")


if __name__ == "__main__":
    main()
