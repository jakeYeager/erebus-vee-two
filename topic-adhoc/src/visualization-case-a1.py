"""
Case A1: Visualization — P-Value Sweep and Distribution Grid

Reads case-a1-results.json and produces:
  - case-a1-pvalue-sweep.png   (line chart, p-value vs bin count)
  - case-a1-distributions.png  (3x4 grid of bar charts)
"""

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = BASE_DIR.parent
DATA_PATH = PROJECT_ROOT / "data" / "global-sets" / "iscgem_global_events.csv"
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a1-results.json"

METRICS = ["solar_secs", "lunar_secs", "midnight_secs"]
BIN_COUNTS = [8, 16, 24, 32]
LUNAR_CYCLE_SECS = 29.53059 * 86400
DAY_SECS = 86400


def is_leap_year(year: int) -> bool:
    """Return True if *year* is a Gregorian leap year."""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def year_length_seconds(year: int) -> int:
    """Return the calendar-year length in seconds for *year*."""
    return 366 * DAY_SECS if is_leap_year(year) else 365 * DAY_SECS


def load_results() -> dict:
    """Load the results JSON."""
    with open(RESULTS_PATH) as f:
        return json.load(f)


def plot_pvalue_sweep(data: dict) -> None:
    """Create the p-value sweep line chart."""
    fig, ax = plt.subplots(figsize=(10, 6))

    styles = {
        "solar_secs": {"color": "steelblue", "linestyle": "-", "marker": "o", "label": "solar_secs"},
        "lunar_secs": {"color": "coral", "linestyle": "-", "marker": "s", "label": "lunar_secs"},
        "midnight_secs": {"color": "gray", "linestyle": "--", "marker": "^", "label": "midnight_secs"},
    }

    x_positions = list(range(len(BIN_COUNTS)))
    x_labels = [str(k) for k in BIN_COUNTS]

    for metric in METRICS:
        p_values = [data["results"][metric][str(k)]["p_value"] for k in BIN_COUNTS]
        style = styles[metric]
        ax.plot(
            x_positions, p_values,
            color=style["color"],
            linestyle=style["linestyle"],
            marker=style["marker"],
            markersize=8,
            linewidth=2,
            label=style["label"],
            zorder=3,
        )
        # Annotate each point
        for i, p in enumerate(p_values):
            label_text = f"{p:.2g}"
            ax.annotate(
                label_text,
                (x_positions[i], p),
                textcoords="offset points",
                xytext=(0, 12),
                ha="center",
                fontsize=8,
            )

    # Reference lines
    ax.axhline(0.05, color="gray", linestyle="--", linewidth=0.8, label=r"$\alpha$ = 0.05")
    ax.axhline(0.004167, color="red", linestyle="--", linewidth=1.0, label=r"Bonferroni $\alpha$")

    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel("Bin Count")
    ax.set_ylabel("p-value (log scale, inverted)")
    ax.set_title("Chi-Square p-value by Bin Count (Phase-Normalized)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a1-pvalue-sweep.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_distributions(data: dict) -> None:
    """Create the 3x4 distribution grid."""
    # Load raw data to recompute bins for plotting
    df = pd.read_csv(DATA_PATH)

    # Compute phases
    year_lengths = df["solaration_year"].apply(year_length_seconds).astype(float)
    solar_phase = np.clip((df["solar_secs"] / year_lengths).values, 0, 1.0 - 1e-9)
    lunar_phase = np.clip((df["lunar_secs"] / LUNAR_CYCLE_SECS).values, 0, 1.0 - 1e-9)
    midnight_phase = np.clip((df["midnight_secs"] / float(DAY_SECS)).values, 0, 1.0 - 1e-9)
    phases = {
        "solar_secs": solar_phase,
        "lunar_secs": lunar_phase,
        "midnight_secs": midnight_phase,
    }

    n = len(df)
    fig, axes = plt.subplots(3, 4, figsize=(16, 10))

    for row, metric in enumerate(METRICS):
        for col, k in enumerate(BIN_COUNTS):
            ax = axes[row, col]
            bin_indices = np.floor(phases[metric] * k).astype(int)
            counts = np.bincount(bin_indices, minlength=k)

            p_val = data["results"][metric][str(k)]["p_value"]

            ax.bar(range(1, k + 1), counts, color="steelblue", edgecolor="white", linewidth=0.3)
            expected = n / k
            ax.axhline(expected, color="red", linestyle="--", linewidth=1.0)

            ax.set_title(f"{metric} -- {k} bins (p={p_val:.2g})", fontsize=8)
            ax.set_xlabel("Bin", fontsize=7)
            if col == 0:
                ax.set_ylabel("Count", fontsize=7)
            ax.tick_params(labelsize=6)

    fig.suptitle("Phase-Normalized Distributions by Metric and Bin Count", fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = OUTPUT_DIR / "case-a1-distributions.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def main() -> None:
    """Generate all visualizations."""
    data = load_results()
    plot_pvalue_sweep(data)
    plot_distributions(data)


if __name__ == "__main__":
    main()
