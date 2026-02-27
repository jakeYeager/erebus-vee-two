"""Case A0 Visualization: Catalog Comparison Charts

Reads case-a0-results.json and produces two PNG images:
1. Magnitude distribution comparison (overlapping histograms)
2. ComCat events by decade and ID prefix source (stacked bar)
"""

import json
import os
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_PATH = os.path.join(BASE_DIR, "output", "case-a0-results.json")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")


def load_results() -> dict[str, Any]:
    """Load the results JSON."""
    with open(RESULTS_PATH, "r") as f:
        return json.load(f)


def plot_magnitude_distributions(results: dict[str, Any]) -> None:
    """Create overlaid magnitude distribution comparison."""
    bins_data = results["magnitude_bins"]
    bin_edges = bins_data["bins"]
    comcat_counts = bins_data["comcat_counts"]
    iscgem_counts = bins_data["iscgem_counts"]

    comcat_n = results["catalogs"]["comcat"]["event_count"]
    iscgem_n = results["catalogs"]["iscgem"]["event_count"]

    comcat_1d_pct = results["magnitude_precision"]["comcat"]["one_decimal"]["pct"]
    iscgem_2d_pct = results["magnitude_precision"]["iscgem"]["two_decimal"]["pct"]

    x = np.array(bin_edges)
    width = 0.04

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(x - width / 2, comcat_counts, width=width, color="steelblue",
           alpha=0.6, label=f"ComCat (n={comcat_n:,})")
    ax.bar(x + width / 2, iscgem_counts, width=width, color="coral",
           alpha=0.6, label=f"ISC-GEM (n={iscgem_n:,})")

    ax.set_xlabel("Magnitude", fontsize=12)
    ax.set_ylabel("Event Count", fontsize=12)
    ax.set_title("Magnitude Distribution: ComCat vs ISC-GEM", fontsize=14)
    ax.legend(fontsize=11)

    annotation = (
        f"ComCat: {comcat_1d_pct}% 1-decimal | "
        f"ISC-GEM: {iscgem_2d_pct}% 2-decimal"
    )
    ax.annotate(annotation, xy=(0.98, 0.95), xycoords="axes fraction",
                ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))

    ax.set_xlim(5.95, 9.65)
    ax.set_xticks(np.arange(6.0, 9.7, 0.5))

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "case-a0-magnitude-distributions.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_comcat_prefix_temporal(results: dict[str, Any]) -> None:
    """Create stacked bar chart of ComCat events by decade and ID prefix."""
    import csv

    # We need per-decade counts for ALL prefix groups, not just iscgem.
    # Reload ComCat data to compute this.
    comcat_path = os.path.join(os.path.dirname(BASE_DIR),
                               results["catalogs"]["comcat"]["file"])
    with open(comcat_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    decades = ["1950s", "1960s", "1970s", "1980s", "1990s", "2000s", "2010s", "2020s"]
    decade_starts = [1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020]
    groups = ["us_native", "iscgem", "other"]

    counts: dict[str, dict[str, int]] = {d: {g: 0 for g in groups} for d in decades}

    for r in rows:
        year = int(float(r["solaration_year"]))
        uid = r["usgs_id"]
        if uid.startswith("iscgem"):
            group = "iscgem"
        elif uid.startswith("us"):
            group = "us_native"
        else:
            group = "other"

        for i in range(len(decade_starts) - 1, -1, -1):
            if year >= decade_starts[i]:
                counts[decades[i]][group] += 1
                break

    # Also handle 1940s events (there might be 1949 records)
    # The spec says X-axis: 1950s through 2020s, so we include 1940s events
    # in the 1950s bucket if they exist, or skip. Actually spec says decades
    # start at 1950s. Let's check for 1949.
    for r in rows:
        year = int(float(r["solaration_year"]))
        if year < 1950:
            uid = r["usgs_id"]
            if uid.startswith("iscgem"):
                group = "iscgem"
            elif uid.startswith("us"):
                group = "us_native"
            else:
                group = "other"
            # Include in 1950s for display purposes
            counts["1950s"][group] += 1

    us_vals = [counts[d]["us_native"] for d in decades]
    isc_vals = [counts[d]["iscgem"] for d in decades]
    other_vals = [counts[d]["other"] for d in decades]

    x = np.arange(len(decades))

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(x, us_vals, color="steelblue", label="us_native")
    ax.bar(x, isc_vals, bottom=us_vals, color="coral", label="iscgem")
    bottoms = [u + i for u, i in zip(us_vals, isc_vals)]
    ax.bar(x, other_vals, bottom=bottoms, color="gray", label="other")

    ax.set_xlabel("Decade", fontsize=12)
    ax.set_ylabel("Event Count", fontsize=12)
    ax.set_title("ComCat Events by Decade and ID Prefix Source", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(decades)
    ax.legend(fontsize=11)

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "case-a0-comcat-prefix-temporal.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main() -> None:
    """Generate all visualizations."""
    results = load_results()
    plot_magnitude_distributions(results)
    plot_comcat_prefix_temporal(results)


if __name__ == "__main__":
    main()
