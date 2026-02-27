"""
Case A0b Visualization: Duplicate Detection and Cross-Catalog Event Accounting

Produces three PNG images from case-a0b-results.json:
  1. ID overlap stacked bar chart
  2. Cross-catalog event accounting grouped bar chart
  3. Unmatched events temporal side-by-side bar chart
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a0b-results.json"


def load_results() -> Dict[str, Any]:
    """Load the results JSON."""
    with open(RESULTS_PATH) as f:
        return json.load(f)


def plot_id_overlap(results: Dict[str, Any]) -> None:
    """Image 1: Stacked bar chart of ID-level accounting."""
    id_ref = results["id_cross_reference"]

    comcat_matched = id_ref["id_matched_in_iscgem"]["count"]
    comcat_unmatched = id_ref["id_unmatched_in_iscgem"]["count"]
    comcat_total = id_ref["comcat_iscgem_prefix_count"]

    iscgem_matched = comcat_matched  # same set of IDs
    iscgem_unmatched = id_ref["iscgem_ids_not_in_comcat"]["count"]
    iscgem_total = 9210

    fig, ax = plt.subplots(figsize=(10, 6))

    categories = [
        f"ComCat iscgem-prefix\n(n={comcat_total:,})",
        f"ISC-GEM full catalog\n(n={iscgem_total:,})",
    ]
    matched_vals = [comcat_matched, iscgem_matched]
    unmatched_vals = [comcat_unmatched, iscgem_unmatched]

    x = np.arange(len(categories))
    width = 0.5

    bars_matched = ax.bar(x, matched_vals, width, label="Matched IDs", color="steelblue")
    bars_unmatched = ax.bar(x, unmatched_vals, width, bottom=matched_vals,
                            label="Unmatched IDs", color="lightgray", edgecolor="gray")

    # Annotate
    for i, (m, u) in enumerate(zip(matched_vals, unmatched_vals)):
        ax.text(i, m / 2, f"{m:,}", ha="center", va="center", fontweight="bold", fontsize=11)
        ax.text(i, m + u / 2, f"{u:,}", ha="center", va="center", fontsize=11, color="dimgray")

    ax.set_ylabel("Event Count")
    ax.set_title("Direct ID Overlap: ComCat iscgem-prefixed Records vs ISC-GEM Catalog")
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a0b-id-overlap.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Saved %s", out_path)


def plot_event_accounting(results: Dict[str, Any]) -> None:
    """Image 2: Grouped bar chart of cross-catalog event accounting."""
    primary = results["cross_catalog_matching"]["primary"]
    strict = results["cross_catalog_matching"]["strict"]

    categories = ["ComCat Only", "Matched (both)", "ISC-GEM Only"]
    primary_vals = [primary["comcat_only"], primary["matched"], primary["iscgem_only"]]
    strict_vals = [strict["comcat_only"], strict["matched"], strict["iscgem_only"]]

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(categories))
    width = 0.35

    bars1 = ax.bar(x - width / 2, primary_vals, width, label="Primary (60s/50km/0.3)", color="steelblue")
    bars2 = ax.bar(x + width / 2, strict_vals, width, label="Strict (30s/25km/0.2)", color="lightsteelblue",
                   edgecolor="steelblue")

    # Annotate
    for bar_group in [bars1, bars2]:
        for bar in bar_group:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, height + 20,
                    f"{int(height):,}", ha="center", va="bottom", fontsize=9)

    ax.set_ylabel("Event Count")
    ax.set_title("Cross-Catalog Event Accounting: ComCat vs ISC-GEM")
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a0b-event-accounting.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Saved %s", out_path)


def plot_unmatched_temporal(results: Dict[str, Any]) -> None:
    """Image 3: Side-by-side bar chart of unmatched events by decade."""
    primary = results["cross_catalog_matching"]["primary"]
    comcat_temporal = primary["comcat_only_temporal"]
    iscgem_temporal = primary["iscgem_only_temporal"]

    # Collect all decades and sort
    all_decades = sorted(set(list(comcat_temporal.keys()) + list(iscgem_temporal.keys())))

    comcat_vals = [comcat_temporal.get(d, 0) for d in all_decades]
    iscgem_vals = [iscgem_temporal.get(d, 0) for d in all_decades]

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(all_decades))
    width = 0.35

    bars1 = ax.bar(x - width / 2, comcat_vals, width, label="ComCat Only", color="steelblue")
    bars2 = ax.bar(x + width / 2, iscgem_vals, width, label="ISC-GEM Only", color="coral")

    # Annotate nonzero bars
    for bar_group in [bars1, bars2]:
        for bar in bar_group:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, height + 2,
                        f"{int(height)}", ha="center", va="bottom", fontsize=8)

    ax.set_ylabel("Event Count")
    ax.set_xlabel("Decade")
    ax.set_title("Unmatched Events by Decade (Primary Tolerance)")
    ax.set_xticks(x)
    ax.set_xticklabels(all_decades)
    ax.legend()

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-a0b-unmatched-temporal.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Saved %s", out_path)


def main() -> None:
    results = load_results()
    plot_id_overlap(results)
    plot_event_accounting(results)
    plot_unmatched_temporal(results)
    logger.info("All visualizations complete.")


if __name__ == "__main__":
    main()
