"""
Case A3.C2: Visualization Script

Generates three figures from the case-a3-c2-results.json output:
  1. case-a3-c2-degradation.png  — chi-square p-value and Cramér's V trajectories
  2. case-a3-c2-interval-decay.png — per-interval z-score decay across removal steps
  3. case-a3-c2-sequence-summary.png — visual table of sequence metrics for major events
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a3-c2-results.json"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Load results
# ---------------------------------------------------------------------------

with open(RESULTS_PATH) as f:
    results = json.load(f)

major_events = results["major_events"]
runs = results["runs"]
seq_metrics = results["sequence_metrics"]

# ---------------------------------------------------------------------------
# Build x-axis tick labels
# ---------------------------------------------------------------------------


def build_tick_labels(major_events: list[dict]) -> list[str]:
    """Build per-step x-axis tick labels (step 0 = 'Baseline', steps 1+ = 'MX.X\\nYYYY').

    Args:
        major_events: List of major event dicts from results JSON.

    Returns:
        List of label strings, one per removal step including baseline.
    """
    labels = ["Baseline"]
    for ev in major_events:
        mag = ev["usgs_mag"]
        year = ev["event_at"][:4]
        labels.append(f"M{mag:.1f}\n{year}")
    return labels


# ---------------------------------------------------------------------------
# Run styling
# ---------------------------------------------------------------------------

RUN_STYLES: dict[str, dict] = {
    "raw_gk": {"color": "steelblue", "linestyle": "-", "label": "Raw + G-K removal"},
    "raw_reas": {"color": "steelblue", "linestyle": "--", "label": "Raw + Reasenberg removal"},
    "mainshock_gk": {"color": "red", "linestyle": "-", "label": "G-K mainshocks only"},
    "mainshock_reas": {"color": "red", "linestyle": "--", "label": "Reasenberg mainshocks only"},
}

RUN_KEYS = ["raw_gk", "raw_reas", "mainshock_gk", "mainshock_reas"]


# ---------------------------------------------------------------------------
# Figure 1: Chi-square degradation trajectory
# ---------------------------------------------------------------------------


def plot_degradation() -> None:
    """Plot chi-square p-value and Cramér's V across removal steps for all four runs."""
    tick_labels = build_tick_labels(major_events)
    x_steps = list(range(len(tick_labels)))

    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(2, 1, figure=fig, hspace=0.15)
    ax_p = fig.add_subplot(gs[0])
    ax_v = fig.add_subplot(gs[1], sharex=ax_p)

    for run_key in RUN_KEYS:
        style = RUN_STYLES[run_key]
        steps = runs[run_key]["steps"]
        p_vals = [s["stats"]["p_chi2_k24"] for s in steps]
        v_vals = [s["stats"]["cramers_v"] for s in steps]

        ax_p.plot(
            x_steps,
            p_vals,
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            marker="o",
            markersize=4,
            label=style["label"],
        )
        ax_v.plot(
            x_steps,
            v_vals,
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            marker="o",
            markersize=4,
            label=style["label"],
        )

    # P-value panel
    ax_p.axhline(0.05, color="gray", linestyle="--", linewidth=1.0, label="p = 0.05 threshold")
    ax_p.set_yscale("log")
    ax_p.set_ylim(1e-12, 1.0)
    ax_p.set_ylabel("Chi-Square p-value (log scale)", fontsize=11)
    ax_p.legend(fontsize=9, loc="upper left")
    ax_p.grid(True, which="both", linestyle=":", alpha=0.5)
    ax_p.set_title(
        "Phased Removal of M≥8.5 Sequences — Chi-Square and Cramér's V Degradation",
        fontsize=13,
        pad=12,
    )
    plt.setp(ax_p.get_xticklabels(), visible=False)

    # Cramér's V panel
    ax_v.set_ylabel("Cramér's V", fontsize=11)
    ax_v.grid(True, linestyle=":", alpha=0.5)
    ax_v.set_xticks(x_steps)
    ax_v.set_xticklabels(tick_labels, fontsize=8)
    ax_v.set_xlabel("Event Removed (Removal Order, Largest First)", fontsize=11)

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-c2-degradation.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Interval z-score decay
# ---------------------------------------------------------------------------


def plot_interval_decay() -> None:
    """Plot A1b interval z-score trajectories across sequential removal steps."""
    tick_labels = build_tick_labels(major_events)
    x_steps = list(range(len(tick_labels)))

    interval_keys = ["interval_1_z", "interval_2_z", "interval_3_z"]
    interval_labels = [
        "Interval 1 (0.167–0.250)",
        "Interval 2 (0.625–0.667)",
        "Interval 3 (0.875–0.917)",
    ]

    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
    fig.suptitle(
        "A1b Interval Elevation Z-Score Through Sequential Removal",
        fontsize=13,
        y=1.01,
    )

    for row_idx, (interval_key, row_label) in enumerate(zip(interval_keys, interval_labels)):
        ax = axes[row_idx]

        for run_key in RUN_KEYS:
            style = RUN_STYLES[run_key]
            steps = runs[run_key]["steps"]
            z_vals = [s["stats"][interval_key] for s in steps]

            ax.plot(
                x_steps,
                z_vals,
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=1.8,
                marker="o",
                markersize=4,
                label=style["label"],
            )

        ax.axhline(1.96, color="gray", linestyle="--", linewidth=1.0, label="z = 1.96 (α=0.05)")
        ax.set_ylim(bottom=0)
        ax.set_ylabel(f"Z-score\n{row_label}", fontsize=9)
        ax.grid(True, linestyle=":", alpha=0.5)
        if row_idx == 0:
            ax.legend(fontsize=8, loc="upper right")

    axes[-1].set_xticks(x_steps)
    axes[-1].set_xticklabels(tick_labels, fontsize=8)
    axes[-1].set_xlabel("Event Removed (Removal Order, Largest First)", fontsize=11)

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-c2-interval-decay.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Sequence summary visual table
# ---------------------------------------------------------------------------

# Canonical names for each major event
CANONICAL_NAMES: dict[str, str] = {
    "iscgem879136": "1960 Valdivia",
    "iscgem7453151": "2004 Sumatra",
    "iscgem869809": "1964 Alaska",
    "iscgem16461282": "2011 Tōhoku",
    "iscgem893648": "1952 Kamchatka",
    "iscgem14340585": "2010 Maule",
    "iscgem859206": "1965 Rat Islands",
    "iscgem895681": "1950 Assam",
    "iscgem7486110": "2005 Nias",
    "iscgem879134": "1960 Valdivia-2",
    "iscgem600860404": "2012 Sumatra OT",
    "iscgem886030": "1957 Andreanof",
    "iscgem873239": "1963 Kuril",
}


def fmt_pct(val: float | None) -> str:
    """Format early_pct as percentage string or N/A.

    Args:
        val: Float value in [0,1] or None.

    Returns:
        Formatted string like '75.0%' or 'N/A'.
    """
    if val is None:
        return "N/A"
    return f"{val * 100:.0f}%"


def fmt_days(val: float) -> str:
    """Format window_days as string.

    Args:
        val: Number of days.

    Returns:
        Formatted string like '1825.3 d' or '—'.
    """
    if val == 0:
        return "—"
    return f"{val:.0f} d"


def plot_sequence_summary() -> None:
    """Render a matplotlib table of sequence metrics for all M>=8.5 events."""
    gk_metrics = {m["usgs_id"]: m for m in seq_metrics["gk"]}
    reas_metrics = {m["usgs_id"]: m for m in seq_metrics["reasenberg"]}

    # Build table rows: sorted by magnitude descending (major_events already sorted)
    col_headers = [
        "Event", "Mag", "Date",
        "G-K\nAftershocks", "G-K\nWindow", "G-K\nEarly%", "G-K\nClass",
        "Reas\nAftershocks", "Reas\nWindow", "Reas\nEarly%", "Reas\nClass",
    ]

    rows = []
    for ev in major_events:
        uid = ev["usgs_id"]
        gk = gk_metrics[uid]
        reas = reas_metrics[uid]
        name = CANONICAL_NAMES.get(uid, uid)
        date_str = ev["event_at"][:10]
        row = [
            name,
            f"{ev['usgs_mag']:.2f}",
            date_str,
            str(gk["aftershock_count"]),
            fmt_days(gk["window_days"]),
            fmt_pct(gk["early_pct"]),
            gk["classification"][:4],
            str(reas["aftershock_count"]),
            fmt_days(reas["window_days"]),
            fmt_pct(reas["early_pct"]),
            reas["classification"][:4],
        ]
        rows.append(row)

    fig, ax = plt.subplots(figsize=(16, 7))
    ax.set_axis_off()
    ax.set_title(
        "Sequence Metrics for M≥8.5 Events — G-K vs. Reasenberg Attribution",
        fontsize=12,
        pad=14,
        fontweight="bold",
    )

    table = ax.table(
        cellText=rows,
        colLabels=col_headers,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.7)

    # Style header row
    n_cols = len(col_headers)
    for col_idx in range(n_cols):
        cell = table[0, col_idx]
        cell.set_facecolor("#2c4a7c")
        cell.set_text_props(color="white", fontweight="bold")

    # Alternate row colors; highlight rows where G-K and Reasenberg classification differ
    for row_idx, ev in enumerate(major_events):
        uid = ev["usgs_id"]
        gk = gk_metrics[uid]
        reas = reas_metrics[uid]
        classifications_differ = gk["classification"] != reas["classification"]

        for col_idx in range(n_cols):
            cell = table[row_idx + 1, col_idx]
            if classifications_differ:
                cell.set_facecolor("#ffe8a0")  # Amber highlight for differing classification
            elif row_idx % 2 == 0:
                cell.set_facecolor("#f2f2f2")
            else:
                cell.set_facecolor("white")

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-c2-sequence-summary.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Generate all three case A3.C2 figures."""
    print("Generating Figure 1: chi-square degradation trajectory...")
    plot_degradation()

    print("Generating Figure 2: interval z-score decay...")
    plot_interval_decay()

    print("Generating Figure 3: sequence summary table...")
    plot_sequence_summary()

    print("All figures generated successfully.")


if __name__ == "__main__":
    main()
