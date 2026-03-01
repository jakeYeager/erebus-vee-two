"""
Case B5 — Visualizations: Solar Declination Rate-of-Change vs. Position Test

Generates three figures:
  1. case-b5-binplots.png       — 2×2 panel: bin distributions for all four variables at k=24
  2. case-b5-variable-ranking.png — horizontal bar chart of Cramér's V by variable at k=24
  3. case-b5-a1b-alignment.png  — 3×3 grid: A1b interval expected values overlaid on distributions
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-b5-results.json"
OUTPUT_DIR = BASE_DIR / "output"
SOLAR_GEO_PATH = BASE_DIR.parent / "data" / "iscgem" / "celestial-geometry" / "solar_geometry_global.csv"

# ── Load data ──────────────────────────────────────────────────────────────────
with open(RESULTS_PATH) as fh:
    results = json.load(fh)

df = pd.read_csv(SOLAR_GEO_PATH)
n = len(df)

actual_ranges = results["actual_ranges"]
variable_stats = results["variable_stats"]
variable_ranking = results["variable_ranking"]
a1b_alignment = results["a1b_alignment"]

# ── Style constants ────────────────────────────────────────────────────────────
STEELBLUE = "steelblue"
EXPECTED_COLOR = "black"
THRESHOLD_COLOR = "darkorange"
A1B_BAND_COLOR = "#cccccc"
EQUINOX_COLOR = "#444444"
DPI = 300
K24 = 24


# ── Helper: annotation string ─────────────────────────────────────────────────

def stat_annotation(var: str, k_key: str = "k24") -> str:
    """Build annotation string for a panel."""
    s = variable_stats[var][k_key]
    chi2 = s["chi2"]
    p = s["p_chi2"]
    v = s["cramer_v"]
    if p == 0.0:
        p_str = "p < 1×10⁻³⁰⁰"
    elif p < 0.0001:
        p_str = f"p = {p:.2e}"
    else:
        p_str = f"p = {p:.4f}"
    return f"χ² = {chi2:.1f}\n{p_str}\nV = {v:.4f}"


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: Four-panel bin distributions at k=24
# ─────────────────────────────────────────────────────────────────────────────

def make_binplots() -> None:
    """Generate 2×2 bin distribution panel figure."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Bin Distributions at k=24: Solar Phase and Geometric Variables",
                 fontsize=14, fontweight="bold", y=1.01)

    panel_order = [
        ("solar_phase",       axes[0, 0]),
        ("solar_declination", axes[0, 1]),
        ("declination_rate",  axes[1, 0]),
        ("earth_sun_distance",axes[1, 1]),
    ]

    # A1b baseline intervals (from Case 3A / Adhoc A1b)
    a1b_intervals = [
        (0.1875, 0.2500, "I₁"),
        (0.5625, 0.7083, "I₂"),
        (0.8542, 0.9375, "I₃"),
    ]

    # Equinox/solstice phase positions
    equinox_phases = [0.0, 0.219, 0.469, 0.719]  # Jan 1, ~Mar 20, ~Jun 21, ~Sep 23
    equinox_labels = ["Jan 1", "Mar eq", "Jun sol", "Sep eq"]

    for var, ax in panel_order:
        counts = np.array(variable_stats[var]["k24"]["bin_counts"])
        k = len(counts)
        expected = n / k
        std_counts = np.std(counts)
        threshold = expected + std_counts
        bin_indices = np.arange(k)

        # Horizontal bars (bins on y-axis, counts on x-axis)
        ax.barh(bin_indices, counts, color=STEELBLUE, height=0.8, alpha=0.85)
        ax.axvline(expected, color=EXPECTED_COLOR, linewidth=1.2,
                   linestyle="--", label=f"Expected ({expected:.0f})")
        ax.axvline(threshold, color=THRESHOLD_COLOR, linewidth=1.0,
                   linestyle=":", label=f"Mean+1SD ({threshold:.0f})")

        # Panel-specific formatting
        if var == "solar_phase":
            # Add A1b interval gray bands
            for lo, hi, label in a1b_intervals:
                lo_bin = lo * k
                hi_bin = hi * k
                ax.axhspan(lo_bin, hi_bin, color=A1B_BAND_COLOR, alpha=0.5, zorder=0)
                ax.text(ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else expected * 1.5,
                        (lo_bin + hi_bin) / 2,
                        label, va="center", ha="right", fontsize=7, color="#555555")
            # Add equinox/solstice lines
            for ep in equinox_phases:
                ax.axhline(ep * k, color=EQUINOX_COLOR, linewidth=0.7,
                           linestyle="-.", alpha=0.5)

            ax.set_yticks(np.arange(0, k, 4))
            ax.set_yticklabels([f"{i/k:.2f}" for i in range(0, k, 4)], fontsize=7)
            ax.set_ylabel("Solar Phase (Julian year fraction)", fontsize=9)
            ax.set_title("Solar Phase (solar_secs / Julian year)", fontsize=10, fontweight="bold")

        elif var == "solar_declination":
            dec_min = actual_ranges["solar_declination"]["min"]
            dec_max = actual_ranges["solar_declination"]["max"]
            tick_positions = [0, k // 4, k // 2, 3 * k // 4, k - 1]
            tick_values = [dec_min + t / (k - 1) * (dec_max - dec_min) for t in tick_positions]
            ax.set_yticks(tick_positions)
            ax.set_yticklabels([f"{v:.1f}°" for v in tick_values], fontsize=7)
            ax.set_ylabel("Solar Declination (degrees)", fontsize=9)
            ax.set_title("Solar Declination Angle", fontsize=10, fontweight="bold")

        elif var == "declination_rate":
            rate_min = actual_ranges["declination_rate"]["min"]
            rate_max = actual_ranges["declination_rate"]["max"]
            tick_positions = [0, k // 4, k // 2, 3 * k // 4, k - 1]
            tick_values = [rate_min + t / (k - 1) * (rate_max - rate_min) for t in tick_positions]
            ax.set_yticks(tick_positions)
            ax.set_yticklabels([f"{v:.3f}" for v in tick_values], fontsize=7)
            ax.set_ylabel("Declination Rate (deg/day)", fontsize=9)
            ax.set_title("Rate of Change of Solar Declination", fontsize=10, fontweight="bold")

        elif var == "earth_sun_distance":
            dist_min = actual_ranges["earth_sun_distance"]["min"]
            dist_max = actual_ranges["earth_sun_distance"]["max"]
            tick_positions = [0, k // 4, k // 2, 3 * k // 4, k - 1]
            tick_values = [dist_min + t / (k - 1) * (dist_max - dist_min) for t in tick_positions]
            ax.set_yticks(tick_positions)
            ax.set_yticklabels([f"{v:.4f}" for v in tick_values], fontsize=7)
            ax.set_ylabel("Earth-Sun Distance (AU)", fontsize=9)
            ax.set_title("Earth-Sun Distance", fontsize=10, fontweight="bold")

        ax.set_xlabel("Event Count", fontsize=9)

        # Annotation
        annot = stat_annotation(var)
        ax.text(0.98, 0.02, annot, transform=ax.transAxes,
                fontsize=7.5, va="bottom", ha="right",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#aaaaaa", alpha=0.85))

        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(axis="x", linestyle=":", alpha=0.4)

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-b5-binplots.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Variable ranking bar chart
# ─────────────────────────────────────────────────────────────────────────────

def make_variable_ranking() -> None:
    """Generate horizontal bar chart of Cramér's V by variable."""
    var_names_raw = [r["variable"] for r in variable_ranking]
    cramer_vs = [r["cramer_v_k24"] for r in variable_ranking]
    p_vals = [r["p_chi2_k24"] for r in variable_ranking]

    # Readable labels
    label_map = {
        "declination_rate": "Declination Rate\n(deg/day)",
        "earth_sun_distance": "Earth-Sun Distance\n(AU)",
        "solar_declination": "Solar Declination\n(degrees)",
        "solar_phase": "Solar Phase\n(annual cycle)",
    }
    labels = [label_map.get(v, v) for v in var_names_raw]

    # Reverse for top-to-bottom ordering
    labels = labels[::-1]
    cramer_vs_rev = cramer_vs[::-1]
    p_vals_rev = p_vals[::-1]
    var_names_rev = var_names_raw[::-1]

    fig, ax = plt.subplots(figsize=(9, 5))
    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, cramer_vs_rev, color=STEELBLUE, height=0.55)

    # Significance annotations
    for i, (p, v_val, bar) in enumerate(zip(p_vals_rev, cramer_vs_rev, bars)):
        sig = "*" if p < 0.05 else "ns"
        sig_color = "darkgreen" if p < 0.05 else "dimgray"
        ax.text(v_val + 0.001, i, f"  {sig}  V={v_val:.4f}",
                va="center", ha="left", fontsize=9, color=sig_color, fontweight="bold")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Cramér's V (k=24)", fontsize=11)
    ax.set_title("Cramér's V by Solar Variable at k=24", fontsize=13, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.4)

    # Legend
    sig_patch = mpatches.Patch(color="darkgreen", label="* significant (p < 0.05)")
    ns_patch = mpatches.Patch(color="dimgray", label="ns = not significant")
    ax.legend(handles=[sig_patch, ns_patch], fontsize=9, loc="lower right")

    # Add note about comparison
    ax.text(0.02, 0.02,
            "Note: solar_phase is cyclic (annual).\nGeometric variables are non-cyclic.\nDirect V comparison requires caution.",
            transform=ax.transAxes, fontsize=7.5, va="bottom", ha="left",
            color="#666666",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f8f8f8",
                      edgecolor="#cccccc", alpha=0.9))

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-b5-variable-ranking.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: A1b interval alignment (3 rows × 3 columns)
# ─────────────────────────────────────────────────────────────────────────────

def make_a1b_alignment() -> None:
    """Generate 3×3 grid showing A1b interval expected values on distributions."""
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    fig.suptitle("A1b Interval Alignment: Expected Physical Values on Geometric Distributions",
                 fontsize=13, fontweight="bold")

    # Column variables (non-cyclic only)
    col_vars = [
        ("solar_declination",  "Solar Declination", "degrees",
         actual_ranges["solar_declination"]["min"],
         actual_ranges["solar_declination"]["max"]),
        ("declination_rate",   "Declination Rate",  "deg/day",
         actual_ranges["declination_rate"]["min"],
         actual_ranges["declination_rate"]["max"]),
        ("earth_sun_distance", "Earth-Sun Distance", "AU",
         actual_ranges["earth_sun_distance"]["min"],
         actual_ranges["earth_sun_distance"]["max"]),
    ]

    # Row intervals
    row_intervals = [
        ("interval_1_phase_0.22", "Interval 1 (phase ≈ 0.22, ~Mar 22)",
         "solar_declination_expected", "declination_rate_expected", "earth_sun_distance_expected"),
        ("interval_2_phase_0.64", "Interval 2 (phase ≈ 0.64, ~Aug 22)",
         "solar_declination_expected", "declination_rate_expected", "earth_sun_distance_expected"),
        ("interval_3_phase_0.90", "Interval 3 (phase ≈ 0.90, ~Nov 24)",
         "solar_declination_expected", "declination_rate_expected", "earth_sun_distance_expected"),
    ]

    expected_keys = [
        ("solar_declination_expected", "declination_rate_expected", "earth_sun_distance_expected")
    ]

    # Build bin count arrays for each variable at k=24
    bin_data = {}
    for var, _, _, vmin, vmax in col_vars:
        counts = np.array(variable_stats[var]["k24"]["bin_counts"])
        bin_edges = np.linspace(vmin, vmax, K24 + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
        bin_data[var] = (counts, bin_centers, vmin, vmax)

    for row_i, (interval_key, row_label, exp_dec_key, exp_rate_key, exp_dist_key) in enumerate(row_intervals):
        alignment = a1b_alignment[interval_key]
        exp_values = {
            "solar_declination": alignment["solar_declination_expected"],
            "declination_rate": alignment["declination_rate_expected"],
            "earth_sun_distance": alignment["earth_sun_distance_expected"],
        }

        for col_i, (var, col_title, units, vmin, vmax) in enumerate(col_vars):
            ax = axes[row_i, col_i]
            counts, bin_centers, vmin_, vmax_ = bin_data[var]
            exp_val = exp_values[var]

            # Histogram as horizontal bar chart
            bar_width = (vmax_ - vmin_) / K24
            ax.bar(bin_centers, counts, width=bar_width * 0.85,
                   color=STEELBLUE, alpha=0.75)

            # Expected value vertical line
            ax.axvline(exp_val, color="red", linewidth=1.8,
                       linestyle="-", label=f"A1b expected: {exp_val:.3g} {units}")

            # Mean line
            expected_count = n / K24
            ax.axhline(expected_count, color="black", linewidth=0.8,
                       linestyle="--", alpha=0.7, label=f"Expected n ({expected_count:.0f})")

            # Column title (top row only)
            if row_i == 0:
                ax.set_title(f"{col_title}\n({units})", fontsize=9, fontweight="bold")

            # Row label (leftmost column only)
            if col_i == 0:
                ax.set_ylabel(row_label + "\nEvent Count", fontsize=7.5, labelpad=4)
            else:
                ax.set_ylabel("Event Count", fontsize=8)

            ax.set_xlabel(f"{units}", fontsize=8)
            ax.legend(fontsize=6.5, loc="upper right")
            ax.grid(axis="y", linestyle=":", alpha=0.4)
            ax.tick_params(labelsize=7)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out_path = OUTPUT_DIR / "case-b5-a1b-alignment.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# ── Entry point ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    make_binplots()
    make_variable_ranking()
    make_a1b_alignment()
    print("All figures generated.")
