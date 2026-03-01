"""Case B3: Tectonic Regime Stratification — Visualization Script.

Generates three PNG figures:
  1. case-b3-binplots.png         — 3-panel bin distributions (thrust, normal,
                                    strike-slip) at k=24
  2. case-b3-mechanism-comparison.png — Cramér's V and Rayleigh R by mechanism
                                        (dual-axis grouped bar chart)
  3. case-b3-coverage-map.png     — global event map colored by mechanism class

cartopy is unavailable; global map uses lat/lon scatter with ocean background.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent  # topic-a2/
DATA_DIR = BASE_DIR.parent / "data" / "iscgem"
FOCAL_DIR = DATA_DIR / "focal-mechanism"
OUTPUT_DIR = BASE_DIR / "output"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("visualization-case-b3")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DPI = 300
SOLAR_YEAR_SECS: float = 31_557_600.0
K = 24

# Mechanism display labels
MECH_LABELS: dict[str, str] = {
    "thrust":      "Thrust",
    "normal":      "Normal",
    "strike_slip": "Strike-Slip",
}
MECHANISM_ORDER: list[str] = ["thrust", "normal", "strike_slip"]

# Colors per mechanism class
MECH_COLORS: dict[str, str] = {
    "thrust":      "#d62728",   # red
    "normal":      "#1f77b4",   # blue
    "strike_slip": "#2ca02c",   # green
    "oblique":     "#ff7f0e",   # orange
    "unmatched":   "#aaaaaa",   # gray
}

# A1b baseline elevated phase intervals
A1B_INTERVALS: list[tuple[float, float]] = [
    (0.1875, 0.25),
    (0.625,  0.656),
    (0.875,  0.917),
]

# Solar calendar reference positions (fraction of year)
EQUINOX_SPRING  = 0.19
SOLSTICE_SUMMER = 0.44
EQUINOX_AUTUMN  = 0.69
SOLSTICE_WINTER = 0.94


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_data() -> tuple[pd.DataFrame, dict]:
    """Load focal join catalog and case-b3-results.json.

    Returns:
        Tuple of (focal join DataFrame with tectonic_class and solar_phase,
                  results dict from JSON).
    """
    # Load focal join
    focal_path = FOCAL_DIR / "focal_join_global.csv"
    df = pd.read_csv(focal_path)

    # Reproduce classify_rake inline for visualization (no import of analysis module)
    def classify_rake(rake: float) -> str:
        """Classify focal mechanism from rake angle."""
        while rake > 180:
            rake -= 360
        while rake <= -180:
            rake += 360
        if 45 <= rake <= 135:
            return "thrust"
        elif -135 <= rake <= -45:
            return "normal"
        elif (-45 < rake <= 45) or (135 < rake <= 180) or (-180 <= rake < -135):
            return "strike_slip"
        else:
            return "oblique"

    df["solar_phase"] = (df["solar_secs"] / SOLAR_YEAR_SECS) % 1.0
    df["tectonic_class"] = "unmatched"

    matched_mask = df["match_confidence"] == "proximity"
    has_mech = matched_mask & df["mechanism"].notna()
    for mech in ["thrust", "normal", "strike_slip"]:
        df.loc[has_mech & (df["mechanism"] == mech), "tectonic_class"] = mech
    df.loc[has_mech & (df["mechanism"] == "oblique"), "tectonic_class"] = "oblique"

    needs_fallback = matched_mask & df["mechanism"].isna() & df["rake"].notna()
    for idx in df.index[needs_fallback]:
        df.at[idx, "tectonic_class"] = classify_rake(df.at[idx, "rake"])

    # Load results JSON
    results_path = OUTPUT_DIR / "case-b3-results.json"
    with open(results_path) as fh:
        results = json.load(fh)

    logger.info("Data loaded: n=%d", len(df))
    return df, results


# ---------------------------------------------------------------------------
# Figure 1: Three-panel bin distributions
# ---------------------------------------------------------------------------
def _compute_threshold(n: int, k: int) -> float:
    """Compute 1-SD elevated-bin threshold.

    Args:
        n: Total event count.
        k: Number of bins.

    Returns:
        Threshold value = expected + sqrt(expected).
    """
    expected = n / k
    return expected + np.sqrt(expected)


def plot_single_mechanism_binplot(
    ax: plt.Axes,
    bin_counts: list[int],
    n: int,
    k: int,
    mech_label: str,
    chi2: float,
    p_chi2: float,
    cramer_v: float,
) -> None:
    """Draw a horizontal bar chart for one mechanism class at k bins.

    Format: steelblue bars, dashed expected-count line, 1-SD threshold,
    A1b baseline interval gray bands, calendar reference dotted lines.

    Args:
        ax: Matplotlib axes.
        bin_counts: Observed bin counts (length k).
        n: Total events in this class.
        k: Number of bins.
        mech_label: Display label for title (e.g., "Thrust").
        chi2: Chi-square statistic.
        p_chi2: Chi-square p-value.
        cramer_v: Cramér's V effect size.
    """
    counts = np.array(bin_counts, dtype=float)
    expected = n / k
    threshold = _compute_threshold(n, k)

    bin_edges = np.linspace(0, 1, k + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = 1.0 / k

    colors = ["orange" if c > threshold else "steelblue" for c in counts]

    ax.barh(
        bin_centers, counts,
        height=bin_width * 0.85,
        color=colors, edgecolor="white", linewidth=0.4,
    )

    # Expected count line
    ax.axvline(
        expected, color="black", linestyle="--", linewidth=1.2,
        label=f"Expected ({expected:.1f})",
    )
    # 1-SD threshold line
    ax.axvline(
        threshold, color="gray", linestyle=":", linewidth=1.0,
        label=f"1-SD threshold ({threshold:.1f})",
    )

    # A1b baseline interval shaded bands
    for ps, pe in A1B_INTERVALS:
        ax.axhspan(ps, pe, color="lightgray", alpha=0.55, zorder=0)

    # Calendar reference lines
    for pos, lbl in [
        (EQUINOX_SPRING,  "Spr Eq"),
        (SOLSTICE_SUMMER, "Sum Sol"),
        (EQUINOX_AUTUMN,  "Aut Eq"),
        (SOLSTICE_WINTER, "Win Sol"),
    ]:
        ax.axhline(pos, color="dimgray", linestyle=":", linewidth=0.7, alpha=0.7)

    ax.set_title(f"{mech_label} (n={n:,})", fontsize=10, fontweight="bold")
    ax.set_xlabel("Event Count", fontsize=9)
    ax.set_ylabel("Solar Phase (fraction of year)", fontsize=9)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0, 1.05, 0.25))
    ax.set_yticklabels(["0.0", "0.25", "0.50", "0.75", "1.0"], fontsize=8)

    # Statistical annotation
    p_str = f"{p_chi2:.2e}" if p_chi2 < 0.001 else f"{p_chi2:.4f}"
    annot = (
        f"n={n:,}\n"
        f"\u03c7\u00b2={chi2:.2f}\n"
        f"p={p_str}\n"
        f"V={cramer_v:.4f}"
    )
    ax.text(
        0.97, 0.97, annot,
        transform=ax.transAxes,
        ha="right", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85),
    )

    # Legend
    legend_handles = [
        mpatches.Patch(color="steelblue", label="Normal bin"),
        mpatches.Patch(color="orange", label=">1 SD elevated"),
        mpatches.Patch(color="lightgray", alpha=0.55, label="A1b baseline interval"),
        mlines.Line2D([0], [0], color="black", linestyle="--", linewidth=1.2,
                      label="Expected count"),
        mlines.Line2D([0], [0], color="gray", linestyle=":", linewidth=1.0,
                      label="1-SD threshold"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="lower right")


def plot_binplots(results: dict) -> None:
    """Generate Figure 1: 1x3 panel bin distributions at k=24.

    Panels: thrust (left), normal (center), strike-slip (right).

    Args:
        results: Loaded results JSON dict.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 9))
    fig.suptitle(
        f"Case B3: Solar Phase Bin Distributions by Tectonic Mechanism (k={K})",
        fontsize=13, fontweight="bold",
    )

    mech_stats = results["mechanism_stats"]
    for col_idx, mech in enumerate(MECHANISM_ORDER):
        ax = axes[col_idx]
        stats = mech_stats[mech]
        k_stats = stats[f"k{K}"]

        plot_single_mechanism_binplot(
            ax=ax,
            bin_counts=k_stats["bin_counts"],
            n=k_stats["n"],
            k=K,
            mech_label=MECH_LABELS[mech],
            chi2=k_stats["chi2"],
            p_chi2=k_stats["p_chi2"],
            cramer_v=k_stats["cramer_v"],
        )

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = OUTPUT_DIR / "case-b3-binplots.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Binplots saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 2: Mechanism comparison bar chart (dual-axis)
# ---------------------------------------------------------------------------
def plot_mechanism_comparison(results: dict) -> None:
    """Generate Figure 2: dual-axis grouped bar chart of Cramér's V and Rayleigh R.

    Left y-axis: Cramér's V (steelblue bars).
    Right y-axis: Rayleigh R (orange bars).
    Significance asterisks where p < 0.05 for chi-square.

    Args:
        results: Loaded results JSON dict.
    """
    mech_stats = results["mechanism_stats"]
    pred = results["prediction_evaluation"]

    labels = [MECH_LABELS[m] for m in MECHANISM_ORDER]
    cramer_vs = [mech_stats[m][f"k{K}"]["cramer_v"] for m in MECHANISM_ORDER]
    rayleigh_rs = [mech_stats[m][f"k{K}"]["rayleigh_R"] for m in MECHANISM_ORDER]
    p_chi2s = [mech_stats[m][f"k{K}"]["p_chi2"] for m in MECHANISM_ORDER]

    x = np.arange(len(labels))
    bar_width = 0.35

    fig, ax1 = plt.subplots(figsize=(9, 6))

    # Cramér's V bars (left axis)
    bars_v = ax1.bar(
        x - bar_width / 2, cramer_vs, bar_width,
        color="steelblue", alpha=0.85, label="Cramér's V",
    )
    ax1.set_xlabel("Tectonic Mechanism", fontsize=11)
    ax1.set_ylabel("Cramér's V", fontsize=11, color="steelblue")
    ax1.tick_params(axis="y", labelcolor="steelblue")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=10)
    ax1.set_ylim(0, max(cramer_vs) * 1.4 + 0.005)

    # Rayleigh R bars (right axis)
    ax2 = ax1.twinx()
    bars_r = ax2.bar(
        x + bar_width / 2, rayleigh_rs, bar_width,
        color="orange", alpha=0.85, label="Rayleigh R",
    )
    ax2.set_ylabel("Rayleigh R", fontsize=11, color="darkorange")
    ax2.tick_params(axis="y", labelcolor="darkorange")
    ax2.set_ylim(0, max(rayleigh_rs) * 1.4 + 0.005)

    # Significance asterisks above Cramér's V bars
    for i, (v, p) in enumerate(zip(cramer_vs, p_chi2s)):
        if p < 0.05:
            ax1.text(
                x[i] - bar_width / 2, v + 0.0005, "*",
                ha="center", va="bottom", fontsize=14, color="black", fontweight="bold",
            )

    # Métivier pattern annotation
    metivier_label = (
        "Métivier pattern (normal >= strike-slip >= thrust): MATCHED"
        if pred["matches_metivier_pattern"]
        else "Métivier pattern (normal >= strike-slip >= thrust): not matched"
    )
    fig.text(
        0.5, 0.01, metivier_label,
        ha="center", va="bottom", fontsize=9,
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.9),
    )

    ax1.set_title(
        "Case B3: Cramér's V and Rayleigh R by Tectonic Mechanism (k=24)\n"
        f"Rank order by V: {' > '.join(MECH_LABELS[m] for m in pred['cramer_v_rank_order'])}",
        fontsize=11, fontweight="bold",
    )

    # Combined legend
    handles = [
        mpatches.Patch(color="steelblue", alpha=0.85, label="Cramér's V (left axis)"),
        mpatches.Patch(color="orange", alpha=0.85, label="Rayleigh R (right axis)"),
        mlines.Line2D([0], [0], marker="*", color="black", linestyle="None",
                      markersize=10, label="p < 0.05 (chi-square)"),
    ]
    ax1.legend(handles=handles, fontsize=9, loc="upper right")

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    out_path = OUTPUT_DIR / "case-b3-mechanism-comparison.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Mechanism comparison chart saved: %s", out_path)


# ---------------------------------------------------------------------------
# Figure 3: Global coverage map
# ---------------------------------------------------------------------------
def plot_coverage_map(df: pd.DataFrame, results: dict) -> None:
    """Generate Figure 3: global event scatter map colored by mechanism class.

    Colors: thrust=red, normal=blue, strike_slip=green, unmatched=gray.
    Point size proportional to magnitude.

    Args:
        df: DataFrame with tectonic_class, latitude, longitude, usgs_mag columns.
        results: Loaded results JSON dict.
    """
    coverage = results["mechanism_stats"]["coverage"]

    class_order = ["thrust", "normal", "strike_slip", "oblique", "unmatched"]
    class_display = {
        "thrust":      f"Thrust (n={coverage['n_thrust']:,})",
        "normal":      f"Normal (n={coverage['n_normal']:,})",
        "strike_slip": f"Strike-Slip (n={coverage['n_strike_slip']:,})",
        "oblique":     f"Oblique (n={coverage['n_oblique']:,})",
        "unmatched":   f"Unmatched (n={coverage['n_unmatched']:,})",
    }
    plot_order = ["unmatched", "oblique", "strike_slip", "normal", "thrust"]

    fig, ax = plt.subplots(figsize=(15, 7))
    ax.set_facecolor("#d6eaf8")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude", fontsize=9)
    ax.set_ylabel("Latitude", fontsize=9)
    ax.grid(color="white", linewidth=0.4, alpha=0.6)
    ax.tick_params(labelsize=8)

    legend_handles = []
    for cls in plot_order:
        subset = df[df["tectonic_class"] == cls]
        if len(subset) == 0:
            continue
        color = MECH_COLORS[cls]
        alpha = 0.25 if cls == "unmatched" else 0.5
        sizes = np.clip((subset["usgs_mag"] - 5.5) ** 2 * 2, 1, 40)
        ax.scatter(
            subset["longitude"], subset["latitude"],
            s=sizes, c=color, alpha=alpha, linewidths=0,
            rasterized=True, zorder=2 if cls == "unmatched" else 3,
        )
        legend_handles.append(mpatches.Patch(
            color=color, alpha=0.8,
            label=class_display[cls],
        ))

    ax.legend(
        handles=legend_handles,
        loc="lower left", fontsize=8,
        framealpha=0.92,
        title="Mechanism Class",
        title_fontsize=8,
    )
    ax.set_title(
        "ISC-GEM Events by Focal Mechanism (GCMT Join)\n"
        "M\u22656.0, 1950\u20132021 | Point size \u221d magnitude",
        fontsize=12, fontweight="bold",
    )

    plt.tight_layout()
    out_path = OUTPUT_DIR / "case-b3-coverage-map.png"
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    logger.info("Coverage map saved: %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Load data and generate all three PNG figures."""
    logger.info("=== Case B3 Visualization ===")

    df, results = load_data()

    logger.info("Generating Figure 1: bin distribution plots...")
    plot_binplots(results)

    logger.info("Generating Figure 2: mechanism comparison chart...")
    plot_mechanism_comparison(results)

    logger.info("Generating Figure 3: global coverage map...")
    plot_coverage_map(df, results)

    logger.info("All figures written to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
