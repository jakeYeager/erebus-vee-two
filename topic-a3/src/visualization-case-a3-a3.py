"""
Case A3.A3: Phase-Concentration Audit — Visualization Script

Generates 5 PNG figures from the A3.A3 results JSON and the full catalog data.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent

RESULTS_PATH = BASE_DIR / "output" / "case-a3-a3-results.json"
CATALOG_PATH = BASE_DIR.parent / "data" / "iscgem" / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH = BASE_DIR.parent / "data" / "iscgem" / "plate-location" / "ocean_class_gshhg_global.csv"

OUTPUT_DIR = BASE_DIR / "output"

# Constants
JULIAN_YEAR_SECS = 31_557_600.0
K = 24
ELEVATED_BINS = [4, 5, 6, 7, 15, 19, 21]
SUPPRESSED_BINS = [2, 8, 10, 11, 12, 13, 16, 18, 22]

# Calendar month labels for bins (k=24, ~2 weeks per bin)
# Bin 0 starts at Jan 1
MONTH_LABELS = [
    "Jan", "", "Feb", "", "Mar", "", "Apr", "",
    "May", "", "Jun", "", "Jul", "", "Aug", "",
    "Sep", "", "Oct", "", "Nov", "", "Dec", "",
]


def load_data() -> tuple[dict, pd.DataFrame]:
    """Load results JSON and catalog DataFrame."""
    with open(RESULTS_PATH) as f:
        results = json.load(f)

    df = pd.read_csv(CATALOG_PATH)
    gshhg = pd.read_csv(GSHHG_PATH)
    df = df.merge(gshhg[["usgs_id", "ocean_class", "dist_to_coast_km"]], on="usgs_id", how="left")
    df["phase"] = (df["solar_secs"] / JULIAN_YEAR_SECS) % 1.0
    df["bin_index"] = (np.floor(df["phase"] * K).astype(int)) % K

    # Tectonic class
    def classify_tectonic(dist: float) -> str:
        if dist <= 50.0:
            return "continental"
        elif dist <= 200.0:
            return "transitional"
        else:
            return "oceanic"

    df["tectonic_class"] = df["dist_to_coast_km"].apply(classify_tectonic)

    # Depth band
    def classify_depth(d: float) -> str:
        if d < 20.0:
            return "shallow"
        elif d < 70.0:
            return "mid_crustal"
        elif d < 300.0:
            return "intermediate"
        else:
            return "deep"

    df["depth_band"] = df["depth"].apply(classify_depth)

    # Sequence role from bin_summary (use results JSON for influence)
    bin_summary = results["bin_summary"]
    bin_influence = {b["bin"]: b["bin_influence"] for b in bin_summary}
    df["chi2_influence"] = df["bin_index"].map(bin_influence)

    # Sequence roles from results
    # Re-build from aftershock files for the sequence_role column
    from pathlib import Path as _Path

    gk_after = pd.read_csv(BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_gk-seq_global.csv")
    reas_after = pd.read_csv(BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_reas-seq_global.csv")
    a1b_after = pd.read_csv(BASE_DIR.parent / "data" / "iscgem" / "declustering-algorithm" / "aftershocks_a1b-seq_global.csv")

    mainshock_with_aftershocks: set[str] = set()
    aftershock_member: set[str] = set()
    for adf in [gk_after, reas_after, a1b_after]:
        mainshock_with_aftershocks.update(adf["parent_id"].dropna().astype(str))
        aftershock_member.update(adf["usgs_id"].dropna().astype(str))

    def classify_role(uid: str) -> str:
        if uid in aftershock_member:
            return "aftershock"
        elif uid in mainshock_with_aftershocks:
            return "mainshock_with_sequence"
        else:
            return "isolated_mainshock"

    df["sequence_role"] = df["usgs_id"].astype(str).apply(classify_role)

    return results, df


# ---------------------------------------------------------------------------
# Figure 1: Signed influence distribution by sequence role
# ---------------------------------------------------------------------------

def figure1_influence_distribution(results: dict, df: pd.DataFrame) -> None:
    """Figure 1: Histogram of chi2_influence per sequence role."""
    logger.info("Generating Figure 1: Influence Distribution...")

    roles = ["isolated_mainshock", "mainshock_with_sequence", "aftershock"]
    role_labels = ["Isolated Mainshock", "Mainshock with Sequence", "Aftershock"]

    # Determine shared x-axis range
    all_influences = df["chi2_influence"].values
    x_min = all_influences.min() - 0.005
    x_max = all_influences.max() + 0.005

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=False)

    for ax, role, label in zip(axes, roles, role_labels):
        subset = df[df["sequence_role"] == role]["chi2_influence"].values
        n_events = len(subset)

        # Compute fractions
        frac_elevated = (subset > 0).mean()
        frac_suppressed = (subset < 0).mean()

        # Get unique influence values (each bin has one value)
        unique_infs = np.unique(subset)

        # Build bar chart using actual bin influence values
        for infl in unique_infs:
            count = (subset == infl).sum()
            if infl > 0:
                color = "coral"
            elif infl < 0:
                color = "steelblue"
            else:
                color = "gray"
            ax.bar(infl, count, width=(x_max - x_min) / (K * 1.5), color=color, edgecolor="white", linewidth=0.5)

        ax.axvline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
        ax.set_xlim(x_min, x_max)
        ax.set_xlabel("Chi-Square Influence", fontsize=10)
        ax.set_ylabel("Event Count", fontsize=10)
        ax.set_title(label, fontsize=11, fontweight="bold")

        # Annotation
        annotation = (
            f"n = {n_events:,}\n"
            f"Elevated: {frac_elevated:.1%}\n"
            f"Suppressed: {frac_suppressed:.1%}"
        )
        ax.text(
            0.97, 0.97, annotation,
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.8),
        )

        # Legend
        handles = [
            mpatches.Patch(color="coral", label="Elevated bin (infl > 0)"),
            mpatches.Patch(color="steelblue", label="Suppressed bin (infl < 0)"),
            mpatches.Patch(color="gray", label="Neutral"),
        ]
        ax.legend(handles=handles, fontsize=8, loc="upper left")

    fig.suptitle("Signed Chi-Square Influence by Sequence Role (A3.A3)", fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()

    out_path = OUTPUT_DIR / "case-a3-a3-influence-distribution.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 1 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Parallel degradation curves
# ---------------------------------------------------------------------------

def figure2_degradation_curves(results: dict) -> None:
    """Figure 2: Two parallel removal curves vs baselines."""
    logger.info("Generating Figure 2: Degradation Curves...")

    deg = results["degradation_curves"]

    # Critical chi2 at p=0.05, dof=K-1=23
    from scipy.stats import chi2 as chi2_dist
    critical_chi2 = chi2_dist.ppf(0.95, df=K - 1)

    # A1b interval bins to track
    A1B_I1 = [4, 5]
    A1B_I2 = [15]
    A1B_I3 = [21]
    STRONG_SUPP = [13]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=False)

    for col, (curve_key, gk_key, rand_key, title_suffix, removal_label) in enumerate([
        ("elevated_removal", "gk_sequential_baseline_elevated", "random_baseline_elevated",
         "Elevated-Bin Removal", "Elevated"),
        ("suppressed_removal", "gk_sequential_baseline_suppressed", "random_baseline_suppressed",
         "Suppressed-Bin Removal", "Suppressed"),
    ]):
        curve_data = deg[curve_key]["curve"]
        gk_data = deg[gk_key]["curve"]
        rand_mean = deg[rand_key]["curve_mean_chi2"]
        rand_p10 = deg[rand_key]["curve_p10_chi2"]
        rand_p90 = deg[rand_key]["curve_p90_chi2"]

        pct_x = [d["pct_catalog_removed"] for d in curve_data]
        chi2_y = [d["chi2"] for d in curve_data]

        gk_pct_x = [d["pct_catalog_removed"] for d in gk_data]
        gk_chi2_y = [d["chi2"] for d in gk_data]

        rand_pct_x = [i / (len(rand_mean) - 1) * pct_x[-1] for i in range(len(rand_mean))]

        persistence_step = deg[curve_key]["signal_persistence_step"]
        persistence_pct = deg[curve_key]["signal_persistence_pct_catalog"]

        # Top panel: chi2 curve
        ax_top = axes[0, col]
        ax_top.fill_between(rand_pct_x, rand_p10, rand_p90, alpha=0.3, color="lightgray", label="Random baseline P10–P90")
        ax_top.plot(rand_pct_x, rand_mean, color="gray", linestyle="-", linewidth=1.0, alpha=0.7, label="Random baseline mean")
        ax_top.plot(gk_pct_x, gk_chi2_y, color="darkorange", linestyle="--", linewidth=1.5, label="G-K sequential baseline")
        ax_top.plot(pct_x, chi2_y, color="steelblue" if col == 0 else "coral", linewidth=2.0, label=f"{removal_label} removal")
        ax_top.axhline(critical_chi2, color="red", linestyle="--", linewidth=1.2, label=f"p=0.05 threshold (χ²={critical_chi2:.1f})")

        # Mark signal persistence
        if persistence_pct is not None and persistence_pct < pct_x[-1]:
            ax_top.axvline(persistence_pct, color="navy", linestyle=":", linewidth=1.5, alpha=0.8)
            ax_top.annotate(
                f"p≥0.05 at\n{persistence_pct:.1f}%",
                xy=(persistence_pct, critical_chi2),
                xytext=(persistence_pct + 1, critical_chi2 * 1.3),
                fontsize=8, color="navy",
                arrowprops=dict(arrowstyle="->", color="navy", lw=1.0),
            )
        else:
            ax_top.annotate(
                "Signal persists\nthrough full removal",
                xy=(pct_x[-1] * 0.5, critical_chi2 * 1.5),
                fontsize=9, color="darkred", ha="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", edgecolor="red", alpha=0.8),
            )

        ax_top.set_title(f"{title_suffix} — Chi-Square Track", fontsize=10, fontweight="bold")
        ax_top.set_ylabel("Chi-Square Statistic", fontsize=9)
        ax_top.set_xlabel("% Catalog Removed", fontsize=9)
        ax_top.legend(fontsize=8, loc="upper right")
        ax_top.grid(True, alpha=0.3)

        # Bottom panel: z-score tracking
        ax_bot = axes[1, col]

        # Aggregate A1b intervals z-scores (average of bins in interval)
        def get_avg_z(curve: list[dict], key: str, bins: list[int]) -> list[float]:
            return [np.mean([d[key].get(str(b), d[key].get(b, 0.0)) for b in bins]) for d in curve]

        pct_x_full = [d["pct_catalog_removed"] for d in curve_data]

        # Z-score for elevated bins
        z_i1 = get_avg_z(curve_data, "elevated_bin_zscores", A1B_I1)
        z_i2 = get_avg_z(curve_data, "elevated_bin_zscores", A1B_I2)
        z_i3 = get_avg_z(curve_data, "elevated_bin_zscores", A1B_I3)
        z_s13 = get_avg_z(curve_data, "suppressed_bin_zscores", STRONG_SUPP)

        ax_bot.plot(pct_x_full, z_i1, color="tomato", linewidth=1.5, label="A1b Interval 1 (bins 4–5)")
        ax_bot.plot(pct_x_full, z_i2, color="orangered", linewidth=1.5, label="A1b Interval 2 (bin 15)")
        ax_bot.plot(pct_x_full, z_i3, color="firebrick", linewidth=1.5, label="A1b Interval 3 (bin 21)")
        ax_bot.plot(pct_x_full, z_s13, color="steelblue", linewidth=1.5, linestyle="--", label="Strongest suppressed (bin 13)")

        ax_bot.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
        ax_bot.axhline(-1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
        ax_bot.axhline(0.0, color="black", linestyle="-", linewidth=0.5, alpha=0.3)

        ax_bot.set_title(f"{title_suffix} — Bin Z-Score Tracking", fontsize=10, fontweight="bold")
        ax_bot.set_ylabel("Z-Score", fontsize=9)
        ax_bot.set_xlabel("% Catalog Removed", fontsize=9)
        ax_bot.legend(fontsize=8, loc="upper right")
        ax_bot.grid(True, alpha=0.3)

    fig.suptitle("Phase-Concentration Audit — Sequential Removal (A3.A3)", fontsize=13, fontweight="bold")
    plt.tight_layout()

    out_path = OUTPUT_DIR / "case-a3-a3-degradation-curves.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 2 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Permutation confidence bands
# ---------------------------------------------------------------------------

def figure3_permutation_tails(results: dict) -> None:
    """Figure 3: Per-bin counts with permutation confidence bands."""
    logger.info("Generating Figure 3: Permutation Confidence Bands...")

    perm = results["permutation_baseline"]
    bin_summary = results["bin_summary"]

    obs_counts = [b["obs"] for b in bin_summary]
    expected = bin_summary[0]["expected"]
    bin_p5 = perm["bin_p5"]
    bin_p95 = perm["bin_p95"]
    sig_elev = perm["sig_elevated_bins_permutation"]
    sig_supp = perm["sig_suppressed_bins_permutation"]

    x = np.arange(K)

    fig, ax = plt.subplots(figsize=(13, 5))

    # Permutation band
    ax.fill_between(x, bin_p5, bin_p95, alpha=0.3, color="lightgray", label="Permutation 5th–95th percentile")

    # Bars
    colors = []
    for b in range(K):
        if b in ELEVATED_BINS:
            colors.append("coral")
        elif b in SUPPRESSED_BINS:
            colors.append("steelblue")
        else:
            colors.append("lightgray")

    bars = ax.bar(x, obs_counts, color=colors, edgecolor="white", linewidth=0.5, alpha=0.85, label="Observed count")

    # Expected line
    ax.axhline(expected, color="black", linestyle="--", linewidth=1.2, alpha=0.7, label=f"Expected ({expected:.0f})")

    # Asterisks for significant bins
    for b in sig_elev:
        ax.text(b, obs_counts[b] + 5, "*", ha="center", va="bottom", color="darkred", fontsize=14, fontweight="bold")

    for b in sig_supp:
        ax.text(b, obs_counts[b] - 5, "*", ha="center", va="top", color="darkblue", fontsize=14, fontweight="bold")

    ax.set_xlabel("Bin Index (0 = Jan 1)", fontsize=10)
    ax.set_ylabel("Event Count", fontsize=10)
    ax.set_title(
        f"Per-Bin Counts vs. Permutation Confidence Bands (A3.A3, n={results['permutation_baseline']['n_permutations']:,} permutations)",
        fontsize=11, fontweight="bold",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([str(i) for i in range(K)], fontsize=8)

    # Secondary x-axis with month labels
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(x)
    ax2.set_xticklabels(MONTH_LABELS, fontsize=8)
    ax2.set_xlabel("Approximate Calendar Month", fontsize=9)

    # Legend
    legend_handles = [
        mpatches.Patch(color="coral", alpha=0.85, label="Elevated bin"),
        mpatches.Patch(color="steelblue", alpha=0.85, label="Suppressed bin"),
        mpatches.Patch(color="lightgray", alpha=0.85, label="Neutral bin"),
        mpatches.Patch(color="lightgray", alpha=0.3, label="Permutation 5th–95th pct."),
        Line2D([0], [0], color="black", linestyle="--", label=f"Expected ({expected:.0f})"),
        Line2D([0], [0], color="darkred", marker="*", linestyle="None", markersize=10, label="Sig. elevated (p<0.05)"),
        Line2D([0], [0], color="darkblue", marker="*", linestyle="None", markersize=10, label="Sig. suppressed (p<0.05)"),
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc="upper left")

    plt.tight_layout()

    out_path = OUTPUT_DIR / "case-a3-a3-permutation-tails.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 3 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 4: Representativeness heatmap
# ---------------------------------------------------------------------------

def figure4_representativeness(results: dict, df: pd.DataFrame) -> None:
    """Figure 4: Tectonic × depth profile of top-influence groups."""
    logger.info("Generating Figure 4: Representativeness Heatmap...")

    repr_data = results["representativeness"]

    tectonic_cats = ["continental", "transitional", "oceanic"]
    depth_cats = ["shallow", "mid_crustal", "intermediate", "deep"]
    n_full = len(df)

    # Signal-bearing stratum
    signal_stratum = df[(df["tectonic_class"] == "continental") | (df["depth_band"] == "mid_crustal")]
    n_stratum = len(signal_stratum)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    groups = [
        ("top_positive_100", "Top-100 Positive Influence"),
        ("top_negative_100", "Top-100 Negative Influence"),
    ]

    dimension_labels = [
        ("tectonic_dist", tectonic_cats, "Tectonic Class"),
        ("depth_dist", depth_cats, "Depth Band"),
    ]

    for col, (group_key, group_label) in enumerate(groups):
        gdata = repr_data[group_key]
        n_group = gdata["n"]

        for row, (dist_key, cats, dim_label) in enumerate(dimension_labels):
            ax = axes[row, col]

            # Group values
            group_vals = np.array([gdata[dist_key].get(c, 0) for c in cats], dtype=float)

            # Full catalog proportions scaled to group size
            full_vals = np.array([(df[dim_label.lower().replace(" ", "_")] if dim_label == "Depth Band"
                                    else df["tectonic_class"]) == c for c in cats], dtype=float)
            # Re-compute directly
            if dist_key == "tectonic_dist":
                full_props = np.array([(df["tectonic_class"] == c).sum() / n_full for c in cats])
                stratum_props = np.array([(signal_stratum["tectonic_class"] == c).sum() / n_stratum for c in cats])
            else:
                full_props = np.array([(df["depth_band"] == c).sum() / n_full for c in cats])
                stratum_props = np.array([(signal_stratum["depth_band"] == c).sum() / n_stratum for c in cats])

            x = np.arange(len(cats))
            width = 0.3

            # Group bars
            ax.bar(x - width / 2, group_vals / n_group, width=width, color="steelblue", alpha=0.85,
                   label=f"Group (n={n_group})")

            # Full catalog bars
            ax.bar(x + width / 2, full_props, width=width, color="lightgray", alpha=0.85,
                   label="Full catalog")

            # Signal stratum outline
            for xi, sp in enumerate(stratum_props):
                ax.bar(xi + width / 2 + width * 0.5, sp, width=width * 0.15,
                       color="none", edgecolor="darkorange", linewidth=2, linestyle="--",
                       label="Signal stratum" if xi == 0 else "")

            ax.set_xticks(x)
            ax.set_xticklabels(cats, rotation=20, ha="right", fontsize=9)
            ax.set_ylabel("Fraction", fontsize=9)
            ax.set_title(f"{group_label}\n{dim_label}", fontsize=9, fontweight="bold")

            # p-value annotation
            if dist_key == "tectonic_dist":
                p_full = gdata.get("p_vs_full_tectonic", gdata.get("p_vs_full", 1.0))
                p_stratum = gdata.get("p_vs_signal_stratum_tectonic", gdata.get("p_vs_signal_stratum", 1.0))
            else:
                p_full = gdata.get("p_vs_full_depth", gdata.get("p_vs_full", 1.0))
                p_stratum = gdata.get("p_vs_signal_stratum_depth", gdata.get("p_vs_signal_stratum", 1.0))

            def sig_label(p: float) -> str:
                if p <= 0.01:
                    return "**"
                elif p <= 0.05:
                    return "*"
                else:
                    return "ns"

            ann = (
                f"vs full: p={p_full:.3f} {sig_label(p_full)}\n"
                f"vs stratum: p={p_stratum:.3f} {sig_label(p_stratum)}"
            )
            ax.text(0.97, 0.97, ann, transform=ax.transAxes, ha="right", va="top",
                    fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                          edgecolor="gray", alpha=0.8))

            ax.legend(fontsize=8, loc="upper left")
            ax.grid(True, alpha=0.2)
            ax.set_ylim(0, max(group_vals.max() / n_group, full_props.max(), stratum_props.max()) * 1.4)

    fig.suptitle("Tectonic and Depth Profile of Top-Influence Event Groups (A3.A3)", fontsize=12, fontweight="bold")
    plt.tight_layout()

    out_path = OUTPUT_DIR / "case-a3-a3-representativeness.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 4 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 5: Sequence membership by bin type
# ---------------------------------------------------------------------------

def figure5_sequence_annotation(results: dict, df: pd.DataFrame) -> None:
    """Figure 5: Stacked bar chart of sequence role by bin type."""
    logger.info("Generating Figure 5: Sequence Annotation...")

    bin_summary = results["bin_summary"]
    bs_map = {b["bin"]: b for b in bin_summary}

    # Define bin type groups
    A1B_I1_BINS = [4, 5]
    A1B_I2_BINS = [15]
    A1B_I3_BINS = [21]
    OTHER_ELEVATED = [b for b in ELEVATED_BINS if b not in A1B_I1_BINS + A1B_I2_BINS + A1B_I3_BINS]
    SUPP_STRONG = [b for b in SUPPRESSED_BINS if bs_map[b]["z"] < -2.0]
    SUPP_MOD = [b for b in SUPPRESSED_BINS if -2.0 <= bs_map[b]["z"] < -1.0]
    NEUTRAL_BINS = [b for b in range(K) if b not in ELEVATED_BINS and b not in SUPPRESSED_BINS]

    group_defs = [
        ("Elevated A1b Interval 1\n(bins 4–5)", A1B_I1_BINS),
        ("Elevated A1b Interval 2\n(bin 15)", A1B_I2_BINS),
        ("Elevated A1b Interval 3\n(bin 21)", A1B_I3_BINS),
        ("Other Elevated", OTHER_ELEVATED),
        ("Suppressed Strong (z<-2)", SUPP_STRONG),
        ("Suppressed Moderate (-2≤z<-1)", SUPP_MOD),
        ("Neutral", NEUTRAL_BINS),
    ]

    role_cols = ["isolated_mainshock", "mainshock_with_sequence", "aftershock"]
    role_colors = ["steelblue", "darkorange", "gray"]
    role_labels = ["Isolated Mainshock", "Mainshock with Sequence", "Aftershock"]

    # Full catalog fractions
    n_full = len(df)
    full_fracs = {r: (df["sequence_role"] == r).sum() / n_full for r in role_cols}

    fig, ax = plt.subplots(figsize=(11, 7))

    y_positions = list(range(len(group_defs)))
    y_labels = [g[0] for g in group_defs]
    bar_height = 0.6

    for yi, (group_name, group_bins) in enumerate(group_defs):
        # Get events in these bins
        mask = df["bin_index"].isin(group_bins)
        group_df = df[mask]
        n_group = len(group_df)
        if n_group == 0:
            continue

        fracs = [
            (group_df["sequence_role"] == r).sum() / n_group for r in role_cols
        ]

        left = 0.0
        for frac, color in zip(fracs, role_colors):
            ax.barh(yi, frac, left=left, height=bar_height, color=color, alpha=0.85, edgecolor="white", linewidth=0.5)
            left += frac

        # Count annotation
        ax.text(1.01, yi, f"n={n_group:,}", va="center", ha="left", fontsize=9)

    # Reference lines for full catalog proportions
    left_ref = 0.0
    for r, color in zip(role_cols, role_colors):
        frac = full_fracs[r]
        ax.axvline(left_ref + frac, color=color, linestyle=":", linewidth=1.5, alpha=0.8)
        left_ref += frac

    ax.set_xlim(0, 1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlabel("Fraction of Events", fontsize=10)
    ax.set_title("Sequence Role by Bin Type — Phase-Concentration Audit (A3.A3)", fontsize=11, fontweight="bold")

    # Legend
    role_handles = [
        mpatches.Patch(color=color, alpha=0.85, label=label)
        for color, label in zip(role_colors, role_labels)
    ]
    ref_handle = Line2D([0], [0], color="black", linestyle=":", linewidth=1.5, label="Full catalog proportion")
    role_handles.append(ref_handle)
    ax.legend(handles=role_handles, fontsize=9, loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)

    plt.tight_layout()

    out_path = OUTPUT_DIR / "case-a3-a3-sequence-annotation.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 5 saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    logger.info("=== Case A3.A3 Visualization ===")
    results, df = load_data()

    figure1_influence_distribution(results, df)
    figure2_degradation_curves(results)
    figure3_permutation_tails(results)
    figure4_representativeness(results, df)
    figure5_sequence_annotation(results, df)

    logger.info("=== All figures generated ===")


if __name__ == "__main__":
    main()
