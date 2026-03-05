"""
Case A3.A2 — Visualization Script

Produces 5 publication-quality figures from the A3.A2 results JSON:
  1. case-a3-a2-baseline-spectrum.png — Baseline Schuster spectrum
  2. case-a3-a2-stratum-mfpa.png — MFPA stratification
  3. case-a3-a2-nh-sh-phase.png — NH/SH phase angle comparison
  4. case-a3-a2-declustering-sensitivity.png — Declustering sensitivity
  5. case-a3-a2-cluster-window.png — Cluster-window sensitivity
"""

import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
)
log = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a3-a2-results.json"
OUTPUT_DIR = BASE_DIR / "output"

# Key period marker definitions
KEY_PERIODS = {
    75.6: ("Period A2.A1 (~75.6d)", "darkorange", ":"),
    91.25: ("Quarter-year (91.25d)", "orange", ":"),
    182.5: ("Half-year (182.5d)", "purple", ":"),
    365.25: ("Annual (365.25d)", "red", ":"),
}
TIDAL_PERIODS = [0.5, 1.0, 14.77, 27.32]

STEELBLUE = "steelblue"
EXPLICIT_PERIOD_LABELS = {
    "quarter_year_756": "75.6d",
    "quarter_year_912": "91.25d",
    "half_year": "182.5d",
    "annual": "365.25d",
}


def load_results() -> dict:
    with open(RESULTS_PATH) as f:
        return json.load(f)


def add_period_markers(ax: plt.Axes, ymin: float, ymax: float) -> None:
    """Add vertical lines for key periods."""
    for period, (label, color, ls) in KEY_PERIODS.items():
        ax.axvline(period, color=color, linestyle=ls, linewidth=1.0, alpha=0.8, label=label)
    for tp in TIDAL_PERIODS:
        ax.axvline(tp, color="lightgray", linestyle=":", linewidth=0.8, alpha=0.6)


def format_schuster_axes(ax: plt.Axes, title: str, annotation: str) -> None:
    """Common axis formatting for Schuster spectrum panels."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Period (days)", fontsize=10)
    ax.set_ylabel("p-value", fontsize=10)
    ax.set_xlim(0.25, 548.0)
    ax.axhline(0.05, color="red", linestyle="--", linewidth=1.0, alpha=0.8, label="p=0.05")
    ax.axhline(0.001, color="red", linestyle=":", linewidth=0.8, alpha=0.6, label="p=0.001")
    ax.set_title(title, fontsize=10)
    ax.text(0.98, 0.04, annotation, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color="gray",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7))


# ---------------------------------------------------------------------------
# Figure 1 — Baseline Schuster spectrum
# ---------------------------------------------------------------------------

def figure1(results: dict) -> None:
    log.info("Figure 1: Baseline Schuster spectrum")
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    fig.suptitle(
        "Cluster-Robust Schuster Spectrum — Full Catalog and Signal-Bearing Stratum (A3.A2)",
        fontsize=12, fontweight="bold", y=0.98,
    )

    datasets = [
        (
            "subtest_1_baseline",
            "schuster_1day",
            "Full Catalog (n=9,210)",
            results["strata_sizes"]["n_full"],
        ),
        (
            "subtest_2_stratified",
            None,  # special — nested
            f"Signal-Bearing Stratum (continental × mid-crustal, n={results['strata_sizes']['n_signal_bearing']})",
            results["strata_sizes"]["n_signal_bearing"],
        ),
    ]

    for ax_idx, (key, subkey, panel_title, n_events) in enumerate(datasets):
        ax = axes[ax_idx]

        if subkey is None:
            s1 = results["subtest_2_stratified"]["signal_bearing"]["schuster_1day"]
        else:
            s1 = results[key][subkey]

        periods = np.array([r["period_days"] for r in s1])
        p_std = np.array([r["p_standard"] for r in s1])
        p_cr = np.array([r["p_cluster_robust"] for r in s1])
        n_clusters_arr = np.array([r["n_clusters"] for r in s1])
        n_clusters_mean = int(np.mean(n_clusters_arr))

        # Compute inflation factor for this panel
        n_std_sig = int(np.sum(p_std < 0.05))
        n_cr_sig = int(np.sum(p_cr < 0.05))
        inf_factor = n_std_sig / n_cr_sig if n_cr_sig > 0 else float("inf")

        add_period_markers(ax, 0, 1)
        ax.plot(periods, p_std, color="gray", linewidth=0.8, alpha=0.7, label="Standard p-value")
        ax.plot(periods, p_cr, color=STEELBLUE, linewidth=1.5, label="Cluster-robust p-value (1d)")

        format_schuster_axes(
            ax, panel_title,
            f"n={n_events:,}  n_clusters≈{n_clusters_mean:,}  inflation={inf_factor:.1f}×",
        )
        ax.legend(fontsize=7, loc="upper right", ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = OUTPUT_DIR / "case-a3-a2-baseline-spectrum.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Figure 2 — MFPA stratification
# ---------------------------------------------------------------------------

def figure2(results: dict) -> None:
    log.info("Figure 2: MFPA stratification")
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    fig.suptitle(
        "MFPA Periodogram — Signal-Bearing vs. Non-Signal-Bearing Stratum (A3.A2)",
        fontsize=12, fontweight="bold", y=0.98,
    )

    strata = [
        ("signal_bearing", "Signal-Bearing Stratum (continental × mid-crustal)"),
        ("non_signal", "Non-Signal-Bearing Remainder"),
    ]

    for ax_idx, (stratum_key, panel_title) in enumerate(strata):
        ax = axes[ax_idx]
        mfpa = results["subtest_2_stratified"][stratum_key]["mfpa"]
        n = results["subtest_2_stratified"][stratum_key]["n"]

        periods = np.array([r["period_days"] for r in mfpa])
        power = np.array([r["power"] for r in mfpa])
        p95 = mfpa[0]["p95_threshold"]
        p99 = mfpa[0]["p99_threshold"]
        fdr_flags = np.array([r["fdr_significant"] for r in mfpa])

        ax.set_xscale("log")
        ax.plot(periods, power, color=STEELBLUE, linewidth=1.2, label="Observed power")
        ax.axhline(p95, color="gray", linestyle="--", linewidth=1.0, label="p95 threshold")
        ax.axhline(p99, color="gray", linestyle=":", linewidth=0.8, label="p99 threshold")

        # Fill where power > p95
        above = power > p95
        ax.fill_between(periods, power, p95, where=above,
                        color="steelblue", alpha=0.2, label="Power > p95")

        # Mark FDR-significant periods
        fdr_periods = periods[fdr_flags]
        fdr_powers = power[fdr_flags]
        if len(fdr_periods) > 0:
            ax.scatter(fdr_periods, fdr_powers + (fdr_powers * 0.08),
                       marker="^", color="orange", s=60, zorder=5, label="FDR significant")
            for fp, fw in zip(fdr_periods, fdr_powers):
                ax.annotate(f"{fp:.1f}d", (fp, fw * 1.12),
                            fontsize=7, ha="center", color="darkorange")

        add_period_markers(ax, 0, 1)

        ax.set_xlabel("Period (days)", fontsize=10)
        ax.set_ylabel("MFPA Power", fontsize=10)
        ax.set_xlim(0.25, 548.0)
        ax.set_title(f"{panel_title} (n={n:,})", fontsize=10)
        ax.legend(fontsize=7, loc="upper right", ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = OUTPUT_DIR / "case-a3-a2-stratum-mfpa.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Figure 3 — NH/SH phase angle comparison
# ---------------------------------------------------------------------------

def figure3(results: dict) -> None:
    log.info("Figure 3: NH/SH phase angle comparison")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("NH/SH Phase Angle and Power Comparison (A3.A2)", fontsize=12, fontweight="bold", y=0.98)

    nh = results["subtest_3_nh_sh"]["nh"]
    sh = results["subtest_3_nh_sh"]["sh"]
    delta_days = results["subtest_3_nh_sh"]["half_year_phase_delta_days"]

    # ---- Top row: phase angle polar-style panels ----
    def phase_panel(ax_polar_host, nh_phase: float, sh_phase: float, period_label: str, period_days: float) -> None:
        """Draw a polar-style phase arrow plot on a standard axes."""
        ax_polar_host.set_aspect("equal")
        ax_polar_host.set_xlim(-1.5, 1.5)
        ax_polar_host.set_ylim(-1.5, 1.5)

        # Draw circle
        theta_circle = np.linspace(0, 2 * np.pi, 300)
        ax_polar_host.plot(np.cos(theta_circle), np.sin(theta_circle), "k-", linewidth=0.8, alpha=0.3)

        # NH arrow (warm red)
        nh_ang = nh_phase * 2 * np.pi
        ax_polar_host.annotate(
            "", xy=(np.cos(nh_ang), np.sin(nh_ang)),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color="crimson", lw=2.0),
        )
        ax_polar_host.text(
            1.15 * np.cos(nh_ang), 1.15 * np.sin(nh_ang),
            f"NH\n{nh_phase:.3f}", ha="center", va="center", fontsize=8, color="crimson"
        )

        # SH arrow (cool blue)
        sh_ang = sh_phase * 2 * np.pi
        ax_polar_host.annotate(
            "", xy=(np.cos(sh_ang), np.sin(sh_ang)),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color="steelblue", lw=2.0),
        )
        ax_polar_host.text(
            1.15 * np.cos(sh_ang), 1.15 * np.sin(sh_ang),
            f"SH\n{sh_phase:.3f}", ha="center", va="center", fontsize=8, color="steelblue"
        )

        # Arc annotations
        ax_polar_host.set_xticks([])
        ax_polar_host.set_yticks([])
        ax_polar_host.spines["top"].set_visible(False)
        ax_polar_host.spines["right"].set_visible(False)
        ax_polar_host.spines["left"].set_visible(False)
        ax_polar_host.spines["bottom"].set_visible(False)

        ax_polar_host.set_title(
            f"{period_label}\n(period={period_days:.2f}d)\nΔ={delta_days:.1f}d | anti-phase=91.25d | B2 pred≈36d",
            fontsize=8
        )

    nh_hy_phase = nh["half_year_dominant_phase"]
    sh_hy_phase = sh["half_year_dominant_phase"]
    nh_qy_phase = nh["quarter_year_912_dominant_phase"]
    sh_qy_phase = sh["quarter_year_912_dominant_phase"]

    phase_panel(axes[0, 0], nh_hy_phase, sh_hy_phase, "Half-Year (182.5d)", 182.5)
    phase_panel(axes[0, 1], nh_qy_phase, sh_qy_phase, "Quarter-Year (91.25d)", 91.25)

    # ---- Bottom row: MFPA power grouped bar charts ----
    ep_keys = ["quarter_year_756", "quarter_year_912", "half_year", "annual"]
    ep_labels = ["75.6d", "91.25d", "182.5d", "365.25d"]

    nh_mfpa = results["subtest_3_nh_sh"]["nh"]["explicit_period_tests"]
    sh_mfpa = results["subtest_3_nh_sh"]["sh"]["explicit_period_tests"]

    nh_powers = [nh_mfpa[k]["power"] for k in ep_keys]
    sh_powers = [sh_mfpa[k]["power"] for k in ep_keys]
    nh_sig = [nh_mfpa[k]["p_mfpa"] for k in ep_keys]
    sh_sig = [sh_mfpa[k]["p_mfpa"] for k in ep_keys]

    # Use p95 threshold from any MFPA result
    p95_val = results["subtest_3_nh_sh"]["nh"]["explicit_period_tests"]["half_year"]["significant_95"]

    # Get actual p95 threshold from mfpa list
    nh_mfpa_list = results["subtest_3_nh_sh"]["nh"]["explicit_period_tests"]["half_year"]
    # We need raw p95 — stored in explicit test:
    # actually stored only in the scan array entries; use a fallback from explicit
    # The explicit period tests store significant_95 bool but not p95 value.
    # Retrieve from the full mfpa scan in subtest_1 (same n_bootstrap, same rng)
    # p95 should be similar; use from subtest_1 baseline
    p95_ref = results["subtest_1_baseline"]["mfpa"][0]["p95_threshold"]

    x = np.arange(len(ep_keys))
    width = 0.35

    for col_idx, (powers_nh, powers_sh, sigs_nh, sigs_sh, legend_title) in enumerate([
        (nh_powers, sh_powers, nh_sig, sh_sig, "NH vs SH MFPA Power"),
    ]):
        ax = axes[1, 0]
        bars_nh = ax.bar(x - width / 2, powers_nh, width, color="crimson", alpha=0.7, label="NH")
        bars_sh = ax.bar(x + width / 2, powers_sh, width, color=STEELBLUE, alpha=0.7, label="SH")
        ax.axhline(p95_ref, color="gray", linestyle="--", linewidth=1.0, label="p95 threshold")

        # Significance annotations
        for i, (pnh, psh, pv_nh, pv_sh) in enumerate(zip(powers_nh, powers_sh, sigs_nh, sigs_sh)):
            sig_nh = "**" if pv_nh < 0.01 else ("*" if pv_nh < 0.05 else "ns")
            sig_sh = "**" if pv_sh < 0.01 else ("*" if pv_sh < 0.05 else "ns")
            ax.text(x[i] - width / 2, pnh + 0.002, sig_nh, ha="center", fontsize=7, color="crimson")
            ax.text(x[i] + width / 2, psh + 0.002, sig_sh, ha="center", fontsize=7, color=STEELBLUE)

        ax.set_xticks(x)
        ax.set_xticklabels(ep_labels)
        ax.set_xlabel("Period", fontsize=9)
        ax.set_ylabel("MFPA Power", fontsize=9)
        ax.set_title("MFPA Power: NH vs SH at Explicit Periods", fontsize=9)
        ax.legend(fontsize=8)

    # Second bottom panel: replicate with just significance flags
    ax2 = axes[1, 1]
    nh_p_mfpa_list = [nh_mfpa[k]["p_mfpa"] for k in ep_keys]
    sh_p_mfpa_list = [sh_mfpa[k]["p_mfpa"] for k in ep_keys]

    ax2.bar(x - width / 2, nh_powers, width, color="crimson", alpha=0.7, label="NH")
    ax2.bar(x + width / 2, sh_powers, width, color=STEELBLUE, alpha=0.7, label="SH")
    ax2.axhline(p95_ref, color="gray", linestyle="--", linewidth=1.0, label="p95 threshold")

    # Phase offset info box
    textstr = (
        f"Half-year Δ = {delta_days:.1f}d\n"
        f"Anti-phase ref: 91.25d\n"
        f"A3.B2 pred: ~36d\n"
        f"Closer to antiphase: {results['subtest_3_nh_sh']['closer_to_antiphase']}"
    )
    ax2.text(0.97, 0.97, textstr, transform=ax2.transAxes,
             ha="right", va="top", fontsize=7,
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    ax2.set_xticks(x)
    ax2.set_xticklabels(ep_labels)
    ax2.set_xlabel("Period", fontsize=9)
    ax2.set_ylabel("MFPA Power", fontsize=9)
    ax2.set_title("MFPA Power: NH vs SH (phase offset summary)", fontsize=9)
    ax2.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = OUTPUT_DIR / "case-a3-a2-nh-sh-phase.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Figure 4 — Declustering sensitivity
# ---------------------------------------------------------------------------

def figure4(results: dict) -> None:
    log.info("Figure 4: Declustering sensitivity")
    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
    fig.suptitle(
        "Declustering Sensitivity — Schuster Spectrum Across Catalog Versions (A3.A2)",
        fontsize=12, fontweight="bold", y=0.98,
    )

    catalogs = [
        ("full", "Full Catalog (n=9,210)"),
        ("gk_mainshocks", "G-K Mainshocks (n=5,883)"),
        ("gk_aftershocks", "G-K Aftershocks (n=3,327)"),
    ]

    # Collect all p_cluster_robust values to find consistent y-axis range
    all_p = []
    for cat_key, _ in catalogs:
        s = results["subtest_4_declustering"][cat_key]["schuster_1day"]
        all_p.extend([r["p_cluster_robust"] for r in s])
    p_min = max(1e-50, min(p for p in all_p if p > 0))
    p_max = 1.5

    mark_periods = [75.6, 91.25, 182.5, 365.25]

    for ax_idx, (cat_key, panel_title) in enumerate(catalogs):
        ax = axes[ax_idx]
        s = results["subtest_4_declustering"][cat_key]["schuster_1day"]
        ep = results["subtest_4_declustering"][cat_key]["explicit_period_tests"]

        periods = np.array([r["period_days"] for r in s])
        p_std = np.array([r["p_standard"] for r in s])
        p_cr = np.array([r["p_cluster_robust"] for r in s])

        add_period_markers(ax, p_min, p_max)
        ax.plot(periods, p_std, color="gray", linewidth=0.8, alpha=0.6, label="Standard p-value")
        ax.plot(periods, p_cr, color=STEELBLUE, linewidth=1.5, label="Cluster-robust p-value (1d)")

        # Mark significance at key explicit periods
        sig_labels = {
            "quarter_year_756": 75.6,
            "quarter_year_912": 91.25,
            "half_year": 182.5,
            "annual": 365.25,
        }
        for ep_key, ep_period in sig_labels.items():
            ep_entry = ep.get(ep_key, {})
            p_cr_val = ep_entry.get("p_cluster_robust", 1.0)
            color = "red" if p_cr_val < 0.05 else "gray"
            marker = "o" if p_cr_val < 0.05 else "o"
            ax.scatter(ep_period, p_cr_val, color=color,
                       edgecolors="darkred" if p_cr_val < 0.05 else "gray",
                       s=50, zorder=6, marker=marker,
                       facecolors=color if p_cr_val < 0.05 else "none")
            ax.annotate(
                f"p={p_cr_val:.3f}",
                (ep_period, max(p_cr_val, p_min * 10)),
                fontsize=6.5, ha="center", va="bottom",
                xytext=(0, 5), textcoords="offset points",
            )

        format_schuster_axes(ax, panel_title, f"n={ep.get('half_year', {}).get('n_clusters', '?')} clusters at 182.5d")
        ax.set_ylim(p_min, p_max)
        ax.legend(fontsize=7, loc="upper right", ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = OUTPUT_DIR / "case-a3-a2-declustering-sensitivity.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Figure 5 — Cluster-window sensitivity
# ---------------------------------------------------------------------------

def figure5(results: dict) -> None:
    log.info("Figure 5: Cluster-window sensitivity")
    fig, ax = plt.subplots(1, 1, figsize=(11, 5))
    fig.suptitle(
        "Cluster-Window Sensitivity: 1-Day vs. 7-Day Threshold (A3.A2, Full Catalog)",
        fontsize=12, fontweight="bold",
    )

    s1 = results["subtest_1_baseline"]["schuster_1day"]
    s7 = results["subtest_1_baseline"]["schuster_7day"]

    periods = np.array([r["period_days"] for r in s1])
    p_std = np.array([r["p_standard"] for r in s1])
    p_cr1 = np.array([r["p_cluster_robust"] for r in s1])
    p_cr7 = np.array([r["p_cluster_robust"] for r in s7])

    add_period_markers(ax, 0, 1)
    ax.plot(periods, p_std, color="gray", linewidth=0.8, alpha=0.7, label="Standard p-value")
    ax.plot(periods, p_cr1, color=STEELBLUE, linewidth=1.5, linestyle="-", label="Cluster-robust 1-day window")
    ax.plot(periods, p_cr7, color=STEELBLUE, linewidth=1.5, linestyle="--", label="Cluster-robust 7-day window")

    inf_1d = results["subtest_1_baseline"]["inflation_factor_1day"]
    inf_7d = results["subtest_1_baseline"]["inflation_factor_7day"]

    annotation = (
        f"1-day inflation: {inf_1d:.1f}×\n"
        f"7-day inflation: {inf_7d:.1f}×\n"
        f"A2.A1 reference: 329×"
    )
    ax.text(0.98, 0.96, annotation, transform=ax.transAxes,
            ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Period (days)", fontsize=10)
    ax.set_ylabel("p-value", fontsize=10)
    ax.set_xlim(0.25, 548.0)
    ax.axhline(0.05, color="red", linestyle="--", linewidth=1.0, alpha=0.8, label="p=0.05")
    ax.legend(fontsize=8, loc="upper left")

    plt.tight_layout()
    out = OUTPUT_DIR / "case-a3-a2-cluster-window.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved %s", out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    results = load_results()
    figure1(results)
    figure2(results)
    figure3(results)
    figure4(results)
    figure5(results)
    log.info("All figures generated.")


if __name__ == "__main__":
    main()
