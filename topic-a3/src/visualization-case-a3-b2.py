"""
Case A3.B2: Hemisphere Stratification Refinement — Visualization

Generates 5 figures:
1. Tectonic × hemisphere signal heatmap (Figure 1)
2. Mid-crustal hemisphere bin distributions (Figure 2)
3. Phase alignment comparison (Figure 3)
4. Declustering sensitivity on tectonic-matched comparison (Figure 4)
5. Interval 1 SH threshold sensitivity (Figure 5)
"""

import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import numpy as np

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a3-b2-results.json"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Load results
# ---------------------------------------------------------------------------
with open(RESULTS_PATH) as fh:
    results = json.load(fh)

subtest1 = results["subtest_1_tectonic_hemisphere"]
subtest2 = results["subtest_2_midcrustal_hemisphere"]
subtest3 = results["subtest_3_phase_alignment"]
subtest4 = results["subtest_4_interval1_threshold"]
catalog_sizes = results["catalog_sizes"]

TECTONIC_CLASSES = ["continental", "transitional", "oceanic"]
CATALOGS = ["full", "gk", "reas"]
CATALOG_LABELS = {"full": "Full", "gk": "G-K", "reas": "Reas"}
HEMISPHERES = ["nh", "sh"]
HEMI_LABELS = {"nh": "NH", "sh": "SH"}

# ---------------------------------------------------------------------------
# Figure 1: Tectonic × hemisphere signal heatmap
# ---------------------------------------------------------------------------

def figure1_tectonic_heatmap() -> None:
    """
    3-row × 6-column heatmap: rows = tectonic class, columns = NH/SH per catalog.
    Cell fill = Cramér's V; white=0, red=max.
    """
    fig, ax = plt.subplots(figsize=(16, 7))
    fig.subplots_adjust(left=0.12, right=0.88, top=0.88, bottom=0.12)

    n_rows = 3  # tectonic classes
    n_cols = 6  # NH/SH for each of 3 catalogs

    # Collect all Cramér's V values for shared color scale
    all_cv: list[float] = []
    for cat in CATALOGS:
        for tclass in TECTONIC_CLASSES:
            for hemi in HEMISPHERES:
                cv = subtest1[cat][tclass][hemi].get("cramers_v")
                if cv is not None:
                    all_cv.append(cv)

    max_cv = max(all_cv) if all_cv else 0.05
    cmap = mcolors.LinearSegmentedColormap.from_list("white_red", ["white", "red"])
    norm = mcolors.Normalize(vmin=0.0, vmax=max_cv)

    # Column headers: Full-NH, Full-SH | GK-NH, GK-SH | Reas-NH, Reas-SH
    col_labels = ["Full NH", "Full SH", "G-K NH", "G-K SH", "Reas NH", "Reas SH"]

    for row_idx, tclass in enumerate(TECTONIC_CLASSES):
        col_idx_base = 0
        for cat in CATALOGS:
            for hemi in HEMISPHERES:
                cell = subtest1[cat][tclass][hemi]
                cv = cell.get("cramers_v")
                p = cell.get("p_chi2")
                chi2 = cell.get("chi2")
                n = cell.get("n", 0)
                low_n = cell.get("low_n", False) or (n < 100)

                x0 = col_idx_base / n_cols
                x1 = (col_idx_base + 1) / n_cols
                y0 = (n_rows - 1 - row_idx) / n_rows
                y1 = (n_rows - row_idx) / n_rows

                # Fill color
                if cv is not None and not low_n:
                    color = cmap(norm(cv))
                else:
                    color = "white"

                rect = plt.Rectangle(
                    (x0, y0), x1 - x0, y1 - y0,
                    transform=ax.transAxes,
                    facecolor=color, edgecolor="gray", linewidth=0.8,
                    zorder=1,
                )
                ax.add_patch(rect)

                # Gray diagonal striping for low-n
                if low_n:
                    for stripe in np.linspace(0, 1, 8):
                        ax.plot(
                            [x0 + stripe * (x1 - x0) * 0.3, x0 + stripe * (x1 - x0)],
                            [y1, y0 + stripe * (y1 - y0)],
                            color="lightgray", linewidth=0.5,
                            transform=ax.transAxes, zorder=2,
                        )

                # Bold border for p < 0.05
                if p is not None and p < 0.05:
                    border_rect = plt.Rectangle(
                        (x0, y0), x1 - x0, y1 - y0,
                        transform=ax.transAxes,
                        facecolor="none", edgecolor="black", linewidth=2.5,
                        zorder=3,
                    )
                    ax.add_patch(border_rect)

                # Text annotations
                cy = (y0 + y1) / 2
                cx = (x0 + x1) / 2

                if low_n:
                    ax.text(cx, cy, f"n={n}\n(low-n)", ha="center", va="center",
                            fontsize=7, transform=ax.transAxes, zorder=4, color="gray")
                elif cv is not None:
                    p_str = f"p={p:.3f}" if p >= 0.001 else f"p={p:.2e}"
                    chi2_str = f"χ²={chi2:.1f}" if chi2 is not None else ""
                    ax.text(cx, cy + 0.04, f"V={cv:.4f}", ha="center", va="center",
                            fontsize=7.5, transform=ax.transAxes, zorder=4, fontweight="bold")
                    ax.text(cx, cy - 0.01, chi2_str, ha="center", va="center",
                            fontsize=6.5, transform=ax.transAxes, zorder=4)
                    ax.text(cx, cy - 0.06, p_str, ha="center", va="center",
                            fontsize=6.5, transform=ax.transAxes, zorder=4)
                    ax.text(cx, cy - 0.11, f"n={n}", ha="center", va="center",
                            fontsize=6, transform=ax.transAxes, zorder=4, color="dimgray")

                col_idx_base += 1

    # Column headers
    for i, label in enumerate(col_labels):
        cx = (i + 0.5) / n_cols
        ax.text(cx, 1.04, label, ha="center", va="bottom", fontsize=9,
                transform=ax.transAxes, fontweight="bold")

    # Row labels (left side)
    for row_idx, tclass in enumerate(TECTONIC_CLASSES):
        cy = (n_rows - 0.5 - row_idx) / n_rows
        ax.text(-0.04, cy, tclass.capitalize(), ha="right", va="center",
                fontsize=10, transform=ax.transAxes, fontweight="bold", rotation=0)

    # Separator lines between catalog triplets
    for sep_x in [2 / n_cols, 4 / n_cols]:
        ax.axvline(sep_x, color="black", linewidth=1.5, zorder=5)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.65])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Cramér's V", fontsize=9)

    ax.set_title(
        "Tectonic-Matched Hemisphere χ² Signal Strength (A3.B2)",
        fontsize=13, fontweight="bold", pad=30,
    )

    out_path = OUTPUT_DIR / "case-a3-b2-tectonic-heatmap.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 1 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Mid-crustal hemisphere bin distributions
# ---------------------------------------------------------------------------

def figure2_midcrustal_binplots() -> None:
    """
    2-row × 3-column grid: rows=NH/SH, columns=full/GK/Reas.
    Horizontal steelblue bar charts with phase interval shading.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Mid-Crustal Solar Phase Distribution by Hemisphere (A3.B2)",
                 fontsize=13, fontweight="bold", y=1.02)

    row_labels = ["Northern Hemisphere (20–70 km)", "Southern Hemisphere (20–70 km)"]
    col_labels = ["Full catalog", "G-K mainshocks", "Reasenberg mainshocks"]
    hemi_keys = ["nh", "sh"]

    # Interval 1 shading: bins 4–5 at k=24
    interval_shade_k24 = {
        "interval_1": [4, 5],
        "interval_2": [15],
        "interval_3": [21],
    }

    for row_idx, hemi in enumerate(hemi_keys):
        for col_idx, cat in enumerate(CATALOGS):
            ax = axes[row_idx][col_idx]
            cell = subtest2[cat][hemi]
            n = cell.get("n", 0)
            k = cell.get("k")
            chi2 = cell.get("chi2")
            p = cell.get("p_chi2")
            cv = cell.get("cramers_v")
            bin_counts = cell.get("bin_counts")

            if bin_counts is None or k is None or n == 0:
                ax.text(0.5, 0.5, "n < 100\n(insufficient)", ha="center", va="center",
                        transform=ax.transAxes, fontsize=10, color="gray")
                ax.set_title(f"{CATALOG_LABELS[cat]}", fontsize=9)
                continue

            bc = np.array(bin_counts)
            expected = n / k
            expected_1sd = expected + np.sqrt(expected)
            bin_centers = [(b + 0.5) / k for b in range(k)]
            bin_phases = np.array(bin_centers)

            # Gray shading for interval bands (at k=24 spacing)
            if k == 24:
                for iname, ibins in interval_shade_k24.items():
                    for b in ibins:
                        if b < k:
                            x_start = b / k
                            x_end = (b + 1) / k
                            ax.axvspan(x_start, x_end, alpha=0.12, color="gray", zorder=0)

            # Horizontal bars — use bin phase as y-axis position
            bar_width = 1.0 / k
            colors = ["steelblue"] * k
            ax.barh(bin_phases, bc, height=bar_width * 0.85,
                    color=colors, edgecolor="white", linewidth=0.3, zorder=2)

            # Expected line
            ax.axvline(expected, color="orange", linestyle="--", linewidth=1.2,
                       label=f"Expected ({expected:.1f})", zorder=3)

            # 1-SD threshold
            ax.axvline(expected_1sd, color="red", linestyle=":", linewidth=1.0,
                       label=f"+1 SD ({expected_1sd:.1f})", zorder=3)

            # Peak bin marker
            peak_bin = int(np.argmax(bc))
            peak_phase_center = (peak_bin + 0.5) / k
            ax.axhline(peak_phase_center, color="navy", linestyle=":",
                       linewidth=1.0, alpha=0.7, zorder=3)

            # Annotations
            p_str = f"p={p:.3f}" if p >= 0.001 else f"p={p:.2e}"
            annot = f"n={n}, k={k}\nχ²={chi2:.2f}\n{p_str}\nV={cv:.4f}"
            ax.text(0.97, 0.03, annot, transform=ax.transAxes, fontsize=7.5,
                    ha="right", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

            ax.set_xlabel("Count", fontsize=8)
            ax.set_ylabel("Phase (solar year)", fontsize=8)
            ax.set_ylim(-0.02, 1.02)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.tick_params(labelsize=7)
            ax.set_title(f"{CATALOG_LABELS[cat]}", fontsize=9)

    # Row and column outer labels
    for row_idx, label in enumerate(row_labels):
        axes[row_idx][0].set_ylabel(f"{label}\n\nPhase", fontsize=8, labelpad=5)
    for col_idx, label in enumerate(col_labels):
        axes[0][col_idx].set_title(label, fontsize=9, fontweight="bold")

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b2-midcrustal-binplots.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 2 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Phase alignment comparison
# ---------------------------------------------------------------------------

def figure3_phase_alignment() -> None:
    """
    Scatter plot: x=NH peak phase, y=SH peak phase. Color by source type.
    Reference lines: y=x (in-phase), y=x±0.5 (anti-phase).
    """
    fig, ax = plt.subplots(figsize=(8, 7))

    pairs = subtest3["pairs"]

    # Color by source type
    source_colors = {
        "subtest1_continental": "#e41a1c",   # red
        "subtest1_transitional": "#4daf4a",  # green
        "subtest1_oceanic": "#984ea3",       # purple
        "subtest2_midcrustal": "#377eb8",    # blue
    }
    source_labels = {
        "subtest1_continental": "Continental (ST1)",
        "subtest1_transitional": "Transitional (ST1)",
        "subtest1_oceanic": "Oceanic (ST1)",
        "subtest2_midcrustal": "Mid-crustal (ST2)",
    }

    # Short labels for points
    label_map_cat = {"full": "Full", "gk": "GK", "reas": "Reas"}
    label_map_src = {
        "subtest1_continental": "Cont",
        "subtest1_transitional": "Trans",
        "subtest1_oceanic": "Oce",
        "subtest2_midcrustal": "Mid",
    }

    # Reference lines: y=x (in-phase)
    x_ref = np.linspace(0, 1, 200)
    ax.plot(x_ref, x_ref, color="gray", linewidth=1.2, linestyle="-",
            label="In-phase (y=x)", zorder=1)
    ax.plot(x_ref, (x_ref + 0.5) % 1.0, color="gray", linewidth=1.0, linestyle="--",
            label="Anti-phase (y=x±0.5)", zorder=1)
    ax.plot(x_ref, (x_ref - 0.5) % 1.0, color="gray", linewidth=1.0, linestyle="--",
            zorder=1)

    # Tolerance bands: ±2 bins (0.083 phase units) around y=x
    ax.fill_between(x_ref, x_ref - 0.083, x_ref + 0.083,
                    alpha=0.15, color="steelblue", label="±2-bin tolerance (in-phase)", zorder=0)
    # Anti-phase band: ±2 bins around y=x±0.5
    ax.fill_between(x_ref, (x_ref + 0.5) % 1.0 - 0.083, (x_ref + 0.5) % 1.0 + 0.083,
                    alpha=0.15, color="orange", label="±2-bin tolerance (anti-phase)", zorder=0)
    ax.fill_between(x_ref, (x_ref - 0.5) % 1.0 - 0.083, (x_ref - 0.5) % 1.0 + 0.083,
                    alpha=0.15, color="orange", zorder=0)

    plotted_sources = set()
    for pair in pairs:
        src = pair["source"]
        cat = pair["catalog"]
        nh_phase = pair["nh_peak_phase"]
        sh_phase = pair["sh_peak_phase"]
        nh_sig = pair["nh_significant"]
        sh_sig = pair["sh_significant"]
        both_sig = nh_sig and sh_sig

        color = source_colors.get(src, "black")
        marker = "o" if both_sig else "o"
        fill = "full" if both_sig else "none"
        mec = color

        label_src = source_labels.get(src, src) if src not in plotted_sources else "_nolegend_"
        if src not in plotted_sources:
            plotted_sources.add(src)

        ax.scatter(
            nh_phase, sh_phase,
            color=color if fill == "full" else "none",
            edgecolors=mec,
            marker=marker,
            s=80,
            linewidths=1.5,
            zorder=4,
            label=label_src,
        )

        # Point label
        short_label = f"{label_map_src.get(src, src)}-{label_map_cat.get(cat, cat)}"
        ax.annotate(
            short_label,
            xy=(nh_phase, sh_phase),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
            color=color,
            zorder=5,
        )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("NH Peak Phase (solar year)", fontsize=10)
    ax.set_ylabel("SH Peak Phase (solar year)", fontsize=10)
    ax.set_title("NH vs. SH Peak Phase — Alignment Test (A3.B2)", fontsize=12, fontweight="bold")

    # Legend
    dominant = subtest3["dominant_alignment"]
    n_pairs = subtest3["n_pairs_evaluated"]
    n_in = subtest3["n_in_phase"]
    n_anti = subtest3["n_anti_phase"]
    n_off = subtest3["n_offset"]
    annot = (f"n pairs={n_pairs}\nIn-phase: {n_in}\nAnti-phase: {n_anti}\n"
             f"Offset: {n_off}\nDominant: {dominant}")
    ax.text(0.02, 0.97, annot, transform=ax.transAxes, fontsize=8,
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85))

    ax.legend(loc="lower right", fontsize=7, framealpha=0.85)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b2-phase-alignment.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 3 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 4: Declustering sensitivity on tectonic-matched comparison
# ---------------------------------------------------------------------------

def figure4_declustering_sensitivity() -> None:
    """
    2-row × 3-column grouped bar chart: rows=NH/SH, columns=tectonic class.
    x=catalog, y=Cramér's V. Asterisks for p<0.05 (*) and p<0.01 (**).
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharey=False)
    fig.suptitle("Declustering Sensitivity — Tectonic × Hemisphere (A3.B2)",
                 fontsize=13, fontweight="bold")

    cat_labels = ["Full", "G-K", "Reasenberg"]
    x_pos = np.arange(len(CATALOGS))
    bar_width = 0.55

    hemi_titles = ["Northern Hemisphere", "Southern Hemisphere"]

    for row_idx, hemi in enumerate(["nh", "sh"]):
        for col_idx, tclass in enumerate(TECTONIC_CLASSES):
            ax = axes[row_idx][col_idx]

            cvs = []
            ps = []
            ns_list = []
            ref_cv = None  # full catalog result

            for cat in CATALOGS:
                cell = subtest1[cat][tclass][hemi]
                cv = cell.get("cramers_v") or 0.0
                p = cell.get("p_chi2")
                n = cell.get("n", 0)
                cvs.append(cv)
                ps.append(p)
                ns_list.append(n)
                if cat == "full":
                    ref_cv = cv

            bars = ax.bar(x_pos, cvs, width=bar_width, color="steelblue",
                          edgecolor="white", linewidth=0.5, zorder=2)

            # Reference line: full catalog Cramér's V
            if ref_cv is not None:
                ax.axhline(ref_cv, color="navy", linewidth=1.2, linestyle="--",
                           alpha=0.6, zorder=3)

            # Annotate bars with p-value significance and n
            for bi, (bar, p_val, n_val) in enumerate(zip(bars, ps, ns_list)):
                if p_val is None:
                    marker = ""
                elif p_val < 0.01:
                    marker = "**"
                elif p_val < 0.05:
                    marker = "*"
                else:
                    marker = ""
                bar_top = bar.get_height()
                if marker:
                    ax.text(bar.get_x() + bar.get_width() / 2, bar_top + 0.0002,
                            marker, ha="center", va="bottom", fontsize=11, color="darkred",
                            fontweight="bold")
                p_label = f"p={p_val:.3f}" if p_val is not None and p_val >= 0.001 else (f"p={p_val:.2e}" if p_val is not None else "n/a")
                ax.text(bar.get_x() + bar.get_width() / 2, bar_top / 2,
                        f"n={n_val}", ha="center", va="center", fontsize=6.5,
                        color="white", fontweight="bold")

            ax.set_xticks(x_pos)
            ax.set_xticklabels(cat_labels, fontsize=8)
            ax.set_ylabel("Cramér's V" if col_idx == 0 else "", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(axis="y", alpha=0.3, linewidth=0.5)

            if row_idx == 0:
                ax.set_title(tclass.capitalize(), fontsize=10, fontweight="bold")
            if col_idx == 0:
                ax.set_ylabel(f"{hemi_titles[row_idx]}\nCramér's V", fontsize=8)

    # Legend
    legend_elements = [
        mpatches.Patch(color="steelblue", label="Cramér's V"),
        plt.Line2D([0], [0], color="navy", linestyle="--", linewidth=1.2, label="Full catalog baseline"),
        mpatches.Patch(color="white", label="* p<0.05,  ** p<0.01", linewidth=0),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3, fontsize=8,
               bbox_to_anchor=(0.5, -0.01), framealpha=0.9)

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b2-declustering-sensitivity.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 4 saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 5: Interval 1 SH threshold sensitivity
# ---------------------------------------------------------------------------

def figure5_interval1_threshold() -> None:
    """
    1-row × 2-column: all SH events and continental SH events only.
    Each panel: x=threshold, y=overlap_fraction; threshold lines; present/absent fill.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle("Interval 1 SH Presence — Threshold Sensitivity (A3.B2)",
                 fontsize=13, fontweight="bold")

    panel_keys = ["sh_all", "sh_continental"]
    panel_titles = ["All SH Events", "Continental SH Events Only"]

    threshold_colors = {
        0.33: "#e41a1c",
        0.40: "#ff7f00",
        0.45: "#4daf4a",
        0.50: "#377eb8",
    }

    for ax, pop_key, panel_title in zip(axes, panel_keys, panel_titles):
        records = subtest4[pop_key]
        thresholds = [r["threshold"] for r in records]
        overlap_fractions = [r["overlap_fraction"] for r in records]
        classifications = [r["interval_1_classification"] for r in records]
        n_sh = records[0]["n_sh"] if records else 0

        # The overlap fraction is constant (does not vary with threshold)
        of_val = overlap_fractions[0] if overlap_fractions else 0.0

        # Threshold dashed lines
        for t, tcolor in threshold_colors.items():
            ax.axhline(t, color=tcolor, linestyle="--", linewidth=1.2,
                       alpha=0.75, label=f"Threshold={t:.2f}")

        # Horizontal line for the actual overlap fraction
        ax.axhline(of_val, color="black", linewidth=2.0, zorder=5, label=f"Overlap fraction={of_val:.3f}")

        # Fill below/above
        x_fill = [0.30, 0.55]
        for rec in records:
            t = rec["threshold"]
            classification = rec["interval_1_classification"]
            fill_color = "steelblue" if classification == "present" else "lightgray"
            # Mark the region around each threshold
            ax.scatter(t, of_val, color=fill_color, s=100, zorder=6,
                       edgecolors="black", linewidth=1.0)
            ax.annotate(
                classification,
                xy=(t, of_val),
                xytext=(0, 12),
                textcoords="offset points",
                ha="center", va="bottom",
                fontsize=9, fontweight="bold",
                color="steelblue" if classification == "present" else "gray",
            )

        ax.set_xlim(0.28, 0.55)
        ax.set_ylim(-0.05, 1.1)
        ax.set_xticks(thresholds)
        ax.set_xticklabels([f"{t:.2f}" for t in thresholds], fontsize=9)
        ax.set_xlabel("Overlap Threshold", fontsize=10)
        ax.set_ylabel("Interval 1 Overlap Fraction", fontsize=10)
        ax.set_title(f"{panel_title}\n(n={n_sh})", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7.5, loc="upper right", framealpha=0.85)
        ax.grid(True, alpha=0.3, linewidth=0.5)

        # Flip point annotation
        flip_key = f"{pop_key}_flip_threshold"
        flip_val = subtest4.get(flip_key)
        if flip_val is not None:
            ax.annotate(
                f"Flip at t={flip_val:.2f}",
                xy=(flip_val, of_val),
                xytext=(0, -25),
                textcoords="offset points",
                ha="center", fontsize=8, color="darkred",
                arrowprops=dict(arrowstyle="->", color="darkred"),
            )
        else:
            ax.text(0.5, 0.08, "No classification flip", transform=ax.transAxes,
                    ha="center", fontsize=8.5, color="dimgray", style="italic")

    fig.tight_layout()
    out_path = OUTPUT_DIR / "case-a3-b2-interval1-threshold.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 5 saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Generate all A3.B2 figures."""
    logger.info("Generating Figure 1: Tectonic heatmap...")
    figure1_tectonic_heatmap()

    logger.info("Generating Figure 2: Mid-crustal bin plots...")
    figure2_midcrustal_binplots()

    logger.info("Generating Figure 3: Phase alignment...")
    figure3_phase_alignment()

    logger.info("Generating Figure 4: Declustering sensitivity...")
    figure4_declustering_sensitivity()

    logger.info("Generating Figure 5: Interval 1 threshold sensitivity...")
    figure5_interval1_threshold()

    logger.info("All figures generated.")


if __name__ == "__main__":
    main()
