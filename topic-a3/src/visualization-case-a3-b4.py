"""
Case A3.B4: Visualizations — Depth × Magnitude Two-Way Stratification with Moho Isolation

Generates:
  - Figure 1: Two-way stratification matrix (p-value heatmap)
  - Figure 2: Cramér's V trends by depth × magnitude band
  - Figure 3: Moho isolation bin distributions (delta = 5, 10, 15 km)
  - Figure 4: Subduction proximity cross-tabulation (mid-crustal)
  - Figure 5: Global Moho depth map with event overlay (Cartopy Robinson)

CRUST1.0 data from https://igppweb.ucsd.edu/~gabi/crust1.html (Laske et al. 2013)
"""

import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = BASE_DIR / "output"
RESULTS_PATH = OUTPUT_DIR / "case-a3-b4-results.json"
EVENTS_PKL   = OUTPUT_DIR / "case-a3-b4-events.pkl"
MOHO_GRID_NPY = OUTPUT_DIR / "case-a3-b4-moho-grid.npy"

# Output PNG paths
FIG1_PATH = OUTPUT_DIR / "case-a3-b4-stratification-matrix.png"
FIG2_PATH = OUTPUT_DIR / "case-a3-b4-cramer-trends.png"
FIG3_PATH = OUTPUT_DIR / "case-a3-b4-moho-isolation.png"
FIG4_PATH = OUTPUT_DIR / "case-a3-b4-subduction-crosstab.png"
FIG5_PATH = OUTPUT_DIR / "case-a3-b4-moho-map.png"

DEPTH_LABELS = [
    "shallow_0-20km",
    "midcrustal_20-70km",
    "intermediate_70-300km",
    "deep_300km+",
]
DEPTH_DISPLAY = ["Shallow\n0–20 km", "Mid-crustal\n20–70 km", "Intermediate\n70–300 km", "Deep\n>300 km"]
MAG_LABELS  = ["m6_6.9", "m7_7.9", "m8_plus"]
MAG_DISPLAY = ["M6.0–6.9", "M7.0–7.9", "M8.0+"]

A2B4_REF = {
    "shallow_0-20km":        {"p_chi2": 0.368, "cramers_v": 0.0187},
    "midcrustal_20-70km":    {"p_chi2": 4.02e-9, "cramers_v": 0.0285},
    "intermediate_70-300km": {"p_chi2": 0.807,  "cramers_v": 0.0268},
    "deep_300km+":            {"p_chi2": 0.626,  "cramers_v": 0.0398},
}


def load_results() -> dict:
    """Load analysis results JSON."""
    with open(RESULTS_PATH) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Figure 1: Two-way stratification matrix
# ---------------------------------------------------------------------------

def figure1_stratification_matrix(results: dict) -> None:
    """P-value heatmap: 4 depth rows × 3 magnitude columns."""
    logger.info("Generating Figure 1: stratification matrix...")

    sm = results["stratification_matrix"]
    db_totals = results["depth_band_totals"]

    n_rows, n_cols = 4, 3
    fig, axes = plt.subplots(
        n_rows, n_cols + 1,
        figsize=(14, 10),
        gridspec_kw={"width_ratios": [1, 1, 1, 0.45]},
    )

    # Color map: white (p=1) → red (p=0)
    cmap = mcolors.LinearSegmentedColormap.from_list("wred", ["white", "red"])
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

    for r, d_label in enumerate(DEPTH_LABELS):
        for c, m_label in enumerate(MAG_LABELS):
            ax = axes[r, c]
            cell = sm[d_label][m_label]
            p_val = cell["p_chi2"]
            n_val = cell["n"]
            k_val = cell["k"]
            low_n = cell["low_n"]

            color = cmap(1.0 - p_val)
            ax.set_facecolor(color)

            # Bold border if significant
            border_lw = 3.0 if p_val < 0.05 else 0.5
            for spine in ax.spines.values():
                spine.set_linewidth(border_lw)
                spine.set_edgecolor("black" if p_val < 0.05 else "#aaaaaa")

            # Text annotations
            p_str = f"{p_val:.2e}"
            ax.text(0.5, 0.72, f"n = {n_val:,}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, fontweight="bold")
            ax.text(0.5, 0.45, f"p = {p_str}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9)
            ax.text(0.5, 0.20, f"k = {k_val}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8, color="gray")

            if low_n:
                ax.text(0.5, 0.85, "LOW-N", transform=ax.transAxes,
                        ha="center", va="center", fontsize=7,
                        color="darkred", alpha=0.6, fontweight="bold")

            ax.set_xticks([])
            ax.set_yticks([])

            # Column headers (top row)
            if r == 0:
                ax.set_title(MAG_DISPLAY[c], fontsize=11, fontweight="bold", pad=6)

        # Row labels (left)
        axes[r, 0].set_ylabel(DEPTH_DISPLAY[r], fontsize=9, labelpad=8)

        # Reference column: A2.B4 totals
        ax_ref = axes[r, n_cols]
        # Use depth band total p from current run
        total_cell = db_totals[d_label]
        p_total = total_cell["p_chi2"]
        ref_p = A2B4_REF[d_label]["p_chi2"]
        color_ref = cmap(1.0 - p_total)
        ax_ref.set_facecolor(color_ref)
        for spine in ax_ref.spines.values():
            spine.set_linewidth(2.0 if p_total < 0.05 else 0.5)
            spine.set_edgecolor("black" if p_total < 0.05 else "#aaaaaa")

        ax_ref.text(0.5, 0.70, f"n={total_cell['n']:,}", transform=ax_ref.transAxes,
                    ha="center", va="center", fontsize=8, fontweight="bold")
        ax_ref.text(0.5, 0.45, f"p={p_total:.2e}", transform=ax_ref.transAxes,
                    ha="center", va="center", fontsize=7)
        ax_ref.text(0.5, 0.22, f"A2.B4\np={ref_p:.2e}", transform=ax_ref.transAxes,
                    ha="center", va="center", fontsize=6, color="#555555")
        ax_ref.set_xticks([])
        ax_ref.set_yticks([])

        if r == 0:
            ax_ref.set_title("All Mag\n(ref A2.B4)", fontsize=9, fontweight="bold", pad=6)

    fig.suptitle(
        "Solar Phase Chi-Square P-Value by Depth × Magnitude",
        fontsize=13, fontweight="bold", y=1.01,
    )

    # Colorbar
    sm_cb = ScalarMappable(cmap=cmap, norm=norm)
    sm_cb.set_array([])
    cbar = fig.colorbar(sm_cb, ax=axes[:, :], orientation="vertical",
                        shrink=0.6, pad=0.02, fraction=0.02)
    cbar.set_label("Chi-square p-value", fontsize=10)
    cbar.ax.invert_yaxis()  # red at bottom (p=0), white at top (p=1)

    plt.tight_layout()
    fig.savefig(FIG1_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 1 saved: {FIG1_PATH}")


# ---------------------------------------------------------------------------
# Figure 2: Cramér's V trends
# ---------------------------------------------------------------------------

def figure2_cramer_trends(results: dict) -> None:
    """Cramér's V across depth bands, one line per magnitude band."""
    logger.info("Generating Figure 2: Cramér's V trends...")

    sm = results["stratification_matrix"]
    db_totals = results["depth_band_totals"]

    fig, ax = plt.subplots(figsize=(10, 6))

    x_pos = np.arange(1, 5)  # ordinal 1–4

    mag_styles = [
        ("m6_6.9",  "M6.0–6.9", "steelblue", "-",  "o"),
        ("m7_7.9",  "M7.0–7.9", "steelblue", "--", "o"),
        ("m8_plus", "M8.0+",    "steelblue", ":",  "o"),
    ]

    for m_label, m_display, color, ls, mk in mag_styles:
        v_vals, p_vals, n_vals, low_ns = [], [], [], []
        for d_label in DEPTH_LABELS:
            cell = sm[d_label][m_label]
            v_vals.append(cell["cramers_v"])
            p_vals.append(cell["p_chi2"])
            n_vals.append(cell["n"])
            low_ns.append(cell["low_n"])

        ax.plot(x_pos, v_vals, color=color, linestyle=ls, linewidth=1.5,
                label=m_display, zorder=2)

        for xi, (v, p, n, low_n) in enumerate(zip(v_vals, p_vals, n_vals, low_ns)):
            if low_n:
                marker_shape = "x"
            elif p < 0.05:
                marker_shape = "o"  # filled
            else:
                marker_shape = "o"  # open

            fill = "full" if (p < 0.05 and not low_n) else "none" if not low_n else "none"
            ax.plot(xi + 1, v,
                    marker=marker_shape,
                    color=color,
                    fillstyle=fill,
                    markersize=9 if not low_n else 8,
                    zorder=3)
            ax.text(xi + 1, v + 0.0008, f"n={n:,}", ha="center", va="bottom",
                    fontsize=7, color="gray")

    # Depth-band totals line (red solid)
    total_v = [db_totals[d]["cramers_v"] for d in DEPTH_LABELS]
    total_p = [db_totals[d]["p_chi2"] for d in DEPTH_LABELS]
    total_n = [db_totals[d]["n"] for d in DEPTH_LABELS]
    ax.plot(x_pos, total_v, color="red", linestyle="-", linewidth=2,
            label="All magnitudes (depth total)", zorder=4)
    for xi, (v, p, n) in enumerate(zip(x_pos, total_v, total_p)):
        fill = "full" if p < 0.05 else "none"
        ax.plot(xi + 1, v, marker="o", color="red", fillstyle=fill, markersize=10, zorder=5)

    # Global reference line — use midcrustal all-mag V from A2.B4
    ax.axhline(0.0285, color="gray", linestyle="--", linewidth=1, alpha=0.6,
               label="A2.B4 global V (mid-crustal all-mag)")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(DEPTH_DISPLAY, fontsize=10)
    ax.set_xlabel("Depth Band", fontsize=11)
    ax.set_ylabel("Cramér's V", fontsize=11)
    ax.set_title("Cramér's V by Depth Band and Magnitude Band", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG2_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 2 saved: {FIG2_PATH}")


# ---------------------------------------------------------------------------
# Figure 3: Moho isolation bin distributions
# ---------------------------------------------------------------------------

def figure3_moho_isolation(results: dict) -> None:
    """3×2 panel: bin distributions and continental vs oceanic comparison per delta."""
    logger.info("Generating Figure 3: Moho isolation distributions...")

    mi = results["moho_isolation"]
    delta_keys = ["delta_5km", "delta_10km", "delta_15km"]
    delta_labels = ["Δ = 5 km", "Δ = 10 km", "Δ = 15 km"]

    # A1b interval shading bins for k=12 (approximate; spec says use gray shaded bands)
    # k=12 interval boundaries: interval_1 ~ bins 2, interval_2 ~ bin 7, interval_3 ~ bin 10
    # For display we use k=12 equivalent of the k=24 intervals
    # interval_1: bins 4-5 at k=24 → bins 2 at k=12 (approx)
    # interval_2: bin 15 at k=24 → bin 7 at k=12
    # interval_3: bin 21 at k=24 → bin 10 at k=12

    fig, axes = plt.subplots(2, 3, figsize=(18, 8))

    for col_idx, (dk, dlabel) in enumerate(zip(delta_keys, delta_labels)):
        data = mi[dk]
        n_total = data["n_total"]
        k_use   = data["k"]
        chi2    = data["chi2"]
        p_val   = data["p_chi2"]
        v_val   = data["cramers_v"]
        bins    = data.get("bin_counts", [])

        # ---- Top row: overall bin distribution ----
        ax_top = axes[0, col_idx]
        if bins:
            expected = n_total / k_use
            x_bins = np.arange(k_use)
            ax_top.barh(x_bins, bins, color="steelblue", height=0.8, alpha=0.8)
            ax_top.axvline(expected, color="black", linestyle="--", linewidth=1.2,
                           label=f"Expected ({expected:.1f})")

            # Gray shaded bands for A1b intervals (approximate for k used)
            if k_use == 24:
                for b in [4, 5]:
                    ax_top.axhspan(b - 0.4, b + 0.4, color="gray", alpha=0.15)
                ax_top.axhspan(15 - 0.4, 15 + 0.4, color="gray", alpha=0.15)
                ax_top.axhspan(21 - 0.4, 21 + 0.4, color="gray", alpha=0.15)
            elif k_use == 16:
                # approximate conversions
                ax_top.axhspan(2 - 0.4, 3 + 0.4, color="gray", alpha=0.15)
                ax_top.axhspan(10 - 0.4, 10 + 0.4, color="gray", alpha=0.15)
                ax_top.axhspan(14 - 0.4, 14 + 0.4, color="gray", alpha=0.15)

        ax_top.set_title(dlabel, fontsize=12, fontweight="bold")
        ax_top.set_xlabel("Count", fontsize=9)
        ax_top.set_ylabel("Phase bin", fontsize=9)
        ax_top.text(0.97, 0.97,
                    f"n={n_total}\nk={k_use}\nχ²={chi2:.2f}\np={p_val:.3e}\nV={v_val:.4f}",
                    transform=ax_top.transAxes, va="top", ha="right",
                    fontsize=8, bbox=dict(boxstyle="round", fc="white", alpha=0.7))
        ax_top.text(0.03, 0.03, "CRUST1.0 (1°×1°)",
                    transform=ax_top.transAxes, fontsize=7, color="gray")

        # ---- Bottom row: continental vs oceanic comparison ----
        ax_bot = axes[1, col_idx]
        cont_data = data.get("continental", {})
        oce_data  = data.get("oceanic", {})

        cont_bins = cont_data.get("bin_counts", [])
        oce_bins  = oce_data.get("bin_counts", [])
        n_cont = cont_data.get("n", 0)
        n_oce  = oce_data.get("n", 0)

        # Normalize by k to get comparable counts per bin
        k_cont = cont_data.get("k", 12)
        k_oce  = oce_data.get("k", 12)

        # Use the smaller k for display or just show raw counts
        k_disp = min(k_cont, k_oce) if cont_bins and oce_bins else 12

        # Trim or pad to k_disp bins
        def trim_bins(b: list, k: int) -> np.ndarray:
            arr = np.array(b, dtype=float)
            if len(arr) == 0:
                return np.zeros(k)
            # Rebin if needed
            if len(arr) == k:
                return arr
            # Just truncate or zero-pad
            result = np.zeros(k)
            result[:min(len(arr), k)] = arr[:min(len(arr), k)]
            return result

        cont_arr = trim_bins(cont_bins, k_disp)
        oce_arr  = trim_bins(oce_bins, k_disp)

        x = np.arange(k_disp)
        width = 0.4
        ax_bot.bar(x - width / 2, cont_arr, width=width, color="red",    alpha=0.7, label=f"Continental (n={n_cont})")
        ax_bot.bar(x + width / 2, oce_arr,  width=width, color="steelblue", alpha=0.7, label=f"Oceanic (n={n_oce})")

        ax_bot.set_xlabel("Phase bin", fontsize=9)
        ax_bot.set_ylabel("Count", fontsize=9)
        ax_bot.legend(fontsize=8)
        ax_bot.text(0.97, 0.97, f"k={k_disp} (display)",
                    transform=ax_bot.transAxes, va="top", ha="right",
                    fontsize=7, color="gray")

    fig.suptitle("Moho-Proximal Solar Phase Distributions by Delta and Tectonic Setting",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG3_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 3 saved: {FIG3_PATH}")


# ---------------------------------------------------------------------------
# Figure 4: Subduction proximity cross-tabulation
# ---------------------------------------------------------------------------

def figure4_subduction_crosstab(results: dict) -> None:
    """2-row: p-value (log) and Cramér's V for near vs far subduction per mag band."""
    logger.info("Generating Figure 4: subduction cross-tabulation...")

    ct = results["subduction_crosstab_midcrustal"]
    display_keys = ["m6_6.9", "m7_7.9", "m8_plus", "all_magnitudes"]
    display_labels = ["M6.0–6.9", "M7.0–7.9", "M8.0+", "All mag"]

    near_p = []
    far_p  = []
    near_v = []
    far_v  = []
    near_n = []
    far_n  = []

    for dk in display_keys:
        cell = ct[dk]
        near_p.append(cell["near_sub"]["p_chi2"])
        far_p.append(cell["far_sub"]["p_chi2"])
        near_v.append(cell["near_sub"]["cramers_v"])
        far_v.append(cell["far_sub"]["cramers_v"])
        near_n.append(cell["near_sub"]["n"])
        far_n.append(cell["far_sub"]["n"])

    x = np.arange(len(display_keys))
    width = 0.35

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Row 1: p-value (log scale)
    near_p_clip = [max(p, 1e-15) for p in near_p]
    far_p_clip  = [max(p, 1e-15) for p in far_p]

    bars1n = ax1.bar(x - width / 2, near_p_clip, width=width,
                     color="orange", alpha=0.8, label="Near-subduction (≤200 km)")
    bars1f = ax1.bar(x + width / 2, far_p_clip,  width=width,
                     color="steelblue", alpha=0.8, label="Far-subduction (>200 km)")
    ax1.set_yscale("log")
    ax1.axhline(0.05, color="black", linestyle="--", linewidth=1.2, label="p = 0.05")
    ax1.set_ylabel("Chi-square p-value (log scale)", fontsize=10)
    ax1.legend(fontsize=9)
    ax1.set_title(
        "Mid-Crustal Solar Signal: Near-Subduction vs. Far-Subduction\n(PB2002 SUB+OCB, ≤200 km)",
        fontsize=12, fontweight="bold"
    )

    # Annotate n above bars
    for i, (nn, nf) in enumerate(zip(near_n, far_n)):
        ax1.text(x[i] - width / 2, near_p_clip[i] * 1.5, f"n={nn}", ha="center", fontsize=8)
        ax1.text(x[i] + width / 2, far_p_clip[i] * 1.5,  f"n={nf}", ha="center", fontsize=8)

    # Row 2: Cramér's V
    ax2.bar(x - width / 2, near_v, width=width,
            color="orange", alpha=0.8, label="Near-subduction")
    ax2.bar(x + width / 2, far_v,  width=width,
            color="steelblue", alpha=0.8, label="Far-subduction")
    # Reference V from A3.B3 baseline (full catalog V)
    ax2.axhline(0.0285, color="gray", linestyle="--", linewidth=1.2,
                label="A2.B4 global V reference (0.0285)")
    ax2.set_ylabel("Cramér's V", fontsize=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(display_labels, fontsize=11)
    ax2.set_xlabel("Magnitude Band (mid-crustal depth, 20–70 km)", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.set_ylim(bottom=0)

    plt.tight_layout()
    fig.savefig(FIG4_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 4 saved: {FIG4_PATH}")


# ---------------------------------------------------------------------------
# Figure 5: Moho depth global map (Cartopy Robinson)
# ---------------------------------------------------------------------------

def figure5_moho_map(results: dict) -> None:
    """Global map: CRUST1.0 Moho depth background + events colored by Moho-proximity."""
    logger.info("Generating Figure 5: global Moho map (Cartopy Robinson)...")

    # Load moho grid and events
    moho_grid = np.load(MOHO_GRID_NPY)  # shape (180, 360)
    df = pd.read_pickle(EVENTS_PKL)

    # Grid axes
    lats = np.arange(89.5, -90.5, -1.0)   # 89.5 → -89.5 (descending)
    lons = np.arange(-179.5, 180.5, 1.0)  # -179.5 → 179.5

    # Subsample every 4th point for performance
    step = 4
    lats_sub = lats[::step]
    lons_sub = lons[::step]
    moho_sub = moho_grid[::step, ::step]

    # Build mesh grid for pcolormesh (need cell edges, not centers)
    lon_edges = np.linspace(-180, 180, moho_sub.shape[1] + 1)
    lat_edges = np.linspace(-90, 90, moho_sub.shape[0] + 1)
    lon_mesh, lat_mesh = np.meshgrid(lon_edges, lat_edges)

    # Flip moho_sub to match lat_edges (ascending lat_edges vs descending lats_sub)
    moho_plot = moho_sub[::-1, :]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()

    # White (shallow) → blue (deep) colormap for Moho depth
    cmap_moho = mcolors.LinearSegmentedColormap.from_list("wblue", ["white", "steelblue"])

    im = ax.pcolormesh(
        lon_mesh, lat_mesh, moho_plot,
        transform=ccrs.PlateCarree(),
        cmap=cmap_moho,
        vmin=3.0, vmax=65.0,
        rasterized=True,
        zorder=1,
    )

    # Coastlines
    ax.add_feature(cfeature.COASTLINE.with_scale("110m"),
                   linewidth=0.4, edgecolor="black", zorder=3)

    # Event overlay: delta=10 km classification
    proximal_mask = df["moho_proximal_10km"].values
    non_prox_mask = ~proximal_mask

    # Only events with valid depth
    has_depth = df["depth"].notna().values
    prox_idx = np.where(proximal_mask & has_depth)[0]
    non_prox_idx = np.where(non_prox_mask & has_depth)[0]

    def event_size(mag_arr: np.ndarray) -> np.ndarray:
        return (mag_arr - 5.5) ** 2 * 2.0

    # Non-proximal (gray, behind)
    if len(non_prox_idx) > 0:
        mags_np = df["usgs_mag"].values[non_prox_idx]
        lons_np = df["longitude"].values[non_prox_idx]
        lats_np = df["latitude"].values[non_prox_idx]
        ax.scatter(
            lons_np, lats_np,
            c="gray", s=event_size(mags_np),
            alpha=0.2, transform=ccrs.PlateCarree(),
            zorder=4, linewidths=0, rasterized=True,
        )

    # Moho-proximal (orange, front)
    if len(prox_idx) > 0:
        mags_pr = df["usgs_mag"].values[prox_idx]
        lons_pr = df["longitude"].values[prox_idx]
        lats_pr = df["latitude"].values[prox_idx]
        ax.scatter(
            lons_pr, lats_pr,
            c="orange", s=event_size(mags_pr),
            alpha=0.7, transform=ccrs.PlateCarree(),
            zorder=5, linewidths=0.2, edgecolors="darkorange", rasterized=True,
        )

    # Legend
    orange_patch = mpatches.Patch(color="orange", alpha=0.8,
                                  label=f"Moho-proximal (Δ≤10 km, n={len(prox_idx):,})")
    gray_patch   = mpatches.Patch(color="gray",   alpha=0.5,
                                  label=f"Other events (n={len(non_prox_idx):,})")
    ax.legend(handles=[orange_patch, gray_patch], loc="lower left",
              fontsize=9, framealpha=0.8)

    # Colorbar for Moho depth
    cbar = fig.colorbar(im, ax=ax, orientation="vertical",
                        shrink=0.5, pad=0.02, fraction=0.02)
    cbar.set_label("CRUST1.0 Moho depth (km)", fontsize=10)

    ax.set_title("CRUST1.0 Moho Depth and Moho-Proximal Events (Δ≤10 km)",
                 fontsize=13, fontweight="bold", pad=10)

    plt.tight_layout()
    fig.savefig(FIG5_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 5 saved: {FIG5_PATH}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Generate all A3.B4 figures."""
    logger.info("Loading results JSON...")
    results = load_results()

    figure1_stratification_matrix(results)
    figure2_cramer_trends(results)
    figure3_moho_isolation(results)
    figure4_subduction_crosstab(results)
    figure5_moho_map(results)

    logger.info("All figures generated successfully.")


if __name__ == "__main__":
    main()
