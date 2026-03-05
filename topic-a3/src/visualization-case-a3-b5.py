"""
Case A3.B5: Visualization Script — Corrected Null-Distribution Geometric Variable Test

Generates 5 figures:
  1. Null distribution correction (corrected vs uniform for non-cyclic variables)
  2. Corrected bin distributions (all 4 variables, k=24 full catalog)
  3. Variable ranking by stratum (grouped bar chart)
  4. Correction impact scatter (chi2_uniform vs chi2_corrected)
  5. A1b interval alignment with solar variables
"""

import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_PATH = BASE_DIR / "output" / "case-a3-b5-results.json"
OUTPUT_DIR = BASE_DIR / "output"

# ---------------------------------------------------------------------------
# Load results
# ---------------------------------------------------------------------------

with open(RESULTS_PATH) as fh:
    results = json.load(fh)

strata = results["strata"]
full_k24 = results["full_multibin"]["k24"]
var_ranking = results["variable_ranking"]
a1b_alignment = results["a1b_alignment"]
null_gen = results["null_generation"]
var_ranges = results["parameters"]["var_ranges"]

VARIABLES = ["solar_phase", "solar_declination", "declination_rate", "earth_sun_distance"]
VAR_LABELS = {
    "solar_phase": "Solar Phase [0, 1)",
    "solar_declination": "Solar Declination (°)",
    "declination_rate": "Declination Rate (°/day)",
    "earth_sun_distance": "Earth-Sun Distance (AU)",
}
VAR_UNITS = {
    "solar_phase": "fraction",
    "solar_declination": "degrees",
    "declination_rate": "°/day",
    "earth_sun_distance": "AU",
}

STRATA_LABELS = {
    "full": "Full Catalog",
    "continental": "Continental",
    "midcrustal": "Mid-Crustal (20–70 km)",
    "continental_midcrustal": "Continental × Mid-Crustal",
}


# ---------------------------------------------------------------------------
# Helper: regenerate analytic null for plotting
# ---------------------------------------------------------------------------

def regenerate_null(k: int, var_ranges_dict: dict) -> dict:
    """Regenerate the analytic null distribution for plotting."""
    n_days = int((2021 - 1950 + 1) * 365.25)
    days = np.arange(0, n_days, dtype=float)
    D = days + (1950 - 2000) * 365.25 - 0.5

    L = np.radians((280.46 + 0.9856474 * D) % 360)
    g = np.radians((357.528 + 0.9856003 * D) % 360)
    lam = L + np.radians(1.915 * np.sin(g) + 0.020 * np.sin(2 * g))

    dec = np.degrees(np.arcsin(np.sin(np.radians(23.439)) * np.sin(lam)))
    dist = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    rate = np.gradient(dec, days)

    result = {}
    for var_name, synthetic_values in [
        ("solar_declination", dec),
        ("declination_rate", rate),
        ("earth_sun_distance", dist),
    ]:
        vmin, vmax = var_ranges_dict[var_name]
        edges = np.linspace(vmin, vmax, k + 1)
        null_counts, _ = np.histogram(synthetic_values, bins=edges)
        total = null_counts.sum()
        result[var_name] = null_counts / total
    return result


# ---------------------------------------------------------------------------
# Figure 1: Null distribution correction
# ---------------------------------------------------------------------------
logger.info("Generating Figure 1: Null distribution correction...")

null_k24 = regenerate_null(24, {
    var: (var_ranges[var][0], var_ranges[var][1])
    for var in ["solar_declination", "declination_rate", "earth_sun_distance"]
})

fig1, axes1 = plt.subplots(1, 3, figsize=(14, 5))
fig1.suptitle("Time-Weighted Null Distributions — Corrected vs. Uniform (A3.B5)", fontsize=13)

non_cyclic_vars = ["solar_declination", "declination_rate", "earth_sun_distance"]

for ax, var in zip(axes1, non_cyclic_vars):
    vmin, vmax = var_ranges[var]
    k = 24
    edges = np.linspace(vmin, vmax, k + 1)
    bin_centers = (edges[:-1] + edges[1:]) / 2.0
    width = (vmax - vmin) / k

    corr_null = null_k24[var]
    uniform_null = np.full(k, 1.0 / k)

    # Plot steelblue bars for corrected null
    ax.bar(bin_centers, corr_null, width=width * 0.85, color="steelblue",
           alpha=0.7, label="Corrected null")

    # Gray dashed line for uniform null
    ax.axhline(1.0 / k, color="gray", linestyle="--", linewidth=1.5, label="Uniform null")

    # Shade deviations: green where corrected > uniform, red where corrected < uniform
    for i in range(k):
        left = edges[i]
        right = edges[i + 1]
        if corr_null[i] > uniform_null[i]:
            ax.fill_betweenx(
                [uniform_null[i], corr_null[i]], left, right,
                color="green", alpha=0.3
            )
        else:
            ax.fill_betweenx(
                [corr_null[i], uniform_null[i]], left, right,
                color="red", alpha=0.3
            )

    max_dev_pct = float(np.abs(corr_null - uniform_null).max() / (1.0 / k) * 100)
    ax.set_title(VAR_LABELS[var], fontsize=10)
    ax.set_xlabel(VAR_UNITS[var])
    ax.set_ylabel("Fraction of time" if ax == axes1[0] else "")
    ax.legend(fontsize=8)
    ax.annotate(
        f"n_synth={null_gen['n_synthetic_points']:,}\nMax dev={max_dev_pct:.1f}%",
        xy=(0.05, 0.93), xycoords="axes fraction", fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)
    )

fig1.tight_layout()
out1 = OUTPUT_DIR / "case-a3-b5-null-distributions.png"
fig1.savefig(out1, dpi=300, bbox_inches="tight")
plt.close(fig1)
logger.info(f"Figure 1 saved to {out1}")


# ---------------------------------------------------------------------------
# Figure 2: Corrected bin distributions
# ---------------------------------------------------------------------------
logger.info("Generating Figure 2: Corrected bin distributions...")

fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10))
fig2.suptitle("Solar Variable Distributions vs. Corrected Null (Full Catalog, k=24)", fontsize=13)

var_positions = [
    ("solar_phase", axes2[0, 0]),
    ("solar_declination", axes2[0, 1]),
    ("declination_rate", axes2[1, 0]),
    ("earth_sun_distance", axes2[1, 1]),
]

# A1b intervals for solar_phase panel
A1B_BINS = {"interval_1": [4, 5], "interval_2": [15], "interval_3": [21]}
EQUINOX_PHASES = [0.0, 0.25, 0.50, 0.75, 1.0]  # spring=0.25, autumn=0.75
SOLSTICE_PHASES = [0.125, 0.625]                 # summer=0.125, winter=0.625

# A1b expected physical values for vertical lines
A1B_EXPECTED = {
    "solar_declination": {
        "interval_1": results["a1b_alignment"]["interval_1"]["solar_declination_expected"],
        "interval_2": results["a1b_alignment"]["interval_2"]["solar_declination_expected"],
        "interval_3": results["a1b_alignment"]["interval_3"]["solar_declination_expected"],
    },
    "declination_rate": {
        "interval_1": results["a1b_alignment"]["interval_1"]["declination_rate_expected"],
        "interval_2": results["a1b_alignment"]["interval_2"]["declination_rate_expected"],
        "interval_3": results["a1b_alignment"]["interval_3"]["declination_rate_expected"],
    },
    "earth_sun_distance": {
        "interval_1": results["a1b_alignment"]["interval_1"]["earth_sun_distance_expected"],
        "interval_2": results["a1b_alignment"]["interval_2"]["earth_sun_distance_expected"],
        "interval_3": results["a1b_alignment"]["interval_3"]["earth_sun_distance_expected"],
    },
}

for var, ax in var_positions:
    stats = full_k24[var]
    k = stats["k"]
    observed = np.array(stats["bin_counts"])
    exp_corr = np.array(stats["expected_corrected"])
    exp_unif = np.array(stats["expected_uniform"])

    if var == "solar_phase":
        bin_centers = (np.arange(k) + 0.5) / k
        width = 1.0 / k
        x_label = "Solar Phase [0, 1)"
    else:
        vmin, vmax = var_ranges[var]
        edges = np.linspace(vmin, vmax, k + 1)
        bin_centers = (edges[:-1] + edges[1:]) / 2.0
        width = (vmax - vmin) / k
        x_label = VAR_LABELS[var]

    # Steelblue bars
    ax.bar(bin_centers, observed, width=width * 0.85, color="steelblue", alpha=0.7,
           label="Observed")

    # Corrected expected: solid line
    ax.plot(bin_centers, exp_corr, "k-", linewidth=1.5, label="Corrected expected")

    # Uniform expected: dashed line
    ax.plot(bin_centers, exp_unif, "r--", linewidth=1.2, label="Uniform expected")

    # 1-SD band around corrected expected
    sd_band = np.sqrt(np.maximum(exp_corr, 1.0))
    ax.fill_between(bin_centers, exp_corr - sd_band, exp_corr + sd_band,
                    color="gray", alpha=0.2)

    if var == "solar_phase":
        # A1b interval shading
        ax.axvspan((4.0) / k, (6.0) / k, alpha=0.12, color="gray", label="A1b intervals")
        ax.axvspan(15.0 / k, 16.0 / k, alpha=0.12, color="gray")
        ax.axvspan(21.0 / k, 22.0 / k, alpha=0.12, color="gray")
        # Equinox markers
        for ep in [0.25, 0.75]:
            ax.axvline(ep, color="green", linestyle=":", linewidth=1.0, alpha=0.7)
        # Solstice markers
        for sp in [0.125, 0.625]:
            ax.axvline(sp, color="orange", linestyle=":", linewidth=1.0, alpha=0.7)
        ax.set_xlim(0, 1)
    else:
        # Vertical dotted red lines at A1b interval expected values
        for ivl_key in ["interval_1", "interval_2", "interval_3"]:
            xval = A1B_EXPECTED[var][ivl_key]
            ax.axvline(xval, color="red", linestyle=":", linewidth=1.0, alpha=0.7)

    # Annotations
    p_corr = stats["p_corrected"]
    v_corr = stats["cramers_v_corrected"]
    chi2_corr = stats["chi2_corrected"]
    chi2_unif = stats["chi2_uniform"]
    p_unif = stats["p_uniform"]

    p_str = f"{p_corr:.2e}" if p_corr < 0.001 else f"{p_corr:.4f}"
    annot = (f"χ²_corr={chi2_corr:.2f}, p={p_str}\n"
             f"V_corr={v_corr:.4f}")
    ax.annotate(annot, xy=(0.98, 0.97), xycoords="axes fraction",
                ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Smaller text for uniform comparison
    p_unif_str = f"{p_unif:.2e}" if p_unif < 0.001 else f"{p_unif:.4f}"
    annot2 = f"χ²_unif={chi2_unif:.2f}, p_unif={p_unif_str}"
    ax.annotate(annot2, xy=(0.98, 0.78), xycoords="axes fraction",
                ha="right", va="top", fontsize=7, color="gray",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7))

    ax.set_xlabel(x_label)
    ax.set_ylabel("Count")
    ax.set_title(VAR_LABELS[var])
    ax.legend(fontsize=7)

fig2.tight_layout()
out2 = OUTPUT_DIR / "case-a3-b5-binplots.png"
fig2.savefig(out2, dpi=300, bbox_inches="tight")
plt.close(fig2)
logger.info(f"Figure 2 saved to {out2}")


# ---------------------------------------------------------------------------
# Figure 3: Variable ranking by stratum
# ---------------------------------------------------------------------------
logger.info("Generating Figure 3: Variable ranking by stratum...")

fig3, axes3 = plt.subplots(2, 2, figsize=(13, 9))
fig3.suptitle("Corrected Cramér's V by Variable and Stratum (A3.B5)", fontsize=13)

stratum_order = ["full", "continental", "midcrustal", "continental_midcrustal"]
ax_positions = [axes3[0, 0], axes3[0, 1], axes3[1, 0], axes3[1, 1]]

var_short_labels = {
    "solar_phase": "solar_phase",
    "solar_declination": "declination",
    "declination_rate": "decl_rate",
    "earth_sun_distance": "e-s dist",
}

for stratum_name, ax in zip(stratum_order, ax_positions):
    stratum_data = strata[stratum_name]
    if stratum_data is None:
        ax.set_title(f"{STRATA_LABELS[stratum_name]} (skipped — n<100)")
        ax.axis("off")
        continue

    n = stratum_data["n"]
    y_positions = np.arange(len(VARIABLES))
    bar_height = 0.35

    corr_vals = []
    unif_vals = []
    sig_markers = []

    for var in VARIABLES:
        vs = stratum_data[var]
        corr_vals.append(vs["cramers_v_corrected"])
        unif_vals.append(vs["cramers_v_uniform"])
        p = vs["p_corrected"]
        if p < 0.01:
            sig_markers.append("**")
        elif p < 0.05:
            sig_markers.append("*")
        else:
            sig_markers.append("ns")

    # Steelblue = corrected, light gray = uniform
    bars_corr = ax.barh(y_positions + bar_height / 2, corr_vals, height=bar_height,
                        color="steelblue", label="Corrected V")
    bars_unif = ax.barh(y_positions - bar_height / 2, unif_vals, height=bar_height,
                        color="lightgray", label="Uniform V")

    # Significance markers on corrected bars
    for i, (val, marker) in enumerate(zip(corr_vals, sig_markers)):
        ax.text(val + 0.0005, i + bar_height / 2, marker,
                va="center", ha="left", fontsize=9, color="black")

    ax.set_yticks(y_positions)
    ax.set_yticklabels([var_short_labels[v] for v in VARIABLES])
    ax.set_xlabel("Cramér's V")
    ax.set_title(f"{STRATA_LABELS[stratum_name]}\n(n={n})", fontsize=10)
    ax.legend(fontsize=8)
    ax.set_xlim(0, max(max(unif_vals) * 1.15, 0.01))

fig3.tight_layout()
out3 = OUTPUT_DIR / "case-a3-b5-variable-ranking.png"
fig3.savefig(out3, dpi=300, bbox_inches="tight")
plt.close(fig3)
logger.info(f"Figure 3 saved to {out3}")


# ---------------------------------------------------------------------------
# Figure 4: Correction impact scatter
# ---------------------------------------------------------------------------
logger.info("Generating Figure 4: Correction impact scatter...")

fig4, ax4 = plt.subplots(figsize=(8, 7))
fig4.suptitle("Chi-Square Before vs. After Null Correction (A3.B5)", fontsize=13)

var_colors = {
    "solar_declination": "steelblue",
    "declination_rate": "darkorange",
    "earth_sun_distance": "green",
    "solar_phase": "red",
}
stratum_markers = {
    "full": "o",
    "continental": "s",
    "midcrustal": "^",
    "continental_midcrustal": "D",
}

all_chi2 = []

for stratum_name, stratum_data in strata.items():
    if stratum_data is None:
        continue
    for var in VARIABLES:
        vs = stratum_data[var]
        chi2_u = vs["chi2_uniform"]
        chi2_c = vs["chi2_corrected"]
        all_chi2.extend([chi2_u, chi2_c])

        ax4.scatter(
            chi2_u, chi2_c,
            color=var_colors[var],
            marker=stratum_markers[stratum_name],
            s=70, alpha=0.8, zorder=3
        )

        # Annotate large deviations
        if chi2_u > 0 and abs(chi2_c - chi2_u) / chi2_u > 0.20:
            label = f"{var[:4]}/{stratum_name[:3]}"
            ax4.annotate(
                label, (chi2_u, chi2_c),
                fontsize=7, xytext=(3, 3), textcoords="offset points"
            )

# Identity line (y=x)
if all_chi2:
    max_val = max(all_chi2) * 1.05
    ax4.plot([0, max_val], [0, max_val], "gray", linestyle="-", linewidth=1.2,
             label="y = x (no change)", zorder=1)
    ax4.set_xlim(0, max_val)
    ax4.set_ylim(0, max_val)

ax4.set_xlabel("χ² Uniform (uncorrected)")
ax4.set_ylabel("χ² Corrected")

# Legend for variables
var_patches = [mpatches.Patch(color=c, label=v) for v, c in var_colors.items()]
# Legend for strata
from matplotlib.lines import Line2D
strata_handles = [
    Line2D([0], [0], marker=m, color="gray", linestyle="None", markersize=8, label=s)
    for s, m in stratum_markers.items()
]
leg1 = ax4.legend(handles=var_patches, loc="upper left", fontsize=8, title="Variable")
ax4.add_artist(leg1)
ax4.legend(handles=strata_handles + [
    Line2D([0], [0], color="gray", linestyle="-", label="y = x")
], loc="lower right", fontsize=8, title="Stratum")

ax4.grid(True, alpha=0.3)

fig4.tight_layout()
out4 = OUTPUT_DIR / "case-a3-b5-correction-delta.png"
fig4.savefig(out4, dpi=300, bbox_inches="tight")
plt.close(fig4)
logger.info(f"Figure 4 saved to {out4}")


# ---------------------------------------------------------------------------
# Figure 5: A1b interval alignment
# ---------------------------------------------------------------------------
logger.info("Generating Figure 5: A1b interval alignment...")

fig5, axes5 = plt.subplots(3, 3, figsize=(13, 10))
fig5.suptitle("A1b Interval Alignment with Solar Variables — Corrected Null (A3.B5)",
              fontsize=13)

intervals = [
    ("interval_1", "Interval 1 (Phase ~0.22, ~Mar 22)"),
    ("interval_2", "Interval 2 (Phase ~0.64, ~Aug 22)"),
    ("interval_3", "Interval 3 (Phase ~0.90, ~Nov 24)"),
]
non_cyclic = ["solar_declination", "declination_rate", "earth_sun_distance"]

for row_idx, (ivl_key, ivl_label) in enumerate(intervals):
    ivl_data = a1b_alignment[ivl_key]

    for col_idx, var in enumerate(non_cyclic):
        ax = axes5[row_idx, col_idx]

        stats = full_k24[var]
        k = stats["k"]
        observed = np.array(stats["bin_counts"])
        exp_corr = np.array(stats["expected_corrected"])

        vmin, vmax = var_ranges[var]
        edges = np.linspace(vmin, vmax, k + 1)
        bin_centers = (edges[:-1] + edges[1:]) / 2.0
        width = (vmax - vmin) / k

        # Normalized ratio: observed / expected_corrected
        ratio = np.where(exp_corr > 0, observed / exp_corr, 1.0)

        ax.bar(bin_centers, ratio, width=width * 0.85, color="steelblue", alpha=0.7)

        # Reference line at 1.0
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)

        # Vertical red dotted line at expected physical value
        expected_val = ivl_data[f"{var}_expected"]
        ax.axvline(expected_val, color="red", linestyle=":", linewidth=1.5)

        # Check if elevated
        bin_idx = ivl_data[f"{var}_bin_idx"]
        is_elevated = ivl_data[f"{var}_elevated"]
        elev_text = "ELEVATED" if is_elevated else "not elevated"
        elev_color = "green" if is_elevated else "gray"

        ax.annotate(
            elev_text, xy=(0.97, 0.97), xycoords="axes fraction",
            ha="right", va="top", fontsize=8, color=elev_color,
            fontweight="bold" if is_elevated else "normal"
        )

        if col_idx == 0:
            ax.set_ylabel(ivl_label.replace("(", "\n("), fontsize=8)
        if row_idx == 0:
            ax.set_title(f"{VAR_LABELS[var]}", fontsize=9)
        if row_idx == 2:
            ax.set_xlabel(VAR_UNITS[var])

        ax.set_xlim(vmin, vmax)

fig5.tight_layout()
out5 = OUTPUT_DIR / "case-a3-b5-a1b-alignment.png"
fig5.savefig(out5, dpi=300, bbox_inches="tight")
plt.close(fig5)
logger.info(f"Figure 5 saved to {out5}")

logger.info("All figures generated successfully.")
