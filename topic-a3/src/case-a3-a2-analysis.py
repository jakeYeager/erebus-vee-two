"""
Case A3.A2 — Stratified Schuster/MFPA Periodicity Audit

Main analysis script. Applies cluster-robust Schuster spectrum and MFPA
to stratified subsets of the ISC-GEM catalog to test four predictions from
prior Tier 2 cases (B2, B3, B4, A3).

Sub-tests:
  1. Baseline replication on full catalog with suppression characterization
  2. Signal-bearing (continental × mid_crustal) vs. non-signal-bearing
  3. NH/SH phase angle comparison for half-year period
  4. Full vs. G-K mainshocks vs. G-K aftershocks declustering sensitivity
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent

import importlib.util as _ilu

def _load_module(name: str, path: Path):
    spec = _ilu.spec_from_file_location(name, path)
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_schuster_mod = _load_module("case_a3_a2_schuster", BASE_DIR / "src" / "case-a3-a2-schuster.py")
_mfpa_mod = _load_module("case_a3_a2_mfpa", BASE_DIR / "src" / "case-a3-a2-mfpa.py")

schuster_spectrum = _schuster_mod.schuster_spectrum
mfpa_spectrum = _mfpa_mod.mfpa_spectrum
apply_fdr_bh = _mfpa_mod.apply_fdr_bh

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
JULIAN_YEAR_SECS = 31_557_600.0
PERIOD_SCAN_DAYS = np.logspace(np.log10(0.25), np.log10(548.0), 200)
EXPLICIT_PERIODS = {
    "quarter_year_756": 75.6,
    "quarter_year_912": 91.25,
    "half_year": 182.5,
    "annual": 365.25,
    "tidal_12h": 0.5,
    "tidal_24h": 1.0,
    "tidal_14d": 14.77,
    "tidal_27d": 27.32,
}
N_BOOTSTRAP = 10_000
CLUSTER_WINDOWS_DAYS = [1.0, 7.0]
RNG_SEED = 42

EPOCH = pd.Timestamp("1950-01-01 00:00:00", tz="UTC")

# A3.A3 permutation-significant suppressed bins (k=24 bins, phase centers)
A3A3_SUPPRESSED_BINS = [2, 8, 10, 11, 12, 13, 16, 18, 22]
A3A3_PERM_SIG_SUPPRESSED_PHASES = [0.104, 0.521, 0.563, 0.688]

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------
DATA_ROOT = BASE_DIR.parent / "data" / "iscgem"
FULL_CATALOG_PATH = DATA_ROOT / "iscgem_global_6-9_1950-2021.csv"
GSHHG_PATH = DATA_ROOT / "plate-location" / "ocean_class_gshhg_global.csv"
GK_MAINSHOCKS_PATH = DATA_ROOT / "declustering-algorithm" / "mainshocks_gk-seq_global.csv"
GK_AFTERSHOCKS_PATH = DATA_ROOT / "declustering-algorithm" / "aftershocks_gk-seq_global.csv"
OUTPUT_PATH = BASE_DIR / "output" / "case-a3-a2-results.json"


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def parse_event_time(df: pd.DataFrame) -> pd.DataFrame:
    """Parse event_at as UTC datetime and compute event_time_days since epoch."""
    df = df.copy()
    df["event_at"] = pd.to_datetime(df["event_at"], utc=True)
    df["event_time_days"] = (df["event_at"] - EPOCH).dt.total_seconds() / 86400.0
    return df


def add_tectonic_class(df: pd.DataFrame) -> pd.DataFrame:
    """Assign tectonic_class from dist_to_coast_km."""
    df = df.copy()
    df["tectonic_class"] = pd.cut(
        df["dist_to_coast_km"],
        bins=[-np.inf, 50.0, 200.0, np.inf],
        labels=["continental", "transitional", "oceanic"],
    )
    return df


def add_depth_band(df: pd.DataFrame) -> pd.DataFrame:
    """Assign depth_band from depth column."""
    df = df.copy()
    df["depth_band"] = pd.cut(
        df["depth"],
        bins=[-np.inf, 20.0, 70.0, 300.0, np.inf],
        labels=["shallow", "mid_crustal", "intermediate", "deep"],
    )
    return df


def explicit_period_results(
    schuster_results: list[dict],
    mfpa_results: list[dict],
    periods: dict,
) -> dict:
    """Extract p-values at the explicit test periods from scan arrays.

    Finds the nearest period in the scan array for each explicit period.
    """
    schuster_arr = np.array([r["period_days"] for r in schuster_results])
    mfpa_arr = np.array([r["period_days"] for r in mfpa_results])

    out = {}
    for key, period in periods.items():
        s_idx = int(np.argmin(np.abs(schuster_arr - period)))
        m_idx = int(np.argmin(np.abs(mfpa_arr - period)))
        sr = schuster_results[s_idx]
        mr = mfpa_results[m_idx]
        out[key] = {
            "period_days": period,
            "p_standard": sr["p_standard"],
            "p_cluster_robust": sr["p_cluster_robust"],
            "n_clusters": sr["n_clusters"],
            "dominant_phase_fraction": sr["dominant_phase_fraction"],
            "power": mr["power"],
            "p_mfpa": mr["p_mfpa"],
            "fdr_significant": mr.get("fdr_significant", False),
            "significant_95": mr["significant_95"],
            "significant_99": mr["significant_99"],
        }
    return out


def run_full_analysis(
    t: np.ndarray,
    dt_cluster_days: float = 1.0,
    label: str = "dataset",
) -> tuple[list[dict], list[dict]]:
    """Run Schuster and MFPA on an event-time array; apply FDR to MFPA."""
    log.info("Running Schuster (%s, dt_cluster=%.1fd, n=%d)", label, dt_cluster_days, len(t))
    s_results = schuster_spectrum(t, PERIOD_SCAN_DAYS, dt_cluster_days=dt_cluster_days, rng_seed=RNG_SEED)

    log.info("Running MFPA (%s, n=%d, n_bootstrap=%d)", label, len(t), N_BOOTSTRAP)
    m_results = mfpa_spectrum(t, PERIOD_SCAN_DAYS, n_bootstrap=N_BOOTSTRAP, rng_seed=RNG_SEED)

    # Apply FDR correction
    p_vals = np.array([r["p_mfpa"] for r in m_results])
    fdr_flags = apply_fdr_bh(p_vals)
    for i, r in enumerate(m_results):
        r["fdr_significant"] = bool(fdr_flags[i])

    return s_results, m_results


def inflation_factor(schuster_results: list[dict]) -> float:
    """Ratio of p_standard < 0.05 count to p_cluster_robust < 0.05 count."""
    n_standard = sum(1 for r in schuster_results if r["p_standard"] < 0.05)
    n_robust = sum(1 for r in schuster_results if r["p_cluster_robust"] < 0.05)
    if n_robust == 0:
        return float("inf")
    return float(n_standard) / float(n_robust)


def circular_distance(a: float, b: float) -> float:
    """Minimum circular distance between two phases in [0, 1)."""
    d = abs(a - b) % 1.0
    return float(min(d, 1.0 - d))


def nearest_suppressed_bin(anti_phase: float, perm_sig_phases: list) -> tuple[int, float]:
    """Find nearest A3.A3 perm-sig suppressed bin to the anti-dominant phase."""
    # Phase centers of A3A3_SUPPRESSED_BINS
    k = 24
    all_bin_centers = [(b + 0.5) / k for b in A3A3_SUPPRESSED_BINS]
    perm_centers = perm_sig_phases

    min_dist = float("inf")
    nearest_bin = -1
    for idx, (bin_idx, center) in enumerate(zip(A3A3_SUPPRESSED_BINS, all_bin_centers)):
        if center not in perm_centers:
            continue
        d = circular_distance(anti_phase, center)
        if d < min_dist:
            min_dist = d
            nearest_bin = bin_idx
    return nearest_bin, min_dist


def schuster_at_period(schuster_results: list[dict], target_period: float) -> dict:
    """Return the Schuster result nearest to target_period."""
    arr = np.array([r["period_days"] for r in schuster_results])
    idx = int(np.argmin(np.abs(arr - target_period)))
    return schuster_results[idx]


def mfpa_at_period(mfpa_results: list[dict], target_period: float) -> dict:
    """Return the MFPA result nearest to target_period."""
    arr = np.array([r["period_days"] for r in mfpa_results])
    idx = int(np.argmin(np.abs(arr - target_period)))
    return mfpa_results[idx]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # -----------------------------------------------------------------------
    # 1. Load and prepare data
    # -----------------------------------------------------------------------
    log.info("Loading full catalog from %s", FULL_CATALOG_PATH)
    df_full = pd.read_csv(FULL_CATALOG_PATH)
    assert len(df_full) == 9210, f"Expected 9,210 rows; got {len(df_full)}"
    df_full = parse_event_time(df_full)
    df_full["phase"] = (df_full["solar_secs"] / JULIAN_YEAR_SECS) % 1.0

    log.info("Loading GSHHG classification from %s", GSHHG_PATH)
    df_gshhg = pd.read_csv(GSHHG_PATH)
    assert len(df_gshhg) == 9210, f"Expected 9,210 GSHHG rows; got {len(df_gshhg)}"

    df_full = df_full.merge(df_gshhg, on="usgs_id", how="left")
    assert df_full["ocean_class"].notna().all(), "NaN values in ocean_class after merge"

    df_full = add_tectonic_class(df_full)
    df_full = add_depth_band(df_full)
    assert df_full["tectonic_class"].notna().all(), "NaN in tectonic_class"

    df_full["hemisphere"] = np.where(df_full["latitude"] >= 0, "NH", "SH")

    # Define strata
    df_signal = df_full[
        (df_full["tectonic_class"] == "continental") & (df_full["depth_band"] == "mid_crustal")
    ]
    df_non_signal = df_full[~df_full.index.isin(df_signal.index)]
    df_nh = df_full[df_full["hemisphere"] == "NH"]
    df_sh = df_full[df_full["hemisphere"] == "SH"]

    # G-K catalogs
    log.info("Loading G-K mainshocks from %s", GK_MAINSHOCKS_PATH)
    df_gk_ms = pd.read_csv(GK_MAINSHOCKS_PATH)
    assert len(df_gk_ms) == 5883, f"Expected 5,883 G-K mainshocks; got {len(df_gk_ms)}"
    df_gk_ms = parse_event_time(df_gk_ms)
    df_gk_ms["phase"] = (df_gk_ms["solar_secs"] / JULIAN_YEAR_SECS) % 1.0

    log.info("Loading G-K aftershocks from %s", GK_AFTERSHOCKS_PATH)
    df_gk_as = pd.read_csv(GK_AFTERSHOCKS_PATH)
    assert len(df_gk_as) == 3327, f"Expected 3,327 G-K aftershocks; got {len(df_gk_as)}"
    df_gk_as = parse_event_time(df_gk_as)
    df_gk_as["phase"] = (df_gk_as["solar_secs"] / JULIAN_YEAR_SECS) % 1.0

    # Log stratum sizes
    n_signal = len(df_signal)
    n_non_signal = len(df_non_signal)
    n_nh = len(df_nh)
    n_sh = len(df_sh)
    n_gk_ms = len(df_gk_ms)
    n_gk_as = len(df_gk_as)

    log.info(
        "Strata: signal_bearing=%d, non_signal=%d, NH=%d, SH=%d, GK_ms=%d, GK_as=%d",
        n_signal, n_non_signal, n_nh, n_sh, n_gk_ms, n_gk_as,
    )

    # -----------------------------------------------------------------------
    # 2. Sub-test 1 — Baseline full catalog
    # -----------------------------------------------------------------------
    log.info("Sub-test 1: baseline full catalog")
    t_full = df_full["event_time_days"].values

    s1_1day, m1_1day = run_full_analysis(t_full, dt_cluster_days=1.0, label="full_1day")
    s1_7day, _ = run_full_analysis(t_full, dt_cluster_days=7.0, label="full_7day")

    # 7-day MFPA (reuse same bootstrap; only need schuster for 7-day)
    # For 7-day, run MFPA separately
    log.info("Running MFPA 7-day baseline (schuster only needed different dt)")
    m1_7day = mfpa_spectrum(t_full, PERIOD_SCAN_DAYS, n_bootstrap=N_BOOTSTRAP, rng_seed=RNG_SEED)
    p_vals_7 = np.array([r["p_mfpa"] for r in m1_7day])
    fdr_7 = apply_fdr_bh(p_vals_7)
    for i, r in enumerate(m1_7day):
        r["fdr_significant"] = bool(fdr_7[i])

    inf_1day = inflation_factor(s1_1day)
    inf_7day = inflation_factor(s1_7day)
    log.info("Inflation factor 1-day=%.1f, 7-day=%.1f (A2.A1 reference: 329x)", inf_1day, inf_7day)

    # Explicit period tests (sub-test 1)
    exp1 = explicit_period_results(s1_1day, m1_1day, EXPLICIT_PERIODS)

    # Top-10 MFPA by power
    sorted_mfpa = sorted(m1_1day, key=lambda r: r["power"], reverse=True)
    top10 = []
    for r in sorted_mfpa[:10]:
        s_match = schuster_at_period(s1_1day, r["period_days"])
        top10.append({
            "period_days": r["period_days"],
            "power": r["power"],
            "p_mfpa": r["p_mfpa"],
            "fdr_significant": r["fdr_significant"],
            "dominant_phase_fraction": s_match["dominant_phase_fraction"],
        })

    # Suppression characterization
    # Half-year dominant and anti-dominant phase
    hy_schuster = schuster_at_period(s1_1day, 182.5)
    hy_dominant = hy_schuster["dominant_phase_fraction"]
    hy_anti = (hy_dominant + 0.5) % 1.0

    nearest_bin, align_dist = nearest_suppressed_bin(hy_anti, A3A3_PERM_SIG_SUPPRESSED_PHASES)
    consistent = bool(align_dist < 1.0 / 48.0)  # half bin width at k=24

    # Naive solstice comparison
    june_solstice_phase = 0.46
    naive_dist = circular_distance(hy_anti, june_solstice_phase)

    # Quarter-year vs half-year power ratio
    qy_mfpa = mfpa_at_period(m1_1day, 75.6)
    qy912_mfpa = mfpa_at_period(m1_1day, 91.25)
    hy_mfpa = mfpa_at_period(m1_1day, 182.5)
    # Use the stronger quarter-year
    qy_power = max(qy_mfpa["power"], qy912_mfpa["power"])
    hy_power = hy_mfpa["power"]
    power_ratio = float(hy_power / qy_power) if qy_power > 0 else float("inf")
    oscillation_symmetric = bool(power_ratio >= 1.0)

    # FDR-significant periods suppression table
    fdr_sig_suppression = []
    for r in m1_1day:
        if r["fdr_significant"]:
            s_match = schuster_at_period(s1_1day, r["period_days"])
            dp = s_match["dominant_phase_fraction"]
            anti = (dp + 0.5) % 1.0
            nb, nd = nearest_suppressed_bin(anti, A3A3_PERM_SIG_SUPPRESSED_PHASES)
            fdr_sig_suppression.append({
                "period_days": r["period_days"],
                "dominant_phase_fraction": dp,
                "anti_dominant_phase": anti,
                "nearest_a3a3_suppressed_bin": nb,
                "alignment_distance": nd,
            })

    suppression_char = {
        "half_year_dominant_phase": float(hy_dominant),
        "half_year_anti_dominant_phase": float(hy_anti),
        "half_year_anti_phase_nearest_a3a3_suppressed_bin": int(nearest_bin),
        "half_year_anti_phase_alignment_distance": float(align_dist),
        "half_year_anti_phase_consistent_with_chi2_suppression": consistent,
        "half_year_naive_solstice_phase_june": june_solstice_phase,
        "half_year_naive_solstice_distance": float(naive_dist),
        "quarter_year_vs_half_year_power_ratio": power_ratio,
        "oscillation_symmetric": oscillation_symmetric,
        "fdr_sig_periods_suppression": fdr_sig_suppression,
    }

    subtest1 = {
        "schuster_1day": s1_1day,
        "schuster_7day": s1_7day,
        "mfpa": m1_1day,
        "explicit_period_tests": exp1,
        "inflation_factor_1day": inf_1day,
        "inflation_factor_7day": inf_7day,
        "top10_mfpa_periods": top10,
        "suppression_characterization": suppression_char,
    }

    # -----------------------------------------------------------------------
    # 3. Sub-test 2 — Signal-bearing vs. non-signal-bearing
    # -----------------------------------------------------------------------
    log.info("Sub-test 2: signal-bearing (n=%d) vs non-signal (n=%d)", n_signal, n_non_signal)
    t_signal = df_signal["event_time_days"].values
    t_non_signal = df_non_signal["event_time_days"].values

    s2_sig, m2_sig = run_full_analysis(t_signal, dt_cluster_days=1.0, label="signal_bearing")
    s2_non, m2_non = run_full_analysis(t_non_signal, dt_cluster_days=1.0, label="non_signal")

    exp2_sig = explicit_period_results(s2_sig, m2_sig, EXPLICIT_PERIODS)
    exp2_non = explicit_period_results(s2_non, m2_non, EXPLICIT_PERIODS)

    subtest2 = {
        "signal_bearing": {
            "n": int(n_signal),
            "schuster_1day": s2_sig,
            "mfpa": m2_sig,
            "explicit_period_tests": exp2_sig,
        },
        "non_signal": {
            "n": int(n_non_signal),
            "schuster_1day": s2_non,
            "mfpa": m2_non,
            "explicit_period_tests": exp2_non,
        },
    }

    # -----------------------------------------------------------------------
    # 4. Sub-test 3 — NH/SH phase angle comparison
    # -----------------------------------------------------------------------
    log.info("Sub-test 3: NH (n=%d) vs SH (n=%d) phase angles", n_nh, n_sh)
    t_nh = df_nh["event_time_days"].values
    t_sh = df_sh["event_time_days"].values

    s3_nh, m3_nh = run_full_analysis(t_nh, dt_cluster_days=1.0, label="NH")
    s3_sh, m3_sh = run_full_analysis(t_sh, dt_cluster_days=1.0, label="SH")

    exp3_nh = explicit_period_results(s3_nh, m3_nh, EXPLICIT_PERIODS)
    exp3_sh = explicit_period_results(s3_sh, m3_sh, EXPLICIT_PERIODS)

    # Half-year dominant phase
    nh_hy = schuster_at_period(s3_nh, 182.5)["dominant_phase_fraction"]
    sh_hy = schuster_at_period(s3_sh, 182.5)["dominant_phase_fraction"]
    nh_qy = schuster_at_period(s3_nh, 91.25)["dominant_phase_fraction"]
    sh_qy = schuster_at_period(s3_sh, 91.25)["dominant_phase_fraction"]

    # Wrapped delta: NH − SH, wrapped to [-0.5, 0.5]
    delta_hy = (nh_hy - sh_hy + 0.5) % 1.0 - 0.5
    delta_hy_days = float(delta_hy * 182.5)

    # Reference comparisons
    ref_antiphase_days = 91.25
    ref_b2_days = 36.0
    closer_to_antiphase = bool(abs(delta_hy_days) > abs(delta_hy_days - ref_antiphase_days))

    subtest3 = {
        "nh": {
            "n": int(n_nh),
            "explicit_period_tests": exp3_nh,
            "half_year_dominant_phase": float(nh_hy),
            "quarter_year_912_dominant_phase": float(nh_qy),
        },
        "sh": {
            "n": int(n_sh),
            "explicit_period_tests": exp3_sh,
            "half_year_dominant_phase": float(sh_hy),
            "quarter_year_912_dominant_phase": float(sh_qy),
        },
        "half_year_phase_delta_fraction": float(delta_hy),
        "half_year_phase_delta_days": delta_hy_days,
        "reference_antiphase_days": ref_antiphase_days,
        "reference_b2_predicted_delta_days": ref_b2_days,
        "closer_to_antiphase": closer_to_antiphase,
    }

    # -----------------------------------------------------------------------
    # 5. Sub-test 4 — Declustering sensitivity
    # -----------------------------------------------------------------------
    log.info("Sub-test 4: declustering sensitivity (full / GK-ms / GK-as)")

    t_gk_ms = df_gk_ms["event_time_days"].values
    t_gk_as = df_gk_as["event_time_days"].values

    s4_full, m4_full = run_full_analysis(t_full, dt_cluster_days=1.0, label="ST4_full")
    s4_gk_ms, m4_gk_ms = run_full_analysis(t_gk_ms, dt_cluster_days=1.0, label="GK_mainshocks")
    s4_gk_as, m4_gk_as = run_full_analysis(t_gk_as, dt_cluster_days=1.0, label="GK_aftershocks")

    exp4_full = explicit_period_results(s4_full, m4_full, EXPLICIT_PERIODS)
    exp4_gk_ms = explicit_period_results(s4_gk_ms, m4_gk_ms, EXPLICIT_PERIODS)
    exp4_gk_as = explicit_period_results(s4_gk_as, m4_gk_as, EXPLICIT_PERIODS)

    # Interval 1 contradiction resolution:
    # If 75.6d period is detected (p_cluster_robust < 0.05) in full catalog
    # but NOT in GK mainshocks, contradiction is explained (aftershock-driven).
    # If it survives in GK mainshocks, contradiction deepens.
    full_qy_sig = exp4_full["quarter_year_756"]["p_cluster_robust"] < 0.05
    ms_qy_sig = exp4_gk_ms["quarter_year_756"]["p_cluster_robust"] < 0.05
    as_qy_sig = exp4_gk_as["quarter_year_756"]["p_cluster_robust"] < 0.05

    interval1_survives_gk_mainshocks = bool(ms_qy_sig)
    # Contradiction resolved = detected in full catalog and aftershocks, NOT in mainshocks
    interval1_contradiction_resolved = bool(
        full_qy_sig and as_qy_sig and not ms_qy_sig
    )

    subtest4 = {
        "full": {
            "n": 9210,
            "schuster_1day": s4_full,
            "mfpa": m4_full,
            "explicit_period_tests": exp4_full,
        },
        "gk_mainshocks": {
            "n": int(n_gk_ms),
            "schuster_1day": s4_gk_ms,
            "mfpa": m4_gk_ms,
            "explicit_period_tests": exp4_gk_ms,
        },
        "gk_aftershocks": {
            "n": int(n_gk_as),
            "schuster_1day": s4_gk_as,
            "mfpa": m4_gk_as,
            "explicit_period_tests": exp4_gk_as,
        },
        "interval1_contradiction_resolved": interval1_contradiction_resolved,
        "interval1_survives_gk_mainshocks": interval1_survives_gk_mainshocks,
    }

    # -----------------------------------------------------------------------
    # 6. Summary
    # -----------------------------------------------------------------------
    # Quarter-year vs half-year: which is dominant?
    qy_top_power = max(
        exp1["quarter_year_756"]["power"],
        exp1["quarter_year_912"]["power"],
    )
    hy_top_power = exp1["half_year"]["power"]
    qy_vs_hy_dominant = "quarter_year" if qy_top_power >= hy_top_power else "half_year"

    half_year_detected_full = bool(
        exp1["half_year"]["p_cluster_robust"] < 0.05
        or exp1["half_year"]["fdr_significant"]
    )

    # Signal-bearing stronger: compare quarter-year power
    sig_qy_power = max(
        exp2_sig["quarter_year_756"]["power"],
        exp2_sig["quarter_year_912"]["power"],
    )
    non_qy_power = max(
        exp2_non["quarter_year_756"]["power"],
        exp2_non["quarter_year_912"]["power"],
    )
    quarter_year_signal_bearing_stronger = bool(sig_qy_power > non_qy_power)

    summary = {
        "quarter_year_signal_bearing_stronger": quarter_year_signal_bearing_stronger,
        "half_year_detected_full_catalog": half_year_detected_full,
        "quarter_year_vs_half_year_dominant": qy_vs_hy_dominant,
        "quarter_year_vs_half_year_power_ratio": power_ratio,
        "oscillation_symmetric": oscillation_symmetric,
        "half_year_trough_consistent_with_chi2_suppression": consistent,
        "nh_sh_phase_delta_days": delta_hy_days,
        "nh_sh_closer_to_antiphase": closer_to_antiphase,
        "interval1_contradiction_resolved": interval1_contradiction_resolved,
        "inflation_factor_1day": inf_1day,
    }

    # -----------------------------------------------------------------------
    # 7. Assemble and write results JSON
    # -----------------------------------------------------------------------
    results = {
        "case": "A3.A2",
        "title": "Stratified Schuster/MFPA Periodicity Audit",
        "parameters": {
            "n_catalog": 9210,
            "julian_year_secs": JULIAN_YEAR_SECS,
            "n_periods_scan": int(len(PERIOD_SCAN_DAYS)),
            "period_range_days": [float(PERIOD_SCAN_DAYS[0]), float(PERIOD_SCAN_DAYS[-1])],
            "explicit_periods": {k: v for k, v in EXPLICIT_PERIODS.items()},
            "n_bootstrap": N_BOOTSTRAP,
            "cluster_windows_days": CLUSTER_WINDOWS_DAYS,
            "rng_seed": RNG_SEED,
        },
        "strata_sizes": {
            "n_full": 9210,
            "n_signal_bearing": int(n_signal),
            "n_non_signal": int(n_non_signal),
            "n_nh": int(n_nh),
            "n_sh": int(n_sh),
            "n_gk_mainshocks": int(n_gk_ms),
            "n_gk_aftershocks": int(n_gk_as),
        },
        "subtest_1_baseline": subtest1,
        "subtest_2_stratified": subtest2,
        "subtest_3_nh_sh": subtest3,
        "subtest_4_declustering": subtest4,
        "summary": summary,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(results, f, indent=2, default=str)

    log.info("Results written to %s", OUTPUT_PATH)

    # Print summary to stdout
    log.info("=== SUMMARY ===")
    log.info("Signal-bearing n=%d, Non-signal n=%d", n_signal, n_non_signal)
    log.info("Inflation factor 1-day: %.1f (A2.A1 reference: 329×)", inf_1day)
    log.info("Inflation factor 7-day: %.1f", inf_7day)
    log.info("Quarter-year vs half-year dominant: %s", qy_vs_hy_dominant)
    log.info("Power ratio (half/quarter): %.4f", power_ratio)
    log.info("Oscillation symmetric: %s", oscillation_symmetric)
    log.info("Half-year trough consistent with A3.A3 suppression: %s", consistent)
    log.info("NH/SH phase delta: %.1f days (antiphase=91.25d, B2=36d)", delta_hy_days)
    log.info("Closer to antiphase: %s", closer_to_antiphase)
    log.info("Interval 1 contradiction resolved: %s", interval1_contradiction_resolved)


if __name__ == "__main__":
    main()
