# A3.B4: Depth × Magnitude Two-Way Stratification with Moho Isolation

**Source:** A2.B4  
**Reference:** A3.B3 (geospatial implementation)


## Instructions

1. Read the source case specifications file to understand how the previous spec file was structured
2. Read the reference case specifications file to understand how geospatial resources were used
3. Read the source case results file to understand the previous results and how to structure a `json` results file
4. Read the expectations for this current case A3.B4
5. Create a specifications file for this case at topic-a3/spec/case-a3-b3-spec.md :
   1. It should satisfy the expectations of this case
   2. Use the data sources for A3.B4 for analysis
   3. It should incapsulate everthing required for an agent to execute the specifications outside of the main context thread

## Source case files

**Specifications file (primary):** @topic-a2/spec/case-b4-spec.md  
**Results file (primary):** @topic-a2/output/case-b4-results.json  
**Specifications reference file (geospatial implementation):** @topic-a3/cases/case-a3-b3-details.md  


## Expectations for A3.B4

**Gap or concern (merged from two sources):**
A2.B4's console summary explicitly flagged a confound: *"Magnitude and depth are correlated (larger events tend to be deeper), so the depth pattern partly reflects the magnitude trend from A3."* The mid-crustal (20–70 km) band dominates the global signal (χ²=85.48, p=4.02×10⁻⁹), but this band also contains the highest concentration of large-magnitude events. The depth result may be largely the A2.A3 magnitude result expressed through depth. The planning-initial question asks about Moho isolation — this is incorporated as a second sub-analysis.

A3.B3 surfaces a deeper confound: the mid-crustal band (20–70 km) is precisely the depth range of slab interface seismicity. A3.B3 found that 65.8% of transitional-zone events are within 200 km of a subduction boundary (GCMT thrust enrichment ratio 1.97), and those thrust events concentrate at 20–50 km depth. The A2.B4 mid-crustal dominance may therefore be a proxy for **subduction geometry** rather than crustal loading depth — widening the original confound from depth↔magnitude to depth↔magnitude↔subduction proximity.

**Intent:** Three-component case:

1. **Two-way stratification sub-test:** Partition events by both depth band and magnitude band simultaneously. For each depth band, test whether the solar-phase signal persists after holding magnitude roughly constant (e.g., only M 6.0–6.9 events within each depth band). This disentangles the depth and magnitude effects and determines whether the mid-crustal signal is independently real or magnitude-driven. Use adaptive k based on cell size: k=24 if n≥500, k=16 if 200≤n<500, k=12 if 100≤n<200; flag n<100 as low-n (report but exclude from trend analysis).

2. **Moho isolation sub-test:** The Mohorovičić discontinuity (~20–35 km depth in continental crust, ~7–10 km in oceanic crust) separates crustal from mantle behavior. Test a narrow band around the Moho as a discrete depth class to determine whether the signal is specifically concentrated at the crust-mantle transition. Moho depth per event is assigned using the GSHHG `ocean_class` field as a tectonic setting proxy: `continental` → continental Moho band 20–35 km; `oceanic` → oceanic Moho band 5–12 km; `transitional` → exclude from Moho sub-test or report separately. This is a meaningful improvement over a fixed 25 km cutoff. See open question on external Moho dataset — a gridded dataset (e.g., CRUST1.0) would supersede the GSHHG proxy if feasible.

3. **Subduction proximity cross-tabulation sub-test:** Within the significant depth×magnitude cells (expected: mid-crustal × M6.0–6.9 and mid-crustal × M7.0–7.9), cross-tabulate with subduction proximity using the B3 PB2002 cKDTree pipeline (re-derived in B4's analysis script; do not load B3's pickle as a cross-topic dependency). This directly tests whether the mid-crustal signal persists across all tectonic settings or is confined to subduction-proximal events, distinguishing depth-driven from geometry-driven signal.

**Resolved questions:**
- *k=24 cell sizes*: Adaptive k scheme resolves this. Projected cell sizes for mid-crustal × M6.0–6.9 ≈ 3,329 (k=24 viable); intermediate × M8.0+ ≈ 52 (low-n flag). See magnitude band breakdown in spec.
- *Tectonic setting for Moho*: GSHHG `ocean_class` provides per-event proxy. Continental class → Moho 20–35 km; oceanic class → Moho 5–12 km.

**Open question (pending due-diligence):**
- Is there an external gridded Moho depth dataset (e.g., CRUST1.0 at 1°×1°, Szwillus et al. 2019, LITHO1.0) that could replace the GSHHG proxy with actual crustal thickness at each epicenter? See research note to be appended below.

## Data sources for A3.B4

- ISC-GEM raw catalog: `data/iscgem/iscgem_global_6-9_1950-2021.csv`
- GSHHG classification: `data/iscgem/plate-location/ocean_class_gshhg_global.csv` — `ocean_class`, `dist_to_coast_km` (Moho proxy)
- GCMT focal mechanism join: `data/iscgem/focal-mechanism/focal_join_global.csv` — mechanism cross-tabulation within depth×magnitude cells
- PB2002 boundary steps: `lib/PB2002_steps.dat` — SUB+OCB segments for subduction proximity (re-derive using B3 pipeline)
- *Moho dataset (TBD)*: **CRUST1.0 recommended** — see Moho dataset research note below

---

## Moho Dataset Due-Diligence Research Note

**Date:** March 4, 2026  
**Question:** Is there an external gridded Moho depth dataset that can replace the GSHHG proxy with actual per-event crustal thickness?

**Answer: Yes. CRUST1.0 is the recommended solution.**

### Evaluated datasets

| Dataset | Resolution | Format | Free | Python access | Uncertainty | Notes |
|---------|-----------|--------|------|--------------|-------------|-------|
| **CRUST1.0** | 1° × 1° | ASCII flat files + NetCDF (IRIS) | Yes | `jrleeman/Crust1.0` — NumPy only, no compiled deps | No per-cell σ | **Recommended** |
| CRUST2.0 | 2° × 2° | ASCII (bundled in Pyrocko) | Yes | `pyrocko.dataset.crust2x2.get_profile(lat, lon)` | No per-cell σ | Coarser; easy if Pyrocko installed |
| LITHO1.0 | 1° tessellated | Native tar + C++; NetCDF (IRIS) | Yes | `litho1pt0` (PyPI) — **requires gfortran** | No per-cell σ | Not viable in pure-Python environment |
| Szwillus 2019 | ~1° (kriged) | CSV (lon, lat, Moho_km, σ_km) | Yes (password on request — Kiel Univ.) | Direct CSV + `scipy.interpolate` | Yes — σ per cell | Only dataset with per-point uncertainty |
| GEMMA | 0.5° × 0.5° | NetCDF | Yes (Polimi GEOlab WPS) | `xarray` or `scipy` | ~3.4 km global RMS | Gravity-constrained, not seismic; caveat near subduction zones |
| MOHV21 | 1° × 1° | Not confirmed public | Unclear | Unknown | Propagated 5-model ensemble | Best theoretical accuracy; data availability TBD |

### Recommended implementation: CRUST1.0

**Source:** Laske, G., Masters, G., Ma, Z., and Pasyanos, M. (2013). Update on CRUST1.0 — A 1-degree Global Model of Earth's Crust. *Geophysical Research Abstracts*, 15, Abstract EGU2013-2658.  
**Download:** https://igppweb.ucsd.edu/~gabi/crust1.html (four ASCII files: `crust1.vp`, `crust1.vs`, `crust1.rho`, `crust1.bnds`) or via IRIS SAGE as NetCDF: https://ds.iris.edu/ds/products/emc-crust10/  
**Python:** https://github.com/jrleeman/Crust1.0 (pure NumPy, no compiled dependencies)

**Implementation sketch for B4 analysis script:**
```python
import numpy as np
from scipy.interpolate import RectBivariateSpline

# Load crust1.bnds: 180 rows × 360 cols × 9 layers
# Layer index 8 (0-based) = Moho depth in km (negative = depth below sea level)
bnds = np.loadtxt("crust1.bnds").reshape(180, 360, 9)
moho_grid = -bnds[:, :, 8]  # positive = depth in km
lats = np.arange(89.5, -90.5, -1.0)   # 89.5 → -89.5°N (cell centers)
lons = np.arange(-179.5, 180.5, 1.0)  # -179.5 → 179.5°E (cell centers)

# Bilinear interpolation surface
moho_interp = RectBivariateSpline(lats[::-1], lons, moho_grid[::-1], kx=1, ky=1)

# Per-event lookup
df["moho_depth_km"] = df.apply(
    lambda r: float(moho_interp(r["latitude"], r["longitude"])), axis=1
)
```

**Impact on Moho isolation sub-test:** CRUST1.0 replaces the GSHHG proxy entirely. Each event gets a continuous, spatially varying Moho depth estimate. The Moho isolation band becomes `abs(depth - moho_depth_km) <= delta` (e.g., delta = 5 km or 10 km), rather than a fixed depth range. This is a materially better approach and should be used in the spec.

**Data file to add to lib/:** Download `crust1.bnds` (the boundary file; ~2 MB uncompressed) to `lib/crust1.bnds`. The other three files (`crust1.vp`, `crust1.vs`, `crust1.rho`) are not needed for Moho depth lookup alone.

**Szwillus 2019 note:** If per-point Moho uncertainty is needed to define a defensible isolation window, Szwillus et al. (2019) (*JGR Solid Earth*, 124, 1626–1652, DOI: 10.1029/2018JB016593) provides a CSV with σ_km per grid cell. Data is available from Kiel University (password on request). This would allow filtering events where σ > delta (i.e., Moho position is too uncertain to classify the event as Moho-proximal). Not required for the initial spec but flagged as a future extension.
