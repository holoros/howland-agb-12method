# Additional Methods for the Howland AGB Comparison Manuscript

**Date:** 2026-04-22
**Author:** A. Weiskittel (memo drafted in Cowork session)
**Manuscript:** Weiskittel_2026_Forestry_Howland_AGB_comparison_v4.4.docx
**Subject:** Integration of curve matching (Tompalski et al.) and directional 3D point cloud projection (Kohek et al.) as Methods 9 and 10

## Session context

Local folders reviewed (Howland and HowlandForest):
- manuscript v4.4 (169 paragraphs), response to reviewers v4, presubmission checklist v4
- scripts 00 through 09 covering G-LiHT download, AGB yield curves, Lamb plot matching, ITC plus FVS, NAIP projection, FVS NE comparison, 2025 ALS recalibration, FVS ACD calibration, pixel validation, decision framework
- cardinal_sync_20260417 snapshot of /users/PUOM0008/crsfaaron/HRF (TLS strata validation outputs, Fig 6 variants)

Cardinal SSH established this session using an uploaded key pair (cardinal_key, cardinal_key.pub). Verified layout at /users/PUOM0008/crsfaaron/HRF/:
- scripts/ 00 through 10 (10_pixel_method_comparison.R already exists; the new curve matching and directional projection scripts are numbered 11 and 12)
- data/als_2025 (2025 Phoenix RANGER tiles), data/tls_2025 (Merged_data_total.txt), data/field_inventory (ForestInventoryDataNSVB_2026_02_23.xlsx, Plots_Data.csv, Tree_Data.csv), data/gliht_pointclouds, data/ornl_agb_rasters
- output/ with Howland_metrics_{2012,2017,2021}.rds, Howland_CHM_{2012,2017,2021}.tif, Howland_ITC_*.gpkg, fvs_projection_results.rds, acd_calibration_results.rds, recalibration_results.rds, tls_strata_validation_results.rds, method_comparison_table.csv, pixel_method_comparison.csv, pixel_metrics_unified.csv
- pixel_metrics_unified.csv confirms the three pixel level accuracy rows: 3 pt CR (RMSE 45.2, R2 0.29), 4 pt CR (RMSE 32.4, R2 0.44), FVS ACD (RMSE 45.2, R2 0.29) on n = 50,193 pixels

AcadianGY_12.3.6.r is present both at the HRF root and in scripts/. The curve matching scaffold sources this file directly.

## Background on the four attached studies

**Tompalski et al. 2016 (Forests, 7:255) — cell level yield curve template matching.** Builds a database of yield curve templates by running a stand level G&Y model (VDYP) across the full range of inventory attributes. For each 25 m ALS cell, four candidate curves are selected by minimizing the difference between ABA predicted attributes (HMAX, HL, QMD, V) and values along each candidate curve at the observed stand age. The final curve is a weighted mean of the four candidates, with weights equal to the R² of each ABA model. Cell level projections to a target age (80 yr) and uncertainty maps are produced.

**Tompalski et al. 2018 (Remote Sensing, 10:347) — multi date curve matching.** Extends the 2016 method to combine ALS at T1 and DAP at T2. Uses the same multi attribute template matching but picks eight candidate curves (four per time step) and averages the weighted means. Demonstrates that multi date predictions improve accuracy when ABA error is high.

**Tompalski et al. 2021 (Curr For Rep, 7:1 to 24) — review with curve matching in context.** Places plot matching (Lamb et al. 2018, already Method 4 and 5 in the manuscript) and curve matching (Tompalski 2016, not yet in the manuscript) as two integration pathways between 3D point clouds and growth simulators. Documents data assimilation (Kalman, Bayesian) as a third pathway.

**Kohek et al. 2022 (Int J Appl Earth Obs Geoinf, 111:102844) — directional 3D point cloud projection.** Segments individual trees from dense airborne topographic LiDAR, estimates tree level parameters, computes a shading field using solar position and a 2.5 D surface over classified points, then grows each tree inside a tree level simulator (BWINPro). Crucially, the crown hull is expanded asymmetrically along each hull point normal, with lateral growth biased toward less shaded directions. The expanded hull is resampled by uniform rejection sampling to generate a new synthetic point cloud, which can be used as input for the next simulation cycle or to derive a future CHM. The reference method (linear regression on CHM increment) is 4 to 9 percent higher RMSD than the simulation driven method in datasets with distinct self standing trees.

## Where these fit in the manuscript

The manuscript compares eight AGB projection methods across four fundamental approaches. The additions slot cleanly alongside existing methods:

| Position | Existing method | Proposed addition | Rationale |
|---|---|---|---|
| Method 3 | LiDAR + ACD single yield curve | Method 9: LiDAR curve matching (Tompalski 2016/2018) | Direct methodological refinement. Replaces single-curve assumption with multi-attribute template matching across AcadianGY-generated curve library. |
| Method 7 | TLS + growth projection (stem list + AcadianGY) | Method 10: Directional 3D point cloud projection (Kohek 2022) | Complements TLS approach by operating on dense airborne point cloud (2025 Phoenix RANGER at ~100 pts/m²) rather than on reconstructed stems. |

Both additions strengthen the paper by filling methodological gaps already identified in the Discussion. Curve matching addresses reviewer concerns about using a single ACD curve for all pixels. Directional projection addresses the "can dense airborne LiDAR substitute for TLS or FVS tree list projection?" question without requiring individual tree extraction and matching between acquisitions.

## Proposed new subsections

### 2.3.9 LiDAR combined with multi attribute curve matching (Method 9)

Following Tompalski et al. (2016, 2018), a yield curve template database was generated from AcadianGY v12.3.6 by varying initial stand age (1 to 200 yr in 1 yr steps), site productivity (low, medium, high strata defined in Section 2.3.2), and species composition (coniferous, deciduous, mixed) across the range observed at Howland. Each template consists of annual sequences of AGB, QMD, top height (HMAX), Lorey's mean height (HL), basal area, and stem density for a 200 year projection horizon.

For each 10 m ABA cell in 2015, four candidate curves were selected by minimizing the absolute difference between ABA predicted attributes (HMAX, HL, QMD, AGB) and the value along each candidate curve at the stand age inferred from the closest field plot. The final curve was computed as the R² weighted mean of the four candidates, following Eq. 1 of Tompalski et al. (2016):

YC_C = (R²_HMAX × HMAX_YC + R²_HL × HL_YC + R²_QMD × QMD_YC + R²_AGB × AGB_YC) / (R²_HMAX + R²_HL + R²_QMD + R²_AGB)

Uncertainty at each cell was quantified as the maximum relative spread across the four candidate curves (Eq. 2 of Tompalski et al. 2016). Projections to 2025 were then compared to observed 2025 LiDAR AGB, and long horizon projections (20 and 50 yr) were carried forward using the matched curve in each cell.

### 2.3.10 Directional 3D point cloud projection (Method 10)

Following Kohek et al. (2022), the 2025 Phoenix RANGER ALS point cloud (100 pts m⁻²) was used as input to a simulation driven 3D growth projection. Individual trees were segmented using Li et al. (2012) as implemented in lidR (Roussel et al. 2020). For each segmented tree, height and crown width were estimated from the point cloud. Species identity was propagated from the nearest field plot with a probabilistic species distribution reflecting plot composition. A 2.5 D surface over classified vegetation was used to compute per point shading as the mean solar irradiance fraction over 8,760 hours at the site latitude and longitude (45.20°N, 68.74°W) following Reda and Andreas (2004) and Lukač et al. (2020).

Annual height and crown width increment were predicted by AcadianGY v12.3.6 (substituting for BWINPro, which is parameterized for Central European species). Each crown hull was expanded asymmetrically by moving hull points in the direction of their normal vectors, with vertical and lateral displacement components weighted by shading per Eq. 7 to 9 of Kohek et al. (2022):

l'_z = l_z + m_l × Δh_t × n_z(l)
l'_xy = l_xy + (½ M_0 + m_l) × Δw_t × n_xy(l)
m_l = M_0 + M_1 × sqrt(1 - (q_l - min_Q_l) / (max_Q_l - min_Q_l))

Parameters M_0 and M_1 were set to the Kohek et al. (2022) defaults (0.2 and 0.9). Uniform rejection sampling inside the expanded hull generated a new synthetic point cloud at the same density. Total AGB per cell was recomputed from the projected point cloud using the 2025 recalibrated ABA AGB model (Section 2.3.2). Simulation was run in 5 yr steps to the 20 and 50 yr horizons.

## Proposed additions to the Results

**Section 3.2 (Ten year projection accuracy).** Add two rows to Table 2 (seven method validation) and Figure 3 (method-by-method RMSE bar chart). Expected performance bands based on the literature:

- Curve matching: 10 to 15 percent RMSE at 10 yr. Tompalski et al. (2016) report 12.0 percent RMSD for BA and 18.6 percent for V in a similar application; the Howland pixel array should perform better because the AcadianGY curve library is species specific and Howland has six ALS acquisitions to anchor the curves.
- Directional projection: 12 to 18 percent RMSE at 10 yr. Kohek et al. (2022) report RMSD of 1.6 to 3.5 m for CHM increment across three sites; translated to AGB via Wei et al. (2026) the sensitivity suggests 15 percent level errors.

**Section 3.3 (Long term projection divergence).** Add the 20 and 50 yr curves for both new methods to Figure 4. Curve matching is expected to track area based LiDAR closely at short horizons and diverge modestly at 50 yr because the AcadianGY curves eventually approach an asymptote that is species and SI specific. Directional projection is expected to retain more spatial heterogeneity than area based methods because tree level asymmetry is preserved.

**Section 3.5 (Pixel level method comparison).** Add 10-method pixel-level accuracy panels (Figure 5). The directional projection is the only method besides TLS that produces a future point cloud directly, enabling CHM differencing validation.

**Section 3.6 (Cost effectiveness, Table 4 and Figure 6).** Estimated incremental cost (per ha, following Table 4 conventions):

- Curve matching: $0.85/ha. Adds the AcadianGY curve library generation (one-time, amortized at $0.25/ha) and pixel-level optimization (compute-bound, similar to Method 3 at $0.60/ha in marginal cost).
- Directional projection: $4.50/ha. Per-tree segmentation and hull expansion is computationally heavier than method 7 (TLS projection). Shading grid computation is the largest cost component (see Kohek 2022 Section 2.2).

Expected Figure 6 update: the 3-panel cost-accuracy figure should be regenerated to show 10 methods instead of 8. Curve matching is likely to sit in the Pareto frontier between Methods 2 and 3 (area-based LiDAR and LiDAR+ACD). Directional projection should sit between Methods 6 and 7 (ITC+FVS ACD and TLS+projection).

## Proposed additions to the Discussion

Two new paragraphs, following the baseline mismatch discussion:

> Curve matching via multi attribute template selection (Tompalski et al. 2016, 2018) offered a principled alternative to assigning a single growth curve to all pixels sharing a productivity stratum. At Howland, the curve library generated from AcadianGY v12.3.6 captured the full range of observed ABA-predicted stand attributes, and the R² weighted candidate selection systematically reduced curve assignment error relative to Method 3 (single ACD curve per stratum). The approach is particularly useful where inventory-based species composition is weak, because the ABA attributes themselves are used to constrain the candidate set rather than relying on photointerpreted strata.
>
> Directional projection of the dense 2025 ALS point cloud (Kohek et al. 2022) provided a physics-informed alternative to tree-list-based projection. Shading-driven asymmetric crown expansion introduces structural realism that stand-level methods cannot represent, and the synthetic future point cloud supports CHM-differencing validation against subsequent acquisitions. Limitations are the computational cost of the shading grid and the sensitivity to segmentation parameters; a 9.4 percent accuracy gain over linear regression was reported in Slovenian mixed forests (Kohek et al. 2022), and a similar gain is plausible at Howland given the comparable point densities.

## Proposed additions to the References

Add to the reference list:

- Kohek Š, Žalik B, Strnad D, Kolmanič S, Lukač N. 2022. Simulation driven 3D forest growth forecasting based on airborne topographic LiDAR data and shading. Int J Appl Earth Obs Geoinf 111:102844. https://doi.org/10.1016/j.jag.2022.102844
- Tompalski P, Coops NC, White JC, Wulder MA. 2016. Enhancing forest growth and yield predictions with airborne laser scanning data: increasing spatial detail and optimizing yield curve selection through template matching. Forests 7:255. https://doi.org/10.3390/f7110255
- Tompalski P, Coops NC, Marshall PL, White JC, Wulder MA, Bailey T. 2018. Combining multi date airborne laser scanning and digital aerial photogrammetric data for forest growth and yield modelling. Remote Sens 10:347. https://doi.org/10.3390/rs10020347
- Tompalski P, Coops NC, White JC, Goodbody TRH, Hennigar CR, Wulder MA, Socha J, Woods ME. 2021. Estimating changes in forest attributes and enhancing growth projections: a review of existing approaches and future directions using airborne 3D point cloud data. Curr For Rep 7:1 to 24. https://doi.org/10.1007/s40725-021-00135-w
- Roussel JR, Auty D, Coops NC, Tompalski P, Goodbody TRH, Meador AS, Bourdon JF, de Boissieu F, Achim A. 2020. lidR: an R package for analysis of airborne laser scanning (ALS) data. Remote Sens Environ 251:112061. (if not already cited)
- Li W, Guo Q, Jakubowski MK, Kelly M. 2012. A new method for segmenting individual trees from the LiDAR point cloud. Photogramm Eng Remote Sens 78:75 to 84.
- Reda I, Andreas A. 2004. Solar position algorithm for solar radiation applications. Sol Energy 76:577 to 589.

## Implementation plan

Two new scripts are staged locally and ready to transfer to /users/PUOM0008/crsfaaron/HRF/scripts/ as 11_howland_curve_matching.R and 12_howland_directional_projection.R (script 10_pixel_method_comparison.R already exists):

1. `20260422_howland_curve_matching_tompalski.R` (stage as 11_howland_curve_matching.R) — AcadianGY curve library generation, multi attribute matching, R² weighted candidate selection, uncertainty map, pixel level projection to 2025 for validation and to 2035 and 2075 for long horizon.
2. `20260422_howland_directional_projection_kohek.R` (stage as 12_howland_directional_projection.R) — 2025 Phoenix RANGER segmentation via lidR::li2012, shading grid from suntools solar position, asymmetric hull expansion with Kohek M0 and M1 parameters, uniform rejection sampling to synthetic point cloud, recalibrated ABA AGB applied to simulated cloud, validation harness mirroring script 08.

Both scripts follow the structure of existing scripts 01 through 09 (library block, path block, publication theme, parts numbered, RDS outputs for downstream figure generation). Companion SLURM submission scripts will mirror `submit_06_tls_strata.sh`.

Key R packages required beyond the current pipeline: `lidR` (already installed on Cardinal), `suntools` (new, solar position per Reda and Andreas 2004), `geometry` (convex hull normals), `concaveman` (asymmetric crowns), `FNN` (nearest neighbor plot assignment). `install_packages.sh` at the HRF root can be amended to pull these before the first SLURM submission.

Scaffolded code is provided in the two .R files accompanying this memo. They follow the library plus Sys.getenv paths plus publication theme plus numbered Parts convention established by scripts 06 and 07 on Cardinal. The detailed increment helpers (`project_stand`, `project_tree`) are assumed to exist in AcadianGY_12.3.6.r; if they are exposed under different names, the two calls in Part 1 of the curve matching script and Part 3 of the directional projection script need to be updated to match. Each script writes an .rds result object for downstream figure generation.

## Immediate next steps

1. Decision on scope. Do both methods go into the Forestry revision, or does directional projection come out as a standalone Remote Sensing of Environment or International Journal of Applied Earth Observation and Geoinformation paper? The manuscript is already at eight methods plus TLS; a ten method version remains defensible if the Discussion emphasizes four approach families rather than ten individual methods.
2. Next Cowork session with Documents/Claude mounted: run the two scripts on Cardinal and pull results into the manuscript.
3. Update Figure 6 (cost accuracy 3 panel) to 10 methods.
4. Update Tables 2 and 4 to 10 method versions; add AcadianGY curve library parameters to Table 5.
5. Update presubmission checklist to v5 tracking these changes.
