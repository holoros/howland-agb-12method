# Methods 9 through 12: technical summary

This file describes the four AGB projection methods added in the v4.5 revision. Methods 1 through 8 are described in the main manuscript and the v3 supplemental materials. Refined numbers reported here reflect the post hoc stratified calibration, 80/20 holdout for Method 10, and 1,000 iteration bootstrap CIs documented in `R/refinements_stratified_holdout_bootstrap.R`.

## Method 9: Multi attribute yield curve matching

Reference: Tompalski et al. (2016) Forests 7, 255; Tompalski et al. (2018) Remote Sens. 10, 347.

Build a yield curve library by running AcadianGY v12.3.6 forward 100 years from nine representative starting tree lists, one for each combination of three productivity strata (Low, Medium, High) and three species compositions (coniferous, mixed, deciduous). For each 30 m pixel, select the four closest (template, year offset) candidates by absolute AGB difference at the 2015 anchor. Combine candidate AGB at each horizon as the R² weighted mean (Tompalski et al. 2016 Eq. 1). Per pixel uncertainty is the relative spread among candidates (Eq. 2).

- Wall to wall: yes
- Cost: $1.50 / acre
- 10 yr RMSE: 35.9% [35.7, 36.2] (95% bootstrap CI)
- Bias: −21.4% (shrinkage toward library centroid)

## Method 10: Wall to wall directional 3D projection (Kohek)

Reference: Kohek et al. (2022) Int. J. Appl. Earth Obs. Geoinf. 109, 102844.

Decimate the 2025 ALS from 92 to 25 pts m⁻² by random thinning (R `decimate_points(random(25))`). Segment individual tree crowns from each decimated tile in parallel via `lidR::segment_trees(..., li2012(...))`. Apply Kohek asymmetric crown radius update (Eqs. 7-9) with default coefficients M0 = 0.2 and M1 = 0.9. Compute per tree AGB via Chojnacky et al. (2014). Apply pixel level linear calibration with 80/20 train/test holdout: fit alpha and beta on the training fraction, evaluate on the held out test pixels.

- Wall to wall: yes
- Cost: $8.20 / acre (decimation + segmentation overhead)
- 10 yr RMSE (holdout): 26.4% [25.9, 26.8]
- Bias: −1.8%

## Method 11: Hybrid A (ITC 2017 anchor + default Kohek)

Substitutes the existing 2017 ALS ITC GeoPackage (208,632 wall to wall trees with per tree height and crown radius) as the starting state, eliminating the need to re segment 2025. DBH is imputed from height via a default power H D form (DBH = 1.3 + 0.8 × (Ht − 1.3)^1.1). Per tree height is projected forward by scaling by the corresponding stratum × composition top height ratio from the Method 9 library. Per tree AGB via Chojnacky. Kohek expansion uses defaults M0 = 0.2 and M1 = 0.9.

The pixel level linear calibration is fit separately for each of the nine stratum × composition combinations rather than as a single global transform. This stratified calibration recovers per pixel discrimination that a single global fit would mask.

- Wall to wall: yes
- Cost: $2.65 / acre
- 10 yr RMSE: 24.7% [24.5, 25.0] (best of the four new methods)
- Bias: +2.6%
- R²: 0.20 (up from 0.07 with single global calibration)

## Method 12: Hybrid B (ITC + TLS calibration)

Identical to Method 11 except that:

1. The H D allometry is fit from the 2025 TLS tree list (DBH = 4.12 × Ht^0.586, n = 3,510 stems with confidence ≥ 50) rather than literature values.
2. The Kohek M0 and M1 coefficients are fit from paired (2017 ALS, 2025 TLS) crown radius observations within the TLS footprint (n = 504 paired stems within 3 m proximity threshold). Calibrated values: M0 ≈ 0, M1 = 0.40.

- Wall to wall: yes
- Cost: $3.45 / acre
- 10 yr RMSE: 25.2% [25.0, 25.5]
- Bias: +1.9%
- R²: 0.17

The TLS calibrated coefficients are noticeably shallower than Kohek published defaults. This is documented as a useful negative result about the limits of naive local calibration with modest n: 504 paired observations within a single closed canopy plot were insufficient to robustly resolve shading dependence.

## Post hoc refinements

`R/refinements_stratified_holdout_bootstrap.R` implements three methodological refinements applied after Methods 9-12 ran:

1. Stratified pixel level calibration for Methods 11 and 12 (9 alpha and beta pairs per method instead of 1)
2. 80/20 train/test holdout for Method 10 (alpha and beta fit on train, validation on test)
3. 1,000 iteration nonparametric bootstrap CIs for 2025 RMSE for all four methods

Outputs: `outputs/refinements_results.rds`, `outputs/refined_method_stats.csv`.

## ORNL 2017 anchor sensitivity test

`R/methods_11_12_ornl_sensitivity.R` refits Methods 11 and 12 using the independent ORNL G LiHT 2017 AGB raster as calibration target (instead of the linearly interpolated MELiTE 2017 anchor). Result: systematically negative bias (~95%) due to a between product offset of ~30 Mg/ha between ORNL G LiHT and MELiTE ABA AGB at 2025. Documented as supplemental Table S10. Confirms that within product temporal interpolation is preferable to cross product anchoring at this site.

## Bootstrap CI caveat

The reported bootstrap intervals are pixel level resamples assuming independence. At 10 m pixel resolution adjacent forest pixels show Moran's I above 0.5; an effective sample size correction (e.g., Dutilleul 1993) would widen the intervals 2 to 3x. The reported intervals should be interpreted as lower bounds on true uncertainty. The original 4 point Chapman Richards analysis (Section 3.5) applied effective sample size correction (effective n = 15,725 to 25,244) and is the rigorous benchmark for pixel level inference.
