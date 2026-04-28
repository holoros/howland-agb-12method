# Changelog

All notable changes to this analysis are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project tracks the manuscript revision pipeline.

## [v4.5] — 2026-04-25

Pre-submission, ready for co-author circulation.

### Added
- Method 9: multi attribute yield curve matching (Tompalski et al. 2016, 2018) wall-to-wall at 50,193 pixels
- Method 10: wall-to-wall directional 3D point cloud projection (Kohek et al. 2022) with 2025 ALS decimated from 92 to 25 pts m⁻²
- Method 11: hybrid A using the legacy 2017 ALS ITC tree list (208,632 trees) as starting state with default Kohek M0 = 0.2 and M1 = 0.9
- Method 12: hybrid B with TLS calibrated H D allometry (DBH = 4.12 × Ht^0.586) and Kohek coefficients fit from 504 paired (2017 ALS, 2025 TLS) crowns
- Stratified pixel level calibration (9 stratum × composition cells) for Methods 11 and 12, replacing single global calibration
- 80/20 train/test holdout for Method 10 (n_train = 34,899; n_test = 9,747)
- 1,000 iteration nonparametric bootstrap CIs for 2025 RMSE for all four new methods
- ORNL 2017 anchor sensitivity test (Table S10): documents that within-product temporal interpolation is preferable to cross-product anchoring
- Tompalski Eq. 1 and Eq. 2, and Kohek Eqs. 7 to 9 displayed explicitly in methods sections 2.9 and 2.10
- Spatial autocorrelation caveat for the bootstrap CIs
- Calibration leakage caveat for Methods 11 and 12 (interpolated 2017 anchor draws 20% weight from 2025 observation)
- Software availability statement with explicit package versions
- Figures 3, 4, 5, 6 regenerated with all 12 methods; Figure 6 now four-panel design (10, 20, 50 yr, average)
- Supplemental figures S11 (bootstrap RMSE distributions), S12 (stratified calibration lines), S13 (ITC + TLS coverage), S14 (residual diagnostics)
- Supplemental tables S7 (per-method parameters), S8 (bootstrap CIs), S9 / S11 (stratified calibration coefficients with SEs), S10 (ORNL sensitivity)
- 6 new references added: Bengtsson 2021, Kohek 2022, Li 2012, Roussel 2020, Tompalski 2016, Tompalski 2018 (all verified via Google Scholar)

### Changed
- Abstract restructured as a coherent 12-method narrative
- Tables 1, 2, and 4 in main manuscript extended with 4 new method rows
- Discussion 12-method mean cost ($2.72/ha) and RMSE means (20.8%, 27.0%, 35.7%) replace the v4.4 8-method values
- Figure 5 caption extended with caveat about simulated R² values for Methods 1, 4, 5, 6, 7, 8

### Fixed
- Tompalski 2018 reference title corrected from "Demonstrating the transferability..." to actual "Combining multi-date airborne laser scanning..."
- Kohek 2022 reference title corrected from "SimForest..." to actual "Simulation driven 3D forest growth forecasting..." with journal corrected from CEA to Int. J. Appl. Earth Obs. Geoinf.
- Self-citation Weiskittel et al. 2026 removed from abstract
- Results P78 changed from "the seven plot level methods" to "the twelve projection methods"

## [v4.4] — 2026-04-18

Pre-revision baseline (8-method comparison). Internal version, not in this repository.

## [v4.3] — 2026-03-26

Internal review version. Not in this repository.
