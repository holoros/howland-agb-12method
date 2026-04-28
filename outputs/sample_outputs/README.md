# Sample outputs

These are the small validation summary CSVs produced by the analysis. They serve as expected output reference for reproducibility checks. Large RDS, TIF, LAS, and PNG outputs are not included; see `docs/REPRODUCE.md` for how to regenerate them.

## Files

| File | Source script | Contents |
|------|---------------|----------|
| `method09_validation.csv`         | R/method09_curve_matching.R | n, bias, bias %, RMSE, RMSE %, R² for Method 9 at 2025 |
| `method10_validation.csv`         | R/method10_directional_kohek.R | Same metrics for Method 10 (full pixel set, before holdout) |
| `method11_validation.csv`         | R/method11_hybrid_itc_anchor.R | Same metrics for Method 11 (before stratified calibration) |
| `method12_validation.csv`         | R/method12_hybrid_itc_tls.R | Same metrics for Method 12 (before stratified calibration) |
| `refinements_summary.csv`         | R/refinements_stratified_holdout_bootstrap.R | Bootstrap RMSE means + 95% CIs for all 4 methods |
| `refined_method_stats.csv`        | R/refinements_stratified_holdout_bootstrap.R | Refined stats: stratified for M11/M12, holdout for M10, with 95% CI bracket |
| `methods_11_12_ornl_validation.csv` | R/methods_11_12_ornl_sensitivity.R | ORNL 2017 anchor sensitivity test validation |
| `methods_11_12_ornl_fits.csv`     | R/methods_11_12_ornl_sensitivity.R | Stratified calibration fits with ORNL anchor (alpha, beta, SEs, R²) |
| `Table_method_summary_12methods.csv` | figures/main_figures.R | Canonical 12-method dataset used for Figs 3-6 |
| `TableS_boot_methods9_12.csv`     | figures/supplemental_figures.R | Bootstrap RMSE summary for FigS5 |
| `TableS_stratified_fits.csv`      | figures/supplemental_figures.R | Stratified calibration fits (MELiTE anchor) for FigS6 and Table S9 |

## Headline numbers (from `refined_method_stats.csv`)

```
method                    bias%   RMSE%   95% CI         R²
9 Curve matching          -21.4   35.9   [35.7, 36.2]   -0.31
10 Pure Kohek (holdout)    -1.8   26.4   [25.9, 26.8]    0.25
11 Hybrid A (stratified)   +2.6   24.7   [24.5, 25.0]    0.20
12 Hybrid B (stratified)   +1.9   25.2   [25.0, 25.5]    0.17
```

## Reproducibility check

After running the full analysis pipeline, validate against these reference CSVs:

```r
# Compare a freshly computed run to the stored reference
ref <- read.csv("outputs/sample_outputs/refined_method_stats.csv")
new <- read.csv("outputs/refinements/refined_method_stats.csv")
all.equal(ref, new, tolerance = 0.01)  # within rounding error
```
