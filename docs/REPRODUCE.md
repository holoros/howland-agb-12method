# Reproduction guide

Step by step instructions to reproduce the v4.5 results for Methods 9 through 12.

## Prerequisites

Before running any script, ensure the prerequisites described in `DATA.md` are in place:

- `output/recalibration_results.rds` (pixel trajectories, p95 model, yield fits, agb_class)
- `output/agb_2025_melite_30m.tif` (validation raster)
- `output/Howland_ITC_2017.gpkg` (208,632 wall to wall trees)
- `output/als_2025_processing/*.las` (26 normalized tiles for Method 10)
- `data/field_inventory/{Plots_Data.csv,Tree_Data.csv}`
- `data/tls_2025/Merged_data_total.txt`
- `data/ornl_agb_rasters/Howland_AGB_2017.tif` (for sensitivity test)
- `AcadianGY_12.3.6.r` (sourced internally; expected at `$HRF_ROOT/`)

Set environment variables before submitting:

```bash
export HRF_ROOT=/users/PUOM0008/crsfaaron/HRF       # adjust to your HPC path
export OUTPUT_DIR=$HRF_ROOT/output
export FIELD_DIR=$HRF_ROOT/data/field_inventory
export TLS_DIR=$HRF_ROOT/data/tls_2025
export NORM_DIR=$HRF_ROOT/output/als_2025_processing
export FIG_DIR=$HRF_ROOT/output/manuscript_figures_v45
```

## Recommended execution order

```bash
# Step 1: Method 9 (build curve library, ~22 minutes)
sbatch slurm/submit_method09_curve_matching.sh
# Outputs: curve_matching_library.rds, curve_matching_results.rds, curve_matching_validation.csv

# Step 2: Methods 11 and 12 (ITC anchored hybrids; fast: <1 minute each)
sbatch slurm/submit_method11_hybrid_a.sh
sbatch slurm/submit_method12_hybrid_b.sh
# Outputs: method11_results.rds, method12_results.rds, method11_validation.csv, method12_validation.csv

# Step 3: Method 10 (decimation + ITC + Kohek; ~40 minutes wall time, requires 25 tiles)
sbatch slurm/submit_method10_kohek.sh
# Outputs: method10_results.rds, method10_segmented_trees_2025.rds, method10_validation.csv

# Step 4: Post hoc refinements (stratified calibration, holdout, bootstrap CIs; <30 sec)
sbatch slurm/submit_refinements.sh
# Outputs: refinements_results.rds, refined_method_stats.csv, refinements_summary.csv

# Step 5: ORNL sensitivity test (independent 2017 anchor; <30 sec)
sbatch slurm/submit_ornl_sensitivity.sh
# Outputs: methods_11_12_ornl_anchor.rds, methods_11_12_ornl_validation.csv,
#          methods_11_12_ornl_fits.csv, FigS14_residuals_methods9_12.png

# Step 6: Generate main figures (Figs 3, 4, 5, 6; <1 minute)
sbatch slurm/submit_main_figures.sh
# Outputs: Fig{3,4,5,6}_*.png and .pdf

# Step 7: Generate supplemental figures (FigS5, S6, S7; <1 minute)
sbatch slurm/submit_supplemental_figures.sh
# Outputs: FigS{5,6,7}_*.png and .pdf
```

## Expected runtime summary

| Step | Method | Wall time | CPU |
|------|--------|-----------|-----|
| 1 | 9 Curve matching | ~22 min | 8 cores |
| 2 | 11 Hybrid A | <1 min | 4 cores |
| 2 | 12 Hybrid B | <1 min | 4 cores |
| 3 | 10 Kohek (decimated) | ~40 min | 16 cores, 8 future workers |
| 4 | Refinements | <30 sec | 4 cores |
| 5 | ORNL sensitivity | <30 sec | 2 cores |
| 6 | Main figures | <1 min | 2 cores |
| 7 | Supplemental figures | <1 min | 2 cores |

Total: about 65 minutes wall clock on Cardinal.

## Expected outputs

After all steps complete, the `outputs/` directory should contain:

```
outputs/
├── method09_results/
│   ├── curve_matching_library.rds
│   ├── curve_matching_results.rds
│   └── curve_matching_validation.csv
├── method10_results/
│   ├── method10_results.rds
│   ├── method10_segmented_trees_2025.rds
│   └── method10_validation.csv
├── method11_results/
│   ├── method11_results.rds
│   └── method11_validation.csv
├── method12_results/
│   ├── method12_results.rds
│   └── method12_validation.csv
├── refinements/
│   ├── refinements_results.rds
│   ├── refined_method_stats.csv
│   └── refinements_summary.csv
├── ornl_sensitivity/
│   ├── methods_11_12_ornl_anchor.rds
│   ├── methods_11_12_ornl_validation.csv
│   └── methods_11_12_ornl_fits.csv
└── figures/
    ├── Fig3_rmse_ranking_12methods.png
    ├── Fig4_rmse_horizons_12methods.png
    ├── Fig5_scatter_12methods.png
    ├── Fig6_cost_accuracy_12methods.png
    ├── FigS5_bootstrap_methods9_12.png
    ├── FigS6_stratified_calibration.png
    ├── FigS7_itc_tls_coverage.png
    └── FigS14_residuals_methods9_12.png
```

## Validation summary expected at the end

The headline 2025 validation numbers should reproduce as:

| # | Method | Bias % | RMSE % | 95% CI | R² |
|---|--------|--------|--------|--------|-----|
| 9 | Curve matching | −21.4 | 35.9 | [35.7, 36.2] | −0.31 |
| 10 | Pure Kohek (holdout) | −1.8 | 26.4 | [25.9, 26.8] | 0.25 |
| 11 | Hybrid A (ITC 2017) | +2.6 | 24.7 | [24.5, 25.0] | 0.20 |
| 12 | Hybrid B (ITC + TLS) | +1.9 | 25.2 | [25.0, 25.5] | 0.17 |

Bootstrap CIs use seed 42; stratified calibration is deterministic; holdout uses seed 2026. Reproducibility within rounding error is expected.

## Adapting to other HPC environments

The SLURM submit scripts target Cardinal modules (`gcc/12.3.0`, `gdal/3.7.3`, `geos/3.12.0`, `proj/9.2.1`, `R/4.4.0`) and account `PUOM0008`. To run on another cluster:

1. Update `module load` lines in all `slurm/*.sh` files to match available modules
2. Replace `--account=PUOM0008` with your allocation
3. Adjust `--time`, `--mem`, `--cpus-per-task` as needed
4. Update environment variables to point at your data layout

To run interactively without SLURM, source the R scripts directly after loading modules and exporting the path environment variables:

```bash
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
export HRF_ROOT=/path/to/your/HRF; export OUTPUT_DIR=$HRF_ROOT/output; ...
Rscript --vanilla R/method11_hybrid_itc_anchor.R
```

## Troubleshooting

- **lidR LAScatalog output template error in Method 10**: ensure `decimate_one()` writes per tile output to the `dec_dir` scratch folder. The script handles this automatically; skip Method 10 setup if `als_2025_decimated_25/` is already populated.

- **Out of memory on Method 10**: reduce `future workers` from 8 to 4 in `R/method10_directional_kohek.R`. 16 GB per worker is the lower limit for 92 to 25 pts m⁻² decimation.

- **Stratified calibration falls back to ratio**: a small group like Medium_deciduous (n=31) may produce a non positive beta. The fallback (alpha = 0, beta = sum(obs) / sum(raw)) is logged and is the intended behavior.

- **Bootstrap CIs differ slightly across runs**: confirm `set.seed(42)` is in effect at the top of `boot_rmse()`. Variation beyond ±0.05 percentage points indicates an environment issue.
