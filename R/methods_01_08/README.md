# Methods 1 through 8: original 8-method analysis

These scripts implement the original eight-method comparison from Weiskittel et al. v4.4 (the Forestry submission baseline). Methods 9 through 12, added in the v4.5 revision, are in the parent `R/` directory.

## Method to script mapping

| # | Method | Primary script | Notes |
|---|--------|----------------|-------|
| 1 | NAIP + ACD curve | `04_howland_naip_projection_validation.R` | NAIP imagery proxy combined with ACD growth rates |
| 2 | Area-based LiDAR (CR + Schumacher) | `01_howland_agb_yield_curves.R`, `06_recalibrate_with_2025als.R` | Chapman-Richards and Schumacher yield curves; 4-point recalibration with 2025 ALS |
| 3 | LiDAR + ACD curve | `01_howland_agb_yield_curves.R` | LiDAR initial state combined with ACD growth rates |
| 4 | ITC + strata averaging | `02_howland_lamb_plot_matching.R`, `03_howland_individual_tree_efi_fvs.R` | Lamb et al. (2018) kNN imputation with TLS-anchored strata averages |
| 5 | ITC + FVS-NE | `05_howland_fvsne_projection_comparison.R` | Stratum-averaged ITC tree list + FVS-NE growth |
| 6 | ITC + FVS-ACD | `05_howland_fvsne_projection_comparison.R`, `07_calibrate_fvsacd.R` | Same ITC inputs + Acadian variant of FVS |
| 7 | TLS + projection | `02_howland_lamb_plot_matching.R` (TLS branch) | 2025 Purdue TLS direct stem list with ACD projection |
| 8 | Field inventory | `08_pixel_validation.R` | 56-plot field inventory baseline with NSVB allometrics |

## Cross-cutting scripts

- `09_decision_framework.R` — assembles the cost-accuracy quadrant analysis used in Figure 6
- `10_pixel_method_comparison.R` — pixel-level comparison across the three wall-to-wall methods (3-pt CR, 4-pt CR, FVS-ACD)
- `generate_manuscript_figures.R` — builds Figures 1 through 13 of the v4.4 manuscript
- `reviewer_response_analyses.R` — analyses added in response to first-round reviewer comments

## Pipeline order

For a fresh reproduction of Methods 1 through 8, run in this order on Cardinal:

```bash
sbatch slurm/submit_methods_01_yield_curves.sh   # Methods 2, 3
sbatch slurm/submit_methods_02_lamb_imputation.sh # Methods 4, 7
sbatch slurm/submit_methods_03_itc_fvs.sh         # Methods 5, 6
sbatch slurm/submit_methods_04_naip.sh            # Method 1
sbatch slurm/submit_methods_05_fvsne_compare.sh   # Methods 5, 6 comparison
sbatch slurm/submit_methods_06_recalibration.sh   # Method 2 4-point recalibration
sbatch slurm/submit_methods_07_fvsacd.sh          # Method 6 calibration
sbatch slurm/submit_methods_08_field_validation.sh # Method 8
sbatch slurm/submit_methods_09_decision_framework.sh # Cost-accuracy framework
sbatch slurm/submit_methods_10_pixel_comparison.sh   # Wall-to-wall comparison
sbatch slurm/submit_methods_main_figures.sh           # Figures 1-13
```

Submit scripts for these are not duplicated in `slurm/`; the canonical SLURM scripts live on Cardinal at `/users/PUOM0008/crsfaaron/HRF/submit_*.sh`. Adapt the v4.5 submit script templates in this repo's `slurm/` directory if you need to re-run on another HPC environment.

## Compute environment

Same as the main `R/` scripts:

- R 4.4.0
- gcc 12.3.0
- gdal 3.7.3, geos 3.12.0, proj 9.2.1
- AcadianGY v12.3.6 (sourced separately at `AcadianGY_12.3.6.r`)
- FVS-NE / FVS-ACD installations (proprietary, not redistributable)

## Caveats

These scripts were developed iteratively across the v4.0 through v4.4 revision cycle and contain Cardinal-specific paths (`/users/PUOM0008/crsfaaron/HRF/`) that need adaptation for other environments. They are provided here for transparency and to enable independent reproduction; for the v4.5 analysis (Methods 9 through 12) the cleaner, environment-aware scripts in the parent `R/` directory are the recommended starting point.
