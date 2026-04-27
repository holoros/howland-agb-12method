# howland-agb-12method

Comparative assessment of twelve methods for projecting forest aboveground biomass at the Howland Research Forest, Maine, USA. Companion analysis code and outputs for Weiskittel et al. (2026), *Forestry*.

## Overview

This repository contains the R analysis code and SLURM submission scripts used to compare twelve aboveground biomass (AGB) projection methods at the Howland Research Forest:

1. NAIP imagery combined with Acadian yield curves
2. Area based LiDAR (Chapman Richards and Schumacher)
3. LiDAR combined with Acadian yield curves
4. ITC segmentation combined with strata averaging
5. ITC combined with FVS NE
6. ITC combined with FVS ACD
7. Terrestrial laser scanning combined with Acadian projection
8. Sample based field inventory
9. Multi attribute yield curve matching (Tompalski et al. 2016, 2018)
10. Wall to wall directional 3D point cloud projection (Kohek et al. 2022) with ALS decimated to 25 pts m⁻²
11. Hybrid A: 2017 ITC anchor with default Kohek expansion coefficients and stratified pixel calibration
12. Hybrid B: ITC + TLS calibration, with TLS fit H D allometry and locally calibrated Kohek M0 and M1 coefficients

Methods 1 through 8 are the original analysis (see Weiskittel et al. 2026 main text and supplementary v3). Methods 9 through 12 are added in the v4.5 revision and are the primary focus of the code in this repository.

## Headline results

| # | Method | Cost ($/ac) | RMSE 10 yr (%) | 95% CI | R² |
|---|--------|------|----------|------------|------|
| 6 | ITC + FVS-ACD (best baseline) | 0.88 | 10.0 | [6.1, 14.3] | — |
| 9 | Curve matching | 1.50 | 35.9 | [35.7, 36.2] | −0.31 |
| 10 | Pure Kohek (holdout) | 8.20 | 26.4 | [25.9, 26.8] | 0.25 |
| 11 | Hybrid A (ITC 2017, MELiTE anchor) | 2.65 | **24.7** | [24.5, 25.0] | 0.20 |
| 12 | Hybrid B (ITC + TLS, MELiTE anchor) | 3.45 | 25.2 | [25.0, 25.5] | 0.17 |

Method 11 is the best of the four new methods. None of them beats the existing ITC + FVS ACD baseline. The directional 3D family adds individual tree resolution that area based methods cannot provide; that is its operational value.

## Repository layout

```
howland-agb-12method/
├── R/                                       # Core analysis scripts (one per method)
│   ├── method09_curve_matching.R            # Method 9 (Tompalski 2016, 2018)
│   ├── method10_directional_kohek.R         # Method 10 (Kohek 2022, decimated ALS)
│   ├── method11_hybrid_itc_anchor.R         # Method 11 (ITC 2017 anchor)
│   ├── method12_hybrid_itc_tls.R            # Method 12 (ITC + TLS calibration)
│   ├── refinements_stratified_holdout_bootstrap.R  # Post hoc refinements
│   └── methods_11_12_ornl_sensitivity.R     # Sensitivity test with ORNL 2017 anchor
├── figures/
│   ├── main_figures.R                       # Figs 3, 4, 5, 6 (manuscript)
│   └── supplemental_figures.R               # Figs S5, S6, S7, S14 (supplemental)
├── slurm/                                   # SLURM submit scripts for OSC Cardinal
├── docs/
│   ├── DATA.md                              # Data sources and access
│   ├── METHODS.md                           # Method by method summary
│   └── REPRODUCE.md                         # Step by step reproduction guide
├── outputs/                                 # Empty placeholder; populated on Cardinal
├── CITATION.cff
├── LICENSE
└── README.md
```

## Compute environment

The analysis was developed and run on the Ohio Supercomputer Center (OSC) Cardinal cluster (RHEL 9, SLURM scheduler) under project allocation PUOM0008.

Required software:

- R 4.4.0
- gcc 12.3.0
- gdal 3.7.3, geos 3.12.0, proj 9.2.1
- AcadianGY v12.3.6 (sourced separately at `AcadianGY_12.3.6.r`)

Required R packages:

- tidyverse 2.0.0
- terra 1.9.1
- sf 1.0.16
- lidR 4.1.2
- FNN 1.1.4
- future 1.34.0
- future.apply 1.11.0
- ggplot2 4.0.2
- patchwork 1.2.0
- ggrepel 0.9.5
- ggtext 0.1.2
- scales 1.3.0

## Quick start

The scripts are written for the Cardinal HPC environment. Adapt paths and module loads as needed for other clusters.

```bash
# 1. Clone
git clone https://github.com/<your-org>/howland-agb-12method.git
cd howland-agb-12method

# 2. Set environment variables to point at your data
export HRF_ROOT=/path/to/HRF
export OUTPUT_DIR=$HRF_ROOT/output
export FIELD_DIR=$HRF_ROOT/data/field_inventory
export TLS_DIR=$HRF_ROOT/data/tls_2025
export NORM_DIR=$HRF_ROOT/output/als_2025_processing

# 3. Submit Method 9 first (curve matching builds the AcadianGY library)
sbatch slurm/submit_method09_curve_matching.sh

# 4. Submit hybrid Methods 11 and 12 (fast: about 30 seconds each)
sbatch slurm/submit_method11_hybrid_a.sh
sbatch slurm/submit_method12_hybrid_b.sh

# 5. Run post hoc refinements (stratified calibration, holdout, bootstrap CIs)
sbatch slurm/submit_refinements.sh

# 6. Run Method 10 (decimation + ITC + Kohek; longer: about 40 minutes)
sbatch slurm/submit_method10_kohek.sh

# 7. Generate figures
sbatch slurm/submit_main_figures.sh
sbatch slurm/submit_supplemental_figures.sh
```

See `docs/REPRODUCE.md` for the full reproduction sequence with expected outputs.

## Data

Raw data are not included in this repository because they are large and partially restricted. See `docs/DATA.md` for sources and access:

- 2012, 2014, 2015, 2017, 2021, 2023, 2025 G LiHT airborne LiDAR (NASA distributed via ORNL DAAC)
- 2025 MELiTE ABA AGB raster (validation reference)
- 2025 Phoenix RANGER ALS (~92 pts m⁻²)
- 2025 RIEGL VZ 400i terrestrial laser scan
- 56 ground inventory plots
- 2017 wall to wall ITC tree list (Howland_ITC_2017.gpkg, 208,632 trees)

## Citation

If you use this code, please cite both this repository and the manuscript:

```
Weiskittel, A. R. (2026). howland-agb-12method: comparative assessment of
  twelve AGB projection methods at Howland Research Forest [Software].
  https://github.com/<your-org>/howland-agb-12method

Weiskittel, A. R., et al. (2026). Comparative assessment of twelve methods
  for projecting forest aboveground biomass at the Howland Research Forest,
  Maine, USA. Forestry, in press.
```

See `CITATION.cff` for machine readable citation metadata.

## License

MIT (see `LICENSE`). Data products distributed under their respective licenses; see `docs/DATA.md`.

## Contact

Aaron R. Weiskittel · University of Maine · aaron.weiskittel@maine.edu

## Acknowledgments

Compute resources provided by the Ohio Supercomputer Center under allocation PUOM0008. ORNL DAAC G LiHT product hosting. NASA G LiHT mission for multi temporal acquisitions. Field inventory and TLS deployment supported by USDA Forest Service Northern Research Station.
