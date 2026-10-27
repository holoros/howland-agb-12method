# Howland reanalysis checks

This script runs the checks requested after the October 22, 2026 review of the Howland AGB projection manuscript. It reads two CSV files exported from the current analysis (version 2 maps, 70 plots) and writes one table per check, so it does not depend on how the original pipeline is organized. It needs Python 3 with numpy, pandas, and scipy.

## Inputs

`cells.csv` has one row per scored 30 m cell, with AGB in Mg C ha⁻¹. Columns are `cell_id`, `x`, `y` (projected cell center in m), `map_2012`, `map_2015`, `map_2017`, `map_2021`, `map_2023`, `map_2025`, `complete_2021` (1 or 0), optional `disturbed` (the >20% loss screen) and `indep_mask` (disturbance from Landsat, Sentinel-2, or the 2024 imagery), and one `pred_<method>_2021` and `pred_<method>_2025` pair per method (for example `pred_2_2021`, `pred_6c_2025`). The script builds the static and constant-increment references itself.

`plots.csv` has `plot_id`, `cell_id` (the cell containing the plot center), `agb_2021`, `agb_2025`, and `calibrated` (1 if the plot calibrated the version 2 maps).

## Run

```
python howland_reanalysis_checks.py --cells cells.csv --plots plots.csv --out results --reps 2000 --refit
```

`--refit` refits the area-based curve and the constant increment inside every bootstrap replicate. `--offset 4.8` fixes the 2025 offset; by default it is estimated per cell from the G-LiHT trend. A run with 2,000 replicates and four block sizes takes several minutes.

## Outputs

`00_summary.md` gives the seed, the constant increment, the offset, and the number of blocks at each size. `01_offset_rescoring.csv` scores every method against the 2025 map with and without the offset. `02_bootstrap_contrasts.csv` gives paired RMSE contrasts against both references at 90, 270, 450, and 540 m, with block counts, intervals, Monte Carlo error of the endpoints, and Holm-adjusted p-values. `03_deciles_by_base_map.csv` gives residuals by decile binned on the 2017, 2015, and 2012 maps. `04_disturbance_sensitivity.csv` rescores without flagged cells. `05_plot_level_scores.csv` scores every method at the independent plots, the calibration plots, and all plots.

## Not covered

The FVS-ACD reruns (ingrowth on in the main run, climate site index of 10 and 14 m, simulated AGB restricted to trees of 10 cm DBH or more) and the refit of NAIP methods 1a to 1d use the AcadianGY run scripts in the Zenodo record and the NAIP pipeline, so they are not part of this script.
