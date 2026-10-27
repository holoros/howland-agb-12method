# 2026-10-26 close-items reproduction (Cardinal job 15894685)

Independent reproduction of the map-only and plot-level quantities of the October 2026 manuscript from public data:
ORNL DAAC 2490 version 2 AGB maps (2012, 2015, 2017, 2021, 2023, 2025), DigiSylva 2.0 plots (Harvard Dataverse),
and the AcadianGY 12.3.6 climate-site-index trajectories of the 2026-10-23 check (`../2026-10-23_v2-reproduction/csi_check.R`).
Input checksums are in `inputs_md5.json`. The study area polygon is not public, so the domain is the convex hull of the
Howland plots buffered by 30 m (224 ha, 1,954 30 m cells), with 0 and 60 m buffers as sensitivity; the 2025 map is
shifted 5 m east and 5 m south as in the manuscript.

What it produces (all in Mg C ha-1, carbon fraction 0.5 assumed for the AcadianGY trajectories):

- `method_scores.csv`: bias, RMSE (% of the observed mean), R2 and change efficiency at 2021 and 2025 for the static and
  constant-increment references, the area-based curve (refit on this domain and with the manuscript parameters), the
  FVS-ACD and FVS-NE curves with the manuscript's Table 3 parameters, the FVS-ACD ingrowth curve, and trajectory
  imputation; the 2025 rows are repeated with the 2025 offset removed at the domain estimate (3.0) and at 3.5, 4.3, 5.1.
- `score_intervals.csv`: 95% block bootstrap intervals (90, 270, 540 m blocks, 1,000 replicates, block counts) for RMSE and
  bias, and for offset-adjusted RMSE with the offset re-estimated in each replicate.
- `contrasts.csv`: paired RMSE contrasts of every method against the static and constant-increment references.
- `three_cornered_hat_all.csv`: three-cornered hat error estimates for all 20 map triplets.
- `plot_predictions.csv`, `plot_scores.csv`: map values and every method's prediction at the plots (3 by 3 window of 10 m
  cells around each plot center) and scores at the 29 plots first measured in 2026, the 49 measured in 2024 and 2025, and all.
- `summary.json`: domain, means, trend, offsets by trend form and buffer, loss fractions, plot against map increments,
  FVS-ACD curve refits at climate site index 10, 12 and 14 m.

Run: `sbatch submit_close_items.slurm` (python3 with rasterio, shapely, scipy, pandas) from a directory holding the six
`Howland_AGB_<year>.tif` files, `plots.tab`, `trees.tab`, and `csi_check/csi_trajectories.csv`.

These values feed Tables A4 and A5 of the manuscript revision of 2026-10-26 and are labelled there as a reproduction on
an approximate domain, not as the manuscript's own domain. The public 2025 layer differs from the layer behind the
manuscript's 2025 numbers (mean 92.2 against 93.5 Mg C ha-1, three-cornered hat error 8.9 to 13.1 against 4.2, offset
3.0 against 4.3), which is an open question for the coauthors.
