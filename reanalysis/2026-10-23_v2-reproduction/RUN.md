# 2026-10-23 version 2 reproduction and climate site index check (Cardinal)

`howland_v2_reproduce.py` reproduces, from the public DAAC 2490 version 2 maps and DigiSylva 2.0 plots, the map means,
G-LiHT rates and trend, project offsets, references, the area-based curve, deciles of observed change binned on 2017,
2015 and 2012, loss screen, map pair agreement, three-cornered hat, paired contrasts, and the plot-level agreement at the
plots first measured in 2026. Outputs: `summary.json`, `contrasts.csv`, `deciles.csv`, `plot_scores.csv`; inputs in `inputs.md5`.

`csi_check.R` (SLURM job 15794617) runs AcadianGY 12.3.6 on the 93 DigiSylva Howland plots at climate site index 10, 12
and 14 m for 50 years; `csi_summary.csv` holds the 8 year gross growth and mortality totals. The per-plot annual
trajectories (`csi_trajectories.csv`, 1.4 MB) are in the Howland project archive and the Zenodo record, not here.

`howland_reanalysis_checks.py` (README in `README_reanalysis_checks.md`) is the turnkey script for the coauthors' own
per-cell and per-plot tables: offset-adjusted rescoring, block bootstrap at 90, 270, 450 and 540 m with Holm
adjustment, deciles binned on 2012, 2015 and 2017, disturbance sensitivity, and plot-level scoring.
