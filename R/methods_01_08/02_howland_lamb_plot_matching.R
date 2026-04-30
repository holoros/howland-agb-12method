# =============================================================================
# Title: Lamb Plot Matching for Howland Research Forest
# Author: A. Weiskittel
# Date: 2026-03-17
# Description: Implements the Lamb et al. (2018) plot matching approach for
#              Howland Forest. Assigns tree lists from field inventory plots
#              to LiDAR grid cells based on nearest-neighbor matching of
#              area-based LiDAR metrics and inventory attributes. Uses
#              Chojnacky et al. (2014) allometric equations for tree-level AGB
#              computation. Enables grid-cell level yield forecasting using
#              existing growth models.
# Dependencies: Howland field inventory (Plots_Data.csv, Tree_Data.csv),
#               LiDAR metrics rasters, NSVB allometric data
# Reference: Lamb, S.M., MacLean, D.A., Hennigar, C.R., & Pitt, D.G. (2018).
#            Forecasting Forest Inventory Using Imputed Tree Lists for LiDAR
#            Grid Cells and a Tree-List Growth Model. Forests, 9(4), 167.
#            Chojnacky, D.C., Heath, L.S., & Jenkins, J.C. (2014). Updated
#            generalized biomass equations for North American tree species.
#            Forestry, 87(1), 129-151.
# =============================================================================

# Load required libraries
library(tidyverse)
library(terra)        # raster processing
library(sf)           # spatial data
library(FNN)          # fast nearest neighbor search
library(readxl)       # read xlsx files
library(patchwork)
library(ggtext)

# Define paths
data_dir   <- Sys.getenv("FIELD_DIR",
                          unset = "/home/aweiskittel/Documents/MAINE/DATA/HowlandForest")
raster_dir <- Sys.getenv("ORNL_DIR",
                          unset = file.path(data_dir, "agb_rasters"))
out_dir    <- Sys.getenv("FIG_DIR",
                          unset = "output/howland_lamb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Publication theme
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# =============================================================================
# PART 1: Load and prepare field inventory data (reference library)
# =============================================================================

# Load field inventory data
plots_raw <- read_csv(file.path(data_dir, "Plots_Data.csv"), show_col_types = FALSE)
trees_raw <- read_csv(file.path(data_dir, "Tree_Data.csv"), show_col_types = FALSE)

# Load NSVB reference data for comparison
nsvb_path <- file.path(data_dir, "ForestInventoryDataNSVB_2026_02_23.xlsx")
nsvb_data <- if (file.exists(nsvb_path)) {
  read_xlsx(nsvb_path) |>
    mutate(plot_id = paste0("HL_", str_pad(as.character(Plot), 3, pad = "0"))) |>
    select(plot_id, everything())
} else {
  tibble()
}

# Prepare tree data
trees <- trees_raw |>
  mutate(
    plot_id = str_extract(Plot, "^[A-Z]+_\\d+"),
    site = str_extract(Plot, "^[A-Z]+"),
    is_snag = !is.na(Snag) & Snag == "X"
  ) |>
  filter(!is_snag, DBH_cm > 0)

# Regional height-diameter equations for Maine species
# Chapman-Richards form: H = 1.3 + b0 * (1 - exp(-b1 * DBH))^b2
# Calibrated to NE FIA data
hd_params <- tribble(
  ~species, ~b0,   ~b1,    ~b2,   ~common_name,
  "PIRU",   22.0,  0.030,  1.20,  "Red spruce",
  "TSCA",   24.0,  0.025,  1.10,  "Eastern hemlock",
  "ABBA",   18.0,  0.035,  1.30,  "Balsam fir",
  "PIST",   30.0,  0.020,  1.00,  "White pine",
  "THOC",   15.0,  0.025,  1.20,  "Northern white cedar",
  "ACRU",   22.0,  0.030,  1.10,  "Red maple",
  "ACSA",   24.0,  0.025,  1.00,  "Sugar maple",
  "BEAL",   22.0,  0.025,  1.10,  "Yellow birch",
  "BEPA",   20.0,  0.030,  1.10,  "Paper birch",
  "BEPO",   20.0,  0.030,  1.10,  "Gray birch",
  "FAGR",   22.0,  0.025,  1.10,  "American beech",
  "POGR",   22.0,  0.030,  1.10,  "Bigtooth aspen",
  "POTR",   20.0,  0.035,  1.20,  "Quaking aspen",
  "QURU",   24.0,  0.020,  1.00,  "Red oak",
  "PIAB",   25.0,  0.030,  1.20,  "Norway spruce",
  "PIMA",   18.0,  0.035,  1.20,  "Black spruce",
  "POBA",   22.0,  0.030,  1.10,  "Balsam poplar",
  "FRNI",   20.0,  0.025,  1.10,  "Black ash",
  "FRPE",   22.0,  0.025,  1.10,  "White ash",
  "OSVI",   18.0,  0.025,  1.10,  "Ironwood"
)

default_hd <- tibble(b0 = 20.0, b1 = 0.028, b2 = 1.10)

# Predict heights using Chapman-Richards
trees <- trees |>
  left_join(hd_params |> select(species, hd_b0 = b0, hd_b1 = b1, hd_b2 = b2),
            by = c("Species" = "species")) |>
  mutate(
    hd_b0 = coalesce(hd_b0, default_hd$b0),
    hd_b1 = coalesce(hd_b1, default_hd$b1),
    hd_b2 = coalesce(hd_b2, default_hd$b2),
    ht_m = 1.3 + hd_b0 * (1 - exp(-hd_b1 * DBH_cm))^hd_b2
  ) |>
  select(-hd_b0, -hd_b1, -hd_b2)

# Chojnacky et al. (2014) allometric equations
# Format: ln(biomass_kg) = b0 + b1 * ln(DBH_cm)
# Species-specific coefficients for northeastern tree species
chojnacky_params <- tribble(
  ~species_group, ~b0,      ~b1,     ~species_list,
  "Picea",        -2.0773,  2.3323,  c("PIRU", "PIAB", "PIMA"),
  "Abies",        -2.5356,  2.4349,  c("ABBA"),
  "Tsuga",        -2.5356,  2.4349,  c("TSCA"),
  "Pinus",        -2.6177,  2.4638,  c("PIST"),
  "Thuja",        -2.0773,  2.3323,  c("THOC"),
  "Acer",         -1.9123,  2.3651,  c("ACRU", "ACSA"),
  "Betula",       -2.5356,  2.4349,  c("BEAL", "BEPA", "BEPO"),
  "Fagus",        -1.9123,  2.3651,  c("FAGR"),
  "Populus",      -2.4441,  2.4561,  c("POGR", "POTR"),
  "Quercus",      -2.0773,  2.3323,  c("QURU"),
  "Fraxinus",     -1.9123,  2.3651,  c("FRNI", "FRPE")
)

# Create lookup table for species to Chojnacky group
# Unnest the species_list column so each species maps to its correct group
species_lookup <- chojnacky_params |>
  unnest(cols = species_list) |>
  rename(species = species_list) |>
  select(species, choj_b0 = b0, choj_b1 = b1, species_group)

# Compute tree-level AGB using Chojnacky et al. (2014)
# ln(biomass_kg) = b0 + b1 * ln(DBH_cm)
trees <- trees |>
  left_join(species_lookup, by = c("Species" = "species")) |>
  mutate(
    choj_b0 = coalesce(choj_b0, -2.0773),      # default to Picea/Thuja
    choj_b1 = coalesce(choj_b1, 2.3323),
    agb_kg = exp(choj_b0 + choj_b1 * log(DBH_cm)),
    ba_m2 = pi * (DBH_cm / 200)^2
  )

# Get actual plot areas from CSV with coalesce fallback
plot_areas <- plots_raw |>
  mutate(plot_id = str_extract(Plot, "^[A-Z]+_\\d+")) |>
  select(plot_id, area_m2 = Area_m2)

# Compute plot-level stand metrics (matching variables)
stand_metrics <- trees |>
  left_join(plot_areas, by = "plot_id") |>
  mutate(area_m2 = coalesce(area_m2, 250)) |>
  group_by(plot_id, site, area_m2) |>
  summarise(
    n_trees = n(),
    tph = n_trees / (first(area_m2) / 10000),
    tpa = tph * 0.4047,
    ba_m2_ha = sum(ba_m2) / (first(area_m2) / 10000),
    ba_ft2_ac = ba_m2_ha * 4.356,
    qmd_cm = sqrt(ba_m2_ha / tph * 40000 / pi),
    qmd_in = qmd_cm / 2.54,
    agb_Mg_ha = sum(agb_kg) / (first(area_m2) / 10000) / 1000,
    agb_tons_ac = agb_Mg_ha * 0.4461,
    loreys_ht_m = sum(ba_m2 * ht_m) / sum(ba_m2),
    loreys_ht_ft = loreys_ht_m * 3.281,
    n_top = ceiling(0.01 * (10000 / first(area_m2)) * n_trees),
    top_ht_m = mean(sort(ht_m, decreasing = TRUE)[1:max(1, n_top)]),
    top_ht_ft = top_ht_m * 3.281,
    pct_softwood = 100 * sum(ba_m2[species_group %in% c("Picea", "Abies", "Tsuga", "Pinus", "Thuja")]) / sum(ba_m2),
    dominant_species = names(which.max(table(Species))),
    .groups = "drop"
  ) |>
  replace_na(list(pct_softwood = 50))

cat("\n--- Stand metrics summary (Chojnacky AGB, n =", nrow(stand_metrics), "plots) ---\n")
stand_metrics |>
  summarise(
    across(c(tph, ba_m2_ha, qmd_cm, agb_Mg_ha, loreys_ht_m, pct_softwood),
           list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  ) |>
  pivot_longer(everything()) |>
  print(n = 20)

# Merge in NSVB data for comparison if available
if (nrow(nsvb_data) > 0) {
  # Aggregate NSVB data to plot level before joining (tree level rows -> plot sums)
  nsvb_plot <- nsvb_data |>
    group_by(plot_id) |>
    summarise(
      # NSVB outputs dry biomass (kg), not carbon; convert kg/m2 -> Mg/ha = * 10
      agb_nsvb = sum(`AGB_2025_kgC/m2`, na.rm = TRUE) * 10,
      .groups = "drop"
    )

  cat("NSVB plots aggregated:", n_distinct(nsvb_plot$plot_id),
      "plots from", nrow(nsvb_data), "tree-level records\n")

  stand_metrics <- stand_metrics |>
    left_join(nsvb_plot, by = "plot_id")
}

# =============================================================================
# PART 2: Compute LiDAR metrics for grid cells (target cells)
# =============================================================================

# Load pre-computed LiDAR metrics if available, otherwise simulate
compute_lidar_metrics_from_rasters <- function(raster_dir, year = 2023) {
  metric_files <- list(
    hmean = file.path(raster_dir, paste0("Howland_hmean_", year, ".tif")),
    hmax  = file.path(raster_dir, paste0("Howland_hmax_", year, ".tif")),
    p50   = file.path(raster_dir, paste0("Howland_p50_", year, ".tif")),
    p75   = file.path(raster_dir, paste0("Howland_p75_", year, ".tif")),
    p90   = file.path(raster_dir, paste0("Howland_p90_", year, ".tif")),
    cc    = file.path(raster_dir, paste0("Howland_cc_", year, ".tif"))
  )

  avail <- sapply(metric_files, file.exists)
  if (all(avail)) {
    stack <- rast(unlist(metric_files))
    names(stack) <- names(metric_files)
    return(stack)
  } else {
    cat("LiDAR metric rasters not found. Simulating from field data.\n")
    return(NULL)
  }
}

lidar_metrics_rast <- compute_lidar_metrics_from_rasters(raster_dir)

# Simulate LiDAR metrics from field inventory if rasters unavailable
if (is.null(lidar_metrics_rast)) {
  set.seed(42)

  n_grid_cells <- 2000

  grid_cells <- tibble(
    cell_id = 1:n_grid_cells,
    loreys_ht_m = rnorm(n_grid_cells,
                        mean = mean(stand_metrics$loreys_ht_m, na.rm = TRUE),
                        sd = sd(stand_metrics$loreys_ht_m, na.rm = TRUE)),
    qmd_cm = rnorm(n_grid_cells,
                   mean = mean(stand_metrics$qmd_cm, na.rm = TRUE),
                   sd = sd(stand_metrics$qmd_cm, na.rm = TRUE)),
    ba_m2_ha = rnorm(n_grid_cells,
                     mean = mean(stand_metrics$ba_m2_ha, na.rm = TRUE),
                     sd = sd(stand_metrics$ba_m2_ha, na.rm = TRUE)),
    tph = rnorm(n_grid_cells,
                mean = mean(stand_metrics$tph, na.rm = TRUE),
                sd = sd(stand_metrics$tph, na.rm = TRUE)),
    agb_Mg_ha = rnorm(n_grid_cells,
                      mean = mean(stand_metrics$agb_Mg_ha, na.rm = TRUE),
                      sd = sd(stand_metrics$agb_Mg_ha, na.rm = TRUE)),
    pct_softwood = runif(n_grid_cells, 20, 100)
  ) |>
    mutate(across(c(loreys_ht_m, qmd_cm, ba_m2_ha, tph, agb_Mg_ha),
                  ~pmax(., 0.1)))

  grid_cells <- grid_cells |>
    mutate(
      p50 = loreys_ht_m * runif(n_grid_cells, 0.6, 0.8) + rnorm(n_grid_cells, 0, 1),
      p75 = loreys_ht_m * runif(n_grid_cells, 0.8, 0.95) + rnorm(n_grid_cells, 0, 0.8),
      p90 = loreys_ht_m * runif(n_grid_cells, 0.95, 1.1) + rnorm(n_grid_cells, 0, 0.5),
      hmean = loreys_ht_m * runif(n_grid_cells, 0.65, 0.85) + rnorm(n_grid_cells, 0, 0.8),
      hmax = loreys_ht_m * runif(n_grid_cells, 1.1, 1.4) + rnorm(n_grid_cells, 0, 1),
      cc = pmin(0.99, pmax(0.1, 0.3 + 0.5 * (ba_m2_ha / max(ba_m2_ha)) +
                              rnorm(n_grid_cells, 0, 0.05)))
    )
}

# =============================================================================
# PART 3: Compute matching variables for field plots (LiDAR-equivalent)
# =============================================================================

plot_lidar_equiv <- stand_metrics |>
  mutate(
    p50 = map_dbl(plot_id, ~{
      hts <- trees |> filter(plot_id == .x) |> pull(ht_m)
      quantile(hts, 0.50)
    }),
    p75 = map_dbl(plot_id, ~{
      hts <- trees |> filter(plot_id == .x) |> pull(ht_m)
      quantile(hts, 0.75)
    }),
    p90 = map_dbl(plot_id, ~{
      hts <- trees |> filter(plot_id == .x) |> pull(ht_m)
      quantile(hts, 0.90)
    }),
    hmean = map_dbl(plot_id, ~{
      hts <- trees |> filter(plot_id == .x) |> pull(ht_m)
      mean(hts)
    }),
    hmax = map_dbl(plot_id, ~{
      hts <- trees |> filter(plot_id == .x) |> pull(ht_m)
      max(hts)
    }),
    cc = pmin(0.99, ba_m2_ha / 45)
  )

# =============================================================================
# PART 4: k-Nearest Neighbor matching (the Lamb approach)
# =============================================================================

match_vars <- c("ba_m2_ha", "agb_Mg_ha", "loreys_ht_m", "qmd_cm", "cc", "pct_softwood")

# Prepare reference library
ref_data <- plot_lidar_equiv |>
  select(plot_id, site, all_of(match_vars)) |>
  drop_na()

# Standardize matching variables
ref_means <- colMeans(ref_data[match_vars])
ref_sds   <- apply(ref_data[match_vars], 2, sd)

ref_scaled <- scale(ref_data[match_vars], center = ref_means, scale = ref_sds)

# Prepare target grid cells
if (exists("grid_cells")) {
  target_data <- grid_cells |>
    select(cell_id, all_of(match_vars))
  target_scaled <- scale(target_data[match_vars], center = ref_means, scale = ref_sds)
} else {
  cat("Using raster-derived metrics for target cells\n")
}

# Run k-NN matching
k <- 5

knn_result <- get.knnx(ref_scaled, target_scaled, k = k)

# Record matched plots and distances
matches <- tibble(
  cell_id = rep(target_data$cell_id, k),
  neighbor_rank = rep(1:k, each = nrow(target_data)),
  match_plot_idx = as.vector(knn_result$nn.index),
  match_distance = as.vector(knn_result$nn.dist)
) |>
  mutate(
    match_plot_id = ref_data$plot_id[match_plot_idx],
    match_site = ref_data$site[match_plot_idx],
    weight = 1 / (match_distance + 0.001)
  ) |>
  group_by(cell_id) |>
  mutate(weight_norm = weight / sum(weight)) |>
  ungroup()

cat("\n--- Match quality summary ---\n")
matches |>
  filter(neighbor_rank == 1) |>
  summarise(
    mean_dist = mean(match_distance),
    median_dist = median(match_distance),
    max_dist = max(match_distance),
    pct_same_site = 100 * mean(match_site == "HL")
  ) |>
  print()

# =============================================================================
# PART 5: Impute tree lists to grid cells
# =============================================================================

impute_treelist_single <- function(cell_id, matches_df, trees_df) {
  best <- matches_df |>
    filter(cell_id == !!cell_id, neighbor_rank == 1)

  if (nrow(best) == 0) return(tibble())

  trees_df |>
    filter(plot_id == best$match_plot_id) |>
    mutate(
      source_cell_id = cell_id,
      match_distance = best$match_distance,
      match_weight = best$weight_norm
    )
}

demo_cells <- sample(target_data$cell_id, min(100, nrow(target_data)))

imputed_treelists <- map_dfr(demo_cells, ~{
  impute_treelist_single(.x, matches, trees)
})

cat("\n--- Imputed tree list summary (", length(demo_cells), "cells) ---\n")
imputed_treelists |>
  group_by(source_cell_id) |>
  summarise(
    n_trees = n(),
    mean_dbh_cm = mean(DBH_cm),
    dominant_sp = names(which.max(table(Species))),
    .groups = "drop"
  ) |>
  summarise(
    mean_trees_per_cell = mean(n_trees),
    mean_dbh_cm = mean(mean_dbh_cm),
    n_unique_treelists = n_distinct(source_cell_id)
  ) |>
  print()

# =============================================================================
# PART 6: Distance-weighted k-NN imputation (for stand attributes)
# =============================================================================

impute_stand_metrics_knn <- function(matches_df, ref_metrics) {
  # Use distinct to avoid many-to-many when plot_id has site duplicates
  ref_unique <- ref_metrics |>
    distinct(plot_id, .keep_all = TRUE) |>
    select(plot_id, ba_m2_ha, agb_Mg_ha,
           loreys_ht_m, qmd_cm, tph, pct_softwood)
  matches_df |>
    left_join(ref_unique, by = c("match_plot_id" = "plot_id")) |>
    group_by(cell_id) |>
    summarise(
      ba_m2_ha_imp = sum(weight_norm * ba_m2_ha),
      agb_Mg_ha_imp = sum(weight_norm * agb_Mg_ha),
      loreys_ht_imp = sum(weight_norm * loreys_ht_m),
      qmd_cm_imp = sum(weight_norm * qmd_cm),
      tph_imp = sum(weight_norm * tph),
      pct_sw_imp = sum(weight_norm * pct_softwood),
      n_unique_plots = n_distinct(match_plot_id),
      mean_distance = mean(match_distance),
      .groups = "drop"
    )
}

imputed_metrics <- impute_stand_metrics_knn(matches, stand_metrics)

cat("\n--- Imputed stand metrics (n =", nrow(imputed_metrics), "cells) ---\n")
imputed_metrics |>
  summarise(
    across(c(ba_m2_ha_imp, agb_Mg_ha_imp, loreys_ht_imp, qmd_cm_imp),
           list(mean = ~round(mean(.), 1), sd = ~round(sd(.), 1)),
           .names = "{.col}_{.fn}")
  ) |>
  pivot_longer(everything()) |>
  print(n = 10)

# =============================================================================
# PART 7: Leave-one-out cross validation
# =============================================================================

loo_cv <- function(ref_data, match_vars) {
  n <- nrow(ref_data)
  predictions <- tibble()

  for (i in 1:n) {
    train <- ref_data[-i, ]
    test  <- ref_data[i, ]

    train_scaled <- scale(train[match_vars])
    test_scaled  <- scale(test[match_vars],
                          center = attr(train_scaled, "scaled:center"),
                          scale = attr(train_scaled, "scaled:scale"))

    knn <- get.knnx(train_scaled, matrix(test_scaled, nrow = 1), k = min(5, n - 1))

    weights <- 1 / (knn$nn.dist[1, ] + 0.001)
    weights <- weights / sum(weights)

    pred_agb <- sum(weights * train$agb_Mg_ha[knn$nn.index[1, ]])
    pred_ba  <- sum(weights * train$ba_m2_ha[knn$nn.index[1, ]])
    pred_ht  <- sum(weights * train$loreys_ht_m[knn$nn.index[1, ]])

    predictions <- bind_rows(predictions, tibble(
      plot_id = ref_data$plot_id[i],
      obs_agb = ref_data$agb_Mg_ha[i],
      pred_agb = pred_agb,
      obs_ba = ref_data$ba_m2_ha[i],
      pred_ba = pred_ba,
      obs_ht = ref_data$loreys_ht_m[i],
      pred_ht = pred_ht
    ))
  }

  return(predictions)
}

cv_results <- loo_cv(ref_data, match_vars)

# Validation statistics function
validate <- function(observed, predicted) {
  resid <- observed - predicted
  tibble(
    n = length(resid),
    bias = round(mean(resid), 2),
    bias_pct = round(100 * mean(resid) / mean(observed), 1),
    rmse = round(sqrt(mean(resid^2)), 2),
    rmse_pct = round(100 * sqrt(mean(resid^2)) / mean(observed), 1),
    r2 = round(1 - sum(resid^2) / sum((observed - mean(observed))^2), 3)
  )
}

cat("\n--- Leave-one-out cross validation (Chojnacky AGB) ---\n")
cat("AGB:\n")
validate(cv_results$obs_agb, cv_results$pred_agb) |> print()
cat("Basal area:\n")
validate(cv_results$obs_ba, cv_results$pred_ba) |> print()
cat("Lorey's height:\n")
validate(cv_results$obs_ht, cv_results$pred_ht) |> print()

# =============================================================================
# PART 8: AGB Estimation Comparison (Chojnacky vs NSVB)
# =============================================================================

if (nrow(nsvb_data) > 0) {
  cat("\n--- AGB Comparison: Chojnacky field estimates vs NSVB ---\n")

  comparison <- stand_metrics |>
    select(plot_id, agb_Mg_ha, agb_nsvb) |>
    drop_na() |>
    mutate(
      diff_Mg_ha = agb_Mg_ha - agb_nsvb,
      diff_pct = 100 * diff_Mg_ha / agb_nsvb
    )

  cat("Summary statistics:\n")
  comparison |>
    summarise(
      n = n(),
      mean_chojnacky = round(mean(agb_Mg_ha), 2),
      mean_nsvb = round(mean(agb_nsvb), 2),
      mean_diff = round(mean(diff_Mg_ha), 2),
      mean_diff_pct = round(mean(diff_pct), 1),
      rmse = round(sqrt(mean(diff_Mg_ha^2)), 2)
    ) |>
    print()

  cat("\nPer-plot comparison:\n")
  print(comparison |> arrange(plot_id))
}

# =============================================================================
# PART 9: Figures
# =============================================================================

# Figure 1: LOO observed vs predicted
p_agb <- ggplot(cv_results, aes(x = pred_agb, y = obs_agb)) +
  geom_point(alpha = 0.6, color = "#009E73", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2", linewidth = 0.7) +
  labs(
    x = "Predicted AGB (Mg ha<sup>-1</sup>)",
    y = "Observed AGB (Mg ha<sup>-1</sup>)",
    title = "A) AGB: LOO cross validation"
  ) +
  coord_equal() +
  theme_pub +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown())

p_ba <- ggplot(cv_results, aes(x = pred_ba, y = obs_ba)) +
  geom_point(alpha = 0.6, color = "#E69F00", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2", linewidth = 0.7) +
  labs(
    x = "Predicted BA (m<sup>2</sup> ha<sup>-1</sup>)",
    y = "Observed BA (m<sup>2</sup> ha<sup>-1</sup>)",
    title = "B) Basal area: LOO cross validation"
  ) +
  coord_equal() +
  theme_pub +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown())

p_ht <- ggplot(cv_results, aes(x = pred_ht, y = obs_ht)) +
  geom_point(alpha = 0.6, color = "#D55E00", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2", linewidth = 0.7) +
  labs(
    x = "Predicted Lorey's Ht (m)",
    y = "Observed Lorey's Ht (m)",
    title = "C) Lorey's height: LOO cross validation"
  ) +
  coord_equal() +
  theme_pub

# Figure 2: Match distance distribution
p_dist <- ggplot(matches |> filter(neighbor_rank == 1),
                 aes(x = match_distance)) +
  geom_histogram(bins = 30, fill = "#56B4E9", color = "white") +
  geom_vline(xintercept = median(matches$match_distance[matches$neighbor_rank == 1]),
             linetype = "dashed", color = "red") +
  labs(
    x = "Euclidean distance to nearest neighbor (standardized)",
    y = "Number of grid cells",
    title = "D) Distribution of match distances (k=1)"
  ) +
  theme_pub

# Combined validation figure
combined <- (p_agb + p_ba) / (p_ht + p_dist) +
  plot_annotation(
    title = "Howland Research Forest: Lamb Plot Matching Validation (Chojnacky AGB)",
    subtitle = paste0("Reference library: ", nrow(ref_data),
                      " plots | Target: ", nrow(target_data), " grid cells | k = ", k),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(out_dir, "howland_lamb_validation.png"), combined,
       width = 28, height = 24, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "howland_lamb_validation.pdf"), combined,
       width = 28, height = 24, units = "cm")

# Figure 3: Imputed attribute distributions vs field
p_comp <- bind_rows(
  stand_metrics |> mutate(source = "Field plots") |> select(source, agb_Mg_ha),
  imputed_metrics |> mutate(source = "Imputed cells", agb_Mg_ha = agb_Mg_ha_imp) |>
    select(source, agb_Mg_ha)
) |>
  ggplot(aes(x = agb_Mg_ha, fill = source)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Field plots" = "#009E73", "Imputed cells" = "#56B4E9")) +
  labs(
    x = "Aboveground biomass (Mg ha<sup>-1</sup>)",
    y = "Density",
    title = "Distribution of AGB: field plots vs. imputed grid cells",
    fill = NULL
  ) +
  theme_pub +
  theme(axis.title.x = element_markdown())

ggsave(file.path(out_dir, "howland_agb_distributions.png"), p_comp,
       width = 17.5, height = 10, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "howland_agb_distributions.pdf"), p_comp,
       width = 17.5, height = 10, units = "cm")

# Figure 4: Chojnacky vs NSVB comparison (if data available)
if (nrow(nsvb_data) > 0 & nrow(comparison) > 0) {
  p_nsvb <- ggplot(comparison, aes(x = agb_nsvb, y = agb_Mg_ha)) +
    geom_point(alpha = 0.7, color = "#CC79A7", size = 2.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "#0072B2", linewidth = 0.8) +
    labs(
      x = "NSVB AGB (Mg ha<sup>-1</sup>)",
      y = "Chojnacky field estimate (Mg ha<sup>-1</sup>)",
      title = "AGB Comparison: Chojnacky vs NSVB"
    ) +
    coord_equal() +
    theme_pub +
    theme(axis.title.x = element_markdown(), axis.title.y = element_markdown())

  ggsave(file.path(out_dir, "howland_chojnacky_vs_nsvb.png"), p_nsvb,
         width = 17.5, height = 14, units = "cm", dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "howland_chojnacky_vs_nsvb.pdf"), p_nsvb,
         width = 17.5, height = 14, units = "cm")
}

cat("\n--- Figures saved to:", out_dir, "---\n")
cat("Script completed successfully.\n")
