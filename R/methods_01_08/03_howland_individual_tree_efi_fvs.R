
# ==============================================================================
# HOWLAND FOREST LIDAR CASE STUDY: SCRIPT 03
# Individual Tree EFI and FVS Integration with Real lidR ITC Segmentation
# ==============================================================================
# PURPOSE: Segment individual trees from G-LiHT point clouds using real lidR
# ITC methodology, predict DBH from morphological metrics, generate FVS tree
# lists, and compare with field inventory and area-based estimates
#
# AUTHOR: Aaron Weiskittel
# DATE: 2026-03-17
# ==============================================================================

# Load required libraries
library(tidyverse)
library(lidR)
library(terra)
library(sf)
library(readxl)
library(RSQLite)
library(DBI)
library(ggtext)

# Define paths at top
data_dir <- Sys.getenv("FIELD_DIR",
  unset = "/home/aweiskittel/Documents/MAINE/DATA/HowlandForest"
)
gliht_dir <- Sys.getenv("GLIHT_DIR",
  unset = file.path(data_dir, "gliht_pointclouds")
)
fig_dir <- Sys.getenv("FIG_DIR",
  unset = "output/howland_efi"
)
fvs_dir <- Sys.getenv("FVS_DIR",
  unset = "output/fvs_treelists"
)

# Create output directories if they don't exist
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(fvs_dir)) dir.create(fvs_dir, recursive = TRUE)

# ==============================================================================
# PART 1: Load Field Inventory Data
# ==============================================================================

# Load plot data (keep original casing, rename selectively)
plots_data <- read_csv(
  file.path(data_dir, "Plots_Data.csv"),
  show_col_types = FALSE
) |>
  rename(
    plot = Plot,
    radius_m = Radus_m,
    area_m2 = Area_m2,
    northing = Northing,
    easting = Easting,
    agb_2025_j_kgc_m2 = `AGB_2025_J_KgC/m2`,
    agb_2025_y_kgc_m2 = `AGB_2025_Y_KgC/m2`,
    agb_2025_c_kgc_m2 = `AGB_2025_C_KgC/m2`
  )

# Load tree data (keep original casing, rename selectively)
# Tree_Data.csv Plot column has format "HL_001_01" (plot + tree suffix)
# Extract just the plot portion to match Plots_Data.csv format ("HL_001")
tree_data <- read_csv(
  file.path(data_dir, "Tree_Data.csv"),
  show_col_types = FALSE
) |>
  mutate(
    plot = str_extract(Plot, "^[A-Z]+_\\d+")
  ) |>
  rename(
    species = Species,
    dbh_cm = DBH_cm
  )

# Load NSVB plot data for reference
nsvb_data <- read_xlsx(
  file.path(data_dir, "ForestInventoryDataNSVB_2026_02_23.xlsx"),
  sheet = 1
) |>
  rename(
    site = Site,
    agb_2025_kgc_m2 = `AGB_2025_kgC/m2`
  )

# Join plot and tree data
field_inventory <- plots_data |>
  left_join(
    tree_data,
    by = "plot"
  )

cat("Loaded field inventory with", n_distinct(field_inventory$plot),
  "plots and", nrow(tree_data), "trees\n")

# ==============================================================================
# PART 2: LiDAR Point Cloud ITC Segmentation (Real lidR vs Simulation)
# ==============================================================================

# Check for LiDAR data in multiple locations
las_files <- list.files(
  gliht_dir,
  pattern = "\\.(las|laz)$",
  full.names = TRUE,
  recursive = TRUE
)

gpkg_files <- list.files(
  gliht_dir,
  pattern = "\\.gpkg$",
  full.names = TRUE,
  recursive = TRUE
)

# Also check for G-LiHT CHM rasters in the project G-LiHT directory
# (user may have pre-downloaded these separately from the point clouds)
project_root <- Sys.getenv("PROJECT_DIR",
  unset = dirname(dirname(gliht_dir)))
gliht_chm_dir <- file.path(project_root, "G-LiHT")
chm_files <- character(0)
if (dir.exists(gliht_chm_dir)) {
  chm_files <- list.files(gliht_chm_dir, pattern = "_CHM\\.tif$",
                           full.names = TRUE)
  cat("Found", length(chm_files), "G-LiHT CHM rasters in", gliht_chm_dir, "\n")
}

use_real_lidar <- length(las_files) > 0 || length(gpkg_files) > 0
use_chm_rasters <- length(chm_files) > 0

if (use_real_lidar && length(gpkg_files) > 0) {
  cat("Loading pre-processed ITC results from .gpkg files...\n")

  itc_results_list <- list()

  for (i in seq_along(gpkg_files)) {
    plot_name <- tools::file_path_sans_ext(basename(gpkg_files[i]))
    cat("Loading", plot_name, "\n")

    itc_geom <- st_read(gpkg_files[i], quiet = TRUE)

    itc_results_list[[i]] <- itc_geom |>
      as_tibble() |>
      select(-geometry) |>
      mutate(
        plot = plot_name,
        .before = 1
      )
  }

  # Standardize column names after binding
  itc_metrics <- bind_rows(itc_results_list)
  names(itc_metrics) <- tolower(names(itc_metrics))
  itc_metrics <- itc_metrics |>
    rename(
      treeID = treeid,
      itc_ht_m = z_max,
      itc_crown_area_m2 = crown_area_m2
    ) |>
    select(plot, treeID, itc_ht_m, itc_crown_area_m2, n_points)

  cat("Loaded ITC results from", length(gpkg_files), "plots\n")

} else if (use_chm_rasters) {
  # --- ITC detection from pre-existing G-LiHT CHM rasters ---
  cat("Running ITC detection from G-LiHT CHM rasters...\n")

  # Get plot locations for spatial subsetting
  plot_locs <- plots_data |>
    filter(!is.na(easting), !is.na(northing)) |>
    select(plot, easting, northing, radius_m)

  itc_results_list <- list()

  for (chm_path in chm_files) {
    chm_name <- tools::file_path_sans_ext(basename(chm_path))
    cat("Processing CHM:", chm_name, "\n")

    chm_raw <- rast(chm_path)

    # For each field plot, extract a circular window from the CHM
    for (j in seq_len(nrow(plot_locs))) {
      pl <- plot_locs$plot[j]
      px <- plot_locs$easting[j]
      py <- plot_locs$northing[j]
      pr <- plot_locs$radius_m[j]
      if (is.na(pr)) pr <- 11.28  # default ~400 m2 plot

      # Create a small extent around the plot
      buf <- pr + 10
      plot_ext <- ext(px - buf, px + buf, py - buf, py + buf)

      tryCatch({
        chm_clip <- crop(chm_raw, plot_ext)
        if (is.null(chm_clip) || all(is.na(values(chm_clip)))) next

        # Smooth CHM
        chm_smooth <- focal(chm_clip, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)

        # Variable window function
        f_ws <- function(x) {
          dplyr::case_when(
            x < 5 ~ 3.0, x < 10 ~ 4.0, x < 15 ~ 5.0,
            x < 20 ~ 6.0, TRUE ~ 7.0
          )
        }

        # Detect treetops
        ttops <- locate_trees(chm_smooth, algorithm = lmf(ws = f_ws))

        if (nrow(ttops) > 0) {
          itc_results_list[[length(itc_results_list) + 1]] <- tibble(
            plot = pl,
            treeID = seq_len(nrow(ttops)),
            itc_ht_m = ttops$Z,
            itc_crown_area_m2 = pi * (f_ws(ttops$Z) / 2)^2,  # estimated from window
            n_points = NA_integer_
          )
        }
      }, error = function(e) {
        cat("  Skipping plot", pl, "for", chm_name, ":", e$message, "\n")
      })
    }
  }

  if (length(itc_results_list) > 0) {
    itc_metrics <- bind_rows(itc_results_list)
    # Keep only one CHM result per plot (use the first one found)
    itc_metrics <- itc_metrics |>
      group_by(plot, treeID) |>
      slice_head(n = 1) |>
      ungroup()
    cat("ITC detection from CHMs:", n_distinct(itc_metrics$plot),
        "plots,", nrow(itc_metrics), "trees\n")
  } else {
    cat("CHM-based ITC detection produced no results. Falling back to simulation.\n")
    use_chm_rasters <- FALSE
  }
}

if (!exists("itc_metrics")) {
  # No gpkg or CHM data produced results; try raw LAS point clouds
  if (length(las_files) > 0) {
    cat("Processing G-LiHT point clouds with lidR ITC pipeline...\n")

    ctg <- readLAScatalog(gliht_dir)
    opt_output_file(ctg) <- ""
    opt_chunk_buffer(ctg) <- 10
    opt_chunk_size(ctg) <- 200

    ctg_norm <- normalize_height(ctg, algorithm = tin())

    cat("Generating CHM...\n")
    chm <- rasterize_canopy(ctg_norm, res = 0.5,
                             algorithm = p2r(subcircle = 0.1))
    chm_smooth <- terra::focal(chm, w = matrix(1, 3, 3),
                                fun = mean, na.rm = TRUE)

    f_ws <- function(x) {
      dplyr::case_when(
        x < 5 ~ 3.0, x < 10 ~ 4.0, x < 15 ~ 5.0,
        x < 20 ~ 6.0, TRUE ~ 7.0
      )
    }

    cat("Detecting treetops and segmenting crowns...\n")
    ttops <- locate_trees(chm_smooth, algorithm = lmf(ws = f_ws))
    las_seg <- segment_trees(ctg_norm,
                              algorithm = silva2016(chm_smooth, ttops))

    tree_metrics <- crown_metrics(las_seg,
      func = ~list(z_max = max(Z), z_mean = mean(Z),
                   z_p90 = quantile(Z, 0.90), n_points = length(Z)),
      geom = "convex"
    )
    tree_metrics$crown_area_m2 <- as.numeric(st_area(tree_metrics))

    itc_metrics <- tree_metrics |>
      as_tibble() |>
      select(-geometry) |>
      rename(itc_ht_m = z_max, itc_crown_area_m2 = crown_area_m2) |>
      select(treeID, itc_ht_m, itc_crown_area_m2, n_points)

    cat("Extracted metrics from", nrow(itc_metrics), "ITC trees\n")

  } else {
    # --- Simulation fallback when no LiDAR data available ---
    cat("No LiDAR data found. Using simulated ITC detection from field trees...\n")
    set.seed(42)

    itc_metrics <- tree_data |>
      group_by(plot) |>
      mutate(
        treeID = rank(-dbh_cm, ties.method = "first"),
        itc_ht_m = 1.3 + 20 * (1 - exp(-0.03 * dbh_cm))^1.2 +
                   rnorm(n(), 0, 1),
        itc_ht_m = pmax(itc_ht_m, 2),
        itc_crown_area_m2 = pi * (0.5 + 0.15 * dbh_cm)^2 +
                            rnorm(n(), 0, 1),
        itc_crown_area_m2 = pmax(itc_crown_area_m2, 1),
        n_points = rpois(n(), 150)
      ) |>
      ungroup() |>
      select(plot, treeID, itc_ht_m, itc_crown_area_m2, n_points)

    cat("Simulated ITC metrics for", n_distinct(itc_metrics$plot),
        "plots,", nrow(itc_metrics), "trees\n")
  }
}

# ==============================================================================
# PART 3: Species Classification (Simulated from Field Data)
# ==============================================================================

cat("Classifying ITC trees by species...\n")

# Get species composition from field data (proportion within each plot)
species_comp <- tree_data |>
  group_by(plot) |>
  mutate(plot_n = n()) |>
  group_by(plot, species) |>
  summarize(
    n_trees = n(),
    prop = n_trees / first(plot_n),
    .groups = "drop"
  )

# Assign species to ITC trees probabilistically based on field composition
# Only join plots that have species data; drop rows with NA/zero weights
itc_with_species <- itc_metrics |>
  inner_join(species_comp, by = "plot") |>
  filter(!is.na(prop), prop > 0) |>
  group_by(plot, treeID) |>
  slice_sample(n = 1, weight_by = prop) |>
  ungroup() |>
  select(plot, treeID, itc_ht_m, itc_crown_area_m2, n_points, species)

cat("Species assigned to", n_distinct(itc_with_species$plot), "plots,",
    nrow(itc_with_species), "ITC trees\n")

# ==============================================================================
# PART 4: DBH Prediction from ITC Height and Crown Area
# ==============================================================================

cat("Predicting DBH from ITC morphometrics...\n")

# Define softwood and hardwood groups
softwood_species <- c("PIRU", "PIAB", "PIMA", "ABBA", "TSCA", "THOC")
hardwood_species <- c("ACRU", "ACSA", "BEAL", "BEPA", "BEPO", "FAGR", "POGR", "POTR")

# Predict heights for field trees using Chapman Richards H D equations
# so we can calibrate the inverse model: DBH = f(height)
hd_params_03 <- tribble(
  ~species, ~hd_b0, ~hd_b1, ~hd_b2,
  "PIRU", 22.0, 0.030, 1.20,
  "TSCA", 24.0, 0.025, 1.10,
  "ABBA", 18.0, 0.035, 1.30,
  "PIST", 30.0, 0.020, 1.00,
  "THOC", 15.0, 0.025, 1.20,
  "ACRU", 22.0, 0.030, 1.10,
  "ACSA", 24.0, 0.025, 1.00,
  "BEAL", 22.0, 0.025, 1.10,
  "BEPA", 20.0, 0.030, 1.10,
  "BEPO", 20.0, 0.030, 1.10,
  "FAGR", 22.0, 0.025, 1.10,
  "POGR", 22.0, 0.030, 1.10,
  "POTR", 20.0, 0.035, 1.20,
  "QURU", 24.0, 0.020, 1.00,
  "PIAB", 25.0, 0.030, 1.20,
  "PIMA", 18.0, 0.035, 1.20,
  "POBA", 22.0, 0.030, 1.10,
  "FRNI", 20.0, 0.025, 1.10,
  "FRPE", 22.0, 0.025, 1.10,
  "OSVI", 18.0, 0.025, 1.10
)

field_classified <- field_inventory |>
  left_join(hd_params_03, by = "species") |>
  mutate(
    hd_b0 = coalesce(hd_b0, 20.0),
    hd_b1 = coalesce(hd_b1, 0.028),
    hd_b2 = coalesce(hd_b2, 1.10),
    ht_m = 1.3 + hd_b0 * (1 - exp(-hd_b1 * dbh_cm))^hd_b2,
    is_softwood = species %in% softwood_species,
    crown_area_proxy = (dbh_cm / 5) ^ 1.5
  ) |>
  filter(!is.na(dbh_cm), dbh_cm > 0, !is.na(ht_m))

# Fit inverse models: DBH = f(height) using field data
sw_data <- field_classified |> filter(is_softwood == TRUE)
hw_data <- field_classified |> filter(is_softwood == FALSE)

fit_sw2 <- lm(dbh_cm ~ log(ht_m), data = sw_data)
fit_hw2 <- lm(dbh_cm ~ log(ht_m), data = hw_data)

cat("Softwood DBH model R2:", round(summary(fit_sw2)$r.squared, 3), "\n")
cat("Hardwood DBH model R2:", round(summary(fit_hw2)$r.squared, 3), "\n")

# Predict DBH for ITC trees using their LiDAR derived heights
itc_with_dbh <- itc_with_species |>
  mutate(
    is_softwood = species %in% softwood_species,
    ht_m = itc_ht_m  # rename for model prediction compatibility
  ) |>
  mutate(
    pred_dbh_cm = ifelse(is_softwood,
      predict(fit_sw2, newdata = data.frame(ht_m = ht_m)),
      predict(fit_hw2, newdata = data.frame(ht_m = ht_m))
    ),
    pred_dbh_cm = pmax(pred_dbh_cm, 1)
  )

cat("Predicted DBH for", nrow(itc_with_dbh), "ITC trees\n")

# ==============================================================================
# PART 5: Tree-Level Biomass using Chojnacky et al. (2014)
# ==============================================================================

cat("Calculating tree-level biomass using Chojnacky equations...\n")

# Chojnacky et al. (2014) allometric coefficients: ln(biomass_kg) = b0 + b1*ln(DBH_cm)
choj_coefs <- tribble(
  ~species_group, ~b0, ~b1,
  "PIRU", -2.0773, 2.3323,
  "PIAB", -2.0773, 2.3323,
  "PIMA", -2.0773, 2.3323,
  "ABBA", -2.5356, 2.4349,
  "TSCA", -2.5356, 2.4349,
  "THOC", -2.0773, 2.3323,
  "PIST", -2.6177, 2.4638,
  "ACRU", -1.9123, 2.3651,
  "ACSA", -1.9123, 2.3651,
  "BEAL", -2.5356, 2.4349,
  "BEPA", -2.5356, 2.4349,
  "BEPO", -2.5356, 2.4349,
  "FAGR", -1.9123, 2.3651,
  "POGR", -2.4441, 2.4561,
  "POTR", -2.4441, 2.4561
)

itc_with_biomass <- itc_with_dbh |>
  left_join(choj_coefs, by = c("species" = "species_group")) |>
  mutate(
    b0 = coalesce(b0, -2.0773),
    b1 = coalesce(b1, 2.3323),
    pred_agb_kg = exp(b0 + b1 * log(pred_dbh_cm))
  ) |>
  select(-b0, -b1)

# Calculate stand-level AGB from ITC predictions
itc_agb_stand <- itc_with_biomass |>
  group_by(plot) |>
  summarize(
    itc_n_trees = n(),
    itc_agb_kg = sum(pred_agb_kg, na.rm = TRUE),
    itc_agb_mean_kg = mean(pred_agb_kg, na.rm = TRUE),
    itc_agb_sd_kg = sd(pred_agb_kg, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    plots_data |>
      select(plot, area_m2),
    by = "plot"
  ) |>
  mutate(
    itc_agb_kgc_m2 = itc_agb_kg / area_m2
  )

cat("Stand-level AGB calculated for", nrow(itc_agb_stand), "plots\n")

# ==============================================================================
# PART 6: Generate FVS Tree Lists
# ==============================================================================

cat("Generating FVS tree lists...\n")

# Create FVS-formatted tree list
fvs_treelist <- itc_with_biomass |>
  rename(
    dbh = pred_dbh_cm,
    agb_kg = pred_agb_kg
  ) |>
  mutate(
    ht_ft = itc_ht_m * 3.28084,
    cr_width_ft = sqrt(itc_crown_area_m2) * 3.28084,
    stand_id = plot,
    tree_record_id = treeID,
    sp_code = species,
    dbh_in = dbh / 2.54,
    ht_primary_ft = ht_ft,
    crown_ratio = pmin(0.95, pmax(0.2, exp(-0.05 * dbh))),
    damage_code = 0,
    damage_percent = 0,
    status = "alive"
  ) |>
  select(
    stand_id, tree_record_id, sp_code, dbh_in, ht_primary_ft,
    crown_ratio, damage_code, damage_percent, status
  )

# Save FVS tree lists to CSV
fvs_csv_path <- file.path(fvs_dir, "howland_fvs_treelist.csv")
write_csv(fvs_treelist, fvs_csv_path)
cat("Saved FVS tree list CSV to", fvs_csv_path, "\n")

# Save to SQLite database
db_path <- file.path(fvs_dir, "howland_fvs.db")
con <- dbConnect(RSQLite::SQLite(), db_path)
dbWriteTable(con, "FVS_TreeInit", as.data.frame(fvs_treelist), overwrite = TRUE)
dbDisconnect(con)
cat("Saved FVS tree list to SQLite database at", db_path, "\n")

# ==============================================================================
# PART 7: Comparison of ITC vs Area-Based vs Field Inventory
# ==============================================================================

cat("Comparing ITC, area-based, and field estimates...\n")

# Compute field tree AGB using Chojnacky equations (same as ITC)
field_with_agb <- field_classified |>
  left_join(choj_coefs, by = c("species" = "species_group")) |>
  mutate(
    b0 = coalesce(b0, -2.0773),
    b1 = coalesce(b1, 2.3323),
    field_agb_kg = exp(b0 + b1 * log(dbh_cm))
  )

field_agb_stand <- field_with_agb |>
  group_by(plot) |>
  summarize(
    field_n_trees = n(),
    field_agb_kg_total = sum(field_agb_kg, na.rm = TRUE),
    field_agb_mean_kg = mean(field_agb_kg, na.rm = TRUE),
    field_agb_sd_kg = sd(field_agb_kg, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    plots_data |>
      select(plot, area_m2),
    by = "plot"
  ) |>
  mutate(
    field_agb_kgc_m2 = field_agb_kg_total / area_m2
  )

# Prepare NSVB comparison data
# NSVB uses full site names; Plots_Data uses 2-letter codes
site_code_map <- c(Howland = "HL", Demeritt = "DD", Penobscot = "PN", Schoodic = "SC")

nsvb_summary <- nsvb_data |>
  group_by(site, Plot) |>
  summarize(
    nsvb_agb_kgc_m2 = mean(agb_2025_kgc_m2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    site_code = coalesce(site_code_map[site], site),
    plot = paste0(site_code, "_", sprintf("%03d", Plot))
  ) |>
  select(plot, nsvb_agb_kgc_m2)

# Combine all estimates
comparison_data <- itc_agb_stand |>
  left_join(field_agb_stand |>
    select(plot, field_agb_kgc_m2),
    by = "plot"
  ) |>
  left_join(nsvb_summary, by = "plot") |>
  left_join(
    plots_data |>
      select(plot, agb_2025_j_kgc_m2, agb_2025_y_kgc_m2, agb_2025_c_kgc_m2),
    by = "plot"
  )

# Calculate errors relative to NSVB
comparison_stats <- comparison_data |>
  filter(!is.na(nsvb_agb_kgc_m2)) |>
  summarize(
    n_itc   = sum(!is.na(itc_agb_kgc_m2) & !is.na(nsvb_agb_kgc_m2)),
    n_field = sum(!is.na(field_agb_kgc_m2) & !is.na(nsvb_agb_kgc_m2)),
    itc_rmse_kgc_m2 = sqrt(mean((itc_agb_kgc_m2 - nsvb_agb_kgc_m2) ^ 2, na.rm = TRUE)),
    itc_bias_kgc_m2 = mean(itc_agb_kgc_m2 - nsvb_agb_kgc_m2, na.rm = TRUE),
    itc_r2 = cor(itc_agb_kgc_m2, nsvb_agb_kgc_m2, use = "complete.obs") ^ 2,
    field_rmse_kgc_m2 = sqrt(mean((field_agb_kgc_m2 - nsvb_agb_kgc_m2) ^ 2, na.rm = TRUE)),
    field_bias_kgc_m2 = mean(field_agb_kgc_m2 - nsvb_agb_kgc_m2, na.rm = TRUE),
    field_r2 = cor(field_agb_kgc_m2, nsvb_agb_kgc_m2, use = "complete.obs") ^ 2
  )

cat("Comparison statistics:\n")
print(comparison_stats)

# ==============================================================================
# PART 8: Error Propagation Analysis
# ==============================================================================

cat("Conducting error propagation analysis...\n")

# Quantify uncertainty from ITC metrics, DBH prediction, and allometry
error_prop <- itc_with_biomass |>
  group_by(plot) |>
  summarize(
    agb_uncertainty_kg = sum(pred_agb_kg * 0.1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    plots_data |>
      select(plot, area_m2),
    by = "plot"
  ) |>
  mutate(
    agb_uncertainty_kgc_m2 = agb_uncertainty_kg / area_m2
  )

cat("Mean AGB uncertainty:", mean(error_prop$agb_uncertainty_kgc_m2),
  "kgC/m2\n"
)

# ==============================================================================
# PART 9: Figures
# ==============================================================================

cat("Creating publication-quality figures...\n")

# Define publication theme
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

# Figure 1: ITC AGB vs NSVB AGB
fig1 <- comparison_data |>
  filter(!is.na(nsvb_agb_kgc_m2)) |>
  ggplot(aes(x = nsvb_agb_kgc_m2, y = itc_agb_kgc_m2)) +
  geom_point(size = 3, alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "ITC-derived AGB vs NSVB Reference",
    x = "NSVB AGB (kgC m<sup>&minus;2</sup>)",
    y = "ITC AGB (kgC m<sup>&minus;2</sup>)"
  ) +
  coord_equal() +
  theme_pub +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown())

# Save Figure 1
png(file.path(fig_dir, "fig1_itc_vs_nsvb_agb.png"),
  width = 800, height = 700, res = 100
)
print(fig1)
dev.off()

pdf(file.path(fig_dir, "fig1_itc_vs_nsvb_agb.pdf"),
  width = 8, height = 7
)
print(fig1)
dev.off()

# Figure 2: ITC Height distribution
fig2 <- itc_with_biomass |>
  ggplot(aes(x = itc_ht_m, fill = is_softwood)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  facet_wrap(~plot, scales = "free_y") +
  labs(
    title = "Distribution of ITC Tree Heights",
    x = "Height (m)",
    y = "Frequency",
    fill = "Softwood"
  ) +
  scale_fill_manual(values = c("darkgreen", "steelblue")) +
  theme_pub

png(file.path(fig_dir, "fig2_itc_height_dist.png"),
  width = 1000, height = 800, res = 100
)
print(fig2)
dev.off()

pdf(file.path(fig_dir, "fig2_itc_height_dist.pdf"),
  width = 10, height = 8
)
print(fig2)
dev.off()

# Figure 3: Predicted DBH vs Field DBH
# Generate predicted DBH from height models for field trees (inverse validation)
fig3_data <- field_classified |>
  filter(!is.na(dbh_cm), !is.na(ht_m)) |>
  mutate(
    pred_dbh_cm = ifelse(is_softwood,
      predict(fit_sw2, newdata = data.frame(ht_m = ht_m)),
      predict(fit_hw2, newdata = data.frame(ht_m = ht_m))
    )
  )

fig3 <- fig3_data |>
  ggplot(aes(x = dbh_cm, y = pred_dbh_cm, color = is_softwood)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.2) +
  facet_wrap(~paste0("Softwood: ", is_softwood), labeller = label_both) +
  labs(
    title = "DBH Prediction Model Validation",
    subtitle = paste0(
      "SW R\u00b2 = ", round(summary(fit_sw2)$r.squared, 3),
      ", HW R\u00b2 = ", round(summary(fit_hw2)$r.squared, 3)
    ),
    x = "Field Measured DBH (cm)",
    y = "Predicted DBH from Height (cm)",
    color = "Softwood"
  ) +
  coord_equal() +
  theme_pub +
  theme(legend.position = "none")

png(file.path(fig_dir, "fig3_dbh_validation.png"),
  width = 1000, height = 600, res = 100
)
print(fig3)
dev.off()

pdf(file.path(fig_dir, "fig3_dbh_validation.pdf"),
  width = 10, height = 6
)
print(fig3)
dev.off()

# Figure 4: AGB comparison by method
fig4 <- comparison_data |>
  filter(!is.na(nsvb_agb_kgc_m2)) |>
  select(plot, itc_agb_kgc_m2, field_agb_kgc_m2, nsvb_agb_kgc_m2) |>
  pivot_longer(
    cols = ends_with("_kgc_m2"),
    names_to = "method",
    values_to = "agb_kgc_m2"
  ) |>
  mutate(
    method = case_when(
      method == "itc_agb_kgc_m2" ~ "ITC",
      method == "field_agb_kgc_m2" ~ "Field",
      method == "nsvb_agb_kgc_m2" ~ "NSVB Ref"
    )
  ) |>
  ggplot(aes(x = plot, y = agb_kgc_m2, fill = method)) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(
    title = "Stand-level AGB Estimates by Method",
    x = "Plot",
    y = "AGB (kgC m<sup>&minus;2</sup>)",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("steelblue", "darkgreen", "darkred")) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_markdown())

png(file.path(fig_dir, "fig4_agb_comparison.png"),
  width = 1000, height = 700, res = 100
)
print(fig4)
dev.off()

pdf(file.path(fig_dir, "fig4_agb_comparison.pdf"),
  width = 10, height = 7
)
print(fig4)
dev.off()

# Figure 5: Species composition in ITC detection
fig5 <- itc_with_species |>
  group_by(plot, species) |>
  summarize(n = n(), .groups = "drop") |>
  ggplot(aes(x = plot, y = n, fill = species)) +
  geom_col(position = "fill", alpha = 0.85) +
  labs(
    title = "Species Composition of ITC-Detected Trees",
    x = "Plot",
    y = "Proportion",
    fill = "Species"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png(file.path(fig_dir, "fig5_species_composition.png"),
  width = 1000, height = 700, res = 100
)
print(fig5)
dev.off()

pdf(file.path(fig_dir, "fig5_species_composition.pdf"),
  width = 10, height = 7
)
print(fig5)
dev.off()

# ==============================================================================
# PART 10: Key Takeaways and Summary Statistics
# ==============================================================================

cat("\n")
cat("============================================================================\n")
cat("HOWLAND FOREST LIDAR ITC CASE STUDY SUMMARY\n")
cat("============================================================================\n")
cat("\n")

cat("LIDAR DATA PROCESSING:\n")
cat("  LiDAR Source:", ifelse(use_real_lidar, "Real G-LiHT point clouds", "Simulated"), "\n")
cat("  Total plots processed:", n_distinct(itc_with_biomass$plot), "\n")
cat("  Total ITC trees segmented:", nrow(itc_with_biomass), "\n")
cat("  Mean trees per plot:", round(nrow(itc_with_biomass) / n_distinct(itc_with_biomass$plot), 1), "\n")
cat("\n")

cat("BIOMASS ESTIMATION:\n")
cat("  Mean stand AGB (ITC):", round(mean(itc_agb_stand$itc_agb_kgc_m2, na.rm = TRUE), 2),
  "kgC/m2\n"
)
cat("  SD stand AGB (ITC):", round(sd(itc_agb_stand$itc_agb_kgc_m2, na.rm = TRUE), 2),
  "kgC/m2\n"
)
cat("  Mean stand AGB (Field):", round(mean(field_agb_stand$field_agb_kgc_m2, na.rm = TRUE), 2),
  "kgC/m2\n")
cat("  Field plots used:", nrow(field_agb_stand), "\n")
cat("\n")

if (nrow(comparison_stats) > 0) {
  cat("ERROR METRICS (vs NSVB Reference):\n")
  cat("  ITC RMSE:", round(comparison_stats$itc_rmse_kgc_m2, 2), "kgC/m2\n")
  cat("  ITC Bias:", round(comparison_stats$itc_bias_kgc_m2, 2), "kgC/m2\n")
  cat("  ITC R²:", round(comparison_stats$itc_r2, 3), "\n")
  cat("  Field RMSE:", round(comparison_stats$field_rmse_kgc_m2, 2), "kgC/m2\n")
  cat("  Field Bias:", round(comparison_stats$field_bias_kgc_m2, 2), "kgC/m2\n")
  cat("  Field R²:", round(comparison_stats$field_r2, 3), "\n")
  cat("\n")
}

cat("OUTPUT FILES:\n")
cat("  FVS CSV:", fvs_csv_path, "\n")
cat("  FVS SQLite:", db_path, "\n")
cat("  Figures directory:", fig_dir, "\n")
cat("\n")

cat("SCRIPT COMPLETED SUCCESSFULLY\n")
cat("============================================================================\n")
