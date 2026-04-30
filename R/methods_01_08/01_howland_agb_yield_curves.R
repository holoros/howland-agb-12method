# =============================================================================
# Title: AGB Yield Curves from Repeated LiDAR at Howland Research Forest
# Author: A. Weiskittel
# Date: 2026-03-17
# Description: Develops area based AGB yield curves using repeated LiDAR
#              acquisitions (2012, 2015, 2017, 2021, 2023, 2025) from the
#              Howland Research Forest in central Maine. Demonstrates that
#              repeated LiDAR enables direct observation of AGB trajectories
#              without traditional plot remeasurement.
#
# Key updates:
#   - Validation uses NSVB allometric system (loaded from xlsx)
#   - Comparison of 4 allometric systems: Jenkins, Young, CRM, NSVB
#   - Stand level growth projections with 20 year forward trajectory
#   - Schumacher yield model as alternative (ln(AGB) = b0 + b1/age)
#   - Placeholder for 2025 AGB raster validation
#
# Dependencies: ORNL DAAC AGB rasters (Wei et al. 2025), field inventory data
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(terra)         # raster processing
library(sf)            # spatial data
library(nlme)          # nonlinear mixed effects
library(patchwork)     # multi panel figures
library(ggtext)        # markdown in axis labels
library(readxl)        # read Excel files

# --- Paths -------------------------------------------------------------------
# Local development paths (used when running interactively on your machine)
# On Cardinal HPC, these are overridden by the master script environment
data_dir   <- Sys.getenv("FIELD_DIR",
                          unset = "/home/aweiskittel/Documents/MAINE/DATA/HowlandForest")
raster_dir <- Sys.getenv("ORNL_DIR",
                          unset = file.path(data_dir, "agb_rasters"))
out_dir    <- Sys.getenv("FIG_DIR",
                          unset = "output/howland_yield")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Publication theme -------------------------------------------------------
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

# Colorblind safe palette
pal_strata <- c("Low productivity" = "#E69F00",
                "Medium productivity" = "#56B4E9",
                "High productivity" = "#009E73")

# =============================================================================
# PART 1: Load and stack multi temporal AGB rasters
# =============================================================================

# The ORNL DAAC dataset (doi: 10.3334/ORNLDAAC/2434) provides AGB maps at
# 10 m resolution for 2012, 2015, 2017, 2021, and 2023. The V2 maps from
# Wei et al. (2026) extend coverage through 2025 using Phoenix RANGER data.
#
# LiDAR sources:
#   2012: G LiHT (leaf on)
#   2015: USGS 3DEP (leaf off)
#   2017: G LiHT (leaf on)
#   2021: G LiHT (leaf on)
#   2023: USGS 3DEP (leaf off)
#   2025: Phoenix RANGER U580 (leaf on, ~100 pts/m2)

lidar_years <- c(2012, 2015, 2017, 2021, 2023, 2025)

# Load rasters into a SpatRaster stack
# Each file is a cloud optimized GeoTIFF with AGB in kgC/m2
agb_files <- paste0(raster_dir, "/Howland_AGB_", lidar_years, ".tif")

# Check which files exist AND are real rasters (not corrupt/empty downloads)
# Files under 1 KB are almost certainly failed downloads (HTML error pages)
avail <- file.exists(agb_files) &
  sapply(agb_files, function(f) {
    if (file.exists(f)) file.info(f)$size > 1000 else FALSE
  })
cat("Available rasters:\n")
data.frame(
  year = lidar_years,
  file = basename(agb_files),
  exists = file.exists(agb_files),
  valid = avail,
  size_bytes = sapply(agb_files, function(f) {
    if (file.exists(f)) file.info(f)$size else 0
  })
) |> print()

# Remove corrupt downloads (< 1 KB) so they can be re-downloaded later
corrupt <- file.exists(agb_files) & !avail
if (any(corrupt)) {
  cat("Removing", sum(corrupt), "corrupt raster files (< 1 KB)\n")
  file.remove(agb_files[corrupt])
}

if (sum(avail) >= 2) {
  # Load rasters individually (extents differ: G-LiHT covers east, 3DEP covers all)
  rast_list <- lapply(agb_files[avail], rast)

  # Crop all to common (intersection) extent and align grids
  common_ext <- ext(rast_list[[1]])
  for (r in rast_list[-1]) {
    common_ext <- intersect(common_ext, ext(r))
  }
  cat("Common extent:", as.vector(common_ext), "\n")

  ref_rast <- crop(rast_list[[1]], common_ext)
  aligned <- lapply(rast_list, function(r) {
    r_crop <- crop(r, common_ext)
    if (!compareGeom(r_crop, ref_rast, stopOnError = FALSE)) {
      r_crop <- resample(r_crop, ref_rast, method = "bilinear")
    }
    r_crop
  })

  agb_stack <- rast(aligned)
  names(agb_stack) <- paste0("AGB_", lidar_years[avail])
  cat("Stacked", nlyr(agb_stack), "rasters (", ncell(agb_stack), "cells)\n")
} else {
  cat("Fewer than 2 valid rasters found. Using simulated data for demonstration.\n")
}

# =============================================================================
# PART 2: Extract pixel level AGB trajectories
# =============================================================================

# Convert the raster stack to a data frame of pixel trajectories.
# Each row = one 10m pixel; columns = AGB at each acquisition year.
# This gives us thousands of "virtual plots" with repeated AGB measurements.

extract_pixel_trajectories <- function(agb_stack, years) {
  # Sample all non NA pixels (or subsample for very large rasters)
  n_cells <- ncell(agb_stack)
  cat("Total cells:", n_cells, "\n")

  # Extract values for all cells
  vals <- values(agb_stack)
  coords <- xyFromCell(agb_stack, 1:n_cells)

  # Build trajectory data frame
  traj <- as_tibble(vals) |>
    mutate(cell_id = 1:n(), x = coords[, 1], y = coords[, 2]) |>
    pivot_longer(
      cols = starts_with("AGB_"),
      names_to = "year_label",
      values_to = "agb_kgC_m2"
    ) |>
    mutate(
      year = as.integer(str_extract(year_label, "\\d{4}")),
      # Convert to Mg/ha: kgC/m2 * 10000 m2/ha / 1000 kg/Mg * 2 (C to biomass)
      agb_Mg_ha = agb_kgC_m2 * 10000 / 1000 * 2
    ) |>
    filter(!is.na(agb_kgC_m2))

  return(traj)
}

# If rasters are available, extract trajectories
# Otherwise, simulate demonstration data based on field plot statistics
if (exists("agb_stack")) {
  traj <- extract_pixel_trajectories(agb_stack, lidar_years[avail])
} else {
  # ---- Simulate demonstration data from field plot AGB distributions ----
  set.seed(42)
  n_pixels <- 5000

  # Use field plot AGB stats to parameterize simulation
  # From Plots_Data.csv: AGB ranges 0 16 kgC/m2, mean ~10 kgC/m2
  # Simulate pixel trajectories using Chapman Richards growth + noise
  sim_pixels <- tibble(
    cell_id = rep(1:n_pixels, each = length(lidar_years)),
    year = rep(lidar_years, n_pixels)
  ) |>
    group_by(cell_id) |>
    mutate(
      # Assign each pixel a "productivity class" (random asymptote)
      Amax = first(rnorm(1, mean = 12, sd = 4)),
      Amax = pmax(Amax, 2),
      # Assign a random "effective age" in 2012
      age_2012 = first(sample(40:140, 1)),
      age = age_2012 + (year - 2012),
      # Chapman Richards growth with pixel specific params
      k = first(rnorm(1, mean = 0.015, sd = 0.003)),
      k = pmax(k, 0.005),
      p = 2.0,
      agb_kgC_m2 = Amax * (1 - exp(-k * age))^p + rnorm(n(), 0, 0.5),
      agb_kgC_m2 = pmax(agb_kgC_m2, 0),
      agb_Mg_ha = agb_kgC_m2 * 20
    ) |>
    ungroup()

  traj <- sim_pixels
  cat("Using simulated pixel trajectories (n =", n_pixels, "pixels)\n")
}

# =============================================================================
# PART 3: Classify pixels into forest type groups
# =============================================================================

# Strategy: Use the 2023 LiDAR metrics or existing forest type classification
# to group pixels. For this demonstration, we classify by AGB level in the
# earliest acquisition year to create "productivity strata."

classify_strata <- function(traj, base_year = NULL) {
  if (is.null(base_year)) base_year <- min(traj$year)

  base_agb <- traj |>
    filter(year == base_year) |>
    select(cell_id, base_agb = agb_kgC_m2)

  # Tercile classification
  breaks <- quantile(base_agb$base_agb, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

  strata <- base_agb |>
    mutate(
      stratum = cut(base_agb,
                    breaks = breaks,
                    labels = c("Low productivity", "Medium productivity", "High productivity"),
                    include.lowest = TRUE)
    ) |>
    select(cell_id, stratum)

  traj |> left_join(strata, by = "cell_id")
}

traj <- classify_strata(traj) |>
  filter(!is.na(stratum))  # Drop pixels that fell outside tercile breaks

# =============================================================================
# PART 4: Fit AGB yield models (Chapman Richards by stratum)
# =============================================================================

# Approach 1: Fit Chapman Richards yield curves by stratum
# AGB = Amax * (1  exp( k * t))^p
# where t = time since some reference (we use year  2000 as a proxy for
# "relative time" since we lack true stand age for every pixel)

# Summarize mean AGB by stratum and year
stratum_means <- traj |>
  group_by(stratum, year) |>
  summarise(
    mean_agb = mean(agb_Mg_ha, na.rm = TRUE),
    sd_agb = sd(agb_Mg_ha, na.rm = TRUE),
    se_agb = sd_agb / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

cat("\n--- Mean AGB by stratum and year ---\n")
print(stratum_means, n = 30)

# Fit Chapman Richards by stratum
# Strategy: try 3 parameter fit first (Amax, k, p). If that fails, fix p = 2.0
# and fit only Amax and k (2 parameter model). This avoids overparameterization
# when the time series is short or covers a narrow age range.
fit_cr_by_stratum <- function(dat) {
  # Use relative time (years since 2000)
  dat <- dat |> mutate(t = year - 2000)

  fits <- dat |>
    group_by(stratum) |>
    group_map(~ {
      d <- .x
      amax_init <- max(d$agb_Mg_ha, na.rm = TRUE) * 1.3

      # Attempt 1: full 3 parameter CR
      fit_3p <- tryCatch({
        nls(
          agb_Mg_ha ~ Amax * (1 - exp(-k * t))^p,
          data = d,
          start = list(Amax = amax_init, k = 0.05, p = 2.0),
          control = nls.control(maxiter = 500, warnOnly = TRUE)
        )
      }, error = function(e) NULL)

      # Attempt 2: fix p = 2.0, fit Amax and k only (more stable)
      fit_2p <- if (is.null(fit_3p)) {
        tryCatch({
          nls(
            agb_Mg_ha ~ Amax * (1 - exp(-k * t))^2,
            data = d,
            start = list(Amax = amax_init, k = 0.05),
            control = nls.control(maxiter = 500, warnOnly = TRUE)
          )
        }, error = function(e) NULL)
      } else NULL

      # Use whichever succeeded
      fit <- if (!is.null(fit_3p)) fit_3p else fit_2p

      if (!is.null(fit)) {
        cf <- coef(fit)
        tibble(
          stratum = .y$stratum,
          Amax = cf["Amax"],
          k = cf["k"],
          p = if ("p" %in% names(cf)) cf["p"] else 2.0,
          rmse = sqrt(mean(residuals(fit)^2)),
          r2 = 1 - sum(residuals(fit)^2) / sum((d$agb_Mg_ha - mean(d$agb_Mg_ha))^2),
          n_params = length(cf)
        )
      } else {
        # Final fallback: linear model
        fit_lm <- lm(agb_Mg_ha ~ t, data = d)
        tibble(
          stratum = .y$stratum,
          Amax = NA_real_, k = NA_real_, p = NA_real_,
          rmse = sqrt(mean(residuals(fit_lm)^2)),
          r2 = summary(fit_lm)$r.squared,
          n_params = 2L
        )
      }
    }) |>
    bind_rows()

  return(fits)
}

# Subsample for model fitting (full pixel data can be very large)
traj_sample <- traj |>
  group_by(stratum) |>
  filter(cell_id %in% sample(unique(cell_id), min(1000, n_distinct(cell_id)))) |>
  ungroup()

cr_fits <- fit_cr_by_stratum(traj_sample)
cat("\n--- Chapman Richards fits by stratum ---\n")
print(cr_fits)

# =============================================================================
# PART 5: Compute PAI (periodic annual increment)
# =============================================================================

# For each pixel, compute AGB change between consecutive LiDAR acquisitions
compute_pai <- function(traj) {
  traj |>
    arrange(cell_id, year) |>
    group_by(cell_id) |>
    mutate(
      agb_prev = lag(agb_Mg_ha),
      year_prev = lag(year),
      interval = year - year_prev,
      pai_Mg_ha_yr = (agb_Mg_ha - agb_prev) / interval
    ) |>
    filter(!is.na(pai_Mg_ha_yr)) |>
    ungroup()
}

pai_data <- compute_pai(traj)

pai_summary <- pai_data |>
  group_by(stratum, year_prev, year) |>
  summarise(
    mean_pai = mean(pai_Mg_ha_yr, na.rm = TRUE),
    sd_pai = sd(pai_Mg_ha_yr, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  mutate(period = paste0(year_prev, "\u2013", year))

cat("\n--- Periodic Annual Increment by stratum and period ---\n")
print(pai_summary)

# =============================================================================
# PART 6: Stand level growth projections
# =============================================================================

# Project the fitted Chapman Richards curves forward 20 years (to 2045)
# using the fitted parameters. Also fit Schumacher yield model as alternative:
# ln(AGB) = b0 + b1/age

project_growth_curves <- function(cr_fits, stratum_means, projection_years = 20) {
  # Get current year
  current_year <- max(stratum_means$year)
  projection_year <- current_year + projection_years

  # Generate future year sequence
  future_years <- seq(current_year, projection_year, by = 1)

  # For each stratum, project using Chapman Richards parameters
  projections <- crossing(
    stratum = unique(cr_fits$stratum),
    year = future_years
  ) |>
    left_join(cr_fits, by = "stratum") |>
    mutate(
      t = year - 2000,
      # Chapman Richards prediction
      agb_cr = Amax * (1 - exp(-k * t))^p,
      # Set to observed value if within historical range
      agb_cr = ifelse(is.na(agb_cr), NA_real_, agb_cr)
    ) |>
    select(stratum, year, agb_cr)

  return(projections)
}

# Generate projections
projections_cr <- project_growth_curves(cr_fits, stratum_means, projection_years = 20)

# Fit Schumacher model as alternative (using synthetic ages)
fit_schumacher <- function(traj) {
  # Estimate stand age from Chapman Richards fit
  # For simplicity, use pixels in each stratum and fit: ln(AGB) = b0 + b1/age
  traj <- traj |>
    filter(agb_Mg_ha > 0) |>   # log(0) is undefined; exclude zero AGB pixels
    mutate(t = year - 2000, age_est = t + 50)

  schumacher_fits <- traj |>
    group_by(stratum) |>
    group_map(~ {
      d <- .x |> filter(is.finite(log(agb_Mg_ha)))
      if (nrow(d) < 3) {
        return(tibble(
          stratum = .y$stratum,
          b0 = NA_real_, b1 = NA_real_,
          rmse = NA_real_, r2 = NA_real_
        ))
      }
      tryCatch({
        fit <- lm(log(agb_Mg_ha) ~ I(1 / age_est), data = d)
        tibble(
          stratum = .y$stratum,
          b0 = coef(fit)[1],
          b1 = coef(fit)[2],
          rmse = sqrt(mean(residuals(fit)^2)),
          r2 = summary(fit)$r.squared
        )
      }, error = function(e) {
        tibble(
          stratum = .y$stratum,
          b0 = NA_real_, b1 = NA_real_,
          rmse = NA_real_, r2 = NA_real_
        )
      })
    }) |>
    bind_rows()

  return(schumacher_fits)
}

schumacher_fits <- fit_schumacher(traj_sample)
cat("\n--- Schumacher fits by stratum ---\n")
print(schumacher_fits)

# Project using Schumacher
projections_schumacher <- crossing(
  stratum = unique(schumacher_fits$stratum),
  year = seq(min(traj$year), max(traj$year) + 20, by = 1)
) |>
  left_join(schumacher_fits, by = "stratum") |>
  mutate(
    t = year - 2000,
    age_est = t + 50,
    agb_schumacher = exp(b0 + b1 / age_est),
    agb_schumacher = ifelse(is.na(agb_schumacher), NA_real_, agb_schumacher)
  ) |>
  select(stratum, year, agb_schumacher)

# Combine observations and projections for visualization
traj_for_plot <- traj |>
  select(stratum, year, agb_Mg_ha) |>
  group_by(stratum, year) |>
  summarise(agb_mean = mean(agb_Mg_ha, na.rm = TRUE),
            agb_sd = sd(agb_Mg_ha, na.rm = TRUE),
            .groups = "drop") |>
  mutate(data_type = "observed")

projections_combined <- projections_cr |>
  left_join(projections_schumacher, by = c("stratum", "year")) |>
  filter(year >= min(traj$year)) |>
  mutate(
    agb_mean = agb_cr,
    agb_sd = NA_real_,
    data_type = "projected_cr"
  ) |>
  select(stratum, year, agb_mean, agb_sd, data_type)

growth_data <- bind_rows(traj_for_plot, projections_combined)

# Save trajectory data for keynote figure generation (master script Part 4)
rds_path <- file.path(Sys.getenv("OUTPUT_DIR", unset = out_dir),
                       "agb_trajectories.rds")
saveRDS(traj, rds_path)
cat("Saved trajectory data to:", rds_path, "\n")

cat("\nGrowth projections to", max(projections_combined$year),
    "computed for all strata.\n")

# =============================================================================
# PART 7: Validate against field plots
# =============================================================================

# Read the NSVB field inventory (primary data source)
nsvb_plots <- read_xlsx(file.path(data_dir, "ForestInventoryDataNSVB_2026_02_23.xlsx"),
                          sheet = 1) |>
  filter(Site == "Howland") |>
  rename(nsvb_agb_kgC_m2 = "AGB_2025_kgC/m2") |>
  mutate(
    plot_id = paste0("HL_", sprintf("%03d", Plot)),
    # NSVB (Westfall et al. 2024) outputs dry biomass (kg), NOT carbon.
    # Convert kg_biomass/m2 -> Mg_biomass/ha: * 10000 / 1000 = * 10
    # (No carbon-to-biomass factor of 2 needed; contrast with Chojnacky which IS carbon)
    nsvb_agb_Mg_ha = nsvb_agb_kgC_m2 * 10
  ) |>
  select(plot_id, northing = Northing, easting = Easting,
         radius_m = Radius_m, area_m2 = Area_m2,
         nsvb_agb_kgC_m2, nsvb_agb_Mg_ha)

# Read the local field data to get all 4 allometric estimates
local_plots <- read_csv(file.path(data_dir, "Plots_Data.csv"),
                         show_col_types = FALSE) |>
  rename(radius_m = Radus_m) |>
  filter(str_detect(Plot, "^HL")) |>
  mutate(
    plot_id = Plot,
    jenkins_agb_kgC_m2 = `AGB_2025_J_KgC/m2`,
    young_agb_kgC_m2 = `AGB_2025_Y_KgC/m2`,
    crm_agb_kgC_m2 = `AGB_2025_C_KgC/m2`,
    jenkins_agb_Mg_ha = jenkins_agb_kgC_m2 * 20,
    young_agb_Mg_ha = young_agb_kgC_m2 * 20,
    crm_agb_Mg_ha = crm_agb_kgC_m2 * 20
  ) |>
  select(plot_id, northing = Northing, easting = Easting, radius_m,
         area_m2 = Area_m2, jenkins_agb_kgC_m2, young_agb_kgC_m2,
         crm_agb_kgC_m2, jenkins_agb_Mg_ha, young_agb_Mg_ha, crm_agb_Mg_ha)

# Merge NSVB with local data
field_validation <- local_plots |>
  left_join(nsvb_plots |> select(plot_id, nsvb_agb_kgC_m2, nsvb_agb_Mg_ha),
            by = "plot_id")

cat("\n--- Field plot AGB summary (Howland, n =", nrow(field_validation), "plots) ---\n")
field_validation |>
  summarise(
    jenkins_mean = round(mean(jenkins_agb_kgC_m2, na.rm = TRUE), 2),
    jenkins_sd = round(sd(jenkins_agb_kgC_m2, na.rm = TRUE), 2),
    young_mean = round(mean(young_agb_kgC_m2, na.rm = TRUE), 2),
    crm_mean = round(mean(crm_agb_kgC_m2, na.rm = TRUE), 2),
    nsvb_mean = round(mean(nsvb_agb_kgC_m2, na.rm = TRUE), 2)
  ) |>
  print()

# =============================================================================
# PART 8: 2025 AGB raster validation placeholder
# =============================================================================

cat("\n========== PART 8: 2025 AGB RASTER VALIDATION PLACEHOLDER ==========\n")

# Look for 2025 AGB raster file
raster_2025_file <- file.path(raster_dir, "Howland_AGB_2025.tif")

if (file.exists(raster_2025_file)) {
  cat("Found 2025 AGB raster. Loading and validating against field inventory...\n")

  # Load the 2025 raster
  agb_2025_raster <- rast(raster_2025_file)
  cat("Raster loaded successfully.\n")

  # Extract AGB values at field plot locations using NSVB coordinates
  if (nrow(field_validation) > 0 & all(!is.na(field_validation$easting))) {
    plot_coords <- field_validation |>
      select(plot_id, x = easting, y = northing) |>
      as.data.frame()

    # Extract raster values at plot coordinates
    extracted_agb_2025 <- terra::extract(agb_2025_raster, plot_coords[, c("x", "y")])

    # Combine with field data
    validation_2025 <- field_validation |>
      mutate(
        raster_agb_kgC_m2 = extracted_agb_2025[, 1],
        raster_agb_Mg_ha = raster_agb_kgC_m2 * 20,
        residual_vs_nsvb = raster_agb_kgC_m2 - nsvb_agb_kgC_m2
      )

    # Compute validation metrics
    rmse_vs_nsvb <- sqrt(mean((validation_2025$raster_agb_kgC_m2 -
                                 validation_2025$nsvb_agb_kgC_m2)^2, na.rm = TRUE))
    bias_vs_nsvb <- mean(validation_2025$raster_agb_kgC_m2 -
                           validation_2025$nsvb_agb_kgC_m2, na.rm = TRUE)
    r2_vs_nsvb <- cor(validation_2025$raster_agb_kgC_m2,
                       validation_2025$nsvb_agb_kgC_m2, use = "complete.obs")^2

    cat("\n--- 2025 AGB Raster vs. NSVB Field Inventory ---\n")
    cat("Number of plots:", nrow(validation_2025), "\n")
    cat("RMSE (kgC/m2):", round(rmse_vs_nsvb, 3), "\n")
    cat("Bias (kgC/m2):", round(bias_vs_nsvb, 3), "\n")
    cat("R squared:", round(r2_vs_nsvb, 3), "\n")

  } else {
    cat("Warning: Field plot coordinates unavailable for extraction.\n")
  }

} else {
  cat("2025 AGB raster not yet available; awaiting from Wei.\n")
  cat("Expected file:", raster_2025_file, "\n")
  cat("As soon as Howland_AGB_2025_10m.tif is placed in the raster directory,\n")
  cat("this validation will automatically execute.\n")
}

# =============================================================================
# PART 9: Figures
# =============================================================================

# --- Figure 1: AGB trajectories by stratum with growth projections ----
p1 <- ggplot() +
  # Historical observations
  geom_line(
    data = stratum_means,
    aes(x = year, y = mean_agb, color = stratum, linetype = "Observed"),
    linewidth = 1.1
  ) +
  geom_point(
    data = stratum_means,
    aes(x = year, y = mean_agb, color = stratum),
    size = 2.5
  ) +
  geom_ribbon(
    data = stratum_means,
    aes(x = year, ymin = mean_agb - 1.96 * se_agb, ymax = mean_agb + 1.96 * se_agb,
        fill = stratum),
    alpha = 0.2
  ) +
  # Projections (Chapman Richards)
  geom_line(
    data = projections_cr |> filter(!is.na(agb_cr)),
    aes(x = year, y = agb_cr, color = stratum, linetype = "Projected CR"),
    linewidth = 1, alpha = 0.8
  ) +
  scale_color_manual(values = pal_strata, name = NULL) +
  scale_fill_manual(values = pal_strata, guide = "none") +
  scale_linetype_manual(values = c("Observed" = "solid", "Projected CR" = "dashed"),
                        name = NULL) +
  coord_cartesian(xlim = c(2010, 2047)) +
  labs(
    x = "Year",
    y = "Aboveground biomass (Mg ha<sup>minus1</sup>)",
    title = "A) LiDAR AGB trajectories and 20 year growth projections"
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown(),
        legend.position = "inside", legend.position.inside = c(0.15, 0.95))

# --- Figure 2: PAI by period and stratum ----
p2 <- ggplot(pai_summary, aes(x = period, y = mean_pai, fill = stratum)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_pai - sd_pai / sqrt(n), ymax = mean_pai + sd_pai / sqrt(n)),
    position = position_dodge(width = 0.7), width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = pal_strata) +
  labs(
    x = "LiDAR acquisition period",
    y = "Periodic annual increment (Mg ha<sup>minus1</sup> yr<sup>minus1</sup>)",
    title = "B) AGB periodic annual increment from repeated LiDAR",
    fill = NULL
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

# --- Figure 3: Field allometric comparison (4 systems) ----
# Create long format for comparison
field_long <- field_validation |>
  select(plot_id, jenkins_agb_kgC_m2, young_agb_kgC_m2,
         crm_agb_kgC_m2, nsvb_agb_kgC_m2) |>
  pivot_longer(
    cols = ends_with("_kgC_m2"),
    names_to = "allometric_system",
    values_to = "agb_kgC_m2"
  ) |>
  mutate(
    allometric_system = case_when(
      allometric_system == "jenkins_agb_kgC_m2" ~ "Jenkins",
      allometric_system == "young_agb_kgC_m2" ~ "Young",
      allometric_system == "crm_agb_kgC_m2" ~ "CRM",
      allometric_system == "nsvb_agb_kgC_m2" ~ "NSVB (preferred)",
      TRUE ~ allometric_system
    ),
    allometric_system = factor(allometric_system,
                                levels = c("Jenkins", "Young", "CRM", "NSVB (preferred)"))
  ) |>
  drop_na(agb_kgC_m2)

p3 <- ggplot(field_long, aes(x = allometric_system, y = agb_kgC_m2, fill = allometric_system)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  scale_fill_manual(
    values = c("Jenkins" = "#E69F00", "Young" = "#56B4E9",
               "CRM" = "#009E73", "NSVB (preferred)" = "#CC79A7"),
    guide = "none"
  ) +
  labs(
    x = "Allometric system",
    y = "Plot level AGB (kgC m<sup>minus2</sup>)",
    title = "C) Comparison of four allometric estimates"
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown())

# --- Figure 4: Chapman Richards vs Schumacher models ----
# Create comparison data
model_comp_data <- tibble(
  year = seq(2012, 2045, by = 1)
) |>
  expand_grid(stratum = unique(cr_fits$stratum)) |>
  left_join(projections_cr, by = c("stratum", "year")) |>
  left_join(projections_schumacher, by = c("stratum", "year"))

p4 <- ggplot() +
  # Chapman Richards
  geom_line(
    data = model_comp_data |> filter(!is.na(agb_cr)),
    aes(x = year, y = agb_cr, color = stratum, linetype = "Chapman Richards"),
    linewidth = 1.1
  ) +
  # Schumacher
  geom_line(
    data = model_comp_data |> filter(!is.na(agb_schumacher)),
    aes(x = year, y = agb_schumacher, color = stratum, linetype = "Schumacher"),
    linewidth = 1.1, alpha = 0.7
  ) +
  # Observed means
  geom_point(
    data = stratum_means,
    aes(x = year, y = mean_agb, color = stratum),
    size = 2
  ) +
  scale_color_manual(values = pal_strata, name = NULL) +
  scale_linetype_manual(values = c("Chapman Richards" = "solid", "Schumacher" = "dashed"),
                        name = "Model") +
  coord_cartesian(xlim = c(2010, 2047)) +
  labs(
    x = "Year",
    y = "Aboveground biomass (Mg ha<sup>minus1</sup>)",
    title = "D) Yield model comparison: Chapman Richards vs. Schumacher"
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown())

# --- Combined multi panel figure ----
combined <- (p1 / p2) | (p3 / p4)
combined <- combined +
  plot_annotation(
    title = "Howland Research Forest: AGB Yield Curves from Repeated LiDAR",
    subtitle = "Data: ORNL DAAC Wei et al. (2025), G LiHT, USGS 3DEP, Phoenix RANGER U580",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(out_dir, "howland_agb_yield_curves.png"), combined,
       width = 35, height = 28, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "howland_agb_yield_curves.pdf"), combined,
       width = 35, height = 28, units = "cm")

cat("\nFigures saved to:", out_dir, "\n")

# =============================================================================
# PART 10: Summary statistics
# =============================================================================

cat("\n========== SUMMARY FOR KEYNOTE ==========\n")
cat("Howland Research Forest, central Maine\n")
cat("Repeated LiDAR acquisitions:", paste(lidar_years, collapse = ", "), "\n")
cat("Total pixels tracked:", n_distinct(traj$cell_id), "\n")
cat("Observation period:", min(traj$year), "to", max(traj$year),
    "(", max(traj$year) - min(traj$year), "years)\n")
cat("\nMean AGB by stratum (most recent year):\n")
traj |>
  filter(year == max(year)) |>
  group_by(stratum) |>
  summarise(
    mean_Mg_ha = round(mean(agb_Mg_ha, na.rm = TRUE), 1),
    sd_Mg_ha = round(sd(agb_Mg_ha, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

cat("\nField calibration plots: n =", nrow(field_validation), "\n")
cat("NSVB AGB range:", round(range(field_validation$nsvb_agb_kgC_m2, na.rm = TRUE), 2),
    "kgC/m2\n")

cat("\nChapman Richards fitted parameters:\n")
print(cr_fits)

cat("\nSchumacher fitted parameters:\n")
print(schumacher_fits)

cat("\nGrowth projection (Chapman Richards) summary:\n")
cat("Projection horizon: 2025 to 2045 (20 years)\n")
projections_cr |>
  filter(year %in% c(2025, 2035, 2045)) |>
  arrange(stratum, year) |>
  print()

cat("\n========== END OF ANALYSIS ==========\n")
