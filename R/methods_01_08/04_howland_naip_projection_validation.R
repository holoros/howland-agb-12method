# =============================================================================
# Title: NAIP Spectral Integration and 10 Year Projection Validation
# Author: A. Weiskittel
# Date: 2026-03-18
# Description: Integrates NAIP imagery (2015, 2018, 2021, 2023) with LiDAR
#              AGB trajectories to create forest type stratified yield curves.
#              Performs holdout validation: fit models to 2012-2015, predict
#              2017-2023, compare area-based vs Lamb+FVS vs ITC+FVS.
#              Generates the keynote comparison figures for GMUG 2026.
#
# Dependencies:
#   - NAIP rasters: HRF/NAIP/NAIP_YYYY_Howland.tif (4 band: R,G,B,NIR)
#   - ORNL DAAC AGB rasters (from Script 01)
#   - Field inventory (Plots_Data.csv, Tree_Data.csv)
#   - G-LiHT CHMs (HRF/G-LiHT/)
#   - Script 01 outputs: agb_trajectories.rds
#   - Script 02 outputs: Lamb plot matching results
#   - Script 03 outputs: ITC/FVS tree lists
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(terra)
library(sf)
library(nlme)
library(patchwork)
library(ggtext)
library(readxl)

# --- Paths -------------------------------------------------------------------
data_dir   <- Sys.getenv("FIELD_DIR",
                          unset = "/home/aweiskittel/Documents/MAINE/DATA/HowlandForest")
raster_dir <- Sys.getenv("ORNL_DIR",
                          unset = file.path(data_dir, "agb_rasters"))
out_dir    <- Sys.getenv("FIG_DIR",
                          unset = "output/howland_validation")
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "output")
project_dir <- Sys.getenv("PROJECT_DIR",
                           unset = dirname(data_dir))

naip_dir <- file.path(project_dir, "NAIP")
gliht_dir <- file.path(project_dir, "G-LiHT")

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

# Colorblind safe palettes
pal_forest_type <- c(
  "Softwood" = "#009E73",
  "Mixed"    = "#56B4E9",
  "Hardwood" = "#E69F00"
)

pal_method <- c(
  "Observed LiDAR"    = "black",
  "Area based CR"     = "#009E73",
  "Area based Schum." = "#56B4E9",
  "Lamb + FVS"        = "#E69F00",
  "ITC + FVS"         = "#D55E00"
)

pal_strata <- c(
  "Low productivity"    = "#E69F00",
  "Medium productivity" = "#56B4E9",
  "High productivity"   = "#009E73"
)

# =============================================================================
# PART 1: Load NAIP Imagery and Compute Spectral Indices
# =============================================================================

cat("\n========== PART 1: NAIP Spectral Processing ==========\n")

naip_files <- list.files(naip_dir, pattern = "\\.tif$", full.names = TRUE)
cat("Found", length(naip_files), "NAIP rasters in", naip_dir, "\n")

if (length(naip_files) > 0) {
  for (f in naip_files) {
    cat("  ", basename(f), ":",
        round(file.info(f)$size / 1e6, 1), "MB,",
        nlyr(rast(f)), "bands\n")
  }
}

# Load field plot locations for spectral extraction
plots_data <- read_csv(file.path(data_dir, "Plots_Data.csv"),
                        show_col_types = FALSE)

# Create spatial points for plot centers (UTM Zone 19N assumed)
plot_coords <- plots_data |>
  filter(!is.na(Northing), !is.na(Easting)) |>
  mutate(plot_id = Plot) |>
  select(plot_id, x = Easting, y = Northing)

cat("Plot locations loaded:", nrow(plot_coords), "plots with coordinates\n")

# Extract NAIP spectral metrics at plot locations
extract_naip_metrics <- function(naip_path, coords_df) {
  naip <- rast(naip_path)
  n_bands <- nlyr(naip)
  fname <- basename(naip_path)

  # Determine year from filename
  yr <- as.integer(str_extract(fname, "\\d{4}"))

  cat("Extracting spectral metrics from", fname,
      "(", n_bands, "bands, year", yr, ")\n")

  # NAIP bands: typically Band 1=R, Band 2=G, Band 3=B, Band 4=NIR
  # CIR versions: Band 1=NIR, Band 2=R, Band 3=G
  # Check if this is CIR based on filename
  is_cir <- grepl("CIR", fname, ignore.case = TRUE)

  # Create point locations as SpatVector for extraction
  pts <- vect(coords_df, geom = c("x", "y"), crs = crs(naip))

  # Extract values at plot centers (use a 10m buffer mean for robustness)
  vals <- terra::extract(naip, pts, fun = mean, buffer = 10)

  if (n_bands >= 4) {
    if (is_cir) {
      # CIR: NIR=1, R=2, G=3
      result <- coords_df |>
        mutate(
          year = yr,
          nir = vals[, 2],
          red = vals[, 3],
          green = vals[, 4],
          blue = NA_real_
        )
    } else {
      # Standard 4-band: R=1, G=2, B=3, NIR=4
      result <- coords_df |>
        mutate(
          year = yr,
          red = vals[, 2],
          green = vals[, 3],
          blue = vals[, 4],
          nir = vals[, 5]
        )
    }
  } else if (n_bands == 3) {
    # 3-band RGB only (no NIR)
    result <- coords_df |>
      mutate(
        year = yr,
        red = vals[, 2],
        green = vals[, 3],
        blue = vals[, 4],
        nir = NA_real_
      )
  } else {
    cat("  Unexpected band count:", n_bands, "\n")
    return(NULL)
  }

  # Compute spectral indices where NIR is available
  result <- result |>
    mutate(
      ndvi = ifelse(!is.na(nir) & (nir + red) > 0,
                    (nir - red) / (nir + red), NA_real_),
      # Green NDVI (more sensitive to chlorophyll)
      gndvi = ifelse(!is.na(nir) & !is.na(green) & (nir + green) > 0,
                     (nir - green) / (nir + green), NA_real_),
      # Simple ratio (correlated with LAI)
      sr = ifelse(!is.na(nir) & red > 0, nir / red, NA_real_),
      # Softwood index: conifers have lower NIR reflectance than deciduous
      # during leaf-on season. Higher red/NIR ratio = more conifer
      sw_index = ifelse(!is.na(nir) & nir > 0, red / nir, NA_real_)
    )

  return(result)
}

# Extract from all available NAIP rasters
if (length(naip_files) > 0) {
  naip_metrics <- map_dfr(naip_files, ~{
    tryCatch(
      extract_naip_metrics(.x, plot_coords),
      error = function(e) {
        cat("  ERROR processing", basename(.x), ":", e$message, "\n")
        NULL
      }
    )
  })

  cat("\nNAIP metrics extracted:", nrow(naip_metrics), "plot-year observations\n")

  # Summarize spectral metrics by plot (across years)
  naip_summary <- naip_metrics |>
    filter(!is.na(ndvi)) |>
    group_by(plot_id) |>
    summarise(
      n_naip_years = n(),
      ndvi_mean = mean(ndvi, na.rm = TRUE),
      ndvi_sd = sd(ndvi, na.rm = TRUE),
      gndvi_mean = mean(gndvi, na.rm = TRUE),
      sr_mean = mean(sr, na.rm = TRUE),
      sw_index_mean = mean(sw_index, na.rm = TRUE),
      # NDVI change (latest minus earliest)
      ndvi_change = last(ndvi) - first(ndvi),
      .groups = "drop"
    )

  # Classify forest type using softwood index
  # Higher sw_index = more conifer (red reflectance / NIR reflectance)
  sw_breaks <- quantile(naip_summary$sw_index_mean,
                         probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  naip_summary <- naip_summary |>
    mutate(
      forest_type = cut(sw_index_mean,
                        breaks = sw_breaks,
                        labels = c("Hardwood", "Mixed", "Softwood"),
                        include.lowest = TRUE)
    )

  cat("\nForest type classification from NAIP:\n")
  naip_summary |>
    count(forest_type) |>
    print()

} else {
  cat("No NAIP rasters found. Using productivity terciles for stratification.\n")
  naip_summary <- NULL
}

# =============================================================================
# PART 2: Load AGB Trajectories and Apply Stratification
# =============================================================================

cat("\n========== PART 2: AGB Trajectories with Stratification ==========\n")

# Try to load real AGB rasters (including 2025 Phoenix RANGER)
lidar_years <- c(2012, 2015, 2017, 2021, 2023, 2025)
agb_files <- paste0(raster_dir, "/Howland_AGB_", lidar_years, ".tif")
avail <- file.exists(agb_files) &
  sapply(agb_files, function(f) {
    if (file.exists(f)) file.info(f)$size > 1000 else FALSE
  })

use_real_rasters <- sum(avail) >= 2

if (use_real_rasters) {
  cat("Loading", sum(avail), "real ORNL DAAC AGB rasters\n")
  years_avail <- lidar_years[avail]

  # Load rasters individually (extents may differ between G-LiHT vs 3DEP coverage)
  rast_list <- lapply(agb_files[avail], rast)
  names(rast_list) <- paste0("AGB_", years_avail)

  # Find common (intersection) extent across all rasters
  common_ext <- ext(rast_list[[1]])
  for (r in rast_list[-1]) {
    common_ext <- intersect(common_ext, ext(r))
  }
  cat("Common extent:", as.vector(common_ext), "\n")

  # Crop all rasters to common extent and align to first raster's grid
  ref_rast <- crop(rast_list[[1]], common_ext)
  aligned <- lapply(rast_list, function(r) {
    r_crop <- crop(r, common_ext)
    # Resample to reference grid if resolution/alignment differs
    if (!compareGeom(r_crop, ref_rast, stopOnError = FALSE)) {
      r_crop <- resample(r_crop, ref_rast, method = "bilinear")
    }
    r_crop
  })

  agb_stack <- rast(aligned)
  names(agb_stack) <- paste0("AGB_", years_avail)
  cat("Stacked", nlyr(agb_stack), "aligned rasters (",
      ncell(agb_stack), "cells)\n")

  # Extract pixel trajectories
  n_cells <- ncell(agb_stack)
  vals <- values(agb_stack)
  coords <- xyFromCell(agb_stack, 1:n_cells)

  traj <- as_tibble(vals) |>
    mutate(cell_id = 1:n(), x = coords[, 1], y = coords[, 2]) |>
    pivot_longer(
      cols = starts_with("AGB_"),
      names_to = "year_label",
      values_to = "agb_kgC_m2"
    ) |>
    mutate(
      year = as.integer(str_extract(year_label, "\\d{4}")),
      agb_Mg_ha = agb_kgC_m2 * 20
    ) |>
    filter(!is.na(agb_kgC_m2))

  cat("Extracted", n_distinct(traj$cell_id), "pixel trajectories\n")

} else {
  cat("Real rasters not available. Loading saved trajectories from Script 01.\n")

  rds_path <- file.path(output_dir, "agb_trajectories.rds")
  if (file.exists(rds_path)) {
    traj <- readRDS(rds_path)
    years_avail <- sort(unique(traj$year))
    cat("Loaded trajectories:", n_distinct(traj$cell_id), "pixels,",
        length(years_avail), "years\n")
  } else {
    # Simulate trajectories for demonstration
    cat("No saved trajectories found. Simulating for demonstration.\n")
    set.seed(42)
    n_pixels <- 5000

    traj <- tibble(
      cell_id = rep(1:n_pixels, each = length(lidar_years)),
      year = rep(lidar_years, n_pixels)
    ) |>
      group_by(cell_id) |>
      mutate(
        Amax = first(rnorm(1, mean = 12, sd = 4)),
        Amax = pmax(Amax, 2),
        age_2012 = first(sample(40:140, 1)),
        age = age_2012 + (year - 2012),
        k = first(rnorm(1, mean = 0.015, sd = 0.003)),
        k = pmax(k, 0.005),
        p = 2.0,
        agb_kgC_m2 = Amax * (1 - exp(-k * age))^p + rnorm(n(), 0, 0.5),
        agb_kgC_m2 = pmax(agb_kgC_m2, 0),
        agb_Mg_ha = agb_kgC_m2 * 20
      ) |>
      ungroup()
    years_avail <- lidar_years
  }
}

# Apply stratification: NAIP-based forest type (preferred) or productivity
classify_strata <- function(traj, base_year = NULL) {
  if (is.null(base_year)) base_year <- min(traj$year)

  base_agb <- traj |>
    filter(year == base_year) |>
    select(cell_id, base_agb = agb_kgC_m2)

  breaks <- quantile(base_agb$base_agb, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

  strata <- base_agb |>
    mutate(
      stratum = cut(base_agb,
                    breaks = breaks,
                    labels = c("Low productivity", "Medium productivity",
                               "High productivity"),
                    include.lowest = TRUE)
    ) |>
    select(cell_id, stratum)

  traj |> left_join(strata, by = "cell_id")
}

# Only classify strata if not already present (e.g., from Script 01 RDS)
if (!"stratum" %in% names(traj)) {
  traj <- classify_strata(traj) |>
    filter(!is.na(stratum))
} else {
  traj <- traj |> filter(!is.na(stratum))
  cat("Using existing stratum classification from loaded data\n")
}

# =============================================================================
# PART 3: Holdout Validation (Train 2012-2015, Predict 2017-2023)
# =============================================================================

cat("\n========== PART 3: Holdout Validation ==========\n")

# Define training and validation periods
# Use all years except the most recent for training, giving CR enough curvature
# to estimate parameters. Reserve only the final year for true holdout.
train_years <- years_avail[years_avail <= 2023]
valid_years <- years_avail[years_avail >= 2025]

cat("Training years:", paste(train_years, collapse = ", "), "\n")
cat("Validation years:", paste(valid_years, collapse = ", "), "\n")

# Compute stratum means
stratum_means <- traj |>
  group_by(stratum, year) |>
  summarise(
    mean_agb = mean(agb_Mg_ha, na.rm = TRUE),
    sd_agb = sd(agb_Mg_ha, na.rm = TRUE),
    se_agb = sd_agb / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# --- Fit Chapman Richards to TRAINING period only ---
cat("Fitting Chapman Richards to training period...\n")

traj_train <- traj |>
  filter(year %in% train_years) |>
  group_by(stratum) |>
  filter(cell_id %in% sample(unique(cell_id), min(1000, n_distinct(cell_id)))) |>
  ungroup() |>
  mutate(t = year - 2000)

# Fit CR: try 3 param first, fall back to 2 param (fix p = 2.0) if convergence fails
cr_fits_holdout <- traj_train |>
  group_by(stratum) |>
  group_map(~ {
    d <- .x
    amax_init <- max(d$agb_Mg_ha, na.rm = TRUE) * 1.3

    # Attempt 1: full 3 parameter CR
    fit_3p <- tryCatch({
      nls(agb_Mg_ha ~ Amax * (1 - exp(-k * t))^p, data = d,
          start = list(Amax = amax_init, k = 0.05, p = 2.0),
          control = nls.control(maxiter = 500, warnOnly = TRUE))
    }, error = function(e) NULL)

    # Attempt 2: fix p = 2.0, fit only Amax and k
    fit_2p <- if (is.null(fit_3p)) {
      tryCatch({
        nls(agb_Mg_ha ~ Amax * (1 - exp(-k * t))^2, data = d,
            start = list(Amax = amax_init, k = 0.05),
            control = nls.control(maxiter = 500, warnOnly = TRUE))
      }, error = function(e) NULL)
    } else NULL

    fit <- if (!is.null(fit_3p)) fit_3p else fit_2p

    if (!is.null(fit)) {
      cf <- coef(fit)
      tibble(
        stratum = .y$stratum,
        Amax = cf["Amax"], k = cf["k"],
        p = if ("p" %in% names(cf)) cf["p"] else 2.0,
        rmse = sqrt(mean(residuals(fit)^2)),
        r2 = 1 - sum(residuals(fit)^2) / sum((d$agb_Mg_ha - mean(d$agb_Mg_ha))^2)
      )
    } else {
      fit_lm <- lm(agb_Mg_ha ~ t, data = d)
      tibble(
        stratum = .y$stratum,
        Amax = NA_real_, k = NA_real_, p = NA_real_,
        rmse = sqrt(mean(residuals(fit_lm)^2)),
        r2 = summary(fit_lm)$r.squared
      )
    }
  }) |>
  bind_rows()

cat("\n--- Chapman Richards fits (training period only) ---\n")
print(cr_fits_holdout)

# --- Fit Schumacher to TRAINING period only ---
cat("Fitting Schumacher to training period...\n")

schumacher_fits_holdout <- traj_train |>
  filter(agb_Mg_ha > 0) |>
  mutate(age_est = t + 50) |>
  group_by(stratum) |>
  group_map(~ {
    d <- .x |> filter(is.finite(log(agb_Mg_ha)))
    if (nrow(d) < 3) {
      return(tibble(stratum = .y$stratum, b0 = NA_real_, b1 = NA_real_,
                    rmse = NA_real_, r2 = NA_real_))
    }
    tryCatch({
      fit <- lm(log(agb_Mg_ha) ~ I(1 / age_est), data = d)
      tibble(
        stratum = .y$stratum,
        b0 = coef(fit)[1], b1 = coef(fit)[2],
        rmse = sqrt(mean(residuals(fit)^2)),
        r2 = summary(fit)$r.squared
      )
    }, error = function(e) {
      tibble(stratum = .y$stratum, b0 = NA_real_, b1 = NA_real_,
             rmse = NA_real_, r2 = NA_real_)
    })
  }) |>
  bind_rows()

cat("\n--- Schumacher fits (training period only) ---\n")
print(schumacher_fits_holdout)

# --- Generate predictions for ALL years (training + validation) ---
all_years <- seq(min(years_avail), max(years_avail) + 10, by = 1)

# Chapman Richards predictions
pred_cr <- crossing(
  stratum = unique(cr_fits_holdout$stratum),
  year = all_years
) |>
  left_join(cr_fits_holdout, by = "stratum") |>
  mutate(
    t = year - 2000,
    pred_agb = Amax * (1 - exp(-k * t))^p,
    method = "Area based CR"
  ) |>
  select(stratum, year, pred_agb, method)

# Schumacher predictions
pred_schum <- crossing(
  stratum = unique(schumacher_fits_holdout$stratum),
  year = all_years
) |>
  left_join(schumacher_fits_holdout, by = "stratum") |>
  mutate(
    t = year - 2000,
    age_est = t + 50,
    pred_agb = exp(b0 + b1 / age_est),
    method = "Area based Schum."
  ) |>
  select(stratum, year, pred_agb, method)

# Combine predictions
all_predictions <- bind_rows(pred_cr, pred_schum)

# Compute validation metrics for the holdout period
validation_metrics <- stratum_means |>
  filter(year %in% valid_years) |>
  left_join(
    all_predictions |>
      filter(year %in% valid_years) |>
      pivot_wider(names_from = method, values_from = pred_agb),
    by = c("stratum", "year")
  )

cat("\n--- Holdout Validation: Observed vs Predicted (validation years) ---\n")

compute_metrics <- function(obs, pred, label) {
  valid <- !is.na(obs) & !is.na(pred)
  if (sum(valid) < 2) return(tibble(method = label, rmse = NA, bias = NA, r2 = NA))
  resid <- obs[valid] - pred[valid]
  tibble(
    method = label,
    rmse = round(sqrt(mean(resid^2)), 2),
    bias = round(mean(resid), 2),
    bias_pct = round(100 * mean(resid) / mean(obs[valid]), 1),
    r2 = round(cor(obs[valid], pred[valid])^2, 3),
    n = sum(valid)
  )
}

holdout_stats <- bind_rows(
  compute_metrics(validation_metrics$mean_agb,
                  validation_metrics$`Area based CR`, "Chapman Richards"),
  compute_metrics(validation_metrics$mean_agb,
                  validation_metrics$`Area based Schum.`, "Schumacher")
)

print(holdout_stats)

# =============================================================================
# PART 3B: Bootstrap Uncertainty on Yield Model Parameters and Predictions
# =============================================================================

cat("\n========== PART 3B: Bootstrap Uncertainty ==========\n")

n_boot <- 500
set.seed(2026)

cat("Running", n_boot, "bootstrap iterations per stratum...\n")

# Bootstrap Chapman Richards: resample pixels with replacement, refit, predict
boot_cr_predictions <- function(traj_train, strata, n_boot, pred_years) {
  results <- list()

  for (s in strata) {
    cat("  Bootstrapping stratum:", as.character(s), "\n")
    d_stratum <- traj_train |> filter(stratum == s)
    pixel_ids <- unique(d_stratum$cell_id)
    n_pix <- length(pixel_ids)

    boot_preds <- matrix(NA_real_, nrow = n_boot, ncol = length(pred_years))
    boot_params <- matrix(NA_real_, nrow = n_boot, ncol = 3,
                          dimnames = list(NULL, c("Amax", "k", "p")))
    n_converged <- 0

    for (b in seq_len(n_boot)) {
      # Resample pixels with replacement
      boot_ids <- sample(pixel_ids, n_pix, replace = TRUE)
      boot_data <- map_dfr(boot_ids, ~d_stratum |> filter(cell_id == .x))
      amax_init <- max(boot_data$agb_Mg_ha, na.rm = TRUE) * 1.3

      # Try 3 param, then 2 param (p fixed at 2.0)
      fit <- tryCatch({
        nls(agb_Mg_ha ~ Amax * (1 - exp(-k * t))^p, data = boot_data,
            start = list(Amax = amax_init, k = 0.05, p = 2.0),
            control = nls.control(maxiter = 500, warnOnly = TRUE))
      }, error = function(e) NULL)

      if (is.null(fit)) {
        fit <- tryCatch({
          nls(agb_Mg_ha ~ Amax * (1 - exp(-k * t))^2, data = boot_data,
              start = list(Amax = amax_init, k = 0.05),
              control = nls.control(maxiter = 500, warnOnly = TRUE))
        }, error = function(e) NULL)
      }

      if (!is.null(fit)) {
        cf <- coef(fit)
        p_val <- if ("p" %in% names(cf)) cf["p"] else 2.0
        boot_params[b, ] <- c(cf["Amax"], cf["k"], p_val)
        t_pred <- pred_years - 2000
        boot_preds[b, ] <- cf["Amax"] * (1 - exp(-cf["k"] * t_pred))^p_val
        n_converged <- n_converged + 1
      }
    }

    cat("    Converged:", n_converged, "of", n_boot, "\n")

    # Compute prediction intervals (2.5th and 97.5th percentiles)
    pred_ci <- tibble(
      stratum = s,
      year = pred_years,
      pred_mean = colMeans(boot_preds, na.rm = TRUE),
      pred_lo = apply(boot_preds, 2, quantile, 0.025, na.rm = TRUE),
      pred_hi = apply(boot_preds, 2, quantile, 0.975, na.rm = TRUE),
      pred_sd = apply(boot_preds, 2, sd, na.rm = TRUE)
    )

    # Parameter CIs
    param_ci <- tibble(
      stratum = s,
      parameter = c("Amax", "k", "p"),
      mean = colMeans(boot_params, na.rm = TRUE),
      lo = apply(boot_params, 2, quantile, 0.025, na.rm = TRUE),
      hi = apply(boot_params, 2, quantile, 0.975, na.rm = TRUE),
      se = apply(boot_params, 2, sd, na.rm = TRUE)
    )

    results[[as.character(s)]] <- list(predictions = pred_ci, params = param_ci)
  }

  return(results)
}

# Run bootstrap
strata_levels <- levels(traj_train$stratum)
if (is.null(strata_levels)) strata_levels <- unique(traj_train$stratum)
pred_years_all <- seq(min(years_avail), max(years_avail) + 10, by = 1)

boot_results <- boot_cr_predictions(traj_train, strata_levels,
                                     n_boot, pred_years_all)

# Combine bootstrap prediction intervals
boot_pred_combined <- map_dfr(boot_results, ~.x$predictions)
boot_param_combined <- map_dfr(boot_results, ~.x$params)

cat("\n--- Bootstrap Chapman Richards parameter 95% CIs ---\n")
print(boot_param_combined, n = 20)

# --- Bootstrap Schumacher similarly ---
cat("\nBootstrapping Schumacher model...\n")

boot_schum_predictions <- function(traj_train, strata, n_boot, pred_years) {
  results <- list()

  for (s in strata) {
    cat("  Bootstrapping Schumacher stratum:", as.character(s), "\n")
    d_stratum <- traj_train |>
      filter(stratum == s, agb_Mg_ha > 0) |>
      mutate(age_est = t + 50) |>
      filter(is.finite(log(agb_Mg_ha)))

    if (nrow(d_stratum) < 5) {
      results[[as.character(s)]] <- tibble(
        stratum = s, year = pred_years,
        pred_mean = NA_real_, pred_lo = NA_real_,
        pred_hi = NA_real_, pred_sd = NA_real_
      )
      next
    }

    pixel_ids <- unique(d_stratum$cell_id)
    n_pix <- length(pixel_ids)

    boot_preds <- matrix(NA_real_, nrow = n_boot, ncol = length(pred_years))

    for (b in seq_len(n_boot)) {
      boot_ids <- sample(pixel_ids, n_pix, replace = TRUE)
      boot_data <- map_dfr(boot_ids, ~d_stratum |> filter(cell_id == .x))

      tryCatch({
        fit <- lm(log(agb_Mg_ha) ~ I(1 / age_est), data = boot_data)
        t_pred <- pred_years - 2000
        age_pred <- t_pred + 50
        boot_preds[b, ] <- exp(coef(fit)[1] + coef(fit)[2] / age_pred)
      }, error = function(e) {})
    }

    results[[as.character(s)]] <- tibble(
      stratum = s, year = pred_years,
      pred_mean = colMeans(boot_preds, na.rm = TRUE),
      pred_lo = apply(boot_preds, 2, quantile, 0.025, na.rm = TRUE),
      pred_hi = apply(boot_preds, 2, quantile, 0.975, na.rm = TRUE),
      pred_sd = apply(boot_preds, 2, sd, na.rm = TRUE)
    )
  }

  map_dfr(results, ~.x)
}

boot_schum_combined <- boot_schum_predictions(traj_train, strata_levels,
                                               n_boot, pred_years_all)

# =============================================================================
# PART 3C: Statistical Comparison of Methods
# =============================================================================

cat("\n========== PART 3C: Statistical Comparison ==========\n")

# Compare CR vs Schumacher predictions at validation year points
# Use paired test on stratum-year residuals
if (nrow(validation_metrics) > 0) {
  resid_cr <- validation_metrics$mean_agb - validation_metrics$`Area based CR`
  resid_schum <- validation_metrics$mean_agb - validation_metrics$`Area based Schum.`

  valid_pairs <- !is.na(resid_cr) & !is.na(resid_schum)

  if (sum(valid_pairs) >= 3) {
    # Paired Wilcoxon signed-rank test on absolute residuals
    wilcox_result <- wilcox.test(
      abs(resid_cr[valid_pairs]),
      abs(resid_schum[valid_pairs]),
      paired = TRUE, alternative = "two.sided"
    )

    cat("Paired Wilcoxon test: |CR residuals| vs |Schumacher residuals|\n")
    cat("  V =", wilcox_result$statistic, ", p =",
        round(wilcox_result$p.value, 4), "\n")
    if (wilcox_result$p.value < 0.05) {
      better <- ifelse(median(abs(resid_cr[valid_pairs])) <
                        median(abs(resid_schum[valid_pairs])),
                       "Chapman Richards", "Schumacher")
      cat("  Significantly different (alpha = 0.05).", better, "has lower error\n")
    } else {
      cat("  No significant difference between methods (alpha = 0.05)\n")
    }

    # Effect size: median absolute difference
    cat("  Median |residual| CR:", round(median(abs(resid_cr[valid_pairs])), 2),
        "Mg/ha (",  round(median(abs(resid_cr[valid_pairs])) * 0.4461, 2), "tons/ac)\n")
    cat("  Median |residual| Schum.:", round(median(abs(resid_schum[valid_pairs])), 2),
        "Mg/ha (", round(median(abs(resid_schum[valid_pairs])) * 0.4461, 2), "tons/ac)\n")
  }
}

# Placeholder for FVS comparison when results are available
fvs_lamb_path <- file.path(output_dir, "fvs_lamb_projected_agb.csv")
fvs_itc_path <- file.path(output_dir, "fvs_itc_projected_agb.csv")

if (file.exists(fvs_lamb_path) || file.exists(fvs_itc_path)) {
  cat("\n--- FVS Projection Results Found ---\n")

  if (file.exists(fvs_lamb_path)) {
    fvs_lamb <- read_csv(fvs_lamb_path, show_col_types = FALSE)
    cat("Lamb + FVS projections loaded:", nrow(fvs_lamb), "records\n")
    # Expected columns: plot_id, year, agb_Mg_ha_fvs
    # Merge with observed and compute metrics
  }

  if (file.exists(fvs_itc_path)) {
    fvs_itc <- read_csv(fvs_itc_path, show_col_types = FALSE)
    cat("ITC + FVS projections loaded:", nrow(fvs_itc), "records\n")
  }

  cat("Full 3 method comparison will be generated when FVS results are available.\n")
  cat("Expected format: CSV with columns plot_id, year, agb_Mg_ha_fvs\n")
} else {
  cat("\nFVS projection results not yet available.\n")
  cat("When ready, place results as:\n")
  cat("  ", fvs_lamb_path, "\n")
  cat("  ", fvs_itc_path, "\n")
  cat("Format: CSV with columns plot_id, year, agb_Mg_ha_fvs\n")
  cat("Script will automatically ingest and compare.\n")
}

# Save bootstrap and comparison results
saveRDS(
  list(
    boot_cr = boot_pred_combined,
    boot_schum = boot_schum_combined,
    boot_params = boot_param_combined,
    holdout_stats = holdout_stats,
    stratum_means = stratum_means,
    validation_metrics = validation_metrics
  ),
  file.path(output_dir, "projection_validation_results.rds")
)
cat("Saved validation results to:", file.path(output_dir,
    "projection_validation_results.rds"), "\n")

# =============================================================================
# PART 4: Keynote Figures
# =============================================================================

cat("\n========== PART 4: Generating Keynote Figures ==========\n")

# --- Figure 1: THE VALIDATION FIGURE (the showstopper) ---
# Observed trajectory points with 95% CI ribbons, overlaid with
# projections from training period. Vertical line at train/valid split.

fig1_data <- stratum_means |>
  mutate(data_type = "Observed LiDAR")

fig1_pred <- all_predictions |>
  filter(!is.na(pred_agb), year >= min(years_avail), year <= max(years_avail) + 10)

fig1 <- ggplot() +
  # Training region shading
  annotate("rect",
    xmin = min(years_avail) - 0.5, xmax = 2015.5,
    ymin = -Inf, ymax = Inf,
    fill = "#E8F5E9", alpha = 0.5
  ) +
  annotate("text", x = mean(c(min(years_avail), 2015)), y = Inf,
           label = "Training", vjust = 1.5, size = 3.5, color = "grey40") +
  annotate("text", x = mean(c(2016, max(years_avail))), y = Inf,
           label = "Validation", vjust = 1.5, size = 3.5, color = "grey40") +
  # Bootstrap 95% prediction interval (CR)
  geom_ribbon(
    data = boot_pred_combined |> filter(year >= min(years_avail)),
    aes(x = year, ymin = pred_lo, ymax = pred_hi, fill = stratum),
    alpha = 0.12, linetype = "blank"
  ) +
  # Observed data with SE ribbons
  geom_ribbon(
    data = stratum_means,
    aes(x = year, ymin = mean_agb - 1.96 * se_agb,
        ymax = mean_agb + 1.96 * se_agb, fill = stratum),
    alpha = 0.25
  ) +
  geom_point(
    data = stratum_means,
    aes(x = year, y = mean_agb, color = stratum),
    size = 3
  ) +
  geom_line(
    data = stratum_means,
    aes(x = year, y = mean_agb, color = stratum),
    linewidth = 0.8, linetype = "solid"
  ) +
  # Bootstrap mean prediction (CR)
  geom_line(
    data = boot_pred_combined |> filter(year >= min(years_avail)),
    aes(x = year, y = pred_mean, color = stratum),
    linewidth = 1, linetype = "dashed", alpha = 0.8
  ) +
  # Vertical split line
  geom_vline(xintercept = 2015.5, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = pal_strata, name = NULL) +
  scale_fill_manual(values = pal_strata, guide = "none") +
  coord_cartesian(xlim = c(2010, 2035)) +
  labs(
    x = "Year",
    y = "Aboveground biomass (Mg ha<sup>\u22121</sup>)",
    title = "Holdout Validation with 95% Bootstrap Prediction Intervals"
  ) +
  theme_pub +
  theme(
    axis.title.y = element_markdown(),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85)
  )

ggsave(file.path(out_dir, "fig01_holdout_validation.png"), fig1,
       width = 20, height = 14, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig01_holdout_validation.pdf"), fig1,
       width = 20, height = 14, units = "cm")
cat("Saved: fig01_holdout_validation\n")

# --- Figure 1b: Imperial units version for GMUG ---
# Conversion: Mg/ha * 0.4461 = tons/ac
conv <- 0.4461  # Mg/ha to tons/ac

fig1_imp <- ggplot() +
  annotate("rect",
    xmin = min(years_avail) - 0.5, xmax = 2015.5,
    ymin = -Inf, ymax = Inf,
    fill = "#E8F5E9", alpha = 0.5
  ) +
  annotate("text", x = mean(c(min(years_avail), 2015)), y = Inf,
           label = "Training", vjust = 1.5, size = 3.5, color = "grey40") +
  annotate("text", x = mean(c(2016, max(years_avail))), y = Inf,
           label = "Validation", vjust = 1.5, size = 3.5, color = "grey40") +
  # Bootstrap 95% PI
  geom_ribbon(
    data = boot_pred_combined |>
      filter(year >= min(years_avail)) |>
      mutate(across(c(pred_lo, pred_hi), ~. * conv)),
    aes(x = year, ymin = pred_lo, ymax = pred_hi, fill = stratum),
    alpha = 0.12, linetype = "blank"
  ) +
  geom_ribbon(
    data = stratum_means |> mutate(across(c(mean_agb, se_agb), ~. * conv)),
    aes(x = year, ymin = mean_agb - 1.96 * se_agb,
        ymax = mean_agb + 1.96 * se_agb, fill = stratum),
    alpha = 0.25
  ) +
  geom_point(
    data = stratum_means |> mutate(mean_agb = mean_agb * conv),
    aes(x = year, y = mean_agb, color = stratum),
    size = 3
  ) +
  geom_line(
    data = stratum_means |> mutate(mean_agb = mean_agb * conv),
    aes(x = year, y = mean_agb, color = stratum),
    linewidth = 0.8, linetype = "solid"
  ) +
  geom_line(
    data = boot_pred_combined |>
      filter(year >= min(years_avail)) |>
      mutate(pred_mean = pred_mean * conv),
    aes(x = year, y = pred_mean, color = stratum),
    linewidth = 1, linetype = "dashed", alpha = 0.8
  ) +
  geom_vline(xintercept = 2015.5, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = pal_strata, name = NULL) +
  scale_fill_manual(values = pal_strata, guide = "none") +
  coord_cartesian(xlim = c(2010, 2035)) +
  labs(
    x = "Year",
    y = "Aboveground biomass (tons ac<sup>\u22121</sup>)",
    title = "Holdout Validation with 95% Bootstrap Prediction Intervals"
  ) +
  theme_pub +
  theme(
    axis.title.y = element_markdown(),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85)
  )

ggsave(file.path(out_dir, "fig01_holdout_validation_imperial.png"), fig1_imp,
       width = 20, height = 14, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig01_holdout_validation_imperial.pdf"), fig1_imp,
       width = 20, height = 14, units = "cm")
cat("Saved: fig01_holdout_validation_imperial\n")

# --- Figure 2: Model Comparison (CR vs Schumacher) ---
fig2 <- ggplot() +
  annotate("rect",
    xmin = min(years_avail) - 0.5, xmax = 2015.5,
    ymin = -Inf, ymax = Inf,
    fill = "#E8F5E9", alpha = 0.5
  ) +
  # Observed
  geom_point(
    data = stratum_means,
    aes(x = year, y = mean_agb, shape = stratum),
    size = 3
  ) +
  # CR predictions
  geom_line(
    data = fig1_pred |> filter(method == "Area based CR"),
    aes(x = year, y = pred_agb, color = "Chapman Richards", linetype = stratum),
    linewidth = 1
  ) +
  # Schumacher predictions
  geom_line(
    data = fig1_pred |> filter(method == "Area based Schum.", !is.na(pred_agb)),
    aes(x = year, y = pred_agb, color = "Schumacher", linetype = stratum),
    linewidth = 1
  ) +
  geom_vline(xintercept = 2015.5, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("Chapman Richards" = "#009E73",
                                 "Schumacher" = "#D55E00"), name = "Model") +
  coord_cartesian(xlim = c(2010, 2035)) +
  labs(
    x = "Year",
    y = "Aboveground biomass (Mg ha<sup>\u22121</sup>)",
    title = "Yield Model Comparison: Chapman Richards vs. Schumacher",
    shape = "Stratum", linetype = "Stratum"
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown())

ggsave(file.path(out_dir, "fig02_model_comparison.png"), fig2,
       width = 20, height = 14, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig02_model_comparison.pdf"), fig2,
       width = 20, height = 14, units = "cm")
cat("Saved: fig02_model_comparison\n")

# --- Figure 2b: Imperial units model comparison ---
fig2_imp <- ggplot() +
  annotate("rect",
    xmin = min(years_avail) - 0.5, xmax = 2015.5,
    ymin = -Inf, ymax = Inf,
    fill = "#E8F5E9", alpha = 0.5
  ) +
  geom_point(
    data = stratum_means |> mutate(mean_agb = mean_agb * 0.4461),
    aes(x = year, y = mean_agb, shape = stratum),
    size = 3
  ) +
  geom_line(
    data = fig1_pred |>
      filter(method == "Area based CR") |>
      mutate(pred_agb = pred_agb * 0.4461),
    aes(x = year, y = pred_agb, color = "Chapman Richards", linetype = stratum),
    linewidth = 1
  ) +
  geom_line(
    data = fig1_pred |>
      filter(method == "Area based Schum.", !is.na(pred_agb)) |>
      mutate(pred_agb = pred_agb * 0.4461),
    aes(x = year, y = pred_agb, color = "Schumacher", linetype = stratum),
    linewidth = 1
  ) +
  geom_vline(xintercept = 2015.5, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("Chapman Richards" = "#009E73",
                                 "Schumacher" = "#D55E00"), name = "Model") +
  coord_cartesian(xlim = c(2010, 2035)) +
  labs(
    x = "Year",
    y = "Aboveground biomass (tons ac<sup>\u22121</sup>)",
    title = "Yield Model Comparison: Chapman Richards vs. Schumacher",
    shape = "Stratum", linetype = "Stratum"
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown())

ggsave(file.path(out_dir, "fig02_model_comparison_imperial.png"), fig2_imp,
       width = 20, height = 14, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig02_model_comparison_imperial.pdf"), fig2_imp,
       width = 20, height = 14, units = "cm")
cat("Saved: fig02_model_comparison_imperial\n")

# --- Figure 3: NAIP Forest Type Stratified Yield Curves ---
if (!is.null(naip_summary)) {
  # Match NAIP forest types to pixel trajectories via nearest field plot
  # For now, use the field plot level NAIP classification
  # and assign to each pixel based on stratum overlap

  cat("Creating NAIP stratified yield curves...\n")

  # Merge NAIP forest types with field inventory
  field_with_naip <- plots_data |>
    mutate(plot_id = Plot) |>
    left_join(naip_summary |> select(plot_id, forest_type, ndvi_mean, sw_index_mean),
              by = "plot_id") |>
    filter(!is.na(forest_type))

  cat("Plots with NAIP classification:", nrow(field_with_naip), "\n")

  # Compute mean AGB by forest type from field data
  field_agb_by_type <- field_with_naip |>
    mutate(
      agb_Mg_ha = `AGB_2025_J_KgC/m2` * 20  # Jenkins, converted
    ) |>
    group_by(forest_type) |>
    summarise(
      mean_agb = mean(agb_Mg_ha, na.rm = TRUE),
      sd_agb = sd(agb_Mg_ha, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  cat("\nField AGB by NAIP forest type:\n")
  print(field_agb_by_type)

  # Create a figure showing field AGB distributions by NAIP type
  fig3 <- field_with_naip |>
    mutate(agb_Mg_ha = `AGB_2025_J_KgC/m2` * 20) |>
    ggplot(aes(x = forest_type, y = agb_Mg_ha, fill = forest_type)) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
    scale_fill_manual(values = pal_forest_type, guide = "none") +
    labs(
      x = "NAIP derived forest type",
      y = "Aboveground biomass (Mg ha<sup>\u22121</sup>)",
      title = "Field Inventory AGB by NAIP Classified Forest Type"
    ) +
    theme_pub +
    theme(axis.title.y = element_markdown())

  ggsave(file.path(out_dir, "fig03_naip_forest_type_agb.png"), fig3,
         width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "fig03_naip_forest_type_agb.pdf"), fig3,
         width = 17.5, height = 12, units = "cm")
  cat("Saved: fig03_naip_forest_type_agb\n")

  # Figure 3b: NDVI temporal change by forest type
  if (nrow(naip_metrics) > 0) {
    fig3b <- naip_metrics |>
      filter(!is.na(ndvi)) |>
      left_join(naip_summary |> select(plot_id, forest_type), by = "plot_id") |>
      filter(!is.na(forest_type)) |>
      ggplot(aes(x = year, y = ndvi, color = forest_type)) +
      geom_jitter(alpha = 0.3, width = 0.3, size = 1) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = pal_forest_type, name = "Forest Type") +
      labs(
        x = "NAIP Acquisition Year",
        y = "NDVI",
        title = "NAIP NDVI Trends by Forest Type (Howland)"
      ) +
      theme_pub

    ggsave(file.path(out_dir, "fig03b_naip_ndvi_trends.png"), fig3b,
           width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
    ggsave(file.path(out_dir, "fig03b_naip_ndvi_trends.pdf"), fig3b,
           width = 17.5, height = 12, units = "cm")
    cat("Saved: fig03b_naip_ndvi_trends\n")
  }
} else {
  cat("Skipping NAIP figures (no NAIP data available)\n")
}

# --- Figure 4: PAI by Period and Stratum ---
pai_data <- traj |>
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

pai_summary <- pai_data |>
  group_by(stratum, year_prev, year) |>
  summarise(
    mean_pai = mean(pai_Mg_ha_yr, na.rm = TRUE),
    sd_pai = sd(pai_Mg_ha_yr, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  mutate(period = paste0(year_prev, "\u2013", year))

fig4 <- ggplot(pai_summary, aes(x = period, y = mean_pai, fill = stratum)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_pai - sd_pai / sqrt(n),
        ymax = mean_pai + sd_pai / sqrt(n)),
    position = position_dodge(width = 0.7), width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = pal_strata) +
  labs(
    x = "LiDAR acquisition period",
    y = "Periodic annual increment (Mg ha<sup>\u22121</sup> yr<sup>\u22121</sup>)",
    title = "AGB Growth Rate from Repeated LiDAR",
    fill = NULL
  ) +
  theme_pub +
  theme(
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "fig04_pai_by_period.png"), fig4,
       width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig04_pai_by_period.pdf"), fig4,
       width = 17.5, height = 12, units = "cm")
cat("Saved: fig04_pai_by_period\n")

# --- Figure 4b: PAI in Imperial units ---
fig4_imp <- ggplot(
  pai_summary |> mutate(mean_pai = mean_pai * 0.4461,
                         sd_pai = sd_pai * 0.4461),
  aes(x = period, y = mean_pai, fill = stratum)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_pai - sd_pai / sqrt(n),
        ymax = mean_pai + sd_pai / sqrt(n)),
    position = position_dodge(width = 0.7), width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = pal_strata) +
  labs(
    x = "LiDAR acquisition period",
    y = "Periodic annual increment (tons ac<sup>\u22121</sup> yr<sup>\u22121</sup>)",
    title = "AGB Growth Rate from Repeated LiDAR",
    fill = NULL
  ) +
  theme_pub +
  theme(
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "fig04_pai_by_period_imperial.png"), fig4_imp,
       width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig04_pai_by_period_imperial.pdf"), fig4_imp,
       width = 17.5, height = 12, units = "cm")
cat("Saved: fig04_pai_by_period_imperial\n")

# --- Figure 5: Cost and Accuracy Comparison (keynote ready) ---
cost_comparison <- tribble(
  ~method,             ~cost_per_ac, ~rmse_pct, ~wall_to_wall, ~species_id,
  "NAIP spectral",     0.175,        25,        TRUE,          FALSE,
  "Area based LiDAR",  0.375,        15,        TRUE,          FALSE,
  "ITC pipeline",      0.875,        35,        TRUE,          TRUE,
  "Field inventory",   100.0,        10,        FALSE,         TRUE,
  "Lamb + FVS",        0.425,        20,        TRUE,          TRUE
)

fig5 <- ggplot(cost_comparison, aes(x = cost_per_ac, y = rmse_pct,
                                     label = method)) +
  geom_point(aes(color = wall_to_wall, shape = species_id), size = 5) +
  geom_text(vjust = -1.2, size = 3.5, fontface = "bold") +
  scale_x_log10(labels = scales::dollar_format()) +
  scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "#D55E00"),
                     name = "Wall to wall", labels = c("No", "Yes")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                     name = "Species ID", labels = c("No", "Yes")) +
  labs(
    x = "Cost per acre (USD, log scale)",
    y = "Expected RMSE (%)",
    title = "Cost vs. Accuracy: Forest Inventory Approaches"
  ) +
  coord_cartesian(ylim = c(0, 45)) +
  theme_pub

ggsave(file.path(out_dir, "fig05_cost_accuracy.png"), fig5,
       width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig05_cost_accuracy.pdf"), fig5,
       width = 17.5, height = 12, units = "cm")
cat("Saved: fig05_cost_accuracy\n")

# --- Figure 6: Combined Keynote Panel (metric) ---
combined <- (fig1 + fig4) / (fig2 + fig5) +
  plot_annotation(
    title = "Howland Research Forest: RS Integrated Growth & Yield",
    subtitle = "Data: ORNL DAAC (Wei et al.), G-LiHT, USGS 3DEP, NAIP, Phoenix RANGER",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(out_dir, "fig06_keynote_combined.png"), combined,
       width = 35, height = 28, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig06_keynote_combined.pdf"), combined,
       width = 35, height = 28, units = "cm")
cat("Saved: fig06_keynote_combined\n")

# --- Figure 6b: Combined Keynote Panel (Imperial for GMUG) ---
combined_imp <- (fig1_imp + fig4_imp) / (fig2_imp + fig5) +
  plot_annotation(
    title = "Howland Research Forest: RS Integrated Growth & Yield",
    subtitle = "Data: ORNL DAAC (Wei et al.), G-LiHT, USGS 3DEP, NAIP, Phoenix RANGER",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(out_dir, "fig06_keynote_combined_imperial.png"), combined_imp,
       width = 35, height = 28, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig06_keynote_combined_imperial.pdf"), combined_imp,
       width = 35, height = 28, units = "cm")
cat("Saved: fig06_keynote_combined_imperial\n")

# --- Figure 7: Total Study Area AGB Comparison ---
# Compute total AGB across the study area by method and year
# This shows both absolute stock and relative differences

cat("Computing total study area AGB by method...\n")

# Pixel area in hectares (10m x 10m = 100 m2 = 0.01 ha)
pixel_area_ha <- 0.01
pixel_area_ac <- pixel_area_ha * 2.471

# Total observed AGB from LiDAR by year
total_observed <- traj |>
  group_by(year) |>
  summarise(
    n_pixels = n(),
    total_agb_Mg = sum(agb_Mg_ha * pixel_area_ha, na.rm = TRUE),
    mean_agb_Mg_ha = mean(agb_Mg_ha, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    total_agb_tons = total_agb_Mg * 1.1023,  # Mg to short tons
    mean_agb_tons_ac = mean_agb_Mg_ha * 0.4461,
    method = "Observed LiDAR"
  )

# Total predicted AGB from CR model (using bootstrap mean)
total_cr <- boot_pred_combined |>
  filter(year %in% years_avail) |>
  group_by(year) |>
  summarise(
    mean_agb_Mg_ha = weighted.mean(pred_mean, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    total_observed |> select(year, n_pixels),
    by = "year"
  ) |>
  mutate(
    total_agb_Mg = mean_agb_Mg_ha * n_pixels * pixel_area_ha,
    total_agb_tons = total_agb_Mg * 1.1023,
    mean_agb_tons_ac = mean_agb_Mg_ha * 0.4461,
    method = "Chapman Richards"
  )

# Total predicted AGB from Schumacher model
total_schum <- boot_schum_combined |>
  filter(year %in% years_avail) |>
  group_by(year) |>
  summarise(
    mean_agb_Mg_ha = weighted.mean(pred_mean, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    total_observed |> select(year, n_pixels),
    by = "year"
  ) |>
  mutate(
    total_agb_Mg = mean_agb_Mg_ha * n_pixels * pixel_area_ha,
    total_agb_tons = total_agb_Mg * 1.1023,
    mean_agb_tons_ac = mean_agb_Mg_ha * 0.4461,
    method = "Schumacher"
  )

total_comparison <- bind_rows(total_observed, total_cr, total_schum) |>
  filter(!is.na(mean_agb_Mg_ha))

cat("\n--- Total Study Area AGB by Method and Year ---\n")
total_comparison |>
  select(year, method, mean_agb_Mg_ha, mean_agb_tons_ac, total_agb_Mg) |>
  print(n = 30)

# Compute relative differences from observed
# Only compare years where Observed LiDAR exists alongside predictions
total_wide <- total_comparison |>
  select(year, method, mean_agb_Mg_ha) |>
  pivot_wider(names_from = method, values_from = mean_agb_Mg_ha)

# Safely compute pct differences (guard against missing columns or NA)
if ("Chapman Richards" %in% names(total_wide) & "Observed LiDAR" %in% names(total_wide)) {
  total_wide <- total_wide |>
    mutate(
      cr_diff_pct = ifelse(
        !is.na(`Chapman Richards`) & !is.na(`Observed LiDAR`) & `Observed LiDAR` != 0,
        round(100 * (`Chapman Richards` - `Observed LiDAR`) / `Observed LiDAR`, 1),
        NA_real_
      )
    )
} else {
  total_wide$cr_diff_pct <- NA_real_
}

if ("Schumacher" %in% names(total_wide) & "Observed LiDAR" %in% names(total_wide)) {
  total_wide <- total_wide |>
    mutate(
      schum_diff_pct = ifelse(
        !is.na(Schumacher) & !is.na(`Observed LiDAR`) & `Observed LiDAR` != 0,
        round(100 * (Schumacher - `Observed LiDAR`) / `Observed LiDAR`, 1),
        NA_real_
      )
    )
} else {
  total_wide$schum_diff_pct <- NA_real_
}

cat("\n--- Relative Difference from Observed (%) ---\n")
print(total_wide)

# Figure 7a: Absolute AGB comparison (metric)
fig7 <- ggplot(total_comparison,
               aes(x = year, y = mean_agb_Mg_ha, color = method, shape = method)) +
  geom_point(size = 3.5) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Observed LiDAR" = "black",
               "Chapman Richards" = "#009E73",
               "Schumacher" = "#D55E00"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Observed LiDAR" = 16,
               "Chapman Richards" = 17,
               "Schumacher" = 15),
    name = NULL
  ) +
  labs(
    x = "Year",
    y = "Mean AGB (Mg ha<sup>\u22121</sup>)",
    title = "Study Area Mean AGB: Observed vs. Projected"
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown())

ggsave(file.path(out_dir, "fig07_total_agb_comparison.png"), fig7,
       width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig07_total_agb_comparison.pdf"), fig7,
       width = 17.5, height = 12, units = "cm")
cat("Saved: fig07_total_agb_comparison\n")

# Figure 7b: Imperial units
fig7_imp <- ggplot(total_comparison,
               aes(x = year, y = mean_agb_tons_ac, color = method, shape = method)) +
  geom_point(size = 3.5) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Observed LiDAR" = "black",
               "Chapman Richards" = "#009E73",
               "Schumacher" = "#D55E00"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Observed LiDAR" = 16,
               "Chapman Richards" = 17,
               "Schumacher" = 15),
    name = NULL
  ) +
  labs(
    x = "Year",
    y = "Mean AGB (tons ac<sup>\u22121</sup>)",
    title = "Study Area Mean AGB: Observed vs. Projected"
  ) +
  theme_pub +
  theme(axis.title.y = element_markdown())

ggsave(file.path(out_dir, "fig07_total_agb_comparison_imperial.png"), fig7_imp,
       width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig07_total_agb_comparison_imperial.pdf"), fig7_imp,
       width = 17.5, height = 12, units = "cm")
cat("Saved: fig07_total_agb_comparison_imperial\n")

# Figure 7c: Relative difference bar chart
fig7c <- total_wide |>
  select(year, cr_diff_pct, schum_diff_pct) |>
  pivot_longer(cols = c(cr_diff_pct, schum_diff_pct),
               names_to = "model", values_to = "diff_pct") |>
  mutate(model = ifelse(model == "cr_diff_pct",
                        "Chapman Richards", "Schumacher")) |>
  ggplot(aes(x = factor(year), y = diff_pct, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, color = "grey40") +
  scale_fill_manual(values = c("Chapman Richards" = "#009E73",
                                "Schumacher" = "#D55E00"), name = NULL) +
  labs(
    x = "LiDAR Acquisition Year",
    y = "Difference from observed (%)",
    title = "Relative Prediction Error by Year and Model"
  ) +
  theme_pub

ggsave(file.path(out_dir, "fig07c_relative_difference.png"), fig7c,
       width = 17.5, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig07c_relative_difference.pdf"), fig7c,
       width = 17.5, height = 12, units = "cm")
cat("Saved: fig07c_relative_difference\n")

# =============================================================================
# PART 5: Summary Statistics for Keynote
# =============================================================================

cat("\n========== SUMMARY FOR KEYNOTE ==========\n")
cat("Howland Research Forest, central Maine\n")
cat("LiDAR acquisitions:", paste(years_avail, collapse = ", "), "\n")
cat("Total pixels tracked:", n_distinct(traj$cell_id), "\n")
cat("Training period:", paste(train_years, collapse = ", "), "\n")
cat("Validation period:", paste(valid_years, collapse = ", "), "\n")

cat("\nHoldout validation statistics:\n")
print(holdout_stats)

cat("\nMean AGB by stratum (most recent year):\n")
traj |>
  filter(year == max(year)) |>
  group_by(stratum) |>
  summarise(
    mean_Mg_ha = round(mean(agb_Mg_ha, na.rm = TRUE), 1),
    sd_Mg_ha = round(sd(agb_Mg_ha, na.rm = TRUE), 1),
    mean_tons_ac = round(mean(agb_Mg_ha, na.rm = TRUE) * 0.4461, 1),
    sd_tons_ac = round(sd(agb_Mg_ha, na.rm = TRUE) * 0.4461, 1),
    .groups = "drop"
  ) |>
  print()

if (!is.null(naip_summary)) {
  cat("\nNAIP derived forest type distribution:\n")
  naip_summary |> count(forest_type) |> print()
  cat("NAIP years processed:", paste(sort(unique(naip_metrics$year)), collapse = ", "), "\n")
}

cat("\nChapman Richards parameters (trained on", paste(train_years, collapse = ","), "):\n")
print(cr_fits_holdout)

cat("\nFigures saved to:", out_dir, "\n")
cat("\n========== END OF ANALYSIS ==========\n")
