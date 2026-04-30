# =============================================================================
# Title: Recalibrate Yield Curves with 2025 UMaine MELiTE ALS
# Author: A. Weiskittel
# Date: 2026-03-26
# Description: Processes raw 2025 ALS point clouds into a 30m AGB raster,
#              adds the 2025 timepoint to the yield curve calibration window,
#              and refits Chapman-Richards models with 4 timepoints (2012-2025).
#              Compares "before" (3-pt linear fallback) vs "after" (4-pt nonlinear)
#              yield curve performance.
# Dependencies: lidR, terra, tidyverse, nlme
#               2025 .las files in data/als_2025/
#               ORNL DAAC AGB rasters in data/ornl_agb_rasters/
#               Pipeline results in output/fvs_projection_results.rds
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(terra)
library(lidR)
library(sf)
library(nlme)

# --- Paths -------------------------------------------------------------------
als_dir     <- "data/als_2025"
ornl_dir    <- Sys.getenv("ORNL_DIR", unset = "data/ornl_agb_rasters")
output_dir  <- "output"
fig_dir     <- "output/manuscript_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Constants ---------------------------------------------------------------
# ORNL DAAC raster conversion: kgC/m2 to Mg biomass/ha
ORNL_FACTOR <- 20  # 10000/1000 * 2.0

# ALS processing parameters
PIXEL_RES   <- 30  # meters, to match ORNL rasters
P95_TO_AGB  <- 2.5 # approximate p95 height to AGB regression coefficient
P95_POWER   <- 1.8 # power term for height to AGB relationship

cat("=== Script 06: Recalibrate with 2025 ALS ===\n")

# =============================================================================
# PART 1: Process 2025 ALS point clouds to 30m CHM and AGB
# =============================================================================
cat("\n--- Part 1: Processing 2025 ALS point clouds ---\n")

las_files <- list.files(als_dir, pattern = "\\.las$", full.names = TRUE)
cat("Found", length(las_files), "LAS tiles\n")

if (length(las_files) == 0) {
  stop("No .las files found in ", als_dir)
}

# Read as catalog for efficient processing
ctg <- readLAScatalog(las_files)
e <- ext(ctg)
cat("Catalog extent:", e[1], e[2], e[3], e[4], "\n")
cat("Point density:", round(density(ctg), 1), "pts/m2\n")

# Normalize heights using TIN ground classification
# Clean up any prior partial run output first
norm_dir <- file.path(output_dir, "als_2025_processing")
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
old_norm <- list.files(norm_dir, pattern = "_norm\\.la[sz]$", full.names = TRUE)
if (length(old_norm) > 0) {
  cat("Removing", length(old_norm), "leftover normalized files from prior run...\n")
  file.remove(old_norm)
}

opt_output_files(ctg) <- paste0(norm_dir, "/{*}_norm")
opt_chunk_size(ctg) <- 0  # process by file so {*} template works
opt_chunk_buffer(ctg) <- 30

cat("Normalizing point cloud heights...\n")
ctg_norm <- normalize_height(ctg, tin())

# Compute p95 height at 30m resolution
opt_output_files(ctg_norm) <- ""
cat("Computing p95 canopy height at", PIXEL_RES, "m resolution...\n")
p95_raster <- pixel_metrics(ctg_norm, ~quantile(Z, 0.95), res = PIXEL_RES)
names(p95_raster) <- "p95_ht_m"

# Load a reference ORNL raster for alignment
ref_path <- file.path(ornl_dir, "Howland_AGB_2015.tif")
if (file.exists(ref_path)) {
  ref_r <- rast(ref_path)
  cat("Resampling to match ORNL raster grid...\n")
  p95_aligned <- resample(p95_raster, ref_r, method = "bilinear")
} else {
  p95_aligned <- p95_raster
  cat("WARNING: No reference raster found, using native grid\n")
}

# Convert p95 height to AGB using area-based regression
# AGB (Mg/ha) = a * p95^b
# These coefficients should be calibrated from plot data; using approximate values
# that match the ORNL AGB relationship at Howland
cat("Converting p95 height to AGB (Mg/ha)...\n")

# Calibrate p95-to-AGB regression using ORNL 2015 raster
if (file.exists(ref_path)) {
  ornl_2015 <- rast(ref_path) * ORNL_FACTOR

  # Sample collocated pixels for calibration
  p95_vals <- values(p95_aligned, na.rm = FALSE)
  agb_vals <- values(ornl_2015, na.rm = FALSE)

  cal_df <- tibble(p95 = as.numeric(p95_vals), agb = as.numeric(agb_vals)) |>
    filter(!is.na(p95), !is.na(agb), p95 > 2, agb > 5)

  cat("Calibrating p95-to-AGB from", nrow(cal_df), "collocated pixels...\n")

  # Fit power model: AGB = a * p95^b
  p95_fit <- nls(agb ~ a * p95^b, data = cal_df,
                 start = list(a = P95_TO_AGB, b = P95_POWER),
                 control = nls.control(maxiter = 200))

  cat("  Fitted: AGB =", round(coef(p95_fit)["a"], 3),
      "* p95^", round(coef(p95_fit)["b"], 3), "\n")
  cat("  R2 =", round(1 - sum(residuals(p95_fit)^2) /
        sum((cal_df$agb - mean(cal_df$agb))^2), 4), "\n")

  # Apply fitted model to 2025 p95
  agb_2025 <- app(p95_aligned, function(x) {
    coef(p95_fit)["a"] * pmax(x, 0)^coef(p95_fit)["b"]
  })
} else {
  # Fallback: approximate conversion
  agb_2025 <- app(p95_aligned, function(x) {
    P95_TO_AGB * pmax(x, 0)^P95_POWER
  })
}

names(agb_2025) <- "agb_Mg_ha_2025"

# Save processed raster
writeRaster(agb_2025, file.path(output_dir, "agb_2025_melite_30m.tif"),
            overwrite = TRUE)
cat("Saved: agb_2025_melite_30m.tif\n")

# Summary statistics
agb_2025_vals <- values(agb_2025, na.rm = TRUE)
cat("2025 MELiTE AGB summary:\n")
cat("  Mean:", round(mean(agb_2025_vals), 1), "Mg/ha\n")
cat("  Median:", round(median(agb_2025_vals), 1), "Mg/ha\n")
cat("  SD:", round(sd(agb_2025_vals), 1), "Mg/ha\n")
cat("  n pixels:", length(agb_2025_vals), "\n")

# =============================================================================
# PART 2: Build 4-timepoint AGB trajectory dataset
# =============================================================================
cat("\n--- Part 2: Building 4-timepoint trajectory ---\n")

# Load existing ORNL rasters
years_ornl <- c(2012, 2014, 2015)
raster_names <- c("Howland_AGB_2012.tif", "Howland_AGB_2014.tif", "Howland_AGB_2015.tif")

pixel_stack <- list()
for (i in seq_along(years_ornl)) {
  rpath <- file.path(ornl_dir, raster_names[i])
  if (file.exists(rpath)) {
    r <- rast(rpath) * ORNL_FACTOR
    pixel_stack[[as.character(years_ornl[i])]] <- r
    cat("Loaded", raster_names[i], ": mean =",
        round(mean(values(r, na.rm = TRUE)), 1), "Mg/ha\n")
  }
}

# Add 2025 MELiTE raster
pixel_stack[["2025"]] <- agb_2025

# Align all rasters to the same grid
ref_grid <- pixel_stack[["2015"]]
for (yr in names(pixel_stack)) {
  if (yr != "2015") {
    pixel_stack[[yr]] <- resample(pixel_stack[[yr]], ref_grid, method = "near")
  }
}

# Extract collocated pixel values across all years
cat("Extracting collocated pixel trajectories...\n")
all_years <- sort(as.numeric(names(pixel_stack)))
agb_matrix <- sapply(as.character(all_years), function(yr) {
  as.numeric(values(pixel_stack[[yr]], na.rm = FALSE))
})
colnames(agb_matrix) <- paste0("agb_", all_years)

pixel_df <- as_tibble(agb_matrix) |>
  mutate(pixel_id = row_number()) |>
  filter(if_all(starts_with("agb_"), ~!is.na(.) & . > 0))

cat("Valid collocated pixels:", nrow(pixel_df), "\n")

# Reshape to long format
pixel_long <- pixel_df |>
  pivot_longer(cols = starts_with("agb_"),
               names_to = "year_col", values_to = "agb_Mg_ha") |>
  mutate(year = as.numeric(str_extract(year_col, "\\d{4}")),
         year_offset = year - min(all_years))

# Compute landscape means
landscape_means <- pixel_long |>
  group_by(year) |>
  summarise(mean_agb = mean(agb_Mg_ha),
            sd_agb = sd(agb_Mg_ha),
            n = n(), .groups = "drop") |>
  mutate(se = sd_agb / sqrt(n))

cat("\nLandscape trajectory (4 timepoints):\n")
print(landscape_means)

# =============================================================================
# PART 3: Refit yield curves with 4 timepoints
# =============================================================================
cat("\n--- Part 3: Refitting yield curves ---\n")

# Stratify pixels by 2015 AGB
pixel_df <- pixel_df |>
  mutate(agb_class = case_when(
    agb_2015 < quantile(agb_2015, 0.33) ~ "Low",
    agb_2015 < quantile(agb_2015, 0.67) ~ "Medium",
    TRUE ~ "High"
  ))

# Reshape for curve fitting
fit_data <- pixel_df |>
  pivot_longer(cols = matches("^agb_\\d{4}$"),
               names_to = "year_col", values_to = "agb") |>
  mutate(year = as.numeric(str_extract(year_col, "\\d{4}")),
         year_offset = year - 2012) |>  # Use 2012 as t=0
  select(pixel_id, agb_class, year, year_offset, agb)

# Fit Chapman-Richards for each stratum
# AGB = A * (1 - exp(-k * t))^p
fit_cr <- function(df, class_name) {
  cat("  Fitting Chapman-Richards for", class_name, "(n =", nrow(df), ")...\n")

  # Starting values from data
  A_start <- max(df$agb) * 1.2
  k_start <- 0.03
  p_start <- 2.0

  # 3-parameter Chapman-Richards
  cr3 <- tryCatch({
    nls(agb ~ A * (1 - exp(-k * year_offset))^p,
        data = df,
        start = list(A = A_start, k = k_start, p = p_start),
        control = nls.control(maxiter = 500, minFactor = 1/4096))
  }, error = function(e) {
    cat("    3-param failed:", e$message, "\n")
    NULL
  })

  if (!is.null(cr3)) {
    ss_res <- sum(residuals(cr3)^2)
    ss_tot <- sum((df$agb - mean(df$agb))^2)
    r2 <- 1 - ss_res / ss_tot
    cat("    3-param CR: A =", round(coef(cr3)["A"], 1),
        ", k =", round(coef(cr3)["k"], 4),
        ", p =", round(coef(cr3)["p"], 2),
        ", R2 =", round(r2, 4), "\n")
    return(list(fit = cr3, type = "Chapman-Richards 3p", r2 = r2))
  }

  # 2-parameter fallback (fix p = 2.0)
  cr2 <- tryCatch({
    nls(agb ~ A * (1 - exp(-k * year_offset))^2,
        data = df,
        start = list(A = A_start, k = k_start),
        control = nls.control(maxiter = 500))
  }, error = function(e) {
    cat("    2-param failed:", e$message, "\n")
    NULL
  })

  if (!is.null(cr2)) {
    ss_res <- sum(residuals(cr2)^2)
    ss_tot <- sum((df$agb - mean(df$agb))^2)
    r2 <- 1 - ss_res / ss_tot
    cat("    2-param CR: A =", round(coef(cr2)["A"], 1),
        ", k =", round(coef(cr2)["k"], 4),
        ", R2 =", round(r2, 4), "\n")
    return(list(fit = cr2, type = "Chapman-Richards 2p", r2 = r2))
  }

  # Linear fallback
  lm_fit <- lm(agb ~ year_offset, data = df)
  r2 <- summary(lm_fit)$r.squared
  cat("    Linear fallback: slope =", round(coef(lm_fit)[2], 3),
      ", R2 =", round(r2, 4), "\n")
  return(list(fit = lm_fit, type = "Linear", r2 = r2))
}

# Fit for each stratum
strata <- c("Low", "Medium", "High")
yield_fits_4pt <- list()

for (cls in strata) {
  cls_data <- fit_data |> filter(agb_class == cls)
  # Subsample for computational efficiency (use 10K pixels per stratum)
  if (nrow(cls_data) > 40000) {
    set.seed(42)
    sampled_pixels <- sample(unique(cls_data$pixel_id), 10000)
    cls_data <- cls_data |> filter(pixel_id %in% sampled_pixels)
  }
  yield_fits_4pt[[cls]] <- fit_cr(cls_data, cls)
}

# =============================================================================
# PART 4: Compare 3-point (original) vs 4-point (recalibrated) predictions
# =============================================================================
cat("\n--- Part 4: Comparing 3-pt vs 4-pt yield curves ---\n")

# Load original pipeline results
res <- readRDS(file.path(output_dir, "fvs_projection_results.rds"))
yield_fits_3pt <- res$yield_curves

# Predict at 2025 validation year for both calibrations
# For 4-pt curves, 2025 is year_offset = 13 (from 2012)
pred_years <- data.frame(year_offset = c(0, 3, 13, 23, 53))  # 2012, 2015, 2025, 2035, 2065

comparison_table <- tibble()

for (cls in strata) {
  # 4-point predictions
  fit4 <- yield_fits_4pt[[cls]]$fit
  pred4 <- tryCatch(
    predict(fit4, newdata = pred_years),
    error = function(e) rep(NA, nrow(pred_years))
  )

  # 3-point predictions (from original pipeline)
  fit3_row <- yield_fits_3pt |> filter(agb_class == cls)
  if (nrow(fit3_row) > 0) {
    fit3 <- fit3_row$fit[[1]]
    # 3-pt uses different initial_agb parameterization
    pred3 <- tryCatch({
      predict(fit3, newdata = data.frame(
        year_offset = pred_years$year_offset,
        initial_agb = mean(fit_data$agb[fit_data$agb_class == cls & fit_data$year_offset == 0])
      ))
    }, error = function(e) rep(NA, nrow(pred_years)))
  } else {
    pred3 <- rep(NA, nrow(pred_years))
  }

  comparison_table <- bind_rows(comparison_table, tibble(
    stratum = cls,
    year = 2012 + pred_years$year_offset,
    pred_3pt = pred3,
    pred_4pt = pred4,
    model_type_4pt = yield_fits_4pt[[cls]]$type
  ))
}

cat("\nComparison (landscape mean):\n")
comp_summary <- comparison_table |>
  group_by(year) |>
  summarise(mean_3pt = mean(pred_3pt, na.rm = TRUE),
            mean_4pt = mean(pred_4pt, na.rm = TRUE),
            .groups = "drop") |>
  mutate(diff = mean_4pt - mean_3pt)
print(comp_summary)

# =============================================================================
# PART 5: Generate comparison figures
# =============================================================================
cat("\n--- Part 5: Generating recalibration figures ---\n")

library(ggtext)
library(patchwork)

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(8, 8, 8, 8),
    plot.tag = element_text(size = 12, face = "bold")
  )

W_DOUBLE <- 17.5
H_STANDARD <- 12

# Figure: Before/after yield curves by stratum
pred_plot_data <- comparison_table |>
  pivot_longer(cols = c(pred_3pt, pred_4pt),
               names_to = "calibration", values_to = "agb") |>
  mutate(calibration = recode(calibration,
    "pred_3pt" = "3-year calibration (2012\u20132015)",
    "pred_4pt" = "4-year calibration (2012\u20132025)"
  ))

fig_recal <- ggplot(pred_plot_data, aes(x = year, y = agb,
                                         color = calibration,
                                         linetype = calibration)) +
  geom_line(linewidth = 1) +
  geom_point(data = landscape_means, aes(x = year, y = mean_agb),
             inherit.aes = FALSE, size = 3, shape = 16) +
  geom_errorbar(data = landscape_means,
                aes(x = year, ymin = mean_agb - 1.96 * se,
                    ymax = mean_agb + 1.96 * se),
                inherit.aes = FALSE, width = 0.5) +
  facet_wrap(~stratum, scales = "free_y") +
  scale_color_manual(values = c("#E69F00", "#0072B2"), name = NULL) +
  scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
  labs(x = "Year",
       y = "AGB (Mg ha<sup>\u22121</sup>)",
       title = "Yield curve recalibration: 3-year vs 4-year calibration window") +
  annotate("text", x = 2025, y = Inf, label = "Validation", vjust = 1.5,
           size = 2.5, color = "grey40") +
  theme_pub +
  theme(axis.title.y = element_markdown())

ggsave(file.path(fig_dir, "Fig9_recalibration_comparison.pdf"), fig_recal,
       width = W_DOUBLE, height = H_STANDARD, units = "cm")
ggsave(file.path(fig_dir, "Fig9_recalibration_comparison.png"), fig_recal,
       width = W_DOUBLE, height = H_STANDARD, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig9\n")

# =============================================================================
# PART 6: Save recalibration results
# =============================================================================
recal_results <- list(
  agb_2025_raster_path = file.path(output_dir, "agb_2025_melite_30m.tif"),
  landscape_means = landscape_means,
  yield_fits_4pt = yield_fits_4pt,
  comparison_table = comparison_table,
  p95_model = if (exists("p95_fit")) coef(p95_fit) else NULL,
  pixel_trajectories = pixel_df
)

saveRDS(recal_results, file.path(output_dir, "recalibration_results.rds"))
cat("\nSaved recalibration results to output/recalibration_results.rds\n")

cat("\n=== Script 06 complete ===\n")
