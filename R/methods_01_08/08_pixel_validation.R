# =============================================================================
# Title: Pixel-Level Validation and Residual Mapping
# Author: A. Weiskittel
# Date: 2026-03-26
# Description: Compares predicted vs observed 2025 AGB at the pixel level
#              (30m resolution) for methods that produce wall-to-wall predictions.
#              Computes spatial RMSE, bias maps, and pixel-level R2. Generates
#              residual heat maps showing where each method over/underpredicts.
# Dependencies: terra, tidyverse
#               fvs_projection_results.rds
#               recalibration_results.rds (from Script 06)
#               ORNL AGB rasters
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(terra)
library(ggtext)
library(patchwork)

# --- Paths -------------------------------------------------------------------
output_dir <- "output"
ornl_dir   <- Sys.getenv("ORNL_DIR", unset = "data/ornl_agb_rasters")
fig_dir    <- "output/manuscript_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Constants ---------------------------------------------------------------
ORNL_FACTOR <- 20  # kgC/m2 to Mg biomass/ha

cat("=== Script 08: Pixel-Level Validation ===\n")

# =============================================================================
# PART 1: Load observed and predicted 2025 rasters
# =============================================================================
cat("\n--- Part 1: Loading rasters ---\n")

# Observed 2025 AGB (from MELiTE processing in Script 06)
obs_2025_path <- file.path(output_dir, "agb_2025_melite_30m.tif")
if (!file.exists(obs_2025_path)) {
  cat("WARNING: MELiTE 2025 raster not found. Using ORNL 2025 if available.\n")
  obs_2025_path <- file.path(ornl_dir, "Howland_AGB_2025.tif")
}

obs_2025 <- rast(obs_2025_path)
if (grepl("ornl", obs_2025_path, ignore.case = TRUE) ||
    grepl("ORNL", obs_2025_path)) {
  obs_2025 <- obs_2025 * ORNL_FACTOR
}
cat("Observed 2025 raster: mean =", round(mean(values(obs_2025, na.rm = TRUE)), 1),
    "Mg/ha\n")

# Load 2015 baseline
agb_2015_path <- file.path(ornl_dir, "Howland_AGB_2015.tif")
agb_2015 <- rast(agb_2015_path) * ORNL_FACTOR
cat("Baseline 2015 raster: mean =", round(mean(values(agb_2015, na.rm = TRUE)), 1),
    "Mg/ha\n")

# Load pipeline results for yield curve predictions
res <- readRDS(file.path(output_dir, "fvs_projection_results.rds"))
yield_curves <- res$yield_curves

# Load recalibration if available
recal_path <- file.path(output_dir, "recalibration_results.rds")
if (file.exists(recal_path)) {
  recal <- readRDS(recal_path)
  has_recal <- TRUE
  cat("Loaded recalibration results\n")
} else {
  has_recal <- FALSE
  cat("No recalibration results (Script 06 not yet run)\n")
}

# =============================================================================
# PART 2: Generate pixel-level predictions for wall-to-wall methods
# =============================================================================
cat("\n--- Part 2: Generating pixel-level predictions ---\n")

# Align rasters
obs_aligned <- resample(obs_2025, agb_2015, method = "near")

# Extract valid pixels
agb_2015_vals <- values(agb_2015, na.rm = FALSE)
obs_2025_vals <- values(obs_aligned, na.rm = FALSE)
valid <- !is.na(agb_2015_vals) & !is.na(obs_2025_vals) &
         agb_2015_vals > 0 & obs_2025_vals > 0

n_valid <- sum(valid)
cat("Valid collocated pixels:", n_valid, "\n")

init_agb <- agb_2015_vals[valid]
obs_agb  <- obs_2025_vals[valid]

# Pre-compute class breaks for yield curve predictions
class_breaks <- yield_curves %>%
  mutate(mean_init = sapply(data, function(d) mean(d$initial_agb[d$year_offset == 0], na.rm = TRUE))) %>%
  select(agb_class, mean_init)

if (nrow(class_breaks) >= 3) {
  low_med_break  <- mean(class_breaks$mean_init[class_breaks$agb_class %in% c("Low", "Medium")])
  med_high_break <- mean(class_breaks$mean_init[class_breaks$agb_class %in% c("Medium", "High")])
} else {
  low_med_break  <- quantile(init_agb, 0.33)
  med_high_break <- quantile(init_agb, 0.67)
}

# Vectorized prediction function
predict_yield_vec <- function(agb_vec, years, yield_fits) {
  pixel_class <- case_when(
    agb_vec < low_med_break  ~ "Low",
    agb_vec < med_high_break ~ "Medium",
    TRUE ~ "High"
  )
  result <- numeric(length(agb_vec))
  for (cls in unique(pixel_class)) {
    idx <- which(pixel_class == cls)
    fit <- yield_fits$fit[yield_fits$agb_class == cls][[1]]
    pred_data <- data.frame(year_offset = rep(years, length(idx)),
                            initial_agb = agb_vec[idx])
    preds <- tryCatch(
      predict(fit, newdata = pred_data),
      error = function(e) agb_vec[idx]
    )
    result[idx] <- preds
  }
  result
}

# Method 1: Area-based CR yield curves (3-point)
cat("Predicting: Area-based CR (3-pt)...\n")
pred_cr3 <- predict_yield_vec(init_agb, 10, yield_curves)

# Method 2: Area-based CR yield curves (4-point, if available)
if (has_recal) {
  cat("Predicting: Area-based CR (4-pt recalibrated)...\n")
  # Use 4-pt yield fits from recalibration
  # Need to adapt predict function for new parameterization
  pred_cr4 <- tryCatch({
    yf4 <- recal$yield_fits_4pt
    pixel_class <- case_when(
      init_agb < low_med_break  ~ "Low",
      init_agb < med_high_break ~ "Medium",
      TRUE ~ "High"
    )
    result <- numeric(length(init_agb))
    for (cls in c("Low", "Medium", "High")) {
      idx <- which(pixel_class == cls)
      fit <- yf4[[cls]]$fit
      pred_data <- data.frame(year_offset = rep(13, length(idx)),  # 2012+13=2025
                              initial_agb = init_agb[idx])
      result[idx] <- tryCatch(
        predict(fit, newdata = pred_data),
        error = function(e) init_agb[idx]
      )
    }
    result
  }, error = function(e) {
    cat("  4-pt prediction failed:", e$message, "\n")
    rep(NA, length(init_agb))
  })
} else {
  pred_cr4 <- rep(NA, length(init_agb))
}

# =============================================================================
# PART 3: Compute pixel-level accuracy metrics
# =============================================================================
cat("\n--- Part 3: Computing accuracy metrics ---\n")

pixel_results <- tibble(
  init_agb = init_agb,
  obs_2025 = obs_agb,
  obs_change = obs_agb - init_agb,
  pred_cr3 = pred_cr3,
  resid_cr3 = pred_cr3 - obs_agb,
  pred_cr4 = pred_cr4,
  resid_cr4 = pred_cr4 - obs_agb
)

# Accuracy metrics function
compute_metrics <- function(pred, obs) {
  valid <- !is.na(pred) & !is.na(obs)
  p <- pred[valid]; o <- obs[valid]
  tibble(
    n = sum(valid),
    bias = mean(p - o),
    mae = mean(abs(p - o)),
    rmse = sqrt(mean((p - o)^2)),
    r2 = cor(p, o)^2,
    bias_pct = 100 * mean(p - o) / mean(o),
    rmse_pct = 100 * sqrt(mean((p - o)^2)) / mean(o)
  )
}

metrics <- bind_rows(
  compute_metrics(pixel_results$pred_cr3, pixel_results$obs_2025) |>
    mutate(method = "3-pt yield curves"),
  compute_metrics(pixel_results$pred_cr4, pixel_results$obs_2025) |>
    mutate(method = "4-pt yield curves (recalibrated)")
) |>
  relocate(method)

cat("\nPixel-level accuracy:\n")
print(metrics)

# =============================================================================
# PART 4: Generate validation figures
# =============================================================================
cat("\n--- Part 4: Generating validation figures ---\n")

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(8, 8, 8, 8),
    plot.tag = element_text(size = 12, face = "bold")
  )

W_DOUBLE <- 17.5
H_STANDARD <- 12

# Subsample for plotting efficiency
set.seed(42)
n_plot <- min(nrow(pixel_results), 50000)
plot_df <- pixel_results |> sample_n(n_plot)

# Panel A: Predicted vs Observed (3-pt)
fig11a <- ggplot(plot_df, aes(x = obs_2025, y = pred_cr3)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(option = "C", name = "Count", trans = "log10") +
  geom_abline(color = "red", linewidth = 0.8) +
  labs(x = "Observed 2025 AGB (Mg ha<sup>\u22121</sup>)",
       y = "Predicted 2025 AGB (Mg ha<sup>\u22121</sup>)",
       title = "A) 3-year calibration") +
  annotate("text", x = 10, y = max(plot_df$pred_cr3, na.rm = TRUE) * 0.9,
           label = paste0("RMSE = ", round(metrics$rmse[1], 1),
                         "\nBias = ", sprintf("%+.1f", metrics$bias[1]),
                         "\nR\u00B2 = ", round(metrics$r2[1], 3)),
           hjust = 0, size = 3, fontface = "bold") +
  coord_equal() +
  theme_pub +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

# Panel B: Predicted vs Observed (4-pt, if available)
if (!all(is.na(plot_df$pred_cr4))) {
  fig11b <- ggplot(plot_df, aes(x = obs_2025, y = pred_cr4)) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(option = "C", name = "Count", trans = "log10") +
    geom_abline(color = "red", linewidth = 0.8) +
    labs(x = "Observed 2025 AGB (Mg ha<sup>\u22121</sup>)",
         y = "Predicted 2025 AGB (Mg ha<sup>\u22121</sup>)",
         title = "B) 4-year calibration (recalibrated)") +
    annotate("text", x = 10, y = max(plot_df$pred_cr4, na.rm = TRUE) * 0.9,
             label = paste0("RMSE = ", round(metrics$rmse[2], 1),
                           "\nBias = ", sprintf("%+.1f", metrics$bias[2]),
                           "\nR\u00B2 = ", round(metrics$r2[2], 3)),
             hjust = 0, size = 3, fontface = "bold") +
    coord_equal() +
    theme_pub +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())

  fig11 <- fig11a + fig11b
} else {
  fig11 <- fig11a
}

ggsave(file.path(fig_dir, "Fig11_pixel_validation.pdf"), fig11,
       width = W_DOUBLE, height = H_STANDARD * 0.85, units = "cm")
ggsave(file.path(fig_dir, "Fig11_pixel_validation.png"), fig11,
       width = W_DOUBLE, height = H_STANDARD * 0.85, units = "cm",
       dpi = 300, bg = "white")
cat("  Saved Fig11\n")

# Panel C: Residual map (spatial)
cat("Generating residual map...\n")
resid_raster <- agb_2015  # template
resid_vals <- rep(NA, ncell(resid_raster))
resid_vals[valid] <- pixel_results$resid_cr3
values(resid_raster) <- resid_vals

resid_df <- as.data.frame(resid_raster, xy = TRUE) |>
  rename(residual = 3) |>
  filter(!is.na(residual))

# Clip extreme residuals for visualization
resid_df <- resid_df |>
  mutate(residual_clipped = pmax(-100, pmin(100, residual)))

fig12 <- ggplot(resid_df, aes(x = x, y = y, fill = residual_clipped)) +
  geom_raster() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Residual\n(Mg ha<sup>\u22121</sup>)",
                       limits = c(-100, 100)) +
  coord_equal() +
  labs(title = "Prediction residuals (3-pt yield curves, 2025)",
       subtitle = "Blue = underprediction, Red = overprediction",
       x = "Easting (m)", y = "Northing (m)") +
  theme_pub +
  theme(legend.title = element_markdown(),
        legend.key.height = unit(1.2, "cm"))

ggsave(file.path(fig_dir, "Fig12_residual_map.pdf"), fig12,
       width = 8.5, height = 10.5, units = "cm")
ggsave(file.path(fig_dir, "Fig12_residual_map.png"), fig12,
       width = 8.5, height = 10.5, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig12\n")

# =============================================================================
# PART 5: Save pixel-level validation results
# =============================================================================
pixel_val_results <- list(
  metrics = metrics,
  pixel_results_sample = plot_df,
  n_total_pixels = n_valid,
  resid_raster_path = file.path(fig_dir, "residual_raster.tif")
)

# Save residual raster
writeRaster(resid_raster, pixel_val_results$resid_raster_path, overwrite = TRUE)

saveRDS(pixel_val_results, file.path(output_dir, "pixel_validation_results.rds"))
cat("\nSaved pixel validation results\n")

cat("\n=== Script 08 complete ===\n")
