# =============================================================================
# Title: FVS-ACD Calibration and Structural Diagnosis
# Author: A. Weiskittel
# Date: 2026-03-26
# Description: Comprehensive calibration of FVS-ACD for Howland AGB projection.
#   Three calibration approaches:
#     A) Simple growth multiplier (1D) - mirrors FVS-NE approach
#     B) Growth x mortality 2D grid search
#     C) Density-dependent transfer: apply ACD growth-density curve to LiDAR pixels
#   Diagnoses why FVS-ACD overpredicts by decomposing annual growth dynamics
#   and showing the baseline mismatch as the dominant error source.
# =============================================================================

library(tidyverse)
library(ggtext)
library(patchwork)
library(terra)

output_dir <- "output"
ornl_dir   <- Sys.getenv("ORNL_DIR", unset = "data/ornl_agb_rasters")
fig_dir    <- "output/manuscript_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ORNL_FACTOR <- 20

cat("=== Script 07: FVS-ACD Calibration and Structural Diagnosis ===\n\n")

# =============================================================================
# PART 1: Load data
# =============================================================================
cat("--- Part 1: Loading data ---\n")

res <- readRDS(file.path(output_dir, "fvs_projection_results.rds"))
acd_trajectory <- res$acd_trajectory

# Observed LiDAR landscape values
ornl_2015 <- rast(file.path(ornl_dir, "Howland_AGB_2015.tif")) * ORNL_FACTOR
obs_2015 <- mean(values(ornl_2015, na.rm = TRUE))

# Check for 2025 MELiTE raster first, fall back to ORNL
obs_2025_path <- file.path(output_dir, "agb_2025_melite_30m.tif")
if (!file.exists(obs_2025_path)) {
  obs_2025_path <- file.path(ornl_dir, "Howland_AGB_2025.tif")
}
ornl_2025 <- rast(obs_2025_path)
if (grepl("ornl|ORNL", obs_2025_path)) ornl_2025 <- ornl_2025 * ORNL_FACTOR
obs_2025 <- mean(values(ornl_2025, na.rm = TRUE))
obs_change <- obs_2025 - obs_2015

cat("Observed LiDAR landscape:\n")
cat("  2015 mean:", round(obs_2015, 1), "Mg/ha\n")
cat("  2025 mean:", round(obs_2025, 1), "Mg/ha\n")
cat("  10yr change:", round(obs_change, 1), "Mg/ha (",
    sprintf("%+.1f%%", 100 * obs_change / obs_2015), ")\n\n")

# Extract per-plot ACD summaries
acd_0  <- acd_trajectory |> filter(year_offset == 0)
acd_10 <- acd_trajectory |> filter(year_offset == 10)

plot_data <- acd_0 |>
  select(PLOT, agb_0 = agb_Mg_ha, ba_0 = ba_m2_ha,
         n0 = n_trees, dbh_0 = mean_dbh_cm) |>
  left_join(
    acd_10 |> select(PLOT, agb_10 = agb_Mg_ha, ba_10 = ba_m2_ha,
                      n10 = n_trees, dbh_10 = mean_dbh_cm),
    by = "PLOT"
  ) |>
  mutate(
    agb_change = agb_10 - agb_0,
    ba_change = ba_10 - ba_0,
    n_change = n10 - n0,
    mortality_pct = 100 * (n0 - n10) / n0,
    rel_change = 100 * agb_change / agb_0
  )

n_plots <- nrow(plot_data)
cat("ACD trajectory: ", n_plots, " plots\n")
cat("  Mean initial AGB:", round(mean(plot_data$agb_0), 1), "Mg/ha\n")
cat("  Mean initial BA:", round(mean(plot_data$ba_0), 1), "m2/ha\n")
cat("  Mean 10yr AGB:", round(mean(plot_data$agb_10), 1), "Mg/ha\n")
cat("  Mean 10yr change:", round(mean(plot_data$agb_change), 1), "Mg/ha\n")
cat("  Plots with negative growth:", sum(plot_data$agb_change < 0), "of", n_plots, "\n")
cat("  Plots with tree mortality:", sum(plot_data$n_change < 0), "of", n_plots, "\n\n")

# =============================================================================
# PART 2: Build annual growth-density relationship from ACD trajectory
# =============================================================================
cat("--- Part 2: Growth-density relationship ---\n")

# Compute annual increments for all plots
annual_inc <- acd_trajectory |>
  arrange(PLOT, year_offset) |>
  group_by(PLOT) |>
  mutate(
    annual_agb_inc = agb_Mg_ha - lag(agb_Mg_ha),
    annual_ba_inc = ba_m2_ha - lag(ba_m2_ha),
    annual_mort = lag(n_trees) - n_trees,
    lag_agb = lag(agb_Mg_ha),
    lag_ba = lag(ba_m2_ha)
  ) |>
  ungroup() |>
  filter(!is.na(annual_agb_inc))

cat("Annual increments computed:", nrow(annual_inc), "plot-years\n")

# Fit growth-density curve: annual AGB increment as f(current AGB)
# This captures how ACD's growth rate changes with stocking level
gd_fit <- loess(annual_agb_inc ~ lag_agb, data = annual_inc, span = 0.5)

# Predict growth rate across the full AGB range
agb_range <- seq(20, 350, by = 5)
pred_inc <- predict(gd_fit, newdata = data.frame(lag_agb = agb_range))

gd_curve <- tibble(agb = agb_range, annual_inc = pred_inc)

# Find the crossover point (where annual increment = 0)
crossover <- approx(pred_inc, agb_range, xout = 0)$y
cat("ACD growth-density crossover (zero increment):",
    round(crossover, 1), "Mg/ha\n")
cat("ACD annual increment at landscape mean (123 Mg/ha):",
    round(predict(gd_fit, newdata = data.frame(lag_agb = 123)), 3), "Mg/ha/yr\n")
cat("ACD annual increment at plot mean (202 Mg/ha):",
    round(predict(gd_fit, newdata = data.frame(lag_agb = 202)), 3), "Mg/ha/yr\n\n")

# =============================================================================
# PART 3: Calibration Method A - Simple growth multiplier (1D)
# =============================================================================
cat("--- Part 3: Method A - Growth multiplier calibration ---\n")

growth_mults <- seq(0.0, 5.0, by = 0.05)
results_A <- tibble(
  multiplier = growth_mults,
  mean_change = sapply(growth_mults, function(m) {
    mean(m * plot_data$agb_change, na.rm = TRUE)
  }),
  mean_2025 = sapply(growth_mults, function(m) {
    mean(plot_data$agb_0 + m * plot_data$agb_change, na.rm = TRUE)
  })
) |>
  mutate(
    abs_error_change = abs(mean_change - obs_change),
    pct_error = 100 * (mean_2025 - obs_2025) / obs_2025
  )

best_A <- results_A |> slice_min(abs_error_change, n = 1)
cat("Method A (growth multiplier):\n")
cat("  Best multiplier:", round(best_A$multiplier, 2), "\n")
cat("  Predicted change:", round(best_A$mean_change, 1), "Mg/ha",
    "(observed:", round(obs_change, 1), ")\n")
cat("  Residual error:", round(best_A$abs_error_change, 1), "Mg/ha\n")
cat("  Landscape % error:", sprintf("%+.1f%%", best_A$pct_error), "\n")

# Can ANY multiplier match observed change?
can_match_A <- min(results_A$abs_error_change) < 1.0
cat("  Can match observed change:", can_match_A, "\n\n")

# =============================================================================
# PART 4: Calibration Method B - Growth x Mortality 2D grid search
# =============================================================================
cat("--- Part 4: Method B - 2D growth x mortality calibration ---\n")

# Decompose ACD net change into gross growth vs mortality loss
# For plots with no mortality (n_change = 0), all change is growth
# For plots with mortality, estimate gross growth assuming surviving trees
# grew at the same rate as the no-mortality cohort

# Annual decomposition
annual_decomp <- annual_inc |>
  mutate(
    trees_died = pmax(0, annual_mort),
    # Approximate mortality AGB loss: assume dead trees had mean AGB
    mort_loss = ifelse(trees_died > 0,
                       trees_died * lag_agb / lag(n_trees),
                       0),
    # Gross growth = net change + mortality loss
    gross_growth = annual_agb_inc + mort_loss
  )

# Summarize over 10-year window per plot
plot_decomp <- annual_decomp |>
  filter(year_offset <= 10) |>
  group_by(PLOT) |>
  summarise(
    agb_0 = first(lag_agb),
    total_gross_growth = sum(gross_growth, na.rm = TRUE),
    total_mort_loss = sum(mort_loss, na.rm = TRUE),
    net_change = sum(annual_agb_inc, na.rm = TRUE),
    .groups = "drop"
  )

cat("10yr decomposition (plot means):\n")
cat("  Gross growth:", round(mean(plot_decomp$total_gross_growth), 1), "Mg/ha\n")
cat("  Mortality loss:", round(mean(plot_decomp$total_mort_loss), 1), "Mg/ha\n")
cat("  Net change:", round(mean(plot_decomp$net_change), 1), "Mg/ha\n\n")

# 2D grid: growth_mult x mortality_mult
g_mults <- seq(0.5, 5.0, by = 0.1)
m_mults <- seq(0.0, 2.0, by = 0.1)

grid_2d <- expand_grid(g_mult = g_mults, m_mult = m_mults) |>
  mutate(
    mean_change = map2_dbl(g_mult, m_mult, function(gm, mm) {
      calibrated_net <- gm * plot_decomp$total_gross_growth -
                        mm * plot_decomp$total_mort_loss
      mean(calibrated_net, na.rm = TRUE)
    }),
    mean_2025 = map2_dbl(g_mult, m_mult, function(gm, mm) {
      calibrated_net <- gm * plot_decomp$total_gross_growth -
                        mm * plot_decomp$total_mort_loss
      mean(plot_decomp$agb_0 + calibrated_net, na.rm = TRUE)
    })
  ) |>
  mutate(
    abs_error = abs(mean_change - obs_change),
    pct_error = 100 * (mean_2025 - obs_2025) / obs_2025
  )

best_B <- grid_2d |> slice_min(abs_error, n = 1) |> slice(1)
cat("Method B (growth x mortality):\n")
cat("  Best growth mult:", round(best_B$g_mult, 1), "\n")
cat("  Best mortality mult:", round(best_B$m_mult, 1), "\n")
cat("  Predicted change:", round(best_B$mean_change, 1), "Mg/ha\n")
cat("  Residual error:", round(best_B$abs_error, 1), "Mg/ha\n")
cat("  Landscape % error:", sprintf("%+.1f%%", best_B$pct_error), "\n\n")

# =============================================================================
# PART 5: Method C - Apply ACD growth-density curve to LiDAR pixel distribution
# =============================================================================
cat("--- Part 5: Method C - Density-dependent landscape transfer ---\n")

# Get LiDAR pixel AGB distribution
pixel_vals <- values(ornl_2015, na.rm = TRUE)
pixel_vals <- pixel_vals[pixel_vals > 0]
n_pixels <- length(pixel_vals)
cat("LiDAR pixels:", n_pixels, "\n")
cat("Pixel AGB range:", round(min(pixel_vals), 1), "to",
    round(max(pixel_vals), 1), "Mg/ha\n")

# Build a fast lookup table from the loess fit for vectorized prediction
# Precompute increments at 0.5 Mg/ha resolution across the full AGB range
lut_agb <- seq(0, 400, by = 0.5)
lut_inc <- predict(gd_fit, newdata = data.frame(lag_agb = lut_agb))
lut_inc[is.na(lut_inc)] <- 0

# Fast vectorized prediction: uses approx() lookup instead of per-pixel loess
predict_vec <- function(agb_vec, years, growth_mult = 1.0) {
  agb <- agb_vec
  for (yr in 1:years) {
    inc <- approx(lut_agb, lut_inc, xout = pmin(pmax(agb, 0), 400),
                  rule = 2)$y
    agb <- agb + growth_mult * inc
    agb <- pmax(agb, 0)
  }
  return(agb)
}

# Apply with default (mult = 1.0)
cat("Applying ACD growth curve to", n_pixels, "pixels (default, vectorized)...\n")
pred_default <- predict_vec(pixel_vals, years = 10, growth_mult = 1.0)
cat("  Default landscape mean 2025:", round(mean(pred_default), 1), "Mg/ha\n")
cat("  Default landscape change:", round(mean(pred_default) - mean(pixel_vals), 1), "Mg/ha\n")

# Grid search over growth multiplier applied to landscape
cat("Grid searching over growth multiplier for landscape application...\n")
c_mults <- seq(0.5, 5.0, by = 0.1)
results_C <- tibble(multiplier = c_mults, mean_2025 = NA_real_, mean_change = NA_real_)

for (i in seq_along(c_mults)) {
  pred_c <- predict_vec(pixel_vals, years = 10, growth_mult = c_mults[i])
  results_C$mean_2025[i] <- mean(pred_c)
  results_C$mean_change[i] <- mean(pred_c) - mean(pixel_vals)
}

results_C <- results_C |>
  mutate(
    abs_error = abs(mean_change - obs_change),
    pct_error = 100 * (mean_2025 - obs_2025) / obs_2025
  )

best_C <- results_C |> slice_min(abs_error, n = 1) |> slice(1)
cat("\nMethod C (landscape transfer):\n")
cat("  Best multiplier:", round(best_C$multiplier, 1), "\n")
cat("  Predicted change:", round(best_C$mean_change, 1), "Mg/ha",
    "(observed:", round(obs_change, 1), ")\n")
cat("  Landscape % error:", sprintf("%+.1f%%", best_C$pct_error), "\n")

# Apply best C multiplier and get pixel-level predictions
pred_C_best <- predict_vec(pixel_vals, years = 10, growth_mult = best_C$multiplier)

# =============================================================================
# PART 6: Summary comparison table
# =============================================================================
cat("\n--- Part 6: Summary comparison ---\n\n")

comparison <- tribble(
  ~Method, ~Description, ~Multiplier, ~Pred_Change, ~Obs_Change, ~Error_pct,
  "Default (plot)", "ACD default from field plots", "1.0",
    round(mean(plot_data$agb_change), 1), round(obs_change, 1),
    sprintf("%+.1f", 100*(mean(plot_data$agb_10)-obs_2025)/obs_2025),
  "A (plot)", "1D growth multiplier, plot-based", as.character(round(best_A$multiplier, 2)),
    round(best_A$mean_change, 1), round(obs_change, 1),
    sprintf("%+.1f", best_A$pct_error),
  "B (plot)", "2D growth x mortality, plot-based",
    paste0("g=", round(best_B$g_mult, 1), ", m=", round(best_B$m_mult, 1)),
    round(best_B$mean_change, 1), round(obs_change, 1),
    sprintf("%+.1f", best_B$pct_error),
  "C (landscape)", "ACD curve on LiDAR pixels", as.character(round(best_C$multiplier, 1)),
    round(best_C$mean_change, 1), round(obs_change, 1),
    sprintf("%+.1f", best_C$pct_error)
)

print(comparison)

# =============================================================================
# PART 7: Generate figures
# =============================================================================
cat("\n--- Part 7: Generating figures ---\n")

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

# --- Panel A: ACD growth-density curve ---
# Show where plots vs landscape sit on this curve
fig10a <- ggplot() +
  geom_point(data = annual_inc,
             aes(x = lag_agb, y = annual_agb_inc),
             alpha = 0.08, size = 0.5, color = "grey60") +
  geom_line(data = gd_curve |> filter(!is.na(annual_inc)),
            aes(x = agb, y = annual_inc),
            linewidth = 1.2, color = "#D55E00") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  # Annotate landscape vs plot means
  geom_vline(xintercept = obs_2015, linetype = "dotted",
             color = "#0072B2", linewidth = 0.8) +
  geom_vline(xintercept = mean(plot_data$agb_0), linetype = "dotted",
             color = "#D55E00", linewidth = 0.8) +
  annotate("text", x = obs_2015 - 2, y = max(gd_curve$annual_inc, na.rm = TRUE) * 0.85,
           label = paste0("Landscape\nmean\n(", round(obs_2015), ")"),
           hjust = 1, size = 2.5, color = "#0072B2", fontface = "bold") +
  annotate("text", x = mean(plot_data$agb_0) + 2,
           y = max(gd_curve$annual_inc, na.rm = TRUE) * 0.85,
           label = paste0("Plot mean\n(", round(mean(plot_data$agb_0)), ")"),
           hjust = 0, size = 2.5, color = "#D55E00", fontface = "bold") +
  # Shade the zone where plots live
  annotate("rect",
           xmin = quantile(plot_data$agb_0, 0.05),
           xmax = quantile(plot_data$agb_0, 0.95),
           ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.08) +
  labs(x = "Current AGB (Mg ha\u207B\u00B9)",
       y = "Annual AGB increment (Mg ha\u207B\u00B9 yr\u207B\u00B9)",
       title = "A) ACD growth-density curve") +
  theme_pub

# --- Panel B: 2D calibration surface ---
fig10b <- ggplot(grid_2d, aes(x = g_mult, y = m_mult, fill = abs_error)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", direction = -1,
                       name = "|Error|\n(Mg/ha)",
                       limits = c(0, max(20, max(grid_2d$abs_error)))) +
  geom_point(data = best_B, aes(x = g_mult, y = m_mult),
             shape = 4, size = 4, stroke = 2, color = "white") +
  geom_point(aes(x = 1.0, y = 1.0),
             shape = 1, size = 3, stroke = 1.5, color = "white") +
  annotate("text", x = 1.0, y = 1.15, label = "Default",
           color = "white", size = 2.5) +
  labs(x = "Growth multiplier",
       y = "Mortality multiplier",
       title = "B) 2D calibration surface") +
  theme_pub

# --- Panel C: Method C - landscape transfer with growth multiplier ---
fig10c <- ggplot(results_C, aes(x = multiplier, y = mean_change)) +
  geom_line(linewidth = 0.8, color = "#0072B2") +
  geom_hline(yintercept = obs_change, linetype = "dashed", color = "#D55E00") +
  geom_point(data = best_C, aes(x = multiplier, y = mean_change),
             size = 3, color = "#0072B2") +
  annotate("text", x = max(c_mults) * 0.7, y = obs_change + 0.5,
           label = paste0("Observed (+", round(obs_change, 1), ")"),
           color = "#D55E00", size = 2.8, fontface = "bold") +
  annotate("text", x = best_C$multiplier + 0.2, y = best_C$mean_change,
           label = paste0("Optimal = ", round(best_C$multiplier, 1)),
           size = 2.8, color = "#0072B2", fontface = "bold", hjust = 0) +
  labs(x = "Growth multiplier",
       y = "Predicted 10yr \u0394AGB (Mg ha\u207B\u00B9)",
       title = "C) ACD curve applied to LiDAR pixels") +
  theme_pub

# --- Panel D: Per-plot comparison: default vs landscape-transferred ---
# Show that the fundamental issue is baseline, not growth rate
plot_comp_d <- tibble(
  source = c(rep("Field plots\n(default ACD)", n_plots),
             rep(paste0("LiDAR pixels\n(ACD curve, mult=",
                        round(best_C$multiplier, 1), ")"),
                 min(5000, n_pixels))),
  init_agb = c(plot_data$agb_0, sample(pixel_vals, min(5000, n_pixels))),
  pred_2025 = c(plot_data$agb_10,
                predict_vec(sample(pixel_vals, min(5000, n_pixels)),
                            years = 10, growth_mult = best_C$multiplier))
) |>
  mutate(change = pred_2025 - init_agb)

fig10d <- ggplot(plot_comp_d, aes(x = init_agb, y = change, color = source)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  geom_hline(yintercept = obs_change, linetype = "dashed", color = "black") +
  annotate("text", x = 10, y = obs_change + 1,
           label = paste0("Observed landscape \u0394AGB = +", round(obs_change, 1)),
           hjust = 0, size = 2.5, fontface = "bold") +
  scale_color_manual(values = c("#D55E00", "#0072B2"), name = NULL) +
  labs(x = "Initial AGB (Mg ha\u207B\u00B9)",
       y = "10yr \u0394AGB (Mg ha\u207B\u00B9)",
       title = "D) Baseline mismatch drives FVS-ACD error") +
  theme_pub +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = "white", color = NA))

# Combine
fig10 <- (fig10a | fig10b) / (fig10c | fig10d) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D")))

ggsave(file.path(fig_dir, "Fig10_acd_calibration.pdf"), fig10,
       width = W_DOUBLE, height = W_DOUBLE * 0.85, units = "cm")
ggsave(file.path(fig_dir, "Fig10_acd_calibration.png"), fig10,
       width = W_DOUBLE, height = W_DOUBLE * 0.85, units = "cm",
       dpi = 300, bg = "white")
cat("  Saved Fig10 (4-panel)\n")

# Also save a standalone growth-density curve figure (for potential supplemental)
fig_gd_standalone <- ggplot() +
  geom_point(data = annual_inc,
             aes(x = lag_agb, y = annual_agb_inc),
             alpha = 0.05, size = 0.5, color = "grey60") +
  geom_line(data = gd_curve |> filter(!is.na(annual_inc)),
            aes(x = agb, y = annual_inc),
            linewidth = 1.2, color = "#D55E00") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = obs_2015, linetype = "dotted",
             color = "#0072B2", linewidth = 0.8) +
  geom_vline(xintercept = mean(plot_data$agb_0), linetype = "dotted",
             color = "#D55E00", linewidth = 0.8) +
  # Add rug for pixel and plot distributions
  geom_rug(data = tibble(x = sample(pixel_vals, 2000)),
           aes(x = x), sides = "b", alpha = 0.05, color = "#0072B2") +
  geom_rug(data = tibble(x = plot_data$agb_0),
           aes(x = x), sides = "t", alpha = 0.3, color = "#D55E00") +
  annotate("text", x = obs_2015 - 2, y = max(gd_curve$annual_inc, na.rm = TRUE) * 0.9,
           label = paste0("LiDAR landscape\nmean = ", round(obs_2015), " Mg/ha"),
           hjust = 1, size = 3, color = "#0072B2", fontface = "bold") +
  annotate("text", x = mean(plot_data$agb_0) + 2,
           y = max(gd_curve$annual_inc, na.rm = TRUE) * 0.9,
           label = paste0("Field plot\nmean = ", round(mean(plot_data$agb_0)), " Mg/ha"),
           hjust = 0, size = 3, color = "#D55E00", fontface = "bold") +
  annotate("text", x = crossover + 5, y = -0.3,
           label = paste0("Zero growth\n= ", round(crossover), " Mg/ha"),
           hjust = 0, size = 2.5, color = "grey40") +
  labs(x = expression("Current AGB (Mg ha"^-1*")"),
       y = expression("Annual AGB increment (Mg ha"^-1*" yr"^-1*")"),
       title = "FVS-ACD growth-density relationship at Howland") +
  theme_pub

ggsave(file.path(fig_dir, "FigS4_acd_growth_density.pdf"), fig_gd_standalone,
       width = 8.5, height = 8.5, units = "cm")
ggsave(file.path(fig_dir, "FigS4_acd_growth_density.png"), fig_gd_standalone,
       width = 8.5, height = 8.5, units = "cm", dpi = 300, bg = "white")
cat("  Saved FigS4 (growth-density standalone)\n")

# =============================================================================
# PART 8: Save results
# =============================================================================
acd_cal_results <- list(
  # Method summaries
  method_A = list(
    best_multiplier = best_A$multiplier,
    pred_change = best_A$mean_change,
    pct_error = best_A$pct_error,
    can_match = can_match_A
  ),
  method_B = list(
    best_growth_mult = best_B$g_mult,
    best_mort_mult = best_B$m_mult,
    pred_change = best_B$mean_change,
    pct_error = best_B$pct_error
  ),
  method_C = list(
    best_multiplier = best_C$multiplier,
    pred_change = best_C$mean_change,
    pct_error = best_C$pct_error
  ),
  # Diagnostic data
  growth_density_curve = gd_curve,
  crossover_agb = crossover,
  plot_decomposition = plot_decomp,
  annual_increments = annual_inc,
  comparison_table = comparison,
  obs_change = obs_change,
  obs_2025 = obs_2025,
  obs_2015 = obs_2015
)

saveRDS(acd_cal_results, file.path(output_dir, "acd_calibration_results.rds"))
cat("\nSaved calibration results\n")

cat("\n=== Script 07 complete ===\n")
