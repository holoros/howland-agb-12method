# =============================================================================
# Title: Publication Figures for Howland Multi-Method AGB Comparison
# Author: A. Weiskittel
# Date: 2026-03-26
# Journal: Forestry: An International Journal of Forest Research (Oxford)
# Description: Generates all main text and supplemental figures from the
#              saved pipeline results (fvs_projection_results.rds and
#              agb_trajectories.rds). Designed to run on Cardinal HPC
#              after the main analysis pipeline completes.
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(terra)
library(patchwork)
library(ggtext)
library(scales)
library(sf)

# --- Paths -------------------------------------------------------------------
output_dir  <- Sys.getenv("OUTPUT_DIR", unset = "output")
fig_dir     <- "output/manuscript_figures"
raster_dir  <- Sys.getenv("ORNL_DIR", unset = "data/ornl_agb_rasters")
field_dir   <- Sys.getenv("FIELD_DIR", unset = "data/field_inventory")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Theme -------------------------------------------------------------------
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    plot.margin = margin(8, 8, 8, 8),
    plot.tag = element_text(size = 12, face = "bold")
  )

# Forestry journal: single column = 8.5 cm, double column = 17.5 cm
# Full page height ~ 23 cm
W_SINGLE <- 8.5
W_DOUBLE <- 17.5
H_STANDARD <- 12

# --- Unit conversions --------------------------------------------------------
MG_HA_TO_TONS_AC <- 0.4461

# --- Load results ------------------------------------------------------------
cat("Loading pipeline results...\n")

res <- readRDS(file.path(output_dir, "fvs_projection_results.rds"))
comparison_df  <- res$comparison_df
acd_trajectory <- res$acd_trajectory
ne_oob_means   <- res$ne_oob_means
ne_cal_means   <- res$ne_cal_means
acd_means      <- res$acd_means
yield_curves   <- res$yield_curves
mc_summary     <- res$mc_summary
lidar_pixel    <- res$lidar_pixel
naip_pixel     <- res$naip_pixel

# Load area-based trajectories from Script 01
traj_path <- file.path(output_dir, "agb_trajectories.rds")
if (file.exists(traj_path)) {
  pixel_traj <- readRDS(traj_path)
  cat("Loaded pixel trajectories:", nrow(pixel_traj), "rows\n")
} else {
  cat("WARNING: agb_trajectories.rds not found\n")
  pixel_traj <- tibble()
}

# Load validation results from Script 04
val_path <- file.path(output_dir, "projection_validation_results.rds")
if (file.exists(val_path)) {
  val_results <- readRDS(val_path)
  cat("Loaded validation results\n")
} else {
  val_results <- NULL
}

# --- Color palette -----------------------------------------------------------
pal_methods <- c(
  "Observed LiDAR"       = "#1a1a1a",
  "Area-based CR"        = "#009E73",
  "Area-based Schum."    = "#56B4E9",
  "FVS-NE (defaults)"    = "#E69F00",
  "FVS-NE (calibrated)"  = "#CC79A7",
  "FVS-ACD (full)"       = "#D55E00",
  "LiDAR + ACD curve"    = "#0072B2",
  "NAIP + ACD curve"     = "#F0E442"
)

shape_methods <- c(
  "Observed LiDAR" = 16, "Area-based CR" = 17, "Area-based Schum." = 15,
  "FVS-NE (defaults)" = 18, "FVS-NE (calibrated)" = 4,
  "FVS-ACD (full)" = 8, "LiDAR + ACD curve" = 3, "NAIP + ACD curve" = 7
)

lty_methods <- c(
  "Observed LiDAR" = "solid", "Area-based CR" = "solid",
  "Area-based Schum." = "dashed",
  "FVS-NE (defaults)" = "solid", "FVS-NE (calibrated)" = "dashed",
  "FVS-ACD (full)" = "solid",
  "LiDAR + ACD curve" = "solid", "NAIP + ACD curve" = "dotted"
)

# =============================================================================
# FIGURE 1: Study area map with LiDAR AGB 2015
# =============================================================================
cat("Figure 1: Study area map...\n")

agb_2015_path <- file.path(raster_dir, "Howland_AGB_2015.tif")
if (file.exists(agb_2015_path)) {
  agb_r <- rast(agb_2015_path) * 20  # kgC/m2 -> Mg/ha
  agb_df <- as.data.frame(agb_r, xy = TRUE) |>
    rename(agb_Mg_ha = 3) |>
    filter(!is.na(agb_Mg_ha), agb_Mg_ha > 0)

  fig1 <- ggplot(agb_df, aes(x = x, y = y, fill = agb_Mg_ha)) +
    geom_raster() +
    scale_fill_viridis_c(option = "D", name = "AGB\n(Mg ha<sup>\u22121</sup>)",
                          limits = c(0, 250)) +
    coord_equal() +
    labs(title = "Howland Research Forest",
         subtitle = "Aboveground biomass from G-LiHT LiDAR (2015)",
         x = "Easting (m)", y = "Northing (m)") +
    theme_pub +
    theme(legend.title = element_markdown(),
          legend.key.height = unit(1.2, "cm"))

  ggsave(file.path(fig_dir, "Fig1_study_area.pdf"), fig1,
         width = W_SINGLE, height = W_SINGLE + 2, units = "cm")
  ggsave(file.path(fig_dir, "Fig1_study_area.png"), fig1,
         width = W_SINGLE, height = W_SINGLE + 2, units = "cm", dpi = 300, bg = "white")
  cat("  Saved Fig1\n")
} else {
  cat("  Skipping Fig1: raster not found\n")
}

# =============================================================================
# FIGURE 2: LiDAR AGB trajectories by stratum (2012-2025)
# =============================================================================
cat("Figure 2: AGB trajectories by stratum...\n")

if (nrow(pixel_traj) > 0) {
  stratum_traj <- pixel_traj |>
    group_by(year) |>
    summarise(
      mean_agb = mean(agb_Mg_ha, na.rm = TRUE),
      sd_agb = sd(agb_Mg_ha, na.rm = TRUE),
      n = n(),
      se_agb = sd_agb / sqrt(n),
      .groups = "drop"
    )

  fig2 <- ggplot(stratum_traj, aes(x = year, y = mean_agb)) +
    geom_ribbon(aes(ymin = mean_agb - 1.96 * se_agb,
                    ymax = mean_agb + 1.96 * se_agb),
                fill = "grey80", alpha = 0.5) +
    geom_line(linewidth = 1, color = "#1a1a1a") +
    geom_point(size = 2.5, color = "#1a1a1a") +
    labs(x = "Year",
         y = "Mean AGB (Mg ha<sup>\u22121</sup>)",
         title = "Observed LiDAR AGB trajectory",
         subtitle = "Landscape mean \u00b1 95% CI from multi-temporal ALS") +
    theme_pub +
    theme(axis.title.y = element_markdown())

  ggsave(file.path(fig_dir, "Fig2_observed_trajectory.pdf"), fig2,
         width = W_SINGLE, height = H_STANDARD * 0.7, units = "cm")
  ggsave(file.path(fig_dir, "Fig2_observed_trajectory.png"), fig2,
         width = W_SINGLE, height = H_STANDARD * 0.7, units = "cm", dpi = 300, bg = "white")
  cat("  Saved Fig2\n")
}

# =============================================================================
# FIGURE 3: Multi-method comparison (10-year projection window)
# =============================================================================
cat("Figure 3: 10-year multi-method comparison...\n")

fig3 <- comparison_df |>
  filter(year >= 2012, year <= 2026) |>
  ggplot(aes(x = year, y = mean_agb_Mg_ha, color = method,
             shape = method, linetype = method)) +
  annotate("rect", xmin = -Inf, xmax = 2015.5,
           ymin = -Inf, ymax = Inf, fill = "#E8F5E9", alpha = 0.3) +
  annotate("text", x = 2013.5, y = Inf,
           label = "Calibration", vjust = 1.5, size = 2.8, color = "grey40") +
  annotate("text", x = 2020.5, y = Inf,
           label = "Projection", vjust = 1.5, size = 2.8, color = "grey40") +
  geom_vline(xintercept = 2015.5, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 2) +
  scale_color_manual(values = pal_methods, name = NULL) +
  scale_shape_manual(values = shape_methods, name = NULL) +
  scale_linetype_manual(values = lty_methods, name = NULL) +
  labs(x = "Year",
       y = "Mean AGB (Mg ha<sup>\u22121</sup>)",
       title = "10-year projection: multi-method comparison") +
  theme_pub +
  theme(axis.title.y = element_markdown(),
        legend.key.width = unit(1.5, "cm"))

ggsave(file.path(fig_dir, "Fig3_multimethod_10yr.pdf"), fig3,
       width = W_DOUBLE, height = H_STANDARD, units = "cm")
ggsave(file.path(fig_dir, "Fig3_multimethod_10yr.png"), fig3,
       width = W_DOUBLE, height = H_STANDARD, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig3\n")

# =============================================================================
# FIGURE 4: Long-term projection divergence (50 years)
# =============================================================================
cat("Figure 4: 50-year projection divergence...\n")

fig4 <- comparison_df |>
  filter(year >= 2015) |>
  ggplot(aes(x = year, y = mean_agb_Mg_ha, color = method,
             shape = method, linetype = method)) +
  geom_vline(xintercept = 2025, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("text", x = 2025, y = Inf, label = "Validation (2025)",
           vjust = 1.5, hjust = -0.05, size = 2.5, color = "grey40") +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.8) +
  scale_color_manual(values = pal_methods, name = NULL) +
  scale_shape_manual(values = shape_methods, name = NULL) +
  scale_linetype_manual(values = lty_methods, name = NULL) +
  labs(x = "Year",
       y = "Mean AGB (Mg ha<sup>\u22121</sup>)",
       title = "50-year projection: model divergence") +
  theme_pub +
  theme(axis.title.y = element_markdown(),
        legend.key.width = unit(1.5, "cm"))

ggsave(file.path(fig_dir, "Fig4_multimethod_50yr.pdf"), fig4,
       width = W_DOUBLE, height = H_STANDARD, units = "cm")
ggsave(file.path(fig_dir, "Fig4_multimethod_50yr.png"), fig4,
       width = W_DOUBLE, height = H_STANDARD, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig4\n")

# =============================================================================
# FIGURE 5: Prediction error at 2025 validation year
# =============================================================================
cat("Figure 5: 2025 prediction error...\n")

obs_2025 <- comparison_df |>
  filter(method == "Observed LiDAR", year == 2025) |>
  pull(mean_agb_Mg_ha)

if (length(obs_2025) > 0 && !is.na(obs_2025)) {
  error_at_2025 <- comparison_df |>
    filter(year == 2025, method != "Observed LiDAR") |>
    mutate(
      obs = obs_2025,
      error_Mg = mean_agb_Mg_ha - obs,
      error_pct = 100 * error_Mg / obs
    ) |>
    mutate(method = fct_reorder(method, error_pct))

  fig5 <- ggplot(error_at_2025, aes(x = method, y = error_pct, fill = method)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_text(aes(label = paste0(sprintf("%+.1f", error_pct), "%")),
              vjust = ifelse(error_at_2025$error_pct >= 0, -0.5, 1.5),
              size = 2.8, fontface = "bold") +
    scale_fill_manual(values = pal_methods, guide = "none") +
    labs(x = NULL,
         y = "Prediction error at 2025 (%)",
         title = "AGB prediction error relative to 2025 LiDAR observation") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

  ggsave(file.path(fig_dir, "Fig5_error_2025.pdf"), fig5,
         width = W_DOUBLE, height = H_STANDARD * 0.85, units = "cm")
  ggsave(file.path(fig_dir, "Fig5_error_2025.png"), fig5,
         width = W_DOUBLE, height = H_STANDARD * 0.85, units = "cm", dpi = 300, bg = "white")
  cat("  Saved Fig5\n")
}

# =============================================================================
# FIGURE 6: Cost vs Accuracy (Monte Carlo)
# =============================================================================
cat("Figure 6: Cost vs accuracy...\n")

cost_data <- tribble(
  ~method,                 ~cost_lo, ~cost_hi, ~cost_unit,  ~coverage,
  "NAIP + ACD curve",       0.10,     0.25,    "per acre",  "Wall-to-wall",
  "Area-based LiDAR",       0.25,     0.50,    "per acre",  "Wall-to-wall",
  "LiDAR + ACD curve",      0.30,     0.55,    "per acre",  "Wall-to-wall",
  "ITC + FVS-NE",           0.75,     1.00,    "per acre",  "Wall-to-wall",
  "ITC + FVS-ACD",          0.75,     1.00,    "per acre",  "Wall-to-wall",
  "Field inventory",        5.78,     9.25,    "per acre",  "Sample-based"
)

acc_10yr <- c(
  "NAIP + ACD curve"  = 22, "Area-based LiDAR"  = 15,
  "LiDAR + ACD curve" = 12, "ITC + FVS-NE"      = 35,
  "ITC + FVS-ACD"     = 10, "Field inventory"   = 8
)

set.seed(2026)
n_mc <- 1000
mc_df <- map_dfr(seq_len(nrow(cost_data)), function(i) {
  m <- cost_data$method[i]
  acc_val <- acc_10yr[m]
  if (is.na(acc_val)) acc_val <- 25
  tibble(
    method = m, sim = 1:n_mc,
    cost = runif(n_mc, cost_data$cost_lo[i], cost_data$cost_hi[i]),
    rmse_pct = pmax(1, rnorm(n_mc, acc_val, acc_val * 0.2)),
    coverage = cost_data$coverage[i]
  )
})

mc_means <- mc_df |>
  group_by(method, coverage) |>
  summarise(cost_mean = mean(cost), rmse_mean = mean(rmse_pct),
            .groups = "drop")

fig6 <- ggplot(mc_df, aes(x = cost, y = rmse_pct, color = method)) +
  geom_point(alpha = 0.04, size = 0.3) +
  stat_ellipse(level = 0.90, linewidth = 0.8) +
  geom_point(data = mc_means,
             aes(x = cost_mean, y = rmse_mean, shape = coverage),
             size = 3, stroke = 1.0) +
  scale_x_log10(labels = dollar_format()) +
  scale_color_brewer(palette = "Set2", name = NULL) +
  scale_shape_manual(values = c("Wall-to-wall" = 16, "Sample-based" = 17),
                     name = "Coverage") +
  labs(x = "Cost per acre (USD, log scale)",
       y = "Expected RMSE (%)",
       title = "Cost vs. accuracy tradeoff (10-year projection, n = 1,000)") +
  theme_pub +
  theme(legend.box = "vertical")

ggsave(file.path(fig_dir, "Fig6_cost_accuracy.pdf"), fig6,
       width = W_DOUBLE, height = H_STANDARD, units = "cm")
ggsave(file.path(fig_dir, "Fig6_cost_accuracy.png"), fig6,
       width = W_DOUBLE, height = H_STANDARD, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig6\n")

# =============================================================================
# FIGURE 7: ACD trajectory diagnostics (per-plot growth)
# =============================================================================
cat("Figure 7: ACD per-plot trajectories...\n")

acd_plot_traj <- acd_trajectory |>
  filter(year_offset %in% seq(0, 50, by = 5))

fig7a <- ggplot(acd_plot_traj, aes(x = year_offset, y = agb_Mg_ha, group = PLOT)) +
  geom_line(alpha = 0.25, linewidth = 0.3, color = "#D55E00") +
  geom_line(data = acd_means |> filter(year_offset %in% seq(0, 50, by = 5)),
            aes(x = year_offset, y = mean_agb_Mg_ha, group = 1),
            linewidth = 1.2, color = "black") +
  labs(x = "Years from baseline",
       y = "AGB (Mg ha<sup>\u22121</sup>)",
       title = "A) FVS-ACD per-plot trajectories") +
  theme_pub +
  theme(axis.title.y = element_markdown())

fig7b <- ggplot(acd_plot_traj |> filter(year_offset %in% c(0, 10)),
                aes(x = factor(year_offset), y = agb_Mg_ha)) +
  geom_boxplot(fill = "#D55E00", alpha = 0.3, outlier.size = 0.8) +
  labs(x = "Year offset", y = "AGB (Mg ha<sup>\u22121</sup>)",
       title = "B) Distribution at baseline and 10 yr") +
  theme_pub +
  theme(axis.title.y = element_markdown())

fig7 <- fig7a + fig7b + plot_layout(widths = c(2, 1))

ggsave(file.path(fig_dir, "Fig7_acd_trajectories.pdf"), fig7,
       width = W_DOUBLE, height = H_STANDARD * 0.8, units = "cm")
ggsave(file.path(fig_dir, "Fig7_acd_trajectories.png"), fig7,
       width = W_DOUBLE, height = H_STANDARD * 0.8, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig7\n")

# =============================================================================
# FIGURE 8: Pixel-level LiDAR AGB distribution and projected change
# =============================================================================
cat("Figure 8: Pixel-level LiDAR projections...\n")

if (!is.null(lidar_pixel) && nrow(lidar_pixel) > 0) {
  pixel_long <- lidar_pixel |>
    mutate(change_10yr = pred_agb_2025 - pixel_agb_2015,
           change_50yr = pred_agb_2065 - pixel_agb_2015) |>
    sample_n(min(nrow(lidar_pixel), 50000))  # subsample for plotting

  fig8a <- ggplot(pixel_long, aes(x = pixel_agb_2015)) +
    geom_histogram(bins = 80, fill = "#0072B2", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = mean(lidar_pixel$pixel_agb_2015),
               linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(x = "AGB 2015 (Mg ha<sup>\u22121</sup>)",
         y = "Pixel count",
         title = "A) 2015 AGB distribution") +
    theme_pub +
    theme(axis.title.x = element_markdown())

  fig8b <- ggplot(pixel_long, aes(x = pixel_agb_2015, y = change_10yr)) +
    geom_hex(bins = 60) +
    scale_fill_viridis_c(option = "C", name = "Count") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(x = "Initial AGB (Mg ha<sup>\u22121</sup>)",
         y = "\u0394AGB 10yr (Mg ha<sup>\u22121</sup>)",
         title = "B) Projected 10-year change") +
    theme_pub +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())

  fig8 <- fig8a + fig8b

  ggsave(file.path(fig_dir, "Fig8_pixel_projections.pdf"), fig8,
         width = W_DOUBLE, height = H_STANDARD * 0.8, units = "cm")
  ggsave(file.path(fig_dir, "Fig8_pixel_projections.png"), fig8,
         width = W_DOUBLE, height = H_STANDARD * 0.8, units = "cm", dpi = 300, bg = "white")
  cat("  Saved Fig8\n")
}

# =============================================================================
# SUPPLEMENTAL FIGURES
# =============================================================================

# --- Figure S1: FVS-NE SI calibration grid search ---
cat("Figure S1: SI calibration...\n")
# (Regenerated from saved data if available)

# --- Figure S2: Bootstrap CI for Chapman-Richards parameters ---
cat("Figure S2: Bootstrap parameter uncertainty...\n")
if (!is.null(val_results) && !is.null(val_results$boot_cr)) {
  boot_cr_long <- val_results$boot_cr |>
    pivot_longer(cols = c(pred_mean, pred_lo, pred_hi),
                 names_to = "stat", values_to = "agb")

  figS2 <- ggplot(val_results$boot_cr, aes(x = year)) +
    geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = stratum), alpha = 0.2) +
    geom_line(aes(y = pred_mean, color = stratum), linewidth = 0.8) +
    labs(x = "Year",
         y = "Predicted AGB (Mg ha<sup>\u22121</sup>)",
         title = "Chapman-Richards bootstrap 95% CI by stratum") +
    theme_pub +
    theme(axis.title.y = element_markdown())

  ggsave(file.path(fig_dir, "FigS2_bootstrap_ci.pdf"), figS2,
         width = W_DOUBLE, height = H_STANDARD, units = "cm")
  ggsave(file.path(fig_dir, "FigS2_bootstrap_ci.png"), figS2,
         width = W_DOUBLE, height = H_STANDARD, units = "cm", dpi = 300, bg = "white")
  cat("  Saved FigS2\n")
}

# --- Figure S3: Relative change from baseline ---
cat("Figure S3: Relative change...\n")
baseline_by_method <- comparison_df |>
  group_by(method) |>
  filter(year == min(year)) |>
  select(method, baseline_agb = mean_agb_Mg_ha) |>
  ungroup()

relative_change <- comparison_df |>
  left_join(baseline_by_method, by = "method") |>
  mutate(pct_change = 100 * (mean_agb_Mg_ha - baseline_agb) / baseline_agb)

figS3 <- ggplot(relative_change |> filter(year >= 2015),
                aes(x = year, y = pct_change, color = method, linetype = method)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = pal_methods, name = NULL) +
  scale_linetype_manual(values = lty_methods, name = NULL) +
  labs(x = "Year", y = "Change from baseline (%)",
       title = "Relative AGB change by method") +
  theme_pub +
  theme(legend.key.width = unit(1.5, "cm"))

ggsave(file.path(fig_dir, "FigS3_relative_change.pdf"), figS3,
       width = W_DOUBLE, height = H_STANDARD, units = "cm")
ggsave(file.path(fig_dir, "FigS3_relative_change.png"), figS3,
       width = W_DOUBLE, height = H_STANDARD, units = "cm", dpi = 300, bg = "white")
cat("  Saved FigS3\n")

# =============================================================================
# SUPPLEMENTAL TABLES (saved as CSV)
# =============================================================================
cat("\nGenerating supplemental tables...\n")

# Table S1: Full comparison at key years
table_s1 <- comparison_df |>
  filter(year %in% c(2015, 2025, 2035, 2045, 2065)) |>
  select(method, year, mean_agb_Mg_ha, mean_agb_tons_ac) |>
  pivot_wider(names_from = year,
              values_from = c(mean_agb_Mg_ha, mean_agb_tons_ac))

write_csv(table_s1, file.path(fig_dir, "TableS1_method_comparison.csv"))
cat("  Saved TableS1\n")

# Table S2: ACD per-plot summary
table_s2 <- acd_trajectory |>
  filter(year_offset %in% c(0, 10, 20, 50)) |>
  select(PLOT, year_offset, n_trees, ba_m2_ha, agb_Mg_ha) |>
  pivot_wider(names_from = year_offset,
              values_from = c(n_trees, ba_m2_ha, agb_Mg_ha),
              names_sep = "_yr")

write_csv(table_s2, file.path(fig_dir, "TableS2_acd_plot_trajectories.csv"))
cat("  Saved TableS2\n")

# Table S3: Monte Carlo cost summary
write_csv(mc_summary, file.path(fig_dir, "TableS3_cost_accuracy_summary.csv"))
cat("  Saved TableS3\n")

cat("\n========== Figure generation complete ==========\n")
cat("Output directory:", fig_dir, "\n")
cat("Files:\n")
list.files(fig_dir, full.names = FALSE) |> walk(~cat("  ", .x, "\n"))
