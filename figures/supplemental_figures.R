# ==============================================================================
# Supplementary figures S5, S6, S7 for Methods 9-12
# S5: Bootstrap CI density plots for Methods 9, 10, 11, 12 RMSE
# S6: Stratum × composition calibration lines (M11 and M12)
# S7: ITC 2017 and TLS 2025 spatial coverage map
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(scales)
  library(patchwork); library(ggrepel); library(sf); library(terra)
})

hrf_root   <- "/users/PUOM0008/crsfaaron/HRF"
output_dir <- file.path(hrf_root, "output")
fig_dir    <- file.path(output_dir, "manuscript_figures_v45")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

r <- readRDS(file.path(output_dir, "refinements_results.rds"))

# ==============================================================================
# FIGURE S5: Bootstrap RMSE distributions for Methods 9-12
# ==============================================================================
cat("Building Figure S5 (bootstrap RMSE densities)...\n")

# Recompute bootstrap samples for the full distribution (not just summary)
# Use the cached pixel level data
set.seed(42)
B <- 1000

boot_samples <- function(obs, pred, B = 1000) {
  n <- length(obs)
  replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    sqrt(mean((pred[idx] - obs[idx])^2))
  })
}

# Method 9
cm <- readRDS(file.path(output_dir, "curve_matching_results.rds"))
recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))
pt <- recal$pixel_trajectories |> as_tibble()
m9 <- cm$snap |> filter(target_year == 2025) |>
  inner_join(pt |> select(pixel_id, agb_2025), by = "pixel_id") |>
  drop_na(agb_pred, agb_2025)
m9_b <- boot_samples(m9$agb_2025, m9$agb_pred, B = B)
m9_rel <- 100 * m9_b / mean(m9$agb_2025)

# Method 10 holdout
m10_t <- r$m10_test
m10_b <- boot_samples(m10_t$agb_2025_obs, m10_t$pred_hold, B = B)
m10_rel <- 100 * m10_b / mean(m10_t$agb_2025_obs)

# Method 11 stratified
m11_p <- r$m11_strat$pred |>
  filter(!is.na(agb_pred_2025_strat), agb_2025 > 0)
m11_b <- boot_samples(m11_p$agb_2025, m11_p$agb_pred_2025_strat, B = B)
m11_rel <- 100 * m11_b / mean(m11_p$agb_2025)

# Method 12 stratified
m12_p <- r$m12_strat$pred |>
  filter(!is.na(agb_pred_2025_strat), agb_2025 > 0)
m12_b <- boot_samples(m12_p$agb_2025, m12_p$agb_pred_2025_strat, B = B)
m12_rel <- 100 * m12_b / mean(m12_p$agb_2025)

boot_df <- bind_rows(
  tibble(method = "9 Curve matching",            rmse_pct = m9_rel),
  tibble(method = "10 Pure Kohek (holdout)",     rmse_pct = m10_rel),
  tibble(method = "11 Hybrid A (stratified)",    rmse_pct = m11_rel),
  tibble(method = "12 Hybrid B (stratified)",    rmse_pct = m12_rel)
) |> mutate(method = factor(method, levels = c(
    "9 Curve matching", "10 Pure Kohek (holdout)",
    "11 Hybrid A (stratified)", "12 Hybrid B (stratified)")))

ci_df <- boot_df |>
  group_by(method) |>
  summarise(med = median(rmse_pct),
             lo  = quantile(rmse_pct, 0.025),
             hi  = quantile(rmse_pct, 0.975),
             .groups = "drop")

pal4 <- c(
  "9 Curve matching"            = "#6E4C9E",
  "10 Pure Kohek (holdout)"     = "#2CA02C",
  "11 Hybrid A (stratified)"    = "#D62728",
  "12 Hybrid B (stratified)"    = "#17BECF")

figS5 <- ggplot(boot_df, aes(rmse_pct, fill = method, color = method)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  geom_vline(data = ci_df, aes(xintercept = med, color = method),
              linetype = "dashed", linewidth = 0.5) +
  geom_text(data = ci_df,
             aes(x = med, y = 0, color = method,
                  label = sprintf("%.1f [%.1f, %.1f]", med, lo, hi)),
             vjust = -0.5, size = 3.0, fontface = "bold",
             show.legend = FALSE) +
  facet_wrap(~ method, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pal4, guide = "none") +
  scale_color_manual(values = pal4, guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  labs(x = "Bootstrap RMSE (% of mean observed 2025 AGB)",
        y = "Density",
        title = "Figure S5. Bootstrap RMSE distributions for Methods 9 through 12",
        subtitle = sprintf("n = %d resamples; median and 95%% CI labeled. Tight intervals reflect large pixel n (>9k each).", B)) +
  theme_pub +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "FigS5_bootstrap_methods9_12.png"), figS5,
        width = 18, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "FigS5_bootstrap_methods9_12.pdf"), figS5,
        width = 18, height = 12, units = "cm")
cat("  Saved FigS5\n")

# ==============================================================================
# FIGURE S6: Stratified calibration lines for Method 11 and 12
# ==============================================================================
cat("Building Figure S6 (stratified calibration)...\n")

m11_fits <- r$m11_strat$fits |>
  mutate(model = "Method 11 (literature H-D, default Kohek)")
m12_fits <- r$m12_strat$fits |>
  mutate(model = "Method 12 (TLS H-D, calibrated Kohek)")

fits_combined <- bind_rows(m11_fits, m12_fits)

# Sample raw vs observed pixel data per group for visualization
# r$m11_strat$pred and r$m12_strat$pred already contain agb_2015 and agb_2025
m11_pred <- r$m11_strat$pred |>
  as_tibble() |>
  mutate(agb_obs_2017 = agb_2015 + (agb_2025 - agb_2015) * 2 / 10,
         model = "Method 11 (literature H-D, default Kohek)") |>
  filter(agb_raw_2017 > 0, agb_obs_2017 > 0, n_trees >= 2)

m12_pred <- r$m12_strat$pred |>
  as_tibble() |>
  mutate(agb_obs_2017 = agb_2015 + (agb_2025 - agb_2015) * 2 / 10,
         model = "Method 12 (TLS H-D, calibrated Kohek)") |>
  filter(agb_raw_2017 > 0, agb_obs_2017 > 0, n_trees >= 2)

# Sample up to N per group for plot legibility
sample_per_group <- function(d, n_per = 200) {
  d |> group_by(group, model) |>
    group_modify(~ slice_sample(.x, n = min(n_per, nrow(.x)))) |>
    ungroup()
}
plot_pts <- bind_rows(sample_per_group(m11_pred), sample_per_group(m12_pred))

# Plot raw vs obs colored by group with regression lines
figS6 <- ggplot(plot_pts, aes(agb_raw_2017, agb_obs_2017, color = group)) +
  geom_point(alpha = 0.20, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_wrap(~ model) +
  scale_color_brewer(palette = "Set1", name = "Stratum × composition") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 350)) +
  labs(x = "Raw ITC-derived AGB (Mg ha⁻¹, before calibration)",
        y = "Observed AGB at 2017 anchor (Mg ha⁻¹)",
        title = "Figure S6. Stratified calibration: 9 fits per method recover per-pixel discrimination",
        subtitle = "Each line is one stratum × composition cell; slopes vary 0.2 to 6.9, intercepts 0 to 196 Mg/ha.") +
  theme_pub +
  theme(legend.position = "right",
         legend.text = element_text(size = 8),
         strip.text = element_text(face = "bold"),
         plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave(file.path(fig_dir, "FigS6_stratified_calibration.png"), figS6,
        width = 22, height = 12, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "FigS6_stratified_calibration.pdf"), figS6,
        width = 22, height = 12, units = "cm")
cat("  Saved FigS6\n")

# ==============================================================================
# FIGURE S7: ITC 2017 + TLS 2025 spatial coverage at Howland
# ==============================================================================
cat("Building Figure S7 (ITC + TLS coverage)...\n")

itc <- st_read(file.path(output_dir, "Howland_ITC_2017.gpkg"), quiet = TRUE) |>
  filter(Z >= 2, Z <= 45)
itc_xy <- st_coordinates(itc)[, 1:2]

# TLS 2025
tls <- read.table(file.path(hrf_root, "data/tls_2025/Merged_data_total.txt"),
                   header = FALSE, skip = 1) |>
  setNames(c("X","Y","Z","TileID","TreeID","GlobalTreeID",
             "Height_m","DBH_cm","conf_loc","conf_DBH"))

# Compute ITC density on a 50 m grid
ref_raster <- rast(file.path(output_dir, "agb_2025_melite_30m.tif"))
ext_r <- ext(ref_raster)
density_grid <- rast(ext = ext_r, resolution = 50, crs = crs(ref_raster))
itc_pts <- vect(itc_xy, type = "points", crs = crs(ref_raster))
density_grid <- rasterize(itc_pts, density_grid, fun = "count")

# Build dataframes for ggplot
density_df <- as.data.frame(density_grid, xy = TRUE) |>
  setNames(c("x", "y", "n_trees")) |>
  filter(!is.na(n_trees))

tls_df <- tibble(X = tls$X, Y = tls$Y)

bbox_tls <- c(xmin = min(tls$X) - 20, xmax = max(tls$X) + 20,
              ymin = min(tls$Y) - 20, ymax = max(tls$Y) + 20)

panel_a <- ggplot() +
  geom_raster(data = density_df,
               aes(x = x, y = y, fill = n_trees)) +
  geom_rect(aes(xmin = bbox_tls["xmin"], xmax = bbox_tls["xmax"],
                 ymin = bbox_tls["ymin"], ymax = bbox_tls["ymax"]),
             fill = NA, color = "red", linewidth = 0.8) +
  scale_fill_viridis_c(name = "ITC trees\n(per 50 m cell)", option = "B",
                        trans = "log10") +
  coord_equal() +
  labs(title = "A. ITC 2017 wall-to-wall density (208,632 trees)",
        subtitle = "Red rectangle = 2025 TLS plot footprint",
        x = "Easting (m, UTM 19N)", y = "Northing (m, UTM 19N)") +
  theme_pub +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"),
         legend.position = "right")

panel_b <- ggplot(tls_df, aes(X, Y)) +
  geom_point(color = "#D62728", size = 0.4, alpha = 0.6) +
  coord_equal(xlim = bbox_tls[c("xmin","xmax")],
               ylim = bbox_tls[c("ymin","ymax")]) +
  labs(title = "B. TLS 2025 tree positions (4,417 trees)",
        subtitle = "Plot footprint, calibration source for Method 12",
        x = "Easting (m, UTM 19N)", y = "Northing (m, UTM 19N)") +
  theme_pub +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))

figS7 <- (panel_a | panel_b) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = "Figure S7. ITC 2017 wall-to-wall coverage and 2025 TLS calibration footprint at Howland",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "FigS7_itc_tls_coverage.png"), figS7,
        width = 22, height = 11, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "FigS7_itc_tls_coverage.pdf"), figS7,
        width = 22, height = 11, units = "cm")
cat("  Saved FigS7\n")

# Save bootstrap CI table (also useful for TableS8)
ci_table <- ci_df |>
  mutate(across(c(med, lo, hi), \(x) round(x, 2)))
write_csv(ci_table, file.path(fig_dir, "TableS_boot_methods9_12.csv"))

# Save fits table for Table S
fits_export <- fits_combined |>
  select(model, group, n, alpha, beta) |>
  mutate(across(c(alpha, beta), \(x) round(x, 3)))
write_csv(fits_export, file.path(fig_dir, "TableS_stratified_fits.csv"))

cat("\n=== Supplementary figures complete ===\n")
