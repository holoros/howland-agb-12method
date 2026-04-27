# ==============================================================================
# Methods 11 and 12 with INDEPENDENT 2017 ORNL ABA AGB anchor
# Replaces the linearly interpolated 2017 reference with the actual 2017 ORNL
# G-LiHT raster, eliminating the data leakage between calibration target and
# 2025 validation reference.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(FNN)
})

hrf_root   <- "/users/PUOM0008/crsfaaron/HRF"
output_dir <- file.path(hrf_root, "output")
ornl_dir   <- file.path(hrf_root, "data/ornl_agb_rasters")

# --- 1. Load reference grid and ORNL 2017 anchor + 2025 target ---------------
ref_raster <- rast(file.path(output_dir, "agb_2025_melite_30m.tif"))
pixel_area_ha <- prod(res(ref_raster)) / 10000

agb_2017_ornl <- rast(file.path(ornl_dir, "Howland_AGB_2017.tif"))
agb_2025_ornl <- rast(file.path(ornl_dir, "Howland_AGB_2025.tif"))
cat("ORNL 2017 raster:", res(agb_2017_ornl), "m,", ncell(agb_2017_ornl), "cells\n")
cat("ORNL 2025 raster:", res(agb_2025_ornl), "m,", ncell(agb_2025_ornl), "cells\n")

# Resample to ref grid if needed
agb_2017_resampled <- if (!compareGeom(agb_2017_ornl, ref_raster, stopOnError = FALSE)) {
  resample(agb_2017_ornl, ref_raster, method = "bilinear")
} else { agb_2017_ornl }
agb_2025_resampled <- if (!compareGeom(agb_2025_ornl, ref_raster, stopOnError = FALSE)) {
  resample(agb_2025_ornl, ref_raster, method = "bilinear")
} else { agb_2025_ornl }

# Per-pixel observed AGB at 2017 (independent anchor) and 2025 (validation)
recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))
pt <- recal$pixel_trajectories |> as_tibble()

# Extract 2017 ORNL values at the 50,193 study pixels
pt$x_coord <- xyFromCell(ref_raster, pt$pixel_id)[, 1]
pt$y_coord <- xyFromCell(ref_raster, pt$pixel_id)[, 2]
pt$agb_2017_ornl <- terra::extract(agb_2017_resampled,
                                     as.matrix(pt[, c("x_coord", "y_coord")]))[, 1]
cat("Pixels with 2017 ORNL AGB:", sum(!is.na(pt$agb_2017_ornl)), "\n")
cat("2017 ORNL AGB range:",
    round(range(pt$agb_2017_ornl, na.rm = TRUE), 1), "Mg/ha\n")
cat("2017 ORNL AGB mean:",
    round(mean(pt$agb_2017_ornl, na.rm = TRUE), 1), "Mg/ha\n")

# --- 2. Load Method 11 and 12 raw predictions (already cached) ----------------
m11 <- readRDS(file.path(output_dir, "method11_results.rds"))
m12 <- readRDS(file.path(output_dir, "method12_results.rds"))
r   <- readRDS(file.path(output_dir, "refinements_results.rds"))

# Pull per-pixel raw AGB predictions and metadata from Method 11/12
m11_pred_full <- m11$pred |> as_tibble()
m12_pred_full <- m12$pred |> as_tibble()

# Pixel-level group (stratum × composition) was previously computed; reuse
# from the refinements results
pix_meta <- r$m11_strat$pred |>
  as_tibble() |>
  select(pixel_id, group)

# --- 3. Refit stratified pixel calibration using INDEPENDENT 2017 ORNL anchor -
fit_stratified <- function(pred_df, pix_meta, anchor) {
  d <- pred_df |>
    inner_join(pix_meta, by = "pixel_id") |>
    inner_join(anchor |> select(pixel_id, agb_2017_ornl),
                by = "pixel_id") |>
    filter(agb_raw_2017 > 0, agb_2017_ornl > 0,
            !is.na(group), n_trees >= 2)
  fits <- d |>
    group_by(group) |>
    summarise(n = n(),
               fit = list(lm(agb_2017_ornl ~ agb_raw_2017,
                              data = pick(agb_2017_ornl, agb_raw_2017))),
               .groups = "drop") |>
    mutate(alpha    = map_dbl(fit, ~ coef(.x)[1]),
            beta     = map_dbl(fit, ~ coef(.x)[2]),
            alpha_se = map_dbl(fit,
                                ~ summary(.x)$coefficients[1, "Std. Error"]),
            beta_se  = map_dbl(fit,
                                ~ summary(.x)$coefficients[2, "Std. Error"]),
            r2       = map_dbl(fit, ~ summary(.x)$r.squared))
  # Fallback if any group has unstable fit
  fits <- fits |>
    mutate(use_alpha = ifelse(beta > 0.1 & beta < 20, alpha, 0),
            use_beta  = ifelse(beta > 0.1 & beta < 20, beta, NA_real_))
  list(fits = fits, data = d)
}

# Anchor is the same for both methods
anchor <- pt |> select(pixel_id, agb_2017_ornl)

cat("\n=== Refitting Method 11 with ORNL 2017 anchor ===\n")
m11_strat_ornl <- fit_stratified(m11_pred_full, pix_meta, anchor)
print(m11_strat_ornl$fits |> select(group, n, alpha, alpha_se, beta, beta_se, r2))

cat("\n=== Refitting Method 12 with ORNL 2017 anchor ===\n")
m12_strat_ornl <- fit_stratified(m12_pred_full, pix_meta, anchor)
print(m12_strat_ornl$fits |> select(group, n, alpha, alpha_se, beta, beta_se, r2))

# --- 4. Apply calibration to 2025 predictions and validate against ORNL 2025 --
apply_and_validate <- function(pred_df, fits, pix_meta, pt) {
  pred <- pred_df |>
    inner_join(pix_meta, by = "pixel_id") |>
    left_join(fits |> select(group, alpha = use_alpha, beta = use_beta),
               by = "group") |>
    mutate(agb_pred_2025 = pmax(alpha + beta * agb_raw_2025, 0)) |>
    inner_join(pt |> select(pixel_id, agb_2025_obs = agb_2025),
                by = "pixel_id") |>
    filter(!is.na(agb_pred_2025), !is.na(agb_2025_obs),
            agb_2025_obs > 0, n_trees >= 1)
  stats <- pred |>
    summarise(n = n(),
               bias = mean(agb_pred_2025 - agb_2025_obs),
               bias_pct = 100 * mean(agb_pred_2025 - agb_2025_obs) /
                          mean(agb_2025_obs),
               rmse = sqrt(mean((agb_pred_2025 - agb_2025_obs)^2)),
               rmse_pct = 100 * sqrt(mean((agb_pred_2025 - agb_2025_obs)^2)) /
                          mean(agb_2025_obs),
               r2 = 1 - sum((agb_pred_2025 - agb_2025_obs)^2) /
                        sum((agb_2025_obs - mean(agb_2025_obs))^2),
               mean_pred = mean(agb_pred_2025),
               mean_obs  = mean(agb_2025_obs))
  list(pred = pred, stats = stats)
}

m11_ornl <- apply_and_validate(m11_pred_full, m11_strat_ornl$fits, pix_meta, pt)
m12_ornl <- apply_and_validate(m12_pred_full, m12_strat_ornl$fits, pix_meta, pt)

cat("\n--- Method 11 with ORNL 2017 anchor (2025 validation) ---\n")
print(m11_ornl$stats)

cat("\n--- Method 12 with ORNL 2017 anchor (2025 validation) ---\n")
print(m12_ornl$stats)

# --- 5. Bootstrap 95% CIs on the new validation ------------------------------
boot_rmse <- function(obs, pred, B = 1000, seed = 42) {
  n <- length(obs); set.seed(seed)
  rmses <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    sqrt(mean((pred[idx] - obs[idx])^2))
  })
  rel <- 100 * rmses / mean(obs)
  tibble(rmse_mean = mean(rmses),
          rmse_lo95 = quantile(rmses, 0.025),
          rmse_hi95 = quantile(rmses, 0.975),
          rmse_pct_mean = mean(rel),
          rmse_pct_lo = quantile(rel, 0.025),
          rmse_pct_hi = quantile(rel, 0.975))
}

m11_boot <- boot_rmse(m11_ornl$pred$agb_2025_obs, m11_ornl$pred$agb_pred_2025)
m12_boot <- boot_rmse(m12_ornl$pred$agb_2025_obs, m12_ornl$pred$agb_pred_2025)
cat("\nBootstrap M11:\n"); print(m11_boot)
cat("Bootstrap M12:\n"); print(m12_boot)

# --- 6. Save all outputs ------------------------------------------------------
saveRDS(list(m11 = m11_ornl, m12 = m12_ornl,
              m11_fits = m11_strat_ornl$fits,
              m12_fits = m12_strat_ornl$fits,
              m11_boot = m11_boot, m12_boot = m12_boot,
              anchor = anchor),
        file.path(output_dir, "methods_11_12_ornl_anchor.rds"))

# Refined stats CSV
refined <- bind_rows(
  bind_cols(tibble(method = "11 Hybrid A (ORNL 2017 anchor)"), m11_ornl$stats,
            tibble(rmse_pct_lo = m11_boot$rmse_pct_lo,
                    rmse_pct_hi = m11_boot$rmse_pct_hi)),
  bind_cols(tibble(method = "12 Hybrid B (ORNL 2017 anchor)"), m12_ornl$stats,
            tibble(rmse_pct_lo = m12_boot$rmse_pct_lo,
                    rmse_pct_hi = m12_boot$rmse_pct_hi))
)
write_csv(refined, file.path(output_dir, "methods_11_12_ornl_validation.csv"))

# Stratified fits with SE for Table S9
fits_export <- bind_rows(
  m11_strat_ornl$fits |> mutate(model = "Method 11 (ORNL 2017 anchor)"),
  m12_strat_ornl$fits |> mutate(model = "Method 12 (ORNL 2017 anchor)")
) |>
  select(model, group, n, alpha, alpha_se, beta, beta_se, r2)
write_csv(fits_export, file.path(output_dir, "methods_11_12_ornl_fits.csv"))

# --- 7. Residual diagnostics (FigS14) ----------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
theme_pub <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# Pull Method 9 and 10 holdout data from refinements
cm <- readRDS(file.path(output_dir, "curve_matching_results.rds"))
m9 <- cm$snap |> filter(target_year == 2025) |>
  inner_join(pt |> select(pixel_id, agb_2025), by = "pixel_id") |>
  transmute(method = "9 Curve matching",
             obs = agb_2025, pred = agb_pred,
             resid = pred - obs)

m10 <- r$m10_test |>
  transmute(method = "10 Pure Kohek (holdout)",
             obs = agb_2025_obs, pred = pred_hold,
             resid = pred - obs)

m11r <- m11_ornl$pred |>
  transmute(method = "11 Hybrid A (ORNL 2017 anchor)",
             obs = agb_2025_obs, pred = agb_pred_2025,
             resid = pred - obs)

m12r <- m12_ornl$pred |>
  transmute(method = "12 Hybrid B (ORNL 2017 anchor)",
             obs = agb_2025_obs, pred = agb_pred_2025,
             resid = pred - obs)

resid_df <- bind_rows(m9, m10, m11r, m12r) |>
  mutate(method = factor(method, levels = c(
    "9 Curve matching", "10 Pure Kohek (holdout)",
    "11 Hybrid A (ORNL 2017 anchor)", "12 Hybrid B (ORNL 2017 anchor)")))

# Panel 1: residuals vs predicted
p_resid <- ggplot(resid_df, aes(pred, resid)) +
  geom_hex(bins = 30, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "blue",
               linewidth = 0.5, span = 0.5) +
  scale_fill_viridis_c(option = "B", trans = "log10") +
  facet_wrap(~ method, nrow = 2) +
  labs(x = "Predicted 2025 AGB (Mg ha^-1)",
        y = "Residual (predicted - observed, Mg ha^-1)",
        title = "Figure S14a. Residual vs. predicted 2025 AGB",
        subtitle = "Hex density; red = zero residual; blue = LOESS smoother") +
  theme_pub +
  theme(strip.text = element_text(face = "bold", size = 9))

# Panel 2: Q-Q plot
qq_df <- resid_df |>
  group_by(method) |>
  mutate(theoretical = qnorm(rank(resid) / (n() + 1))) |>
  ungroup()

p_qq <- ggplot(qq_df, aes(theoretical, scale(resid)[, 1])) +
  geom_point(alpha = 0.05, size = 0.3) +
  geom_abline(color = "red", linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ method, nrow = 1) +
  labs(x = "Theoretical normal quantile",
        y = "Standardized residual",
        title = "Figure S14b. Normal Q-Q of residuals",
        subtitle = "Heavy tails (departure from red line) indicate non Gaussian residual structure") +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-6, 6)) +
  theme_pub +
  theme(strip.text = element_text(face = "bold", size = 9))

figS14 <- p_resid / p_qq +
  plot_layout(heights = c(1.6, 1)) +
  plot_annotation(theme = theme(
    plot.title = element_text(face = "bold", size = 13)))

fig_dir <- file.path(output_dir, "manuscript_figures_v45")
ggsave(file.path(fig_dir, "FigS14_residuals_methods9_12.png"), figS14,
        width = 24, height = 22, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "FigS14_residuals_methods9_12.pdf"), figS14,
        width = 24, height = 22, units = "cm")
cat("\nSaved FigS14\n")

cat("\n=== Methods 11/12 ORNL anchor rerun + FigS14 complete ===\n")
