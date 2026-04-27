# ==============================================================================
# POST-HOC REFINEMENTS for Methods 9, 10, 11, 12
# P1: Stratified calibration for Methods 11 and 12 (by stratum × composition)
# P2: 80/20 train/test holdout for Method 10
# P3: Bootstrap 95% CIs for RMSE of all four methods
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(FNN)
})

hrf_root   <- "/users/PUOM0008/crsfaaron/HRF"
output_dir <- file.path(hrf_root, "output")

# --- Shared inputs -----------------------------------------------------------
recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))
pt <- recal$pixel_trajectories |> as_tibble()
ref_raster <- rast(file.path(output_dir, "agb_2025_melite_30m.tif"))

# Build pixel level species composition from ITC 2017 trees
itc <- st_read(file.path(output_dir, "Howland_ITC_2017.gpkg"), quiet = TRUE) |>
  filter(Z >= 2, Z <= 45)
itc_xy <- st_coordinates(itc)[, 1:2]
itc_df <- itc |> st_drop_geometry() |>
  mutate(X = itc_xy[,1], Y = itc_xy[,2]) |>
  mutate(pixel_id = cellFromXY(ref_raster,
                                 as.matrix(data.frame(X,Y)))) |>
  filter(!is.na(pixel_id))

# Species comp from field plots
plots <- read_csv(file.path(hrf_root, "data/field_inventory/Plots_Data.csv"),
                   show_col_types = FALSE) |>
  rename(plot = Plot, northing = Northing, easting = Easting)
trees_fld <- read_csv(file.path(hrf_root, "data/field_inventory/Tree_Data.csv"),
                       show_col_types = FALSE) |>
  rename(plot_tree = Plot, species = Species, dbh_cm = DBH_cm, snag = Snag) |>
  mutate(plot = stringr::str_sub(plot_tree, 1, 6))
pct_con <- trees_fld |>
  filter(is.na(snag) | snag == 0) |>
  group_by(plot) |>
  summarise(pct_con = mean(species %in% c("PIRU","ABBA","PIMA","TSCA",
                                            "PIST","PIGL","THOC","LALA")),
            .groups = "drop")
plots <- plots |> left_join(pct_con, by = "plot") |>
  mutate(spp_comp = case_when(pct_con >= 0.75 ~ "coniferous",
                                pct_con <= 0.25 ~ "deciduous",
                                TRUE ~ "mixed"))
plot_xy <- plots |> filter(!is.na(easting)) |>
  select(easting, northing) |> as.matrix()
nnp <- get.knnx(plot_xy, as.matrix(itc_df[, c("X","Y")]), k = 1)
itc_df$spp_comp <- plots |> filter(!is.na(easting)) |>
  pull(spp_comp) |> (\(v) v[nnp$nn.index[,1]])()

# Pixel level dominant composition (mode)
pix_comp <- itc_df |>
  group_by(pixel_id) |>
  summarise(spp_comp = names(sort(table(spp_comp), decreasing = TRUE))[1],
            .groups = "drop")

# Pixel level stratum from observations
pix_stratum <- pt |> select(pixel_id, stratum3 = agb_class)

pix_meta <- pix_stratum |>
  inner_join(pix_comp, by = "pixel_id") |>
  mutate(group = paste(stratum3, spp_comp, sep = "_"))

# --- Stratified calibration helper -------------------------------------------
stratified_calibrate <- function(res_name, raw_2017_col = "agb_raw_2017",
                                  anchor_year_for_obs = 2017) {
  res <- readRDS(file.path(output_dir, paste0(res_name, "_results.rds")))
  pred <- res$pred |> as_tibble() |>
    inner_join(pix_meta, by = "pixel_id") |>
    inner_join(pt |> select(pixel_id, agb_2015, agb_2025),
                by = "pixel_id")

  # Observed 2017 via linear interpolation between 2015 and 2025
  pred <- pred |>
    mutate(agb_obs_2017 = agb_2015 + (agb_2025 - agb_2015) * 2 / 10)

  # Fit one calibration per group (stratum × composition)
  # Constrain beta > 0 and reasonable alpha
  fits <- pred |>
    filter(.data[[raw_2017_col]] > 0, agb_obs_2017 > 0, n_trees >= 2) |>
    group_by(group) |>
    summarise(
      n = n(),
      fit = list(lm(agb_obs_2017 ~ .data[[raw_2017_col]], data = cur_data())),
      .groups = "drop"
    ) |>
    mutate(
      alpha = map_dbl(fit, ~ coef(.x)[1]),
      beta  = map_dbl(fit, ~ coef(.x)[2]),
      # Fallback to ratio through origin if fit unstable
      beta_fallback = map_dbl(fit, ~ {
        d <- .x$model; sum(d[[1]]) / sum(d[[2]])
      })
    )
  # Where beta is non-positive or extreme, fall back to ratio
  fits <- fits |>
    mutate(alpha = ifelse(is.finite(alpha) & beta > 0.1 & beta < 20,
                            alpha, 0),
           beta  = ifelse(is.finite(beta) & beta > 0.1 & beta < 20,
                            beta, beta_fallback))

  cat("\nStratified calibration for", res_name, ":\n")
  print(fits |> select(group, n, alpha, beta))

  # Apply to 2025 predictions
  pred <- pred |> left_join(fits |> select(group, alpha, beta),
                              by = "group")
  raw_2025 <- "agb_raw_2025"
  pred$agb_pred_2025_strat <- pmax(pred$alpha + pred$beta *
                                     pred[[raw_2025]], 0)

  # Validate
  val <- pred |>
    filter(!is.na(agb_pred_2025_strat), agb_2025 > 0, n_trees >= 1)
  stats <- val |>
    summarise(
      n = n(),
      bias = mean(agb_pred_2025_strat - agb_2025),
      bias_pct = 100 * mean(agb_pred_2025_strat - agb_2025) /
                  mean(agb_2025),
      rmse = sqrt(mean((agb_pred_2025_strat - agb_2025)^2)),
      rmse_pct = 100 * sqrt(mean((agb_pred_2025_strat - agb_2025)^2)) /
                  mean(agb_2025),
      r2 = 1 - sum((agb_pred_2025_strat - agb_2025)^2) /
               sum((agb_2025 - mean(agb_2025))^2)
    )
  cat("Stratified validation:\n"); print(stats)
  list(pred = pred, stats = stats, fits = fits)
}

cat("\n==== P1: Stratified calibration for Methods 11 and 12 ====\n")
m11_strat <- stratified_calibrate("method11")
m12_strat <- stratified_calibrate("method12")

# --- P2: 80/20 train/test for Method 10 --------------------------------------
cat("\n==== P2: Honest 80/20 holdout for Method 10 ====\n")
m10 <- readRDS(file.path(output_dir, "method10_results.rds"))
m10_val <- m10$val_pixels |> as_tibble()

set.seed(2026)
m10_val$split <- ifelse(runif(nrow(m10_val)) < 0.8, "train", "test")
train <- m10_val |> filter(split == "train", agb_raw_2025 > 0,
                             agb_2025_obs > 0, n_trees >= 2)
test  <- m10_val |> filter(split == "test",  agb_raw_2025 > 0,
                             agb_2025_obs > 0)

cat("Train n:", nrow(train), " Test n:", nrow(test), "\n")
cal_fit <- lm(agb_2025_obs ~ agb_raw_2025, data = train)
alpha_h <- coef(cal_fit)[1]; beta_h <- coef(cal_fit)[2]
cat("Holdout alpha =", round(alpha_h,1), " beta =", round(beta_h,3), "\n")

test$pred_hold <- pmax(alpha_h + beta_h * test$agb_raw_2025, 0)
m10_hold_stats <- test |>
  summarise(n = n(),
            bias = mean(pred_hold - agb_2025_obs),
            bias_pct = 100 * mean(pred_hold - agb_2025_obs) /
                        mean(agb_2025_obs),
            rmse = sqrt(mean((pred_hold - agb_2025_obs)^2)),
            rmse_pct = 100 * sqrt(mean((pred_hold - agb_2025_obs)^2)) /
                        mean(agb_2025_obs),
            r2 = 1 - sum((pred_hold - agb_2025_obs)^2) /
                     sum((agb_2025_obs - mean(agb_2025_obs))^2))
cat("Holdout test stats:\n"); print(m10_hold_stats)

# --- P3: Bootstrap 95% CIs for RMSE of all four methods ----------------------
cat("\n==== P3: Bootstrap 95% CIs for RMSE (1000 iterations) ====\n")

boot_rmse <- function(obs, pred, B = 1000, seed = 42) {
  stopifnot(length(obs) == length(pred))
  n <- length(obs)
  set.seed(seed)
  rmses <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    sqrt(mean((pred[idx] - obs[idx])^2))
  })
  rel <- 100 * rmses / mean(obs)
  tibble(
    rmse_mean   = mean(rmses),
    rmse_lo95   = quantile(rmses, 0.025),
    rmse_hi95   = quantile(rmses, 0.975),
    rmse_pct_mean = mean(rel),
    rmse_pct_lo   = quantile(rel, 0.025),
    rmse_pct_hi   = quantile(rel, 0.975)
  )
}

# Method 9 — curve matching
cm <- readRDS(file.path(output_dir, "curve_matching_results.rds"))
cm25 <- cm$snap |> filter(target_year == 2025) |>
  inner_join(pt |> select(pixel_id, agb_2025), by = "pixel_id") |>
  drop_na(agb_pred, agb_2025)
m9_boot <- boot_rmse(cm25$agb_2025, cm25$agb_pred)

# Method 10 — use holdout test
m10_boot <- boot_rmse(test$agb_2025_obs, test$pred_hold)

# Method 11 — stratified calibration result
m11_boot <- boot_rmse(m11_strat$pred$agb_2025,
                        m11_strat$pred$agb_pred_2025_strat)

# Method 12 — stratified calibration result
m12_boot <- boot_rmse(m12_strat$pred$agb_2025,
                        m12_strat$pred$agb_pred_2025_strat)

summary_tbl <- bind_rows(
  tibble(method = "9 Curve matching",            m9_boot),
  tibble(method = "10 Pure Kohek (holdout)",     m10_boot),
  tibble(method = "11 Hybrid A (stratified)",    m11_boot),
  tibble(method = "12 Hybrid B (stratified)",    m12_boot)
)
cat("\nBootstrap RMSE summary:\n"); print(summary_tbl)

write_csv(summary_tbl, file.path(output_dir, "refinements_summary.csv"))

# Full refined validation stats
refined_stats <- bind_rows(
  tibble(method = "9 Curve matching",
         n = nrow(cm25),
         bias = mean(cm25$agb_pred - cm25$agb_2025),
         bias_pct = 100 * mean(cm25$agb_pred - cm25$agb_2025) /
                    mean(cm25$agb_2025),
         rmse = sqrt(mean((cm25$agb_pred - cm25$agb_2025)^2)),
         rmse_pct = 100 * sqrt(mean((cm25$agb_pred - cm25$agb_2025)^2)) /
                    mean(cm25$agb_2025),
         r2 = 1 - sum((cm25$agb_pred - cm25$agb_2025)^2) /
              sum((cm25$agb_2025 - mean(cm25$agb_2025))^2)),
  bind_cols(tibble(method = "10 Pure Kohek (holdout)"), m10_hold_stats),
  bind_cols(tibble(method = "11 Hybrid A (stratified)"), m11_strat$stats),
  bind_cols(tibble(method = "12 Hybrid B (stratified)"), m12_strat$stats)
) |>
  left_join(summary_tbl |> select(method, rmse_pct_lo, rmse_pct_hi),
             by = "method")

cat("\n==== REFINED METHOD STATS ====\n")
print(refined_stats)
write_csv(refined_stats, file.path(output_dir, "refined_method_stats.csv"))

saveRDS(list(m9_boot = m9_boot, m10_boot = m10_boot,
              m10_hold = m10_hold_stats,
              m10_test = test,
              m11_strat = m11_strat, m12_strat = m12_strat,
              refined_stats = refined_stats),
        file.path(output_dir, "refinements_results.rds"))

cat("\n=== Refinements complete ===\n")
