# ==============================================================================
# HOWLAND FOREST: METHOD 12
# Directional 3D projection HYBRID B:
#   ITC 2017 anchor + TLS-calibrated H-D and Kohek coefficients
# ==============================================================================
# Same pipeline as Method 11 but substitutes local TLS-derived H-D allometry
# and fits M0, M1 from paired (2017 ALS, 2025 TLS) crown radii where trees
# can be matched by proximity. Falls back to Method 11 defaults where TLS
# data are insufficient.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(FNN)
})

hrf_root   <- Sys.getenv("HRF_ROOT",   unset = "/users/PUOM0008/crsfaaron/HRF")
output_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path(hrf_root, "output"))
field_dir  <- Sys.getenv("FIELD_DIR",  unset = file.path(hrf_root, "data/field_inventory"))
tls_dir    <- Sys.getenv("TLS_DIR",
                          unset = file.path(hrf_root, "data/tls_2025"))

ANCHOR_YR <- 2017
HORIZONS  <- c(2025, 2035, 2045, 2075)
# Matching threshold for 2017 ALS to 2025 TLS (generous since trees grow
# and move very little)
MATCH_TOL_M <- 3.0

cat("=== Method 12: ITC + TLS calibrated directional projection ===\n")

# --- Reference grid and pixel observations -----------------------------------
ref_raster <- rast(file.path(output_dir, "agb_2025_melite_30m.tif"))
pixel_area_ha <- prod(res(ref_raster)) / 10000

recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))
pt <- recal$pixel_trajectories |> as_tibble()
cell_xy <- xyFromCell(ref_raster, pt$pixel_id)
pt$x <- cell_xy[, 1]; pt$y <- cell_xy[, 2]

# --- ITC 2017 ---------------------------------------------------------------
itc <- st_read(file.path(output_dir, "Howland_ITC_2017.gpkg"), quiet = TRUE)
itc <- itc |> filter(Z >= 2, Z <= 45)
itc_xy <- st_coordinates(itc)[, 1:2]
itc_df <- itc |>
  st_drop_geometry() |>
  mutate(X = itc_xy[,1], Y = itc_xy[,2]) |>
  mutate(pixel_id = cellFromXY(ref_raster,
                                 as.matrix(data.frame(X,Y)))) |>
  filter(!is.na(pixel_id)) |>
  inner_join(pt |> select(pixel_id, agb_class, agb_2012,
                            agb_2015, agb_2025),
              by = "pixel_id")
cat("ITC trees inside Howland study pixels:", nrow(itc_df), "\n")

# --- TLS 2025 tree list -----------------------------------------------------
tls <- read.table(file.path(tls_dir, "Merged_data_total.txt"),
                   header = FALSE, skip = 1) |>
  setNames(c("X","Y","Z","TileID","TreeID","GlobalTreeID",
             "Height_m","DBH_cm","conf_loc","conf_DBH"))
tls <- tls |>
  filter(Height_m >= 2, Height_m <= 45,
          DBH_cm > 0, DBH_cm <= 150,
          conf_DBH >= 50)
cat("TLS trees with DBH:", nrow(tls), "\n")

# --- Step 1: TLS calibrated H-D allometry -----------------------------------
# Fit DBH = a * Ht^b via nls
hd_fit <- tryCatch(
  nls(DBH_cm ~ a * Height_m^b, data = tls,
      start = list(a = 1.0, b = 1.2)),
  error = function(e) NULL)
if (!is.null(hd_fit)) {
  hd_a <- coef(hd_fit)["a"]; hd_b <- coef(hd_fit)["b"]
  cat("TLS H-D power fit: DBH =", round(hd_a, 3), "* Ht^",
      round(hd_b, 3), "\n")
} else {
  hd_a <- 0.35; hd_b <- 1.30
  cat("TLS H-D fit failed; using default a =", hd_a, "b =", hd_b, "\n")
}

itc_df$dbh_2017 <- pmax(hd_a * itc_df$Z^hd_b, 1)

# --- Step 2: match 2017 ALS to 2025 TLS to calibrate Kohek M0/M1 ------------
tls_xy <- as.matrix(tls[, c("X", "Y")])
itc_in_tls_bbox <- itc_df |>
  filter(X >= min(tls$X) - 5, X <= max(tls$X) + 5,
          Y >= min(tls$Y) - 5, Y <= max(tls$Y) + 5)
cat("ITC trees in TLS footprint:", nrow(itc_in_tls_bbox), "\n")

if (nrow(itc_in_tls_bbox) > 10) {
  nn <- get.knnx(tls_xy,
                  as.matrix(itc_in_tls_bbox[, c("X","Y")]), k = 1)
  paired <- itc_in_tls_bbox |>
    mutate(tls_dist = nn$nn.dist[,1],
            tls_idx  = nn$nn.index[,1]) |>
    filter(tls_dist <= MATCH_TOL_M)
  paired$tls_dbh_2025  <- tls$DBH_cm[paired$tls_idx]
  paired$tls_ht_2025   <- tls$Height_m[paired$tls_idx]
  # TLS crown radius approximated from DBH via species-agnostic
  # allometric (Condes 2005 for mixed forests)
  paired$tls_cr_2025   <- 0.4 * paired$tls_dbh_2025^0.6
  # 2017 ALS crown radius already in paired$crown_radius_m
  paired <- paired |>
    mutate(growth_ratio = tls_cr_2025 / crown_radius_m,
            ht_ratio     = tls_ht_2025 / Z) |>
    filter(is.finite(growth_ratio), growth_ratio > 0.5,
            growth_ratio < 3)
  cat("Calibration pairs:", nrow(paired), "\n")

  if (nrow(paired) >= 20) {
    # Fit Kohek: growth_ratio ~ 1 + (M0 + M1 * shade) * 8yr
    # Approximate shade from 5 m neighborhood height already in ITC
    r_ht <- rast(ext(ref_raster), resolution = 2, crs = crs(ref_raster))
    pts <- vect(as.matrix(itc_df[, c("X","Y")]), type = "points",
                 atts = data.frame(z = itc_df$Z),
                 crs = crs(ref_raster))
    r_ht <- rasterize(pts, r_ht, field = "z", fun = "max")
    r_fmax <- focal(r_ht, w = matrix(1, 3, 3), fun = "max", na.rm = TRUE)
    paired$neigh_max <- terra::extract(r_fmax,
      as.matrix(paired[, c("X","Y")]))[, 1]
    paired$neigh_max[is.na(paired$neigh_max)] <- paired$Z[is.na(paired$neigh_max)]
    paired$shade <- pmin(pmax((paired$neigh_max - paired$Z) /
                               pmax(paired$Z, 5), 0), 1)

    # (growth_ratio - 1) / 8 = M0 + M1 * shade
    paired$y_lhs <- pmax((paired$growth_ratio - 1) / 8, -0.1)
    kohek_fit <- lm(y_lhs ~ shade, data = paired)
    M0_fit <- coef(kohek_fit)[1]
    M1_fit <- coef(kohek_fit)[2]
    cat("TLS calibrated Kohek: M0 =", round(M0_fit, 3),
        ", M1 =", round(M1_fit, 3), "\n")
    if (!is.finite(M0_fit) || M0_fit < 0) M0_fit <- 0.02
    if (!is.finite(M1_fit))               M1_fit <- 0.05
  } else {
    cat("Too few pairs; using defaults\n")
    M0_fit <- 0.02; M1_fit <- 0.05
  }
} else {
  M0_fit <- 0.02; M1_fit <- 0.05
  paired <- NULL
}

# --- Step 3: height and DBH projection (as in Method 11) --------------------
cm <- readRDS(file.path(output_dir, "curve_matching_library.rds")) |>
  as_tibble()

plots <- read_csv(file.path(field_dir, "Plots_Data.csv"),
                   show_col_types = FALSE) |>
  rename(plot = Plot, northing = Northing, easting = Easting)
trees_fld <- read_csv(file.path(field_dir, "Tree_Data.csv"),
                       show_col_types = FALSE) |>
  rename(plot_tree = Plot, species = Species, dbh_cm = DBH_cm,
          snag = Snag) |>
  mutate(plot = stringr::str_sub(plot_tree, 1, 6))
pct_con <- trees_fld |>
  filter(is.na(snag) | snag == 0) |>
  group_by(plot) |>
  summarise(pct_con = mean(species %in% c("PIRU","ABBA","PIMA","TSCA",
                                            "PIST","PIGL","THOC","LALA")),
            .groups = "drop")
plots <- plots |>
  left_join(pct_con, by = "plot") |>
  mutate(spp_comp = case_when(
    pct_con >= 0.75 ~ "coniferous",
    pct_con <= 0.25 ~ "deciduous",
    TRUE            ~ "mixed"))
plot_coords <- plots |> filter(!is.na(easting)) |>
  select(easting, northing) |> as.matrix()
nnp <- get.knnx(plot_coords,
                 as.matrix(itc_df[, c("X","Y")]), k = 1)
itc_df$spp_comp <- plots |> filter(!is.na(easting)) |>
  pull(spp_comp) |>
  (\(v) v[nnp$nn.index[,1]])()
itc_df$stratum3    <- itc_df$agb_class
itc_df$template_id <- paste(itc_df$stratum3, itc_df$spp_comp, sep = "_")

anchor_offset  <- ANCHOR_YR - 2015
horizon_offset <- HORIZONS   - 2015
tpl <- cm |>
  filter(year_offset %in% c(anchor_offset, horizon_offset)) |>
  group_by(template_id, year_offset) |>
  summarise(top_ht = mean(top_ht_m, na.rm = TRUE),
            .groups = "drop") |>
  pivot_wider(id_cols = template_id,
              names_from = year_offset,
              values_from = top_ht,
              names_prefix = "topht_")

itc_df <- itc_df |> left_join(tpl, by = "template_id")
anchor_col <- paste0("topht_", anchor_offset)
for (h in HORIZONS) {
  offs <- h - 2015
  scale_col <- paste0("topht_", offs)
  itc_df[[paste0("ht_", h)]] <- pmin(
    itc_df$Z * itc_df[[scale_col]] / itc_df[[anchor_col]],
    45)
  itc_df[[paste0("dbh_", h)]] <- pmax(hd_a *
    itc_df[[paste0("ht_", h)]]^hd_b, 1)
}

# --- Per tree AGB (Chojnacky) -----------------------------------------------
itc_df$agb_kg_2017 <- exp(-2.0773 + 2.3323 * log(itc_df$dbh_2017))
for (h in HORIZONS) {
  itc_df[[paste0("agb_kg_", h)]] <- exp(-2.0773 + 2.3323 *
    log(itc_df[[paste0("dbh_", h)]]))
}

# --- Pixel sum and calibration ----------------------------------------------
raw <- itc_df |>
  as_tibble() |>
  group_by(pixel_id) |>
  summarise(agb_raw_2017 = sum(agb_kg_2017, na.rm = TRUE) / 1000 / pixel_area_ha,
            agb_raw_2025 = sum(agb_kg_2025, na.rm = TRUE) / 1000 / pixel_area_ha,
            agb_raw_2035 = sum(agb_kg_2035, na.rm = TRUE) / 1000 / pixel_area_ha,
            agb_raw_2045 = sum(agb_kg_2045, na.rm = TRUE) / 1000 / pixel_area_ha,
            agb_raw_2075 = sum(agb_kg_2075, na.rm = TRUE) / 1000 / pixel_area_ha,
            n_trees      = n(),
            .groups = "drop")

cal_df <- raw |>
  inner_join(pt |> select(pixel_id, agb_2015, agb_2025), by = "pixel_id") |>
  mutate(agb_obs_2017 = agb_2015 + (agb_2025 - agb_2015) * 2 / 10) |>
  drop_na(agb_raw_2017, agb_obs_2017) |>
  filter(agb_raw_2017 > 0, agb_obs_2017 > 0, n_trees >= 2)
cal_fit <- lm(agb_obs_2017 ~ agb_raw_2017, data = cal_df)
alpha <- coef(cal_fit)[1]; beta <- coef(cal_fit)[2]
if (!is.finite(beta) || beta <= 0) {
  beta <- sum(cal_df$agb_obs_2017) / sum(cal_df$agb_raw_2017)
  alpha <- 0
}
cat("Pixel calibration: alpha =", round(alpha, 1),
    ", beta =", round(beta, 3), "\n")

pred <- raw
for (h in HORIZONS) {
  src <- paste0("agb_raw_", h); dst <- paste0("agb_pred_", h)
  pred[[dst]] <- pmax(alpha + beta * pred[[src]], 0)
}

val <- pred |>
  inner_join(pt |> select(pixel_id, agb_2025_obs = agb_2025),
              by = "pixel_id") |>
  drop_na(agb_pred_2025, agb_2025_obs) |>
  filter(n_trees >= 1, agb_2025_obs > 0)
val_stats <- val |>
  summarise(n = n(),
            bias = mean(agb_pred_2025 - agb_2025_obs),
            bias_pct = 100 * mean(agb_pred_2025 - agb_2025_obs) /
                        mean(agb_2025_obs),
            rmse = sqrt(mean((agb_pred_2025 - agb_2025_obs)^2)),
            rmse_pct = 100 * sqrt(mean((agb_pred_2025 - agb_2025_obs)^2)) /
                        mean(agb_2025_obs),
            r2 = 1 - sum((agb_pred_2025 - agb_2025_obs)^2) /
                      sum((agb_2025_obs - mean(agb_2025_obs))^2))
cat("\n--- Method 12 validation (2025) ---\n")
print(val_stats)
write_csv(val_stats, file.path(output_dir, "method12_validation.csv"))
saveRDS(list(pred = pred, val = val_stats, val_pixels = val,
              alpha = alpha, beta = beta,
              hd_a = hd_a, hd_b = hd_b,
              M0_tls = M0_fit, M1_tls = M1_fit,
              paired = paired),
        file.path(output_dir, "method12_results.rds"))

cat("\n=== Method 12 complete ===\n")
