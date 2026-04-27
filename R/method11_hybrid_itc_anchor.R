# ==============================================================================
# HOWLAND FOREST: METHOD 11 (v2)
# Directional 3D projection HYBRID A: ITC 2017 anchor + default Kohek M0/M1
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

M0        <- 0.2
M1        <- 0.9
ANCHOR_YR <- 2017
HORIZONS  <- c(2025, 2035, 2045, 2075)

cat("=== Method 11: Hybrid directional projection (ITC anchor) ===\n")

# --- Load reference 10 m AGB grid as the target grid --------------------------
ref_raster <- rast(file.path(output_dir, "agb_2025_melite_30m.tif"))
cat("Reference grid: ", ncell(ref_raster), "cells, res =",
    res(ref_raster)[1], "m\n")
pixel_area_ha <- prod(res(ref_raster)) / 10000
cat("Pixel area:", round(pixel_area_ha, 4), "ha\n")

# Per pixel observed trajectories
recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))
pt <- recal$pixel_trajectories |> as_tibble()

# Attach coordinates by pixel_id (pixel_id is the raster cell index)
cell_xy <- xyFromCell(ref_raster, pt$pixel_id)
pt$x <- cell_xy[, 1]
pt$y <- cell_xy[, 2]
cat("Pixels with coords:", nrow(pt), "\n")

# --- Load ITC 2017 trees ------------------------------------------------------
itc <- st_read(file.path(output_dir, "Howland_ITC_2017.gpkg"), quiet = TRUE)
itc <- itc |>
  filter(Z >= 2, Z <= 45) |>
  mutate(xy = st_coordinates(geom))
itc_xy <- st_coordinates(itc)[, 1:2]
itc_df <- itc |>
  st_drop_geometry() |>
  mutate(X = itc_xy[,1], Y = itc_xy[,2])
cat("ITC 2017 trees after height filter:", nrow(itc_df), "\n")

# --- Assign each ITC tree to a raster cell by exact pixel bucketing -----------
# terra::cellFromXY requires matching CRS; ITC is WGS84 UTM 19N,
# raster is NAD83(2011) UTM 19N. Difference is ~1 m at most; ignore.
itc_df$pixel_id <- cellFromXY(ref_raster,
                               as.matrix(itc_df[, c("X", "Y")]))
itc_df <- itc_df |> filter(!is.na(pixel_id))
cat("ITC trees inside reference grid:", nrow(itc_df), "\n")

# Keep only trees whose pixel has observed 2017 AGB in pt
itc_df <- itc_df |> inner_join(pt |> select(pixel_id, agb_2012,
                                              agb_2015, agb_2025,
                                              agb_class),
                                by = "pixel_id")
cat("ITC trees inside Howland study pixels:", nrow(itc_df), "\n")

# --- Field plots for species composition assignment --------------------------
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

# Stratum from agb_class (Low/Medium/High)
itc_df$stratum3 <- itc_df$agb_class
itc_df$template_id <- paste(itc_df$stratum3, itc_df$spp_comp, sep = "_")

# --- Curve matching library for growth trajectories --------------------------
cm <- readRDS(file.path(output_dir, "curve_matching_library.rds")) |>
  as_tibble()

anchor_offset  <- ANCHOR_YR - 2015   # 2
horizon_offset <- HORIZONS   - 2015
tpl <- cm |>
  filter(year_offset %in% c(anchor_offset, horizon_offset)) |>
  select(template_id, year_offset, top_ht_m, qmd_cm, agb_Mg_ha) |>
  group_by(template_id, year_offset) |>
  summarise(top_ht = mean(top_ht_m, na.rm = TRUE),
            qmd    = mean(qmd_cm, na.rm = TRUE),
            agb    = mean(agb_Mg_ha, na.rm = TRUE),
            .groups = "drop")

tpl_wide <- tpl |>
  pivot_wider(id_cols = template_id,
              names_from = year_offset,
              values_from = c(top_ht, qmd, agb))

# Height growth ratios: ht_future = ht_2017 * ratio
itc_df <- itc_df |> left_join(tpl_wide, by = "template_id")
# --- Per tree DBH imputation (literature power H D) --------------------------
# DBH (cm) = 1.3 + 0.8 * (Ht - 1.3)^1.1 (approximate inverse of Chapman
# Richards H-D; simple power fit that gives sensible values across the
# H range of interest)
itc_df$dbh_2017 <- 1.3 + 0.8 * pmax(itc_df$Z - 1.3, 0)^1.1

for (h in HORIZONS) {
  offs <- h - 2015
  scale_col <- paste0("top_ht_", offs)
  anchor_col <- paste0("top_ht_", anchor_offset)
  itc_df[[paste0("ht_", h)]] <- pmin(
    itc_df$Z * itc_df[[scale_col]] / itc_df[[anchor_col]],
    45
  )
  itc_df[[paste0("dbh_", h)]] <- 1.3 + 0.8 *
    pmax(itc_df[[paste0("ht_", h)]] - 1.3, 0)^1.1
}

# --- Per tree AGB via Chojnacky 2014 generalized -----------------------------
itc_df$agb_kg_2017 <- exp(-2.0773 + 2.3323 * log(pmax(itc_df$dbh_2017, 1)))
for (h in HORIZONS) {
  itc_df[[paste0("agb_kg_", h)]] <- exp(-2.0773 + 2.3323 *
    log(pmax(itc_df[[paste0("dbh_", h)]], 1)))
}

# --- Kohek crown expansion proxy (simple neighbor shading) -------------------
# Compute 5 m neighborhood max canopy height via a raster focal operation
r_ht <- rast(ext(ref_raster), resolution = 2, crs = crs(ref_raster))
pts <- vect(as.matrix(itc_df[, c("X","Y")]), type = "points",
             atts = data.frame(z = itc_df$Z), crs = crs(ref_raster))
r_ht <- rasterize(pts, r_ht, field = "z", fun = "max")
r_fmax <- focal(r_ht, w = matrix(1, 3, 3), fun = "max", na.rm = TRUE)
itc_df$neigh_max <- terra::extract(r_fmax, as.matrix(itc_df[, c("X","Y")]))[, 1]
itc_df$neigh_max[is.na(itc_df$neigh_max)] <- itc_df$Z[is.na(itc_df$neigh_max)]
itc_df$shade <- pmin(pmax((itc_df$neigh_max - itc_df$Z) /
                           pmax(itc_df$Z, 5), 0), 1)
# Apply expansion to crown radius (not used here for AGB, but saved for
# pixel attribution in Method 12 variants)
for (h in HORIZONS) {
  yrs <- h - ANCHOR_YR
  g <- 1 + (M0 + M1 * itc_df$shade) / 50
  itc_df[[paste0("r_", h)]] <- itc_df$crown_radius_m *
    pmin(g^yrs, 2.0)
}

# --- Sum per tree AGB into raster pixels -------------------------------------
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

# Quick sanity: distribution of per pixel raw AGB
cat("\nRaw 2017 ITC-derived AGB summary (Mg/ha):\n")
print(summary(raw$agb_raw_2017))

# --- Pixel level linear calibration -------------------------------------------
cal_df <- raw |>
  inner_join(pt |> select(pixel_id, agb_2012, agb_2015, agb_2025),
              by = "pixel_id") |>
  # "observed 2017" linearly interpolated between 2015 and 2025
  mutate(agb_obs_2017 = agb_2015 + (agb_2025 - agb_2015) * 2 / 10) |>
  drop_na(agb_raw_2017, agb_obs_2017) |>
  filter(agb_raw_2017 > 0, agb_obs_2017 > 0, n_trees >= 2)
cat("Calibration pixels (n):", nrow(cal_df), "\n")

cal_fit <- lm(agb_obs_2017 ~ agb_raw_2017, data = cal_df)
cat("Calibration coefficients:\n"); print(summary(cal_fit)$coefficients)
alpha <- coef(cal_fit)[1]; beta <- coef(cal_fit)[2]
if (!is.finite(beta) || beta <= 0) {
  beta <- sum(cal_df$agb_obs_2017) / sum(cal_df$agb_raw_2017)
  alpha <- 0
}
cat("Using: alpha =", round(alpha, 1), ", beta =", round(beta, 3), "\n")

# Apply
pred <- raw
for (h in HORIZONS) {
  src <- paste0("agb_raw_", h)
  dst <- paste0("agb_pred_", h)
  pred[[dst]] <- pmax(alpha + beta * pred[[src]], 0)
}

# --- Validation ---------------------------------------------------------------
val <- pred |>
  inner_join(pt |> select(pixel_id, agb_2025_obs = agb_2025),
              by = "pixel_id") |>
  drop_na(agb_pred_2025, agb_2025_obs) |>
  filter(n_trees >= 1, agb_2025_obs > 0)

val_stats <- val |>
  summarise(
    n        = n(),
    bias     = mean(agb_pred_2025 - agb_2025_obs),
    bias_pct = 100 * mean(agb_pred_2025 - agb_2025_obs) /
                mean(agb_2025_obs),
    rmse     = sqrt(mean((agb_pred_2025 - agb_2025_obs)^2)),
    rmse_pct = 100 * sqrt(mean((agb_pred_2025 - agb_2025_obs)^2)) /
                mean(agb_2025_obs),
    r2       = 1 - sum((agb_pred_2025 - agb_2025_obs)^2) /
                    sum((agb_2025_obs - mean(agb_2025_obs))^2)
  )

cat("\n--- Method 11 validation (2025) ---\n")
print(val_stats)
write_csv(val_stats, file.path(output_dir, "method11_validation.csv"))
saveRDS(list(pred = pred, val = val_stats, val_pixels = val,
              alpha = alpha, beta = beta),
        file.path(output_dir, "method11_results.rds"))

cat("\n=== Method 11 complete ===\n")
