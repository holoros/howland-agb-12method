# ==============================================================================
# HOWLAND FOREST: METHOD 10
# Wall-to-wall pure Kohek (2022) directional projection with 2025 ALS
# decimated to 25 pts/m²
# ==============================================================================
# Decimates the 26 normalized 2025 ALS tiles from ~92 pts/m² to 25 pts/m²
# using random thinning, writes decimated copies to a scratch folder, then
# runs the same ITC segmentation + directional expansion pipeline as
# Method 11 but anchored on freshly segmented 2025 trees (instead of the
# 2017 ITC). Parallelized across tiles via future.apply.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lidR)
  library(sf)
  library(terra)
  library(FNN)
  library(future)
  library(future.apply)
})

hrf_root   <- Sys.getenv("HRF_ROOT",   unset = "/users/PUOM0008/crsfaaron/HRF")
output_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path(hrf_root, "output"))
field_dir  <- Sys.getenv("FIELD_DIR",  unset = file.path(hrf_root, "data/field_inventory"))
norm_dir   <- Sys.getenv("NORM_DIR",
                          unset = file.path(output_dir, "als_2025_processing"))
dec_dir    <- file.path(output_dir, "als_2025_decimated_25")
dir.create(dec_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_DENSITY <- 25
M0             <- 0.2
M1             <- 0.9
ITC_HMIN_M     <- 3.0
HORIZONS       <- c(2035, 2045, 2075)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
future::plan(multisession, workers = min(n_cores, 8))
lidR::set_lidr_threads(2)  # avoid oversubscription inside each worker

cat("=== Method 10: Wall-to-wall Kohek, ALS decimated to 25 pts/m² ===\n")
cat("future workers:", min(n_cores, 8),
    " lidR threads per worker:", lidR::get_lidr_threads(), "\n")

# --- Step 1: Decimate tiles if not already done ------------------------------
las_files <- list.files(norm_dir, pattern = "_norm\\.la[sz]$",
                         full.names = TRUE)
cat("Source tiles:", length(las_files), "\n")

decimate_one <- function(f) {
  out <- file.path(dec_dir, gsub("_norm\\.", "_dec25.", basename(f)))
  if (file.exists(out)) return(out)
  las <- readLAS(f, filter = paste0("-drop_z_below ", ITC_HMIN_M))
  if (is.empty(las)) return(NA)
  las_d <- decimate_points(las, random(TARGET_DENSITY))
  writeLAS(las_d, out)
  out
}

dec_start <- Sys.time()
dec_files <- future.apply::future_lapply(las_files, decimate_one,
  future.seed = TRUE, future.packages = c("lidR"))
dec_files <- unlist(dec_files)
dec_files <- dec_files[!is.na(dec_files)]
cat("Decimation elapsed:",
    round(as.numeric(difftime(Sys.time(), dec_start, units = "secs"))), "s\n")
cat("Decimated tiles:", length(dec_files), "\n")

# --- Step 2: Segment trees from each decimated tile --------------------------
alg_li <- li2012(dt1 = 1.5, dt2 = 2.0, R = 2.0,
                  Zu = 15, speed_up = 10, hmin = ITC_HMIN_M)

segment_one_tile <- function(f) {
  las <- readLAS(f, filter = paste0("-drop_z_below ", ITC_HMIN_M))
  if (is.empty(las)) return(NULL)
  las <- segment_trees(las, algorithm = alg_li)
  tm <- crown_metrics(las, func = .stdtreemetrics, geom = "point")
  if (is.null(tm) || nrow(tm) == 0) return(NULL)
  xy <- st_coordinates(tm)[, 1:2]
  data.frame(
    treeID         = tm$treeID,
    Z              = tm$Z,
    crown_radius_m = sqrt(pmax(tm$convhull_area, 1) / pi),
    X              = xy[,1],
    Y              = xy[,2]
  )
}

seg_start <- Sys.time()
trees_list <- future.apply::future_lapply(dec_files, segment_one_tile,
  future.seed = TRUE, future.packages = c("lidR", "sf"))
trees_list <- Filter(Negate(is.null), trees_list)
seg <- do.call(rbind, trees_list)
cat("Segmentation elapsed:",
    round(as.numeric(difftime(Sys.time(), seg_start, units = "secs"))), "s\n")
cat("Total trees segmented:", nrow(seg), "\n")
saveRDS(seg, file.path(output_dir, "method10_segmented_trees_2025.rds"))

future::plan(sequential)  # release workers for memory
cat("Cleared future workers\n")

# --- Step 3: Attribute trees to reference 10 m grid -------------------------
ref_raster <- rast(file.path(output_dir, "agb_2025_melite_30m.tif"))
pixel_area_ha <- prod(res(ref_raster)) / 10000

seg$pixel_id <- cellFromXY(ref_raster, as.matrix(seg[, c("X","Y")]))
seg <- seg |> filter(!is.na(pixel_id), Z >= ITC_HMIN_M, Z <= 45)

recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))
pt <- recal$pixel_trajectories |> as_tibble()
seg <- seg |> inner_join(pt |> select(pixel_id, agb_2015,
                                        agb_2025, agb_class),
                          by = "pixel_id")
cat("Trees inside Howland pixels:", nrow(seg), "\n")

# --- Step 4: Species assignment, DBH imputation, projection ------------------
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
nnp <- get.knnx(plot_coords, as.matrix(seg[, c("X","Y")]), k = 1)
seg$spp_comp <- plots |> filter(!is.na(easting)) |>
  pull(spp_comp) |>
  (\(v) v[nnp$nn.index[,1]])()
seg$stratum3    <- seg$agb_class
seg$template_id <- paste(seg$stratum3, seg$spp_comp, sep = "_")

# DBH from Z via default power H-D
hd_a <- 0.35; hd_b <- 1.30
seg$dbh_2025 <- pmax(hd_a * seg$Z^hd_b, 1)
seg$agb_kg_2025 <- exp(-2.0773 + 2.3323 * log(seg$dbh_2025))

# Project forward using curve library top_ht growth
cm <- readRDS(file.path(output_dir, "curve_matching_library.rds")) |>
  as_tibble()
anchor_offset  <- 10   # 2025 - 2015
horizon_offset <- HORIZONS - 2015
tpl <- cm |>
  filter(year_offset %in% c(anchor_offset, horizon_offset)) |>
  group_by(template_id, year_offset) |>
  summarise(top_ht = mean(top_ht_m, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(id_cols = template_id, names_from = year_offset,
              values_from = top_ht, names_prefix = "topht_")
seg <- seg |> left_join(tpl, by = "template_id")
anchor_col <- paste0("topht_", anchor_offset)
for (h in HORIZONS) {
  offs <- h - 2015
  scale_col <- paste0("topht_", offs)
  seg[[paste0("ht_", h)]] <- pmin(seg$Z * seg[[scale_col]] / seg[[anchor_col]],
                                    45)
  seg[[paste0("dbh_", h)]] <- pmax(hd_a * seg[[paste0("ht_", h)]]^hd_b, 1)
  seg[[paste0("agb_kg_", h)]] <- exp(-2.0773 + 2.3323 *
    log(seg[[paste0("dbh_", h)]]))
}

# --- Step 5: Kohek expansion proxy (shade via neighbor max) ------------------
r_ht <- rast(ext(ref_raster), resolution = 2, crs = crs(ref_raster))
pts <- vect(as.matrix(seg[, c("X","Y")]), type = "points",
             atts = data.frame(z = seg$Z), crs = crs(ref_raster))
r_ht <- rasterize(pts, r_ht, field = "z", fun = "max")
r_fmax <- focal(r_ht, w = matrix(1,3,3), fun = "max", na.rm = TRUE)
seg$neigh_max <- terra::extract(r_fmax, as.matrix(seg[, c("X","Y")]))[,1]
seg$neigh_max[is.na(seg$neigh_max)] <- seg$Z[is.na(seg$neigh_max)]
seg$shade <- pmin(pmax((seg$neigh_max - seg$Z) / pmax(seg$Z, 5), 0), 1)

for (h in HORIZONS) {
  yrs <- h - 2025
  g <- 1 + (M0 + M1 * seg$shade) / 50
  seg[[paste0("r_", h)]] <- seg$crown_radius_m * pmin(g^yrs, 2.0)
}

# --- Step 6: Pixel sum + calibration -----------------------------------------
raw <- seg |>
  group_by(pixel_id) |>
  summarise(agb_raw_2025 = sum(agb_kg_2025) / 1000 / pixel_area_ha,
            agb_raw_2035 = sum(agb_kg_2035) / 1000 / pixel_area_ha,
            agb_raw_2045 = sum(agb_kg_2045) / 1000 / pixel_area_ha,
            agb_raw_2075 = sum(agb_kg_2075) / 1000 / pixel_area_ha,
            n_trees      = n(), .groups = "drop")

cal_df <- raw |>
  inner_join(pt |> select(pixel_id, agb_2025), by = "pixel_id") |>
  drop_na(agb_raw_2025, agb_2025) |>
  filter(agb_raw_2025 > 0, agb_2025 > 0, n_trees >= 2)
cal_fit <- lm(agb_2025 ~ agb_raw_2025, data = cal_df)
alpha <- coef(cal_fit)[1]; beta <- coef(cal_fit)[2]
if (!is.finite(beta) || beta <= 0) {
  beta <- sum(cal_df$agb_2025) / sum(cal_df$agb_raw_2025); alpha <- 0
}
cat("Calibration: alpha =", round(alpha,1), ", beta =", round(beta,3), "\n")

pred <- raw
for (h in c(2025, HORIZONS)) {
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
cat("\n--- Method 10 validation (2025 self fit) ---\n")
print(val_stats)
write_csv(val_stats, file.path(output_dir, "method10_validation.csv"))
saveRDS(list(pred = pred, val = val_stats, val_pixels = val,
              alpha = alpha, beta = beta),
        file.path(output_dir, "method10_results.rds"))

cat("\n=== Method 10 complete ===\n")
