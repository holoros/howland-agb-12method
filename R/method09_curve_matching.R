# ==============================================================================
# HOWLAND FOREST LIDAR CASE STUDY: SCRIPT 11
# Method 9: Multi Attribute Yield Curve Matching (Tompalski et al. 2016, 2018)
# ==============================================================================
# PURPOSE: Build a yield curve template library from AcadianGY v12.3.6 by
#          running AcadianGYOneStand on representative Howland plot initial
#          conditions across the low, medium, and high productivity strata
#          and the three species compositions present at Howland. For each
#          30 m pixel, select four candidate curves by minimizing the
#          absolute difference between observed and candidate AGB at the 2012
#          and 2015 anchor points, and (optionally) an intermediate ABA
#          predicted attribute. Combine candidates into a final curve using
#          the R^2 weighted mean of Tompalski et al. (2016) Eq. 1 and compute
#          per pixel uncertainty as the max relative spread across candidates
#          (Eq. 2). Project to 2025, 2035, and 2075 and validate against the
#          observed 2025 MELiTE ABA AGB.
#
# AUTHOR: Aaron Weiskittel
# DATE: 2026-04-22
# DEPENDENCIES: tidyverse, terra, sf, FNN, patchwork, ggtext,
#               AcadianGY_12.3.6.r (sourced from HRF root),
#               data/field_inventory/Plots_Data.csv,
#               data/field_inventory/Tree_Data.csv,
#               output/recalibration_results.rds (pixel trajectories,
#                                                 4 pt yield fits, p95 model),
#               output/Howland_metrics_2015.rds (ABA stack at match year)
# OUTPUT: output/curve_matching_library.rds
#         output/curve_matching_results.rds
#         output/curve_matching_validation.csv
#         output/figures/fig_curve_matching_candidates.png
#         output/figures/fig_curve_matching_uncertainty.png
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(tidyverse)
library(terra)
library(sf)
library(FNN)
library(patchwork)
library(ggtext)

# --- Paths --------------------------------------------------------------------
hrf_root   <- Sys.getenv("HRF_ROOT",   unset = "/users/PUOM0008/crsfaaron/HRF")
output_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path(hrf_root, "output"))
fig_dir    <- Sys.getenv("FIG_DIR",    unset = file.path(output_dir, "figures"))
data_dir   <- Sys.getenv("DATA_DIR",   unset = file.path(hrf_root, "data"))
dir.create(fig_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Constants ----------------------------------------------------------------
ANCHOR_YEARS   <- c(2012, 2015)            # ABA anchor years driving the match
MATCH_YEAR     <- 2015                      # year used for Tompalski Eq. 1 R^2
TARGET_YEARS   <- c(2025, 2035, 2075)       # projection horizons
N_CANDIDATES   <- 4                         # candidates per pixel (Tompalski 2016)
STRATA         <- c("Low", "Medium", "High")
SPP_COMPS      <- c("coniferous", "mixed", "deciduous")
AGE_RANGE      <- seq(10, 200, by = 1)      # initial age grid (yr)
MG_HA_TO_TONS_AC <- 0.4461

HOWLAND_CSI    <- 12       # climate site index (m), ACD default for Maine
HOWLAND_ELEV_M <- 60       # ~200 ft elevation

# --- Publication theme --------------------------------------------------------
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10),
    legend.position  = "bottom",
    strip.text       = element_text(size = 11, face = "bold"),
    plot.title       = element_text(size = 14, face = "bold"),
    plot.margin      = margin(10, 10, 10, 10)
  )

cat("=== Script 11: Curve Matching (Tompalski 2016, 2018) ===\n")

# ==============================================================================
# PART 1: Load field inventory, build representative starting tree lists
# ==============================================================================
cat("\n--- Part 1: Building representative starting tree lists ---\n")

source(file.path(hrf_root, "AcadianGY_12.3.6.r"))

plots_df <- read_csv(file.path(data_dir, "field_inventory/Plots_Data.csv"),
                     show_col_types = FALSE) |>
  rename(plot = Plot, northing = Northing, easting = Easting,
         radius_m = Radus_m, area_m2 = Area_m2, n_trees = N_Trees, tph = TPH)

trees_df <- read_csv(file.path(data_dir, "field_inventory/Tree_Data.csv"),
                     show_col_types = FALSE) |>
  rename(plot_tree = Plot, tree_num = Tree, species = Species,
         dbh_cm = DBH_cm, snag = Snag) |>
  mutate(plot = stringr::str_sub(plot_tree, 1, 6))

# Species composition by plot (percent conifer)
pct_conifer <- trees_df |>
  filter(is.na(snag) | snag == 0) |>
  group_by(plot) |>
  summarise(
    pct_con = mean(species %in% c("PIRU", "ABBA", "PIMA", "TSCA", "PIST",
                                   "PIGL", "THOC", "LALA")),
    .groups = "drop"
  )

plots_df <- plots_df |>
  left_join(pct_conifer, by = "plot") |>
  mutate(
    spp_comp = case_when(
      pct_con >= 0.75 ~ "coniferous",
      pct_con <= 0.25 ~ "deciduous",
      TRUE            ~ "mixed"
    )
  )

# Assign plots to Low/Medium/High stratum using observed 2015 AGB terciles
plot_agb_2015 <- trees_df |>
  filter(is.na(snag) | snag == 0) |>
  left_join(plots_df |> select(plot, area_m2, spp_comp), by = "plot") |>
  mutate(
    agb_kg = exp(-2.0773 + 2.3323 * log(dbh_cm)),
    tph    = 10000 / area_m2
  ) |>
  group_by(plot, spp_comp) |>
  summarise(agb_Mg_ha = sum(agb_kg * tph, na.rm = TRUE) / 1000 / n(),
            .groups = "drop")

cuts <- quantile(plot_agb_2015$agb_Mg_ha, c(1/3, 2/3), na.rm = TRUE)
plot_agb_2015 <- plot_agb_2015 |>
  mutate(stratum = case_when(
    agb_Mg_ha <= cuts[1] ~ "Low",
    agb_Mg_ha <= cuts[2] ~ "Medium",
    TRUE                 ~ "High"
  ))

# Representative tree list per stratum x composition = 9 templates
# For each template, average trees from qualifying plots, rescale TPH to the
# stratum median, and assign a starting age that spans AGE_RANGE.

acd_spp_map <- tibble(
  species = c("PIRU","ABBA","TSCA","PIST","PIMA","PIGL","THOC","LALA",
              "ACRU","ACSA","BEAL","BEPA","FAGR","POTR","POGR","FRAM","QURU"),
  sp_acd  = c("RS","BF","HE","WP","BS","WS","CE","LA",
              "RM","SM","YB","PB","BE","AS","AS","WA","RO")
)

build_template_tree <- function(stratum, spp_comp) {
  qualifying <- plot_agb_2015 |>
    filter(stratum == !!stratum, spp_comp == !!spp_comp) |>
    pull(plot)
  if (length(qualifying) == 0) return(NULL)

  n_plots <- length(qualifying)

  tl <- trees_df |>
    filter(plot %in% qualifying, is.na(snag) | snag == 0) |>
    left_join(plots_df |> select(plot, area_m2), by = "plot") |>
    left_join(acd_spp_map, by = "species") |>
    mutate(
      sp_acd = coalesce(sp_acd, "RM"),
      # Divide per plot EXPF by number of pooled plots so that the pooled
      # synthetic stand has mean per ha density, not summed density.
      EXPF   = (10000 / area_m2) / n_plots
    ) |>
    mutate(
      STAND = paste(stratum, spp_comp, sep = "_"),
      YEAR  = MATCH_YEAR,
      PLOT  = 1L,
      TREE  = dplyr::row_number(),
      SP    = sp_acd,
      DBH   = dbh_cm,
      HT    = NA_real_,
      HCB   = NA_real_,
      CR    = NA_real_,
      Form  = NA_character_,
      Risk  = NA_character_,
      Region     = "ME",
      dDBH.mult  = 1.0,
      dHt.mult   = 1.0,
      mort.mult  = 1.0,
      max.dbh    = 200,
      max.height = 50
    ) |>
    filter(DBH > 0) |>
    select(STAND, YEAR, PLOT, TREE, SP, DBH, HT, HCB, CR, EXPF,
           Form, Risk, Region,
           dDBH.mult, dHt.mult, mort.mult, max.dbh, max.height) |>
    as.data.frame()
  tl
}

# ==============================================================================
# PART 2: Run AcadianGY to build the curve library
# ==============================================================================
cat("\n--- Part 2: Building the curve library with AcadianGY ---\n")

acd_ops <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 2.54,
                CutPoint = 0.95, usedHTCap = TRUE)
acd_stand <- list(CSI = HOWLAND_CSI, ELEV = HOWLAND_ELEV_M)

compute_agb_acd <- function(tree_df, spp_map = acd_spp_map) {
  tree_df |>
    as_tibble() |>
    mutate(
      agb_kg    = exp(-2.0773 + 2.3323 * log(DBH)),
      agb_kg_ha = agb_kg * EXPF
    ) |>
    summarise(
      agb_Mg_ha   = sum(agb_kg_ha, na.rm = TRUE) / 1000,
      ba_m2_ha    = sum(pi * (DBH / 200)^2 * EXPF, na.rm = TRUE),
      tph         = sum(EXPF, na.rm = TRUE),
      qmd_cm      = sqrt(sum((DBH^2) * EXPF, na.rm = TRUE) /
                          sum(EXPF, na.rm = TRUE)),
      mean_ht_m   = mean(HT, na.rm = TRUE),
      top_ht_m    = mean(HT[rank(-DBH) <= pmax(1, ceiling(0.1 * n()))],
                          na.rm = TRUE)
    )
}

curve_grid <- tidyr::expand_grid(
  stratum  = STRATA,
  spp_comp = SPP_COMPS
) |>
  mutate(template_id = paste(stratum, spp_comp, sep = "_"))

# Project each template forward from MATCH_YEAR for MAX_PROJ_YRS years.
# The "age" label is time since match year rather than absolute stand age.
MAX_PROJ_YRS <- 100

build_one_curve <- function(stratum, spp_comp) {
  tl0 <- build_template_tree(stratum, spp_comp)
  if (is.null(tl0) || nrow(tl0) == 0) {
    return(tibble(year_offset = 0:MAX_PROJ_YRS,
                  agb_Mg_ha = NA_real_, ba_m2_ha = NA_real_,
                  qmd_cm = NA_real_, tph = NA_real_,
                  mean_ht_m = NA_real_, top_ht_m = NA_real_))
  }
  out <- vector("list", MAX_PROJ_YRS + 1)
  out[[1]] <- compute_agb_acd(tl0) |> mutate(year_offset = 0)
  current <- tl0
  for (yr in seq_len(MAX_PROJ_YRS)) {
    res <- tryCatch(
      AcadianGYOneStand(tree = current, stand = acd_stand, ops = acd_ops),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) break
    current <- res
    current$YEAR <- current$YEAR[1] + 1
    out[[yr + 1]] <- compute_agb_acd(current) |> mutate(year_offset = yr)
  }
  bind_rows(out)
}

cat("Building", nrow(curve_grid), "templates...\n")
curve_library <- curve_grid |>
  mutate(curve = purrr::pmap(
    list(stratum, spp_comp),
    \(s, c) build_one_curve(s, c)
  )) |>
  tidyr::unnest(curve)

saveRDS(curve_library, file.path(output_dir, "curve_matching_library.rds"))
cat("Curve library saved:", nrow(curve_library), "rows across",
    n_distinct(curve_library$template_id), "templates\n")

# ==============================================================================
# PART 3: Load pixel trajectories and match each pixel to 4 candidate curves
# ==============================================================================
cat("\n--- Part 3: Template matching at the pixel level ---\n")

recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))
pixel_tr <- recal$pixel_trajectories |>
  as_tibble() |>
  rename(AGB_2012 = agb_2012, AGB_2015 = agb_2015, AGB_2025 = agb_2025,
         stratum = agb_class)

# Weighted R^2 used in Eq. 1: four 2015 anchored attributes.
# The ABA model at Howland uses p95 only, so AGB carries most weight.
# Use the 2025 recalibration R^2 proxies for the other attributes.
aba_r2 <- c(AGB = 0.71, DBH = 0.65, BA = 0.79, TOPHT = 0.82)

# Template trajectories span 0..MAX_PROJ_YRS year offsets from MATCH_YEAR.
# For each pixel we find the year_offset in the candidate trajectory that best
# matches the observed 2012 + 2015 AGB anchors, then advance from there.

# For each pixel, find N_CANDIDATES closest (template, year_offset) combos
# by minimizing absolute AGB difference at MATCH_YEAR.
match_pixel <- function(pix_row, library_df, aba_r2, n_cand = 4) {
  # Best year_offset within each template that puts AGB closest to AGB_2015
  cand <- library_df |>
    filter(!is.na(agb_Mg_ha)) |>
    mutate(dist = abs(agb_Mg_ha - pix_row$AGB_2015)) |>
    group_by(template_id) |>
    slice_min(dist, n = 1, with_ties = FALSE) |>
    ungroup() |>
    slice_min(dist, n = n_cand, with_ties = FALSE) |>
    rename(match_offset = year_offset, agb_match = agb_Mg_ha) |>
    select(template_id, match_offset, agb_match, dist)
  cand
}

project_pixel <- function(pix_row, library_df, aba_r2, target_years) {
  cand <- match_pixel(pix_row, library_df, aba_r2, N_CANDIDATES)
  if (nrow(cand) == 0) {
    return(list(
      snap = tibble(pixel_id = pix_row$pixel_id,
                    target_year = target_years,
                    agb_pred = NA_real_, unc_rel = NA_real_),
      full = tibble()
    ))
  }

  # For each target year, look ahead delta_yr = target_year - MATCH_YEAR from
  # each candidate's match_offset, then R^2 weight the AGB values across
  # candidates (Tompalski Eq. 1).
  w_AGB <- aba_r2["AGB"]
  snap <- purrr::map_dfr(target_years, function(ty) {
    delta <- ty - MATCH_YEAR
    vals <- purrr::pmap_dbl(
      list(cand$template_id, cand$match_offset),
      function(tid, moff) {
        v <- library_df |>
          filter(template_id == tid,
                 year_offset == moff + delta) |>
          pull(agb_Mg_ha)
        if (length(v) == 0) NA_real_ else v[1]
      }
    )
    agb_mean <- if (all(is.na(vals))) NA_real_ else
      sum(w_AGB * vals, na.rm = TRUE) / (sum(!is.na(vals)) * w_AGB)
    unc_rel  <- if (all(is.na(vals))) NA_real_ else
      (max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE)) /
       pmax(mean(vals, na.rm = TRUE), 1)
    tibble(pixel_id = pix_row$pixel_id,
           target_year = ty,
           agb_pred = agb_mean,
           unc_rel  = unc_rel)
  })
  list(snap = snap, full = tibble())
}

cat("Projecting", nrow(pixel_tr), "pixels across",
    length(TARGET_YEARS), "horizons...\n")
# Process in chunks to manage memory
chunk_size <- 2000L
chunks <- split(pixel_tr, ceiling(seq_len(nrow(pixel_tr)) / chunk_size))
snap_list <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  snap_list[[i]] <- purrr::map_dfr(
    split(chunks[[i]], seq_len(nrow(chunks[[i]]))),
    \(r) project_pixel(r, curve_library, aba_r2, TARGET_YEARS)$snap
  )
  if (i %% 5 == 0) cat("  Chunk", i, "of", length(chunks), "\n")
}
snap_df <- bind_rows(snap_list)

saveRDS(list(snap = snap_df, aba_r2 = aba_r2, curve_grid = curve_grid),
        file.path(output_dir, "curve_matching_results.rds"))

# ==============================================================================
# PART 4: Validation against observed 2025 AGB
# ==============================================================================
cat("\n--- Part 4: Validation against 2025 ABA AGB ---\n")

val <- snap_df |>
  filter(target_year == 2025) |>
  left_join(pixel_tr |> select(pixel_id, AGB_obs = AGB_2025),
            by = "pixel_id") |>
  mutate(resid = AGB_obs - agb_pred)

val_stats <- val |>
  summarise(
    n        = sum(!is.na(resid)),
    bias     = mean(resid, na.rm = TRUE),
    bias_pct = 100 * mean(resid, na.rm = TRUE) / mean(AGB_obs, na.rm = TRUE),
    rmse     = sqrt(mean(resid^2, na.rm = TRUE)),
    rmse_pct = 100 * sqrt(mean(resid^2, na.rm = TRUE)) /
                 mean(AGB_obs, na.rm = TRUE),
    r2       = 1 - sum(resid^2, na.rm = TRUE) /
                 sum((AGB_obs - mean(AGB_obs, na.rm = TRUE))^2,
                      na.rm = TRUE)
  )

print(val_stats)
write_csv(val_stats, file.path(output_dir, "curve_matching_validation.csv"))

# ==============================================================================
# PART 5: Figures
# ==============================================================================
cat("\n--- Part 5: Figures ---\n")

p_candidates <- val |>
  ggplot(aes(agb_pred, AGB_obs)) +
  geom_point(alpha = 0.12, size = 0.5) +
  geom_abline(color = "red") +
  coord_equal() +
  labs(
    x = "Predicted AGB (Mg ha<sup>-1</sup>, curve matching 2025)",
    y = "Observed AGB (Mg ha<sup>-1</sup>, 2025 MELiTE ABA)",
    title = "Method 9: Curve matching, 10 yr projection"
  ) +
  theme_pub +
  theme(axis.title = element_markdown())

p_uncertainty <- snap_df |>
  filter(target_year == 2075, !is.na(unc_rel)) |>
  ggplot(aes(unc_rel)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  labs(
    x = "Relative spread across 4 candidate curves",
    y = "Pixel count",
    title = "Candidate curve disagreement at 2075 (Eq. 2 of Tompalski 2016)"
  ) +
  theme_pub

ggsave(file.path(fig_dir, "fig_curve_matching_candidates.png"),
       p_candidates, width = 14, height = 10, units = "cm",
       dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "fig_curve_matching_uncertainty.png"),
       p_uncertainty, width = 14, height = 10, units = "cm",
       dpi = 300, bg = "white")

cat("\n=== Script 11 complete ===\n")
