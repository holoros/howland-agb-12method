# =============================================================================
# Title: FVS-NE and FVS-ACD Growth Projection with Multi-Method Comparison
# Author: A. Weiskittel
# Date: 2026-03-22
# Description: Complete rewrite using the ACTUAL AcadianGY_12.3.6.r code for
#              FVS-ACD projections, plus FVS-NE in two forms (out-of-box and
#              locally calibrated). Builds stand-level yield curves from
#              corrected projections and applies them to LiDAR/NAIP pixels.
#              Includes enhanced Monte Carlo for short (10yr) and long (20+yr)
#              cost/uncertainty analysis.
#
#              Methods compared:
#                (1) Area-based Chapman Richards yield curves (Script 01/04)
#                (2) Area-based Schumacher yield model (Script 01/04)
#                (3) FVS-NE individual tree projection (out-of-box defaults)
#                (4) FVS-NE individual tree projection (locally calibrated)
#                (5) FVS-ACD individual tree projection (full AcadianGY v12.3.5)
#                (6) NAIP + FVS-ACD yield curve (Aaron's preferred approach)
#                (7) LiDAR + FVS-ACD yield curve
#                (8) Observed LiDAR AGB (ORNL DAAC + Phoenix RANGER)
#
# References:
#   Teck, R.M. & Hilt, D.E. (1991). Individual-tree diameter growth model
#     for the Northeastern United States. USFS Res. Paper NE-649.
#   FVS Staff (2008, rev. 2025). Northeast (NE) Variant Overview.
#   Kuehne, C., Russell, M., Weiskittel, A. & Kershaw Jr, J. (2020).
#     Comparing strategies for representing individual-tree secondary growth
#     in mixed-species stands. For. Ecol. Manage. 459:117823.
#   Chen, C., Rijal, B. & Weiskittel, A. (in review). Comparative assessment
#     of time-explicit, state-space and simultaneous models for stand-level
#     volume growth and yield predictions.
#   Chojnacky, D.C., Heath, L.S. & Jenkins, J.C. (2014). Updated generalized
#     biomass equations for North American tree species. Forestry 87(1):129-151.
#
# Dependencies: Scripts 01-04 outputs, field inventory data, AcadianGY_12.3.6.r
# =============================================================================

# --- Libraries ---------------------------------------------------------------
# plyr dependency removed: AcadianGY_12.3.6.r modernized to use dplyr (2026-03-24)
library(tidyverse)
library(terra)
library(patchwork)
library(ggtext)
library(readxl)

# --- Paths -------------------------------------------------------------------
data_dir    <- Sys.getenv("FIELD_DIR",
                           unset = "/home/aweiskittel/Documents/MAINE/DATA/HowlandForest")
raster_dir  <- Sys.getenv("ORNL_DIR",
                           unset = file.path(data_dir, "agb_rasters"))
output_dir  <- Sys.getenv("OUTPUT_DIR", unset = "output")
fig_dir     <- Sys.getenv("FIG_DIR", unset = "output/howland_validation")
project_dir <- Sys.getenv("PROJECT_DIR", unset = dirname(data_dir))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Source AcadianGY_12.3.6.r -----------------------------------------------
# Look for the ACD source in several locations
acd_paths <- c(
  file.path(data_dir, "AcadianGY_12.3.6.r"),
  file.path(data_dir, "scripts", "AcadianGY_12.3.6.r"),
  file.path(dirname(data_dir), "AcadianGY_12.3.6.r"),
  file.path(project_dir, "AcadianGY_12.3.6.r"),
  # Cardinal HPC: script and ACD source sit together in scripts/
  file.path(getwd(), "AcadianGY_12.3.6.r"),
  file.path(getwd(), "scripts", "AcadianGY_12.3.6.r"),
  "AcadianGY_12.3.6.r"
)
acd_source <- acd_paths[file.exists(acd_paths)][1]
if (is.na(acd_source)) {
  stop("Cannot find AcadianGY_12.3.6.r. Place it in: ", data_dir)
}
cat("Sourcing AcadianGY from:", acd_source, "\n")
source(acd_source)
cat("AcadianGY version:", AcadianVersionTag, "\n")

# --- Publication theme -------------------------------------------------------
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Unit conversion constants
MG_HA_TO_TONS_AC <- 0.4461
CM_TO_IN <- 1 / 2.54
M_TO_FT <- 3.28084
M2_HA_TO_FT2_AC <- 4.35596

# =============================================================================
# PART 1: Species Mapping and Data Preparation
# =============================================================================

cat("\n========== PART 1: Data Preparation ==========\n")

# 4-letter ITIS codes to ACD 2-letter codes
spp_map <- tribble(
  ~species_4, ~sp_acd,
  "ABBA",     "BF",    # balsam fir
  "PIRU",     "RS",    # red spruce
  "TSCA",     "EH",    # eastern hemlock
  "PIST",     "WP",    # white pine
  "THOC",     "WC",    # northern white cedar
  "ACRU",     "RM",    # red maple
  "BEPA",     "PB",    # paper birch
  "FAGR",     "AB",    # American beech
  "QURU",     "RO",    # red oak
  "POTR",     "QA",    # quaking aspen
  "POGR",     "BT",    # bigtooth aspen
  "BEPO",     "GB",    # gray birch
  "BEAL",     "YB",    # yellow birch
  "FRNI",     "BA",    # black ash
  "ACSA",     "SM",    # sugar maple
  "PIAB",     "NS",    # Norway spruce
  "OSVI",     "HH",    # eastern hophornbeam
  "FRPE",     "GA",    # green ash
  "POBA",     "BP",    # balsam poplar
  "PIMA",     "BS",    # black spruce
  "ACPL",     "OH",    # Norway maple -> other HW
  "UNK",      "99"     # unknown
)

# Chojnacky allometric coefficients (4-letter codes for compatibility)
choj_coefs <- tribble(
  ~species_4, ~choj_b0, ~choj_b1,
  "PIRU",     -2.0773,  2.3323,
  "PIAB",     -2.0773,  2.3323,
  "PIMA",     -2.0773,  2.3323,
  "ABBA",     -2.5356,  2.4349,
  "TSCA",     -2.5356,  2.4349,
  "THOC",     -2.0773,  2.3323,
  "PIST",     -2.6177,  2.4638,
  "ACRU",     -1.9123,  2.3651,
  "ACSA",     -1.9123,  2.3651,
  "BEAL",     -2.5356,  2.4349,
  "BEPA",     -2.5356,  2.4349,
  "BEPO",     -2.5356,  2.4349,
  "FAGR",     -1.9123,  2.3651,
  "POGR",     -2.4441,  2.4561,
  "POTR",     -2.4441,  2.4561,
  "QURU",     -2.0773,  2.3323,
  "FRNI",     -1.9123,  2.3651,
  "FRPE",     -1.9123,  2.3651,
  "OSVI",     -1.9123,  2.3651,
  "POBA",     -2.4441,  2.4561,
  "ACPL",     -1.9123,  2.3651,
  "UNK",      -2.0773,  2.3323
)

# Load field inventory
tree_data <- read_csv(file.path(data_dir, "Tree_Data.csv"),
                       show_col_types = FALSE) |>
  rename(species_4 = Species, dbh_cm = DBH_cm) |>
  mutate(
    plot = str_extract(Plot, "^[A-Z]+_\\d+"),
    tree_num = as.integer(str_extract(Plot, "\\d+$")),
    is_snag = !is.na(Snag) & Snag == "X"
  )

plots_data <- read_csv(file.path(data_dir, "Plots_Data.csv"),
                        show_col_types = FALSE) |>
  rename(plot = Plot, area_m2 = Area_m2)

# Filter to live Howland trees with valid DBH
hl_trees <- tree_data |>
  filter(str_starts(plot, "HL"), !is_snag, !is.na(dbh_cm), dbh_cm > 0) |>
  left_join(plots_data |> select(plot, area_m2), by = "plot") |>
  filter(!is.na(area_m2)) |>
  left_join(spp_map, by = "species_4") |>
  mutate(sp_acd = coalesce(sp_acd, "99"))

cat("Howland live tree list:", n_distinct(hl_trees$plot), "plots,",
    nrow(hl_trees), "trees\n")
cat("Species composition:\n")
hl_trees |> count(species_4, sp_acd, sort = TRUE) |> print(n = 25)

# Howland site parameters
HOWLAND_ELEV_M  <- 60    # ~200 ft elevation
HOWLAND_CSI     <- 12    # climate site index (m), ACD default for Maine
HOWLAND_YEAR    <- 2015  # field inventory date

# =============================================================================
# PART 2: FVS-ACD Projection via AcadianGY_12.3.6.r (Full Model)
# =============================================================================

cat("\n========== PART 2: FVS-ACD Projection (Full AcadianGY v12.3.5) ==========\n")

# Format tree list for AcadianGYOneStand
# Required columns: STAND, YEAR, PLOT, TREE, SP, DBH, HT, EXPF
# HT can be NA (will be predicted from CSI-based H-D model)
# EXPF = expansion factor (trees per hectare per stem)

format_for_acd <- function(trees_df, year = 2015) {
  trees_df |>
    mutate(
      STAND  = "HOWLAND",
      YEAR   = year,
      PLOT   = plot,
      TREE   = tree_num,
      SP     = sp_acd,
      DBH    = dbh_cm,
      HT     = NA_real_,        # let ACD predict from CSI
      HCB    = NA_real_,        # let ACD predict from model
      CR     = NA_real_,        # will be computed internally
      EXPF   = 10000 / area_m2, # trees/ha per stem
      Form   = NA_character_,
      Risk   = NA_character_,
      Region = "ME",
      # Calibration multipliers (1.0 = no adjustment from defaults)
      dDBH.mult  = 1.0,
      dHt.mult   = 1.0,
      mort.mult  = 1.0,
      max.dbh    = 200,         # max DBH cap in cm (generous default)
      max.height = 50           # max height cap in m (generous default)
    ) |>
    filter(DBH > 0) |>
    select(STAND, YEAR, PLOT, TREE, SP, DBH, HT, HCB, CR, EXPF,
           Form, Risk, Region,
           dDBH.mult, dHt.mult, mort.mult, max.dbh, max.height) |>
    as.data.frame()
}

acd_tree_input <- format_for_acd(hl_trees, year = HOWLAND_YEAR)
cat("ACD input: ", nrow(acd_tree_input), "trees across",
    n_distinct(acd_tree_input$PLOT), "plots\n")

# ACD projection options
acd_ops <- list(
  verbose   = FALSE,
  INGROWTH  = "Y",
  MinDBH    = 2.54,  # minimum DBH in cm (1 inch)
  CutPoint  = 0.95,
  usedHTCap = TRUE
)

acd_stand <- list(
  CSI  = HOWLAND_CSI,
  ELEV = HOWLAND_ELEV_M
)

# Function to compute AGB from a tree dataframe (ACD format)
compute_agb_acd <- function(tree_df, choj = choj_coefs, spp_mapping = spp_map) {
  tree_df |>
    as_tibble() |>
    left_join(spp_mapping, by = c("SP" = "sp_acd")) |>
    left_join(choj, by = "species_4") |>
    mutate(
      choj_b0 = coalesce(choj_b0, -2.0773),
      choj_b1 = coalesce(choj_b1, 2.3323),
      agb_kg = exp(choj_b0 + choj_b1 * log(DBH)),
      agb_kg_ha = agb_kg * EXPF
    ) |>
    group_by(PLOT) |>
    summarise(
      n_trees = n(),
      ba_m2_ha = sum(pi * (DBH/200)^2 * EXPF, na.rm = TRUE),
      agb_Mg_ha = sum(agb_kg_ha, na.rm = TRUE) / 1000,
      mean_dbh_cm = mean(DBH, na.rm = TRUE),
      .groups = "drop"
    )
}

# Run FVS-ACD for 50 years (1 year at a time)
n_years_long <- 50
acd_results <- list()

# Year 0 (initial)
acd_results[[1]] <- compute_agb_acd(acd_tree_input) |>
  mutate(year_offset = 0, year = HOWLAND_YEAR)

current_trees_acd <- acd_tree_input
cat("Projecting FVS-ACD for", n_years_long, "years...\n")

for (yr in seq_len(n_years_long)) {
  if (yr %% 5 == 0) cat("  Year", yr, "of", n_years_long, "\n")

  # AcadianGYOneStand grows trees for one cycle (1 year)
  tryCatch({
    current_trees_acd <- AcadianGYOneStand(
      tree  = current_trees_acd,
      stand = acd_stand,
      ops   = acd_ops
    )
    # Update YEAR for next iteration
    current_trees_acd$YEAR <- HOWLAND_YEAR + yr
  }, error = function(e) {
    cat("  ACD error at year", yr, ":", conditionMessage(e), "\n")
    cat("  Continuing with last successful tree list\n")
  })

  acd_results[[yr + 1]] <- compute_agb_acd(current_trees_acd) |>
    mutate(year_offset = yr, year = HOWLAND_YEAR + yr)
}

acd_trajectory <- bind_rows(acd_results)

# Summarize across plots
acd_means <- acd_trajectory |>
  group_by(year_offset, year) |>
  summarise(
    mean_agb_Mg_ha = mean(agb_Mg_ha, na.rm = TRUE),
    sd_agb_Mg_ha   = sd(agb_Mg_ha, na.rm = TRUE),
    mean_ba_m2_ha  = mean(ba_m2_ha, na.rm = TRUE),
    n_plots        = n_distinct(PLOT),
    .groups = "drop"
  ) |>
  mutate(
    mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
    method = "FVS-ACD (full)"
  )

cat("\nFVS-ACD (full AcadianGY) projection summary:\n")
cat("  Year 0 AGB:", round(acd_means$mean_agb_Mg_ha[1], 1), "Mg/ha\n")
cat("  Year 10 AGB:", round(acd_means$mean_agb_Mg_ha[acd_means$year_offset == 10], 1), "Mg/ha\n")
cat("  Year 20 AGB:", round(acd_means$mean_agb_Mg_ha[acd_means$year_offset == 20], 1), "Mg/ha\n")
cat("  Year 30 AGB:", round(acd_means$mean_agb_Mg_ha[acd_means$year_offset == 30], 1), "Mg/ha\n")
cat("  Year 50 AGB:", round(acd_means$mean_agb_Mg_ha[acd_means$year_offset == 50], 1), "Mg/ha\n")
cat("  10yr change:", round(acd_means$mean_agb_Mg_ha[acd_means$year_offset == 10] -
    acd_means$mean_agb_Mg_ha[1], 1), "Mg/ha\n")

# --- ACD diagnostics: per plot trajectory for debugging ---
cat("\n--- ACD Diagnostics: per plot AGB at years 0, 5, 10 ---\n")
acd_diag <- acd_trajectory |>
  filter(year_offset %in% c(0, 5, 10)) |>
  select(PLOT, year_offset, n_trees, ba_m2_ha, agb_Mg_ha, mean_dbh_cm) |>
  pivot_wider(
    names_from = year_offset,
    values_from = c(n_trees, ba_m2_ha, agb_Mg_ha, mean_dbh_cm),
    names_sep = "_yr"
  )
print(acd_diag, n = 60)

# Which plots lost the most AGB?
acd_change <- acd_trajectory |>
  filter(year_offset %in% c(0, 10)) |>
  select(PLOT, year_offset, agb_Mg_ha, n_trees, ba_m2_ha) |>
  pivot_wider(names_from = year_offset, values_from = c(agb_Mg_ha, n_trees, ba_m2_ha),
              names_sep = "_yr") |>
  mutate(agb_change = agb_Mg_ha_yr10 - agb_Mg_ha_yr0,
         trees_change = n_trees_yr10 - n_trees_yr0,
         ba_change = ba_m2_ha_yr10 - ba_m2_ha_yr0) |>
  arrange(agb_change)
cat("\n--- ACD: Plots ranked by 10yr AGB change ---\n")
print(acd_change, n = 60)

# Tree count trajectory (stand level)
cat("\n--- ACD: Mean trees/plot and BA over time ---\n")
acd_trajectory |>
  filter(year_offset %in% seq(0, 50, by = 5)) |>
  group_by(year_offset, year) |>
  summarise(
    mean_n_trees = mean(n_trees, na.rm = TRUE),
    mean_ba = round(mean(ba_m2_ha, na.rm = TRUE), 1),
    mean_agb = round(mean(agb_Mg_ha, na.rm = TRUE), 1),
    mean_dbh = round(mean(mean_dbh_cm, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print(n = 20)

# =============================================================================
# PART 3: FVS-NE Projection (Two Versions)
# =============================================================================

cat("\n========== PART 3: FVS-NE Projection (Out-of-Box + Calibrated) ==========\n")

# --- 3A: Coefficients from FVS-NE Overview (Sept 2025) ---
# Large tree diameter growth: {4.7.1.1-3}
fvsne_di_coefs <- tribble(
  ~species_4, ~spp_group, ~B1,        ~B2,       ~B3,        ~default_si_ft,
  "ABBA",      1,         0.0008829,  0.0602785, 0.012785,   52,
  "PIRU",      4,         0.0008236,  0.0549439, 0.011942,   50,
  "PIAB",      4,         0.0008236,  0.0549439, 0.011942,   50,
  "PIMA",      4,         0.0008236,  0.0549439, 0.011942,   50,
  "PIST",      9,         0.0009050,  0.0517297, 0.012329,   65,
  "THOC",      9,         0.0009050,  0.0517297, 0.012329,   45,
  "TSCA",     10,         0.0008737,  0.0940538, 0.009149,   52,
  "ACRU",     12,         0.0007906,  0.0651982, 0.016191,   60,
  "ACSA",     13,         0.0007439,  0.0706905, 0.016240,   60,
  "BEAL",     14,         0.0006668,  0.0768212, 0.019046,   60,
  "BEPA",     15,         0.0009766,  0.0832328, 0.023978,   60,
  "BEPO",     15,         0.0009766,  0.0832328, 0.023978,   55,
  "FAGR",     17,         0.0006911,  0.0730441, 0.013029,   60,
  "POGR",     20,         0.0011885,  0.0920050, 0.016877,   60,
  "POTR",     19,         0.0008815,  0.1419212, 0.019904,   60,
  "QURU",     25,         0.0008920,  0.0979702, 0.018024,   65,
  "FRNI",     18,         0.0008992,  0.0925395, 0.015004,   55,
  "FRPE",     18,         0.0008992,  0.0925395, 0.015004,   60,
  "OSVI",     27,         0.0009567,  0.1038458, 0.020653,   50,
  "POBA",     20,         0.0011885,  0.0920050, 0.016877,   55
)

# Background mortality coefficients (Table 5.0.1)
fvsne_mort_coefs <- tribble(
  ~species_4, ~p0,       ~p1,
  "ABBA",     5.1676998, -0.0077681,
  "PIRU",     5.1676998, -0.0077681,
  "PIAB",     5.1676998, -0.0077681,
  "PIMA",     5.1676998, -0.0077681,
  "PIST",     5.5876999, -0.0053480,
  "THOC",     5.1676998, -0.0077681,
  "TSCA",     5.1676998, -0.0077681,
  "ACRU",     5.1676998, -0.0077681,
  "ACSA",     5.1676998, -0.0077681,
  "BEAL",     5.9617000, -0.0340128,
  "BEPA",     5.9617000, -0.0340128,
  "BEPO",     5.9617000, -0.0340128,
  "FAGR",     5.1676998, -0.0077681,
  "POGR",     5.9617000, -0.0340128,
  "POTR",     5.9617000, -0.0340128,
  "QURU",     5.9617000, -0.0340128,
  "FRNI",     5.1676998, -0.0077681,
  "FRPE",     5.1676998, -0.0077681,
  "OSVI",     5.1676998, -0.0077681,
  "POBA",     5.9617000, -0.0340128
)

# Height-diameter parameters (Chapman-Richards, NE calibration)
hd_params <- tribble(
  ~species_4, ~hd_a, ~hd_b, ~hd_c,
  "PIRU",     22.0,  0.030, 1.20,
  "TSCA",     24.0,  0.025, 1.10,
  "ABBA",     18.0,  0.035, 1.30,
  "PIST",     30.0,  0.020, 1.00,
  "THOC",     15.0,  0.025, 1.20,
  "ACRU",     22.0,  0.030, 1.10,
  "ACSA",     24.0,  0.025, 1.00,
  "BEAL",     22.0,  0.025, 1.10,
  "BEPA",     20.0,  0.030, 1.10,
  "BEPO",     20.0,  0.030, 1.10,
  "FAGR",     22.0,  0.025, 1.10,
  "POGR",     22.0,  0.030, 1.10,
  "POTR",     20.0,  0.035, 1.20,
  "QURU",     24.0,  0.020, 1.00,
  "PIAB",     25.0,  0.030, 1.20,
  "PIMA",     18.0,  0.035, 1.20,
  "POBA",     22.0,  0.030, 1.10,
  "FRNI",     20.0,  0.025, 1.10,
  "FRPE",     22.0,  0.025, 1.10,
  "OSVI",     18.0,  0.025, 1.10
)

# --- 3B: FVS-NE annual growth engine ---
fvsne_grow_one_year <- function(trees, di_coefs, mort_coefs, hd_coefs,
                                 si_multiplier = 1.0) {
  if (nrow(trees) == 0) return(trees)

  trees <- trees |>
    mutate(
      dbh_in = dbh_cm * CM_TO_IN,
      ba_tree_ft2 = 0.005454 * dbh_in^2
    )

  plot_area_ac <- first(trees$plot_area_m2) / 4046.86
  trees <- trees |>
    arrange(desc(dbh_in)) |>
    mutate(
      cum_ba = cumsum(ba_tree_ft2),
      bal_ft2ac = (cum_ba - ba_tree_ft2) / plot_area_ac
    )

  trees <- trees |>
    left_join(di_coefs |> select(species_4, B1, B2, B3, default_si_ft),
              by = "species_4") |>
    left_join(mort_coefs |> select(species_4, p0, p1), by = "species_4") |>
    mutate(
      B1 = coalesce(B1, 0.0008236),
      B2 = coalesce(B2, 0.0550000),
      B3 = coalesce(B3, 0.012000),
      default_si_ft = coalesce(default_si_ft, 55),
      p0 = coalesce(p0, 5.1677),
      p1 = coalesce(p1, -0.0078)
    )

  # Apply site index multiplier for calibration
  trees <- trees |>
    mutate(
      effective_si = default_si_ft * si_multiplier,
      potbag = B1 * effective_si * (1 - exp(-B2 * dbh_in)) * 0.7,
      gmod = pmax(0.5, exp(-B3 * bal_ft2ac)),
      abag = pmax(potbag * gmod, 0),
      new_ba_ft2 = ba_tree_ft2 + abag,
      dgrow_in = pmax(sqrt(new_ba_ft2 / 0.005454) - dbh_in, 0),
      dbh_in_new = dbh_in + dgrow_in,
      dbh_cm_new = dbh_in_new / CM_TO_IN
    )

  # Background mortality
  trees <- trees |>
    mutate(
      ri = (1 / (1 + exp(p0 + p1 * dbh_in))) * 0.5,
      p_mort_annual = 1 - (1 - ri)^1,
      alive = runif(n()) > p_mort_annual
    )

  # Density-dependent mortality: when stand SDI > 55% of max SDI
  # SDI per acre = Σ TPA_i × (DBH_in_i / 10)^1.605
  tpa_per_tree <- 1 / plot_area_ac  # each tree record = 1/plot_area trees per acre
  stand_ba_ft2ac <- sum(trees$ba_tree_ft2) / plot_area_ac
  stand_sdi <- sum(tpa_per_tree * (trees$dbh_in / 10)^1.605)
  max_sdi <- 400  # approximate max SDI for NE spruce-fir
  if (stand_sdi > 0.55 * max_sdi) {
    # additional mortality proportional to excess density
    excess_ratio <- (stand_sdi - 0.55 * max_sdi) / (0.45 * max_sdi)
    excess_ratio <- min(excess_ratio, 1.0)
    # kill smallest trees preferentially
    trees <- trees |>
      mutate(
        dd_mort_prob = excess_ratio * 0.05 * (1 - dbh_in / max(dbh_in)),
        alive = alive & (runif(n()) > dd_mort_prob)
      )
  }

  # Height update from H-D curve
  trees <- trees |>
    left_join(hd_coefs, by = "species_4") |>
    mutate(
      hd_a = coalesce(hd_a, 20.0),
      hd_b = coalesce(hd_b, 0.028),
      hd_c = coalesce(hd_c, 1.10),
      ht_m_new = 1.3 + hd_a * (1 - exp(-hd_b * dbh_cm_new))^hd_c
    )

  trees |>
    filter(alive) |>
    mutate(dbh_cm = dbh_cm_new, ht_m = ht_m_new) |>
    select(any_of(c("tree_id", "treeID", "plot", "species_4", "dbh_cm", "ht_m",
                     "plot_area_m2")))
}

# --- 3C: Generic projection wrapper for FVS-NE ---
project_treelist_ne <- function(trees, n_years, di_coefs, mort_coefs, hd_coefs,
                                 si_multiplier = 1.0, choj = choj_coefs) {
  results <- list()

  # Initial AGB
  trees_agb <- trees |>
    left_join(choj, by = "species_4") |>
    mutate(
      choj_b0 = coalesce(choj_b0, -2.0773),
      choj_b1 = coalesce(choj_b1, 2.3323),
      agb_kg = exp(choj_b0 + choj_b1 * log(dbh_cm))
    )

  results[[1]] <- trees_agb |>
    group_by(plot) |>
    summarise(
      year_offset = 0, n_trees = n(), plot_area_m2 = first(plot_area_m2),
      agb_kg_total = sum(agb_kg, na.rm = TRUE),
      mean_dbh_cm = mean(dbh_cm, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(agb_Mg_ha = agb_kg_total / plot_area_m2 * 10,
           agb_tons_ac = agb_Mg_ha * MG_HA_TO_TONS_AC)

  current_trees <- trees

  for (yr in seq_len(n_years)) {
    if (yr %% 5 == 0) cat("  Year", yr, "of", n_years, "\n")
    current_trees <- fvsne_grow_one_year(
      current_trees, di_coefs, mort_coefs, hd_coefs,
      si_multiplier = si_multiplier
    )

    trees_agb <- current_trees |>
      left_join(choj, by = "species_4") |>
      mutate(
        choj_b0 = coalesce(choj_b0, -2.0773),
        choj_b1 = coalesce(choj_b1, 2.3323),
        agb_kg = exp(choj_b0 + choj_b1 * log(dbh_cm))
      )

    results[[yr + 1]] <- trees_agb |>
      group_by(plot) |>
      summarise(
        year_offset = yr, n_trees = n(), plot_area_m2 = first(plot_area_m2),
        agb_kg_total = sum(agb_kg, na.rm = TRUE),
        mean_dbh_cm = mean(dbh_cm, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(agb_Mg_ha = agb_kg_total / plot_area_m2 * 10,
             agb_tons_ac = agb_Mg_ha * MG_HA_TO_TONS_AC)
  }

  bind_rows(results)
}

# Prepare NE tree list
ne_trees <- hl_trees |>
  transmute(
    plot,
    species_4,
    dbh_cm,
    ht_m = NA_real_,
    plot_area_m2 = area_m2
  ) |>
  filter(!is.na(dbh_cm), dbh_cm > 0)

# --- Version 1: FVS-NE Out-of-Box (default SI) ---
cat("\nProjecting FVS-NE Out-of-Box (30 years, default SI)...\n")
set.seed(2026)
ne_oob_projection <- project_treelist_ne(
  ne_trees, n_years = n_years_long,
  di_coefs = fvsne_di_coefs, mort_coefs = fvsne_mort_coefs,
  hd_coefs = hd_params, si_multiplier = 1.0
)

ne_oob_means <- ne_oob_projection |>
  group_by(year_offset) |>
  summarise(mean_agb_Mg_ha = mean(agb_Mg_ha, na.rm = TRUE),
            sd_agb_Mg_ha = sd(agb_Mg_ha, na.rm = TRUE),
            n_plots = n_distinct(plot), .groups = "drop") |>
  mutate(year = HOWLAND_YEAR + year_offset,
         mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
         method = "FVS-NE (defaults)")

# --- Version 2: FVS-NE Locally Calibrated ---
# Back-solve for SI multiplier using observed 2012-2025 AGB trajectory
# Load observed AGB from Script 01 outputs
obs_rds_path <- file.path(output_dir, "agb_trajectories.rds")
obs_agb_2025 <- NA_real_
obs_agb_2015 <- NA_real_

if (file.exists(obs_rds_path)) {
  traj <- readRDS(obs_rds_path)
  obs_summary <- traj |>
    group_by(year) |>
    summarise(mean_agb = mean(agb_Mg_ha, na.rm = TRUE), .groups = "drop")
  obs_agb_2025 <- obs_summary |> filter(year == 2025) |> pull(mean_agb)
  obs_agb_2015 <- obs_summary |> filter(year == 2015) |> pull(mean_agb)
  if (length(obs_agb_2025) == 0) obs_agb_2025 <- NA_real_
  if (length(obs_agb_2015) == 0) obs_agb_2015 <- NA_real_
  cat("Observed AGB 2015:", round(obs_agb_2015, 1), "Mg/ha\n")
  cat("Observed AGB 2025:", round(obs_agb_2025, 1), "Mg/ha\n")
}

# Calibration: find SI multiplier that matches observed 10yr change
# Use binary search
cat("\nCalibrating FVS-NE site index multiplier...\n")

ne_oob_10yr <- ne_oob_means |> filter(year_offset == 10) |> pull(mean_agb_Mg_ha)
ne_oob_0yr  <- ne_oob_means |> filter(year_offset == 0)  |> pull(mean_agb_Mg_ha)

if (!is.na(obs_agb_2025) && !is.na(obs_agb_2015)) {
  observed_change <- obs_agb_2025 - obs_agb_2015
  # The OOB change using field-plot-derived starting AGB
  oob_change <- ne_oob_10yr - ne_oob_0yr
  # Simple proportional calibration: SI_mult = observed_change / oob_change
  # But we need to handle the case where changes are very different
  # Use a grid search for robustness
  si_grid <- seq(0.2, 1.5, by = 0.05)
  grid_results <- map_dbl(si_grid, function(mult) {
    set.seed(2026)
    proj <- project_treelist_ne(
      ne_trees, n_years = 10,
      di_coefs = fvsne_di_coefs, mort_coefs = fvsne_mort_coefs,
      hd_coefs = hd_params, si_multiplier = mult
    )
    proj_means <- proj |>
      group_by(year_offset) |>
      summarise(m = mean(agb_Mg_ha, na.rm = TRUE), .groups = "drop")
    pred_10yr <- proj_means |> filter(year_offset == 10) |> pull(m)
    pred_0yr  <- proj_means |> filter(year_offset == 0)  |> pull(m)
    pred_change <- pred_10yr - pred_0yr
    abs(pred_change - observed_change)
  })

  best_mult <- si_grid[which.min(grid_results)]
  cat("  Best SI multiplier:", best_mult, "\n")
  cat("  Observed 10yr AGB change:", round(observed_change, 1), "Mg/ha\n")
  cat("  NOTE: LiDAR baseline AGB =", round(obs_agb_2015, 1),
      "Mg/ha vs field plot baseline =", round(ne_oob_0yr, 1), "Mg/ha\n")
  cat("  Absolute 2025 error will reflect this baseline mismatch.\n")
  if (best_mult >= max(si_grid) - 0.05) {
    cat("  WARNING: SI multiplier hit grid ceiling (", max(si_grid),
        "). Consider extending si_grid range.\n")
  }
} else {
  # Default: reduce SI by 40% (common for NE spruce-fir)
  best_mult <- 0.60
  cat("  No observed trajectory available; using default SI multiplier:", best_mult, "\n")
}

cat("Projecting FVS-NE Calibrated (30 years, SI mult =", best_mult, ")...\n")
set.seed(2026)
ne_cal_projection <- project_treelist_ne(
  ne_trees, n_years = n_years_long,
  di_coefs = fvsne_di_coefs, mort_coefs = fvsne_mort_coefs,
  hd_coefs = hd_params, si_multiplier = best_mult
)

ne_cal_means <- ne_cal_projection |>
  group_by(year_offset) |>
  summarise(mean_agb_Mg_ha = mean(agb_Mg_ha, na.rm = TRUE),
            sd_agb_Mg_ha = sd(agb_Mg_ha, na.rm = TRUE),
            n_plots = n_distinct(plot), .groups = "drop") |>
  mutate(year = HOWLAND_YEAR + year_offset,
         mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
         method = "FVS-NE (calibrated)")

cat("\nFVS-NE Out-of-Box: Year 0 =", round(ne_oob_means$mean_agb_Mg_ha[1], 1),
    "-> Year 10 =", round(ne_oob_10yr, 1), "Mg/ha\n")
cat("FVS-NE Calibrated: Year 0 =", round(ne_cal_means$mean_agb_Mg_ha[1], 1),
    "-> Year 10 =", round(ne_cal_means$mean_agb_Mg_ha[ne_cal_means$year_offset == 10], 1),
    "Mg/ha\n")

# =============================================================================
# PART 4: Build Yield Curves from FVS-ACD Projections
# =============================================================================

cat("\n========== PART 4: Yield Curves from FVS-ACD ==========\n")

# Stratify plots by initial AGB class
plot_initial <- acd_trajectory |>
  filter(year_offset == 0) |>
  mutate(
    agb_class = case_when(
      agb_Mg_ha < quantile(agb_Mg_ha, 0.33) ~ "Low",
      agb_Mg_ha < quantile(agb_Mg_ha, 0.67) ~ "Medium",
      TRUE ~ "High"
    )
  ) |>
  select(PLOT, agb_class, initial_agb = agb_Mg_ha)

# Join class to full trajectory
acd_with_class <- acd_trajectory |>
  left_join(plot_initial |> select(PLOT, agb_class, initial_agb), by = "PLOT")

# Fit Chapman-Richards yield curves by class
# AGB(t) = AGB_init + Amax * (1 - exp(-k * t))^p
# where t = years from baseline

yield_curves <- acd_with_class |>
  group_by(agb_class) |>
  nest() |>
  mutate(
    fit = map(data, function(d) {
      tryCatch({
        nls(agb_Mg_ha ~ initial_agb + Amax * (1 - exp(-k * year_offset))^p,
            data = d,
            start = list(Amax = 30, k = 0.05, p = 1.5),
            control = nls.control(maxiter = 200, warnOnly = TRUE))
      }, error = function(e) {
        # Fallback: linear model
        lm(agb_Mg_ha ~ year_offset + initial_agb, data = d)
      })
    }),
    params = map(fit, ~coef(.x))
  )

cat("Yield curve parameters by AGB class:\n")
for (i in seq_len(nrow(yield_curves))) {
  cat("  ", yield_curves$agb_class[i], ":",
      paste(names(yield_curves$params[[i]]),
            round(yield_curves$params[[i]], 4), sep = "=", collapse = ", "), "\n")
}

# Pre-compute class break thresholds (ONCE, not per pixel)
class_breaks <- plot_initial |>
  group_by(agb_class) |>
  summarise(mean_agb = mean(initial_agb), .groups = "drop")

low_med_break  <- mean(class_breaks$mean_agb[class_breaks$agb_class %in% c("Low", "Medium")])
med_high_break <- mean(class_breaks$mean_agb[class_breaks$agb_class %in% c("Medium", "High")])
cat("  Class breaks: Low <", round(low_med_break, 1),
    "< Medium <", round(med_high_break, 1), "< High\n")

# Vectorized prediction function for a vector of initial AGB values at a single time
predict_yield_vec <- function(agb_vec, years, yield_fits) {
  # Assign all pixels to classes at once
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

# Keep scalar version as fallback for NAIP proxy
predict_yield <- function(initial_agb, years, yield_fits) {
  agb_class <- case_when(
    initial_agb < low_med_break  ~ "Low",
    initial_agb < med_high_break ~ "Medium",
    TRUE ~ "High"
  )
  fit <- yield_fits$fit[yield_fits$agb_class == agb_class][[1]]
  pred_data <- data.frame(year_offset = years, initial_agb = initial_agb)
  tryCatch(
    predict(fit, newdata = pred_data),
    error = function(e) rep(initial_agb, length(years))
  )
}

# =============================================================================
# PART 5: Apply Yield Curves to LiDAR and NAIP Pixels
# =============================================================================

cat("\n========== PART 5: Pixel-Level Projections (LiDAR + NAIP) ==========\n")

# Load 2015 AGB raster (3DEP baseline) for pixel-level projection
agb_2015_path <- file.path(raster_dir, "Howland_AGB_2015.tif")
agb_2025_path <- file.path(raster_dir, "Howland_AGB_2025.tif")

lidar_pixel_results <- NULL
naip_pixel_results  <- NULL

if (file.exists(agb_2015_path)) {
  agb_2015_r <- rast(agb_2015_path)
  # Raster is in kgC/m2; convert to Mg biomass/ha:
  #   kgC/m2 * 10000 m2/ha / 1000 kg/Mg * 2 (C to total biomass) = * 20
  agb_vals_2015 <- values(agb_2015_r, na.rm = FALSE) * 20
  valid_mask <- !is.na(agb_vals_2015) & agb_vals_2015 > 0

  n_valid <- sum(valid_mask)
  cat("LiDAR 2015 raster:", n_valid, "valid pixels\n")

  # Project each pixel forward using yield curves (vectorized)
  if (n_valid > 0) {
    agb_valid <- agb_vals_2015[valid_mask]

    # Predict at 10, 20, and 50 years (vectorized for speed)
    cat("  Projecting", length(agb_valid), "pixels at 10, 20, 50 yr...\n")
    proj_10yr <- predict_yield_vec(agb_valid, 10, yield_curves)
    proj_20yr <- predict_yield_vec(agb_valid, 20, yield_curves)
    proj_50yr <- predict_yield_vec(agb_valid, 50, yield_curves)

    lidar_pixel_results <- tibble(
      pixel_agb_2015 = agb_valid,
      pred_agb_2025  = proj_10yr,
      pred_agb_2035  = proj_20yr,
      pred_agb_2065  = proj_50yr
    )

    cat("  LiDAR+ACD yield curve: mean 2015 =", round(mean(agb_valid), 1),
        "-> mean 2025 =", round(mean(proj_10yr), 1),
        "-> mean 2035 =", round(mean(proj_20yr), 1),
        "-> mean 2065 =", round(mean(proj_50yr), 1), "Mg/ha\n")
  }
} else {
  cat("2015 LiDAR raster not found at:", agb_2015_path, "\n")
}

# NAIP-based projection (using Script 04's spectral model predictions)
naip_agb_path <- file.path(output_dir, "naip_agb_predictions.rds")
if (file.exists(naip_agb_path)) {
  naip_agb <- readRDS(naip_agb_path)
  cat("NAIP AGB predictions loaded:", nrow(naip_agb), "pixels\n")

  naip_proj_10yr <- predict_yield_vec(naip_agb$pred_agb, 10, yield_curves)
  naip_proj_20yr <- predict_yield_vec(naip_agb$pred_agb, 20, yield_curves)
  naip_proj_50yr <- predict_yield_vec(naip_agb$pred_agb, 50, yield_curves)

  naip_pixel_results <- tibble(
    pixel_agb_naip = naip_agb$pred_agb,
    pred_agb_2025  = naip_proj_10yr,
    pred_agb_2035  = naip_proj_20yr,
    pred_agb_2065  = naip_proj_50yr
  )

  cat("  NAIP+ACD yield curve: mean NAIP =", round(mean(naip_agb$pred_agb), 1),
      "-> mean 2025 =", round(mean(naip_proj_10yr), 1), "Mg/ha\n")
} else {
  cat("NAIP AGB predictions not found. Using plot-level proxy.\n")
  # Create proxy NAIP results from plot data with added noise
  proxy_agb <- pmax(plot_initial$initial_agb + rnorm(nrow(plot_initial), 0, 15), 5)
  naip_pixel_results <- tibble(
    pixel_agb_naip = proxy_agb,
    pred_agb_2025  = predict_yield_vec(proxy_agb, 10, yield_curves),
    pred_agb_2035  = predict_yield_vec(proxy_agb, 20, yield_curves),
    pred_agb_2065  = predict_yield_vec(proxy_agb, 50, yield_curves)
  )
}

# =============================================================================
# PART 6: Assemble All Methods for Comparison
# =============================================================================

cat("\n========== PART 6: Multi-Method Comparison ==========\n")

# Load area-based predictions from Script 04
val_rds_path <- file.path(output_dir, "projection_validation_results.rds")
stratum_means <- tibble()
boot_cr <- NULL
boot_schum <- NULL

if (file.exists(val_rds_path)) {
  val_results <- readRDS(val_rds_path)
  stratum_means <- val_results$stratum_means
  boot_cr <- val_results$boot_cr
  boot_schum <- val_results$boot_schum
  cat("Loaded area-based predictions from Script 04\n")
} else {
  rds_path <- file.path(output_dir, "agb_trajectories.rds")
  if (file.exists(rds_path)) {
    traj <- readRDS(rds_path)
    stratum_means <- traj |>
      group_by(year) |>
      summarise(mean_agb = mean(agb_Mg_ha, na.rm = TRUE), .groups = "drop")
  }
}

# Observed LiDAR means
area_based_means <- stratum_means |>
  group_by(year) |>
  summarise(mean_agb_Mg_ha = weighted.mean(mean_agb, na.rm = TRUE),
            .groups = "drop") |>
  mutate(mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
         method = "Observed LiDAR")

# Area-based CR
cr_means <- if (!is.null(boot_cr)) {
  boot_cr |>
    group_by(year) |>
    summarise(mean_agb_Mg_ha = weighted.mean(pred_mean, na.rm = TRUE),
              .groups = "drop") |>
    mutate(mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
           method = "Area-based CR")
} else { tibble() }

# Area-based Schumacher
schum_means <- if (!is.null(boot_schum)) {
  boot_schum |>
    group_by(year) |>
    summarise(mean_agb_Mg_ha = weighted.mean(pred_mean, na.rm = TRUE),
              .groups = "drop") |>
    mutate(mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
           method = "Area-based Schum.")
} else { tibble() }

# LiDAR + ACD yield curve means
lidar_yc_means <- tibble()
if (!is.null(lidar_pixel_results)) {
  lidar_yc_means <- tibble(
    year = c(2015, 2025, 2035, 2065),
    mean_agb_Mg_ha = c(
      mean(lidar_pixel_results$pixel_agb_2015),
      mean(lidar_pixel_results$pred_agb_2025),
      mean(lidar_pixel_results$pred_agb_2035),
      mean(lidar_pixel_results$pred_agb_2065)
    )
  ) |>
    mutate(mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
           method = "LiDAR + ACD curve")
}

# NAIP + ACD yield curve means
naip_yc_means <- tibble()
if (!is.null(naip_pixel_results)) {
  naip_yc_means <- tibble(
    year = c(2015, 2025, 2035, 2065),
    mean_agb_Mg_ha = c(
      mean(naip_pixel_results$pixel_agb_naip),
      mean(naip_pixel_results$pred_agb_2025),
      mean(naip_pixel_results$pred_agb_2035),
      mean(naip_pixel_results$pred_agb_2065)
    )
  ) |>
    mutate(mean_agb_tons_ac = mean_agb_Mg_ha * MG_HA_TO_TONS_AC,
           method = "NAIP + ACD curve")
}

# Combine all methods
comparison_df <- bind_rows(
  area_based_means,
  cr_means,
  schum_means,
  ne_oob_means |> select(year, mean_agb_Mg_ha, mean_agb_tons_ac, method),
  ne_cal_means |> select(year, mean_agb_Mg_ha, mean_agb_tons_ac, method),
  acd_means |> select(year, mean_agb_Mg_ha, mean_agb_tons_ac, method),
  lidar_yc_means,
  naip_yc_means
) |>
  filter(!is.na(mean_agb_Mg_ha))

cat("\n--- Multi-Method Comparison at Key Years ---\n")
comparison_df |>
  filter(year %in% c(2015, 2025, 2035, 2045, 2055, 2065)) |>
  pivot_wider(names_from = method, values_from = c(mean_agb_Mg_ha, mean_agb_tons_ac)) |>
  print()

# =============================================================================
# PART 7: Keynote Figures
# =============================================================================

cat("\n========== PART 7: Generating Figures ==========\n")

# Color palette for all methods
pal_methods <- c(
  "Observed LiDAR"       = "black",
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

# --- Figure A: Multi-Method Comparison (metric, 10yr) ---
fig_comp_10yr <- comparison_df |>
  filter(year >= 2012, year <= 2026) |>
  ggplot(aes(x = year, y = mean_agb_Mg_ha, color = method, shape = method)) +
  annotate("rect", xmin = -Inf, xmax = 2015.5,
           ymin = -Inf, ymax = Inf, fill = "#E8F5E9", alpha = 0.4) +
  annotate("text", x = 2013, y = Inf,
           label = "Calibration", vjust = 1.5, size = 3.5, color = "grey40") +
  annotate("text", x = 2020, y = Inf,
           label = "Projection", vjust = 1.5, size = 3.5, color = "grey40") +
  geom_vline(xintercept = 2015.5, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = pal_methods, name = NULL) +
  scale_shape_manual(values = shape_methods, name = NULL) +
  labs(x = "Year",
       y = "Mean AGB (Mg ha<sup>\u22121</sup>)",
       title = "10 Year Projection: Multi-Method Comparison") +
  theme_pub + theme(axis.title.y = element_markdown())

ggsave(file.path(fig_dir, "fig_multimethod_10yr.png"), fig_comp_10yr,
       width = 24, height = 16, units = "cm", dpi = 300, bg = "white")
cat("Saved: fig_multimethod_10yr.png\n")

# --- Figure B: Long-term projection (30yr) ---
fig_comp_long <- comparison_df |>
  filter(year >= 2015) |>
  ggplot(aes(x = year, y = mean_agb_Mg_ha, color = method, shape = method)) +
  geom_vline(xintercept = 2025, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  annotate("text", x = 2025, y = Inf, label = "Validation year (2025)",
           vjust = 1.5, hjust = -0.05, size = 3, color = "grey40") +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal_methods, name = NULL) +
  scale_shape_manual(values = shape_methods, name = NULL) +
  labs(x = "Year",
       y = "Mean AGB (Mg ha<sup>\u22121</sup>)",
       title = "50 Year Projection: Growth Model Divergence") +
  theme_pub + theme(axis.title.y = element_markdown())

ggsave(file.path(fig_dir, "fig_multimethod_30yr.png"), fig_comp_long,
       width = 24, height = 16, units = "cm", dpi = 300, bg = "white")
cat("Saved: fig_multimethod_30yr.png\n")

# --- Figure C: Error at 2025 bar chart ---
if (!is.na(obs_agb_2025) && length(obs_agb_2025) > 0) {
  obs_2025 <- obs_agb_2025

  error_at_2025 <- comparison_df |>
    filter(year == 2025, method != "Observed LiDAR") |>
    mutate(
      obs_2025 = obs_2025,
      error_Mg_ha = mean_agb_Mg_ha - obs_2025,
      error_pct = 100 * error_Mg_ha / obs_2025,
      error_tons_ac = error_Mg_ha * MG_HA_TO_TONS_AC
    )

  cat("\n--- 2025 Prediction Error ---\n")
  error_at_2025 |>
    select(method, mean_agb_Mg_ha, obs_2025, error_Mg_ha, error_pct) |>
    print()

  fig_error <- ggplot(error_at_2025,
                       aes(x = reorder(method, error_pct), y = error_pct,
                           fill = method)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_text(aes(label = paste0(round(error_pct, 1), "%")),
              vjust = ifelse(error_at_2025$error_pct >= 0, -0.5, 1.5),
              size = 4, fontface = "bold") +
    scale_fill_manual(values = pal_methods, guide = "none") +
    labs(x = NULL, y = "Prediction error at 2025 (%)",
         title = "2025 AGB Prediction Error: Projected vs. Phoenix RANGER") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(fig_dir, "fig_error_at_2025.png"), fig_error,
         width = 24, height = 14, units = "cm", dpi = 300, bg = "white")
  cat("Saved: fig_error_at_2025.png\n")
}

# --- Figure D: Relative change from baseline ---
baseline_by_method <- comparison_df |>
  group_by(method) |>
  filter(year == min(year)) |>
  select(method, baseline_agb = mean_agb_Mg_ha) |>
  ungroup()

relative_change <- comparison_df |>
  left_join(baseline_by_method, by = "method") |>
  mutate(pct_change = 100 * (mean_agb_Mg_ha - baseline_agb) / baseline_agb)

fig_relative <- ggplot(relative_change |> filter(year >= 2015),
                        aes(x = year, y = pct_change,
                            color = method, shape = method)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal_methods, name = NULL) +
  scale_shape_manual(values = shape_methods, name = NULL) +
  labs(x = "Year", y = "Change from baseline (%)",
       title = "Relative AGB Change by Method (50 year)") +
  theme_pub

ggsave(file.path(fig_dir, "fig_relative_change.png"), fig_relative,
       width = 24, height = 16, units = "cm", dpi = 300, bg = "white")
cat("Saved: fig_relative_change.png\n")

# =============================================================================
# PART 8: Enhanced Monte Carlo Financial Simulation
# =============================================================================

cat("\n========== PART 8: Enhanced Monte Carlo (10yr + 20yr) ==========\n")

cost_data <- tribble(
  ~method,                 ~cost_lo, ~cost_hi, ~cost_unit,  ~coverage,
  "NAIP + ACD curve",       0.10,     0.25,    "per acre",  "Wall to wall",
  "Area-based LiDAR",       0.25,     0.50,    "per acre",  "Wall to wall",
  "LiDAR + ACD curve",      0.30,     0.55,    "per acre",  "Wall to wall",
  "ITC + FVS-NE",           0.75,     1.00,    "per acre",  "Wall to wall",
  "ITC + FVS-ACD",          0.75,     1.00,    "per acre",  "Wall to wall",
  # Field inventory: 1 plot per 5-10 acres, ~$5-8k total for 865 ac study area
  "Field inventory",        5.78,     9.25,    "per acre",  "Sample based"
)

# Accuracy estimates (RMSE % of mean AGB)
acc_10yr <- c(
  "NAIP + ACD curve"  = 22,
  "Area-based LiDAR"  = 15,
  "LiDAR + ACD curve" = 12,
  "ITC + FVS-NE"      = 35,
  "ITC + FVS-ACD"     = 10,
  "Field inventory"   = 8
)

# Long-term accuracy degrades differently by method
acc_20yr <- c(
  "NAIP + ACD curve"  = 28,
  "Area-based LiDAR"  = 25,
  "LiDAR + ACD curve" = 16,
  "ITC + FVS-NE"      = 50,
  "ITC + FVS-ACD"     = 14,
  "Field inventory"   = 12
)

# 50-year accuracy: extrapolation uncertainty grows substantially
acc_50yr <- c(
  "NAIP + ACD curve"  = 40,
  "Area-based LiDAR"  = 38,
  "LiDAR + ACD curve" = 22,
  "ITC + FVS-NE"      = 70,
  "ITC + FVS-ACD"     = 18,
  "Field inventory"   = 15
)

# Update accuracy from pipeline results where available
if (!is.na(obs_agb_2025) && length(obs_agb_2025) > 0) {
  # FVS-ACD full model error at 2025
  acd_pred_2025 <- acd_means |> filter(year == 2025) |> pull(mean_agb_Mg_ha)
  if (length(acd_pred_2025) > 0 && !is.na(acd_pred_2025)) {
    acd_err <- 100 * abs(acd_pred_2025 - obs_agb_2025) / obs_agb_2025
    acc_10yr["ITC + FVS-ACD"] <- acd_err
    cat("FVS-ACD 2025 error:", round(acd_err, 1), "% (updated from pipeline)\n")
  }

  # FVS-NE OOB error
  ne_pred_2025 <- ne_oob_means |> filter(year == 2025) |> pull(mean_agb_Mg_ha)
  if (length(ne_pred_2025) > 0 && !is.na(ne_pred_2025)) {
    ne_err <- 100 * abs(ne_pred_2025 - obs_agb_2025) / obs_agb_2025
    acc_10yr["ITC + FVS-NE"] <- ne_err
    cat("FVS-NE OOB 2025 error:", round(ne_err, 1), "% (updated from pipeline)\n")
  }
}

cat("Final 10yr accuracy (RMSE%):\n")
print(acc_10yr)
cat("Final 20yr accuracy (RMSE%):\n")
print(acc_20yr)
cat("Final 50yr accuracy (RMSE%):\n")
print(acc_50yr)

# Monte Carlo simulation
set.seed(2026)
n_mc <- 1000

run_mc <- function(acc_vec, horizon_label) {
  map_dfr(seq_len(nrow(cost_data)), function(i) {
    m <- cost_data$method[i]
    acc_val <- acc_vec[m]
    if (is.na(acc_val)) acc_val <- 25  # default
    tibble(
      method = m, sim = 1:n_mc, horizon = horizon_label,
      cost = runif(n_mc, cost_data$cost_lo[i], cost_data$cost_hi[i]),
      rmse_pct = pmax(1, rnorm(n_mc, acc_val, acc_val * 0.2)),
      coverage = cost_data$coverage[i]
    )
  })
}

mc_10yr <- run_mc(acc_10yr, "10 year")
mc_20yr <- run_mc(acc_20yr, "20 year")
mc_50yr <- run_mc(acc_50yr, "50 year")
mc_all  <- bind_rows(mc_10yr, mc_20yr, mc_50yr)

mc_summary <- mc_all |>
  group_by(method, coverage, horizon) |>
  summarise(
    cost_mean = mean(cost), rmse_mean = mean(rmse_pct),
    rmse_lo = quantile(rmse_pct, 0.025), rmse_hi = quantile(rmse_pct, 0.975),
    .groups = "drop"
  )

# --- Figure E: Cost vs Accuracy (faceted by horizon) ---
fig_mc <- ggplot(mc_all, aes(x = cost, y = rmse_pct, color = method)) +
  geom_point(alpha = 0.05, size = 0.3) +
  stat_ellipse(level = 0.90, linewidth = 1.0) +
  geom_point(data = mc_summary,
             aes(x = cost_mean, y = rmse_mean, shape = coverage),
             size = 4, stroke = 1.2) +
  facet_wrap(~horizon, scales = "free_y") +
  scale_x_log10(labels = scales::dollar_format()) +
  scale_color_brewer(palette = "Set2", name = NULL) +
  scale_shape_manual(values = c("Wall to wall" = 16, "Sample based" = 17),
                     name = "Coverage") +
  labs(x = "Cost per acre (USD, log scale)", y = "Expected RMSE (%)",
       title = "Cost vs. Accuracy: 10, 20, and 50 Year Projections (n = 1,000)") +
  theme_pub + theme(legend.box = "vertical")

ggsave(file.path(fig_dir, "fig_monte_carlo_enhanced.png"), fig_mc,
       width = 28, height = 16, units = "cm", dpi = 300, bg = "white")
cat("Saved: fig_monte_carlo_enhanced.png\n")

# --- Figure F: 4-Panel Keynote Combined ---
if (exists("fig_comp_10yr") && exists("fig_relative") &&
    exists("fig_error") && exists("fig_mc")) {
  fig_keynote <- (fig_comp_10yr + fig_error) / (fig_relative + fig_mc) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Howland Research Forest: Multi-Method AGB Projection Comparison",
      tag_levels = "A"
    ) &
    theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))

  ggsave(file.path(fig_dir, "fig_keynote_combined.png"), fig_keynote,
         width = 36, height = 28, units = "cm", dpi = 300, bg = "white")
  cat("Saved: fig_keynote_combined.png\n")
}

# Howland study area cost estimates
howland_acres <- 865
total_cost <- cost_data |>
  mutate(
    total_low = ifelse(cost_unit == "per acre",
                       cost_lo * howland_acres, cost_lo * 150),
    total_high = ifelse(cost_unit == "per acre",
                        cost_hi * howland_acres, cost_hi * 150)
  )

cat("\n--- Total Cost for Howland (", howland_acres, "acres) ---\n")
total_cost |>
  select(method, cost_lo, cost_hi, cost_unit, total_low, total_high) |>
  mutate(total_low = scales::dollar(total_low),
         total_high = scales::dollar(total_high)) |>
  print()

# =============================================================================
# PART 9: Save All Results
# =============================================================================

cat("\n========== PART 9: Saving Results ==========\n")

# Save comparison data for downstream use
saveRDS(
  list(
    comparison_df  = comparison_df,
    acd_trajectory = acd_trajectory,
    ne_oob_means   = ne_oob_means,
    ne_cal_means   = ne_cal_means,
    acd_means      = acd_means,
    yield_curves   = yield_curves,
    mc_summary     = mc_summary,
    si_multiplier  = best_mult,
    howland_csi    = HOWLAND_CSI,
    lidar_pixel    = lidar_pixel_results,
    naip_pixel     = naip_pixel_results
  ),
  file = file.path(output_dir, "fvs_projection_results.rds")
)
cat("Saved: fvs_projection_results.rds\n")

# Save comparison table as CSV
write_csv(comparison_df, file.path(output_dir, "method_comparison_table.csv"))
cat("Saved: method_comparison_table.csv\n")

cat("\n========== Script 05 Complete ==========\n")
cat("Methods compared:", n_distinct(comparison_df$method), "\n")
cat("Projection horizon: 50 years\n")
cat("FVS-ACD version:", AcadianVersionTag, "\n")
cat("FVS-NE SI multiplier:", best_mult, "\n")
cat("Figures saved to:", fig_dir, "\n")
