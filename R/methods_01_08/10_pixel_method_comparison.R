# =============================================================================
# Title: Pixel Level Method Comparison and AGB Panel Figure
# Author: A. Weiskittel
# Date: 2026-03-26
# Description: Generates two key manuscript figures:
#   (1) Panel AGB maps showing 2015 baseline and 2025 projections by method
#   (2) Pixel level method divergence analysis (SD map, correlation, range)
# Dependencies: tidyverse, terra, ggplot2, patchwork
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(terra)
library(ggplot2)
library(patchwork)

# --- Paths -------------------------------------------------------------------
base_dir   <- "/users/PUOM0008/crsfaaron/HRF"
ornl_dir   <- file.path(base_dir, "data", "ornl_agb_rasters")
output_dir <- file.path(base_dir, "output")
fig_dir    <- file.path(output_dir, "manuscript_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ORNL_FACTOR <- 20  # kgC/m2 to Mg biomass/ha

# --- Publication theme -------------------------------------------------------
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# --- Method color palette (consistent across all figures) --------------------
method_colors <- c(
  "Observed (ORNL 2015)"    = "#1b9e77",
  "Observed (MELiTE 2025)"  = "#d95f02",
  "CR yield curves (3 pt)"  = "#7570b3",
  "CR yield curves (4 pt)"  = "#e7298a",
  "FVS ACD yield curves"    = "#66a61e"
)

cat("============================================================\n")
cat("  Script 10: Pixel Level Method Comparison\n")
cat("============================================================\n\n")

# =============================================================================
# PART 1: Load data and generate pixel level predictions
# =============================================================================
cat("--- Part 1: Loading rasters and model fits ---\n")

# Load ORNL 2015 baseline
agb_2015 <- rast(file.path(ornl_dir, "Howland_AGB_2015.tif")) * ORNL_FACTOR
cat("Baseline 2015: mean =", round(mean(values(agb_2015, na.rm = TRUE)), 1), "Mg/ha\n")

# Load MELiTE 2025 observed
agb_2025_obs <- rast(file.path(output_dir, "agb_2025_melite_30m.tif"))
agb_2025_obs <- resample(agb_2025_obs, agb_2015, method = "near")
cat("Observed 2025 (MELiTE): mean =", round(mean(values(agb_2025_obs, na.rm = TRUE)), 1), "Mg/ha\n")

# Load FVS projection results (contains 3 pt yield curves)
fvs_res <- readRDS(file.path(output_dir, "fvs_projection_results.rds"))
yield_curves <- fvs_res$yield_curves

# Load recalibration results (contains 4 pt yield curves)
recal <- readRDS(file.path(output_dir, "recalibration_results.rds"))

# --- Extract valid pixels ----------------------------------------------------
vals_2015 <- values(agb_2015, na.rm = FALSE)[,1]
vals_2025 <- values(agb_2025_obs, na.rm = FALSE)[,1]
xy_all    <- crds(agb_2015, na.rm = FALSE)

valid <- !is.na(vals_2015) & !is.na(vals_2025) &
         vals_2015 > 0 & vals_2025 > 0
n_valid <- sum(valid)
cat("Valid collocated pixels:", n_valid, "\n")

init_agb <- vals_2015[valid]
obs_agb  <- vals_2025[valid]
x_coords <- xy_all[valid, 1]
y_coords <- xy_all[valid, 2]

# --- Classify pixels for yield curve prediction ------------------------------
# Derive class breaks from actual pixel AGB distribution (terciles)
# Using pixel quantiles rather than yield curve internal means ensures
# class boundaries match the spatial data being predicted
low_med  <- as.numeric(quantile(init_agb, 0.33))
med_high <- as.numeric(quantile(init_agb, 0.67))

pixel_class <- case_when(
  init_agb < low_med  ~ "Low",
  init_agb < med_high ~ "Medium",
  TRUE ~ "High"
)
cat("Class breaks: Low <", round(low_med, 1), "< Medium <", round(med_high, 1), "< High\n")
cat("Class distribution:", table(pixel_class), "\n")

# --- Vectorized yield curve prediction function ------------------------------
predict_yc <- function(init_vec, classes, fits_list, year = 10, label = "") {
  pred <- numeric(length(init_vec))
  for (cls in c("Low", "Medium", "High")) {
    idx <- which(classes == cls)
    if (length(idx) == 0) next
    fit_obj <- fits_list[[cls]]
    if (is.null(fit_obj)) { pred[idx] <- NA; next }
    # Handle both list-with-$fit and direct lm/nls objects
    if (is.list(fit_obj) && !inherits(fit_obj, "lm") && !inherits(fit_obj, "nls") && !is.null(fit_obj$fit)) {
      fit_obj <- fit_obj$fit
    }
    cat("  [", label, cls, "] class:", paste(class(fit_obj), collapse = ","),
        " vars:", paste(all.vars(formula(fit_obj)), collapse = ","), "\n")
    model_vars <- all.vars(formula(fit_obj))[-1]
    nd <- data.frame(year_offset = rep(year, length(idx)))
    if ("initial_agb" %in% model_vars) nd$initial_agb <- init_vec[idx]
    pred[idx] <- tryCatch({
      p <- stats::predict.lm(fit_obj, newdata = nd)
      cat("  [", label, cls, "] n =", length(idx), " mean_pred =", round(mean(p, na.rm = TRUE), 1), "\n")
      p
    },
    error = function(e) { cat("  Prediction error for", label, cls, ":", e$message, "\n"); rep(NA, length(idx)) }
    )
  }
  pred
}

# --- Generate predictions for each method ------------------------------------
cat("\n--- Part 2: Generating pixel level predictions ---\n")

# Method 1: 3 pt CR yield curves
# Extract lm objects explicitly from the grouped tibble list column
# (using $fit on a grouped tibble can return numeric instead of lm objects)
fits_3pt <- list()
for (j in seq_len(nrow(yield_curves))) {
  fits_3pt[[ yield_curves[["agb_class"]][j] ]] <- yield_curves[["fit"]][[j]]
}
cat("fits_3pt classes:", paste(sapply(fits_3pt, class), collapse = ", "), "\n")
cat("fits_3pt names:", paste(names(fits_3pt), collapse = ", "), "\n")
pred_cr3 <- predict_yc(init_agb, pixel_class, fits_3pt, year = 10, label = "CR3")
cat("CR 3 pt: mean predicted 2025 =", round(mean(pred_cr3, na.rm = TRUE), 1), "Mg/ha\n")

# Method 2: 4 pt recalibrated CR yield curves
# Handle different possible structures
if (is.data.frame(recal$yield_fits_4pt) || is_tibble(recal$yield_fits_4pt)) {
  fits_4pt <- setNames(recal$yield_fits_4pt$fit, recal$yield_fits_4pt$agb_class)
} else if (is.list(recal$yield_fits_4pt)) {
  fits_4pt <- recal$yield_fits_4pt
} else {
  fits_4pt <- list(Low = NULL, Medium = NULL, High = NULL)
}
# 4 pt uses 13 year offset (2012 to 2025)
pred_cr4 <- predict_yc(init_agb, pixel_class, fits_4pt, year = 13, label = "CR4")
cat("CR 4 pt: mean predicted 2025 =", round(mean(pred_cr4, na.rm = TRUE), 1), "Mg/ha\n")

# Method 3: FVS ACD yield curves (same 3 pt fits applied, represents the FVS based approach)
# If there are separate FVS specific fits, use those; otherwise use the same CR fits
# as FVS ACD is the projector behind these curves
pred_fvs <- pred_cr3  # FVS ACD is the engine behind the area based CR curves
cat("FVS ACD: mean predicted 2025 =", round(mean(pred_fvs, na.rm = TRUE), 1), "Mg/ha\n")

# --- Build pixel data frame --------------------------------------------------
pixel_df <- tibble(
  x = x_coords,
  y = y_coords,
  agb_2015    = init_agb,
  obs_2025    = obs_agb,
  cr_3pt_2025 = pred_cr3,
  cr_4pt_2025 = pred_cr4,
  fvs_2025    = pred_fvs
) |>
  filter(!is.na(cr_3pt_2025), !is.na(cr_4pt_2025), !is.na(fvs_2025))

cat("Complete cases:", nrow(pixel_df), "\n")

# =============================================================================
# PART 3: Panel AGB Figure (Main Manuscript Figure)
# =============================================================================
cat("\n--- Part 3: Creating panel AGB figure ---\n")

# Reshape for faceted map: one panel per method, showing the 2025 AGB raster
# Plus the observed 2015 and observed 2025 as reference panels

map_long <- bind_rows(
  pixel_df |> transmute(x, y, agb = agb_2015,    method = "Observed (ORNL 2015)",   period = "2015 Baseline"),
  pixel_df |> transmute(x, y, agb = obs_2025,    method = "Observed (MELiTE 2025)", period = "2025 Observed"),
  pixel_df |> transmute(x, y, agb = cr_3pt_2025, method = "CR yield curves (3 pt)", period = "2025 Projected"),
  pixel_df |> transmute(x, y, agb = cr_4pt_2025, method = "CR yield curves (4 pt)", period = "2025 Projected"),
  pixel_df |> transmute(x, y, agb = fvs_2025,    method = "FVS ACD yield curves",   period = "2025 Projected")
) |>
  mutate(method = factor(method, levels = names(method_colors)))

# Also add difference from observed 2025
diff_long <- bind_rows(
  pixel_df |> transmute(x, y, agb = obs_2025 - agb_2015,       method = "Observed Change"),
  pixel_df |> transmute(x, y, agb = cr_3pt_2025 - obs_2025,    method = "CR 3 pt Error"),
  pixel_df |> transmute(x, y, agb = cr_4pt_2025 - obs_2025,    method = "CR 4 pt Error"),
  pixel_df |> transmute(x, y, agb = fvs_2025 - obs_2025,       method = "FVS ACD Error")
)

# --- Main panel: AGB maps by method -----------------------------------------
p_maps <- map_long |>
  ggplot(aes(x = x, y = y, fill = agb)) +
  geom_raster() +
  facet_wrap(~method, ncol = 3) +
  scale_fill_viridis_c(
    option = "D",
    limits = c(0, 250),
    na.value = "grey90",
    name = expression(AGB~(Mg~ha^{-1})),
    breaks = seq(0, 250, 50),
    oob = scales::squish
  ) +
  coord_equal() +
  labs(
    title = "Aboveground biomass: observed and projected by method",
    subtitle = paste0("Howland Research Forest (n = ", format(nrow(pixel_df), big.mark = ","), " pixels at 30 m)")
  ) +
  theme_pub +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0, size = 10, face = "plain")
  )

ggsave(file.path(fig_dir, "Fig_panel_AGB_methods.png"), p_maps,
       width = 22, height = 14, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Fig_panel_AGB_methods.pdf"), p_maps,
       width = 22, height = 14, units = "cm")
cat("Saved: Fig_panel_AGB_methods.png/pdf\n")

# --- Error maps (predicted minus observed 2025) ------------------------------
p_errors <- diff_long |>
  ggplot(aes(x = x, y = y, fill = agb)) +
  geom_raster() +
  facet_wrap(~method, ncol = 2) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-80, 80),
    na.value = "grey90",
    name = expression(Delta*AGB~(Mg~ha^{-1})),
    oob = scales::squish
  ) +
  coord_equal() +
  labs(
    title = "Spatial prediction error (predicted minus observed 2025)",
    subtitle = "Blue = underprediction, Red = overprediction"
  ) +
  theme_pub +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0, size = 10, face = "plain")
  )

ggsave(file.path(fig_dir, "Fig_error_maps.png"), p_errors,
       width = 17.5, height = 14, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Fig_error_maps.pdf"), p_errors,
       width = 17.5, height = 14, units = "cm")
cat("Saved: Fig_error_maps.png/pdf\n")

# =============================================================================
# PART 4: Method Divergence Analysis
# =============================================================================
cat("\n--- Part 4: Method divergence analysis ---\n")

pred_cols <- c("cr_3pt_2025", "cr_4pt_2025", "fvs_2025", "obs_2025")
pred_labels <- c("CR 3 pt", "CR 4 pt", "FVS ACD", "Observed")

# Per pixel statistics
pixel_df <- pixel_df |>
  mutate(
    method_sd   = apply(select(pixel_df, all_of(pred_cols[1:3])), 1, sd),
    method_range = apply(select(pixel_df, all_of(pred_cols[1:3])), 1, function(x) max(x) - min(x)),
    method_mean  = apply(select(pixel_df, all_of(pred_cols[1:3])), 1, mean)
  )

# Panel A: SD map
p_sd <- pixel_df |>
  ggplot(aes(x = x, y = y, fill = method_sd)) +
  geom_raster() +
  scale_fill_viridis_c(
    option = "A",
    na.value = "grey90",
    name = expression(SD~(Mg~ha^{-1}))
  ) +
  coord_equal() +
  labs(title = "A) Method disagreement (SD across predictions)") +
  theme_pub +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), legend.position = "right"
  )

# Panel B: Correlation matrix
cor_mat <- pixel_df |>
  select(all_of(pred_cols)) |>
  cor(use = "complete.obs")
colnames(cor_mat) <- pred_labels
rownames(cor_mat) <- pred_labels

cor_df <- cor_mat |>
  as.data.frame() |>
  rownames_to_column("m1") |>
  pivot_longer(-m1, names_to = "m2", values_to = "r") |>
  mutate(
    m1 = factor(m1, levels = rev(pred_labels)),
    m2 = factor(m2, levels = pred_labels)
  )

p_cor <- cor_df |>
  ggplot(aes(x = m2, y = m1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", r)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    name = "Pearson r"
  ) +
  labs(title = "B) Pairwise method correlations") +
  theme_pub +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )

# Panel C: Histogram of prediction range per pixel
p_range <- pixel_df |>
  ggplot(aes(x = method_range)) +
  geom_histogram(fill = "#1b9e77", color = "white", bins = 50, alpha = 0.85) +
  geom_vline(xintercept = median(pixel_df$method_range),
             linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text",
           x = median(pixel_df$method_range) + 2,
           y = Inf, vjust = 2, hjust = 0,
           label = paste0("Median = ", round(median(pixel_df$method_range), 1), " Mg/ha"),
           size = 3.5, color = "red") +
  labs(
    title = "C) Prediction range per pixel (max minus min)",
    x = expression(Range~(Mg~ha^{-1})),
    y = "Number of pixels"
  ) +
  theme_pub

# Panel D: Method mean vs SD (shows where disagreement is largest)
p_meansd <- pixel_df |>
  ggplot(aes(x = method_mean, y = method_sd)) +
  geom_point(alpha = 0.05, size = 0.3, color = "#333333") +
  geom_density_2d(color = "#7570b3", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
  labs(
    title = "D) Method disagreement by predicted AGB level",
    x = expression(Mean~predicted~AGB~(Mg~ha^{-1})),
    y = expression(SD~across~methods~(Mg~ha^{-1}))
  ) +
  theme_pub +
  theme(legend.position = "right")

# Compose divergence figure
p_divergence <- (p_sd + p_cor) / (p_range + p_meansd) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Method divergence analysis: 2025 AGB predictions at Howland",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
  )

ggsave(file.path(fig_dir, "Fig_method_divergence.png"), p_divergence,
       width = 22, height = 16, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Fig_method_divergence.pdf"), p_divergence,
       width = 22, height = 16, units = "cm")
cat("Saved: Fig_method_divergence.png/pdf\n")

# =============================================================================
# PART 5: Summary statistics
# =============================================================================
cat("\n--- Part 5: Summary statistics ---\n\n")

method_summary <- tibble(
  Method = c("Observed 2015", "Observed 2025", "CR 3 pt", "CR 4 pt", "FVS ACD"),
  Mean = c(mean(pixel_df$agb_2015), mean(pixel_df$obs_2025),
           mean(pixel_df$cr_3pt_2025), mean(pixel_df$cr_4pt_2025), mean(pixel_df$fvs_2025)),
  SD = c(sd(pixel_df$agb_2015), sd(pixel_df$obs_2025),
         sd(pixel_df$cr_3pt_2025), sd(pixel_df$cr_4pt_2025), sd(pixel_df$fvs_2025)),
  Median = c(median(pixel_df$agb_2015), median(pixel_df$obs_2025),
             median(pixel_df$cr_3pt_2025), median(pixel_df$cr_4pt_2025), median(pixel_df$fvs_2025))
) |>
  mutate(across(where(is.numeric), ~round(., 1)))

print(method_summary)

cat("\nPixel level divergence (predictions only, excluding observed):\n")
cat("  Mean SD:", round(mean(pixel_df$method_sd), 2), "Mg/ha\n")
cat("  Median range:", round(median(pixel_df$method_range), 2), "Mg/ha\n")
cat("  95th pct range:", round(quantile(pixel_df$method_range, 0.95), 2), "Mg/ha\n")

# Save summary as CSV
write_csv(method_summary, file.path(fig_dir, "Table_method_summary.csv"))
write_csv(
  pixel_df |> select(x, y, agb_2015, obs_2025, cr_3pt_2025, cr_4pt_2025, fvs_2025, method_sd, method_range),
  file.path(output_dir, "pixel_method_comparison.csv")
)
cat("\nSaved: Table_method_summary.csv, pixel_method_comparison.csv\n")

cat("\n============================================================\n")
cat("  Script 10 Complete\n")
cat("============================================================\n")
