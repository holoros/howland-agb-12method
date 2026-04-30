# =============================================================================
# Title: Practitioner Decision Framework for AGB Projection Methods
# Author: A. Weiskittel
# Date: 2026-03-26
# Description: Creates a decision tree figure (Fig 9 or Fig 13) showing which
#              AGB projection method to use based on budget, spatial coverage
#              needs, projection horizon, and available data. Designed as a
#              practical takeaway for forest managers and carbon analysts.
# Dependencies: ggplot2, ggtext
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggtext)
library(patchwork)

# --- Paths -------------------------------------------------------------------
fig_dir <- "output/manuscript_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Theme -------------------------------------------------------------------
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(8, 8, 8, 8),
    plot.tag = element_text(size = 12, face = "bold")
  )

W_DOUBLE <- 17.5
H_STANDARD <- 12

cat("=== Script 09: Decision Framework Figure ===\n")

# =============================================================================
# PART 1: Method comparison summary table (data for the figure)
# =============================================================================

method_summary <- tribble(
  ~method, ~category, ~cost_per_acre, ~error_10yr_pct, ~spatial, ~horizon_max, ~data_needed, ~complexity,
  "Area-based CR",          "Remote sensing",  0.38, 2.4,  "Wall-to-wall", 50, "Multi-temporal ALS (3+ years)", "Low",
  "Area-based Schumacher",  "Remote sensing",  0.38, 4.6,  "Wall-to-wall", 50, "Multi-temporal ALS (3+ years)", "Low",
  "LiDAR + ACD curves",     "Hybrid",          0.43, 26.4, "Wall-to-wall", 50, "Single ALS + field plots",      "Medium",
  "NAIP + ACD curves",      "Hybrid",          0.18, 46.7, "Wall-to-wall", 50, "NAIP imagery + field plots",    "Medium",
  "FVS-NE (defaults)",       "Process model",   7.52, 24.2, "Plot-based",   100, "Field inventory",               "High",
  "FVS-NE (calibrated)",     "Process model",   7.52, 32.0, "Plot-based",   100, "Field inventory + ALS",         "High",
  "FVS-ACD (full)",          "Process model",   7.52, 44.9, "Plot-based",   100, "Field inventory",               "High"
)

# =============================================================================
# PART 2: Cost vs accuracy bubble chart with decision zones
# =============================================================================

# Decision zones based on practitioner use cases
zone_df <- tribble(
  ~xmin, ~xmax, ~ymin, ~ymax, ~label, ~fill,
  0.05,  0.30,  0,     15,    "Screening &\nrapid assessment", "#E8F5E9",
  0.20,  0.60,  0,     10,    "Operational\ncarbon monitoring", "#E3F2FD",
  5.00,  10.0,  0,     50,    "Research &\nprocess understanding", "#FFF3E0"
)

fig_bubble <- ggplot(method_summary,
                     aes(x = cost_per_acre, y = error_10yr_pct)) +
  # Decision zones
  geom_rect(data = zone_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, alpha = 0.25) +
  geom_text(data = zone_df,
            aes(x = sqrt(xmin * xmax), y = ymax * 0.85, label = label),
            inherit.aes = FALSE, size = 2.5, color = "grey30",
            fontface = "italic", lineheight = 0.9) +
  scale_fill_identity() +
  # Method points
  geom_point(aes(color = category, shape = spatial), size = 4, alpha = 0.85) +
  geom_text(aes(label = method), size = 2.3, vjust = -1.2,
            check_overlap = TRUE) +
  scale_x_log10(labels = scales::dollar_format(),
                breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10)) +
  scale_color_manual(values = c("Remote sensing" = "#009E73",
                                "Hybrid" = "#0072B2",
                                "Process model" = "#D55E00"),
                     name = "Approach") +
  scale_shape_manual(values = c("Wall-to-wall" = 16, "Plot-based" = 17),
                     name = "Coverage") +
  labs(x = "Cost per acre (USD, log scale)",
       y = "Absolute 10-year prediction error (%)",
       title = "Method selection guide: cost, accuracy, and spatial coverage") +
  theme_pub +
  theme(legend.box = "vertical",
        panel.grid.major = element_line(color = "grey90"))

# =============================================================================
# PART 3: Decision flowchart (text-based for the manuscript)
# =============================================================================

# Create a structured decision matrix
decision_matrix <- tribble(
  ~scenario,                              ~recommended_method,  ~rationale,
  "Budget < $0.25/ac, screening",          "NAIP + yield curves", "Lowest cost, adequate for initial assessment",
  "Budget $0.25-1/ac, carbon monitoring",  "Area-based CR",       "Best accuracy/cost ratio, wall-to-wall",
  "Budget > $5/ac, species-level detail",  "FVS-NE or FVS-ACD",  "Individual tree dynamics, species projections",
  "Projection < 20 years",                 "Area-based CR",       "Empirically constrained, low divergence",
  "Projection > 20 years",                 "FVS + periodic recal","Process models needed but must recalibrate",
  "Wall-to-wall map required",             "LiDAR-based methods", "Only remote sensing provides full coverage",
  "Disturbance expected",                  "FVS variants",        "Process models handle harvest/mortality events"
)

# Create a formatted decision table figure
decision_long <- decision_matrix |>
  mutate(row_id = row_number(),
         scenario = str_wrap(scenario, 30),
         rationale = str_wrap(rationale, 40))

fig_decision <- ggplot(decision_long, aes(y = fct_rev(factor(row_id)))) +
  geom_tile(aes(x = 1), fill = "grey95", height = 0.9, width = 1.8) +
  geom_tile(aes(x = 2.5), fill = "#E3F2FD", height = 0.9, width = 1.5) +
  geom_tile(aes(x = 4.2), fill = "grey95", height = 0.9, width = 2) +
  geom_text(aes(x = 1, label = scenario), size = 2.5, lineheight = 0.85) +
  geom_text(aes(x = 2.5, label = recommended_method),
            size = 2.5, fontface = "bold", color = "#0072B2") +
  geom_text(aes(x = 4.2, label = rationale), size = 2.2, lineheight = 0.85) +
  # Headers
  annotate("text", x = 1, y = nrow(decision_long) + 0.7,
           label = "Scenario", fontface = "bold", size = 3) +
  annotate("text", x = 2.5, y = nrow(decision_long) + 0.7,
           label = "Recommended", fontface = "bold", size = 3) +
  annotate("text", x = 4.2, y = nrow(decision_long) + 0.7,
           label = "Rationale", fontface = "bold", size = 3) +
  scale_x_continuous(limits = c(0, 5.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(15, 10, 10, 10),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) +
  labs(title = "Practitioner decision guide")

# Combined figure
fig_combined <- fig_bubble / fig_decision + plot_layout(heights = c(3, 2))

ggsave(file.path(fig_dir, "Fig13_decision_framework.pdf"), fig_combined,
       width = W_DOUBLE, height = 22, units = "cm")
ggsave(file.path(fig_dir, "Fig13_decision_framework.png"), fig_combined,
       width = W_DOUBLE, height = 22, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig13 (decision framework)\n")

# Also save just the bubble chart
ggsave(file.path(fig_dir, "Fig13a_cost_accuracy_zones.pdf"), fig_bubble,
       width = W_DOUBLE, height = H_STANDARD, units = "cm")
ggsave(file.path(fig_dir, "Fig13a_cost_accuracy_zones.png"), fig_bubble,
       width = W_DOUBLE, height = H_STANDARD, units = "cm", dpi = 300, bg = "white")
cat("  Saved Fig13a (bubble chart only)\n")

cat("\n=== Script 09 complete ===\n")
