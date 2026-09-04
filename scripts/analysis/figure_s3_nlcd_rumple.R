# =============================================================================
# figure_s3_nlcd_rumple.R
#
# Purpose: Plot rumple by NLCD class.
#
# Inputs:
#   outputs/insitu_covre_and_lidar_all_sites.csv
#
# Outputs:
#   outputs/nlcd_invasion_structure.csv
#   outputs/nlcd_scatter_rumple_rann.png
# =============================================================================


# =============================================================================
# Libraries
# =============================================================================
library(tidyverse)
library(patchwork)
library(cowplot)


# =============================================================================
# Structural data
# =============================================================================
### Read in full structural data
str_full <- read.csv("outputs/insitu_covre_and_lidar_all_sites.csv", header = T) |>
  mutate(rer = nspp_exotic / nspp_total)   # RRNN (same as data_preperation.R)

### List metrics
structural_metrics <- c("mean.max.canopy.ht.aop", "max.canopy.ht.aop", "rumple.aop",
                        "deepgap.fraction.aop", "cover.fraction.aop", "top.rugosity.aop",
                        "vert.sd.aop", "vertCV.aop", "sd.sd.aop", "entropy.aop",
                        "GFP.AOP.aop", "VAI.AOP.aop", "VCI.AOP.aop", "q25.aop", "q50.aop")

invasion_metrics <- c("shannon_exotic", "nspp_exotic", "rel_cover_exotic", "cover_exotic",
                      "shannon_native", "nspp_native", "rel_cover_native", "cover_native",
                      "nspp_notexotic", "rel_cover_notexotic", "shannon_total", "nspp_total",
                      "rer")


# =============================================================================
# USER INPUT
# =============================================================================
# Structural attribute to show (x axis)
all_str <- "rumple.aop"

# Inv attirbute to show (y axis)
# RANN; "rer" for RRNN
all_inv <- "rel_cover_exotic"

# color by site or nlcd class
# "site" or "nlcd"
all_colour <- "nlcd"               

# filter sites
# empty -> all sites
all_sites <- c()                  

# filter nlcd classes
# empty -> all; same names as nlcd_show
all_nlcd <- c()                  

# Show the overall trendline
# black overall lm
all_trend_overall <- TRUE          

# Trend per nlcd class
# one lm per site or NLCD class
all_trend_group   <- TRUE          


# =============================================================================
# ~ C O L O R S ~
# =============================================================================
### Stable site colors
sites_all <- sort(unique(str_full$site))
str_full <- str_full |>
  mutate(site = factor(site, levels = sites_all))

site_colours <- setNames(scales::hue_pal()(length(sites_all)), sites_all)

# NLCD class colours from MLRC
nlcd_colours <- c(
  "Evergreen forest" = "#1C5F2C",
  "Mixed forest" = "#B5C58E",
  "Deciduous forest" = "#68AB5F",
  "Woody wetlands" = "#6BAED6",
  "Shrub/scrub" = "#CCB879",
  "Dwarf scrub" = "#D1D182",
  "Grassland/herbaceous" = "#C4B84A",
  "Pasture/hay" = "#DCD939"
)


# =============================================================================
# Structural metric + inv metrics table by NLCD class
# =============================================================================
nlcd_mean_col <- paste0("mean_", all_str)
nlcd_sd_col   <- paste0("sd_", all_str)

nlcd_base <- str_full |>
  mutate(
    nlcd_class = str_trim(str_extract(NLCD_plot_des, "^[^:]+")),
    is_invaded = as.integer(invaded == "invaded")
  ) |>
  filter(!is.na(NLCD_code_plot), !is.na(nlcd_class), nlcd_class != "")

extra_nlcd <- setdiff(unique(nlcd_base$nlcd_class), names(nlcd_colours))
if (length(extra_nlcd) > 0) {
  nlcd_colours <- c(nlcd_colours, setNames(rep("grey70", length(extra_nlcd)), extra_nlcd))
}
nlcd_base <- nlcd_base |>
  mutate(nlcd_class = factor(nlcd_class, levels = names(nlcd_colours)))

nlcd_labels <- nlcd_base |>
  group_by(NLCD_code_plot) |>
  summarise(nlcd_class = first(nlcd_class), .groups = "drop")

plot_n <- nlcd_base |>
  group_by(NLCD_code_plot) |>
  summarise(
    n_plots = n_distinct(plotID),
    sd_structural = sd(.data[[all_str]], na.rm = TRUE),
    .groups = "drop"
  )

nlcd_by_year <- nlcd_base |>
  group_by(NLCD_code_plot, year) |>
  summarise(
    n = n(),
    nspp = mean(nspp_exotic, na.rm = TRUE),
    i = mean(is_invaded, na.rm = TRUE),
    ir = mean(rer, na.rm = TRUE),
    ii = mean(rel_cover_exotic, na.rm = TRUE),
    str_m = mean(.data[[all_str]], na.rm = TRUE),
    .groups = "drop"
  )

nlcd_table <- nlcd_by_year |>
  group_by(NLCD_code_plot) |>
  summarise(
    mean_exotic_spp = mean(nspp, na.rm = TRUE),
    percent_invaded_plots = mean(i, na.rm = TRUE),
    mean_RRNN = mean(ir, na.rm = TRUE),
    mean_RANN = mean(ii, na.rm = TRUE),
    mean_structural = mean(str_m, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(nlcd_labels, by = "NLCD_code_plot") |>
  left_join(plot_n, by = "NLCD_code_plot") |>
  drop_na(nlcd_class) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  select(
    `NLCD class` = nlcd_class,
    `Mean number of non-native species` = mean_exotic_spp,
    `Percent of invaded plots` = percent_invaded_plots,
    `Mean RRNN` = mean_RRNN,
    `Mean RANN` = mean_RANN,
    n_plots,
    mean_structural,
    sd_structural
  ) |>
  rename(
    `Number of plots` = n_plots,
    !!nlcd_mean_col := mean_structural,
    !!nlcd_sd_col := sd_structural
  ) |>
  arrange(`NLCD class`)

print(nlcd_table)
write_csv(nlcd_table, "outputs/nlcd_invasion_structure.csv")


# =============================================================================
# Two-panel scatter - Rumple vs. RANN
#   Left: all plots (all NLCD classes)
#   Right: forested plots only (Deciduous, Evergreen, Mixed, Woody wetlands)
# =============================================================================
## Define forest classes
forested_classes <- c("Deciduous forest", "Evergreen forest",
                      "Mixed forest", "Woody wetlands")

rumple_rann_base <- str_full |>
  mutate(
    nlcd_class = factor(
      str_trim(str_extract(NLCD_plot_des, "^[^:]+")),
      levels = names(nlcd_colours)
    ),
    RANN = rel_cover_exotic,
    Rumple = rumple.aop
  ) |>
  # Remove NAs for these plots
  drop_na(Rumple, RANN, nlcd_class)

# Compute shared y-axis limits across BOTH panels so axes are fixed
rann_max <- max(rumple_rann_base$RANN, na.rm = TRUE)
shared_ylim <- c(0, rann_max * 1.03)

# shared panel builder - y limits fixed; no legend
make_rumple_panel <- function(df, panel_title, y_lab = "RANN") {
  ggplot(df, aes(x = Rumple, y = RANN)) +
    # Points
    geom_point(
      aes(fill = nlcd_class),
      shape = 21, colour = "black", alpha = 0.5, stroke = 0.4, size = 2.5
    ) +
    # Per-NLCD trend lines
    geom_smooth(
      aes(colour = nlcd_class, fill = nlcd_class),
      method = "lm", se = TRUE, linewidth = 1, alpha = 0.20,
      show.legend = FALSE
    ) +
    # # Overall black trend line
    # geom_smooth(
    #   method = "lm", se = TRUE, colour = "black",
    #   linewidth = 1.5, alpha = 0.3, show.legend = FALSE
    # ) +
    scale_colour_manual(values = nlcd_colours, drop = TRUE, guide = "none") +
    scale_fill_manual(
      values = nlcd_colours, name = "NLCD class", drop = TRUE,
      guide  = guide_legend(override.aes = list(
        shape = 21, colour = "black", size = 3.5, alpha = 1, stroke = 0.5
      ))
    ) +
    coord_cartesian(ylim = shared_ylim) +
    labs(title = panel_title, x = "Rumple", y = y_lab) +
    theme_bw(base_size = 15)
}

# Build panels
p_rumple_all <- make_rumple_panel(rumple_rann_base, "A. All NLCD classes")
p_rumple_forest <- make_rumple_panel(
  filter(rumple_rann_base, nlcd_class %in% forested_classes),
  panel_title = "B. Forested NLCD classes only",
  y_lab       = NULL
)

# Extract legend from the left panel
legend_grob <- cowplot::get_legend(p_rumple_all)

# Strip legends from both panels, assemble with legend on the right
p_left  <- p_rumple_all    + theme(legend.position = "none")
p_right <- p_rumple_forest + theme(legend.position = "none")

p_panels <- cowplot::plot_grid(p_left, p_right,
                               nrow = 1, align = "h", axis = "tb")

# Title
title_row <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Rumple vs. relative abundance of non-natives (RANN)",
    fontface = "bold", x = 0.01, hjust = 0, size = 13
  )

# Final plot
p_rumple_rann_final <- cowplot::plot_grid(
  title_row,
  cowplot::plot_grid(p_panels, legend_grob, nrow = 1, rel_widths = c(1, 0.18)),
  ncol = 1, rel_heights = c(0.07, 1)
)

# Save
ggsave("outputs/nlcd_scatter_rumple_rann.png", p_rumple_rann_final,
       width = 14, height = 6, dpi = 150)
