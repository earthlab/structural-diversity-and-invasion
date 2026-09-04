# =============================================================================
# boxS2_dela_caseStudy.R
#
# Purpose: DELA case study to classify invasive species by plant functional type
#   (growth habit) and compare forest structural metrics between plots dominated
#   by different invasive functional types. Also looks at native richness and
#   NLCD habitat class.
#
# Inputs:
#   data/input_cover_data.Rda                        NEON 1m2 cover, DELA only
#   data/USDA_PLANTS_growthhabit.csv                 USDA growth-habit lookup
#   outputs/insitu_covre_and_lidar_all_sites.csv     structural + invasion metrics
#
# Outputs:
#   outputs/dela_functype_structural.png    structural metrics by functional type
#   outputs/dela_functype_community.png     native richness, RANN, NLCD by type
#   outputs/dela_functype_scatter.png       figure in box S2
# =============================================================================


# =============================================================================
# Libraries
# =============================================================================
library(tidyverse)
library(patchwork)
library(neonPlantEcology)


# =============================================================================
# ~ C O L O R S ~
# =============================================================================
# Note: we're removing grass/forbs because there are only two plots with those
# PFTs as the dominant invasive
functype_colours <- c(
  "Tree" = "#9C4E86",
  "Shrub" = "#167984",
  "Vine" = "#5AAD6D"
)


# =============================================================================
# USDA PLANTS Growth Habit Data
# =============================================================================
usda_raw <- read_csv(
  "data/USDA_PLANTS_growthhabit.csv",
  #locale = locale(encoding = "latin1"),
  show_col_types = FALSE
)

# Manual overrides for NEON taxon codes not found in USDA (code mismatches):
# MIVI  = Microstegium vimineum (Japanese stiltgrass) --> Forb/herb
# DUCHE = Duchesnea indica (mock strawberry) --> Forb/herb
# DUIN  = likely Duchesnea indica alternate code --> Forb/herb
# ASTER = genus-level Aster sp. --> Forb/herb NEON has this listed as invasive??
manual_habits <- tribble(
  ~taxonID,  ~growth_habit,
  "MIVI", "Forb/herb",
  "DUCHE", "Forb/herb",
  "DUIN", "Forb/herb",
  "ASTER", "Forb/herb",
  # LISI is listed as a tree but USDA PLANTS online has it as a shrub/tree
  # and it functions more like a tall shrub in DELA
  "LISI", "Shrub"
)

usda <- usda_raw |>
  rename(taxonID = `Accepted Symbol`) |>
  # LYJA has two entries: Forb/herb and Vine
  # Restricting to vine for the invasion type at DELA
  filter(!(taxonID == "LYJA" & growth_habit == "Forb/herb")) |>
  # Remove any taxon codes that will be overridden by manual_habits
  filter(!taxonID %in% manual_habits$taxonID) |>
  bind_rows(manual_habits)


# =============================================================================
# NEON 1m2 data PPPC
# =============================================================================
# Add you NEON API token here
### DO NOT COMMIT WITH THIS STILL IN THE CODE
token <- ""

root <- paste0(getwd(), "/")

outfile <- paste0(root, "data/input_cover_data.Rda")

data <- read.csv(paste0(root, "data/NEON_sites_dates_for_cover.csv")) |>
  filter(siteID == "DELA")

if(!file.exists(outfile)){
  # takes a while
  input_data <- lapply(data$siteID, function(site) {
    neonPlantEcology::npe_download(sites = site, token = token)
  })
}else{
  load(outfile)
}

save(input_data, file = outfile)

# Years available in the structural data for DELA
structural_years <- c(2015L, 2016L, 2017L)

# Join to growth habit data
cover_exotic <- input_data[[1]][["div_1m2Data"]] |>
  filter(
    divDataType == "plantSpecies",
    nativeStatusCode == "I",
    !is.na(taxonID), taxonID != ""
  ) |>
  mutate(
    year = as.integer(format(as.Date(endDate), "%Y")),
    percentCover = as.numeric(percentCover)
  ) |>
  filter(year %in% structural_years, !is.na(percentCover)) |>
  left_join(usda, by = "taxonID")


### Determine dominant invader

# Sum cover by growth habit within each plot/year
dom_type_yr <- cover_exotic |>
  filter(!is.na(growth_habit)) |>
  group_by(plotID, year, growth_habit) |>
  summarise(total_cover = sum(percentCover, na.rm = TRUE), .groups = "drop") |>
  group_by(plotID, year) |>
  slice_max(total_cover, n = 1, with_ties = FALSE) |>
  ungroup() |>
  rename(dom_functype = growth_habit)

# Pick the dominant invader growth habit for each plot across the years sampled
dom_type_plot <- dom_type_yr |>
  group_by(plotID) |>
  summarise(
    dom_functype = names(sort(table(dom_functype), decreasing = TRUE))[1],
    .groups = "drop"
  )


# =============================================================================
# Structural metrics
# =============================================================================
structural_metrics <- c(
  "entropy.aop",
  "rumple.aop",
  "VAI.AOP.aop",
  "deepgap.fraction.aop",
  "cover.fraction.aop",
  "mean.max.canopy.ht.aop"
)

metric_labels <- c(
  entropy.aop = "Canopy Entropy",
  rumple.aop = "Rumple",
  VAI.AOP.aop = "Vegetation Area Index",
  deepgap.fraction.aop = "Deep Gap Fraction",
  cover.fraction.aop = "Cover Fraction",
  mean.max.canopy.ht.aop = "Mean Max Canopy Ht (m)"
)

# Filter to DELA
str_dela <- read.csv("outputs/insitu_covre_and_lidar_all_sites.csv") |>
  filter(site == "DELA") |>
  mutate(
    year = as.integer(year),
    nlcd_class = str_trim(str_extract(NLCD_plot_des, "^[^:]+"))
  )

# Average structural and invasion metrics across years per plot
str_plot <- str_dela |>
  group_by(plotID) |>
  summarise(
    across(all_of(structural_metrics), \(x) mean(x, na.rm = TRUE)),
    nspp_native = mean(nspp_native, na.rm = TRUE),
    nspp_exotic = mean(nspp_exotic, na.rm = TRUE),
    rel_cover_exotic = mean(rel_cover_exotic, na.rm = TRUE),
    nlcd_class = first(nlcd_class),
    .groups = "drop"
  )


# =============================================================================
# Final analysis df
# =============================================================================
analysis_df <- str_plot |>
  inner_join(dom_type_plot, by = "plotID") |>
  mutate(
    display_group = factor(
      dom_functype,
      levels = c("Tree", "Shrub", "Vine",
                 "Forb/herb", "Subshrub")
    )
  ) |>
  # Restrict to Woody Wetland plots dominated by the three comparable groups.
  # This removes 8 Vine-Evergreen plots (different habitat) and 2 Forb/herb
  # plots (n too small for comparison), leaving 18 plots.
  filter(
    grepl("wetland", nlcd_class, ignore.case = TRUE),
    display_group %in% c("Tree", "Shrub", "Vine")
  ) |>
  mutate(display_group = droplevels(display_group))

n_matched <- nrow(analysis_df)
n_total   <- nrow(str_plot)

### Top inv species per plot
top_sp <- cover_exotic |>
  filter(!is.na(growth_habit)) |>
  group_by(plotID, taxonID, growth_habit) |>
  summarise(total_cover = sum(percentCover, na.rm = TRUE), .groups = "drop") |>
  group_by(plotID) |>
  slice_max(total_cover, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(growth_habit, plotID)


# =============================================================================
# Figure 1: Structural characteristics across growth habits
# =============================================================================
present_types <- levels(droplevels(analysis_df$display_group))

str_long <- analysis_df |>
  pivot_longer(
    all_of(structural_metrics),
    names_to  = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric_label = metric_labels[metric],
    metric_label = factor(metric_label, levels = unname(metric_labels))
  )

# n labels for legend
n_labels <- analysis_df |>
  count(display_group) |>
  mutate(label = paste0("n=", n))

p_struct <- ggplot(str_long,
                   aes(x = display_group, y = value, fill = display_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5) +
  geom_jitter(
    shape = 21, width = 0.18, height = 0, alpha = 0.8,
    colour = "black", stroke = 0.4, size = 2.5
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = functype_colours,
    breaks = present_types,
    name = "Dominant invasive\nfunctional type"
  ) +
  labs(
    title = "DELA: forest structural metrics by dominant invasive plant functional type",
    subtitle = sprintf(
      "Each point = one plot (mean across survey years 2015-2017); n = %d plots total",
      n_matched
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  )

# ggsave("outputs/dela_functype_structural.png", p_struct,
#        width = 12, height = 8, dpi = 150)


# =============================================================================
# Figure 2: Native species and RANN per growth habit
# =============================================================================
# Panel A: native species richness
pA <- ggplot(analysis_df,
             aes(x = display_group, y = nspp_native, fill = display_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5) +
  geom_jitter(
    shape = 21, width = 0.18, height = 0.1, alpha = 0.8,
    colour = "black", stroke = 0.4, size = 2.5
  ) +
  geom_text(
    data = n_labels, aes(x = display_group, label = label),
    y = -Inf, vjust = -0.3, size = 3.5, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = functype_colours, breaks = present_types, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +
  labs(
    title = "A. Native species richness",
    x = NULL,
    y = "No. native species"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel B: relative exotic cover (RANN)
pB <- ggplot(analysis_df,
             aes(x = display_group, y = rel_cover_exotic, fill = display_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5) +
  geom_jitter(
    shape = 21, width = 0.18, height = 0, alpha = 0.8,
    colour = "black", stroke = 0.4, size = 2.5
  ) +
  scale_fill_manual(values = functype_colours, breaks = present_types, guide = "none") +
  labs(
    title = "B. Relative non-native cover (RANN)",
    x = NULL,
    y = "RANN"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_community <- (pA | pB) +
  plot_annotation(
    title = "DELA: native community context by dominant invasive functional type",
    subtitle = sprintf("n = %d Woody Wetland plots (mean across 2015-2017)", n_matched)
  )

# ggsave("outputs/dela_functype_community.png", p_community,
#        width = 9, height = 5, dpi = 150)


# =============================================================================
# Figure 3: final plot comparison of RANN and structural attributes per GH
# =============================================================================
# Four focal metrics that capture distinct structural axes
scatter_metrics <- c(
  #"entropy.aop",
  "rumple.aop",
  "deepgap.fraction.aop"
  #"mean.max.canopy.ht.aop"
)
scatter_labels <- c(
  #entropy.aop = "Canopy Entropy",
  rumple.aop = "Rumple",
  deepgap.fraction.aop = "Deep Gap Fraction"
  #mean.max.canopy.ht.aop = "Mean Max Canopy Height"
)

scatter_long <- analysis_df |>
  pivot_longer(
    all_of(scatter_metrics),
    names_to  = "metric",
    values_to = "str_val"
  ) |>
  mutate(
    rann = rel_cover_exotic,
    metric_label = factor(scatter_labels[metric],
                          levels = unname(scatter_labels))
  )

# Final plot + trend lines
p_scatter <- ggplot(scatter_long,
                    aes(x = str_val, y = rann)) +
  geom_smooth(
    aes(colour = display_group, fill = after_scale(colour)),
    method = "lm", se = TRUE,
    linewidth = 1, alpha = 0.15,
    show.legend = FALSE
  ) +
  # Points: border = functional type, fill = native richness gradient
  geom_point(
    aes(fill = nspp_native, colour = display_group),
    shape = 21, size = 3.5, stroke = 2, alpha = 0.9
  ) +
  facet_wrap(~ metric_label, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  # Border colour and SE ribbon fill both driven by functype_colours
  scale_colour_manual(
    values = functype_colours,
    breaks = levels(analysis_df$display_group),
    name = "Dominant invasive\nfunctional type"
  ) +
  # Native richness scale
  scale_fill_gradient(
    low = "white",
    high = "black",
    name = "Native species\nrichness"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(shape = 21, fill = "grey80",
                          size = 4, stroke = 1.2, alpha = 1),
      order = 1
    ),
    fill = guide_colorbar(order = 2, barheight = 6)
  ) +
  labs(
    title = "DELA: structural metrics vs. relative abundance of non-natives (RANN)",
    x = NULL,
    y = "RANN"
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(size = 15),
    legend.position = "right"
  )

ggsave("outputs/dela_functype_scatter.png", p_scatter,
       width = 11, height = 5, dpi = 150)