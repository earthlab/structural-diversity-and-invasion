# make a table of basic information about invasion rates
library(tidyverse)

# invasion table ===============================================================

dd <- read_csv("data/data_with_climate_norms.csv") %>%
  dplyr::mutate(folded_aspect = topomicro::folded_aspect(aspect),
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect,
                rer = (nspp_exotic/nspp_total)*100) 

glimpse(dd)
unique(dd$invaded)

invasion_by_site <- dd |>
  mutate(invaded = ifelse(invaded == 'invaded', 1, 0)) |>
  group_by(site, year) |>
  summarise(n = n(),
            nspp = mean(nspp_exotic),
            i = sum(invaded),
            ir = mean(rer),
            ii = mean(rel_cover_exotic)) |>
  ungroup()  |>
  mutate(i = i/n) |>
  group_by(site) |>
  summarise(mean_exotic_spp = mean(nspp),
            percent_invaded_plots = mean(i) * 100,
            mean_ir = mean(ir),
            mean_ii = mean(ii) * 100) |>
  ungroup()
write_csv(invasion_by_site |> mutate_if(is.numeric, round, 2),
          'output/table_X_invasion_by_site.csv')

invasion_by_site |>
  mutate(i_class = case_when(
    percent_invaded_plots == 100 ~ "100",
    percent_invaded_plots == 0 ~ "0",
    percent_invaded_plots > 0 & percent_invaded_plots < 50 ~ '0-50',
    percent_invaded_plots >= 50 & percent_invaded_plots < 100 ~ '50-100'
  )) |>
  pull(i_class) |>
  table()

# invasion by nlcd class =======================================================

ddd <- read_csv('data_quality_control/forest_div_quality_control_v1.csv') |>
  dplyr::select(NLCD_plot_des_main, NLCD_code_plot) |>
  unique()

dd <- read_csv("data/data_with_climate_norms.csv") |>
  left_join(ddd)


glimpse(dd)
unique(dd$invaded)

invasion_by_nlcd <- dd |>
  group_by(NLCD_plot_des_main, year) |>
  summarise(n = n(),
            nspp = mean(nspp_exotic),
            i = sum(i_cat),
            ir = mean(rer) *100,
            ii = mean(rel_cover_exotic)) |>
  ungroup()  |>
  mutate(i = i/n) |>
  group_by(NLCD_plot_des_main) |>
  summarise(mean_exotic_spp = mean(nspp),
            percent_invaded_plots = mean(i) * 100,
            relative_richness_non_native = mean(ir),#n = n(),
            relative_abundance_non_native = mean(ii) * 100) |>
  ungroup() |>
  na.omit()


plot_counts <- dd |>
  group_by(NLCD_plot_des_main) |>
  summarise(n = n())

write_csv(invasion_by_nlcd |> mutate_if(is.numeric, round, 2) |> left_join(plot_counts), 'output/invasion_by_nlcd.csv')

# table of structural metrics ==================================================

dd |>
  dplyr::select()

invasion_by_plot <- dd |>
  mutate(invaded = ifelse(invaded == 'invaded', 1, 0)) |>
  group_by(plotID, year) |>
  summarise(n = n(),
            i = sum(invaded),
            ir = mean(rer),
            ii = mean(rel_cover_exotic)) |>
  ungroup()  |>
  mutate(i = i/n) |>
  group_by(plotID) |>
  summarise(i = mean(i),
            ir = mean(ir),
            ii = mean(ii)) |>
  ungroup()


invasion_by_plot |>
  mutate(i_class = case_when(
    i == 1 ~ "100",
    i == 0 ~ "0",
    i > 0 & i < .5 ~ '0-50',
    i >= 0.5 & i < 1 ~ '50-100'
  )) |>
  pull(i_class) |>
  table()

# structure metrics by site ====================================================
dd |>
  dplyr::select(9:23#, -sd.sd.aop
                ) %>%
  pivot_longer(cols = names(.)) |>
  dplyr::group_by(name) |>
  summarise(mean = mean(value, na.rm=T),
            sd = sd(value, na.rm=T)) |>
  ungroup() |>
  mutate_if(is.numeric, signif, 2) |>
  mutate(name = str_remove_all(name, ".aop")|> str_remove_all('.AOP')) |>
  print() |>
  write_csv('output/structural_metric_table.csv')



