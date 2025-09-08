# alternate plant cover metrics extraction
library(tidyverse)
# devtools::install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
input_data_file <- list.files("data", pattern = "input_cover_data.Rda", full.names = TRUE)
load(input_data_file)

outfile_plot <- paste0("data/plot_site_summaries.csv")
previous_output <- read_csv("output/structural_metrics_by_plot_old.csv")

monthyear <- previous_output %>%
  group_by(siteID) %>%
  reframe(monthyear = unique(monthyear)) %>%
  ungroup() %>%
  mutate(year = str_sub(monthyear, 1,4))

if(!file.exists(outfile_plot)){
  plot_level <- lapply(input_data,
                       function(x)neonPlantEcology::get_diversity_info(x, scale = "plot",betadiversity = TRUE))%>%
    bind_rows()
  
  site_level <- lapply(input_data,
                       function(x)neonPlantEcology::get_diversity_info(x, scale = "site",betadiversity = TRUE)) %>%
    bind_rows()
  site_plot <- bind_rows(plot_level, site_level)%>%
    right_join(monthyear, by = c("site" = "siteID", "year"))
  
  write_csv(site_plot, file = outfile_plot)
}else{site_plot <- read_csv(outfile_plot)}



sites_n_dates <- read_csv(paste0("data/NEON_sites_dates_for_cover.csv")) %>%
  mutate_at(2:4, function(x) as.numeric(str_sub(x,1,4))) 

res <- list()
for(i in 1:length(unique(sites_n_dates$siteID))){
  ss <-  unique(sites_n_dates$siteID)[i]
  yearz <- sites_n_dates %>%
    filter(siteID == ss) %>%
    dplyr::select(-siteID) %>%
    as.matrix() %>%
    as.numeric()
  yearz <- yearz[!is.na(yearz)]
  
  res[[i]] <- filter(site_plot, site == ss, year %in% yearz)
}


lidar_data <- read_csv("output/lidar_structural_metrics.csv")  %>%
  dplyr::mutate(year = str_sub(monthyear,1,4) %>% as.numeric(),
                year_plot = str_c(year, plotID))%>%
  dplyr::select(-monthyear, -sitemonthyear, -siteID, -year)
plot_yrz <- unique(lidar_data$year_plot)

site_plot_filtered <- bind_rows(res) %>%
  dplyr::mutate(year_plot = str_c(year, plotID)) %>%
  filter(year_plot %in% plot_yrz) %>%
  dplyr::select(-subplotID, -scale)

spl_out <- site_plot_filtered  %>%
  left_join(lidar_data, by = c("year_plot", "plotID"))

write_csv(spl_out,"output/plants_lidar_locations.csv")
# veg_types <- read_csv('data/field-sites.csv')
