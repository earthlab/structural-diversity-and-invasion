# adam's analysis - data prep
library(tidyverse)
library(ranger)
library(randomForest)
library(pdp)
library(topomicro)
library(sf)
library(terra)
d_old <- read_csv('data_quality_control/forest_div_quality_control_v1.csv') |>
  dplyr::filter(!site %in% c('SJER', 'ONAQ', 'BLAN', "STEI"))|>
  dplyr::select(plotID,starts_with("X_TSC"), -ends_with('median.1') , year, site) |>
  pivot_longer(-c(plotID, year, site)) |>
  mutate(tsc_year = str_extract(name, '\\d{4}')) |>
  dplyr::filter(year == tsc_year) |>
  dplyr::mutate(value = ifelse(value == 0, 99999, value)) |>
  dplyr::select(plotID, year, days_since_change = value, site) |>
  unique();d_old

dd <- read_csv("output/insitu_and_lidar_all_in_one_final_5_31_2023.csv") %>% 
  dplyr::mutate(folded_aspect = folded_aspect(aspect),
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect,
                slope_fa = slope*folded_aspect,
                rer = nspp_exotic/nspp_total,
                i_cat = ifelse(nspp_exotic >0, 1, 0)) |>
  dplyr::filter(!site %in% c('SJER', 'ONAQ', 'BLAN', "STEI")) |>
  dplyr::rename(rumple = "rumple.aop", deepgap.fraction = "deepgap.fraction.aop",
                mean.max.canopy.ht = mean.max.canopy.ht.aop,
                max.canopy.ht = max.canopy.ht.aop, minelev = minElev,
                cover.fraction = "cover.fraction.aop", top.rugosity = "top.rugosity.aop",
                vert.sd = "vert.sd.aop",  vertcv = "vertCV.aop", sd.sd = "sd.sd.aop",
                entropy = "entropy.aop", gfp = "GFP.AOP.aop", vai = "VAI.AOP.aop",
                vci = "VCI.AOP.aop", q25= "q25.aop", q50 = "q50.aop") 

dd|>
  dplyr::left_join(d_old)
dd_ak <- filter(dd, site %in% c("HEAL", "DEJU", "BONA"))
dd_cus <- filter(dd, !site %in% c("HEAL", "DEJU", "BONA")) |>
  dplyr::left_join(d_old)


plotids_ak <- dd_ak$plotID %>% unique()
plotids_cus <- dd_cus$plotID %>% unique()



# maybe separate slope and folded aspect
# maybe get rid of .aop
centroids_cus <- read_csv("data/All_NEON_TOS_Plot_Centroids_V8.csv") %>% 
  filter(plotID %in% plotids_cus, subtype == 'basePlot') %>%
  dplyr::select(latitude, longitude, plotID, domainID) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

centroids_ak <- read_csv("data/All_NEON_TOS_Plot_Centroids_V8.csv") %>% 
  filter(plotID %in% plotids_ak, subtype == 'basePlot') %>%
  dplyr::select(latitude, longitude, plotID, domainID) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
# 30 year normals of ppt, mean temp, and VPD min and max (resolution 800m) from https://prism.oregonstate.edu/normals/
bils <-list.files("data", recursive = T,
                        pattern = "annual_bil.bil$",
                        full.names = TRUE) %>%
  terra::rast() 

bils_ak <- list.files("data", recursive = T,
                    pattern = "14.txt$",
                    full.names = TRUE) %>%
  terra::rast() 


terra::crs(bils_ak) <- 'epsg:6397'

norms <- centroids %>%
  mutate(terra::extract(bils, vect(.), ID=F)) %>%
  rename(map = 4, mat = 5, tmax=6, tmin=7) %>%
  st_set_geometry(NULL) %>%
  unique()

caa <- centroids_ak %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2],
         x = ifelse(x < 0, 360 +x)) |>
  st_set_geometry(NULL) |>
  st_as_sf(coords =c('x', 'y'))

plot(bils_ak[[1]]);plot(caa[0], add=T)

norms_ak <- caa %>%
  mutate(terra::extract(bils_ak, vect(.), ID=F)) %>%
  rename(map = 4, mat = 5, tmax=6, tmin=7) %>%
  st_set_geometry(NULL) %>%
  unique()
# nrow(norms)
# nrow(dd)
# norms[15,]
# dd[81,]

d_cus <- dd_cus %>%
  left_join(norms)

d_ak <- dd_ak |>
  left_join(norms_ak)

bind_rows(d_cus, d_ak) |> #summary()
  write_csv("data/data_with_climate_norms.csv")

 