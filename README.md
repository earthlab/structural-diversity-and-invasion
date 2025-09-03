<!-- FILLME:START -->
# Structural Diversity and Invasion
This repository downloads, cleans, and analyzes National Ecological Observatory Network (NEON) lidar and plot level plant cover data to look at the relationship between forest structural diversity and resistance to plant invasion across 26 NEON forested sites.

## Summary workflow

### 1. Download NEON plot level data for all date site combinations  
a. Downloads data for all date site combinations NEON_sites_dates_for_cover.csv and All_NEON_TOS_Plot_Centroids_V8.csv from DDDD. The files are in compressed .rar format. 

b. Runs above files using all date site combinations to make a combined dataset. The resultant is a csv file that contains all plot level information including vegetation, topography, soil, and coordinates of plot centers. We use *helper.R* for coordinate conversion to match with the reference systems of other support files. 
