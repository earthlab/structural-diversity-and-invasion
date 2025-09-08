<!-- FILLME:START -->
# Structural Diversity and Invasion
This repository downloads, cleans, and analyzes National Ecological Observatory Network (NEON) lidar and plot level plant cover data to look at the relationship between forest structural diversity and resistance to plant invasion across 26 NEON forested sites.

## Summary workflow

### 1. Download NEON plot level cover data for all date site combinations  
a. Downloads data for all date site combinations NEON_sites_dates_for_cover.csv and All_NEON_TOS_Plot_Centroids_V8.csv from DDDD. The files are in compressed .rar format. We selected only forest NEON sites for further processing in this study. 

### 2. Combine plot vegetation data with site specific data
a. Runs above files using all date site combinations to make a combined dataset. The resultant is a csv file that contains all plot level information including vegetation, topography, soil, and coordinates of plot centers. We use *helper.R* for coordinate conversion to match with the reference systems of other support files. 

b. Runs above files using all date site combinations to make a combined dataset. The resultant is a csv file that contains all plot level information including vegetation, topography, soil, and coordinates of plot centers. We use *helper.R* for coordinate conversion to match with the reference systems of other support files. 

### 3. Download NEON plot level lidar point cloud data

Downloads lidar data from NEON for all selected site-date combination from step 1. To avoid the boundary effects, we clipped lidar point clouds 200m x 200m (plot size is 20m x 20m). This results lidar point cloud per plot per each data-site combination in .laz format. Files have not uploaded here due to size restriction (~2.63 GB in total). 

### 4. Structural diversity metric extraction 

### 5. Analysis on forest structural diversity on invasion
