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

Downloads lidar data from [NEON](https://data.neonscience.org/data-products/DP1.30003.001) from NEON for all site–date combinations identified in Step 1. To minimize boundary effects, each point cloud was clipped to a 200 m × 200 m area (surrounding a 20 m × 20 m plot). This produced one LiDAR point cloud (in .laz format) for each plot within each site–date combination. Due to size constraints (~2.63 GB in total), the files are not included in this repository.

### 4. Structural diversity metric extraction 
Extracts a series of structural divesity metrics using pointcloud and point cloud derived canopy height meodels for each plot based on [LaRue et al., 2019](https://iopscience.iop.org/article/10.1088/1748-9326/ab49bb). Description of structural diversity metrics derived and references can be found [here](docs/Structural-diversity-metrics.pdf)

### 5. Analysis on forest structural diversity on invasion
Analyses the effects of structural diversity and site-specific conditions on forest invasion, using invasion metrics derived from plot-level plant cover data. The results demonstrate that forest structural diversity can serve as a useful indicator of invasion potential across 26 NEON forest sites distributed throughout the conterminous United States.
