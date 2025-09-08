# alternate download of cover data

# neondiveRsity ================
# devtools::install_github("admahood/neonPlantEcology")
library(neonPlantEcology)


root <- paste0(getwd(), "/")

outfile <- paste0(root, "data/input_cover_data.Rda")

data <- read.csv(paste0(root, "data/NEON_sites_dates_for_cover.csv"))

if(!file.exits(outfile)){
  # takes a while
  input_data <- lapply(data$siteID,neonPlantEcology::download_plant_div)
}else{
  load(outfile)
}

save(input_data, file = outfile)
