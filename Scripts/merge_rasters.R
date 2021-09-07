library(raster)
library(sf)

load("/projects/nc57/Chapter_3/SpatialData/RL_shp_prepro.RData")
rm(RL_shp_Aus)

#Create an empty aus raster, same dimensions etc, fill with lowest number/zero
aus_raster <- raster(nrows = 3954,
                     ncols = 4859,
                     xmn = 113.15,
                     xmx = 153.6417,
                     ymn = -43.64167,
                     ymx = -10.69167,
                     crs = "+proj=longlat +datum=WGS84 +no_defs",
                     resolution = 0.008333333)
values(aus_raster) <- 0.00001
aus_raster <- mask(aus_raster, Aus_Coast)

species_path <- "/home/dcla0008/nc57_scratch/Zonation/Input_data/species"
species_area_path <- "/home/dcla0008/nc57_scratch/Zonation/Input_data/species_area"

models <- list.files("/home/dcla0008/nc57_scratch/Zonation/Input_data/maxent_output", full.names = T)

#merge extended maxent with empty raster
lapply(models[1:length(models)], function(i) {
  
  r <- raster(i)
  r_new <- merge(r, aus_raster)
  i <- gsub("/home/dcla0008/nc57_scratch/Zonation/Input_data/maxent_output/", "", i)
  i <- gsub("_pred.tif", "", i)
  writeRaster(r_new, filename = paste0(species_path,"/", i, "_pred.tif"), overwrite = T)
  writeRaster(r_new, filename = paste0(species_area_path,"/", i, "_pred.tif"), overwrite = T)
  
})
