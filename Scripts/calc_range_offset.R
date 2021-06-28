############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "calc_range_offset.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "calc_range_offset.R"),
              local_path = file.path("Scripts", "calc_range_offset.R"), overwrite = T)
##################################################################################################
                              # Creating bias file (with bossMaps) #
                                              # TEST #
#Could potentially alter Pin for each species based off information from cleaned occurrence records
#Could join the prop_within to shapefile

#selecting species
Ap <- RL_shp_Aus %>% filter(BINOMIAL == "Acanthiza pusilla")

#convert to sp
Ap_sp <- as(Ap, Class = "Spatial")

#elevation raster
#getData("alt", country = "AUS", mask = T, path = file.path("SpatialData", "Raster", "Elevation"))
Aus_elev <- raster(file.path("SpatialData", "Raster", "Elevation", "Aus_msk_alt.gri"))
#resolution is in degrees. 0.008333333 degrees = 0.5 arc min ~ 1km
Aus_elev <- mask(Aus_elev, Aus_Coast)
Aus_elev <- crop(Aus_elev, Aus_Coast)

#crop elevation to range map
Aus_elev_Ap <- crop(Aus_elev, extent(Ap_sp))

#Bioclim rasters
# Aus1 <- getData("worldclim", var = "bio", lon = 110, lat = -25, res = 0.5, path = file.path("SpatialData", "Raster", "Worldclim"))
# Aus2 <- getData("worldclim", var = "bio", lon = 130, lat = -25, res = 0.5, path = file.path("SpatialData", "Raster", "Worldclim"))
# Aus3 <- getData("worldclim", var = "bio", lon = 150, lat = -25, res = 0.5, path = file.path("SpatialData", "Raster", "Worldclim"))
# Aus4 <- getData("worldclim", var = "bio", lon = 110, lat = -30, res = 0.5, path = file.path("SpatialData", "Raster", "Worldclim"))
# Aus5 <- getData("worldclim", var = "bio", lon = 130, lat = -30, res = 0.5, path = file.path("SpatialData", "Raster", "Worldclim"))
# Aus6 <- getData("worldclim", var = "bio", lon = 150, lat = -30, res = 0.5, path = file.path("SpatialData", "Raster", "Worldclim"))
# Aus_top <- merge(Aus1, Aus2)
# Aus_top <- merge(Aus_top, Aus3)
# Aus_bot <- merge(Aus4, Aus5)
# Aus_bot <- merge(Aus_bot, Aus6)
# Aus_bio <- merge(Aus_top, Aus_bot)
# rm(Aus1,Aus2,Aus3,Aus4,Aus5,Aus6,Aus_top,Aus_bot)
# Aus_bio <- mask(Aus_bio, Aus_Coast)
# Aus_bio <- crop(Aus_bio, Aus_Coast)
# writeRaster(Aus_bio, filename = file.path("SpatialData", "Raster", "Worldclim", "Aus_bio.grd"), bandorder = "BIL", overwrite = T)
Aus_bio <- stack(file.path("SpatialData", "Raster", "Worldclim", "Aus_bio.gri"))

#crop wordclim to range map
Aus_bio_Ap <- crop(Aus_bio, extent(Ap_sp))

#Calculate distance-to-range
#range and domain need to be the same extent (thus same projection)
start_time <- Sys.time()
rdist <- rangeDist(
              range = Ap_sp,
              domain = Aus_bio_Ap,
              domainkm = 1000,
              mask = F,
              fact = 4
              )
end_time <- Sys.time()
end_time - start_time

#load rdist
load(file.path("SpatialData", "Raster", "rangeDist", "Ap_rdist.RData"))

#Crop "environmental data" to this new domain
Aus_bio_domain <- crop(Aus_bio, rdist)

# Mask pixels with no environmenta data (over ocean, etc.)
rdist <- mask(rdist, Aus_bio_domain[[1]])
names(rdist) <- "rangeDist"

#Evaluate curves
rates <- checkRates(rdist, skew = 0.2)

#Calculate frequency table of distances
dists <- freq(rdist, useNA = "no", digits = 2)

#Create offset
expert <- rangeOffset(
  rdist,
  parms = c(rate = 0.215, skew = 0.2, shift = 0, prob = 0.9), #need to decide on parameter values. If 0.9, rate should be higher
  dists = dists,
  doNormalize = T,
  verbose = T,
  doWriteRaster = T,
  filename = file.path("SpatialData", "Raster", "Ap_offset.tif"), #in a loop will use e.g. paste0(i, "_offset.tif")
  overwrite = T,
  datatype = "FLT4S"
)

#Can extend offset raster using extend(). Prob needed for zonation.
#example
expert <- extend(expert, Aus_Coast)
