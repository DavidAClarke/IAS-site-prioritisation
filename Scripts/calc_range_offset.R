############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "calc_range_offset.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "calc_range_offset.R"),
              local_path = file.path("Scripts", "calc_range_offset"), overwrite = T)
##################################################################################################
                              # Creating bias file (with bossMaps) #
                                              # TEST #
#selecting species
Ap <- RL_Aus_shp %>% filter(BINOMIAL == "Acanthiza pusilla")

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

#Calculate distance-to-range
#range and domain need to be the same extent (thus same projection)
rdist <- rangeDist(
              range = Ap_sp,
              domain = Aus_elev_Ap,
              domainkm = 1000,
              mask = F,
              fact = 2
              )

#Crop "environmental data" to this new domain
elev_domain = crop(Aus_elev_Ap, rdist)

# Mask pixels with no environmenta data (over ocean, etc.)
rdist <- mask(rdist, elev_domain)
names(rdist) <- "rangeDist"

#Evaluate curves
rates <- checkRates(rdist, skew = 0.2)
rates_posdif <- rates %>% filter(pdif >0 & prob >= 0.5) #row 15 here looks good?
