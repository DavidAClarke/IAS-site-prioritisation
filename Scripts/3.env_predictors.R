                  #Obtaining environmental data for use in bias file creation and SDMs
#Australia elevation
#getData("alt", country = "AUS", mask = T, path = file.path("SpatialData", "Raster", "Elevation"))
Aus_elev <- raster(file.path("SpatialData", "Raster", "Elevation", "Aus_msk_alt.gri"))
Aus_elev <- mask(Aus_elev, Aus_Coast)
Aus_elev <- crop(Aus_elev, Aus_Coast)

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

#Make sure any NAs propagate through the layers
Aus_bio <- check.env(Aus_bio)

#Determine multicollinearity of predictors
raster.cor.plot(Aus_bio)

#Pick variables to retain based on plots
# layers 6, 2, 4, 10? (maybe also 1, 12?)
Aus_bio <- Aus_bio[[c("", "", "")]]

#Combine
env_predictors <- stack(Aus_bio, Aus_elev)
