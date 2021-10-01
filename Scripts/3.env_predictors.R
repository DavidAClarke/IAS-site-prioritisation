                  #Obtaining environmental data for use in bias file creation and SDMs
#Australia elevation
#getData("alt", country = "AUS", mask = T, path = file.path("SpatialData", "Raster", "Elevation"))
#Aus_elev <- raster(file.path("SpatialData", "Raster", "Elevation", "Aus_msk_alt.gri"))
Aus_elev <- raster("E:/SpatialData/Raster/Elevation/Aus_msk_alt.gri") #external drive
Aus_elev <- crop(Aus_elev, Aus_Coast)
Aus_elev <- mask(Aus_elev, Aus_Coast)

#for resampling vegetation layer
Aus_elev_proj <- projectRaster(from = Aus_elev, res = 1000, crs = "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")


#Bioclim rasters. Getting resolution of 0.5 for enttirety of Australia
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
#1 = Annual Mean Temperature
#2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#3 = Isothermality (BIO2/BIO7) (×100)

#5 = Max Temperature of Warmest Month
#6 = Min Temperature of Coldest Month
#7 = Temperature annual range
#14 = Precipitation of Driest Month
#13 = Precipitation of Wettest Month
Aus_bio_min <- Aus_bio[[c("layer.1", "layer.2", "layer.3", "layer.4", "layer.5","layer.6","layer.7", 
                          "layer.12","layer.13", "layer.14")]]


#Australian Vegetation
veg <- raster(file.path("SpatialData", "GRID_NVIS6_0_AUST_EXT_MVG", "aus6_0e_mvg", "w001000.adf"))
veg <- resample(veg, Aus_elev_proj, method = "ngb")
veg <- projectRaster(from = veg, to = Aus_elev, method = "ngb")
veg <- mask(veg, Aus_Coast)
veg <- as.factor(veg)
rat <- levels(veg)[[1]]
rat$Name <- c("Rainforests and Vine Thickets", "Eucalypt Tall Open Forests",
              "Eucalypt Open Forests", "Eucalypt Low Open Forests", "Eucalypt Woodlands",
              "Acacia Forests and Woodlands", "Callitris Forests and Woodlands",
              "Casuarina Forests and Woodlands", "Melaleuca Forests and Woodlands",
              "Other Forests and Woodlands", "Eucalypt Open Woodlands", "Tropical Eucalypt Woodlands/Grasslands",
              "Acacia Open Woodlands", "Mallee Woodlands and Shrublands", "Low Closed Forests and Tall Closed Shrublands",
             "Acacia Shrublands", "Other Shrublands", "Heathlands", "Tussock Grasslands",
              "Hummock Grasslands", "Other Grasslands, Herblands, Sedgelands and Rushlands",
              "Chenopod Shrublands, Samphire Shrublands and Forblands", "Mangroves",
              "Inland Aquatic - freshwater, salt lakes, lagoons", "Cleared, non-native vegetation, buildings",
              "Unclassified native vegetation", "Naturally bare - sand, rock, claypan, mudflat",
              "Sea and estuaries", "Regrowth, modified native vegetation", "Unclassified forest",
              "Other Open Woodlands", "Unknown")
levels(veg) <- rat
rm(Aus_elev_proj)
writeRaster(veg, filename = "F:/SpatialData/Raster/Aus_veg.grd", overwrite = T)
Aus_veg <- raster(file.path("SpatialData", "Raster", "Aus_veg.gri"))

# reclassify
Aus_veg[Aus_veg %in% c(1:4, 30)] <- 1 #1, 2,3,4, 30 - Open forests/Rainforests and Vine Thickets
Aus_veg[Aus_veg %in% c(5:13, 31:32)] <- 2 #5:13, 31,32 - Woodlands
Aus_veg[Aus_veg %in% c(14:17)] <- 3 #14:17 - Shrublands
Aus_veg[Aus_veg == 18] <- 4 #18 - Heathlands
Aus_veg[Aus_veg %in% c(19:22)] <- 5 #19:22 - Grasslands
Aus_veg[Aus_veg %in% c(23,24,28)] <- 6 #23,24, 28 - Aquatic
Aus_veg[Aus_veg %in% c(26,29)] <- 7 #29, 26 - Native vegetation
Aus_veg[Aus_veg %in% c(25,27, 33)] <- 8 #25, 27 - Disturbed/bare/unknown

# set new RAT for reclassfied raster
Aus_veg <- raster::ratify(Aus_veg)
rat <- data.frame(
  ID = 1:8,
  landcover = c("Open forests/Rainforests and Vine Thickets", 
                "Woodlands",
                "Shrublands",
                "Heathlands",
                "Grasslands",
                "Aquatic",
                "Native vegetation",
                "Disturbed/bare/unknown")
)
levels(Aus_veg) <- rat

#Calculate proportion of cover
Aus_veg_agg <- lapply(unique(Aus_veg), function(land_class) {

  aggregate(Aus_veg, fact = 5, fun = function(vals, na.rm){

    sum(vals == land_class, na.rm = na.rm)/length(vals)
  })

})
Aus_veg_agg <- stack(Aus_veg_agg)
names(Aus_veg_agg) <- rat$landcover
Aus_veg_agg <- mask(Aus_veg_agg, Aus_Coast)
Aus_veg_agg <- resample(Aus_veg_agg, Aus_bio_2.5 ,method = "bilinear")

#Combine
env_predictors <- stack(Aus_bio_min, Aus_elev, Aus_veg_agg)
writeRaster(env_predictors, filename = file.path("SpatialData", "Raster", "env_predictors.grd"))
rm(Aus_bio, Aus_bio_min, Aus_elev, Aus_veg_agg)
