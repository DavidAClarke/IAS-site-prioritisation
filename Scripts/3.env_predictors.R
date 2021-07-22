                  #Obtaining environmental data for use in bias file creation and SDMs
#Australia elevation
#getData("alt", country = "AUS", mask = T, path = file.path("SpatialData", "Raster", "Elevation"))
Aus_elev <- raster(file.path("SpatialData", "Raster", "Elevation", "Aus_msk_alt.gri"))
#Aus_elev <- raster("E:/SpatialData/Raster/Elevation/Aus_msk_alt.gri") #external drive
Aus_elev <- crop(Aus_elev, Aus_Coast)
Aus_elev <- mask(Aus_elev, Aus_Coast)
#for resampling veg
#Aus_elev_proj <- projectRaster(from = Aus_elev, res = 1000, crs = "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")


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
#raster.cor.plot(Aus_bio)

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
#Determine multicollinearity of predictors
#raster.cor.plot(Aus_bio_min)


#Global Lakes and Wetlands
GLWD_3 <- raster(file.path("SpatialData", "Raster", "glwd_3", "w001001.adf"))
GLWD_3 <- raster("C:/Users/dcla0008/Documents/GLWD-level3/glwd_3/w001001.adf") #work computer 2
crs(GLWD_3) <- crs(Aus_Coast)
GLWD_3 <- crop(GLWD_3, Aus_Coast)
rat <- levels(GLWD_3)[[1]]
rat$Type <- c("Lake","Reservoir", "River", "Freshwater Marsh, Floodplain",
              "Swamp Forest, Flooded Forest", 
              "Coastal Wetland (incl. Mangrove, Estuary, Delta, Lagoon)", 
              "Pan, Brackish/Saline Wetland", "Bog, Fen, Mire (Peatland)", 
              "Intermittent Wetland/Lake", "50-100% Wetland", "25-50% Wetland", 
              "Wetland Compex (0-25% Wetland)")
levels(GLWD_3) <- rat
GLWD_3 <- deratify(GLWD_3, "Type")
crs(GLWD_3) <- crs(Aus_Coast)

#Create distance to water layer (I think works)
GLWD <- spex::polygonize(GLWD_3) #convert to polygon
GLWD_proj <- st_transform(GLWD, 3577) #transform for dealing better with units
GLWD_grid <- st_make_grid(GLWD_proj, cellsize = 5000, what = "centers") #cellsize influences file size
GLWD_grid <- st_intersection(GLWD_grid, GLWD_proj)
GLWD_proj <- st_cast(GLWD_proj, "MULTILINESTRING")
save(GLWD_proj, file = file.path("SpatialData", "GLWD_proj.RData"))
save(GLWD_grid, file = file.path("SpatialData", "GLWD_grid.RData"))
dist <- st_distance(GLWD_proj, GLWD_grid) #big file!
df <- data.frame(dist = as.vector(dist)/1000, st_coordinates(GLWD_grid))
ext <- extent(as(GLWD_grid, "Spatial"))
r <- raster(resolution = 1000, ext = ext, crs = "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0
+units=m +no_defs")
dist_sf <- st_as_sf(df, coords = c("X","Y")) %>% st_set_crs(3577)
dist_raster <- rasterize(dist_sf, r, "dist", fun = "mean") #also try fasterize
writeRaster(dist_raster, filename = file.path("SpatialData", "dist_to_water.grd"), overwrite = T)

#Vegetation
# veg <- raster(file.path("SpatialData", "GRID_NVIS6_0_AUST_EXT_MVG", "aus6_0e_mvg", "w001000.adf"))
# veg <- resample(veg, Aus_elev_proj, method = "ngb")
# veg <- projectRaster(from = veg, to = Aus_elev, method = "ngb")
# veg <- mask(veg, Aus_Coast)
# veg <- as.factor(veg)
# rat <- levels(veg)[[1]]
# rat$Name <- c("Rainforests and Vine Thickets", "Eucalypt Tall Open Forests",
#               "Eucalypt Open Forests", "Eucalypt Low Open Forests", "Eucalypt Woodlands",
#               "Acacia Forests and Woodlands", "Callitris Forests and Woodlands",
#               "Casuarina Forests and Woodlands", "Melaleuca Forests and Woodlands",
#               "Other Forests and Woodlands", "Eucalypt Open Woodlands", "Tropical Eucalypt Woodlands/Grasslands",
#               "Acacia Open Woodlands", "Mallee Woodlands and Shrublands", "Low Closed Forests and Tall Closed Shrublands",
#              " Acacia Shrublands", "Other Shrublands", "Heathlands", "Tussock Grasslands", 
#               "Hummock Grasslands", "Other Grasslands, Herblands, Sedgelands and Rushlands", 
#               "Chenopod Shrublands, Samphire Shrublands and Forblands", "Mangroves",
#               "Inland Aquatic - freshwater, salt lakes, lagoons", "Cleared, non-native vegetation, buildings",
#               "Unclassified native vegetation", "Naturally bare - sand, rock, claypan, mudflat",
#               "Sea and estuaries", "Regrowth, modified native vegetation", "Unclassified forest",
#               "Other Open Woodlands", "Mallee Open Woodlands and Sparse Mallee Shrublands",
#               "Unknown/no data")
# levels(veg) <- rat
# #veg <- deratify(veg, "Name")
# rm(Aus_elev_proj)
#maybe need to separate each class into individual rasters, then stack them
#Here is an example of how to do it:
# New_LULC_2015 <- raster()
# 
# values(New_LULC_2015) <- sample(x = 1:15,size = ncell(New_LULC_2015),replace = T)
# 
# rat <- data.frame(
#   ID= 1:15,
#   LandCover = c("Evergreen","Deciduous","Mixed Forest",
#                 "grass","Degraded Forest","Scrubland",
#                 "Dry Grassland","Wet Grassland","Plantation",
#                 "Settlement","High Elevation Plantation","Coastal",
#                 "Barren Land","Cropland","Wetlands")
# )
# 
# levels(New_LULC_2015)<- rat
# 
# ls() # to see which objects are created in the workspace
# ## [1] "New_LULC_2015" "rat" 
# 
# # create class names without space
# rat$class_names <- gsub(pattern = ' ' ,replacement = '_',x = rat$LandCover)
# 
# for(i in 1:15){
#   assign(rat$class_names[i],mask(New_LULC_2015,New_LULC_2015 != i, maskvalue=1))
# }


#writeRaster(veg, filename = "F:/SpatialData/Raster/Aus_veg.grd", overwrite = T)
Aus_veg <- raster(file.path("SpatialData", "Raster", "Aus_veg.gri"))

#Surface Hydrology
#gdb_path <- file.path("SpatialData", "SurfaceHydrologyPolygonsNational.gdb")
#ogrListLayers(gdb_path)
# HP <- readOGR(gdb_path, "HydroPolys")
# HP_sf <- st_as_sf(HP)
# rm(HP)
# HP_sf <- st_transform(HP_sf, 4326)
# HP_sf <- HP_sf %>%
#   dplyr::select(FEATURETYPE, TYPE, SHAPE_Length, SHAPE_Area, geometry)
# #
# FeatureType <- HP_sf %>%
#   st_drop_geometry() %>%
#   distinct(FEATURETYPE = FEATURETYPE, .keep_all = T) %>%
#   arrange(TYPE) %>%
#   dplyr::select(FEATURETYPE)
# FeatureType <- as.vector(FeatureType$FEATURETYPE)
# #
# HP_sf <- st_make_valid(HP_sf)
# HP_sf <- st_crop(HP_sf, Aus_Coast)
# HP_sf <- st_cast(HP_sf, "MULTIPOLYGON")
# HP_rst <- fasterize(HP_sf, Aus_elev, field = "TYPE")
# HP_rst <- mask(HP_rst, Aus_Coast)
# #maybe need to separate each class into individual rasters, then stack them
# writeRaster(HP_rst, filename = "E:/SpatialData/Raster/Aus_Hyd.grd", overwrite = T)
HP_rst <- raster(file.path("SpatialData", "Raster", "Aus_Hyd.gri"))

#Combine
env_predictors <- stack(Aus_bio_min, Aus_elev, Aus_veg, HP_rst)
rm(Aus_bio, Aus_bio_min, Aus_elev, Aus_veg, HP_rst)
