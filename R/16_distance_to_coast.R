#Get map of Aus coast
#source("R/02_aus_coast.R") if not loaded

pm <- raster(file.path("Zonation", "species_CAZ_proj", "species_CAZ", 
                       "species_CAZ_out", 
                       "species_CAZ.CAZ_E.rank.compressed.tif"))

#Transform to 3577
Aus_Coast_proj <- st_transform(Aus_Coast, 3577)

#create the fishnet
grid <- st_make_grid(Aus_Coast_proj, cellsize = 5000, what = "centers")

#only extract the points in the limits of Australia
grid <- st_intersection(grid, Aus_Coast_proj)

#transform Australia from polygon shape to line
Aus_Coast_proj <- st_cast(Aus_Coast_proj, "MULTILINESTRING")

#calculation of the distance between the coast and our points
dist <- st_distance(Aus_Coast_proj, grid)

#create a data.frame with the distance and the coordinates of the points
df <- data.frame(dist = as.vector(dist)/1000,
                 st_coordinates(grid))


#Create empty aus raster the same specs
ext <- extent(as(grid, "Spatial"))
r <- raster(resolution = 5000, ext = ext, crs = crs(Aus_Coast_proj))

#convert the points to a spatial object class sf
dist_sf <- st_as_sf(df, coords = c("X", "Y")) %>%
  st_set_crs(3577)

#create the distance raster
dist_raster <- rasterize(dist_sf, r, "dist", fun = mean) #r is empty raster

#Re-project back to WGS84
dist_raster_2 <- projectRaster(dist_raster,pm, res = res(pm), 
                               crs = crs(Aus_Coast))

#export the raster
writeRaster(dist_raster_2, file = file.path("SpatialData", "Raster", 
                                            "dist_aus_coast.tif"), 
            format = "GTiff", overwrite = TRUE)
