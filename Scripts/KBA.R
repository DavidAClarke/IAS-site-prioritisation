## Loading and pre-processing of Key Biodiversity Areas (KBA) ##

#Load spatial data
KBA_shp <- st_read("SpatialData/Vector/KBA2019_no_overlap/KBA2019_no_overlap.shp")

#Preparing attribute data
KBA_shp <- KBA_shp %>%
  dplyr::select(Name, STATE, THREAT_STA, IBA_CRITER, HECTARES, geometry) %>%
  filter(STATE != "Other Territories")

#Crop to Aus Coast
#KBA_shp <- st_crop(st_make_valid(KBA_shp), Aus_Coast)

#Projections (GDA94 = 3577)
KBA_shp_proj <- st_transform(KBA_shp, 3577)

#Remove unprojected shapefile
rm(KBA_shp)

##Intersection with projected Australian coast layer
Aus_KBA <- st_intersection(st_buffer(Aus_Coast_proj, dist = 0), st_buffer(KBA_shp_proj, dist = 0))

#Create raster template
#rst_template <- raster(resolution = 1000,
#                       crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#                       ext = extent(Aus_Coast_proj)) #this extent may change. Need to cater to file with smallest extent

#Rasterize projected KBA shapefile
#rst_Aus_KBA <- rasterize(st_zm(Aus_KBA), rst_template, field = 1) #1 = presence of KBA

#Write raster
#dir.create(file.path("SpatialData", "Input_zonation")) #this part may get moved to general analysis script
#writeRaster(rst_Aus_KBA, filename = file.path("SpatialData", "Input_zonation", "KBA.tif"), format = "GTiff")


#Map of threat status
tmap_mode("plot")
KBA_threatSta <- tm_shape(Aus_Coast_proj) +
  tm_fill("white") +
  tm_borders("black", lwd = 0.5) +
  tm_shape(Aus_KBA) +
  tm_polygons("THREAT_STA", palette = "-viridis", alpha = 0.8)