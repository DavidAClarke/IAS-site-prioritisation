############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "ECNES.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts","ECNES.R"),
              local_path = file.path("Scripts", "ECNES.R"))

########################################## Load data ###########################################
ECNES_shp <- st_read("SpatialData/Vector/ECnes_public_jun2020.shp/shapefile/ECnes_public_jun2020.shp")
ECNES_shp_min <- ECNES_shp %>%
  dplyr::select(community, epbc, pres_rank, SHAPE_Area, geometry) %>%
  filter(pres_rank == 2)
rm(ECNES_shp)

####################################### Spatial pre-processing ###############################
ECNES_shp_min_proj <- st_transform(ECNES_shp_min, 3577)

Aus_ECNES <- st_intersection(st_buffer(ECNES_shp_min_proj, dist = 0), st_buffer(Aus_Coast_proj, dist = 0))
