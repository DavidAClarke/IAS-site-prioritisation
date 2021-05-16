# Ramsar and upstream catchments

# Import spatial data
Ramsar.shp <- st_read("SpatialData/Vector/ramsar_wetlands/ramsar_wetlands.shp")
Upstream.shp <- st_read("SpatialData/Vector/ramsar_upstrm_catch/ramsar_upstrm_catch.shp")

# Write Ramsar csv
Ramsar <- st_drop_geometry(Ramsar.shp)
write.csv(Ramsar, file = "SpatialData/Vector/ramsar_wetlands/ramsar_wetlands.csv")

# Load Ramsar data (now with WEIGHTS)
Ramsar_crit <- read_csv("SpatialData/Vector/ramsar_wetlands/ramsar_wetlands.csv")
Ramsar_crit <- Ramsar_crit %>% dplyr::select("RAMSAR_NAM","CRITERIA", "%CRITERIA", "WEIGHT")

# Join criteria and weight data to shapefile
Ramsar.shp <- left_join(Ramsar.shp, Ramsar_crit, by = "RAMSAR_NAM")
