############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "CAPAD.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"),
              local_path = file.path("Scripts", "CAPAD.R"))

####################################### CAPAD Protected Areas ######################################
#Load spatial data
CAPAD.shp <- st_read("SpatialData/Vector/CAPAD2018_terrestrial/CAPAD2018_terrestrial.shp")

CAPAD.shp <- CAPAD.shp %>%
  dplyr::select(`NAME`, `IUCN`, `EPBC`, `SHAPE_AREA`, `geometry`) %>%
  filter(IUCN != "III" & IUCN != "VI" & IUCN != "NA" & IUCN != "NAS")

# Re-project
CAPAD_proj <- st_transform(CAPAD.shp, 3577)

# Intersect with coast
Aus_CAPAD <- st_intersection(st_buffer(CAPAD_proj, dist = 0), st_buffer(Aus_Coast_proj, dist = 0))

# Alter feature names
Aus_CAPAD$NAME <- gsub(" |/", "_", Aus_CAPAD$NAME)

# Number and area of categories (entire shapefile)
IUCN_sum <- CAPAD.shp %>% 
  st_drop_geometry() %>%
  select(IUCN, SHAPE_AREA) %>%
  group_by(IUCN) %>%
  summarise(cats = n(), Area = sum(SHAPE_AREA))