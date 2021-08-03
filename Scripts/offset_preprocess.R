#Load spatial data
Coast_shp <- st_read("/projects/nc57/Chapter_3/SpatialData/Vector/coastal-gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp")

#Only keep mainland and Tasmania
Aus_Coast <- Coast_shp %>% dplyr::filter(id == 6 | id == 32)
Aus_Coast <- st_make_valid(Aus_Coast)

#Removing original shapefile
rm(Coast_shp)

#Elevation
Aus_elev <- raster("/projects/nc57/Chapter_3/SpatialData/Raster/Elevation/Aus_elev_2.gri")
Aus_elev <- crop(Aus_elev, Aus_Coast)
Aus_elev <- mask(Aus_elev, Aus_Coast)

#Worldclim
Aus_bio <- stack("/projects/nc57/Chapter_3/SpatialData/Raster/Aus_bio.gri")
Aus_bio_min <- Aus_bio[[c("layer.1", "layer.2", "layer.3", "layer.4", "layer.5","layer.6","layer.7", 
                          "layer.12","layer.13", "layer.14")]]

#Vegetation
Aus_veg <- raster("/projects/nc57/Chapter_3/SpatialData/Raster/Aus_veg.gri")
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


Aus_veg_agg <- lapply(unique(Aus_veg), function(land_class) {
  
  aggregate(Aus_veg, fact = 5, fun = function(vals, na.rm){
    
    sum(vals == land_class, na.rm = na.rm)/length(vals)
  })
  
})
Aus_veg_agg <- stack(Aus_veg_agg)
names(Aus_veg_agg) <- rat$landcover
Aus_veg_agg <- mask(Aus_veg_agg, Aus_Coast)
Aus_veg_agg <- resample(Aus_veg_agg, Aus_elev ,method = "bilinear")
extent(Aus_veg_agg) <- extent(Aus_elev)

#Combine
env_predictors <- stack(Aus_bio_min, Aus_elev, Aus_veg_agg)
rm(Aus_bio, Aus_bio_min, Aus_elev, Aus_veg_agg, Aus_veg)
gc()

##########################################################################################################################################
RL_shp_0 <- st_read("/projects/nc57/Chapter_3/SpatialData/Vector/redlist_species_data_18052021/data_0.shp")
RL_shp_1 <- st_read("/projects/nc57/Chapter_3/SpatialData/Vector/redlist_species_data_18052021/data_1.shp")
RL_shp_2 <- st_read("/projects/nc57/Chapter_3/SpatialData/Vector/redlist_species_data_18052021/data_2.shp")
RL_shp <- rbind(RL_shp_0,RL_shp_1,RL_shp_2)
rm(RL_shp_0,RL_shp_1,RL_shp_2)

Ranges <- read.csv("/projects/nc57/Chapter_3/SpeciesData/RangeProportions.csv")
Ranges_Aus <- Ranges %>% dplyr::filter(NAME_0 == "Australia")
speciesNamesRange <- unique(Ranges_Aus$scientificName)
RL_shp_Aus <- RL_shp %>% dplyr::filter(BINOMIAL %in% speciesNamesRange)
rm(RL_shp)
gc()

RL_Aus_assess <- read.csv("/projects/nc57/Chapter_3/SpeciesData/redlist_species_data_18052021/assessments.csv")
RL_Aus_threats <- read.csv("/projects/nc57/Chapter_3/SpeciesData/redlist_species_data_18052021/threats.csv")
RL_Aus_taxon <- read.csv("/projects/nc57/Chapter_3/SpeciesData/redlist_species_data_18052021/taxonomy.csv")
RL_Aus_habitats <- read.csv("/projects/nc57/Chapter_3/SpeciesData/redlist_species_data_18052021/habitats.csv")

RL_info <- RL_Aus_assess %>% 
  left_join(RL_Aus_threats, by = "scientificName") %>%
  left_join(RL_Aus_taxon, by = "scientificName") %>%
  left_join(RL_Aus_habitats, by = "scientificName") %>%
  dplyr::select(scientificName, orderName, className, phylumName, kingdomName, 
                redlistCategory, redlistCriteria, yearPublished,code.x, name.x, 
                stressCode, stressName, ias, scope, name.y) %>%
  dplyr::mutate(className = as.factor(className)) %>%
  dplyr::mutate(orderName = as.factor(orderName)) %>%
  dplyr:: mutate(phylumName = as.factor(phylumName)) %>%
  dplyr:: mutate(kingdomName = as.factor(kingdomName)) %>%
  dplyr::filter(scientificName %in% speciesNamesRange) %>% 
  dplyr::filter(redlistCategory != "Lower Risk/near threatened", 
                redlistCategory != "Lower Risk/conservation dependent",
                redlistCategory != "Lower Risk/least concern") %>%
  dplyr::filter(!str_detect(name.y, "Marine")) %>% #this should exclude MOST marine species
  dplyr::filter(!str_detect(name.y, "Inter-Reef Rubble Substrate")) %>%
  dplyr::filter(!str_detect(name.y, "Inter-Reef Soft Substrate" )) %>%
  dplyr::filter(!str_detect(name.y, "Outer Reef Channel")) %>%
  dplyr:: filter(!str_detect(name.y, "Back Slope")) %>%
  dplyr::filter(!str_detect(name.y, "Foreslope")) %>%
  dplyr::filter(!str_detect(name.y, "Lagoon")) %>%
  dplyr::mutate(redlistCategory = factor(redlistCategory, 
                                         levels = c("Data Deficient", "Least Concern", "Near Threatened",
                                                    "Vulnerable", "Endangered", "Critically Endangered", "Extinct"),
                                         labels = c("DD", "LC", "NT", "VU", "EN", "CR", "EX")))
rm(RL_Aus_assess, RL_Aus_taxon, RL_Aus_threats, RL_Aus_habitats)
RL_info_names <- unique(RL_info$scientificName)

RL_shp_Aus <- RL_shp_Aus %>%
  filter(BINOMIAL %in% RL_info_names)
TaxInfo <- read.csv("/projects/nc57/Chapter_3/SpeciesData/TaxInfo_shp.csv")
TaxInfo <- as_tibble(TaxInfo)
TaxInfo_names <- TaxInfo %>% 
  distinct(verbatimScientificName = verbatimScientificName, .keep_all = T) %>%
  rename(acceptedName = "scientificName") %>%
  rename(BINOMIAL = "verbatimScientificName") %>%
  dplyr::select(BINOMIAL, acceptedName)
RL_shp_Aus <- RL_shp_Aus %>%
  left_join(TaxInfo_names, by = "BINOMIAL")
rm(TaxInfo_names, TaxInfo)
gc()

RL_shp_Aus[57,16] <- "Onthophagus froggattellus" #Onthophagus bicornis
RL_shp_Aus[3591,16] <- "Austrelaps ramsayi" #Austrelaps superbus
RL_shp_Aus[3649,16] <- "Gehyra montium" #Gehyra punctata
RL_shp_Aus[3950,16] <- "Eremiascincus phantasmus" #Eremiascincus fasciolatus
RL_shp_Aus[c(6716,6717),16] <- "Chelon planiceps" #Planiliza planiceps
RL_shp_Aus[c(6750,6751),16] <- "Ozimops planiceps"
RL_shp_Aus <- st_make_valid(RL_shp_Aus)
RL_shp_Aus <- RL_shp_Aus %>% 
  drop_na(acceptedName) %>% #removing unresolved names
  dplyr::select(BINOMIAL, PRESENCE, ORIGIN, SEASONAL, LEGEND, acceptedName) %>%
  st_intersection(Aus_Coast)
speciesNames_range <- unique(RL_shp_Aus$acceptedName)
speciesNames_range <- speciesNames_range[1:2] #include range
