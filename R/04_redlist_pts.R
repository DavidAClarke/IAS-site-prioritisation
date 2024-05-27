################## Load Red List spatial and assessment data ###################

#Load spatial data
#RL_Aus_shp <- st_read("SpatialData/Vector/redlist_species_data_Aus/data_0.shp")
load(file.path("SpatialData","Vector", "lvl_1.RData"))
lvl_1_sf <- st_as_sf(lvl_1)

#Load plant point data
RL_plants <- read.csv(file.path("SpatialData", "Vector", 
                                "redlist_species_data_plantpoints", 
                                "points_data.csv"))
RL_plants <- RL_plants %>%
  dplyr::select("binomial", "presence", "origin", "seasonal", "legend", 
                "longitude", "latitude") %>%
  mutate(Country = map.where(x = RL_plants$longitude, 
                             y = RL_plants$latitude)) %>%
  rename(Long = "longitude") %>%
  rename(Lat = "latitude") %>%
  rename(scientificName = "binomial")

RL_plants <- xy_match(RL_plants, RL_plants_new, lvl_1, Long, Lat)
RL_plants <- RL_plants %>% drop_na(NAME_0) %>% filter(NAME_0 == "Australia")
write.csv(RL_plants, file.path("SpatialData", "Vector", 
                               "redlist_species_data_plantpoints", 
                               "points_data_filtered.csv"))
RL_plants <- read.csv(file.path("SpatialData", "Vector", 
                                "redlist_species_data_plantpoints", 
                                "points_data_filtered.csv"))
# RL_plants <- subset(RL_plants, with(RL_plants, unsplit(table(scientificName), 
#                                                        scientificName)) >= 5)
# RL_plants_Aus <- RL_plants %>% 
#   filter(Country == "Australia" | Country == "Australia:Tasmania")

################################################################################
# #For potential Red List paper
# RL_plants_Aus <- filter(RL_plants, grepl("^Australia", Country))
# Aus_species <- unique(RL_plants_Aus$scientificName)
# #Removing species that have no Australian distribution
# Aus_RL_species <- RL_plants %>% filter(scientificName %in% Aus_species)
################################################################################

# RL_plants_Aus <- subset(RL_plants_Aus,
#                         with(RL_plants_Aus, unsplit(table(scientificName), 
#                                                     scientificName)) >= 5)
# RL_plants_points <- RL_plants %>%
#   dplyr::select("Long", "Lat")
# RL_plants_points.sp <- SpatialPointsDataFrame(RL_plants_points,
#                                               data = RL_plants,
#                                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# RL_plants_points_sf <- st_as_sf(RL_plants_points.sp)
#rm(RL_plants_points.sp,RL_plants_points,RL_plants)


#need Red List assessment and threats data. Join to RL spatial layer
RL_Aus_assess <- read.csv(file.path("SpeciesData", 
                                    "redlist_species_data_18052021", 
                                    "assessments.csv"))
RL_Aus_threats <- read.csv(file.path("SpeciesData", 
                                     "redlist_species_data_18052021", 
                                     "threats.csv"))
RL_Aus_taxon <- read.csv(file.path("SpeciesData", 
                                   "redlist_species_data_18052021", 
                                   "taxonomy.csv"))
RL_Aus_habitats <- read.csv(file.path("SpeciesData", 
                                      "redlist_species_data_18052021", 
                                      "habitats.csv"))

############################# Data Preparation #################################
#Joining taxonomic information to red list assessment data
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
  dplyr::filter(!str_detect(name.y, "Marine")) %>% #exclude MOST marine species
  dplyr::filter(!str_detect(name.y, "Inter-Reef Rubble Substrate")) %>%
  dplyr::filter(!str_detect(name.y, "Inter-Reef Soft Substrate" )) %>%
  dplyr::filter(!str_detect(name.y, "Outer Reef Channel")) %>%
  dplyr:: filter(!str_detect(name.y, "Back Slope")) %>%
  dplyr::filter(!str_detect(name.y, "Foreslope")) %>%
  dplyr::filter(!str_detect(name.y, "Lagoon")) %>%
  dplyr::mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern", 
                                             "Near Threatened","Vulnerable", 
                                             "Endangered", 
                                             "Critically Endangered","Extinct"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", "CR", 
                                             "EX")))
rm(RL_Aus_assess, RL_Aus_taxon, RL_Aus_threats, RL_Aus_habitats)
RL_info_names <- unique(RL_info$scientificName)

####################### Taxonomic harmonization ################################
##Red List plants
species_names_plants <- unique(RL_plants$scientificName)
TaxInfo_plants <- traitdataform::get_gbif_taxonomy(species_names_plants)
TaxInfo_plants <- as_tibble(TaxInfo_plants)
TaxInfo_plants_names <- TaxInfo_plants %>%
 rename(acceptedName = "scientificName") %>%
 rename(scientificName = "verbatimScientificName") %>%
 dplyr::select(scientificName, acceptedName)
write.csv(TaxInfo_plants_names, file = file.path("SpeciesData", 
                                                 "TaxInfo_plants_names.csv"))
TaxInfo_plants_names <- read.csv(file.path("SpeciesData", 
                                           "TaxInfo_plants_names.csv"))
RL_plants <- RL_plants %>%
  left_join(TaxInfo_plants_names, by = "scientificName")
write.csv(RL_plants, file = file.path("SpatialData", "Vector", 
                                      "redlist_species_data_plantpoints",
                                      "RL_plants_points.csv"))
RL_plants <- read.csv(file.path("SpatialData", "Vector", 
                                "redlist_species_data_plantpoints",
                                "RL_plants_points.csv"))

#Which species have unresolved taxonomies?
which(is.na(RL_plants$acceptedName == TRUE))

