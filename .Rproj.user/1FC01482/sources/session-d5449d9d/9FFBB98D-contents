###################### Load Red List spatial and assessment data ###############

#Updated downloaded shapefiles. 
RL_shp_0 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_0.shp")
RL_shp_1 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_1.shp")
RL_shp_2 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_2.shp")
RL_shp <- rbind(RL_shp_0,RL_shp_1,RL_shp_2)
rm(RL_shp_0,RL_shp_1,RL_shp_2)

#Calculate area and proportion of total area for each species x country combination
#lvl_1_sf_moll <- st_transform(lvl_1_sf, crs = "+proj=moll") #Mollweide projection
#RL_shp_moll <- st_transform(RL_shp, crs = "+proj=moll") #Mollweide projection
##need to edit function to remove RL_shp_moll from function
#RL_shp_area <- calc_area(RL_shp_moll, "BINOMIAL", lvl_1_sf_moll) 
#write.csv(RL_shp_area, file.path("SpeciesData", "RangeProportions.csv"))

#Which species have Australian ranges?
Ranges <- read.csv(file.path("SpeciesData", "RangeProportions.csv"))
Ranges_Aus <- Ranges %>% dplyr::filter(NAME_0 == "Australia")
speciesNamesRange <- unique(Ranges_Aus$scientificName)
RL_shp_Aus <- RL_shp %>% dplyr::filter(BINOMIAL %in% speciesNamesRange)
rm(RL_shp)
#length(unique(RL_shp_Aus$BINOMIAL))

#More recent download - broader scope
RL_Aus_assess <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "assessments.csv"))
RL_Aus_threats <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "threats.csv"))
RL_Aus_taxon <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "taxonomy.csv"))
RL_Aus_habitats <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "habitats.csv"))

####################################### Data Preparation ####################################
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

######################### Taxonomic harmonization ##############################
##Shapefile
RL_shp_Aus <- RL_shp_Aus %>%
  dplyr::filter(BINOMIAL %in% RL_info_names)
# species_names <- unique(RL_shp_Aus$BINOMIAL)
# TaxInfo <- traitdataform::get_gbif_taxonomy(species_names) 
# write.csv(TaxInfo, file = file.path("SpeciesData", "TaxInfo_shp.csv"))
TaxInfo <- read.csv(file.path("SpeciesData", "TaxInfo_shp.csv"))
TaxInfo <- as_tibble(TaxInfo)
TaxInfo_names <- TaxInfo %>% 
  distinct(verbatimScientificName = verbatimScientificName, .keep_all = T) %>%
  rename(acceptedName = "scientificName") %>%
  rename(BINOMIAL = "verbatimScientificName") %>%
  dplyr::select(BINOMIAL, acceptedName)
RL_shp_Aus <- RL_shp_Aus %>%
  left_join(TaxInfo_names, by = "BINOMIAL")
rm(TaxInfo_names, TaxInfo)

#Which species have unresolved taxonomies?
#which(is.na(RL_shp_Aus$acceptedName == TRUE))
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
#write.csv(speciesNames_range, file = file.path("SpeciesData", "RL_shp_Aus.csv"))
# length(unique(RL_shp_Aus$acceptedName))
# RL_Aus_Range <- st_drop_geometry(RL_shp_Aus)
# write.csv(RL_Aus_Range, file = file.path("SpeciesData", "RL_Aus_Range.csv"))

