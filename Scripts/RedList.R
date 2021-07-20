############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "RedList.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "RedList.R"),
              local_path = file.path("Scripts", "RedList.R"), overwrite = T)

############################### Load Red List spatial and assessment data #########################

#Load spatial data
#RL_Aus_shp <- st_read("SpatialData/Vector/redlist_species_data_Aus/data_0.shp")
load(file.path("SpatialData","Vector", "lvl_1.RData"))
lvl_1_sf <- st_as_sf(lvl_1)


#Updated downloaded shapefiles. Do I combine them all? How do they differ?
#RL_shp_0 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_0.shp")
RL_shp_0 <- st_read("F:/SpatialData/Vector/redlist_species_data_18052021/data_0.shp") #external drive
#RL_shp_1 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_1.shp")
RL_shp_1 <- st_read("F:/SpatialData/Vector/redlist_species_data_18052021/data_1.shp") #external drive
#RL_shp_2 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_2.shp")
RL_shp_2 <- st_read("F:/SpatialData/Vector/redlist_species_data_18052021/data_2.shp") #external drive
RL_shp <- rbind(RL_shp_0,RL_shp_1,RL_shp_2)
rm(RL_shp_0,RL_shp_1,RL_shp_2)

#Calculate area and proportion of total area for each species x country combination
#lvl_1_sf_moll <- st_transform(lvl_1_sf, crs = "+proj=moll") #Mollweide projection
#RL_shp_moll <- st_transform(RL_shp, crs = "+proj=moll") #Mollweide projection
#RL_shp_area <- calc_area(RL_shp_moll, "BINOMIAL", lvl_1_sf_moll) #need to edit function to remove RL_shp_moll from function
#write.csv(RL_shp_area, file.path("SpeciesData", "RangeProportions.csv"))

#Which species have Australian ranges?
Ranges <- read.csv(file.path("SpeciesData", "RangeProportions.csv"))
Ranges_Aus <- Ranges %>% dplyr::filter(NAME_0 == "Australia")
speciesNamesRange <- unique(Ranges_Aus$scientificName)
RL_shp_Aus <- RL_shp %>% dplyr::filter(BINOMIAL %in% speciesNamesRange)
rm(RL_shp)
#length(unique(RL_shp_Aus$BINOMIAL))

#Load plant point data
RL_plants <- read.csv(file.path("SpatialData", "Vector", "redlist_species_data_plantpoints", "points_data.csv"))
RL_plants <- RL_plants %>%
  dplyr::select("binomial", "presence", "origin", "seasonal", "legend", "longitude", "latitude") %>%
  mutate(Country = map.where(x = RL_plants$longitude, y = RL_plants$latitude)) %>%
  rename(Long = "longitude") %>%
  rename(Lat = "latitude") %>%
  rename(scientificName = "binomial")
RL_plants <- xy_match(RL_plants, RL_plants_new, lvl_1)
RL_plants <- RL_plants %>% drop_na(NAME_0) %>% filter(NAME_0 == "Australia")
write.csv(RL_plants, file.path("SpatialData", "Vector", "redlist_species_data_plantpoints", "points_data_filtered.csv"))
RL_plants <- read.csv(file.path("SpatialData", "Vector", "redlist_species_data_plantpoints", "points_data_filtered.csv"))
#RL_plants <- subset(RL_plants, with(RL_plants, unsplit(table(scientificName), scientificName)) >= 5)
#RL_plants_Aus <- RL_plants %>% filter(Country == "Australia" | Country == "Australia:Tasmania")

###########################################################################################
# #For potential Red List paper
# RL_plants_Aus <- filter(RL_plants, grepl("^Australia", Country))
# Aus_species <- unique(RL_plants_Aus$scientificName)
# #Removing species that have no Australian distribution
# Aus_RL_species <- RL_plants %>% filter(scientificName %in% Aus_species)
##############################################################################################

# RL_plants_Aus <- subset(RL_plants_Aus, with(RL_plants_Aus, unsplit(table(scientificName), scientificName)) >= 5)
# RL_plants_points <- RL_plants %>%
#   dplyr::select("Long", "Lat")
# RL_plants_points.sp <- SpatialPointsDataFrame(RL_plants_points,
#                                               data = RL_plants,
#                                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# RL_plants_points_sf <- st_as_sf(RL_plants_points.sp)
#rm(RL_plants_points.sp,RL_plants_points,RL_plants)



#Also need to load Red List assessment and threats data. Join to RL spatial layer #~#
#RL_Aus_assess <- read_excel(file.path("SpeciesData", "redlist_species_data", "assessments.xlsx"), sheet = "assessments")
#RL_Aus_threats <- read_excel(file.path("SpeciesData", "redlist_species_data", "threats.xlsx"), sheet = "threats")

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
  dplyr::mutate(phylumName = as.factor(phylumName)) %>%
  dplyr::mutate(kingdomName = as.factor(kingdomName)) %>%
  dplyr::filter(redlistCategory != "Lower Risk/near threatened", 
         redlistCategory != "Lower Risk/conservation dependent",
         redlistCategory != "Lower Risk/least concern") %>%
  dplyr::filter(!str_detect(name.y, "Marine")) %>% #this should exclude MOST marine species
  dplyr::mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern", "Near Threatened",
                                             "Vulnerable", "Endangered", "Critically Endangered", "Extinct"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", "CR", "EX")))
rm(RL_Aus_assess, RL_Aus_taxon, RL_Aus_threats, RL_Aus_habitats)

############################################ Taxonomic harmonization #####################################
##Shapefile
# species_names <- unique(RL_shp_Aus$BINOMIAL)
# TaxInfo <- traitdataform::get_gbif_taxonomy(species_names) #function from traitdataform package
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

#Which species have unresolved taxonomies?
#which(is.na(RL_shp_Aus$acceptedName == TRUE))
RL_shp_Aus[57,16] <- "Onthophagus froggattellus" #Onthophagus bicornis
RL_shp_Aus[3686,16] <- "Austrelaps ramsayi" #Austrelaps superbus
RL_shp_Aus[3744,16] <- "Gehyra montium" #Gehyra punctata
RL_shp_Aus[4047,16] <- "Eremiascincus phantasmus" #Eremiascincus fasciolatus
RL_shp_Aus[c(6932,6933),16] <- "Chelon planiceps" #Planiliza planiceps
RL_shp_Aus[c(6970,6971),16] <- "Ozimops planiceps"
RL_shp_Aus <- RL_shp_Aus %>% drop_na(acceptedName) #removing unresolved names
# length(unique(RL_shp_Aus$acceptedName))
# RL_Aus_Range <- st_drop_geometry(RL_shp_Aus)
# write.csv(RL_Aus_Range, file = file.path("SpeciesData", "RL_Aus_Range.csv"))

##Assessment information
#species_names_assess <- unique(RL_info$scientificName)
# TaxInfo_assess <- traitdataform::get_gbif_taxonomy(species_names_assess)
# TaxInfo_assess <- as_tibble(TaxInfo_assess)
# TaxInfo_assess_names <- TaxInfo_assess %>% 
#   rename(acceptedName = "scientificName") %>%
#   rename(scientificName = "verbatimScientificName") %>%
#   dplyr::select(scientificName, acceptedName)
#write.csv(TaxInfo_assess_names, file = file.path("SpeciesData", "TaxInfo_assess_names.csv"))
# TaxInfo_assess_names <- read.csv(file.path("SpeciesData", "TaxInfo_assess_names.csv"))
# RL_info <- RL_info %>%
#   left_join(TaxInfo_assess_names, by = "scientificName")

##Red List plants
# species_names_plants <- unique(RL_plants$scientificName)
# TaxInfo_plants <- traitdataform::get_gbif_taxonomy(species_names_plants)
# TaxInfo_plants <- as_tibble(TaxInfo_plants)
# TaxInfo_plants_names <- TaxInfo_plants %>% 
#  rename(acceptedName = "scientificName") %>%
#  rename(scientificName = "verbatimScientificName") %>%
#  dplyr::select(scientificName, acceptedName)
# write.csv(TaxInfo_plants_names, file = file.path("SpeciesData", "TaxInfo_plants_names.csv"))
# TaxInfo_plants_names <- read.csv(file.path("SpeciesData", "TaxInfo_plants_names.csv"))
# RL_plants <- RL_plants %>% 
#   left_join(TaxInfo_plants_names, by = "scientificName")
# write.csv(RL_plants, file = file.path("SpatialData", "Vector", "redlist_species_data_plantpoints","RL_plants_points.csv"))
RL_plants <- read.csv(file.path("SpatialData", "Vector", "redlist_species_data_plantpoints","RL_plants_points.csv"))

#Which species have unresolved taxonomies?
which(is.na(RL_plants$acceptedName == TRUE))


#also, jusdt noticed that Sus scrofa (feral pig?) is in the shapefile
#maybe run the species list against a list of known introduced species in Australia
#GRIIS? Is there a "comprehensive" list of introduced species in Australia?


############################### Downloading GBIF data (rgibf and Taxize) #####################################

#Getting number of GBIF occurrences
gbif_counts <- function(keys) {
  res.out <- lapply(keys[1:length(keys)], function(i) {
       occs <- occ_count(taxonKey = i,
                         basisOfRecord = c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION'),
                         georeferenced = T, 
                         country = "AU")
       usagekey <- i
       info <- cbind(usagekey, occs)
   })
  res.out <- do.call(rbind, res.out)

return(data.frame(res.out))
}

counts <- gbif_counts(gbif_taxon_keys)
matches <- left_join(matches, counts, by = "usagekey")
matches <- matches %>% distinct(original_sciname = original_sciname, .keep_all = T)
matches[5073,27] <- 0 # row for "Ozimops planiceps"
#write.csv(matches, file = file.path("SpeciesData", "all_matches.csv"))
matches <- read.csv(file.path("SpeciesData", "all_matches.csv"))

#Request download
res <- occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("basisOfRecord", c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION', 'PRESERVED_SPECIMEN')),
  pred("country", "AU"),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email)
occ_download_meta(res) #check the download status

#0299151-200613084148143 - this is the link for my download
#The following first downloads to machine and then imports into R
#occ_download_get("0299151-200613084148143") %>% occ_download_import()
#maybe assign to object e.g.
#dat <- occ_download_get("0299151-200613084148143", path = file.path("SpatialData", "Vector")) %>% 
#occ_download_import(fill = FALSE, quote = "", na.strings = c("", NA))
#The fill = FALSE and quote = "" are in anticipation of errors
#Then only select the most relevant columns
#I may issues as it is a BIG file (use VM)

#Importing GBIF data
#first get column names
my_cols <- fread(file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "0299151-200613084148143.csv"), nrows = 0)
my_cols_min <- c("taxonKey","family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince",
             "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode")
dataDir <- path.expand("~/Chapter_3/Chapter_3/SpatialData/Vector/SpeciesOccurrences/GBIF")
dataFiles <- dir(dataDir, pattern = "csv$", full.names = T)
gbif_occs <- fread(cmd = paste("grep", "Acanthiza.pusilla",dataFiles), col.names = names(my_cols)) #works!
gbif_occs_min <- gbif_occs %>% dplyr::select(my_cols_min)
#write.csv(gbif_occs_min, file = file.path("SpatialData", "Vector","SpeciesOccurrences", "GBIF", "gbif_occs_min.csv"))
Ap_gbif_occs <- read.csv(file.path("SpatialData", "Vector","SpeciesOccurrences", "GBIF", "gbif_occs_min.csv"))

#Summary information for Red List categories
RL_assess_sum <- RL_info %>%
  distinct(scientificName = scientificName, .keep_all = T) %>%
  group_by(className, redlistCategory) %>%
  count(redlistCategory)

#Summary information for Red List categories of IAS threatened species
RL_IAS_assess_sum <- RL_info %>%
  filter(code == "8.1.1" | code == "8.1.2" | code == "8.2.2") %>%
  distinct(scientificName = scientificName, .keep_all = T) %>%
  group_by(className, redlistCategory) %>%
  count(redlistCategory)

###################################### Spatial operations ##################################
#Creating convex polygons for plants
coordinates(RL_plants_Aus) <- ~ Long + Lat
Aus_plant_mcp <- mcp(RL_plants_Aus[,1], percent = 95)
Aus_plant_mcp_sf <- st_as_sf(Aus_plant_mcp)
Aus_plant_mcp_sf <- Aus_plant_mcp_sf %>% rename(scientificName = "id")
Aus_plant_mcp_sf <- merge(RL_info, Aus_plant_mcp_sf, by = "scientificName")
Aus_plant_mcp_sf <- st_as_sf(Aus_plant_mcp_sf)
st_crs(Aus_plant_mcp_sf) <- st_crs(Aus_Coast)
Aus_plant_mcp_sf_IAS <- Aus_plant_mcp_sf %>% filter(code == "8.1.1" | code == "8.1.2" | code == "8.2.2")



#Combining spatial and ancillary information 
RL_shp <- RL_shp %>%
  dplyr::select(BINOMIAL, PRESENCE, ORIGIN, SEASONAL, LEGEND) %>%
  rename(scientificName = "BINOMIAL")
new_bb <- c(113, -44, 154, -9)
names(new_bb) <- c("xmin", "ymin", "xmax", "ymax")
RL_Aus_shp <- st_crop(RL_shp, new_bb) #bringing extent down before projection
RL_Aus_shp <- merge(RL_info, RL_Aus_shp, by = "scientificName")
RL_Aus_shp <- st_as_sf(RL_Aus_shp) #now ready for further spatial analysis
RL_Aus_shp_IAS <- RL_Aus_shp %>% filter(code == "8.1.1" | code == "8.1.2" | code == "8.2.2")

#Re-project spatial layers
RL_Aus_shp_proj <- st_transform(RL_Aus_shp, 3577)
Aus_plant_mcp_sf_proj <- st_transform(Aus_plant_mcp_sf, 3577)
RL_Aus_shp_IAS_proj <- st_transform(RL_Aus_shp_IAS, 3577)
Aus_plant_mcp_sf_IAS_proj <- st_transform(Aus_plant_mcp_sf_IAS, 3577)

#Breaking it up
RL_Aus_shp_proj_1 <- RL_Aus_shp_proj[1:3823,]
RL_Aus_shp_proj_2 <- RL_Aus_shp_proj[3824:7642,]
RL_Aus_shp_proj_3 <- RL_Aus_shp_proj[7643:11470,]
RL_Aus_shp_proj_4 <- RL_Aus_shp_proj[11471:15284,]

#Intersection with coastal layer
Aus_RL_1 <- st_intersection(st_buffer(Aus_Coast_proj, dist = 0), st_buffer(RL_Aus_shp_proj_1, dist = 0))
Aus_RL_IAS <- st_intersection(st_buffer(Aus_Coast_proj, dist = 0), st_buffer(RL_Aus_shp_IAS_proj, dist = 0))
Aus_RL_plants <- st_intersection(st_buffer(Aus_Coast_proj, dist = 0), st_buffer(Aus_plant_mcp_sf_proj, dist = 0))

#Remove unnecessary files
rm(RL_Aus_shp)
rm(RL_Aus_shp_IAS)
rm(RL_Aus_shp_proj)
rm(RL_Aus_shp_IAS_proj)

#Species in each Red List shapefile
scientificNames <- unique(Aus_RL$scientificName)
scientificNames_IAS <- unique(Aus_RL_IAS$scientificName)

####################################### Calculating Red List Indices ##############################
#Red List Index for all species
#Worst = 0, Best = 1
#Based on single (latest) assessment
RLItotal <- RL_Aus_assess_min %>%
  mutate(redlistCategory = as.character(redlistCategory)) %>%
  summarise(class = "Total", RLI = rli(redlistCategory, boot = T, runs = 10000), species = length(redlistCategory))

#Red List Index for each taxonomic class
#Worst = 0, Best = 1
#Based on single (latest) assessment
RLIclass <- RL_Aus_assess_min %>% 
  group_by(class) %>% 
  mutate(redlistCategory = as.character(redlistCategory)) %>%
  summarise(RLI = rli(redlistCategory, boot = T, runs = 10000), species = length(redlistCategory))

#Joining all and class RLIs
RLI <- as_tibble(round(rbind(RLIclass$RLI[,1:3], RLItotal$RLI[,1:3]),3))
RLI <- RLI %>%
  mutate(class = c("Actinopterygii","Amphibia", "Arachnida","Aves","Bivalvia",
                   "Bryopsida","Chondrichthyes","Gastropoda", "Insecta","Liliopsida",
                   "Magnoliopsida","Malacostraca","Mammalia","Pinopsida", "Reptilia", 
                   "Sarcopterygii","Total"))

#Red List Index for all species threatened by IAS
#Worst = 0, Best = 1
#Based on single (latest) assessment
IAS_RLItotal <- RL_Aus_join_IAS %>%
  mutate(redlistCategory = as.character(redlistCategory)) %>%
  summarise(class = "Total", RLI = rli(redlistCategory, boot = T, runs = 10000), species = length(redlistCategory))

#Red List Index for each taxonomic class with IAS threats
#Worst = 0, Best = 1
#Based on single (latest) assessment
IAS_RLIclass <- RL_Aus_join_IAS %>% 
  group_by(class) %>% 
  mutate(redlistCategory = as.character(redlistCategory)) %>%
  summarise(RLI = rli(redlistCategory, boot = T, runs = 10000), species = length(redlistCategory))

#Joining all and class RLIs
IAS_RLI <- as_tibble(round(rbind(IAS_RLIclass$RLI[,1:3], IAS_RLItotal$RLI[,1:3]),3))
IAS_RLI <- IAS_RLI %>%
           mutate(class = c("Actinopterygii","Amphibia", "Aves","Insecta","Mammalia",
                            "Pinopsida","Reptilia","Sarcopterygii","Total"))



####################################### Figures ##########################################

#Red List plot
mycols <- c("#808080", "#008000", "#ADFF2F", "#FFFF00", "#FFA500", "#FF0000") #IUCN colours

#All data & all threats
All_species <- ggplot(RL_Aus_assess_sum, aes(x = redlistCategory, y = n, fill = redlistCategory)) +
                geom_bar(stat = "identity") +
                facet_wrap(~class) +
                theme_bw() +
                theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.text = element_text(size = 12), 
                 axis.title.x = element_text(size = 14),
                 axis.title.y = element_text(size = 14),
                 strip.text.x = element_text(size = 12)) +
               scale_fill_manual(values = mycols) +
               scale_y_continuous(trans = "log1p") +
               ylab("Number of species (log + 1)") +
               xlab("IUCN Red List category")

#IAS as Threat
IAS_species <- ggplot(RL_Aus_join_sum, aes(x = redlistCategory, y = n, fill = redlistCategory)) +
                geom_bar(stat = "identity") +
                facet_wrap(~class) +
                theme_bw() +
                theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.text = element_text(size = 12), 
                 axis.title.x = element_text(size = 14),
                 axis.title.y = element_text(size = 14),
                 strip.text.x = element_text(size = 12)) +
               scale_fill_manual(values = mycols) +
               scale_y_continuous(trans = "log1p") +
               ylab("Number of species (log + 1)") +
               xlab("IUCN Red List category")

#Red List Index plots (maybe include both in one plot, separated by colour)
#IAS threatened and all species
RLI_plot <- ggplot(RLI, aes(x=Median, y=reorder(class, -Median))) + #y=reorder(class, -Median) for decreasing median
              geom_point(stat="identity") +
              geom_errorbar(aes(xmin=LowCL, xmax=UpCL), width=.2) +
              geom_point(stat="identity", data = IAS_RLI, colour = "red") + #Total species RLI
              geom_errorbar(aes(xmin=LowCL, xmax=UpCL), data = IAS_RLI, width=.3, colour = "red") + #Total species RLI
              theme_bw() +
              theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size = 12), 
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                strip.text.x = element_text(size = 12),
                legend.position = "none") +
              scale_x_continuous(limits = c(0,1))+
              ylab("Taxonomic class") +
              xlab(paste("Red List Index", "\U00B1", "95% CI"))
