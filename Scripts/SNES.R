############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "SNES.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "SNES.R"),
              local_path = file.path("Scripts", "SNES.R"), overwrite = T)

############################### Load SNES spatial and assessment data #########################

#Load spatial data
SNES_shp <- st_read("SpatialData/Vector/snes_public_grids_08Aug2019_shapefile/snes_species_combined.shp")
SNES_shp_min <- SNES_shp %>% 
  drop_na(threatened) %>% 
  filter(pres_rank == "2" & threatened != "Extinct in the wild") %>% # 1 = species may occur; 2 = species likely occur
  dplyr::select(sci_name, threatened, tax_group, tax_class) %>%
  rename(scientificName = "sci_name") %>%
  rename(threatStatus = "threatened") %>%
  rename(class = "tax_class") %>%
  rename(group = "tax_group") %>%
  mutate(threatStatus = factor(threatStatus, 
                               levels = c("Conservation Dependent","Vulnerable", 
                                          "Endangered", "Critically Endangered"),
                               labels = c("CD", "VU", "EN", "CR"))) %>%
  mutate(class = as.factor(class)) %>%
  mutate(group = as.factor(group))
rm(SNES_shp)
  
    
#Load assessment data
SNES <- read_csv(file.path("Data", "AUS_ThreatenedSpeciesList_04092020.csv")) 
SNES_min <- SNES %>%
  dplyr::select(`Scientific Name`, `Threatened status`, `Class` ,`Profile`) %>% #add URLs column
  rename(scientificName = "Scientific Name") %>%
  rename(threatStatus = "Threatened status") %>%
  rename(class = "Class") %>%
  mutate(threatStatus = as.factor(threatStatus)) %>%
  mutate(class = as.factor(class)) %>%
  filter(threatStatus != "Extinct" & threatStatus != "Extinct in the wild")
rm(SNES)
####################################### Summary information ##################################
#Shapefile information
SNES_threat_sum <- SNES_shp_min %>%
  st_drop_geometry() %>%
  group_by(group, threatStatus) %>% 
  count(threatStatus)

#Assessment information


#################################### Obtaining threat information ############################
# First determine if a given SPRAT page has a "Threats" section

#NOTE: the SPRAT for Zyzomys pedunculatus shows that this first function is not very useful.
#      it has threat information, including on invasives, under a different heading.
#      perhaps only second function (threats_info) is to be used

# two arguments: SPRAT urls and associated species names
#Because of some issues with timeouts I am doing this piecewise
#ThreatsPresent <- Threat_present(SNES_min$Profile,SNES_min$scientificName, target_path = "websites")
#write.csv(ThreatsPresent_1800, file = file.path("Data", "ThreatsPresent.csv"))

# Load and join (rbind) different ThreatsPresent_ csv's. Name the joined files ThreatsPresent for next functions
# ThreatsPresent_200 <- read.csv(file.path("Data", "ThreatsPresent_200.csv"))
# ThreatsPresent_400 <- read.csv(file.path("Data", "ThreatsPresent_400.csv"))
# ThreatsPresent_600 <- read.csv(file.path("Data", "ThreatsPresent_600.csv"))
# ThreatsPresent_800 <- read.csv(file.path("Data", "ThreatsPresent_800.csv"))
# ThreatsPresent_1000 <- read.csv(file.path("Data", "ThreatsPresent_1000.csv"))
# ThreatsPresent_1200 <- read.csv(file.path("Data", "ThreatsPresent_1200.csv"))
# ThreatsPresent_1400 <- read.csv(file.path("Data", "ThreatsPresent_1400.csv"))
# ThreatsPresent_1600 <- read.csv(file.path("Data", "ThreatsPresent_1600.csv"))
# ThreatsPresent_1600_ex <- read.csv(file.path("Data", "ThreatsPresent_1600_ex.csv"))
# ThreatsPresent_1800 <- read.csv(file.path("Data", "ThreatsPresent_1800.csv"))
# ThreatsPresent <- rbind(ThreatsPresent_200, ThreatsPresent_400, ThreatsPresent_600,
#                         ThreatsPresent_800, ThreatsPresent_1000, ThreatsPresent_1200,
#                         ThreatsPresent_1400, ThreatsPresent_1600, ThreatsPresent_1600_ex,
#                         ThreatsPresent_1800)
# write.csv(ThreatsPresent, file = file.path("Data", "ThreatsPresent.csv"))


# The species pages identified as having a Threats section will then be mined for keywords (another function)
#ThreatsPresent <- read.csv(file.path("Data", "ThreatsPresent.csv"))
#Threats_Y <- ThreatsPresent %>% filter(ThreatPresent == 1)
#ThreatsInfo <- Threat_info(SNES_min$Profile,SNES_min$scientificName, target_folder = "websites")
#write.csv(ThreatsInfo, file = file.path("Data", "ThreatsInfo.csv"))

# Load and join (rbind) different ThreatsInfo_ csv's. 
# ThreatsInfo_200 <- read.csv(file.path("Data", "ThreatsInfo_200.csv"))
# ThreatsInfo_400 <- read.csv(file.path("Data", "ThreatsInfo_400.csv"))
# ThreatsInfo_600 <- read.csv(file.path("Data", "ThreatsInfo_600.csv"))
# ThreatsInfo_800 <- read.csv(file.path("Data", "ThreatsInfo_800.csv"))
# ThreatsInfo_1000 <- read.csv(file.path("Data", "ThreatsInfo_1000.csv"))
# ThreatsInfo_1200 <- read.csv(file.path("Data", "ThreatsInfo_1200.csv"))
# ThreatsInfo_1400 <- read.csv(file.path("Data", "ThreatsInfo_1400.csv"))
# ThreatsInfo_1600 <- read.csv(file.path("Data", "ThreatsInfo_1600.csv"))
# ThreatsInfo_1800 <- read.csv(file.path("Data", "ThreatsInfo_1800.csv"))
# ThreatsInfo <- rbind(ThreatsInfo_200, ThreatsInfo_400, ThreatsInfo_600,
#                      ThreatsInfo_800, ThreatsInfo_1000, ThreatsInfo_1200,
#                      ThreatsInfo_1400, ThreatsInfo_1600, ThreatsInfo_1800)
# write.csv(ThreatsInfo, file = file.path("Data", "ThreatsInfo.csv"))
ThreatsInfo <- read.csv(file = file.path("Data", "ThreatsInfo.csv"))

#Try to combine the next few things into one function
#Keywords
IAS_keywords <- c("feral", "exotic", "invasive", "non-native")

# Determining if SPRAT information contains information on IAS
# This technically works.
species_keywords <- sapply(IAS_keywords[1:length(IAS_keywords)], function(i) {
  
  if_else(str_detect(ThreatsInfo$ThreatInfo, i) == T, 1, 0)
  
})

IAS_threats <- if_else(rowSums(species_keywords) > 0, 1, 0)

ThreatsInfo <- ThreatsInfo %>% 
  mutate(IAS_threats = IAS_threats) %>% 
  dplyr::select(scientificNames, IAS_threats) %>%
  rename(scientificName = "scientificNames")

# Join ThreatsInfo to SNES_min to get "class" information for the species threatened by IAS
SNES_min <- left_join(SNES_min, ThreatsInfo, by = "scientificName")

#Combining spatial and ancillary information 
SNES_shp <- merge(ThreatsInfo, SNES_shp_min, by = "scientificName")
SNES_shp <- st_as_sf(SNES_shp) #now ready for further spatial analysis
SNES_shp$scientificName <- gsub(" |/", "_", SNES_shp$scientificName)

#Re-project spatial layers
SNES_shp_proj <- st_transform(SNES_shp, 3577)

#Intersection with coastal layer
Aus_SNES <- st_intersection(st_buffer(Aus_Coast_proj, dist = 0), st_buffer(SNES_shp_proj, dist = 0))

#Subset of IAS threatened species
Aus_SNES_IAS <- Aus_SNES %>% filter(IAS_threats == 1)

#Remove unecessary files
rm(SNES_shp_min)
rm(SNES_shp_proj)
rm(SNES_shp)

########################################### Figures #########################################
#Should manually place the legend where a box would be if there were more groups
All_SNES_species <- ggplot(SNES_threat_sum, aes(x = threatStatus, y = n, fill = threatStatus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~group) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  #scale_fill_manual(values = mycols) +
  scale_y_continuous(trans = "log1p") +
  ylab("Number of species (log+1)") +
  xlab("Federal threat status")