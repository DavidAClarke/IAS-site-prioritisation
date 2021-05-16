############################### Load Red List spatial and assessment data #########################

#Load spatial data
RL_Aus_shp <- st_read("SpatialData/Vector/redlist_species_data_Aus/data_0.shp")

#Also need to load Red List assessment and threats data. Join to RL spatial layer #~#
RL_Aus_assess <- read_excel(file.path("Data", "redlist_species_data", "assessments.xlsx"), sheet = "assessments")
RL_Aus_threats <- read_excel(file.path("Data", "redlist_species_data", "threats.xlsx"), sheet = "threats")


#################################### Preparing assessment data #############################
#Getting higher level taxonomic information
#Information for some species was not retrieved during function.
scinames <- unique(RL_Aus_assess$scientificName)
RL_class <- GetTax(scinames) #Custom function
RL_class <- RL_class %>%
  mutate(canonicalName = as.character(canonicalName)) %>%
  dplyr::rename(scientificName = "canonicalName") %>%
  mutate(class = stringr::str_replace(class, "Ancylastrum cumingianus", "Gastropoda")) %>%
  mutate(order = stringr::str_replace(order, "Gastropoda", "Lymnaeida")) %>%
  mutate(class = stringr::str_replace(class, "Notomys robustus", "Mammalia")) %>%
  mutate(class = stringr::str_replace(class, "Threskiornis moluccus", "Aves")) %>%
  mutate(order = stringr::str_replace(order, "Aves", "Pelecaniformes")) %>%
  mutate(class = stringr::str_replace(class, "Ctenophorus graafi", "Reptilia")) %>%
  mutate(order = stringr::str_replace(order, "Ctenophorus graafi", "Squamata")) %>%
  mutate(class = stringr::str_replace(class, "Ctenophorus infans", "Reptilia")) %>%
  mutate(order = stringr::str_replace(order, "Ctenophorus infans", "Squamata")) %>%
  mutate(class = stringr::str_replace(class, "Ctenophorus slateri", "Reptilia")) %>%
  mutate(order = stringr::str_replace(order, "Ctenophorus slateri", "Squamata")) %>%
  mutate(class = stringr::str_replace(class, "Nyctimystes infrafrenatus", "Amphibia")) %>%
  mutate(order = stringr::str_replace(order, "Nyctimystes infrafrenatus", "Anura")) %>%
  mutate(order = stringr::str_replace(order, "Gelochelidon macrotarsa", "Charadriiformes")) %>%
  mutate(class = stringr::str_replace(class, "Gelochelidon macrotarsa", "Aves")) %>%
  mutate(class = stringr::str_replace(class, "Herennia oz", "Arachnida")) %>%
  mutate(class = stringr::str_replace(class, "Nephila edulis", "Arachnida")) %>%
  mutate(class = stringr::str_replace(class, "Elasmobranchii", "Chondrichthyes"))

#Joining taxonomic information to red list assessment data
RL_Aus_assess <- left_join(RL_Aus_assess, RL_class, by = "scientificName")

#Preparation of assessment data
RL_Aus_assess_min <- RL_Aus_assess %>%
  dplyr::select(scientificName, order, class, redlistCategory, redlistCriteria, yearPublished) %>%
  mutate(class = as.factor(class)) %>%
  mutate(order = as.factor(order))%>%
  filter(redlistCategory != "Extinct", redlistCategory != "Lower Risk/conservation dependent") %>%
  mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern", "Near Threatened",
                                             "Vulnerable", "Endangered", "Critically Endangered"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", "CR")))

#Summary information for Red List categories
RL_Aus_assess_sum <- RL_Aus_assess_min %>%
  group_by(class, redlistCategory) %>%
  count(redlistCategory)

#Joining threat data with taxonomic information
#Note: threat data is smaller than assessment data
RL_Aus_threats <- left_join(RL_Aus_threats, RL_class, by = "scientificName")

#Subsetting species with IAS as known threat
RL_Aus_IASthreats <- RL_Aus_threats %>%
  dplyr::select(scientificName, order, class, code, name, stressCode, stressName, ias, scope) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(scientificName = as.character(scientificName)) %>%
  mutate(ias = as.character(ias)) %>%
  filter(code == "8.1.1" | code == "8.1.2" | code == "8.2.2")

RL_Aus_IASthreats_distinct <- RL_Aus_IASthreats %>%
  distinct(scientificName = scientificName, .keep_all = T)

scinames_IAS <- unique(RL_Aus_IASthreats$scientificName)

#Summary IAS threat information
RL_Aus_IASthreats_sum <- RL_Aus_IASthreats %>%
  group_by(class, name) %>%
  count(name)

#Joining assessment and IAS threat data (loses 13 extinct species)
RL_Aus_join_IAS <- RL_Aus_assess_min %>% 
               left_join(RL_Aus_IASthreats_distinct, by = c("scientificName", "order", "class")) %>%
               filter(scientificName %in% scinames_IAS)

##Summary information for Red List categories threatened by IAS
RL_Aus_join_IAS_sum <- RL_Aus_join_IAS %>%
  group_by(class, redlistCategory) %>%
  count(redlistCategory)

###################################### Spatial operations ##################################

#Combining spatial and ancillary information 
RL_Aus_shp <- RL_Aus_shp %>%
  dplyr::select(BINOMIAL, PRESENCE, ORIGIN, SEASONAL, LEGEND) %>%
  rename(scientificName = "BINOMIAL")
RL_Aus_shp <- st_crop(RL_Aus_shp, Aus_Coast) #bringing extent down before projection
RL_Aus_shp_IAS <- merge(RL_Aus_join_IAS, RL_Aus_shp, by = "scientificName")
RL_Aus_shp_IAS <- st_as_sf(RL_Aus_shp_IAS) #now ready for further spatial analysis

#Re-project spatial layers
RL_Aus_shp_proj <- st_transform(RL_Aus_shp, 3577)
RL_Aus_shp_IAS_proj <- st_transform(RL_Aus_shp_IAS, 3577)

#Intersection with coastal layer
Aus_RL <- st_intersection(st_buffer(Aus_Coast_proj, dist = 0), st_buffer(RL_Aus_shp_proj, dist = 0))
Aus_RL_IAS <- st_intersection(st_buffer(Aus_Coast_proj, dist = 0), st_buffer(RL_Aus_shp_IAS_proj, dist = 0))

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
