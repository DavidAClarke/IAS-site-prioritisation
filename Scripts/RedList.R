############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "RedList.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "RedList.R"),
              local_path = file.path("Scripts", "RedList.R"), overwrite = T)

############################### Load Red List spatial and assessment data #########################

# For downloading on virtual machine
# RL_files <- c("data_0.cpg", "data_0.dbf", "data_0.prj", "data_0.shp", "data_0.shx", "data_1.cpg", "data_1.dbf",
#  "data_1.prj", "data_1.shp", "data_1.shx", "data_2.cpg", "data_2.dbf", "data_2.prj", "data_2.shp",
#  "data_2.shx")
# 
# # Download each file
# lapply(RL_files[1:length(RL_files)], function(i) {
#   
#   drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "SpatialData", "Vector", "redlist_species_data_18052021", i),
#                 local_path = file.path("SpatialData", "Vector", "redlist_species_data_18052021", i))
#   
# })

#Load spatial data
#RL_Aus_shp <- st_read("SpatialData/Vector/redlist_species_data_Aus/data_0.shp")

#Updated downloaded shapefiles. Do I combine them all? How do they differ?
RL_shp_0 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_0.shp")
RL_shp_1 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_1.shp")
RL_shp_2 <- st_read("SpatialData/Vector/redlist_species_data_18052021/data_2.shp")
RL_shp <- rbind(RL_shp_0,RL_shp_1,RL_shp_2)
rm(RL_shp_0,RL_shp_1,RL_shp_2)

#Save new combined shapefile
st_write(RL_shp, 
         dsn = file.path("SpatialData", "Vector", "redlist_species_data_18052021"),
         layer = "RL_shp.shp",
         driver = "ESRI Shapefile")

#Load plant point data
RL_plants <- read.csv(file.path("SpatialData", "Vector", "redlist_species_data_plantpoints", "points_data.csv"))
RL_plants <- RL_plants %>%
  dplyr::select("binomial", "presence", "origin", "seasonal", "legend", "longitude", "latitude")
RL_plants_points <- RL_plants %>%
  dplyr::select("longitude", "latitude")
RL_plants_points.sp <- SpatialPointsDataFrame(RL_plants_points,
                                              data = RL_plants,
                                              proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
RL_plants_points_sf <- st_as_sf(RL_plants_points.sp)
rm(RL_plants_points.sp,RL_plants_points,RL_plants)



#Also need to load Red List assessment and threats data. Join to RL spatial layer #~#
#RL_Aus_assess <- read_excel(file.path("SpeciesData", "redlist_species_data", "assessments.xlsx"), sheet = "assessments")
#RL_Aus_threats <- read_excel(file.path("SpeciesData", "redlist_species_data", "threats.xlsx"), sheet = "threats")

#More recent download - broader scope
RL_Aus_assess <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "assessments.csv"))
RL_Aus_threats <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "threats.csv"))
RL_Aus_taxon <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "taxonomy.csv"))

####################################### Data Preparation ####################################
#Joining taxonomic information to red list assessment data
RL_info <- RL_Aus_assess %>% 
  left_join(RL_Aus_threats, by = "scientificName") %>%
  left_join(RL_Aus_taxon, by = "scientificName") %>%
  dplyr::select(scientificName, orderName, className, phylumName, kingdomName, 
                redlistCategory, redlistCriteria, yearPublished,code, name, 
                stressCode, stressName, ias, scope) %>%
  mutate(className = as.factor(className)) %>%
  mutate(orderName = as.factor(orderName)) %>%
  mutate(phylumName = as.factor(phylumName)) %>%
  mutate(kingdomName = as.factor(kingdomName)) %>%
  filter(redlistCategory != "Lower Risk/near threatened", 
         redlistCategory != "Lower Risk/conservation dependent",
         redlistCategory != "Lower Risk/least concern") %>%
  mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern", "Near Threatened",
                                             "Vulnerable", "Endangered", "Critically Endangered", "Extinct"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", "CR", "EX")))
rm(RL_Aus_assess, RL_Aus_taxon, RL_Aus_threats)

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
plant_mcp <- mcp(RL_plants_points.sp[1:60,1], percent = 95)

#Might be good to overlay this with something like GADM to get country information
#That would be necessary to calculate proportion of Australia to overall range
#Combining spatial and ancillary information 
RL_Aus_shp <- RL_shp_0 %>%
  dplyr::select(BINOMIAL, PRESENCE, ORIGIN, SEASONAL, LEGEND) %>%
  rename(scientificName = "BINOMIAL")
new_bb <- c(113, -44, 154, -9)
names(new_bb) <- c("xmin", "ymin", "xmax", "ymax")
RL_Aus_shp <- st_crop(RL_Aus_shp, new_bb) #bringing extent down before projection
RL_Aus_shp <- merge(RL_info, RL_Aus_shp, by = "scientificName")
RL_Aus_shp <- st_as_sf(RL_Aus_shp) #now ready for further spatial analysis
RL_Aus_shp_IAS <- RL_Aus_shp %>% filter(code == "8.1.1" | code == "8.1.2" | code == "8.2.2")

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
