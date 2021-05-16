####################################### Red List Ecosystems #####################################

############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "RedListEco.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"),
              local_path = file.path("Scripts", "RedListEco.R"))

################################# Load spatial data #############################################
# Load IBRA (sub)regions
IBRA_subreg.shp <- st_read("SpatialData/Vector/IBRA7_subregions/ibra7_subregions.shp")


####################################### Create Ecosystem layers #################################
#Alpine snow patch herbfields (Endangered)
ASPH <- IBRA_subreg.shp %>% 
  filter(SUB_NAME_7 == "Snowy Mountains" | SUB_NAME_7 == "Victorian Alps") %>%
  mutate(Ecosystem = "Alpine snow patch herbfields") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Endangered")

#Eastern Stirling Range Montane Heath and Thicket (Critically Endangered)
ESRMHT <- IBRA_subreg.shp %>% 
  filter(SUB_NAME_7 == "Fitzgerald") %>%
  mutate(Ecosystem = "Eastern Stirling Range Montane Heath and Thicket") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Critically Endangered")

#Mountain ash forest (Critically Endangered)
MAF <- IBRA_subreg.shp %>% 
  filter(SUB_NAME_7 == "Highlands-Southern Fall") %>%
  mutate(Ecosystem = "Mountain ash forest") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Critically Endangered")

#Busselton ironstone shrublands (Critically Endangered)
BIS <- IBRA_subreg.shp %>% 
  filter(SUB_NAME_7 == "Perth") %>%
  mutate(Ecosystem = "Busselton ironstone shrublands") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Critically Endangered")

#Coastal lowland rainforest (Endangered)
CLR <- IBRA_subreg.shp %>% 
  filter(SUB_NAME_7 == "Herbert" | SUB_NAME_7 == "Tully" |  SUB_NAME_7 == "Innisfail" | SUB_NAME_7 == "Macalister" |
         SUB_NAME_7 == "Daintree-Bloomfield") %>%
  mutate(Ecosystem = "Coastal lowland rainforest") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Endangered")

#Lake Eyre Basin (Least Concern)
LEB <- IBRA_subreg.shp %>% 
  filter(SUB_NAME_7 == "Claude River Downs" |
         SUB_NAME_7 == "Carnarvon Ranges" |
         SUB_NAME_7 == "Southern Downs" |
         SUB_NAME_7 == "Barrier Range" |
         SUB_NAME_7 == "Barrier Range Outwash" |
         SUB_NAME_7 == "Bimbowrie" |
         SUB_NAME_7 == "Curnamona" |
         SUB_NAME_7 == "Atartinga" |
         SUB_NAME_7 == "Dulcie" |
         SUB_NAME_7 == "Everard Block" |
         REG_NAME_7 == "Channel Country" |
         REG_NAME_7 == "Desert Uplands" |
         REG_NAME_7 == "Finke" |
         SUB_NAME_7 == "Broughton" |
         SUB_NAME_7 == "Olary Spur" |
         SUB_NAME_7 == "Southern Flinders" |
         SUB_NAME_7 == "Northern Flinders" |
         SUB_NAME_7 == "Central Flinders" |
         SUB_NAME_7 == "Torrens" |
         REG_NAME_7 == "MacDonnell Ranges" |
         REG_NAME_7 == "Mitchell Grass Downs" |
         REG_NAME_7 == "Mount Isa Inlier" |
         SUB_NAME_7 == "West Warrego" |
         SUB_NAME_7 == "Northern Uplands" |
         SUB_NAME_7 == "West Bulloo" |
         REG_NAME_7 == "Simpson Strzelecki Dunefields" |
         REG_NAME_7 == "Stony Plains" |
         SUB_NAME_7 == "Sandover") %>%
  mutate(Ecosystem = "Lake Eyre Basin") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Least Concern")

#Cumberland Plain Woodland (Critically Endangered)
CPW <- IBRA_subreg.shp %>% 
  filter(SUB_NAME_7 == "Yengo" | SUB_NAME_7 == "Cumberland" |  SUB_NAME_7 == "Burragorang" | SUB_NAME_7 == "Sydney Cataract" |
         SUB_NAME_7 == "Pittwater") %>%
  mutate(Ecosystem = "Cumberland Plain Woodland") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Critically Endangered")

#Georgina gidgee woodlands (Vulnerable)
GGW <- IBRA_subreg.shp %>% filter(SUB_NAME_7 == "Yuendumu" | SUB_NAME_7 == "Atartinga" |
                                  SUB_NAME_7 == "Mount Chapple" | SUB_NAME_7 == "Dulcie" |
                                    SUB_NAME_7 == "Toko Plains" | SUB_NAME_7 == "Sturt Stony Desert" |
                                    SUB_NAME_7 == "Diamantina-Eyre" |
                                    SUB_NAME_7 == "Cooper-Diamantina Plains" |
                                    SUB_NAME_7 == "Ashburton Range" | SUB_NAME_7 == "Davenport" |
                                    SUB_NAME_7 == "Barkly" | SUB_NAME_7 == "Henbury" |
                                    SUB_NAME_7 == "Finke River" | SUB_NAME_7 == "Tieyon" |
                                    SUB_NAME_7 == "Pedirka" | SUB_NAME_7 == "Mackay" |
                                    SUB_NAME_7 == "Ehrenberg" | SUB_NAME_7 == "Amedeus" |
                                    SUB_NAME_7 == "Lake Bennett" | SUB_NAME_7 == "Lake Lewis" |
                                    SUB_NAME_7 == "McDonnell" | SUB_NAME_7 == "Watarrka" |
                                    SUB_NAME_7 == "Hartz Range" | SUB_NAME_7 == "Mitchell Grass Downs P1" |
                                    SUB_NAME_7 == "Barkly Tableland" | SUB_NAME_7 == "Georgina Limestone" |
                                    SUB_NAME_7 == "Southwestern Plateaus and Floodouts" |
                                    SUB_NAME_7 == "Andado" | SUB_NAME_7 == "Simpson Desert" |
                                    SUB_NAME_7 == "Dieri" | SUB_NAME_7 == "Breakaways" |
                                    SUB_NAME_7 == "Oodnadatta" | SUB_NAME_7 == "Macumba" |
                                    SUB_NAME_7 == "Witjira" | SUB_NAME_7 == "Wycliffe" |
                                    SUB_NAME_7 == "Sandover") %>%
  mutate(Ecosystem = "Georgina gidgee woodlands") %>%
  dplyr::select(Ecosystem, SUB_NAME_7, REG_NAME_7, SQ_KM, geometry) %>%
  mutate(Status = "Vulnerable")

#combine ecosystems?
Ecosystems <- rbind(ASPH, BIS, CLR, CPW, ESRMHT, GGW, LEB, MAF)

#################################################### Spatial operations ###########################################################
#Transform
IBRA_proj <- st_transform(IBRA_subreg.shp, 3577)

#Intersect with coast
IBRA_proj <- st_intersection(st_buffer(IBRA_proj, dist = 0), st_buffer(Aus_Coast_proj, dist = 0))

#Change extent (?)

#Add Red list status attribute