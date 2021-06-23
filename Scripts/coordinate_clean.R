############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "coordinate_clean.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "coordinate_clean.R"),
              local_path = file.path("Scripts", "coordinate_clean.R"), overwrite = T)

################################# Cleaning occurrence data #################################

#Create function to do this for any species
#What is needed?
#Data
#1. Species occurrences from both GBIF and ALA (in GBIF format)
#2. species range map

#Packages
#1. coordinateCleaner
#2. sp
#3. countrycode
#4. tidyverse

# remove records without coordinates
Ap_gbif_occs <- Ap_gbif_occs %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

#convert country code from ISO2c to ISO3c
Ap_gbif_occs$countryCode <-  countrycode(Ap_gbif_occs$countryCode, origin =  'iso2c', destination = 'iso3c')

#flag problems
flags <- clean_coordinates(x = Ap_gbif_occs,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("centroids", "equal","gbif", "institutions",
                                     "zeros", "countries", "duplicates")) # most test are on by default

#Exclude problematic records
Ap_gbif_occs_cl <- Ap_gbif_occs[flags$.summary,]

#The flagged records
Ap_gbif_occs_fl <- Ap_gbif_occs[!flags$.summary,]

#Remove records with low coordinate precision
hist(Ap_gbif_occs_cl$coordinateUncertaintyInMeters / 1000, breaks = 20)
Ap_gbif_occs_cl <- Ap_gbif_occs_cl %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters)) %>% #uncertainty less than 100m
  filter(basisOfRecord == "HUMAN_OBSERVATION" |
           basisOfRecord == "OBSERVATION" |
           basisOfRecord == "MACHINE_OBSERVATION")

######################################## Separate script/function? ############################
#Should be based on the combination of GBIF and ALA records

#Determine what proportion of cleaned records occur within IUCN range map
#selecting species
Ap <- RL_shp_Aus %>% filter(BINOMIAL == "Acanthiza pusilla")

#convert to sp
Ap_sp <- as(Ap, Class = "Spatial")

#Convert occurrences to spatial layer
coordinates(Ap_gbif_occs_cl) <- ~ decimalLongitude + decimalLatitude

#Ensure the same CRS
proj4string(Ap_gbif_occs_cl) <- proj4string(Ap_sp)

#Get occurrences that occur within range map
within_range <- Ap_gbif_occs_cl[Ap_sp,]
prop_within <- nrow(within_range)/nrow(Ap_gbif_occs_cl)

#[For function?] create data frame with information
#example
Ap_info <- data.frame(c(scientificName = "Acanthiza pusilla", 
                        occs_in_range = nrow(within_range),
                        prop_in_range = prop_within)) 
