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
#5. data.table

#Species list
speciesNames <- read.csv(file.path("SpeciesData", "all_accepted_names.csv"))
speciesNames <- speciesNames$x

#Path
path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "0310609-200613084148143.csv"))

#Import data (turn into function)
data_import_n_clean <- function(sppList, path) {
  
  sppList <- gsub(" ", ".", sppList)
  
  my_cols <- fread(cmd = paste("head -n 1", path))
  
  my_cols_min <- c("family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince", "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode")
  
  res.out <- lapply(sppList[1:length(sppList)], function(i){
    
    print(i)
    
    tryCatch({
      gbif_occs <- fread(cmd = paste("grep", i ,path), col.names = names(my_cols), quote = "", na.strings = c("",NA))
      if(is_empty(gbif_occs)) stop("No occurrence records")
      gbif_occs <- gbif_occs %>% dplyr::select(my_cols_min)
      
      #think its here where I include cleaning stuff
      
      # remove records without coordinates
      gbif_occs <- gbif_occs %>%
        filter(!is.na(decimalLongitude)) %>%
        filter(!is.na(decimalLatitude))
      
      #convert country code from ISO2c to ISO3c
      gbif_occs$countryCode <-  countrycode(gbif_occs$countryCode, origin =  'iso2c', destination = 'iso3c')
      
      #flag problems
      flags <- clean_coordinates(x = gbif_occs,
                                 lon = "decimalLongitude",
                                 lat = "decimalLatitude",
                                 countries = "countryCode",
                                 species = "species",
                                 tests = c("centroids", "equal","gbif", "institutions",
                                           "zeros", "duplicates"))
      
      #Exclude problematic records
      gbif_occs_cl <- gbif_occs[flags$.summary,]
      
      #The flagged records
      gbif_occs_fl <- gbif_occs[!flags$.summary,]
      
      #Excluding records with more than 500m uncertainty
      gbif_occs_cl <- gbif_occs_cl %>%
        filter(coordinateUncertaintyInMeters <= 500 | is.na(coordinateUncertaintyInMeters))
      
      #write csv of cleaned records
      write.csv(gbif_occs_cl, file = file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "Cleaned", paste0(i,"_cl.csv")))
      
      #write csv of flagged records
      write.csv(gbif_occs_fl, file = file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "Flagged", paste0(i,"_fl.csv")))
      
    }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
    
  })
  
  #res.out <- do.call(rbind, res.out)
  #return(res.out)
}


############################################################################################################
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
Ap_info <- data.frame(scientificName = "Acanthiza pusilla", 
                      occs_in_range = nrow(within_range),
                      prop_in_range = prop_within)
