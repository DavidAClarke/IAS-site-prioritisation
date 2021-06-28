############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "occurrence_records.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "occurrence_records.R"),
              local_path = file.path("Scripts", "occurrence_records.R"), overwrite = T)

################################ Downloading/importing occurrence records ###########################
                                         #Both GBIF and ALA

############################################# GBIF #########################################
# library(rgbif)

#GBIF credentials (remove credentials when publishing scripts)
user <- "davidclarke" 
pwd <- "!8208GBIF153"
email <- "david.clarke1@monash.edu"

#Get accepted names for use in download
speciesNames <- unique(c(RL_shp_Aus$acceptedName, RL_plants$acceptedName))
write.csv(speciesNames, file = file.path("SpeciesData", "all_accepted_names.csv"))
speciesNames <- read.csv(file.path("SpeciesData", "all_accepted_names.csv"))

#Match the names
# gbif_taxon_keys <- speciesNames %>% #speciesNames = combination of harmonized range map and plant occurrence species (see RedList.R)
#   taxize::get_gbifid_(method = "backbone") %>%
#   imap(~ .x %>% mutate(original_sciname = .y)) %>%
#   bind_rows() %T>%
#   write.csv(file = file.path("SpeciesData", "all_matches.csv")) %>% #can then examine the data
#   filter(matchtype == "EXACT" & status == "ACCEPTED") %>%
#   pull(usagekey)
# matches <- read.csv(file.path("SpeciesData", "all_matches.csv"))
# setdiff(speciesNames, matches$canonicalname) #matches (therefore GBIF) doesn't contain "Ozimops planiceps". 

#Request download
# res <- occ_download(
#   pred_in("taxonKey", gbif_taxon_keys),
#   pred_in("basisOfRecord", c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION', 'PRESERVED_SPECIMEN')),
#   pred("country", "AU"),
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SIMPLE_CSV",
#   user=user,pwd=pwd,email=email)

#Download data
#occ_download_get("0299151-200613084148143", path = file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF"))
occ_download_get("0310609-200613084148143", path = file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF"))


#first get column names
# my_cols <- fread(file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "0310609-200613084148143.csv"), nrows = 0)
# my_cols_min <- c("family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince",
#                  "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode")
# path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "0310609-200613084148143.csv"))
# gbif_occs <- fread(cmd = paste("grep", "Acanthiza.pusilla", path), col.names = names(my_cols), quote = "", na.strings = c("",NA)) #works!
# gbif_occs_min <- gbif_occs %>% dplyr::select(my_cols_min)



################################################ ALA #########################################
#Downloading ALA records
library(galah)

#Setting configuration
ala_config(cache_directory = file.path("SpatialData", "Vector", "SpeciesOccurrences", "ALA"), 
           email = "david.clarke1@monash.edu", 
           verbose = T)

#Species list
speciesList <- speciesNames$x #if speciesNames is read in
 
# #Obtaining number of ALA occurrence records (change to galah functionality)
get_ALA_occ_num <- function(species_names) {

  res.out <- lapply(species_names[1:length(species_names)], function(i){

    occ_num <- ala_counts(taxa = i, type = "record")

    occ_num <- data.frame(species = i, occ_num = occ_num)

  })

  res.out <- do.call(rbind, res.out)
  return(res.out)
}

#Get number of ALA occurrences per species
# ALA_occ_nums <- get_ALA_occ_num(speciesNames)
# ALA_occ_nums <- ALA_occ_nums %>%
#   filter(occ_num > 0)
# speciesNames_upd <- unique(ALA_occ_nums$species)
# write.csv(speciesNames_upd, file = file.path("SpatialData", "Vector", "SpeciesOccurrences", "ALA", "speciesNames_upd.csv"))
speciesNames_upd <- read.csv(file.path("SpatialData", "Vector", "SpeciesOccurrences", "ALA", "speciesNames_upd.csv"))
speciesNames_upd <- as.character(speciesNames_upd$x)

#Getting ALA occurrences (change to galah functionality)
get_ALA_occs <- function(species_names){
  
  #make empty data.frame
  df <- data.frame(family = character(), species = character(), taxonRank = character(), scientificName = character(),
                   decimalLongitude = numeric(), decimalLatitude = numeric(), countryCode = character(),
                   stateProvince = character(), coordinateUncertaintyInMeters = numeric(),
                   coordinatePrecision = numeric(), basisOfRecord = character(), year = integer(), individualCount = numeric(),
                   institutionCode = character())
  
  #write empty df to csv
  write.csv(df, file = file.path("SpatialData", "Vector", "SpeciesOccurrences", "ALA", "ALA_occs.csv"))
  
  res.out <- lapply(species_names[1:length(species_names)], function(i){
    
    #Only really need this chunk. Just input a vector of names
    occs <- ala_occurrences(taxa = i,  #this will be "taxon = i" in loop
                      filters = select_filters(
                        basisOfRecord = c("HumanObservation", "MachineObservation", "PreservedSpecimen")),
                      columns = select_columns(
                        "family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince",
                        "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode"),
                      mint_doi = T
                      )
                        
    occs <- occs %>%
      mutate(countryCode = "AU") %>%
      mutate(basisOfRecord = factor(basisOfRecord)) %>%
      mutate(basisOfRecord = recode(basisOfRecord,
                                    HumanObservation = "HUMAN_OBSERVATION",
                                    MachineObservation = "MACHINE_OBSERVATION",
                                     PreservedSpecimen = "PRESERVED_SPECIMEN"))
    
    #append occs to csv (using write.table)
    write.table(occs, file = file.path("SpatialData", "Vector", "SpeciesOccurrences", "ALA", "ALA_occs.csv"), append = T,sep = ",",
                col.names = F)
    })
  
  res.out <- do.call(rbind,res.out)
  return(res.out)
  
}

ALA_occs <- get_ALA_occs(speciesNames_upd)

