############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "occurrence_records.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "occurrence_records.R"),
              local_path = file.path("Scripts", "occurrence_records.R"), overwrite = T)

######################################## Downloading occurrence records ###########################
#This should contain both GBIF and ALA

#Downloading ALA records
library(ALA4R)

#Species list
speciesNames <- unique(c(RL_shp_Aus$acceptedName, RL_plants$acceptedName))

#Following should occur for every species
#Example: Ozimops planiceps
Op <- occurrences(taxon = "Ozimops planiceps",  #this will be "taxon = i" in loop
                  fields = c("family",
                             "species",
                             "rank",
                             "accepted_name_usage",
                             "longitude",
                             "latitude",
                             "country_code",
                             "state",
                             "coordinate_uncertainty",
                             "coordinate_precision",
                             "basis_of_record",
                             "year",
                             "individual_count",
                             "institution_code"),
                  email = "david.clarke1@monash.edu",
                  qa = "all")

#Just get the data
Op <- Op$data[,1:14]

#Rename point record fields to match GBIF
Op <- Op %>% 
  rename(decimalLongitude = "longitude") %>%
  rename(decimalLatitude = "latitude") %>%
  rename(coordinateUncertaintyInMeters = "coordinateUncertaintyInMetres") %>%
  mutate(countryCode = "AU") %>%
  mutate(basisOfRecord = factor(basisOfRecord)) %>%
  mutate(basisOfRecord = recode(basisOfRecord,
                                HumanObservation = "HUMAN_OBSERVATION",
                                MachineObservation = "MACHINE_OBSERVATION",
                                ))
