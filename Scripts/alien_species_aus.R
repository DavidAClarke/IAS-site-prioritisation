#Alien species in Australia
aliens_Aus <- read.csv("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Data general/hseebens-sTwist_Workflow-7b32cf8/Output/SInAS_AlienSpeciesDB_2.3.1.csv", sep = " ")
aliens_Aus_taxalist <- read.csv("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Data general/hseebens-sTwist_Workflow-7b32cf8/Output/SInAS_AlienSpeciesDB_2.3.1_FullTaxaList.csv", sep = " ")

#Unique names
#alien_names <- unique(aliens_Aus_taxalist$species)

#Further harmonization (does this get me the same as aliens_Aus_taxalist$species?)
# AlienTaxInfo <- traitdataform::get_gbif_taxonomy(alien_names) #function from traitdataform package
# write.csv(AlienTaxInfo, file = file.path("SpeciesData", "AlienTaxInfo.csv"))
AlienTaxInfo <- read.csv(file.path("SpeciesData", "AlienTaxInfo.csv"))
AlienTaxInfo <- as_tibble(AlienTaxInfo)
AlienTaxInfo_names <- AlienTaxInfo %>% 
  distinct(verbatimScientificName = verbatimScientificName, .keep_all = T) %>%
  rename(acceptedName = "scientificName") %>%
  rename(Taxon = "verbatimScientificName") %>%
  dplyr::select(Taxon, acceptedName)
aliens_Aus <- aliens_Aus %>%
  left_join(AlienTaxInfo_names, by = "Taxon")

#Match the names
# acc_alien_names <- unique(aliens_Aus$acceptedName)
# gbif_taxon_keys <- acc_alien_names %>% #speciesNames = combination of harmonized range map and plant occurrence species (see RedList.R)
#   taxize::get_gbifid_(method = "backbone") %>%
#   imap(~ .x %>% mutate(original_sciname = .y)) %>%
#   bind_rows() %T>%
#   write.csv(file = file.path("SpeciesData", "alien_all_matches.csv")) %>% #can then examine the data
#   filter(matchtype == "EXACT" & status == "ACCEPTED") %>%
#   pull(usagekey)

alien_matches <- read.csv(file = file.path("SpeciesData", "alien_all_matches.csv"))
gbif_taxon_keys <- alien_matches %>% 
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% 
  pull(usagekey)

#GBIF credentials (remove credentials when publishing scripts)
user <- "davidclarke" 
pwd <- "!8208GBIF153"
email <- "david.clarke1@monash.edu"

#Request download
#For these I should probably get all occurrence records from the world, enabling the IAS SDM if required.
res <- occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("basisOfRecord", c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION', 'PRESERVED_SPECIMEN')),
  #pred("country", "AU"),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email)