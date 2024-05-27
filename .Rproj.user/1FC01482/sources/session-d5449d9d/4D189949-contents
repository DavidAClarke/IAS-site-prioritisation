#Load zonation input species
spp_list <- read.csv(file.path("Zonation", "specs.csv"))
spp_list <- spp_list %>% rename(acceptedName = "x")

redlist_status <- read.csv(file.path("SpeciesData", 
                                     "redlist_species_data_18052021", 
                                     "assessments.csv"))
redlist_taxonomy <- read.csv(file.path("SpeciesData", 
                                       "redlist_species_data_18052021", 
                                       "taxonomy.csv"))
threats <- read.csv(file.path("SpeciesData", 
                              "redlist_species_data_18052021", "threats.csv"))
TaxInfo <- read.csv(file.path("SpeciesData", "TaxInfo_assess_names.csv"))

spp_list_info <- spp_list %>%
  left_join(TaxInfo, by = "acceptedName") %>%
  left_join(redlist_status, by= "scientificName") %>%
  left_join(redlist_taxonomy, by = "scientificName") %>%
  left_join(threats, by = "scientificName") %>%
  dplyr::distinct(acceptedName = acceptedName, .keep_all = T) %>%
  dplyr::select(acceptedName, className, redlistCategory, code, ias) %>%
  dplyr::mutate(weight = NA) %>%
  drop_na(className)

for(i in 1:length(spp_list_info$acceptedName)) {
  
  if(spp_list_info[i,]$redlistCategory == "Least Concern") {
    
    spp_list_info[i,]$weight <- 1
    
  }
  if(spp_list_info[i,]$redlistCategory == "Near Threatened") {
    
    spp_list_info[i,]$weight <- 2
    
  }
  if(spp_list_info[i,]$redlistCategory == "Vulnerable") {
    
    spp_list_info[i,]$weight <- 4
    
  }
  if(spp_list_info[i,]$redlistCategory == "Endangered") {
    
    spp_list_info[i,]$weight <- 6
    
  }
  if(spp_list_info[i,]$redlistCategory == "Critically Endangered") {
    
    spp_list_info[i,]$weight <- 8
    
  }
  if(spp_list_info[i,]$redlistCategory == "Data Deficient") {
    
    spp_list_info[i,]$weight <- 2
    
  }}

spp_list_info <- spp_list_info %>%
  drop_na(weight) #removing extinct species

for(i in 1:length(spp_list_info$acceptedName)) {
  
  if(is.na(spp_list_info[i,]$code)) 
    
    next
  
  if(spp_list_info[i,]$code == "8.1.1" || spp_list_info[i,]$code == "8.1.2"){
    
    spp_list_info[i,]$weight <- spp_list_info[i,]$weight *2
  }}

write.csv(spp_list_info,file.path("Zonation", "spp_list.csv"))
