                                          #Maxent via ENMTools
#######################################################################################################################
options(java.parameters = "-Xmx8g" )
library(tidyverse)
library(sf)
library(ENMTools)

#Add paths to necessary data
env_path <- "E:/Maxent/Env_layers"
occ_path <- "E:/Maxent/Occurrences"
bias_path <- "E:/Maxent/Bias_layers"
zonation_data_path <- "E:/Zonation/Input_data/maxent_output"

source("Scripts/2.aus_coast.R")

#Get species list to loop over (speciesNames_range?) (using plants for testing)
RL_plants <- read.csv("E:/SpatialData/Vector/redlist_species_data_plantpoints/RL_plants_points.csv")
RL_plants <- subset(RL_plants, with(RL_plants, unsplit(table(acceptedName), acceptedName)) >= 5)
speciesNames_points <- unique(RL_plants$acceptedName)
rm(RL_plants)

#Species list will ultimately just be all species

#Removing points with no environmental data
rm_occs <- function(env, pts) {
  
  oc <- extract(env, pts[,2:3])
  oc <- as.data.frame(oc)
  pts <- cbind(pts, oc)
  pts <- na.omit(pts)
  pts <- pts %>% dplyr::select(scientificName, Longitude, Latitude)
  return(pts)
  
}

#Create empty vectors for model evaluations
mod_evals <- data.frame(Species = character(), 
                        n = numeric(), 
                        np = numeric(), 
                        auc = numeric(), 
                        ECE = numeric(),
                        MCE = numeric(),
                        Boyce = numeric())
write.csv(mod_evals, file = "E:/Maxent/mod_evals.csv")

start_time <- Sys.time()
res.out <- lapply(speciesNames_points[1:length(speciesNames_points)], function(i){
  
  print(i)
  i <- gsub(" ", "_", i)
  
  tryCatch({
  #First look for occurrence records
  occs_files <- list.files(file.path(occ_path, i), pattern = ".csv")
  if(is_empty(occs_files)) stop("No occurrence records")
  occs_cl <- read.csv(file.path(occ_path, i, paste0(i, "_cl.csv")))
  
  if(nrow(occs_cl) < 15) stop("Less than 15 cleaned/non-duplicated occurrence records")
  occs_cl <- occs_cl %>%
    dplyr::select(scientificName, decimalLongitude,decimalLatitude) %>%
    dplyr::rename(Longitude = "decimalLongitude") %>%
    dplyr::rename(Latitude = "decimalLatitude")
  
  #Look for bias file
  bias_files <- list.files(file.path(bias_path, i), pattern = ".gri")
  if(is_empty(bias_files)) stop("No bias file")
  bias <- raster(file.path(bias_path, i, paste0(i, "_bias.gri")))
  
  #Read in environmental data
  Env_range <- stack(file.path(env_path, i, paste0(i, "_env.gri")))
  #Env and bias should be the same extent, resolution, etc. 
  
  #Remove occurrences that aren't within Env_range
  occs_cl <- rm_occs(Env_range, occs_cl)
  
  #Create enmtools.species object
  species <- enmtools.species(species.name = i, 
                              presence.points = occs_cl) #must be named Longitude and Latitude
  
  #Run maxent model
  maxent_mod <- enmtools.maxent(species = species,
                               env = Env_range[[1:11]], #for now just use the first 11
                               nback = 2000,
                               test.prop = 0.2, 
                               bias = bias,
                               args = "outputformat=raw")
  
  #Model calibration
  mod_cal <- enmtools.calibrate(maxent_mod)
  
  #Model stats
  n <- nrow(occs_cl)
  np <- maxent_mod$test.evaluation@np
  auc <- maxent_mod$test.evaluation@auc
  ECE <- mod_cal$ECE
  MCE <- mod_cal$MCE
  Boyce <- mod_cal$continuous.boyce$Spearman.cor
  evals <- cbind(i, n, np, auc, ECE, MCE, Boyce)
  write.table(evals, file = "E:/Maxent/mod_evals.csv", append = T, sep = ",", col.names = F)
  
  #Write predicted distribution layer for input into zonation
  suit_layer <- maxent_mod$suitability
  suit_layer <- extend(suit_layer, Aus_Coast)
  writeRaster(suit_layer, filename = file.path(zonation_data_path, paste0(i, "_pred.tif")), overwrite = T)
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
})
end_time <- Sys.time()
end_time - start_time