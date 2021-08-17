#Maxent via ENMTools
#######################################################################################################################
options(java.parameters = "-Xmx16g" ) 
library(tidyverse)
library(sf)
library(ENMTools)
library(rJava)

#Load data
load("/projects/nc57/Chapter_3/SpatialData/RL_shp_prepro.RData")
rm(RL_shp_Aus)

#Species
speciesNames_range <- speciesNames_range[251:500] #match respective offset group 

#Add paths to necessary data
env_path <- "/home/dcla0008/nc57_scratch/Maxent/Env_layers"
bias_path <- "/home/dcla0008/nc57_scratch/Maxent/Bias_layers"
occ_path <- "/home/dcla0008/nc57_scratch/Maxent/Occurrences"
zonation_data_path <- "/home/dcla0008/nc57_scratch/Zonation/Input_data/maxent_output"

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
write.csv(mod_evals, file = "/home/dcla0008/nc57_scratch/Maxent/mod_evals2.csv")


res.out <- lapply(speciesNames_range[1:length(speciesNames_range)], function(i){
  
  print(i)
  i <- gsub(" ", "_", i)
  
  tryCatch({
    #First look for occurrence records
    occs_files <- list.files(paste0(occ_path,"/", i), pattern = ".csv")
    if(is_empty(occs_files)) stop("No occurrence records")
    occs_cl <- read.csv(paste0(occ_path,"/", i, "/", i, "_cl.csv"))
    
    if(nrow(occs_cl) < 15) stop("Less than 15 cleaned/non-duplicated occurrence records")
    occs_cl <- occs_cl %>%
      dplyr::select(scientificName, decimalLongitude,decimalLatitude) %>%
      dplyr::rename(Longitude = "decimalLongitude") %>%
      dplyr::rename(Latitude = "decimalLatitude")
    
    #Look for bias file
    bias_files <- list.files(paste0(bias_path,"/", i), pattern = ".gri")
    if(is_empty(bias_files)) stop("No bias file")
    bias <- raster(paste0(bias_path,"/", i, "/",i, "_bias.gri"))
    
    #Read in environmental data
    Env_range <- stack(paste0(env_path,"/", i, "/",i, "_env.gri"))
    Env_range <- check.env(Env_range)
    #Env and bias should be the same extent, resolution, etc. 
    
    #Remove occurrences that aren't within Env_range
    occs_cl <- rm_occs(Env_range, occs_cl)
    
    #Create enmtools.species object
    species <- enmtools.species(species.name = i, 
                                presence.points = occs_cl) #must be named Longitude and Latitude
    
    #Run maxent model
    maxent_mod <- enmtools.maxent(species = species,
                                  env = Env_range, #for now just use the first 11
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
    write.table(evals, file = "/home/dcla0008/nc57_scratch/Maxent/mod_evals2.csv", append = T, sep = ",", col.names = F)
    
    #Write predicted distribution layer for input into zonation
    suit_layer <- maxent_mod$suitability
    suit_layer_ext <- extend(suit_layer, Aus_Coast)
    writeRaster(suit_layer_ext, filename = paste0(zonation_data_path, "/",i, "_pred.tif"), overwrite = T)
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
})
