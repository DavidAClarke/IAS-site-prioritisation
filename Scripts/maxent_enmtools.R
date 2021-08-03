                                          #Maxent via ENMTools
#Environmental data needs to be cropped to species
#use environmental data cropped during offset/bias
#example
Range <- RL_shp_Aus %>% dplyr::filter(acceptedName == "Onthophagus manya")
Range <- as(Range, Class = "Spatial")
Env <- env_predictors
Env_Range <- crop(Env, Range)

#get some points from gbif
Om_points <- occ_data(scientificName = "Onthophagus manya", country = "AU", hasCoordinate = T, hasGeospatialIssue = F)
Om_points <- Om_points$data
Om_points <- Om_points %>%
  dplyr::select(scientificName, decimalLongitude,decimalLatitude) %>%
  rename(Longitude = "decimalLongitude") %>%
  rename(Latitude = "decimalLatitude")

#Writing to folder
write.csv(Om_points, file.path("Maxent", "Occurrences", "Onthophagus_manya", "Onthophagus_manya_pres.csv"), row.names = F)

#Offset rasters are not written with crs, therefore, need to be given after reading in
Om_bias <- raster(file.path("Maxent", "Bias_layers", "Onthophagus_manya", "Onthophagus_manya_bias.asc"))
crs(Om_bias) <- crs(Env_Range) #maybe not needed if saved as .grd

#Create enmtools.species object
species <- enmtools.species(species.name = "Onthophagus manya",
                            presence.points = Om_points) #must be named Longitude and Latitude

#Run maxent model
#try adding additional arguments e.g. output type, also adding report = file.path()
start_time <- Sys.time()
Om_maxent <- enmtools.maxent(species = species,
                             env = Env_Range,
                             test.prop = 0.2, #also try with "block"
                             bias = Om_bias,
                             args = c("outputformat=raw", 
                                      "replicates=10", 
                                      "replicatetype=crossvalidate"))
end_time <- Sys.time()
end_time - start_time

env.evaluate(species, Om_maxent, Env_Range, test.eval = T)


#Add paths to necessary data
env_path <- 
occurrence_path <- 
bias_path <- 

#Get species list to loop over (speciesNames_range?)

lapply()
  
  #First look for occurrence records
  occs_files <- list.files()
  if(is_empty(occs_files)) stop("No occurrence records")
  occs_cl <- read.csv()
  if(nrows(occs_cl) < 15) stop("Less than 15 cleaned/non-duplicated occurrence records")
  occs_cl <- occs_cl %>%
    dplyr::select(scientificName, decimalLongitude,decimalLatitude) %>%
    rename(Longitude = "decimalLongitude") %>%
    rename(Latitude = "decimalLatitude")
  
  #Look for bias file
  bias_files <- list.files()
  if(is_empty(bias_files)) stop("No bias file")
  bias <- raster() #does it have a crs after reading in?
  
  #Read in environmental data
  Env_range <- raster()
  #Make sure env and bias are the same extent, resolution, etc. Should be but check anyway
  
  #Mask occurrences that aren't within Env_range?
  
  #Create enmtools.species object
  species <- enmtools.species(species.name = i, 
                              presence.points = occs_cl) #must be named Longitude and Latitude
  
  #Run maxent model
  maxent_mod <- enmtools.maxent(species = species,
                               env = Env_Range,
                               test.prop = 0.2, #also try with "block"
                               bias = Om_bias,
                               args = c("outputformat=raw", 
                                        "replicates=10", 
                                        "replicatetype=crossvalidate"))
  
  #Evaluate
  env.evaluate(species, maxent_mod, Env_Range, test.eval = T)
  
  #Write predicted distribution layer for input into zonation
  zonation_data_path <- 
  writeRaster()
