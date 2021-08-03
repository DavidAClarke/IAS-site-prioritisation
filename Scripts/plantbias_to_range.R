                                    #Converting plant bias files to range maps
#Get species list (from RL plant points)
RL_plants <- read.csv("E:/SpatialData/Vector/redlist_species_data_plantpoints/RL_plants_points.csv")
RL_plants <- unique(RL_plants$acceptedName)

#Set bias path
Bias_path <- "E:/Maxent/Bias_layers"

#cutoff function
rc <- function(x) {
       
    ifelse(x <  0.005, 0,
    ifelse(x >=  0.005, 1, NA)) }

lapply(RL_plants[1:length(RL_plants)], function(i) {
  
  print(i)
  
  i <- gsub(" ", "_", i)
  
  #Load in bias file
  tryCatch({
  bias_present <- list.files(file.path(Bias_path, i), pattern = ".gri$")
  if(is_empty(bias_present) == T) stop("No bias files")
  bias <- raster(file.path(Bias_path, i, paste0(i,"_bias.gri")))
  
  #Convert to binary
  bias_bin <- calc(bias, fun=rc)

  #Extend the extent to match range maps
  bias_bin <- extend(bias_bin, RL_shp_Aus)
  bias_bin <- mask(bias_bin, Aus_Coast)

  #Write raster to Input_data for zonation
  writeRaster(bias_bin, filename = file.path(Bias_path, i, paste0(i,"_bias_bin.tif")), overwrite = T)
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
})
