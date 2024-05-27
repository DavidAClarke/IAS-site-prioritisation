############################## plant point bias files ##########################

#Load data
env_predictors <- stack("/projects/nc57/Chapter_3/SpatialData/Raster/env_predictors.gri")
load("/projects/nc57/Chapter_3/SpatialData/RL_shp_prepro.RData")
rm(RL_shp_Aus)

#Only species names required are those with Red List plant points
RL_plants <- read.csv("/projects/nc57/Chapter_3/SpatialData/Vector/RL_plants_points.csv")

#Need remove species with few/one points e.g. "Scleria tessellata" has only one point
RL_plants <- subset(RL_plants, 
                    with(RL_plants, unsplit(table(acceptedName), 
                                            acceptedName)) >= 5)
speciesNames_points <- unique(RL_plants$acceptedName)

#Paths to local drive
env_path <- "/home/dcla0008/nc57_scratch/Maxent/Env_layers"
bias_path <- "/home/dcla0008/nc57_scratch/Maxent/Bias_layers"
occ_path <- "/home/dcla0008/nc57_scratch/Maxent/Occurrences"
gbif_path <- "/projects/nc57/Chapter_3/SpatialData/Vector/SpeciesOccurrences/0310609-200613084148143.csv"
ala_path <- "/projects/nc57/Chapter_3/SpatialData/Vector/SpeciesOccurrences/ALA_occs.csv"


lapply(speciesNames_points[1:length(speciesNames_points)], function(i){

  print(i)
  
  i <- gsub(" ", ".", i)
  
  my_gbif_cols <- fread(cmd = paste("head -n 1", gbif_path))
  my_ala_cols <- fread(cmd = paste("head -n 1", ala_path))
  my_ala_cols <- as.character(my_ala_cols[1,])
  
  my_cols_min <- c("family","species", "taxonRank","scientificName",
                   "decimalLongitude","decimalLatitude","countryCode",
                   "stateProvince", "coordinateUncertaintyInMeters", 
                   "coordinatePrecision","basisOfRecord","year", 
                   "individualCount", "institutionCode")
  
  
  tryCatch({
    gbif_occs <- fread(cmd = paste("grep", i ,gbif_path), 
                       col.names = names(my_gbif_cols), 
                       colClasses = "character",
                       quote = "", 
                       na.strings = c("",NA))
    ala_occs <- fread(cmd = paste("grep", i ,ala_path), 
                      col.names = my_ala_cols, 
                      colClasses = "character",
                      na.strings = c("",NA))
    if(is_empty(gbif_occs) & is_empty(ala_occs)) stop("No occurrence records")
    if(is_empty(gbif_occs) & is_empty(ala_occs) == F){
      print("No GBIF but there are ALA records")
      ala_occs <- ala_occs %>% 
        dplyr::select(my_cols_min)
      occs <- ala_occs %>%
        dplyr::mutate(decimalLongitude = as.numeric(decimalLongitude)) %>%
        dplyr::mutate(decimalLatitude = as.numeric(decimalLatitude)) %>%
        dplyr::mutate(coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)) %>%
        dplyr::mutate(coordinatePrecision = as.numeric(coordinatePrecision)) %>%
        dplyr::mutate(individualCount = as.integer(individualCount))
    } else 
      if(is_empty(gbif_occs) == F & is_empty(ala_occs)){
        print("No ALA but there are GBIF records")
        gbif_occs <- gbif_occs %>% 
          dplyr::select(my_cols_min) 
        occs <- gbif_occs %>%
          dplyr::mutate(decimalLongitude = as.numeric(decimalLongitude)) %>%
          dplyr::mutate(decimalLatitude = as.numeric(decimalLatitude)) %>%
          dplyr::mutate(coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)) %>%
          dplyr::mutate(coordinatePrecision = as.numeric(coordinatePrecision)) %>%
          dplyr::mutate(individualCount = as.integer(individualCount))
      } else 
        if(is_empty(gbif_occs) == F & is_empty(ala_occs) == F){
          print("Both GBIF and ALA records present")
          gbif_occs <- gbif_occs %>% 
            dplyr::select(my_cols_min)
          ala_occs <- ala_occs %>% 
            dplyr::select(my_cols_min)
          occs <- rbind(gbif_occs, ala_occs)
          occs <- occs %>%
            dplyr::mutate(decimalLongitude = as.numeric(decimalLongitude)) %>%
            dplyr::mutate(decimalLatitude = as.numeric(decimalLatitude)) %>%
            dplyr::mutate(coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)) %>%
            dplyr::mutate(coordinatePrecision = as.numeric(coordinatePrecision)) %>%
            dplyr::mutate(individualCount = as.integer(individualCount))
        }
    
    
    # remove records without coordinates
    occs <- occs %>%
      dplyr::filter(!is.na(decimalLongitude)) %>%
      dplyr::filter(!is.na(decimalLatitude))
    
    #convert country code from ISO2c to ISO3c
    occs$countryCode <-  countrycode(occs$countryCode, 
                                     origin =  'iso2c', destination = 'iso3c')
    
    #flag problems
    flags <- clean_coordinates(x = occs,
                               lon = "decimalLongitude",
                               lat = "decimalLatitude",
                               countries = "countryCode",
                               species = "species",
                               tests = c("centroids", "equal","gbif", 
                                         "institutions","zeros", "countries", 
                                         "duplicates"))
    
    #Exclude problematic records
    occs_cl <- occs[flags$.summary,]
    
    #Remove points that don't fall within Country
    Coast_sp <- as(Aus_Coast, Class = "Spatial")
    occs_cl <- xy_match(occs_cl, occs_cl_n, Coast_sp)
    occs_cl <- occs_cl %>%
      drop_na(id)
    
    #The flagged records
    #occs_fl <- occs[!flags$.summary,]
    
    #Excluding records with more than 500m uncertainty
    #occs_cl <- occs_cl %>%
    #  filter(coordinateUncertaintyInMeters <= 500 | is.na(coordinateUncertaintyInMeters))
    
    #write csv of cleaned records
    i <- gsub("\\.", "_", i)
    #write.csv(occs_cl, file = file.path(occ_path,i, paste0(i,"_cl.csv")))
    write.csv(occs_cl, file = paste0(occ_path,"/",i,"/",i,"_cl.csv"))
    
    #Choose species
    i <- gsub("_", " ", i)
    species <- RL_plants %>% dplyr::filter(acceptedName == i)
    species <- SpatialPointsDataFrame(species[,7:8], species, 
                                      proj4string = crs(env_predictors))
    
    #predictors <- mask(predictors, Aus_Coast) #taken care of in env_predictors
    predictors <- crop(env_predictors, species) #use species occurrence extent
    i <- gsub(" ", "_", i)
    writeRaster(predictors, 
                filename = paste0(env_path,"/", i,"/",i, "_env.grd"), 
                overwrite = T)
    #model.extent <- extent(min(occs$decimalLongitude)-10,
      #max(occs$decimalLongitude)+10,
      #min(occs$decimalLatitude)-10,
      #max(occs$decimalLatitude)+10)

    #rasterise
    species_ras <- rasterize(species, predictors, 1)
    species_ras <- mask(species_ras, Aus_Coast) #%>% crop(Aus_Coast)

    #Get 2D kernal density estimate
    presences <- which(values(species_ras) == 1)
    
    pres.locs <- coordinates(species_ras)[presences, ]
    if(is_empty(pres.locs)) stop("pres.locs is empty")
    if(quantile(pres.locs)[2]-quantile(pres.locs)[3] == 0) {
      pres.locs + sample(c(-0.0001, 0.0001), 
                         size = length(pres.locs), 
                         replace = T)}
    dens <- kde2d(pres.locs[,1], pres.locs[,2], 
                  n = c(nrow(species_ras), ncol(species_ras)))
    dens.ras <- raster(dens)
    bias <- resample(dens.ras, predictors)
    path <- paste0(bias_path,"/", i)
    writeRaster(bias, filename = paste0(path,"/",i, "_bias.grd"), 
                overwrite = T) #change
    }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  })

