###################################################### Calculating range offset - cluster edition ###############################

#                                                 MAKE SURE ALL DIRECTORY PATHS ARE CORRECT FOR CLUSTER

################################################################ Libraries #################################################
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(bossMaps)
library(data.table)
library(CoordinateCleaner)
library(countrycode)

#################################################################XY function###########################################
xy_match <- function(initial, pre, Layer) {
  pre <- initial
  
  coordinates(pre) <- ~ decimalLongitude + decimalLatitude
  crs(pre) <- crs(Layer)
  
  final <- cbind(initial, over(pre, Layer))
  return(final)
}

###########################################################Spatial information
#Load workspace (Coast shapefile and species list)
load("/projects/nc57/Chapter_3/SpatialData/RL_shp_prepro.RData")
RL_shp_Aus <- RL_shp_Aus[-c(1155, 1228, 1528, 2358, 5014, 5181, 5327),]
env_predictors <- stack("/projects/nc57/Chapter_3/SpatialData/Raster/env_predictors.gri")

#Paths to local drive
env_path <- "/home/dcla0008/nc57_scratch/Maxent/Env_layers"
bias_path <- "/home/dcla0008/nc57_scratch/Maxent/Bias_layers"
occ_path <- "/home/dcla0008/nc57_scratch/Maxent/Occurrences"
gbif_path <- "/projects/nc57/Chapter_3/SpatialData/Vector/SpeciesOccurrences/0310609-200613084148143.csv"
ala_path <- "/projects/nc57/Chapter_3/SpatialData/Vector/SpeciesOccurrences/ALA_occs.csv"

res.out <- lapply(speciesNames_range[1:length(speciesNames_range)], function(i) {
  
  print(i)
  
  i <- gsub(" ", ".", i)
  
  my_gbif_cols <- fread(cmd = paste("head -n 1", gbif_path))
  my_ala_cols <- fread(cmd = paste("head -n 1", ala_path))
  my_ala_cols <- as.character(my_ala_cols[1,])
  
  my_cols_min <- c("family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince", "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode")
  
  
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
    occs$countryCode <-  countrycode(occs$countryCode, origin =  'iso2c', destination = 'iso3c')
    
    #flag problems
    flags <- clean_coordinates(x = occs,
                               lon = "decimalLongitude",
                               lat = "decimalLatitude",
                               countries = "countryCode",
                               species = "species",
                               tests = c("centroids", "equal","gbif", "institutions",
                                         "zeros", "countries", "duplicates"))
    
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
    
    #write csv of flagged records
    #write.csv(gbif_occs_fl, file = file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "Flagged", paste0(i,"_fl.csv")))
    
    #}, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
    
    #tryCatch({
    
    i <- gsub("_", " ", i)
    
    #selecting species
    Range <- RL_shp_Aus %>% dplyr::filter(acceptedName == i)
    
    #Checking geometry type
    # if(st_geometry_type(Range) == "GEOMETRYCOLLECTION") {
    #   
    #   Range <- st_cast(Range)
    # }
    Range <- st_buffer(Range, dist = 1000)
    
    #convert to sp
    Range <- as(Range, Class = "Spatial")
    
    #Domain
    #environmental data or occurrence extent
    # occ_extent <- extent(min(occs_cl$decimalLongitude),
    #                      max(occs_cl$decimalLongitude),
    #                      min(occs_cl$decimalLatitude),
    #                      max(occs_cl$decimalLatitude))
    # 
    # if(extent(Range) > occ_extent) {
    #     
    Env_Range <- crop(env_predictors[[1]], Range)
    #Env_Range <- raster::stack(Env_Range)
    #} else {
    
    #   Range@bbox <- as.matrix(extent(occ_extent))
    #   Env_Range <- crop(env_predictors, occ_extent)
    #   Env_Range <- raster::stack(Env_Range)
    # }
    
    #Calculate distance-to-range
    #range and domain need to be the same extent (thus same projection)
    start_time <- Sys.time()
    
    rdist <- rangeDist(
      range = Range, 
      domain = Env_Range, #I wonder if this should be based on occurrence records?
      domainkm = 1000,
      mask = F,
      fact = 6
    )
    end_time <- Sys.time()
    end_time - start_time
    
    #Crop "environmental data" to this new domain
    Env_domain <- crop(env_predictors, rdist)
    #write environmental layer
    i <- gsub(" ", "_", i)
    writeRaster(Env_domain, filename = paste0(env_path,"/",i,"/",i, "_env.grd"), overwrite = T)
    
    # Mask pixels with no environmenta data (over ocean, etc.)
    rdist <- mask(rdist, Env_domain[[1]])
    names(rdist) <- "rangeDist"
    
    #Evaluate curves
    #rates <- checkRates(rdist, skew = 0.2)
    
    #Calculate frequency table of distances
    dists <- freq(rdist, useNA = "no", digits = 2)
    
    #Create offset
    #bias path
    expert <- rangeOffset(
      rdist,
      parms = c(rate = 1.0, skew = 0.2, shift = 0, prob = 0.9),
      dists = dists,
      doNormalize = T,
      verbose = T,
      doWriteRaster = T,
      filename = paste0(bias_path,"/",i,"/",i, "_bias.grd"), 
      overwrite = T,
      datatype = "FLT4S"
    )
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  #, write.error.dump.file = T, 
  #write.error.dump.folder = "C:/Users/dcla0008/Documents/Maxent")
  
  gc()
})
res.out
