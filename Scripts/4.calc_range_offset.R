############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "4.calc_range_offset.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "4.calc_range_offset.R"),
              local_path = file.path("Scripts", "4.calc_range_offset.R"), overwrite = T)
##################################################################################################
                              # Creating bias file (with bossMaps) #

#species names required are only those from the shapefile
speciesNames_range <- read.csv(file.path("SpeciesData", "RL_shp_Aus.csv"))
speciesNames_range <- unique(speciesNames_range$x)



#Paths to local drive (VM D drive)
env_path <- "D:/Chapter_3/Chapter_3/Maxent/Env_layers"
bias_path <- "D:/Chapter_3/Chapter_3/Maxent/Bias_layers"
occ_path <- "D:/Chapter_3/Chapter_3/Maxent/Occurrences"
gbif_path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "0310609-200613084148143.csv"))
ala_path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "SpeciesOccurrences", "ALA", "ALA_occs.csv"))

lapply(speciesNames_range[1:length(speciesNames_range)], function(i) {
  
  print(i)
  
  i <- gsub(" ", ".", i)
  
  my_gbif_cols <- fread(cmd = paste("head -n 1", gbif_path))
  my_ala_cols <- fread(cmd = paste("head -n 1", ala_path))
  my_ala_cols <- as.character(my_ala_cols[1,])
  
  my_cols_min <- c("family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince", "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode")
  
  
  tryCatch({
    gbif_occs <- fread(cmd = paste("grep", i ,gbif_path), col.names = names(my_gbif_cols), quote = "", na.strings = c("",NA))
    gbif_occs <- gbif_occs %>% dplyr::select(my_cols_min)
    ala_occs <- fread(cmd = paste("grep", i ,ala_path), col.names = my_ala_cols, na.strings = c("",NA))
    ala_occs <- ala_occs %>% dplyr::select(my_cols_min)
    occs <- rbind(gbif_occs, ala_occs)
    if(is_empty(occs)) stop("No occurrence records")
    
    
    #think its here where I include cleaning stuff
    
    # remove records without coordinates
    occs <- occs %>%
      filter(!is.na(decimalLongitude)) %>%
      filter(!is.na(decimalLatitude))
    
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
    occs_fl <- occs[!flags$.summary,]
    
    #Excluding records with more than 500m uncertainty
    #occs_cl <- occs_cl %>%
    #  filter(coordinateUncertaintyInMeters <= 500 | is.na(coordinateUncertaintyInMeters))
    
    #write csv of cleaned records
    i <- gsub("\\.", "_", i)
    write.csv(occs_cl, file = file.path(occ_path,i, paste0(i,"_cl.csv")))
    
    #write csv of flagged records
    #write.csv(gbif_occs_fl, file = file.path("SpatialData", "Vector", "SpeciesOccurrences", "GBIF", "Flagged", paste0(i,"_fl.csv")))
    
  #}, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  
  #tryCatch({
  
  i <- gsub("_", " ", i)
    
  #selecting species
  Range <- RL_shp_Aus %>% dplyr::filter(acceptedName == i)

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
    Env_Range <- crop(env_predictors, Range)
    Env_Range <- raster::stack(Env_Range)
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
              fact = 4
              )
  end_time <- Sys.time()
  end_time - start_time

  #Crop "environmental data" to this new domain
  Env_domain <- crop(env_predictors, rdist)
  #write environmental layer
  i <- gsub(" ", "_", i)
  writeRaster(Env_domain, filename = file.path(env_path, i, paste0(i, "_env.grd")), overwrite = T)
  
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
    filename = file.path(bias_path,i, paste0(i, "_bias.grd")), 
    overwrite = T,
    datatype = "FLT4S"
    )
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  #, write.error.dump.file = T, 
  #write.error.dump.folder = "C:/Users/dcla0008/Documents/Maxent")
  
  gc()

})

#Can extend offset raster using extend(). Prob needed for zonation.
#Actually may be required for predicted SDM (bias file if not making sdm)
#example
expert <- extend(expert, Aus_Coast)
