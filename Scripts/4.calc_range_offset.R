############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "calc_range_offset.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "calc_range_offset.R"),
              local_path = file.path("Scripts", "calc_range_offset.R"), overwrite = T)
##################################################################################################
                              # Creating bias file (with bossMaps) #

#species names required are only those from the shapefile
RL_Aus_Range <- read.csv(file.path("SpeciesData", "RL_Aus_Range.csv"))
speciesNames_range <- unique(RL_Aus_Range$acceptedName)
                                             
lapply(speciesNames_range[1:length(speciesNames_range)], function(i) {
  
  print(i)
  
  #selecting species
  Range <- RL_shp_Aus %>% dplyr::filter(acceptedName == i)

  #convert to sp
  Range <- as(Range, Class = "Spatial")

  #Domain
  #environmental data
  Env <- Aus_elev #will be combination of climate and elevation
  Env_Range <- crop(Env, Range)

  #Calculate distance-to-range
  #range and domain need to be the same extent (thus same projection)
  start_time <- Sys.time()
  rdist <- rangeDist(
              range = Range, 
              domain = Env_Range, #I wonder if this should be based on occurrence records?
              #model.extent <- extent(min(occs$decimalLongitude)-10,
                                    #max(occs$decimalLongitude)+10,
                                    #min(occs$decimalLatitude)-10,
                                    #max(occs$decimalLatitude)+10)
              domainkm = 1000,
              mask = F,
              fact = 6
              )
  end_time <- Sys.time()
  end_time - start_time

  #Crop "environmental data" to this new domain
  Env_domain <- crop(Env, rdist)

  # Mask pixels with no environmenta data (over ocean, etc.)
  rdist <- mask(rdist, Env_domain[[1]])
  names(rdist) <- "rangeDist"

  #Evaluate curves
  #rates <- checkRates(rdist, skew = 0.2)

  #Calculate frequency table of distances
  dists <- freq(rdist, useNA = "no", digits = 2)

  #Create offset
  #bias path
  i <- gsub(" ", "_", i)
  path <- file.path(bias_path, i)
  expert <- rangeOffset(
    rdist,
    parms = c(rate = 0.215, skew = 0.2, shift = 0, prob = 0.9),
    dists = dists,
    doNormalize = T,
    verbose = T,
    doWriteRaster = T,
    filename = file.path(path, paste0(i, "_bias.asc")), 
    overwrite = T,
    datatype = "FLT4S"
    )

})

#Can extend offset raster using extend(). Prob needed for zonation.
#Actually may be required for predicted SDM (bias file if not making sdm)
#example
expert <- extend(expert, Aus_Coast)
