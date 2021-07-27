                               # Plant bias files
#Only species names required are those with Red List plant points
RL_plants <- read.csv(file.path("SpatialData", "Vector", "redlist_species_data_plantpoints","RL_plants_points.csv"))
#Need remove species with few/one points e.g. "Scleria tessellata" has only one point
RL_plants <- subset(RL_plants, with(RL_plants, unsplit(table(acceptedName), acceptedName)) >= 5)
speciesNames_points <- unique(RL_plants$acceptedName)

### I am up to using "Croton choristadenius": it has "only finite values are allowed in 'lims'" issue with kde2d

start_time <- Sys.time()
lapply(speciesNames_points[1:length(speciesNames_points)], function(i){

    print(i)
    tryCatch({
    #Choose species
    species <- RL_plants %>% dplyr::filter(acceptedName == i)
    species <- SpatialPointsDataFrame(species[,7:8], species, proj4string = crs(env_predictors))
    
    #predictors <- mask(predictors, Aus_Coast) #already taken care of in env_predictors
    predictors <- crop(env_predictors, species) #maybe use species occurrence extent
    i <- gsub(" ", "_", i)
    writeRaster(predictors, filename = file.path(env_path, i, paste0(i, "_env.grd")), overwrite = T)
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
    if(quantile(pres.locs)[2]-quantile(pres.locs)[3] == 0) {pres.locs + sample(c(-0.0001, 0.0001), size = length(pres.locs), replace = T)}
    dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(species_ras), ncol(species_ras)))
    dens.ras <- raster(dens)
    bias <- resample(dens.ras, predictors)
    path <- file.path(bias_path, i)
    writeRaster(bias, filename = file.path(path, paste0(i, "_bias.grd")), overwrite = T) #change
    }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  })
end_time <- Sys.time()
end_time - start_time
