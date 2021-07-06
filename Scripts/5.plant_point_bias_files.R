                               # Plant bias files
#Only species names required are those with Red List plant points
RL_plants <- read.csv(file.path("SpatialData", "Vector", "redlist_species_data_plantpoints","RL_plants_points.csv"))
#Need remove species with few/one points e.g. "Scleria tessellata" has only one point
RL_plants <- subset(RL_plants, with(RL_plants, unsplit(table(acceptedName), acceptedName)) >= 5)
speciesNames_points <- unique(RL_plants$acceptedName)

lapply(speciesNames_points[1:length(speciesNames_points)], function(i){

    print(i)
  
    #Choose species
    species <- RL_plants %>% dplyr::filter(acceptedName == i)
    species <- SpatialPointsDataFrame(species[,7:8], species, proj4string = crs(predictors))
    #predictors <- mask(predictors, Aus_Coast) #already taken care of in env_predictors
    predictors <- crop(predictors, species) #maybe use species occurrence extent
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
    dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(species_ras), ncol(species_ras)))
    dens.ras <- raster(dens)
    bias <- resample(dens.ras, predictors)
    i <- gsub(" ", "_", i)
    path <- file.path(bias_path, i)
    writeRaster(bias, filename = file.path(path, paste0(i, "_bias.asc")), overwrite = T) #change
  })

