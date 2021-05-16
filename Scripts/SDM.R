# Distribution modelling

#try to use letsR (e.g. lets.presab.points())
#see here: https://rmacroecology.netlify.app/2018/01/23/a-guide-to-transform-species-shapefiles-into-a-presence-absence-matrix-based-on-a-user-defined-grid-system/

# Species to be modelled
#spec_names <- c("Vespa velutina", "Bombus impatiens", "Lymantria dispar", "Harmonia axyridis", 
#                "Solenopsis invicta", "Wasmannia auropunctata", "Glischrochilus quadrisignatus",
#                "Digitonthophagus gazella", "Icerya purchasi", "Pheidole megacephala")

# For testing: Pheidole megacephala
Pmeg <- "Pheidole megacephala"

#Obtaining GBIF occurrences
#spec_names_occs <- Get_gbif_Occ(spec_names)
Pmeg_occs <- Get_gbif_Occ(Pmeg)

#Convert to a spatial object
Pmeg_occs_sp <- SpatialPointsDataFrame(coords=cbind(Pmeg_occs$decimalLongitude, Pmeg_occs$decimalLatitude),proj4string = CRS(projection(Coast_shp)), data=Pmeg_occs)
Pmeg_occs_sf <- st_as_sf(Pmeg_occs.sp)

#Occurrence data cleaning

  #look into CoordinateCleaner package 
  #(https://ropensci.github.io/CoordinateCleaner/articles/Cleaning_GBIF_data_with_CoordinateCleaner.html)

#Modelling
  
  #look into biomod2 package 