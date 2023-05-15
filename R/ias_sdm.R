# IAS SDM - Setup largely inlfuenced by Polaina et al 2020

#Required libraries
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(data.table)
library(rgbif)
library(biomod2)
library(countrycode)
library(CoordinateCleaner)
library(usdm)
#library(rgdal)

#Required initial data
#bio <- getData("worldclim", var = "bio", res = 2.5, path = "E:/SpatialData/Raster/Worldclim")
#Personal PC = F:, original work PC = E:
bio_names <- list.files(path = "C:/Users/david/Documents/PhD/Chapter_3/SpatialData/Raster/Worldclim/wc2-5", pattern = ".bil$", full.names = T)
bio <- stack(bio_names)
source("Scripts/2.aus_coast.R") #if not already loaded

#Need an empty reference raster (global)
ext <- extent(bio)
crs <- "+proj=longlat +datum=WGS84 +no_defs"
res <- 0.04166667
vals <- 0
ref_raster <- raster(ext = ext, crs = crs, resolution = res, vals = vals)

#Need an Australian raster that is a subset of reference raster (Regional)
Aus_raster <- crop(ref_raster, extent(Aus_Coast))

## Create directories and Set all required paths
#Gloal models
#dir.create(file.path("SpatialData", "IAS_distributions", "IAS_global"))
global_model_path <- paste0("C:/Users/david/Documents/PhD/Chapter_3/SpatialData/IAS_distributions/IAS_global")

#Regional models
#dir.create(file.path("SpatialData", "IAS_distributions", "IAS_regional"))
regional_model_path <- paste0("C:/Users/david/Documents/PhD/Chapter_3/SpatialData/IAS_distributions/IAS_regional")

#Maxent path
maxent_jar_path <- "C:/Users/david/Documents/PhD/Chapter_3/maxent.jar" #could change. No spaces allowed

# #Insect species to model (these all worked)
# spp_list <- c("Digitonthophagus gazella", "Pheidole megacephala",
#               "Vespula germanica", "Tetramorium bicarinatum", "Paratrechina longicornis")

#Additional species to model (some may not work)
spp_list <- c("Apis mellifera", "Solenopsis geminata", "Monomorium floricola", "Monomorium destructor", 
              "Linepithema humile", "Vespula vulgaris", "Polistes chinensis antennalis", 
              "Megachile rotundata", "Bombus terrestris", "Wasmannia auropunctata", "Apis cerana", 
              "Solenopsis invicta", "Heteronychus arator", "Digitonthophagus gazella", "Pheidole megacephala",
              "Vespula germanica", "Tetramorium bicarinatum", "Paratrechina longicornis")
#need to download Polistes separately

#Download occurrence data
#GBIF credentials (remove credentials when publishing scripts)
user <- "davidclarke"
pwd <- "GBIFN#8208!"
email <- "david.clarke@latrobe.edu.au"

#Match the names
gbif_taxon_keys <- spp_list %>%
  taxize::get_gbifid_(method = "backbone") %>%
  imap(~ .x %>% mutate(original_sciname = .y)) %>%
  bind_rows() %T>%
  write.csv(file = file.path("SpeciesData", "IAS_all_matches.csv")) %>% #can then examine the data
  dplyr::filter(matchtype == "EXACT" & status == "ACCEPTED") %>%
  pull(usagekey)

#Request download
res <- occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("basisOfRecord", c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION', 'PRESERVED_SPECIMEN')),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email)

#Download data
occ_path <- "C:/Users/david/Documents/PhD/Chapter_3/SpatialData/Vector/IAS_Occurrences"
#occ_download_get("0325811-200613084148143", path = occ_path) #initial species
#occ_download_get("0376270-210914110416597", path = occ_path) #additional species
occ_download_get("0381185-210914110416597", path = occ_path)

#Occurrence data path
#data_path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "IAS_Occurrences", "0325811-200613084148143.csv"))
#data_path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "IAS_Occurrences", "0376270-210914110416597.csv"))#additional species
data_path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "IAS_Occurrences", "0381185-210914110416597.csv"))

#Add decimal point for use in fread()
spp_list <- gsub(" ", ".", spp_list)

#Obtain just column names of occurrence dataset
my_cols <- fread(cmd = paste("head -n 1", data_path))

#Columns to get subset
my_cols_min <- c("family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince", "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode")

#Convert occurrences to shapefiles
lapply(spp_list[1:length(spp_list)], function(i) {
  
  #Load data
  gbif_occs <- fread(cmd = paste("grep", i ,data_path), 
                     col.names = names(my_cols), 
                     quote = "", 
                     na.strings = c("",NA))
  
  gbif_occs <- gbif_occs %>% 
    dplyr::select(my_cols_min)
  
  #remove records without coordinates (just in case)
  gbif_occs <- gbif_occs %>%
    dplyr::filter(!is.na(decimalLongitude)) %>%
    dplyr::filter(!is.na(decimalLatitude))
  
  #convert country code from ISO2c to ISO3c
  gbif_occs$countryCode <-  countrycode(gbif_occs$countryCode, 
                                        origin =  'iso2c', 
                                        destination = 'iso3c')
  
  #Convert to spatial object
  gbif_occs_shp <- SpatialPointsDataFrame(coords=cbind(gbif_occs$decimalLongitude, gbif_occs$decimalLatitude),
                                         proj4string = CRS(projection(ref_raster)),
                                         data=as.data.frame(gbif_occs))
  gbif_occs_shp <- st_as_sf(gbif_occs_shp)
  
  #Remove . in species name for writing to shapefile
  i <- gsub("\\.", "_", i)
  
  #Save as shapefile
  st_write(gbif_occs_shp, 
           dsn = file.path(occ_path, paste0(i,"_data.shp")), 
           append = F)
           
})

#Clean occurrence data
radius_uncert <- 2500 #radius of cell size: 2.5 arc min very ~ 5km
  
lapply(spp_list[1:length(spp_list)], function(i) {
  
  #Repalce decimal in species name
  i <- gsub("\\.", "_", i)
  
  #Load shapefile
  #occ_shp <- readOGR(file.path(occ_path,paste0(i,"_data.shp")))
  occ_shp <- st_read(file.path(occ_path,paste0(i,"_data.shp")))
  
  #Removing those with uncertainty > radius_uncert and NAs
  occ_shp_cert <- subset(occ_shp, crdnUIM<= radius_uncert)
  occ_shp_cert <- as_Spatial(occ_shp_cert)
  
  
  #flag problems
  flags <- clean_coordinates(x = occ_shp_cert,
                             lon = "dcmlLng",
                             lat = "dcmlLtt",
                             countries = "cntryCd",
                             species = "species", 
                             tests = c("centroids", "equal","gbif", "institutions",
                                       "zeros"))
  
  #Exclude problematic records
  occ_shp_cert_cl <- occ_shp_cert[flags$.summary,]
  
  #The flagged records
  occ_shp_cert_fl <- occ_shp_cert[!flags$.summary,]
  
  ## Getting unique values per grid cell
  coord_data <- occ_shp_cert_cl[,c("dcmlLng","dcmlLtt")]
  extract1 <- raster::extract(ref_raster, coord_data, cellnumbers = T)
  coord_raster_ref_sp <- as.data.frame(coordinates(ref_raster)[extract1[,1],])
  df1 <- as.data.frame(coord_raster_ref_sp) # OCURRENCES
  df2 <- df1[!duplicated(df1[,c('x','y')]),]
  df2 <- na.omit(df2) # PRESENCES
  
  # Converting to a spatial object - OCCURRENCES
  occ.sp <- SpatialPointsDataFrame(coords=cbind(df1$x, df1$y),
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                                   data=as.data.frame(cbind(df1$x, df1$y)))
  names(occ.sp) <- c("x","y")
  occ_sf <- st_as_sf(occ.sp)
  
  # Converting to a spatial object - PRESENCES
  pres.sp <- SpatialPointsDataFrame(coords=cbind(df2$x, df2$y),
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                                    data=as.data.frame(cbind(df2$x, df2$y)))
  names(pres.sp) <- c("x","y")
  pres_sf <- st_as_sf(pres.sp)
  
  #Save as shapefile
  st_write(occ_sf, 
           dsn = file.path(occ_path,paste0(i,"_occ.shp")),
           append = F)
  
  st_write(pres_sf, 
           dsn = file.path(occ_path,paste0(i,"_pres.shp")),
           append = F)
  
  
  #Keeping the NAs as well as those with "certain" data
  occ_shp_na <- occ_shp[is.na(occ_shp$crdnUIM),]
  occ_shp_na <- as_Spatial(occ_shp_na)
  
  #flag problems
  flags <- clean_coordinates(x = occ_shp_na,
                             lon = "dcmlLng",
                             lat = "dcmlLtt",
                             countries = "cntryCd",
                             species = "species",
                             tests = c("centroids", "equal","gbif", "institutions",
                                       "zeros"))
  
  #Exclude problematic records according to the previous function:
  occ_shp_na_cl <- occ_shp_na[flags$.summary,]
  
  #The flagged records
  occ_shp_na_fl <- occ_shp_na[!flags$.summary,]
  
  #Joining to the previous set (cert) to save it
  occ_shp_cl <- rbind(occ_shp_cert_cl, occ_shp_na_cl)
  
  #Getting unique values per grid cell
  coord_data <- occ_shp_cl[,c("dcmlLng","dcmlLtt")]
  extract1 <- raster::extract(ref_raster,coord_data,cellnumbers = T)
  coord_raster_ref_sp <- as.data.frame(coordinates(ref_raster)[extract1[,1],])
  df1 <- as.data.frame(coord_raster_ref_sp) # OCCURRENCES
  df2 <- df1[!duplicated(df1[,c('x','y')]),]
  df2 <- na.omit(df2) ## PRESENCES
  
  #Converting to a spatial object - OCURRENCES
  occ.sp2 <- SpatialPointsDataFrame(coords=cbind(df1$x, df1$y),
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                                    data=as.data.frame(cbind(df1$x, df1$y)))
  names(occ.sp2) <- c("x","y")
  occ_sf_2 <- st_as_sf(occ.sp2)
  
  #Converting to a spatial object - PRESENCES
  pres.sp2 <- SpatialPointsDataFrame(coords=cbind(df2$x, df2$y),
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                                     data=as.data.frame(cbind(df2$x, df2$y)))
  names(pres.sp2) <- c("x","y")
  pres_sf_2 <- st_as_sf(pres.sp2)
  
  #Save as shapefile
  st_write(occ_sf_2, 
           dsn = file.path(occ_path,paste0(i,"_occ_all.shp")),
           append = F)
  
  st_write(pres_sf_2, 
           dsn = file.path(occ_path,paste0(i,"_pres_all.shp")),
           append = F)
  
})
###########################################################################################
#PREDICTORS 

#Climate data (2.5 minutes of a degree)
#Reducing collinearity
nocorrvar <- vifstep(bio, th = 4)
predictors <- as.character(nocorrvar@results$Variables)
clim_sub <- raster::subset(bio,c(predictors))

#Set parameters for the model
name_model_folder <- "IAS_global"
model_class <- "global"

#Individual models
my_models <- c("GLM","GAM","MAXENT.Phillips","FDA","GBM") # Algorithms to fit the models   
my_pa_strategy <- c("random") # Strategy to select PAs   
my_pa_dist_min <- 0 # Distance where to start throwing PAs  
my_runs <- 4 # Number of PAs runs
my_n_pa <- 20000 # Number of pseudoabsences, when not dependent of no. presences
my_proj_name <- "global"  # Name for the output projections

lapply(spp_list[1:length(spp_list)], function(i){
  
  #Remove spaces in species name
  i <- gsub("\\.", "_", i)
  
  #Load in presence shapefile
  #sp_pres <- readOGR(file.path(occ_path,paste0(i,"_pres_all.shp")))
  sp_pres <- st_read(file.path(occ_path,paste0(i,"_pres_all.shp")))
  n_mod <- dim(sp_pres)[1]
  
  #Pseudoabsences
  n_pa <- my_n_pa
  
                                ## Format data for 'BIOMOD'
  #Presences
  mypres <- rep(1,dim(sp_pres)[1])
  
  #Presences coordinates
  sp_pres <- as_Spatial(sp_pres)
  mycoord <- coordinates(sp_pres)
  
  #Clean data
  clean_data <- BIOMOD_FormatingData(resp.var=mypres, 
                                     expl.var=clim_sub, 
                                     resp.xy=mycoord,
                                     resp.name=i, 
                                     PA.nb.rep=my_runs,
                                     PA.nb.absences=n_pa,
                                     PA.strategy = my_pa_strategy, 
                                     PA.dist.min=my_pa_dist_min)
  
  #Settings for the model algorithms
  model_opt <- BIOMOD_ModelingOptions(GLM=list(type='quadratic', interaction.level=0),
                                      GBM=list(n.trees=1000),
                                      MAXENT.Phillips = list(path_to_maxent.jar = maxent_jar_path, product=FALSE),
                                      GAM=list(k=3))
  
  #Path for models
  #must be full path (back to e.g. C:, and have no spaces)
  #global_model_path <- paste0("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/", name_model_folder)
  setwd(global_model_path)
  
  #Run the models
  #Separate
  all_models <- BIOMOD_Modeling(data=clean_data, 
                                models=my_models, 
                                models.options=model_opt, 
                                NbRunEval = 4, 
                                DataSplit = 70, 
                                VarImport = 3, 
                                do.full.models = F, 
                                modeling.id = "global")
  
  #Project the models
  #Separate
  models_proj_current <- BIOMOD_Projection(modeling.output = all_models,
                                           new.env = clim_sub,
                                           proj.name = my_proj_name ,
                                           binary.meth = "TSS",
                                           do.stack = FALSE )
  
})

#model_path <- paste0("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/", name_model_folder)

#Ensemble model
lapply(spp_list[1:length(spp_list)], function(i) {
  
  setwd(global_model_path)
  
  #i <- gsub(".", "_", i)
  
  #Load separate models
  model_name <- file.path(global_model_path, i, paste0(i,".",model_class,".models.out"))
  model_i <- load(model_name)
  
  #Load separate projection models
  model_proj <- file.path(global_model_path, i, paste0("proj_", model_class), paste0(i,".",model_class,".projection.out"))
  model_pi <- load(model_proj)
  
  # Fit ensemble model
  all_test <- get_evaluations(get(model_i))
  all_tss <- all_test[2, 1, 1:5, 1:4, 1:3]
  all_tss[is.na(all_tss)] <- 0 #if there are any NAs the following won't work properly
  
  if(length(all_tss[all_tss>=0.7]) != 0) {
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), 
                                                  em.by='all',
                                                  eval.metric = 'TSS', 
                                                  eval.metric.quality.threshold = 0.7, 
                                                  models.eval.meth = c('KAPPA','TSS','ROC'),
                                                  prob.mean=FALSE, 
                                                  prob.cv=FALSE, 
                                                  committee.averaging=TRUE, 
                                                  prob.mean.weight = FALSE, 
                                                  VarImport = 3)
  }
  
  else {
    val_tss <- stats::quantile(all_tss, probs = 0.9,names=FALSE,na.rm=TRUE)
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), 
                                                  em.by='all',
                                                  eval.metric = 'TSS', 
                                                  eval.metric.quality.threshold = val_tss, 
                                                  models.eval.meth = c('KAPPA','TSS','ROC'), 
                                                  prob.mean=FALSE, 
                                                  prob.cv=FALSE,
                                                  committee.averaging=TRUE, 
                                                  prob.mean.weight = FALSE, 
                                                  VarImport = 3)
  }
  
  # Project ensemble model
  models_ensemble_proj_current <- BIOMOD_EnsembleForecasting(EM.output = all_ensemble_model,
                                                             projection.output = get(model_pi),
                                                             binary.meth = "TSS",
                                                             do.stack = FALSE )
})


#Accuracy
lapply(spp_list[1:length(spp_list)], function(i) {
  
  #model_path <- paste0("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/", name_model_folder)
  setwd(global_model_path)
  
  #Load separated models
  model_name <- file.path(global_model_path, i, paste0(i,".",model_class,".models.out"))
  model_i <- load(model_name)
  
  #Load unique ensemble model
  model_ensem <- file.path(global_model_path, i, paste0(i,".",model_class,"ensemble.models.out"))
  model_ei <- load(model_ensem)
  
  # Get evaluations of ensemble models to decide which one to choose:
  all_ensemble_models_scores <- get_evaluations(get(model_ei))
  write.table(all_ensemble_models_scores, file = file.path(global_model_path, i, paste0(i,"_eval_ensemble.csv")))
  
  # Get variables importance for all models
  all_models_var_import <- get_variables_importance(get(model_i))
  
  # Calculate the mean of variable importance by algorithm
  table_importance <- apply(all_models_var_import, c(1,2), mean, na.rm=TRUE)
  write.table(table_importance, file = file.path(global_model_path, i, paste0(i,"_importance_var.csv")))
  
})
  
                                            ###Australian predictions
#worldclim - use already loaded "bio"
Aus_bio_2.5 <- crop(bio, Aus_raster)

#Elevation
Aus_elev <- raster("C:/Users/david/Documents/PhD/Chapter_3/SpatialData/Raster/Elevation/AUS_msk_alt.gri")
Aus_elev <- crop(Aus_elev, Aus_Coast)
Aus_elev <- mask(Aus_elev, Aus_Coast)
Aus_elev_2.5 <- resample(Aus_elev, Aus_bio_2.5[[1]], method = "bilinear") 

#Vegetation
Aus_veg <- raster("C:/Users/david/Documents/PhD/Chapter_3/SpatialData/Raster/Aus_veg.gri")

# reclassify
Aus_veg[Aus_veg %in% c(1:4, 30)] <- 1 #1, 2,3,4, 30 - Open forests/Rainforests and Vine Thickets
Aus_veg[Aus_veg %in% c(5:13, 31:32)] <- 2 #5:13, 31,32 - Woodlands
Aus_veg[Aus_veg %in% c(14:17)] <- 3 #14:17 - Shrublands
Aus_veg[Aus_veg == 18] <- 4 #18 - Heathlands
Aus_veg[Aus_veg %in% c(19:22)] <- 5 #19:22 - Grasslands
Aus_veg[Aus_veg %in% c(23,24,28)] <- 6 #23,24, 28 - Aquatic
Aus_veg[Aus_veg %in% c(26,29)] <- 7 #29, 26 - Native vegetation
Aus_veg[Aus_veg %in% c(25,27, 33)] <- 8 #25, 27 - Disturbed/bare/unknown

# set new RAT for reclassfied raster
Aus_veg <- raster::ratify(Aus_veg)
rat <- data.frame(
  ID = 1:8,
  landcover = c("Open forests/Rainforests and Vine Thickets", 
                "Woodlands",
                "Shrublands",
                "Heathlands",
                "Grasslands",
                "Aquatic",
                "Native vegetation",
                "Disturbed/bare/unknown")
)
levels(Aus_veg) <- rat
Aus_veg_2.5 <- resample(Aus_veg, Aus_bio_2.5[[1]], method = "ngb")
Aus_veg_2.5 <- mask(Aus_veg_2.5, Aus_bio_2.5[[1]])

#Human footprint
impact <- raster("C:/Users/david/Documents/PhD/Chapter_3/SpatialData/Raster/wildareas-v3-2009-human-footprint.tif")
Aus_Coast_moll <- st_transform(Aus_Coast, crs = "+proj=moll")
aus_impact <- crop(impact, Aus_Coast_moll)
aus_impact <- projectRaster(from = aus_impact, res = 0.008333333, crs = "+proj=longlat +datum=WGS84 +no_defs")
aus_impact <- mask(aus_impact, Aus_Coast)
aus_impact_2.5 <- resample(aus_impact, Aus_bio_2.5[[1]])

#Stacking all
all_predictors <- stack(Aus_bio_2.5, Aus_elev_2.5, Aus_veg_2.5, aus_impact_2.5)

#Reducing collinearity
reg_nocorrvar <- vifstep(all_predictors, th = 4)
reg_predictors <- as.character(reg_nocorrvar@results$Variables)
pred_au_sub <- raster::subset(all_predictors,c(reg_predictors))

###Set parameters for regional model
name_model_folder <- "IAS_regional" 
model_class <- "regional" 
regional_model_path <- "C:/Users/david/Documents/PhD/Chapter_3/SpatialData/IAS_distributions/IAS_regional"

#Individual models
my_models <- c("GLM","GAM","MAXENT.Phillips","FDA","GBM") # Algorithms to fit the models   
my_pa_strategy <- c("random") # Strategy to select PAs   
my_pa_dist_min <- 0 # Distance where to start throwing PAs  
my_runs <- 4 # Number of PAs runs
my_n_pa <- 5000 # Number of pseudoabsences, when not dependent of no. presences
my_proj_name <- "regional"


lapply(spp_list[1:length(spp_list)], function(i) {
  
  #Remove decimal in species name
  i <- gsub("\\.", "_", i)
  
  #Load presences
  #sp_pres <- readOGR(file.path(occ_path,paste0(i,"_pres.shp"))) #try with pres_all.shp
  sp_pres <- st_read(file.path(occ_path,paste0(i,"_pres_all.shp"))) #try with pres_all.shp.
  #used all for M. floricola and M. destructor
  
  #Remove decimal in species name
  i <- gsub("_", ".", i)
  
  #Subset of Australia
  sp_pres_au <- rasterize(sp_pres, Aus_raster, 1)
  n_mod <- summary(as.factor(sp_pres_au@data@values))[1]
  
  #Pseudoabsences
  # Select number of pseudoabsences
  n_pa <- my_n_pa
  
  # Weights based on the global model:
  all_pa <- raster(file.path(global_model_path, i, "proj_global", "individual_projections", paste0(i, "_EMcaByTSS_mergedAlgo_mergedRun_mergedData.grd")))
  pa_au <- raster::crop(all_pa,Aus_raster)
  pa_au <- raster::mask(pa_au,Aus_raster)
  pa_au_1 <- pa_au/1000 # conversion from 0 to 1
  weights_pa <- 1/(1+((pa_au_1/(pa_au_1-1))^2)) # inverse logistic transformation
  
  # Format data for 'BIOMOD'
  # Subset of only presences
  subset_1 <- as.data.frame(rasterToPoints(sp_pres_au,fun=function(x){x==1}))
  
  # Presences
  mypres <- rep(1,dim(subset_1)[1])
  
  # Presences coordinates
  mycoord <- subset_1[,c("x","y")]
  
  # Clean data
  clean_data <- BIOMOD_FormatingData(resp.var=mypres, 
                                     expl.var=pred_au_sub, 
                                     resp.xy=mycoord,
                                     resp.name=i, 
                                     PA.nb.rep=my_runs,
                                     PA.nb.absences=n_pa,
                                     PA.strategy = my_pa_strategy, 
                                     PA.dist.min=my_pa_dist_min)
  
  # Adding "weights" to the data
  # Get the presences + PA dataset
  my_pres_PA_df <- data.frame(clean_data@coord, 
                              obs = clean_data@data.species, 
                              clean_data@PA)
  
  # Add the weight vector to the presences
  pres <- subset(my_pres_PA_df, obs==1)
  pres$yweights <- 1
  
  # Weight of absences based on global model
  abs <- subset(my_pres_PA_df, is.na(obs)==TRUE)
  weight_values_ab <- raster::extract(weights_pa,cbind(abs$x,abs$y))
  abs$yweights <- weight_values_ab
  
  new_my_pres_PA_df <- rbind(pres,abs)
  
  # Settings for the model algorithms
  model_opt <- BIOMOD_ModelingOptions(GLM=list(type='quadratic', interaction.level=0),
                                      GBM=list(n.trees=1000),
                                      MAXENT.Phillips = list(path_to_maxent.jar = maxent_jar_path,product=FALSE), #change path
                                      GAM=list(k=3))
  
  #Set path for model folder
  #model_path <- paste0("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/", name_model_folder)
  setwd(regional_model_path)
  
  # Run the models
  # Separate
  all_models <- BIOMOD_Modeling(data=clean_data, 
                                models=my_models, 
                                models.options=model_opt, 
                                NbRunEval = 4, 
                                DataSplit = 70, 
                                VarImport = 3, 
                                do.full.models = F, 
                                modeling.id = "regional", 
                                Yweights=new_my_pres_PA_df$yweights)
  
  # Project the models
  # Separate
  models_proj_current <- BIOMOD_Projection(modeling.output = all_models,
                                           new.env = pred_au_sub,
                                           proj.name = my_proj_name ,
                                           binary.meth = "TSS",
                                           do.stack = FALSE )
})


#Ensemble model
lapply(spp_list[1:length(spp_list)], function(i) {
  
  #Load separate models
  model_name <- file.path(regional_model_path, i, paste0(i,".",model_class,".models.out"))
  model_i <- load(model_name)
  
  #Load separate projection models
  model_proj <- file.path(regional_model_path, i, paste0("proj_", model_class), paste0(i,".",model_class,".projection.out"))
  model_pi <- load(model_proj)
  
  # Fit ensemble model
  all_test <- get_evaluations(get(model_i))
  all_tss <- all_test[2, 1, 1:5, 1:4, 1:3]
  all_tss[is.na(all_tss)] <- 0 #if there are any NAs the following won't work properly
  
  if (length(all_tss[all_tss>=0.7]) != 0) {
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), 
                                                  em.by='all',
                                                  eval.metric = 'TSS', 
                                                  eval.metric.quality.threshold = 0.7, 
                                                  models.eval.meth = c('KAPPA','TSS','ROC'),
                                                  prob.mean=FALSE, 
                                                  prob.cv=FALSE, 
                                                  committee.averaging=TRUE, 
                                                  prob.mean.weight = FALSE, 
                                                  VarImport = 3)
  }
  
  else {
    val_tss <- stats::quantile(all_tss, probs = 0.9,names=FALSE,na.rm=TRUE)
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), 
                                                  em.by='all',
                                                  eval.metric = 'TSS', 
                                                  eval.metric.quality.threshold = val_tss, 
                                                  models.eval.meth = c('KAPPA','TSS','ROC'), 
                                                  prob.mean=FALSE, 
                                                  prob.cv=FALSE,
                                                  committee.averaging=TRUE, 
                                                  prob.mean.weight = FALSE, 
                                                  VarImport = 3)
  }
  
  # Project ensemble model
  models_ensemble_proj_current <- BIOMOD_EnsembleForecasting(EM.output = all_ensemble_model,
                                                             projection.output = get(model_pi),
                                                             binary.meth = "TSS",
                                                             do.stack = FALSE )
  
})

#Accuracy
lapply(spp_list[1:length(spp_list)], function(i) {
  
  #model_path <- paste0("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/", name_model_folder)
  setwd(regional_model_path)
  
  #Load separated models
  model_name <- file.path(regional_model_path, i, paste0(i,".",model_class,".models.out"))
  model_i <- load(model_name)
  
  #Load unique ensemble model
  model_ensem <- file.path(regional_model_path, i, paste0(i,".",model_class,"ensemble.models.out"))
  model_ei <- load(model_ensem)
  
  # Get evaluations of ensemble models to decide which one to choose:
  all_ensemble_models_scores <- get_evaluations(get(model_ei))
  write.table(all_ensemble_models_scores, file = file.path(regional_model_path, i, paste0(i,"_eval_ensemble.csv")))
  
  # Get variables importance for all models
  all_models_var_import <- get_variables_importance(get(model_i))
  
  # Calculate the mean of variable importance by algorithm
  table_importance <- apply(all_models_var_import, c(1,2), mean, na.rm=TRUE)
  write.table(table_importance, file = file.path(regional_model_path, i, paste0(i,"_importance_var.csv")))
  
})


### Results
# Here you decide which results you want the output from. Maybe I do it for both global and regional
name_model_folder <- "IAS regional"
dataset_type <- "_pres" # Options: _pres/ _pres_all ## from GBIF
global <- "no"
model_class <- "regional"
raster_to_use <- Aus_raster 
#model_path <- paste0("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/", name_model_folder)
setwd(regional_model_path)

## Sample size (n)
# empty vector to store n for all species
col_n <- c()

#lapply(spp_list[1:length(spp_list)], function(i) {
for(i in spp_list){
  
  if (global=="no"){raster1 <- Aus_raster}
  else {raster1 <- ref_raster}
  
  #replace decimal
  i <- gsub("\\.", "_", i)
  
  # GBIF points
  #sp_pres <- readOGR(file.path(occ_path,paste0(i,dataset_type,".shp")))
  sp_pres <- st_read(file.path(occ_path,paste0(i,dataset_type,".shp")))
  
  # rasterize
  points_mod_r <- rasterize(sp_pres,raster1,1,background=0)
  
  n_mod <- summary(as.factor(points_mod_r@data@values))["1"]
  
  # Create a vector with all the species:
  col_n <- c(col_n,n_mod)
  }
#}) 

## TSS, Sensitivity, Specificity, Cut-off (for binary transformation)
# empty vectors to store information for all species
col_tss <- c()
col_sens <- c()
col_spec <- c()
col_cut <- c()

for(i in spp_list){
  
  # File with all the info. per species:
  model_eval <- read.csv(file = file.path(regional_model_path, i, paste0(i,"_eval_ensemble.csv")), sep = " ")
  
  # selecting the parameters that I want
  tss_mod <- round(as.numeric(as.character(model_eval[2,1])),2)
  sens_mod <- round(as.numeric(as.character(model_eval[2,3])),2)
  spec_mod <- round(as.numeric(as.character(model_eval[2,4])),2)
  cut_mod <- round(as.numeric(as.character(model_eval[2,2]))/1000,2)
  
  # Create vectors with all the species:
  col_tss <- c(col_tss, tss_mod)
  col_sens <- c(col_sens,sens_mod)
  col_spec <- c(col_spec,spec_mod)
  col_cut <- c(col_cut, cut_mod)
  
}


## Final table with info for all: 
table_a <- as.data.frame(cbind(as.character(spp_list),col_n,col_tss,col_sens,col_spec,col_cut))
names(table_a) <- c("Species","n","TSS","Sensitivity","Specificity","Cut-off binary")

# Write in the directory:
write.table(table_a, 
            file = file.path(regional_model_path, paste0("Table_",name_model_folder,"_accuracy.txt")),
            row.names = FALSE)



