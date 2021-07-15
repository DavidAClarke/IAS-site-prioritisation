# IAS SDM - Based off Polaina et al 2020

#Required libraries
library(sf)
library(sp)
library(raster)
library(data.table)
library(rgbif)
library(biomod2)
library(countrycode)
library(CoordinateCleaner)
library(usdm)

#Required initial data
#bio <- getData("worldclim", var = "bio", res = 2.5, path = "E:/SpatialData/Raster/Worldclim")
bio_names <- list.files(path = "E:/SpatialData/Raster/Worldclim/wc2-5", pattern = ".bil$", full.names = T)
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

#IAS insect species list
spp_list <- c("Vespa velutina", "Bombus impatiens", "Lymantria dispar", "Harmonia axyridis", 
              "Solenopsis invicta", "Wasmannia auropunctata", "Glischrochilus quadrisignatus",
              "Digitonthophagus gazella", "Icerya purchasi", "Pheidole megacephala")

#Download occurrence data
#GBIF credentials (remove credentials when publishing scripts)
# user <- "davidclarke" 
# pwd <- "!8208GBIF153"
# email <- "david.clarke1@monash.edu"

#Match the names
# gbif_taxon_keys <- spp_list %>% 
#   taxize::get_gbifid_(method = "backbone") %>%
#   imap(~ .x %>% mutate(original_sciname = .y)) %>%
#   bind_rows() %T>%
#   write.csv(file = file.path("SpeciesData", "IAS_all_matches.csv")) %>% #can then examine the data
#   dplyr::filter(matchtype == "EXACT" & status == "ACCEPTED") %>%
#   pull(usagekey)

#Request download
# res <- occ_download(
#   pred_in("taxonKey", gbif_taxon_keys),
#   pred_in("basisOfRecord", c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION', 'PRESERVED_SPECIMEN')),
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SIMPLE_CSV",
#   user=user,pwd=pwd,email=email)

#Download data
occ_path <- file.path("SpatialData", "Vector", "IAS_Occurrences")
#occ_download_get("0324249-200613084148143", path = occ_path)

#Occurrence data path
data_path <- paste0(getwd(), "/", file.path("SpatialData", "Vector", "IAS_Occurrences", "0324249-200613084148143.csv"))

#Maybe combine loading and cleaning functions
#Load occurrence data, convert to spatial object, write to shapefile
spp_list <- gsub(" ", ".", spp_list)

my_cols <- fread(cmd = paste("head -n 1", data_path))

my_cols_min <- c("family","species", "taxonRank","scientificName","decimalLongitude","decimalLatitude","countryCode","stateProvince", "coordinateUncertaintyInMeters", "coordinatePrecision","basisOfRecord","year", "individualCount", "institutionCode")

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
           dsn = file.path(occ_path, paste0(i,"_data.shp")))
           
})

#Clean occurrence data

radius_uncert <- 2500 #radius of cell size: 2.5 arc min very ~ 5km
  
lapply(spp_list[1:length(spp_list)], function(i) {
  
  #Repalce decimal in species name
  i <- gsub("\\.", "_", i)
  
  #Load shapefile
  occ_shp <- readOGR(file.path(occ_path,paste0(i,"_data.shp")))
  
  #Removing those with uncertainty > radius_uncert and NAs
  occ_shp_cert <- subset(occ_shp, crdnUIM<= radius_uncert)
  
  
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
  
  # Converting to a spatial object - OCURRENCES
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
           dsn = file.path(occ_path,paste0(i,"_occ.shp")))
  
  st_write(pres_sf, 
           dsn = file.path(occ_path,paste0(i,"_pres.shp")))
  
  
  #Keeping the NAs as well as those with "certain" data
  occ_shp_na <- occ_shp[is.na(occ_shp$crdnUIM),]
  
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
           dsn = file.path(occ_path,paste0(i,"_occ_all.shp")))
  
  st_write(pres_sf_2, 
           dsn = file.path(occ_path,paste0(i,"_pres_all.shp")))
  
})
###########################################################################################
#PREDICTORS - Polaina et al 2020 only use climatic variables for global model

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
my_runs <- 3 # Number of PAs runs
my_n_pa <- 20000 # Number of pseudoabsences, when not dependent of no. presences
my_proj_name <- "global"  # Name for the output projections

lapply(spp_list[1:length(spp_list)], function(i){
  
  #Remove spaces in species name
  i <- gsub("//.", "_", i)
  
  #Load in presence shapefile
  sp_pres <- readOGR(file.path(occ_path,paste0(i,"_pres_all.shp")))
  n_mod <- dim(sp_pres)[1]
  
  #Pseudoabsences
  n_pa <- my_n_pa
  
                                ## Format data for 'BIOMOD'
  #Presences
  mypres <- rep(1,dim(sp_pres)[1])
  
  #Presences coordinates
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
  maxent_jar_path <- "C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/maxent.jar" #could change. No spaces allowed
  model_opt <- BIOMOD_ModelingOptions(GLM=list(type='quadratic', interaction.level=0),
                                      GBM=list(n.trees=1000),
                                      MAXENT.Phillips = list(path_to_maxent.jar = maxent_jar_path, product=FALSE),
                                      GAM=list(k=3))
  
  #Path for models
  #must be full path (back to e.g. C:, and have no spaces)
  model_path <- paste0("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/", name_model_folder)
  setwd(model_path)
  
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

#Ensemble model
start_time <- Sys.time()
lapply(spp_list[1:length(spp_list)], function(i) {
  
  #i <- gsub(".", "_", i)
  
  #Load separate models
  model_name <- file.path(model_path, i, paste0(i,".",model_class,".models.out"))
  model_i <- load(model_name)
  
  #Load separate projection models
  model_proj <- file.path(model_path, i, paste0("proj_", model_class), paste0(i,".",model_class,".projection.out"))
  model_pi <- load(model_proj)
  
  # Fit ensemble model
  all_test <- get_evaluations(get(model_i))
  all_tss <- all_test[2, 1, 1:5, 1:4, 1:3]
  
  if (all_tss>=0.7) {
    all_ensemble_model <- BIOMOD_EnsembleModeling(modeling.output = get(model_i), 
                                                  em.by='all',
                                                  eval.metric = 'TSS', 
                                                  eval.metric.quality.threshold = 0.7, 
                                                  models.eval.meth = c('KAPPA','TSS','ROC'),
                                                  prob.mean=FALSE, 
                                                  prob.cv=TRUE, 
                                                  committee.averaging=TRUE, 
                                                  prob.mean.weight = TRUE, 
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
                                                  prob.cv=TRUE,
                                                  committee.averaging=TRUE, 
                                                  prob.mean.weight = TRUE, 
                                                  VarImport = 3)
  }
  
  # Project ensemble model
  models_ensemble_proj_current <- BIOMOD_EnsembleForecasting(EM.output = all_ensemble_model,
                                                             projection.output = get(model_pi),
                                                             binary.meth = "TSS",
                                                             do.stack = FALSE )
})
end_time <- Sys.time()
end_time - start_time
  
  

