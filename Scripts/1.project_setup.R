                                          #Project setup
##Libraries
library(tidyverse)
library(readxl)
library(sf)
library(sp)
library(maps)
library(raster)
library(MASS)
library(bossMaps)
#library(dismo) #not to run but to evaluate maxent SDMs
library(rJava)
library(rgbif)
library(CoordinateCleaner)
#library(galah)
#library(zonator)
library(biomod2)
library(rgeos)
#library(ENMTools)
library(data.table)
library(countryCode)
library(rgdal)
library(usdm)

##required to be installed
#"taxize"
#"traitdataform"


##Load xy_match function
xy_match <- function(initial, pre, Layer, Long, Lat) {
  pre <- initial
  
  coordinates(pre) <- ~ Long + Lat
  crs(pre) <- crs(Layer)
  
  final <- cbind(initial, over(pre, Layer))
  return(final)
}

##Load initially required data (other data can maybe come in later)
source("Scripts/2.aus_coast.R")
source("Scripts/3.env_predictors.R")
source("Scripts/4.redlist_shp.R")
source("Scripts/4.redlist_pts.R")
source("Scripts/5.calc_range_offset.R")
source("Scripts/6.plant_point_bias_files.R")

##Creating folder structure for running Maxent
#Specify maxent path
#dir.create("Maxent") #direct the installation of maxent to here.
maxent_path <- file.path("Maxent")

#Create necessary folders
dir.create(file.path(maxent_path, "Env_layers")) #species specific 
dir.create(file.path(maxent_path, "Bias_layers")) #both range offsets and density maps
dir.create(file.path(maxent_path, "Occurrences")) #cleaned occurrence records
dir.create(file.path(maxent_path, "Outputs")) #maxent outputs

#Create necessary subfolders for each species
#speciesNames = final species set
lapply(speciesNames[1:length(speciesNames)], function(i){
  
  i <- gsub(" ", "_", i)
  
  dir.create(file.path(maxent_path, "Env_layers", i))
  dir.create(file.path(maxent_path, "Bias_layers", i))
  dir.create(file.path(maxent_path, "Occurrences", i))
  dir.create(file.path(maxent_path, "Outputs", i))
  
})

#Create paths
env_path <- file.path(maxent_path, "Env_layers")
bias_path <- file.path(maxent_path, "Bias_layers")
occ_path <- file.path(maxent_path, "Occurrences")


##Creating folder structure for running Zonation
#Create path to a new project
dir.create("Zonation")
zonation_path <- file.path("Zonation")
dir.create(project_path, "Input_rasters")