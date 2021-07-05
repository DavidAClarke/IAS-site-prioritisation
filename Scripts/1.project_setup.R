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
library(dismo) #not to run but to evaluate maxent SDMs
library(rJava)
library(rgbif)
library(CoordinateCleaner)
library(galah)
library(zonator)
library(biomod2)

##required to be installed
#"taxize"
#"traitdataform"

##Load initially required data (other data can maybe come in later)
source("Scripts/Coast.R")
source("Scripts/env_predictors.R")
#source("Scripts/RedList.R") #needs to be tidied up
speciesNames <- read.csv(file.path("SpeciesData", "all_accepted_names.csv"))
speciesNames <- speciesNames$x


##Creating folder structure for running Maxent
#Specify maxent path
dir.create("Maxent") #direct the installation of maxent to here.
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
  
  dir.create(file.path(maxent_path, "env_layers", i))
  dir.create(file.path(maxent_path, "bias_layers", i))
  dir.create(file.path(maxent_path, "occurrences", i))
  dir.create(file.path(maxent_path, "outputs", i))
  
})

##Creating folder structure for running Zonation
#Create path to a new project
dir.create("Zonation")
zonation_path <- file.path("Zonation")
dir.create(project_path, "Input_rasters")