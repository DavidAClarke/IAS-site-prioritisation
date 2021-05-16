############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "runAnalysis.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "runAnalysis.R"),
              local_path = file.path("Scripts", "runAnalysis.R"))


################################ Main script for running the analyses #########################

# Clear the environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(readxl)
library(rvest)
library(purrr)
library(sf)
library(tmap)
library(raster)
library(gdalUtils)
library(rgbif)
library(CoordinateCleaner)
library(taxize)
library(red)
library(letsR)
library(biomod2)
library(zonator)
library(rdrop2)

# Load bespoke functions
source("Scripts/Functions.R")

                              # Cleaning and pre-processing of input data #
# Australian coast
source("Scripts/Coast.R")

# IUCN Red List
source("Scripts/RedList.R")

# Species of National Environmental Significance
source("Scripts/SNES.R")

# Key Biodiversity Areas
source("Scripts/KBA.R")

# Ramsar and upstream catchments
source("Scripts/RamsUp.R")

# Protected Areas
source("Scripts/CAPAD.R")

# Red List of Ecosystems
source("Scripts/RedListEco.R")

                                    # Summary statistics of input data #


                                        # Create feature rasters #
#Create folder to send raster files for input into Zonation
#maybe have folders fo the different data e.g. KBA, RL_ALL, PA, etc, instead of what I currently have planned
#then use file.copy to copy rasters into the RL_IAS, RL_ALL,etc folders
#dir.create(file.path("SpatialData", "Input_zonation", "RL_IAS"))
#dir.create(file.path("SpatialData", "Input_zonation", "RL_ALL"))
#dir.create(file.path("SpatialData", "Input_zonation", "SNES_IAS"))
#dir.create(file.path("SpatialData", "Input_zonation", "SNES_ALL"))

#Key Biodiversity Areas (DONE)
#feature_rst(Aus_KBA, "RL_IAS", "Name")
feature_rst(Aus_KBA, "RL_ALL", "Name")
feature_rst(Aus_KBA, "SNES_IAS", "Name")
feature_rst(Aus_KBA, "SNES_ALL", "Name")

#Red List species
#feature_rst(Aus_RL_IAS, "RL_IAS", "scientificName")
feature_rst(Aus_RL, "RL_ALL", "scientificName")

#Ramsar and Upstream catchments
#feature_rst(Aus_Rams, "RL_IAS", "RAMSAR_NAM")
#feature_rst(Aus_Ups, "RL_IAS", "RAM_NAME")
feature_rst(Aus_Rams, "RL_ALL", "RAMSAR_NAM")
feature_rst(Aus_Ups, "RL_ALL", "RAM_NAME")
feature_rst(Aus_Rams, "SNES_IAS", "RAMSAR_NAM")
feature_rst(Aus_Ups, "SNES_IAS", "RAM_NAME")
feature_rst(Aus_Rams, "SNES_ALL", "RAMSAR_NAM")
feature_rst(Aus_Ups, "SNES_ALL", "RAM_NAME")

#Red List Ecosystems
feature_rst(Aus_Ecosystems, "RL_IAS", "Ecosystem")

#Protected areas
feature_rst(Aus_CAPAD, "PAs", "NAME")

#SNES
feature_rst(SNES_shp, "SNES", "scientificName")

                                           # Running Zonation #
