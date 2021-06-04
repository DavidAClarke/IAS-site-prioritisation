############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "runAnalysis.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts","runAnalysis.R"),
              local_path = file.path("Scripts", "runAnalysis.R"), overwrite = T)


################################ Main script for running the analyses #########################

# Clear the environment
rm(list = ls())

# Load required libraries
library(rdrop2)
library(tidyverse)
library(readxl)
library(rvest)
library(purrr)
library(sf)
library(sp)
library(tmap)
library(maps)
library(raster)
library(gdalUtils)
library(rgbif)
library(CoordinateCleaner)
library(taxize)
library(traitdataform)
library(red)
library(letsR)
library(biomod2)
library(zonator)
library(adehabitatHR)

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
feature_rst(Aus_KBA, "KBAs", "Name")

#Red List species
feature_rst(Aus_RL_IAS, "RedList_species_IAS", "scientificName")
feature_rst(Aus_RL, "RedList_species_all", "scientificName")

#Ramsar and Upstream catchments
feature_rst(Aus_Rams, "Ramsar", "RAMSAR_NAM")
feature_rst(Aus_Ups, "Upstream", "RAM_NAME")

#Red List Ecosystems
feature_rst(Aus_Ecosystems, "RedList_Ecosystems", "Ecosystem")

#Protected areas
feature_rst(Aus_CAPAD, "PAs", "NAME")

#SNES
feature_rst(Aus_SNES, "SNES_all", "scientificName")
feature_rst(Aus_SNES_IAS, "SNES_IAS", "scientificName")

                                  # Copy files into folders for zonation #

                                           # Running Zonation #
