################################################################################
## Script name: 00_main_menu.R
################################################################################

## Load libraries
pkgs <- c("tidyverse", "sp", "sf", "raster", "data.table", "CoordinateCleaner",
          "countrycode", "MASS", "bossMaps", "rJava", "ENMTools", "zonator",
          "taxize", "biomod2", "usdm", "rgbif", "galah", "spatstat", "stars", 
          "ggpubr","cowplot", "RColorBrewer", "PNWColors")
lapply(pkgs, require, character.only = T)


## File paths
shared_data <- here("C:/Users/name/Documents/postdoc/projects/shared_data/")

## User defined functions
source("R/01_functions.R")

# Look at all scripts between 01 & 17

## All results
source("R/17_results.R")