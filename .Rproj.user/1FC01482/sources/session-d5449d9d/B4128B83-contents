#Spatial layer plot (https://www.urbandemographics.org/post/figures-map-layers-r/)
##################################################################################
library(tidyverse)
library(sf)
library(raster)
library(stars)
library(ggnewscale)
library(zonator)

#######################################################################################
#Functions to tilt sf
rotate_data <- function(data, x_add = 0, y_add = 0) {
  
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  
  rotate_matrix <- function(x){ 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geometry = .$geometry * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}

rotate_data_geom <- function(data, x_add = 0, y_add = 0) {
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  
  rotate_matrix <- function(x) { 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geom = .$geom * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}
################################### Priority maps ######################
#Load priority rank rasters

#Convert raster to stars (st_as_stars())
CAZ_st <- st_as_stars(CAZ_var_ras)
CAZ_wgt_st <- st_as_stars(CAZ_wgt_var_ras)
CAZ_wgt_KBA_st <- st_as_stars(CAZ_wgt_KBA_var_ras)

CAZ_area_st <- st_as_stars(CAZ_area_var_ras)
CAZ_area_wgt_st <- st_as_stars(CAZ_area_wgt_var_ras)
CAZ_area_wgt_KBA_st <- st_as_stars(CAZ_area_wgt_KBA_var_ras)

#Convert to sf (st_as_sf())
CAZ_sf <- st_as_sf(CAZ_st)
CAZ_wgt_sf <- st_as_sf(CAZ_wgt_st)
CAZ_wgt_KBA_sf <- st_as_sf(CAZ_wgt_KBA_st)

CAZ_area_sf <- st_as_sf(CAZ_area_st)
CAZ_area_wgt_sf <- st_as_sf(CAZ_area_wgt_st)
CAZ_area_wgt_KBA_sf <- st_as_sf(CAZ_area_wgt_KBA_st)

x <- 182
color <- "black"
leg <- zlegend("spectral")
#xaxis goes from 180 to 280; extend to 300 to include text next to map
##########################################################################################
#Make plot
#NOTE: first map made will be on the bottom
speciesMaps <- ggplot() +
  geom_sf(
    data = 
      CAZ_wgt_KBA_sf %>% 
      rotate_data(),
    aes(fill=species_wgt_CAZ_KBA.CAZ_ME.rank.compressed), 
    color=NA, 
    show.legend = FALSE
  ) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values) +
  annotate("text", 
           label='Species + weights\n+ KBA', 
           x=x, y= -75.0, hjust = 0, color=color, size = 10) +
  new_scale_fill() + 
  new_scale_color() +
  geom_sf(
    data = 
      CAZ_wgt_sf %>% 
      rotate_data(y_add = 35),
    aes(fill=species_wgt_CAZ.CAZ_E.rank.compressed), 
    color=NA, 
    show.legend = FALSE
  ) +
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values) +
  annotate("text", 
           label='Species + weights', #\n puts it on the next line
           x=x, y= -35.0, hjust = 0, color=color, size = 10) +
  new_scale_fill() + 
  new_scale_color() +
  geom_sf(
    data = 
      CAZ_sf %>% 
      rotate_data(y_add = 70),
    aes(fill=species_CAZ.CAZ_E.rank.compressed), 
    color=NA, 
    show.legend = FALSE
  ) +
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values) +
  annotate("text", 
           label='Species', #\n puts it on the next line
           x=x, y= 0.0, hjust = 0, color=color, size = 10) +
  theme_void()
###########################################################################
species_areaMaps <- ggplot() +
  geom_sf(
    data = 
      CAZ_area_wgt_KBA_sf %>% 
      rotate_data(),
    aes(fill=species_area_wgt_CAZ_KBA.CAZ_E.rank.compressed), 
    color=NA, 
    show.legend = FALSE
  ) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values) +
  annotate("text", 
           label='Species + area\n+ weights + KBA', 
           x=x, y= -75.0, hjust = 0, color=color, size = 10) +
  new_scale_fill() + 
  new_scale_color() +
  geom_sf(
    data = 
      CAZ_area_wgt_sf %>% 
      rotate_data(y_add = 35),
    aes(fill=species_area_wgt_CAZ.CAZ_E.rank.compressed), 
    color=NA, 
    show.legend = FALSE
  ) +
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values) +
  annotate("text", 
           label='Species +area\n+ weights', #\n puts it on the next line
           x=x, y= -35.0, hjust = 0, color=color, size = 10) +
  new_scale_fill() + 
  new_scale_color() +
  geom_sf(
    data = 
      CAZ_area_sf %>% 
      rotate_data(y_add = 70),
    aes(fill=species_area_CAZ.CAZ_E.rank.compressed), 
    color=NA, 
    show.legend = FALSE
  ) +
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values) +
  annotate("text", 
           label='Species + area', #\n puts it on the next line
           x=x, y= 0.0, hjust = 0, color=color, size = 10) +
  theme_void()

