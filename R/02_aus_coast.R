                        ## Loading and pre-processing of coastal shapefile ##

#Load spatial data
Coast_shp <- st_read("E:/SpatialData/Vector/coastal-gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp")


#Only keep mainland and Tasmania
Aus_Coast <- Coast_shp %>% dplyr::filter(id == 6 | id == 32)
Aus_Coast <- st_make_valid(Aus_Coast)

#Removing original shapefile
rm(Coast_shp)

#Projecting (GDA94 = 3577)
#Aus_Coast_proj <- st_transform(Aus_Coast, 3577)

#Removing unprojected shapefile
#rm(Aus_Coast)
