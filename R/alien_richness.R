library(letsR)

# Alien species in Australia

#Load and combine cleaned occurrences
# files <- list.files("G:/Alien_Aus_GBIF", pattern = "final", full.names = T)
# 
# final <- data.frame(num = numeric(),
#                     species = character(),
#                     decimalLongitude = numeric(),
#                     decimalLatitude = numeric())
# write.csv(final, file = "E:/Alien_Aus_GBIF/Final.csv")
# 
# lapply(files[1:length(files)], function(i) {
#   
#   f <- read.csv(i)
#   write.table(f, file = "E:/Alien_Aus_GBIF/Final.csv", append = T, col.names = F, sep = ",")
#   
# })

#Load points
alien_points <- read.csv("G:/Alien_Aus_GBIF/Final.csv")
occ_mat <- cbind(alien_points$decimalLongitude, alien_points$decimalLatitude)
alien_presab <- lets.presab.points(xy = occ_mat,
                                   species = alien_points$species,
                                   xmn = 113.1667,
                                   xmx = 153.625,
                                   ymn = -43.625,
                                   ymx = -10.70833,
                                   resol = 0.5,#0.04166667
                                   count = T)
r <- alien_presab$Richness_Raster
r <- resample(r, PM_bin)
r <- mask(r, PM_bin)
r <- round(r)
writeRaster(r, filename = file.path("SpatialData", "Raster", "alien_richness.tif"), overwrite = T)
plot(r, col = rev(bpy.colors(60)))

#or
rdf <- raster::as.data.frame(r, xy = T)
rdf <- na.omit(rdf)
ggplot() +
  geom_raster(data = rdf , aes(x = x, y = y, fill = alien_richness)) +
  scale_fill_viridis_c(option = "C") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())
