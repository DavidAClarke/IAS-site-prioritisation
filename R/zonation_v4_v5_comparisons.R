### Testing script

mycols <- grDevices::colorRampPalette(c("#EE602D", "#FFCF26", "#01AEE6"))

z_colors_spectral <- list(values=c(0.0, 0.2, 0.5, 0.75, 0.9, 0.95, 0.98, 1.0),
                          labels=c("0.00-0.20", "0.20-0.50", "0.50-0.75",
                                   "0.75-0.90", "0.90-0.95", "0.95-0.98",
                                   "0.98-1.00"),
                          colors=c("#2b83ba", "#80bfab", "#c7e8ad", "#ffffbf",
                                   "#fdc980", "#f07c4a", "#d7191c"))

species_v4 <- rast(here("zonation", "species_CAZ_proj", "species_CAZ", "species_CAZ_out","species_CAZ.CAZ_E.rank.compressed.tif"))
species_v5 <- rast(here("zonation", "species_scenarios", "species_equal", "output", "rankmap.tif"))
# plot(species_v4, col = mycols(100))
# plot(species_v5, col = mycols(100))

rank_plot <- function(rankmap){
  
  ras_st <- stars::st_as_stars(rankmap)
  ras_sf <- st_as_sf(ras_st)
  fill_ras_sf <- st_drop_geometry(ras_sf)
  ggplot()+
    geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = z_colors_spectral$colors,
                         values = z_colors_spectral$values,
                         name = "Site\nsensitivity",
                         breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
    theme_bw() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank())
  
}


ras_st <- stars::st_as_stars(species_v4)
ras_sf <- st_as_sf(ras_st)
fill_ras_sf <- st_drop_geometry(ras_sf)
ggplot()+
  geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = z_colors_spectral$colors,
                       values = z_colors_spectral$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ras_st <- stars::st_as_stars(species_v5)
ras_sf <- st_as_sf(ras_st)
fill_ras_sf <- st_drop_geometry(ras_sf)
ggplot()+
  geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = z_colors_spectral$colors,
                       values = z_colors_spectral$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

layerCor(c(species_v4, species_v5), "cor", asSample = F)

species_area_v4 <- rast(here("zonation", "species_area_CAZ_proj", "species_area_CAZ", "species_area_CAZ_out","species_area_CAZ.CAZ_E.rank.compressed.tif"))
species_area_v5 <- rast(here("zonation", "species_area_scenarios", "species_area_equal", "output", "rankmap.tif"))

ras_st <- stars::st_as_stars(species_area_v4)
ras_sf <- st_as_sf(ras_st)
fill_ras_sf <- st_drop_geometry(ras_sf)
ggplot()+
  geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = z_colors_spectral$colors,
                       values = z_colors_spectral$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ras_st <- stars::st_as_stars(species_area_v5)
ras_sf <- st_as_sf(ras_st)
fill_ras_sf <- st_drop_geometry(ras_sf)
ggplot()+
  geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = z_colors_spectral$colors,
                       values = z_colors_spectral$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

layerCor(c(species_area_v4, species_area_v5), "cor", asSample = F)
################################## Weighted ####################################
species_w_v4 <- rast(here("zonation", "species_wgt_CAZ_proj", "species_wgt_CAZ", 
                          "species_wgt_CAZ_out",
                          "species_wgt_CAZ.CAZ_E.rank.compressed.tif"))
species_w_v5 <- rast(here("zonation", "species_scenarios", "species_weight", 
                          "output", "rankmap.tif"))

ras_st <- stars::st_as_stars(species_w_v4)
ras_sf <- st_as_sf(ras_st)
fill_ras_sf <- st_drop_geometry(ras_sf)
ggplot()+
  geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = z_colors_spectral$colors,
                       values = z_colors_spectral$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ras_st <- stars::st_as_stars(species_w_v5)
ras_sf <- st_as_sf(ras_st)
fill_ras_sf <- st_drop_geometry(ras_sf)
ggplot()+
  geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = z_colors_spectral$colors,
                       values = z_colors_spectral$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

layerCor(c(species_w_v4, species_w_v5), "cor", asSample = F)



species_area_w_v5 <- rast(here("zonation", "species_area_scenarios", "species_area_weight", 
                          "output", "rankmap.tif"))
ras_st <- stars::st_as_stars(species_area_w_v5)
ras_sf <- st_as_sf(ras_st)
fill_ras_sf <- st_drop_geometry(ras_sf)
ggplot()+
  geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = z_colors_spectral$colors,
                       values = z_colors_spectral$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())



rank_diff <- function(rankmap1, rankmap2){
  
  species_only_diff <- rankmap1 - rankmap2
  coolwarm_hcl <- colorspace::diverging_hcl(11,h = c(250, 10), c = 100, 
                                            l = c(37, 88), power = c(0.7, 1.7))
  species_only_diff_sf <- species_only_diff %>%
    stars::st_as_stars() %>%
    st_as_sf()
  
  ggplot()+
    geom_sf(data = species_only_diff_sf, aes(fill=rankmap), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = rev(coolwarm_hcl),
                         #values = leg$values,
                         name = "Sensitivity\ndifference",
                         breaks = c(-0.5, 0.0, 0.5)) +
    theme_bw() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank())
  
}

r <- rast(here("Eucalyptus_marginata_pred.tif"))
y <- subst(r, NA, 0)
yy <- mask(resample(y, species_area_v5), species_area_v5)