### Testing script

pkgs <- c("tidyverse", "sf", "terra", "stars", "here", "SSIMmap")
lapply(pkgs, require, character.only = T)

source("R/01_new_functions.R")

#mycols <- grDevices::colorRampPalette(c("#EE602D", "#FFCF26", "#01AEE6"))

z_colors_spectral <- list(values=c(0.0, 0.2, 0.5, 0.75, 0.9, 0.95, 0.98, 1.0),
                          labels=c("0.00-0.20", "0.20-0.50", "0.50-0.75",
                                   "0.75-0.90", "0.90-0.95", "0.95-0.98",
                                   "0.98-1.00"),
                          colors=c("#2b83ba", "#80bfab", "#c7e8ad", "#ffffbf",
                                   "#fdc980", "#f07c4a", "#d7191c"))

species_path <- here(dirname(here()), "data", "zonation", "species_scenarios")
species_area_path <- here(dirname(here()), "data", "zonation", "species_area_scenarios")

species_scenarios <- c("species_equal", "species_weight", "species_scheme_1", 
                       "species_scheme_2", "species_scheme_3", "species_scheme_4",
                       "species_equal_KBA", "species_weight_KBA", "species_scheme_1_KBA", 
                       "species_scheme_2_KBA", "species_scheme_3_KBA", "species_scheme_4_KBA",
                       "species_random")

species_area_scenarios <- c("species_area_equal", "species_area_weight", "species_area_scheme_1", 
                       "species_area_scheme_2", 
                       "species_area_scheme_3", "species_area_scheme_4",
                       "species_area_equal_KBA", "species_area_weight_KBA", "species_area_scheme_1_KBA", 
                       "species_area_scheme_2_KBA", "species_area_scheme_3_KBA", "species_area_scheme_4_KBA",
                       "species_area_random")

# species_fig_list <- list()
species_rank_stack <- rast()

for(sp in species_scenarios){
  
  ind <- which(species_scenarios == sp)
  r <- rast(here(species_path, sp, "output", "rankmap.tif"))
  names(r) <- sp
  species_rank_stack <- c(species_rank_stack, r)
  # species_fig_list[[ind]] <- rank_plot(r)
  # ggsave(plot = species_fig_list[[ind]],
  #        filename = paste0(sp, "_rankmap.pdf"),
  #        device = cairo_pdf,
  #        dpi = 300,
  #        path = here(dirname(here()), "figures")) #dirname() lets you go one folder up
  
}

#rm(species_fig_list)

# species_area_fig_list <- list()
species_area_rank_stack <- rast()

for(sp in species_area_scenarios){
  
  ind <- which(species_area_scenarios == sp)
  r <- rast(here(species_area_path, sp, "output", "rankmap.tif"))
  names(r) <- sp
  species_area_rank_stack <- c(species_area_rank_stack, r)
  # species_area_fig_list[[ind]] <- rank_plot(r)
  # ggsave(plot = species_area_fig_list[[ind]],
  #        filename = paste0(sp, "_rankmap.pdf"),
  #        device = cairo_pdf,
  #        dpi = 300,
  #        path = here(dirname(here()), "figures")) #dirname() lets you go one folder up
  
}

full_rank_stack <- c(species_rank_stack, species_area_rank_stack)

priority_cors <- ras_cor(full_rank_stack)
colnames(priority_cors) <- names(full_rank_stack)
rownames(priority_cors) <- names(full_rank_stack)
write.csv(priority_cors, file = here(dirname(here()), "data", "priority_cors.csv"))

## Structural similarity among sensitive sites
# High similarity between KBA and non-KBA equivalents because only highest sensitive fraction differs
ssims <- ssim(full_rank_stack)
for(i in 1:length(ssims)){
  
  write.csv(ssims[[i]], file = here(dirname(here()), "data", "ssim", paste0(names(ssims)[i], ".csv")))
  
}

ssim_mat <- read.csv(here(dirname(here()), "data", "ssim", "ssim_mat.csv"), row.names = 1)
ssim_mat <- as.matrix(ssim_mat)
ssim_mat[upper.tri(ssim_mat)] <- t(ssim_mat)[upper.tri(ssim_mat)] #or make upper tri jaccard

ssim_df <- as.data.frame(as.table(ssim_mat)) %>%
  drop_na(Freq)
g1 <- ggplot(ssim_df, aes(Var1, Var2, fill= Freq)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdPu")

g2 <- ggplot(ssim_df, aes(Var1, Var2, fill= Freq)) + 
  geom_tile()
  coord_flip() 

############################### Alternative ####################################
# I could have sim on one tri and jaccard on the other tri
library(ComplexHeatmap)
  
col1 = colorRamp2(c(0, 1), c("navy", "orange"))
col2 = colorRamp2(c(0, 1), c("purple", "lightblue"))

ht1 <- Heatmap(ssim_mat, 
               rect_gp = gpar(type = "none"), 
               col = col1,
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_row_names = F,
               cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i >= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                  }
                })

ht2 <- Heatmap(ssim_mat, 
               rect_gp = gpar(type = "none"), 
               col = col2,
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_column_names = F,
               cell_fun = function(j, i, x, y, w, h, fill) {
                if(i <= j) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                }
              })

draw(ht1 + ht2, ht_gap = unit(-200, "mm"))
