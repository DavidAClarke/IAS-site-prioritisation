### Testing script

pkgs <- c("tidyverse", "sf", "terra", "stars", "here", "SSIMmap", "raster",
          "ComplexHeatmap", "flextable", "spatstat")
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
  #        path = here(dirname(here()), "figures"))
  
}

#rm(species_fig_list)

#species_area_fig_list <- list()
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
  #        path = here(dirname(here()), "figures")) 
  
}

full_rank_stack <- c(species_rank_stack, species_area_rank_stack)
rm(species_rank_stack, species_area_rank_stack, r)
writeRaster(full_rank_stack, here(dirname(here()), "data", "zonation", "full_rank_stack.tif"))
full_rank_stack <- rast(here(dirname(here()), "data", "zonation", "full_rank_stack.tif"))

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

## SSI heatmaps (average index values)
ssim_mat <- read.csv(here(dirname(here()), "data", "ssim", "ssim_mat.csv"), row.names = 1)
ssim_heat(ssim_mat, pal = "RdPu", nm = "SSIM", lines = F) 

sim_mat <- read.csv(here(dirname(here()), "data", "ssim", "sim_mat.csv"), row.names = 1)
ssim_heat(sim_mat, pal = "GnBu", nm = "SIM", lines = F) 

siv_mat <- read.csv(here(dirname(here()), "data", "ssim", "siv_mat.csv"), row.names = 1)
ssim_heat(siv_mat, pal = "Reds", nm = "SIV", lines = F) 

sip_mat <- read.csv(here(dirname(here()), "data", "ssim", "sip_mat.csv"), row.names = 1)
ssim_heat(sip_mat, pal = "YlOrBr", nm = "SIP", lines = F) 



############################### Alternative ####################################
# I could have sim on one tri and jaccard on the other tri
#ssim_mat[upper.tri(ssim_mat)] <- t(ssim_mat)[upper.tri(ssim_mat)] #or make upper tri jaccard

sip_mmat <- as.matrix(sip_mat)
colnames(sip_mmat) <- gsub("_", " ", colnames(sip_mmat))
jac2.5 <- read.csv(file = here(dirname(here()), "data", "jaccard", "jaccard_two.csv"), row.names = 1)
jac2.5[jac2.5 == "-"] <- NA
jac2.5 <- apply(jac2.5, 2, as.numeric)
rownames(jac2.5) <- colnames(jac2.5)
rownames(jac2.5) <- gsub("_", " ", rownames(jac2.5))
  
col1 <- circlize::colorRamp2(c(0, 0.5, 0.99,1), c("#ef476f", "#ffd166", "#26547c","black"))
col2 <- circlize::colorRamp2(c(0, 0.5, 0.99,1), c("#C04000", "white","#008080","black"))

ht1 <- Heatmap(sip_mmat, 
               name = "SIP",
               rect_gp = gpar(type = "none"), 
               col = col1,
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_row_names = T,
               column_names_gp = gpar(fontsize = 10),
               cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i >= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    grid.text(sprintf("%.2f", sip_mmat[i, j]), x, y, gp = gpar(fontsize = 10))
                  }
                })

ht2 <- Heatmap(jac2.5, 
               name = "Jaccard",
               rect_gp = gpar(type = "none"), 
               col = col2,
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_column_names = F,
               show_row_names = T,
               row_names_side = "right",
               row_names_gp = gpar(fontsize = 10),
               cell_fun = function(j, i, x, y, w, h, fill) {
                if(i <= j) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                  grid.text(sprintf("%.2f", jac2.5[i, j]), x, y, gp = gpar(fontsize = 10))
                }
              })

draw(ht1 + ht2, ht_gap = unit(-247.5, "mm"))
