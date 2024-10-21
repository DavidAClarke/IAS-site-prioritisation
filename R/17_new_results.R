#Results

################################################################################
##Summary information for all species used (n = 5113)
source("R/13_redlist_summary.R")
sum_info
IAS_threats
IAS_species

#All data & all threats
All_species

#Species threatened by IAS
IAS_threatend_species
IAS_threatend_species_stack

#Red List index - IAS impact
RLI_IAS

################################## Zonation results ############################
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

species_area_scenarios <- c("species_area_equal", "species_area_weight", 
                            "species_area_scheme_1", "species_area_scheme_2", 
                            "species_area_scheme_3", "species_area_scheme_4",
                            "species_area_equal_KBA", "species_area_weight_KBA", 
                            "species_area_scheme_1_KBA", "species_area_scheme_2_KBA", 
                            "species_area_scheme_3_KBA", "species_area_scheme_4_KBA",
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

## Jaccard similarities
source("R/13_jaccard_similarities.R")

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

## Create heatmap with SIP on lower tri and jaccard (top 2%) on upper tri
sip_mmat <- as.matrix(sip_mat)
colnames(sip_mmat) <- gsub("_", " ", colnames(sip_mmat))
jac2.5 <- read.csv(file = here(dirname(here()), "data", "jaccard", "jaccard_two.csv"), row.names = 1)
jac2.5[jac2.5 == "-"] <- NA
jac2.5 <- apply(jac2.5, 2, as.numeric)
rownames(jac2.5) <- colnames(jac2.5)
rownames(jac2.5) <- gsub("_", " ", rownames(jac2.5))

col1 <- circlize::colorRamp2(c(0, 0.5, 0.99,1), 
                             c("#ef476f", "#ffd166", "#26547c","black"))
col2 <- circlize::colorRamp2(c(0, 0.5, 0.99,1), 
                             c("#C04000", "white","#008080","black"))

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

draw(ht1 + ht2, ht_gap = unit(-230, "mm"))

###################### Biodiversity feature performance ########################
# ##Rename groups (may need to do two for the variants with and without areas)
# species_only <- c(CAZ_var, CAZ_wgt_var, CAZ_wgt_KBA_inv_var, RAN_var)
# species_area <- c(CAZ_area_var, CAZ_area_wgt_var, CAZ_area_wgt_KBA_inv_var, 
#                   RAN_area_var)
# 
# species_only_groups <- c("1" = "Invertebrate", "2" = "Fish", "3" = "Plant", 
#                          "4" = "Reptile",
#                          "6" = "Mammal", "7" = "Amphibian", "8" = "Fungi", 
#                          "9" = "Bird")
# 
# species_area_groups <- c("1" = "Invertebrate", "2" = "Fish", "3" = "Plant", 
#                          "4" = "Reptile",
#                          "6" = "Mammal", "7" = "Amphibian", "8" = "Fungi", 
#                          "9" = "Bird",
#                          "10" = "Community", "11" = "Ecosystem", 
#                          "12" = "Ramsar", "13" = "Upstream")
# 
# #species + weights
# perf_1 <- performance_plot(CAZ_wgt_var, species_only_groups, "Set2")
# 
# #species + weights + KBA
# perf_2 <- performance_plot(CAZ_wgt_KBA_inv_var, species_only_groups, "Set2")
#   
# #species + area + weights
# perf_3 <- performance_plot(CAZ_area_wgt_var, species_area_groups, "Set3")
#   
# #species + area + weights + KBA
# perf_4 <- performance_plot(CAZ_area_wgt_KBA_inv_var, species_area_groups, "Set3")
# 
# 
# ggarrange(perf_1, perf_2, common.legend = T, ncol = 2, nrow = 1, labels = c("A","B"))
# ggarrange(perf_3, perf_4, common.legend = T, ncol = 2, nrow = 1, labels = c("A","B"))


########################################################################
##Rank rasters





######################################################################
##Plotting rank rasters


# 
# #Comparing effects of feature weighting  - rank differences
# #Red means weighted gave higher priority over unweighted, 
# #blue means it gave lower priority
# p8 <- rank_diff(CAZ_wgt_var_ras,CAZ_var_ras)
# p9 <- rank_diff(CAZ_area_wgt_var_ras,CAZ_area_var_ras)
#   
# p7 <- ggarrange(p1,p4,p2,p5, 
#                 common.legend = T, 
#                 ncol = 2, nrow = 2,
#                 labels = c("A","C","B","D"))
# 
# p10 <- ggarrange(p8,p9, 
#                 common.legend = T, 
#                 ncol = 2, nrow = 2,
#                 labels = c("E","F"))
# 
# ggarrange(p7,p10, 
#           common.legend = F, 
#           ncol = 1, nrow = 2)
################################################################################
#Susceptible sites
regional_model_path <- here(dirname(here()), "data", "IAS_distributions", "IAS_regional")
spp_list <- c("Apis mellifera",  "Monomorium floricola",
              "Monomorium destructor","Linepithema humile", "Vespula vulgaris",
              "Bombus terrestris", "Heteronychus arator",
              "Digitonthophagus gazella", "Pheidole megacephala",
              "Vespula germanica","Tetramorium bicarinatum",
              "Paratrechina longicornis")

## Prepare for priority sites results
susceptible_site_prep <- lapply(spp_list[1:length(spp_list)], function(i){
  
  susc_site_prep(i, full_rank_stack)
  
})

## Overlap among all IAS
ias_stack <- rast()

for(i in 1:length(susceptible_site_prep)){
  
  r <- susceptible_site_prep[[i]][[2]]
  ias_stack <- c(ias_stack, r)
  
}

ias_rich <- sum(ias_stack)
ias_rich_one <- ias_rich
ias_rich_one[ias_rich_one != 1] <- NA

ias_sum_sp <- resample(ias_rich, full_rank_stack[[2]], method = "near")
ias_sum_sp_ar <- resample(ias_rich, full_rank_stack[[15]], method = "near")

ras_sf <- st_as_stars(ias_sum_sp_ar) %>%
  st_as_sf()

fill_ras_sf <- st_drop_geometry(ras_sf)

cl <- colorRampPalette(c("#e69b99","#2c6184"))

ias_map <- ggplot()+
  geom_sf(data = ras_sf, 
          aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) +
  scale_fill_gradientn(colours = cl(12),
                       name = "Alien\nrichness",
                       breaks = seq(0,11)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

non.na <- as.numeric(global(ias_sum_sp, fun = "notNA")) #number of non NA cells
no.cells <- c()
perc.total <- c()

for(i in 0:11){
  
  c <- length(which(values(ias_sum_sp == i)))
  no.cells <- c(no.cells, c)
  p <- (c/non.na)*100
  perc.total <- c(perc.total, p)
  
}

richness <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
totaltable <- data.frame("Richness" = richness, 
                         "Number of cells" = no.cells, 
                         "Percent total" = round(perc.total,2))

pretty_table <- formattable::formattable(totaltable,
                                         align = c("c", "c", "c"),
                                         list(Richness = formattable::color_tile("#e69b99","#2c6184")))

## IAS SDM evaluations
eval_scores <- read.table(here(regional_model_path, 
                               "Table_IAS_regional_accuracy.txt"), 
                          header = T)
pretty_eval_table <- formattable::formattable(eval_scores,
                                              align = c("l", "c", "c", "c", "c"))


####################### Priority sites - species + weights #####################
species_list <- list()
vals_list <- list()

for(i in 1:length(susceptible_site_prep)){
  for(j in 1:nlyr(full_rank_stack)){
    for(k in 1:nlyr(full_rank_stack)){
      if(j == k){
      
        vals <- get_msk_vals(full_rank_stack[[j]], 
                           susceptible_site_prep[[i]][[1]][[k]])
        vals_list[[j]] <- vals
      
    }
   }
  }
  species_list[[i]] <- vals_list
}

names(species_list) <- spp_list #spp_list is from susceptible_site_prep.R


#Get cell values for Kolmogorov-smirnoff tests
#Could use violin plots to make visual comparisons.
#x axis = sensitivities, y = ias insects; 3 violins per sp?
#maybe need two plots; a species one and a species + areas one.
#simply use geom_violin() instead of geom_boxplot()
#could also lump all together instead of separate species to see combined results

#AM_vals <- c(species_list[[1]][[2]], species_list[[1]][[8]])
vals <- c()
nms <- c()
type <- c()
code <- c()

for(i in 1:length(species_list)){
  
  spv <- unlist(species_list[[i]][1:length(species_list[[i]])])
  vals <- c(vals, spv)
  
  spn <- rep(spp_list[i], length(spv))
  nms <- c(nms, spn)
  
  cc <- paste0(substr(unlist(str_split(spp_list[i], " ")), start = 1, stop = 1), collapse = "")
  spc <- rep(cc, length(spv))
  code <- c(code, spc)
  
  for(j in 1:length(species_list[[i]])){
    
    spt <- rep(names(full_rank_stack)[j], length(species_list[[i]][[j]]))
    type <- c(type, spt)
    
  }
  
}

df <- data.frame(nms, code, type, vals)
write.csv(df, file = here(dirname(here()), "data", "priority_site_vals.csv"))

# df_1 <- df %>% filter(nms == "Pheidole megacephala" |
#                         nms == "Vespula germanica" |
#                         nms == "Digitonthophagus gazella" |
#                         nms == "Paratrechina longicornis" |
#                         nms == "Tetramorium bicarinatum" |
#                         nms == "Apis mellifera")
# 
# df_2 <- df %>% filter(nms == "Monomorium floricola" |
#                         nms == "Monomorium destructor" |
#                         nms == "Linepithema humile" |
#                         nms == "Vespula vulgaris" |
#                         nms == "Bombus terrestris" |
#                         nms == "Heteronychus arator")
# 
# val_plot_1 <- ggplot(df_1, aes(x = vals, y = nms, fill = type)) +
#   geom_violin(draw_quantiles = 0.5, 
#             adjust = 0.2, #changes the smoothness; lower is more faithful to the data
#             scale = "width") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text(size = 12), 
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_blank(),
#         strip.text.x = element_text(size = 12),
#         axis.text.y = element_text(face = "italic"),
#         legend.position = "top") +
#   scale_fill_manual(values = c("#66B2FF", "#FFB266"),
#                     labels = c("KBA mask", "species + weight"),
#                     name = "Scenarios") +
#   scale_y_discrete(limits = rev)
# 
# val_plot_2 <- ggplot(df_2, aes(x = vals, y = nms, fill = type)) +
#   geom_violin(draw_quantiles = 0.5, 
#               adjust = 0.2, #changes the smoothness; lower is more faithful to the data
#               scale = "width") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text(size = 12), 
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_blank(),
#         strip.text.x = element_text(size = 12),
#         axis.text.y = element_text(face = "italic"),
#         legend.position = "top") +
#   scale_fill_manual(values = c("#66B2FF", "#FFB266"),
#                     labels = c("KBA mask", "species + weight"),
#                     name = "Scenarios") +
#   scale_y_discrete(limits = rev)
# 
# ggarrange(val_plot_1, val_plot_2, ncol = 1, nrow = 2, common.legend = T)


# #Multi-panel figures
# PM_spec_comb <- Figure("Pheidole megacephala", CAZ_wgt_var_ras, PM_spec_bin)
# VG_spec_comb <- Figure("Vespula germanica", CAZ_wgt_var_ras, PM_spec_bin)
# DG_spec_comb <- Figure("Digitonthophagus gazella", CAZ_wgt_var_ras, DG_spec_bin)
# PL_spec_comb <- Figure("Paratrechina longicornis", CAZ_wgt_var_ras, PL_spec_bin)
# TB_spec_comb <- Figure("Tetramorium bicarinatum", CAZ_wgt_var_ras, TB_spec_bin)
# AM_spec_comb <- Figure("Apis mellifera", CAZ_wgt_var_ras, AM_spec_bin)
# MF_spec_comb <- Figure("Monomorium floricola", CAZ_wgt_var_ras, MF_spec_bin)
# MD_spec_comb <- Figure("Monomorium destructor", CAZ_wgt_var_ras, MD_spec_bin)
# LH_spec_comb <- Figure("Linepithema humile", CAZ_wgt_var_ras, LH_spec_bin)
# VV_spec_comb <- Figure("Vespula vulgaris", CAZ_wgt_var_ras, VV_spec_bin)
# #MR_comb <- Figure("Megachile rotundata", CAZ_wgt_var_ras, MR_bin)
# BT_spec_comb <- Figure("Bombus terrestris", CAZ_wgt_var_ras, BT_spec_bin)
# HA_spec_comb <- Figure("Heteronychus arator", CAZ_wgt_var_ras, HA_spec_bin)


#Difference between species + weight IAS
lapply(c(species_scenarios, species_area_scenarios), FUN = function(i){
  
  df_min <- df %>% filter(type == i)
  
  g <- ggstatsplot::ggbetweenstats(data = df_min, 
                              y = vals, 
                              x = code, 
                              type = "nonparametric",
                              pairwise.display = "non-significant",
                              p.adjust.method = "bonferroni",
                              package = "awtools",
                              palette = "bpalette",
                              xlab = "Insect species",
                              ylab = "Priority site distribution",
                              ggtheme = ggplot2::theme_bw())
  
  ggsave(plot = g,
         filename = paste0(i, "_spec_comp.pdf"),
         device = cairo_pdf,
         dpi = 300,
         width = 11.69,
         height = 8.27,
         units = "in",
         path = here(dirname(here()), "figures"))
  
  
})


#Proportion difference (KBA vs no KBA vs Random) in number of top sensitive sites
#Top two (i.e. >= 0.98 sensitivity)
props <- c(1, 0.98, 0.95, 0.90, 0.75, 0.50, 0.25, 0.00)
types <- unique(df$type)
sp_props <- c()
sp_name <- c()
tyt <- c()

for(ss in spp_list){
    for(tt in types){
      
      st <- df %>% filter(nms == ss) %>% filter(type == tt)
      mp <- multi_props(st$vals, props)
      sp_props <- c(sp_props,mp)
      sp <- rep(ss, length(mp))
      sp_name <- c(sp_name, sp)
      ty <- rep(tt, length(mp))
      tyt <- c(tyt, ty)
    
  }
}

sens_props <- rep(props, 312)
props_df <- data.frame(sp_name, tyt, sp_props, sens_props)
write.csv(props_df, file = here(dirname(here()), "data", "props_df.csv"))

# Total <- as.data.frame(rbind(PM_spec_prop,VG_spec_prop,DG_spec_prop,
#                              TB_spec_prop,PL_spec_prop,AM_spec_prop,
#                              MF_spec_prop,MD_spec_prop,LH_spec_prop,
#                              VV_spec_prop,BT_spec_prop,HA_spec_prop,
#            PM_KBA_prop, VG_KBA_prop, DG_KBA_prop, TB_KBA_prop, PL_KBA_prop,
#            AM_KBA_prop,MF_KBA_prop,MD_KBA_prop,LH_KBA_prop,VV_KBA_prop,
#            BT_KBA_prop,HA_KBA_prop,
#            PM_RAN_prop, VG_RAN_prop, DG_RAN_prop, TB_RAN_prop, PL_RAN_prop,
#            AM_RAN_prop,MF_RAN_prop,MD_RAN_prop,LH_RAN_prop,VV_RAN_prop,
#            BT_RAN_prop,HA_RAN_prop))
# Type <- c(rep("Unmasked", 12), rep("Masked", 12), rep("Random", 12))
# Species <- rep(c("P. megacephala", "V. germanica", "D. gazella",
#                  "T. bicarinatum", "P. longicornis", "A. mellifera",
#                  "M. floricola", "M. destructor", "L. humile",
#                  "V. vulgaris", "B. terrestris", "H. arator"), 3)
# Total <- cbind(Species, Type,Total)
# Total <- Total %>%
#   as_tibble() %>%
#   mutate(Species = factor(Species)) %>%
#   mutate(Type = factor(Type)) %>%
#   pivot_longer(cols = !Species & !Type, 
#                names_to = "SiteSensitivity", 
#                values_to = "DistributionCoverage") %>%
#   mutate(SiteSensitivity = as.double(SiteSensitivity))
# #cols <- c("#b409a7","#0db02f", "#5bb4f2","#040200","#9c5200")
# Total_min <- Total %>%
#   dplyr::filter(SiteSensitivity == 1 | SiteSensitivity == 0.98)
# df <- data.frame(xmin = c(1,0.98,0.95,0.9,0.75,0.5,0.25),
#                  xmax = c(0.98,0.95,0.9,0.75,0.5,0.25,0),
#                  ymin = rep(-Inf,7),
#                  ymax = rep(Inf,7))
# p <- ggplot(Total, 
#        aes(x = SiteSensitivity, 
#            y = DistributionCoverage, 
#            colour = Species,
#            group = interaction(Species,Type))) +
#   geom_line(stat = "identity", aes(linetype = Type), size = 0.8) +
#   scale_linetype_manual(values = c("solid", "dotted", "longdash")) +
#   scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) + 
#   ylab("Distribution Coverage (%)") +
#   xlab("Site Sensitivity") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text(size = 12), 
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = "bold")) +
#   scale_x_continuous(expand = c(0,0),
#                      limits = c(1,0),
#                      trans = "reverse") +
#   scale_y_continuous(expand = c(0,0))
# p2 <-  p + 
#     geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin=ymin, ymax=ymax),
#               fill = rev(leg$colors), alpha = 0.2, inherit.aes = F) +
#     font("legend.text", face = "italic")
# 
# #For inset
# ins <- ggplot(Total_min, 
#          aes(x = SiteSensitivity, 
#              y = DistributionCoverage, 
#              colour = Species,
#              group = interaction(Species,Type))) +
#     geom_line(stat = "identity", aes(linetype = Type), size = 0.8) +
#     scale_linetype_manual(values = c("solid", "dotted", "longdash")) +
#     scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) +
#     ylab("Distribution Coverage (%)") +
#     xlab("Site Sensitivity") +
#     theme_bw() +
#     theme(axis.line = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = alpha(leg$colors[7],0.2)),
#           axis.text = element_blank(), 
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           legend.text = element_blank(),
#           legend.title = element_blank(),
#           legend.position = "none",
#           axis.ticks = element_blank()) +
#     scale_x_continuous(expand = c(0,0),
#                        limits = c(1,0.98),
#                        trans = "reverse") +
#     scale_y_continuous(expand = c(0,0))
# 
# ################## Priority sites - species + area + weights ###################
# 
# #Get cell values for Kolmogorov-smirnoff tests
# PM_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, PM_area_bin)
# PM_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, PM_area_KBA_bin)
# PM_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, PM_area_RAN_bin)
# 
# VG_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, VG_area_bin)
# VG_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, VG_area_KBA_bin)
# VG_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, VG_area_RAN_bin)
# 
# DG_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, DG_area_bin)
# DG_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, DG_area_KBA_bin)
# DG_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, DG_area_RAN_bin)
# 
# TB_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, TB_area_bin)
# TB_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, TB_area_KBA_bin)
# TB_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, TB_area_RAN_bin)
# 
# PL_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, PL_area_bin)
# PL_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, PL_area_KBA_bin)
# PL_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, PL_area_RAN_bin)
# 
# AM_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, AM_area_bin)
# AM_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, AM_area_KBA_bin)
# AM_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, AM_area_RAN_bin)
# 
# MF_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MF_area_bin)
# MF_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MF_area_KBA_bin)
# MF_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, MF_area_RAN_bin)
# 
# MD_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MD_area_bin)
# MD_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MD_area_KBA_bin)
# MD_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, MD_area_RAN_bin)
# 
# LH_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, LH_area_bin)
# LH_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, LH_area_KBA_bin)
# LH_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, LH_area_RAN_bin)
# 
# VV_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, VV_area_bin)
# VV_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, VV_area_KBA_bin)
# VV_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, VV_area_RAN_bin)
# 
# # MR_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MR_bin)
# # MR_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MR_bin_KBA)
# # MR_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, MR_bin_RAN)
# 
# BT_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, BT_area_bin)
# BT_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, BT_area_KBA_bin)
# BT_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, BT_area_RAN_bin)
# 
# HA_area_vals <- get_msk_vals(CAZ_area_wgt_var_ras, HA_area_bin)
# HA_area_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, HA_area_KBA_bin)
# HA_area_RAN_vals <- get_msk_vals(RAN_area_var_ras, HA_area_RAN_bin)
# 
# #Multi-panel figures
# # PM_area_comb <- Figure("Pheidole megacephala", CAZ_area_wgt_var_ras, PM_area_bin)
# # VG_area_comb <- Figure("Vespula germanica", CAZ_area_wgt_var_ras, VG_area_bin)
# # DG_area_comb <- Figure("Digitonthophagus gazella", CAZ_area_wgt_var_ras, DG_area_bin)
# # PL_area_comb <- Figure("Paratrechina longicornis", CAZ_area_wgt_var_ras, PL_area_bin)
# # TB_area_comb <- Figure("Tetramorium bicarinatum", CAZ_area_wgt_var_ras, TB_area_bin)
# # AM_area_comb <- Figure("Apis mellifera", CAZ_area_wgt_var_ras, AM_area_bin)
# # MF_area_comb <- Figure("Monomorium floricola", CAZ_area_wgt_var_ras, MF_area_bin)
# # MD_area_comb <- Figure("Monomorium destructor", CAZ_area_wgt_var_ras, MD_area_bin)
# # LH_area_comb <- Figure("Linepithema humile", CAZ_area_wgt_var_ras, LH_area_bin)
# # VV_area_comb <- Figure("Vespula vulgaris", CAZ_area_wgt_var_ras, VV_area_bin)
# # #MR_comb <- Figure("Megachile rotundata", CAZ_area_wgt_var_ras, MR_bin)
# # BT_area_comb <- Figure("Bombus terrestris", CAZ_area_wgt_var_ras, BT_area_bin)
# # HA_area_comb <- Figure("Heteronychus arator", CAZ_area_wgt_var_ras, HA_area_bin)
# 
# #Kolmogorov-smirnoff tests - KBA vs no KBA
# #a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution
# #presence of ties a result of rounding
# PM_KBA_area <- wilcox.test(PM_area_KBA_vals, PM_area_vals, alternative = "two.sided")
# VG_KBA_area <- wilcox.test(VG_area_KBA_vals, VG_area_vals, alternative = "two.sided")
# DG_KBA_area <- wilcox.test(DG_area_KBA_vals, DG_area_vals, alternative = "two.sided")
# TB_KBA_area <- wilcox.test(TB_area_KBA_vals, TB_area_vals, alternative = "two.sided")
# PL_KBA_area <- wilcox.test(PL_area_KBA_vals, PL_area_vals, alternative = "two.sided")
# AM_KBA_area <- wilcox.test(AM_area_KBA_vals, AM_area_vals, alternative = "two.sided")
# MF_KBA_area <- wilcox.test(MF_area_KBA_vals, MF_area_vals, alternative = "two.sided")
# MD_KBA_area <- wilcox.test(MD_area_KBA_vals, MD_area_vals, alternative = "two.sided")
# LH_KBA_area <- wilcox.test(LH_area_KBA_vals, LH_area_vals, alternative = "two.sided")
# VV_KBA_area <- wilcox.test(VV_area_KBA_vals, VV_area_vals, alternative = "two.sided")
# #ks.test(MR_KBA_vals, MR_vals, alternative = "two.sided")
# BT_KBA_area <- wilcox.test(BT_area_KBA_vals, BT_area_vals, alternative = "two.sided")
# HA_KBA_area <- wilcox.test(HA_area_KBA_vals, HA_area_vals, alternative = "two.sided")
# 
# 
# #Difference between species + area + weight IAS
# area_vals <- c(PM_area_vals,VG_area_vals,DG_area_vals,TB_area_vals,PL_area_vals,
#                AM_area_vals,MF_area_vals,MD_area_vals,LH_area_vals,VV_area_vals, 
#                BT_area_vals,HA_area_vals)
# 
# spec_names <- c(rep("Pm", length(PM_area_vals)),
#                 rep("Vg", length(VG_area_vals)),
#                 rep("Dg", length(DG_area_vals)),
#                 rep("Tb", length(TB_area_vals)),
#                 rep("Pl", length(PL_area_vals)),
#                 rep("Am", length(AM_area_vals)),
#                 rep("Mf", length(MF_area_vals)),
#                 rep("Md", length(MD_area_vals)),
#                 rep("Lh", length(LH_area_vals)),
#                 rep("Vv", length(VV_area_vals)),
#                 rep("Bt", length(BT_area_vals)),
#                 rep("Ha", length(HA_area_vals)))
# 
# vals <- data.frame(spec_names, area_vals)
# kruskal.test(area_vals ~ spec_names, data = vals)
# FSA::dunnTest(area_vals ~ spec_names, data = vals, method = "bonferroni")
# 
# ggstatsplot::ggbetweenstats(data = vals, 
#                             y = area_vals, 
#                             x = spec_names, 
#                             type = "nonparametric",
#                             pairwise.display = "non-significant",
#                             p.adjust.method = "bonferroni",
#                             package = "awtools",
#                             palette = "bpalette",
#                             xlab = "Insect species",
#                             ylab = "Priority site distribution",
#                             ggtheme = ggplot2::theme_bw())

#Difference between species + area + weight IAS
# PM_VG_area <- ks.test(PM_area_vals, VG_area_vals, alternative = "two.sided")
# PM_DG_area <- ks.test(PM_area_vals, DG_area_vals, alternative = "two.sided")
# PM_TB_area <- ks.test(PM_area_vals, TB_area_vals, alternative = "two.sided")
# PM_PL_area <- ks.test(PM_area_vals, PL_area_vals, alternative = "two.sided")
# PM_AM_area <- ks.test(PM_area_vals, AM_area_vals, alternative = "two.sided")
# PM_MF_area <- ks.test(PM_area_vals, MF_area_vals, alternative = "two.sided")
# PM_MD_area <- ks.test(PM_area_vals, MD_area_vals, alternative = "two.sided")
# PM_LH_area <- ks.test(PM_area_vals, LH_area_vals, alternative = "two.sided")
# PM_VV_area <- ks.test(PM_area_vals, VV_area_vals, alternative = "two.sided")
# #PM_MR <- ks.test(PM_vals, MR_vals, alternative = "two.sided")
# PM_BT_area <- ks.test(PM_area_vals, BT_area_vals, alternative = "two.sided")
# PM_HA_area <- ks.test(PM_area_vals, HA_area_vals, alternative = "two.sided")
# 
# VG_DG_area <- ks.test(VG_area_vals, DG_area_vals, alternative = "two.sided")
# VG_TB_area <- ks.test(VG_area_vals, TB_area_vals, alternative = "two.sided")
# VG_PL_area <- ks.test(VG_area_vals, PL_area_vals, alternative = "two.sided")
# VG_AM_area <- ks.test(VG_area_vals, AM_area_vals, alternative = "two.sided")
# VG_MF_area <- ks.test(VG_area_vals, MF_area_vals, alternative = "two.sided")
# VG_MD_area <- ks.test(VG_area_vals, MD_area_vals, alternative = "two.sided")
# VG_LH_area <- ks.test(VG_area_vals, LH_area_vals, alternative = "two.sided")
# VG_VV_area <- ks.test(VG_area_vals, VV_area_vals, alternative = "two.sided")
# #VG_MR <- ks.test(VG_vals, MR_vals, alternative = "two.sided")
# VG_BT_area <- ks.test(VG_area_vals, BT_area_vals, alternative = "two.sided")
# VG_HA_area <- ks.test(VG_area_vals, HA_area_vals, alternative = "two.sided")
# 
# DG_TB_area <- ks.test(DG_area_vals, TB_area_vals, alternative = "two.sided")
# DG_PL_area <- ks.test(DG_area_vals, PL_area_vals, alternative = "two.sided")
# DG_AM_area <- ks.test(DG_area_vals, AM_area_vals, alternative = "two.sided")
# DG_MF_area <- ks.test(DG_area_vals, MF_area_vals, alternative = "two.sided")
# DG_MD_area <- ks.test(DG_area_vals, MD_area_vals, alternative = "two.sided")
# DG_LH_area <- ks.test(DG_area_vals, LH_area_vals, alternative = "two.sided")
# DG_VV_area <- ks.test(DG_area_vals, VV_area_vals, alternative = "two.sided")
# #DG_MR <- ks.test(DG_vals, MR_vals, alternative = "two.sided")
# DG_BT_area <- ks.test(DG_area_vals, BT_area_vals, alternative = "two.sided")
# DG_HA_area <- ks.test(DG_area_vals, HA_area_vals, alternative = "two.sided")
# 
# TB_PL_area <- ks.test(TB_area_vals, PL_area_vals, alternative = "two.sided")
# TB_AM_area <- ks.test(TB_area_vals, AM_area_vals, alternative = "two.sided")
# TB_MF_area <- ks.test(TB_area_vals, MF_area_vals, alternative = "two.sided")
# TB_MD_area <- ks.test(TB_area_vals, MD_area_vals, alternative = "two.sided")
# TB_LH_area <- ks.test(TB_area_vals, LH_area_vals, alternative = "two.sided")
# TB_VV_area <- ks.test(TB_area_vals, VV_area_vals, alternative = "two.sided")
# #TB_MR <- ks.test(TB_vals, MR_vals, alternative = "two.sided")
# TB_BT_area <- ks.test(TB_area_vals, BT_area_vals, alternative = "two.sided")
# TB_HA_area <- ks.test(TB_area_vals, HA_area_vals, alternative = "two.sided")
# 
# PL_AM_area <- ks.test(PL_area_vals, AM_area_vals, alternative = "two.sided")
# PL_MF_area <- ks.test(PL_area_vals, MF_area_vals, alternative = "two.sided")
# PL_MD_area <- ks.test(PL_area_vals, MD_area_vals, alternative = "two.sided")
# PL_LH_area <- ks.test(PL_area_vals, LH_area_vals, alternative = "two.sided")
# PL_VV_area <- ks.test(PL_area_vals, VV_area_vals, alternative = "two.sided")
# #PL_MR <- ks.test(PL_vals, MR_vals, alternative = "two.sided")
# PL_BT_area <- ks.test(PL_area_vals, BT_area_vals, alternative = "two.sided")
# PL_HA_area <- ks.test(PL_area_vals, HA_area_vals, alternative = "two.sided")
# 
# AM_MF_area <- ks.test(AM_area_vals, MF_area_vals, alternative = "two.sided")
# AM_MD_area <- ks.test(AM_area_vals, MD_area_vals, alternative = "two.sided")
# AM_LH_area <- ks.test(AM_area_vals, LH_area_vals, alternative = "two.sided")
# AM_VV_area <- ks.test(AM_area_vals, VV_area_vals, alternative = "two.sided")
# #AM_MR <- ks.test(AM_vals, MR_vals, alternative = "two.sided")
# AM_BT_area <- ks.test(AM_area_vals, BT_area_vals, alternative = "two.sided")
# AM_HA_area <- ks.test(AM_area_vals, HA_area_vals, alternative = "two.sided")
# 
# MF_MD_area <- ks.test(MF_area_vals, MD_area_vals, alternative = "two.sided")
# MF_LH_area <- ks.test(MF_area_vals, LH_area_vals, alternative = "two.sided")
# MF_VV_area <- ks.test(MF_area_vals, VV_area_vals, alternative = "two.sided")
# #MF_MR <- ks.test(MF_vals, MR_vals, alternative = "two.sided")
# MF_BT_area <- ks.test(MF_area_vals, BT_area_vals, alternative = "two.sided")
# MF_HA_area <- ks.test(MF_area_vals, HA_area_vals, alternative = "two.sided")
# 
# MD_LH_area <- ks.test(MD_area_vals, LH_area_vals, alternative = "two.sided")
# MD_VV_area <- ks.test(MD_area_vals, VV_area_vals, alternative = "two.sided")
# #MD_MR <- ks.test(MD_vals, MR_vals, alternative = "two.sided")
# MD_BT_area <- ks.test(MD_area_vals, BT_area_vals, alternative = "two.sided")
# MD_HA_area <- ks.test(MD_area_vals, HA_area_vals, alternative = "two.sided")
# 
# LH_VV_area <- ks.test(LH_area_vals, VV_area_vals, alternative = "two.sided")
# #LH_MR <- ks.test(LH_vals, MR_vals, alternative = "two.sided")
# LH_BT_area <- ks.test(LH_area_vals, BT_area_vals, alternative = "two.sided")
# LH_HA_area <- ks.test(LH_area_vals, HA_area_vals, alternative = "two.sided")
# 
# #VV_MR <- ks.test(VV_vals, MR_vals, alternative = "two.sided")
# VV_BT_area <- ks.test(VV_area_vals, BT_area_vals, alternative = "two.sided")
# VV_HA_area <- ks.test(VV_area_vals, HA_area_vals, alternative = "two.sided")
# 
# #MR_BT <- ks.test(MR_vals, BT_vals, alternative = "two.sided")
# #MR_HA <- ks.test(MR_vals, HA_vals, alternative = "two.sided")
# 
# BT_HA_area <- ks.test(BT_area_vals, HA_area_vals, alternative = "two.sided")

# v <- c(0,PM_VG_area$statistic,PM_DG_area$statistic,PM_TB_area$statistic,PM_PL_area$statistic,PM_AM_area$statistic,PM_MF_area$statistic,PM_MD_area$statistic,PM_LH_area$statistic,PM_VV_area$statistic,PM_BT_area$statistic,PM_HA_area$statistic,
#        0,0,VG_DG_area$statistic,VG_TB_area$statistic,VG_PL_area$statistic,VG_AM_area$statistic,VG_MF_area$statistic,VG_MD_area$statistic,VG_LH_area$statistic,VG_VV_area$statistic,VG_BT_area$statistic,VG_HA_area$statistic,
#        0,0,0,DG_TB_area$statistic,DG_PL_area$statistic,DG_AM_area$statistic,DG_MF_area$statistic,DG_MD_area$statistic,DG_LH_area$statistic,DG_VV_area$statistic,DG_BT_area$statistic,DG_HA_area$statistic,
#        0,0,0,0,TB_PL_area$statistic,TB_AM_area$statistic,TB_MF_area$statistic,TB_MD_area$statistic,TB_LH_area$statistic,TB_VV_area$statistic,TB_BT_area$statistic,TB_HA_area$statistic,
#        0,0,0,0,0,PL_AM_area$statistic,PL_MF_area$statistic,PL_MD_area$statistic,PL_LH_area$statistic,PL_VV_area$statistic,PL_BT_area$statistic,PL_HA_area$statistic,
#        0,0,0,0,0,0,AM_MF_area$statistic,AM_MD_area$statistic,AM_LH_area$statistic,AM_VV_area$statistic,AM_BT_area$statistic,AM_HA_area$statistic,
#        0,0,0,0,0,0,0,MF_MD_area$statistic,MF_LH_area$statistic,MF_VV_area$statistic,MF_BT_area$statistic,MF_HA_area$statistic,
#        0,0,0,0,0,0,0,0,MD_LH_area$statistic,MD_VV_area$statistic,MD_BT_area$statistic,MD_HA_area$statistic,
#        0,0,0,0,0,0,0,0,0,LH_VV_area$statistic,LH_BT_area$statistic,LH_HA_area$statistic,
#        0,0,0,0,0,0,0,0,0,0,VV_BT_area$statistic,VV_HA_area$statistic,
#        0,0,0,0,0,0,0,0,0,0,0,BT_HA_area$statistic,
#        0,0,0,0,0,0,0,0,0,0,0,0)
# 
# tm <- matrix(v, nrow = 12, ncol = 12)
# rownames(tm) <- c(":italic(Pheidole~~megacephala)", 
#                   ":italic(Vespula~~germanica)", 
#                   ":italic(Digitonthophagus~~gazella)", 
#                   ":italic(Tetramorium~~bicarinatum)", 
#                   ":italic(Paratrechina~~longicornis)",
#                   ":italic(Apis~~mellifera)",
#                   ":italic(Monomorium~~floricola)",
#                   ":italic(Monomorium~~destructor)",
#                   ":italic(Linepithema~~humile)",
#                   ":italic(Vespula~~vulgaris)",
#                   ":italic(Bombus~~terrestris)",
#                   ":italic(Heteronychus~~arator)")
# colnames(tm) <- c("Pm", "Vg", "Dg", "Tb", "Pl", "Am", "Mf", "Md", "Lh", "Vv", "Bt", "Ha")
# corrplot::corrplot(tm, type = "lower", method = "color", 
#                    cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
#                    tl.srt = 0, tl.col = "black", tl.cex = 0.8, addCoef.col = "black",
#                    mar = c(0,0,0,0))

# #Alternative
# rv <- rev(v)
# tmr <- matrix(rv, nrow = 12, ncol = 12)
# rownames(tmr) <- rev(c(":italic(Pheidole~~megacephala)", 
#                   ":italic(Vespula~~germanica)", 
#                   ":italic(Digitonthophagus~~gazella)", 
#                   ":italic(Tetramorium~~bicarinatum)", 
#                   ":italic(Paratrechina~~longicornis)",
#                   ":italic(Apis~~mellifera)",
#                   ":italic(Monomorium~~floricola)",
#                   ":italic(Monomorium~~destructor)",
#                   ":italic(Linepithema~~humile)",
#                   ":italic(Vespula~~vulgaris)",
#                   ":italic(Bombus~~terrestris)",
#                   ":italic(Heteronychus~~arator)"))
# colnames(tm) <- rev(c("Pm", "Vg", "Dg", "Tb", "Pl", "Am", "Mf", "Md", "Lh", "Vv", "Bt", "Ha"))
# corrplot::corrplot(tm, type = "upper", method = "color", 
#                    cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
#                    tl.srt = 0, tl.col = "black", tl.cex = 0.8, addCoef.col = "black",
#                    mar = c(0,0,0,0))

#Proportion difference (KBA vs no KBA) in number of top sensitive sites
#Top two (i.e. >= 0.98 sensitivity)
# props <- c(1, 0.98, 0.95, 0.90, 0.75, 0.50, 0.25, 0.00)
# diffs <- c()
# 
# PM_area_prop <- multi_props(PM_area_vals, props)
# PM_area_KBA_prop <- multi_props(PM_KBA_vals, props)
# PM_area_RAN_prop <- multi_props(PM_area_RAN_vals, props)
# 
# VG_area_prop <- multi_props(VG_area_vals, props)
# VG_area_KBA_prop <- multi_props(VG_area_KBA_vals, props)
# VG_area_RAN_prop <- multi_props(VG_area_RAN_vals, props)
# 
# DG_area_prop <- multi_props(DG_area_vals, props)
# DG_area_KBA_prop <- multi_props(DG_area_KBA_vals, props)
# DG_area_RAN_prop <- multi_props(DG_area_RAN_vals, props)
# 
# TB_area_prop <- multi_props(TB_area_vals, props)
# TB_area_KBA_prop <- multi_props(TB_area_KBA_vals, props)
# TB_area_RAN_prop <- multi_props(TB_area_RAN_vals, props)
# 
# PL_area_prop <- multi_props(PL_area_vals, props)
# PL_area_KBA_prop <- multi_props(PL_area_KBA_vals, props)
# PL_area_RAN_prop <- multi_props(PL_area_RAN_vals, props)
# 
# AM_area_prop <- multi_props(AM_area_vals, props)
# AM_area_KBA_prop <- multi_props(AM_area_KBA_vals, props)
# AM_area_RAN_prop <- multi_props(AM_area_RAN_vals, props)
# 
# MF_area_prop <- multi_props(MF_area_vals, props)
# MF_area_KBA_prop <- multi_props(MF_area_KBA_vals, props)
# MF_area_RAN_prop <- multi_props(MF_area_RAN_vals, props)
# 
# MD_area_prop <- multi_props(MD_area_vals, props)
# MD_area_KBA_prop <- multi_props(MD_area_KBA_vals, props)
# MD_area_RAN_prop <- multi_props(MD_area_RAN_vals, props)
# 
# LH_area_prop <- multi_props(LH_area_vals, props)
# LH_area_KBA_prop <- multi_props(LH_area_KBA_vals, props)
# LH_area_RAN_prop <- multi_props(LH_area_RAN_vals, props)
# 
# VV_area_prop <- multi_props(VV_area_vals, props)
# VV_area_KBA_prop <- multi_props(VV_area_KBA_vals, props)
# VV_area_RAN_prop <- multi_props(VV_area_RAN_vals, props)
# 
# # MR_prop <- multi_props(MR_vals, props)
# # MR_KBA_prop <- multi_props(MR_KBA_vals, props)
# # MR_RAN_prop <- multi_props(MR_RAN_vals, props)
# 
# BT_area_prop <- multi_props(BT_area_vals, props)
# BT_area_KBA_prop <- multi_props(BT_area_KBA_vals, props)
# BT_area_RAN_prop <- multi_props(BT_area_RAN_vals, props)
# 
# HA_area_prop <- multi_props(HA_area_vals, props)
# HA_area_KBA_prop <- multi_props(HA_area_KBA_vals, props)
# HA_area_RAN_prop <- multi_props(HA_area_RAN_vals, props)
# 
# Total <- as.data.frame(rbind(PM_area_prop,VG_area_prop,DG_area_prop,
#                              TB_area_prop,PL_area_prop,AM_area_prop,
#                              MF_area_prop,MD_area_prop,LH_area_prop,
#                              VV_area_prop,BT_area_prop,HA_area_prop,
#                              PM_area_KBA_prop, VG_area_KBA_prop, 
#                              DG_area_KBA_prop, TB_area_KBA_prop, 
#                              PL_area_KBA_prop,AM_area_KBA_prop,MF_area_KBA_prop,
#                              MD_area_KBA_prop,LH_area_KBA_prop,VV_area_KBA_prop,
#                              BT_area_KBA_prop,HA_area_KBA_prop,
#                              PM_area_RAN_prop, VG_area_RAN_prop, 
#                              DG_area_RAN_prop, TB_area_RAN_prop, 
#                              PL_area_RAN_prop,AM_area_RAN_prop,MF_area_RAN_prop,
#                              MD_area_RAN_prop,LH_area_RAN_prop,VV_area_RAN_prop,
#                              BT_area_RAN_prop,HA_area_RAN_prop))
# Type <- c(rep("Unmasked", 12), rep("Masked", 12), rep("Random", 12))
# Species <- rep(c("P. megacephala", "V. germanica", "D. gazella",
#                  "T. bicarinatum", "P. longicornis", "A. mellifera",
#                  "M. floricola", "M. destructor", "L. humile",
#                  "V. vulgaris", "B. terrestris", "H. arator"), 3)
# Total <- cbind(Species, Type,Total)
# Total <- Total %>%
#   as_tibble() %>%
#   mutate(Species = factor(Species)) %>%
#   mutate(Type = factor(Type)) %>%
#   pivot_longer(cols = !Species & !Type, 
#                names_to = "SiteSensitivity", 
#                values_to = "DistributionCoverage") %>%
#   mutate(SiteSensitivity = as.double(SiteSensitivity))
# #cols <- c("#b409a7","#0db02f", "#5bb4f2","#040200","#9c5200")
# Total_min <- Total %>%
#   dplyr::filter(SiteSensitivity == 1 | SiteSensitivity == 0.98)
# df <- data.frame(xmin = c(1,0.98,0.95,0.9,0.75,0.5,0.25),
#                  xmax = c(0.98,0.95,0.9,0.75,0.5,0.25,0),
#                  ymin = rep(-Inf,7),
#                  ymax = rep(Inf,7))
# p <- ggplot(Total, 
#             aes(x = SiteSensitivity, 
#                 y = DistributionCoverage, 
#                 colour = Species,
#                 group = interaction(Species,Type))) +
#   geom_line(stat = "identity", aes(linetype = Type), size = 0.8) +
#   scale_linetype_manual(values = c("solid", "dotted", "longdash")) +
#   scale_colour_manual(values = pnw_palette("Starfish", 12, "continuous")) +
#   ylab("Distribution Coverage (%)") +
#   xlab("Site Sensitivity") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text(size = 12), 
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12, face = "bold")) +
#   scale_x_continuous(expand = c(0,0),
#                      limits = c(1,0),
#                      trans = "reverse") +
#   scale_y_continuous(expand = c(0,0))
# p2 <-  p + 
#   geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin=ymin, ymax=ymax),
#             fill = rev(leg$colors), alpha = 0.2, inherit.aes = F) +
#   font("legend.text", face = "italic")
# 
# #For inset
# ins <- ggplot(Total_min, 
#               aes(x = SiteSensitivity, 
#                   y = DistributionCoverage, 
#                   colour = Species,
#                   group = interaction(Species,Type))) +
#   geom_line(stat = "identity", aes(linetype = Type),size = 0.8) +
#   scale_linetype_manual(values = c("solid", "dotted", "longdash")) +
#   #scale_size_manual(values = c(1,1,1)) +
#   scale_colour_manual(values = pnw_palette("Starfish", 12)) +
#   ylab("Distribution Coverage (%)") +
#   xlab("Site Sensitivity") +
#   theme_bw() +
#   theme(axis.line = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = alpha(leg$colors[7],0.2)),
#         axis.text = element_blank(), 
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.text = element_blank(),
#         legend.title = element_blank(),
#         legend.position = "none",
#         axis.ticks = element_blank()) +
#   scale_x_continuous(expand = c(0,0),
#                      limits = c(1,0.98),
#                      trans = "reverse") +
#   scale_y_continuous(expand = c(0,0))

#An accompanying figure could be showing numerical difference between non-kba and kba
#Take differences, e.g. Here negative values means higher distribution coverage for KBAs
# pm_diff <- PM_area_prop - PM_area_KBA_prop
# vg_diff <- VG_area_prop - VG_area_KBA_prop
# 
# props_df <- read.csv(here(dirname(here()), "data", "props_df.csv"))
# 
# props_df_min <- props_df %>% 
#   filter(tyt == "species_equal" | tyt == "species_equal_KBA")
# 
# #Then just follow through like before but using these values
# Total_diff <- data.frame(rbind(pm_diff, vg_diff))
# species_diff <- c("Pheidole megacephala", "Vespula germanica")
# Total_diff <- cbind(species_diff, Total_diff)
# Total_diff <- Total_diff %>%
#   as_tibble() %>%
#   mutate(species = factor(species)) %>%
#   pivot_longer(cols = !species, 
#                names_to = "SiteSensitivity", 
#                values_to = "ScenarioDifference") %>%
#   mutate(SiteSensitivity = as.double(SiteSensitivity))
# 
# ggplot(Total_diff, 
#        aes(x = SiteSensitivity, 
#            y = ScenarioDifference, 
#            colour = species_diff)) +
#   geom_line(stat = "identity") +
#   scale_x_continuous(expand = c(0,0),
#                      limits = c(1,0),
#                      trans = "reverse") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text(size = 12), 
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12, face = "bold")) +
#   geom_hline(yintercept = 0) +
#   geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin=ymin, ymax=ymax),
#             fill = rev(leg$colors), alpha = 0.2, inherit.aes = F) +
#   font("legend.text", face = "italic")


## Distance to coast
#Prior to doing any kind of distance work, need to project everything to GDA94
full_rank_stack_proj <- project(full_rank_stack, y = "epsg:3112", res = 5000)

## For each variant, produce plot and calculate correlation, between top fraction sites and distance to coast
# ant <- full_rank_stack[[1]]
# ant[is.na(ant)] <- 9999
# ant[ant < 9999] <- NA
# d <- distance(ant)
# writeRaster(d, here(dirname(here()), "data", "dist_to_coast.tif"))

dist_coast <- rast(here(dirname(here()), "data", "dist_to_coast.tif"))
dist_coast <- mask(dist_coast, full_rank_stack[[1]])
dist_coast_proj <- project(dist_coast, y = "epsg:3112", res = 5000)

dists <- values(dist_coast_proj)

#extract site priority at each distance
dist_coord <- as.data.frame(crds(dist_coast_proj))


priority_dist <- data.frame(lapply(full_rank_stack_proj, function(i){
  
  priority <- terra::extract(i, dist_coord, ID = F)
  
  priority <- squeeze(priority[[1]])
  
  names(priority) <- names(i)
  
  return(priority)
  
}), dist_coord)

colnames(priority_dist) <- c(names(full_rank_stack_proj), 
                             "longitude", "latitude")

# Looking for spatial pattern of highest sensitive sites
df_new <- spat_priority_dist(priority_dist, 26)

## Get polygon for cluster tests
my_pol <- dist_coast_proj
my_pol[my_pol >= 0] <- 1
my_pol <- as.polygons(my_pol)

## Create empty list and populate with clustering test
res_list <- list()
plot_list <- list()
for(i in 1:(length(df_new) -2)){
  
  res <- clus_fun(df_new, my_pol, i)
  res_list[[i]] <-  res[[1]]
  plot_list[[i]] <-  res[[2]]
 
}
names(res_list) <- names(df_new)[-c(27,28)]
names(plot_list) <- names(df_new)[-c(27,28)]

## Density plots
#plot(density(plot_list[[1]], sigma = 50000))

## Convert cluster results into data frame
ce_res <- lapply(1:length(res_list), FUN = function(i){
  
  cedf <- as.data.frame(t(as.data.frame(unlist(res_list[[i]]))))
  
})

ce_res <- do.call(rbind, ce_res)
rownames(ce_res) <- names(res_list)
ce_res <- ce_res %>% 
  mutate(statistic.R = round(as.numeric(statistic.R), 2)) %>%
  mutate(p.value = round(as.numeric(p.value), 2)) %>%
  mutate(scenario = rownames(ce_res)) %>%
  relocate(scenario)

ft <- flextable(ce_res)
ft <- autofit(ft, add_w = 0, add_h = 0)
save_as_docx("Table S1" = ft, path = here(dirname(here()), "table_s1.docx"))

## Look at correlation between site sensitivity and distance to the coast
site_coast_cors <- apply(priority_dist, 2, FUN = function(i)
  
  cor.test(i, dists[!is.na(dists)], method = "spearman")
  
)

## Convert correlation results into data frame
sc_res <- lapply(1:(length(site_coast_cors)-2), FUN = function(i){
  
  scdf <- as.data.frame(t(as.data.frame(unlist(site_coast_cors[[i]]))))
  
})

sc_res <- do.call(rbind, sc_res)
rownames(sc_res) <- names(res_list)
