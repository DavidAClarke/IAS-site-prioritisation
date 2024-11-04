#Results

################################################################################
##Summary information for all species used (n = 5113)
source("R/15_redlist_summary.R")
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

species_fig_list <- list()
species_rank_stack <- rast()

for(sp in species_scenarios){
  
  ind <- which(species_scenarios == sp)
  r <- rast(here(species_path, sp, "output", "rankmap.tif"))
  names(r) <- sp
  species_rank_stack <- c(species_rank_stack, r)
  species_fig_list[[ind]] <- rank_plot(r) +

    ggtitle(sp) +

    theme(plot.title = element_text(face = "italic", size = 10, hjust = 0.5))
  # ggsave(plot = species_fig_list[[ind]],
  #        filename = paste0(sp, "_rankmap.pdf"),
  #        device = cairo_pdf,
  #        dpi = 300,
  #        path = here(dirname(here()), "figures"))
  
}

ggpubr::ggarrange(plotlist = species_fig_list[-13], nrow = 4, ncol = 3, common.legend = T)

species_area_fig_list <- list()
species_area_rank_stack <- rast()

for(sp in species_area_scenarios){
  
  ind <- which(species_area_scenarios == sp)
  r <- rast(here(species_area_path, sp, "output", "rankmap.tif"))
  names(r) <- sp
  species_area_rank_stack <- c(species_area_rank_stack, r)
  species_area_fig_list[[ind]] <- rank_plot(r) +

    ggtitle(sp) +

    theme(plot.title = element_text(face = "italic", size = 10, hjust = 0.5))
  # ggsave(plot = species_area_fig_list[[ind]],
  #        filename = paste0(sp, "_rankmap.pdf"),
  #        device = cairo_pdf,
  #        dpi = 300,
  #        path = here(dirname(here()), "figures")) 
  
}

ggpubr::ggarrange(plotlist = species_area_fig_list[-13], nrow = 4, ncol = 3, common.legend = T)

full_rank_stack <- c(species_rank_stack, species_area_rank_stack)
rm(species_rank_stack, species_area_rank_stack, r)
writeRaster(full_rank_stack, here(dirname(here()), "data", "zonation", "full_rank_stack.tif"))
full_rank_stack <- rast(here(dirname(here()), "data", "zonation", "full_rank_stack.tif"))

## Correlations
priority_cors <- ras_cor(full_rank_stack)
colnames(priority_cors) <- names(full_rank_stack)
rownames(priority_cors) <- names(full_rank_stack)
write.csv(priority_cors, file = here(dirname(here()), "data", "priority_cors.csv"))

priority_cors <- read.csv(here(dirname(here()), "data", "priority_cors.csv"), row.names = 1)

corrplot::corrplot(round(as.matrix(priority_cors),2), 
                   type = "lower", 
                   method = "color",
                   cl.pos = "n", 
                   tl.pos = "l",
                   col=colorRampPalette(brewer.pal(n=9, name="Blues"))(15),
                   tl.srt = 0, 
                   tl.col = "black", 
                   tl.cex = 0.8, 
                   number.cex = 0.8,
                   addCoef.col = "white",
                   mar = c(0,0,0,0))

## Jaccard similarities
source("R/14_jaccard_similarities.R")

## Structural similarity among sensitive sites
# High similarity between KBA and non-KBA equivalents because only highest sensitive fraction differs
ssims <- ssim(full_rank_stack)

for(i in 1:length(ssims)){
  
  write.csv(ssims[[i]], file = here(dirname(here()), "data", "ssim", 
                                    paste0(names(ssims)[i], ".csv")))
  
}
################################################################################
## SSIM maps
r <- rast(here(dirname(here()), 
            "data", "ssim", "species_area_weight-species_area_weight_KBA.tif"))

r1 <- st_as_stars(r[[2]]) %>%
    st_as_sf()

r1v <- st_drop_geometry(r1)
  
gsim <- ggplot() +
      geom_sf(data = r1, 
            aes(fill=r1v[,1]), 
            color=NA, 
            show.legend = T) +
      colorspace::scale_fill_continuous_sequential("ag_sunset",
                                                    name = "SIM") +
      theme_bw() +
      theme(axis.line = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank())

r2 <- st_as_stars(r[[3]]) %>%
  st_as_sf()

r2v <- st_drop_geometry(r2)

gsiv <- ggplot() +
  geom_sf(data = r2, 
          aes(fill=r2v[,1]), 
          color=NA, 
          show.legend = T) +
  colorspace::scale_fill_continuous_sequential("batlow",
                                               name = "SIV") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

r3 <- st_as_stars(r[[4]]) %>%
  st_as_sf()

r3v <- st_drop_geometry(r3)

gsip <- ggplot() +
  geom_sf(data = r3, 
          aes(fill=r3v[,1]), 
          color=NA, 
          show.legend = T) +
  colorspace::scale_fill_continuous_diverging("Blue-Red",
                                               name = "SIP") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

r4 <- st_as_stars(full_rank_stack[[15]] - full_rank_stack[[21]]) %>%
  st_as_sf()

r4v <- st_drop_geometry(r4)

gdif <- ggplot() +
  geom_sf(data = r4, 
          aes(fill=r4v[,1]), 
          color=NA, 
          show.legend = T) +
  colorspace::scale_fill_continuous_diverging("Purple-Green",
                                              name = "DIFF") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

g1 <- rank_plot(full_rank_stack[[15]])
g2 <- rank_plot(full_rank_stack[[21]])

ggpubr::ggarrange(g1, gdif, gsim, 
                  g2, gsiv, gsip, 
                  nrow = 2, ncol = 3, 
                  common.legend = F,
                  labels = "AUTO")
################################################################################
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


################################ Susceptible sites #############################
regional_model_path <- here(dirname(here()), "data", "IAS_distributions", "IAS_regional")

spp_list <- c("Apis mellifera",  "Monomorium floricola",
              "Monomorium destructor","Linepithema humile", "Vespula vulgaris",
              "Bombus terrestris", "Heteronychus arator",
              "Digitonthophagus gazella", "Pheidole megacephala",
              "Vespula germanica","Tetramorium bicarinatum",
              "Paratrechina longicornis")

## Create IAS distribution plots
ias_list <- list()

for(sp in spp_list){

  ind <- which(spp_list == sp)

  ias_list[[ind]] <- IAS_plot(sp) +

    ggtitle(sp) +

    theme(plot.title = element_text(face = "italic", size = 12, hjust = 0.5))

  spn <- gsub(" ", "_", sp)

  ggsave(plot = ias_list[[ind]],
         filename = paste0(spn, "_dist.pdf"),
         device = cairo_pdf,
         dpi = 300,
         path = here(dirname(here()), "figures"))

}

ggpubr::ggarrange(plotlist = ias_list, nrow = 4, ncol = 3, common.legend = T)

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

ggsave(plot = ias_map,
       filename = "ias_map.pdf",
       device = cairo_pdf,
       dpi = 300,
       path = here(dirname(here()), "figures"))

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

total_ft <- autofit(flextable(totaltable)) %>%
  bg(j = "Richness", bg = colorRampPalette(c("#e69b99","#2c6184"))(12), part = "body")

total_ft2 <- labelizor(x = total_ft,
                       part = "header",
                       labels = c("Number.of.cells" = "No. cells",
                                  "Percent.total" = "Total (%)"))

save_as_docx(total_ft2, path = here(dirname(here()), "total_ft2.docx"))
save_as_image(total_ft2, path = here(dirname(here()), "total_ft2.svg"), res = 500)

## IAS SDM evaluations
eval_scores <- read.table(here(regional_model_path, 
                               "Table_IAS_regional_accuracy.txt"), 
                          header = T)

eval_ft <- autofit(flextable(eval_scores))

eval_ft2 <- labelizor(x = eval_ft,
                       part = "header",
                       labels = c("Cut.off.binary" = "Binary cutoff"))

save_as_docx(eval_ft2, path = here(dirname(here()), "eval_scores.docx"))
save_as_image(eval_ft2, path = here(dirname(here()), "eval_scores.svg"), res = 500)

############################### Priority sites  ################################
## Priority maps
pri_maps <- list()

for(i in seq_along(susceptible_site_prep)){
  
  sp <- spp_list[[i]]
  
  pri_maps[[i]] <- priority_map(full_rank_stack[["species_area_weight"]], 
                                susceptible_site_prep[[i]][[2]]) +
    
    ggtitle(sp) +
    
    theme(plot.title = element_text(face = "italic", size = 12, hjust = 0.5))
  
}

ggpubr::ggarrange(plotlist = pri_maps, nrow = 4, ncol = 3, common.legend = T)


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

ps_df <- read.csv(here(dirname(here()), "data", "priority_site_vals.csv"))

# Want to loop over species and scenarios,comparing KBA and non-KBA scenarios
# for each species. First species, then scenarios
# Separate out KBA ones
scenarios <- c(species_scenarios, species_area_scenarios)
KBA <- scenarios[str_detect(scenarios, "KBA")]
nonKBA <- scenarios[str_detect(scenarios, "KBA", negate = T)]
nonKBA <- nonKBA[str_detect(nonKBA, "random", negate = T)]
scenarios <- c(nonKBA, KBA)

sp <- c()
wcs <- c()
wcp <- c()
ty1 <- c()
ty2 <- c()

for(i in spp_list){
  
  df_min <- ps_df %>% filter(nms == i)
  
  for(j in seq_along(nonKBA)){
    
    df <- df_min %>% filter(type == nonKBA[j] | type == KBA[j])
    
    v1 <- df %>% filter(type == nonKBA[j]) %>% pull(vals)
    v2 <- df %>% filter(type == KBA[j]) %>% pull(vals)
    
    wc <- wilcox.test(v1, v2, alternative = "two.sided")
    wcs <- c(wcs, wc$statistic)
    wcp <- c(wcp, wc$p.value)
    sp <- c(sp, i)
    ty1 <- c(ty1, nonKBA[j])
    ty2 <- c(ty2, KBA[j])
    
    
  }
}

wc_df <- data.frame(sp, ty1, ty2, wcs, wcp)
write.csv(wc_df, here(dirname(here()), "data", "wilcox_df.csv"))

################################################################################
## Calculate median site sensitivity for each species X scenario

spec <- c()
scen <- c()
medsen <- c()

for(i in spp_list){
  
  df_min <- ps_df %>% filter(nms == i) 
  
  for(sc in scenarios){
    
    vs <- df_min %>% filter(type == sc) %>% pull(vals)
    spec <- c(spec, i)
    scen <- c(scen, sc)
    medsen <- c(medsen, median(vs))
    
   }
}

med_df <- data.frame(spec, scen, medsen)
write.csv(med_df, here(dirname(here()), "data", "medsen_df.csv"))

bpalette <- c('#c62828','#f44336','#9c27b0','#673ab7','#3f51b5','#2196f3',
             '#29b6f6','#006064','#009688','#4caf50','#8bc34a','#ffeb3b')


med_df <- med_df %>% mutate(scen = gsub("_", " ", scen))

ggplot(med_df, aes(x = spec, y = medsen, fill = spec)) +
  
  geom_boxplot() +
  
  scale_fill_manual(values = bpalette) +
  
  theme_bw() +
  
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(), 
        legend.position = "none") +
  
  ylab("Median value of site sensitivity") +
  
  xlab("Introduced species")

## Compare scenarios
ggplot(med_df, aes(x = scen, y = medsen, fill = scen)) +
  
  geom_boxplot() +
  
  #scale_fill_manual(values = bpalette) +
  
  theme_bw() +
  
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(), 
        legend.position = "none") +
  
  ylab("Median value of site sensitivity") +
  
  xlab("Sensitive site identification scenario")

################################################################################
for(i in seq_along(nonKBA)){
  
  sp <- gsub("_", " ", nonKBA[i])
  
  scen_one <- ps_df %>% filter(type == nonKBA[i] | type == KBA[i])
  
  g <- ggplot(scen_one, aes(x = vals, y = reorder(nms, desc(nms)), fill = type, height = after_stat(density))) +
    ggridges::geom_density_ridges(stat = "density", alpha = 0.7, scale = 1) +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.text.y = element_text(face = "italic"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          legend.position = "top") +
    scale_x_continuous(expand = c(0.01,0.01)) +
    scale_fill_manual(values = c("#FDDB27FF", "#00B1D2FF"),
                      labels = c(sp, "KBA mask"),
                      name = "Scenarios") +
    ylab("Introduced species") +
    xlab("Site sensitivity")
  
  ggsave(plot = g,
         filename = paste0(nonKBA[i], "_violin.pdf"),
         device = cairo_pdf,
         dpi = 300,
         width = 8.27,
         height = 11.69,
         units = "in",
         path = here(dirname(here()), "figures"))
  
  
}


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
sc_res$scenario <- names(res_list)
sc_res <- sc_res %>% 
  relocate(scenario) %>%
  dplyr::select(scenario, estimate.rho, p.value) %>%
  mutate(across(2:3, as.numeric)) %>%
  mutate(across(2:3, round,2))

ft <- flextable(sc_res)
ft <- autofit(ft, add_w = 0, add_h = 0)
save_as_docx("Table S3" = ft, path = here(dirname(here()), "table_s3.docx"))
