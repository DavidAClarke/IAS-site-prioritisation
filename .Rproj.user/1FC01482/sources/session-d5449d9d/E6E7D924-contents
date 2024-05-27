#Results pertaining to sensitive site identification
###################################### Zonation results ############################
##Load in different projects
#Core Area Zonation
CAZ <- load_zproject(file.path("Zonation", "species_CAZ_proj"))
CAZ_wgt <- load_zproject(file.path("Zonation", "species_wgt_CAZ_proj"))
#CAZ_wgt_KBA <- load_zproject(file.path("Zonation", "species_wgt_CAZ_KBA_proj"))
CAZ_wgt_KBA_inv <- load_zproject(file.path("Zonation", "species_wgt_CAZ_KBA_proj_inv"))
CAZ_area <- load_zproject(file.path("Zonation", "species_area_CAZ_proj"))
CAZ_area_wgt <- load_zproject(file.path("Zonation", "species_area_wgt_CAZ_proj"))
#CAZ_area_wgt_KBA <- load_zproject(file.path("Zonation", "species_area_wgt_CAZ_KBA_proj"))
CAZ_area_wgt_KBA_inv <- load_zproject(file.path("Zonation", "species_area_wgt_CAZ_KBA_proj_inv"))

#Random
RAN <- load_zproject(file.path("Zonation", "species_RAN_proj"))
RAN_wgt <- load_zproject(file.path("Zonation", "species_wgt_RAN_proj"))
#RAN_wgt_KBA <- load_zproject(file.path("Zonation", "species_wgt_RAN_KBA_proj"))
RAN_area <- load_zproject(file.path("Zonation", "species_area_RAN_proj"))
RAN_area_wgt <- load_zproject(file.path("Zonation", "species_area_wgt_RAN_proj"))

##Get all the variants
#Core Area Zonation
CAZ_var <- get_variant(CAZ, 1)
CAZ_wgt_var <- get_variant(CAZ_wgt, 1)
#CAZ_wgt_KBA_var <- get_variant(CAZ_wgt_KBA, 1)
CAZ_wgt_KBA_inv_var <- get_variant(CAZ_wgt_KBA_inv, 1)
CAZ_area_var <- get_variant(CAZ_area, 1)
CAZ_area_wgt_var <- get_variant(CAZ_area_wgt, 1)
#CAZ_area_wgt_KBA_var <- get_variant(CAZ_area_wgt_KBA, 1)
CAZ_area_wgt_KBA_inv_var <- get_variant(CAZ_area_wgt_KBA_inv, 1)

#Random
RAN_var <- get_variant(RAN, 1)
RAN_wgt_var <- get_variant(RAN_wgt, 1)
#RAN_wgt_KBA_var <- get_variant(RAN_wgt_KBA, 1)
RAN_area_var <- get_variant(RAN_area, 1)
RAN_area_wgt_var <- get_variant(RAN_area_wgt, 1)

#############################################################################
##Rename groups (may need to do two for the variants with and without areas)
species_only <- c(CAZ_var, CAZ_wgt_var, CAZ_wgt_KBA_inv_var, RAN_var, RAN_wgt_var)
species_area <- c(CAZ_area_var, CAZ_area_wgt_var, CAZ_area_wgt_KBA_inv_var, RAN_area_var, RAN_area_wgt_var)

species_only_groups <- c("1" = "Invertebrate", "2" = "Fish", "3" = "Plant", "4" = "Reptile",
                         "6" = "Mammal", "7" = "Amphibian", "8" = "Fungi", "9" = "Bird")

species_area_groups <- c("1" = "Invertebrate", "2" = "Fish", "3" = "Plant", "4" = "Reptile",
                         "6" = "Mammal", "7" = "Amphibian", "8" = "Fungi", "9" = "Bird",
                         "10" = "Community", "11" = "Ecosystem", "12" = "Ramsar", "13" = "Upstream")


#species + weights
groupnames(CAZ_wgt_var) <- species_only_groups
lost.levels <- seq(0,1, by = 0.05)
results.caz <- results(CAZ_wgt_var)
perf <- performance(results.caz, lost.levels,melted = TRUE, groups = T)
perf <- na.omit(perf)
perf_wmean <- perf %>% dplyr::filter(str_detect(feature, "w.mean"))
perf_wmean$feature <- as.factor(perf_wmean$feature)
perf_min <- perf %>% dplyr::filter(str_detect(feature, "min."))
perf_min$feature <- as.factor(perf_min$feature)
feature.groups <- c("Amphibian", "Bird", "Fish", "Fungi", "Invertebrate", "Mammal", "Plant", "Reptile")
feat.cols <- RColorBrewer::brewer.pal(8, "Dark2")
names(feat.cols) <- levels(perf_wmean$feature)
colScale <- scale_colour_manual(name = "Feature",
                                values = feat.cols,
                                labels = feature.groups)
perf_1 <- ggplot(data = perf, mapping = aes(x = pr.lost, y = perf.levels, colour = feature)) +
  geom_line(data = perf_wmean, size = 1) +
  colScale +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  ylab("Distribution remaining") +
  xlab("Proportion of landscape lost")

#species + weights + KBA
groupnames(CAZ_wgt_KBA_inv_var) <- species_only_groups
lost.levels <- seq(0,1, by = 0.05)
KBA_results.caz <- results(CAZ_wgt_KBA_inv_var)
perf_KBA <- performance(KBA_results.caz, lost.levels,melted = TRUE, groups = T)
perf_KBA <- na.omit(perf_KBA)
perf_KBA_wmean <- perf_KBA %>% dplyr::filter(str_detect(feature, "w.mean"))
perf_KBA_wmean$feature <- as.factor(perf_wmean$feature)
#perf_min <- perf %>% dplyr::filter(str_detect(feature, "min."))
#perf_min$feature <- as.factor(perf_min$feature)
feature.groups <- c("Amphibian", "Bird", "Fish", "Fungi", "Invertebrate", "Mammal", "Plant", "Reptile")
feat.cols <- RColorBrewer::brewer.pal(8, "Dark2")
names(feat.cols) <- levels(perf_wmean$feature)
colScale <- scale_colour_manual(name = "Feature",
                                values = feat.cols,
                                labels = feature.groups)
perf_2 <- ggplot(data = perf_KBA, mapping = aes(x = pr.lost, y = perf.levels, colour = feature)) +
  geom_line(data = perf_KBA_wmean, size = 1) +
  colScale +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  ylab("Distribution remaining") +
  xlab("Proportion of landscape lost")

ggarrange(perf_1, perf_2, common.legend = T, ncol = 2, nrow = 1, labels = c("A","B"))

#Species + area +weights
groupnames(CAZ_area_wgt_var) <- species_area_groups
lost.levels <- seq(0,1, by = 0.05)
results_area.caz <- results(CAZ_area_wgt_var)
perf_area <- performance(results_area.caz, lost.levels,melted = TRUE, groups = T)
perf_area <- na.omit(perf_area)
perf_area_wmean <- perf_area %>% dplyr::filter(str_detect(feature, "w.mean"))
perf_area_wmean$feature <- as.factor(perf_area_wmean$feature)
#perf_min <- perf %>% dplyr::filter(str_detect(feature, "min."))
#perf_min$feature <- as.factor(perf_min$feature)
feature_area.groups <- c("Amphibian", "Bird", "Community", "Ecosystem", "Fish", "Fungi", "Invertebrate", 
                         "Mammal", "Plant", "Ramsar", "Reptile", "Upstream")
feat_area.cols <- RColorBrewer::brewer.pal(12, "Set3")
names(feat_area.cols) <- levels(perf_area_wmean$feature)
colScale <- scale_colour_manual(name = "Feature",
                                values = feat_area.cols,
                                labels = feature_area.groups)
perf_3 <- ggplot(data = perf_area, mapping = aes(x = pr.lost, y = perf.levels, colour = feature)) +
  geom_line(data = perf_area_wmean, size = 1) +
  colScale +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  ylab("Distribution remaining") +
  xlab("Proportion of landscape lost")

#Species + area +weights + KBA
groupnames(CAZ_area_wgt_KBA_inv_var) <- species_area_groups
lost.levels <- seq(0,1, by = 0.05)
results_area_KBA.caz <- results(CAZ_area_wgt_KBA_inv_var)
perf_area_KBA <- performance(results_area_KBA.caz, lost.levels,melted = TRUE, groups = T)
perf_area_KBA <- na.omit(perf_area_KBA)
perf_area_KBA_wmean <- perf_area_KBA %>% dplyr::filter(str_detect(feature, "w.mean"))
perf_area_KBA_wmean$feature <- as.factor(perf_area_KBA_wmean$feature)
#perf_min <- perf %>% dplyr::filter(str_detect(feature, "min."))
#perf_min$feature <- as.factor(perf_min$feature)
feature_area.groups <- c("Amphibian", "Bird", "Community", "Ecosystem", "Fish", "Fungi", "Invertebrate", 
                         "Mammal", "Plant", "Ramsar", "Reptile", "Upstream")
feat_area.cols <- RColorBrewer::brewer.pal(12, "Set3")
names(feat_area.cols) <- levels(perf_area_wmean$feature)
colScale <- scale_colour_manual(name = "Feature",
                                values = feat_area.cols,
                                labels = feature_area.groups)
perf_4 <- ggplot(data = perf_area_KBA, mapping = aes(x = pr.lost, y = perf.levels, colour = feature)) +
  geom_line(data = perf_area_KBA_wmean, size = 1) +
  colScale +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  ylab("Distribution remaining") +
  xlab("Proportion of landscape lost")

ggarrange(perf_3, perf_4, common.legend = T, ncol = 2, nrow = 1, labels = c("A","B"))
########################################################################
##Rank rasters
#species only - no weights
CAZ_var_ras <- rank_raster(CAZ_var)
RAN_var_ras <- rank_raster(RAN_var)

#species only - with weights
CAZ_wgt_var_ras <- rank_raster(CAZ_wgt_var)
RAN_wgt_var_ras <- rank_raster(RAN_wgt_var)

#species only - with weights and KBA mask
#CAZ_wgt_KBA_var_ras <- rank_raster(CAZ_wgt_KBA_var)
RAN_wgt_KBA_var_ras <- rank_raster(RAN_wgt_KBA_var)

#species only - with weights and inverted KBA mask
CAZ_wgt_KBA_inv_var_ras <- rank_raster(CAZ_wgt_KBA_inv_var)

#species and areas - no weights
CAZ_area_var_ras <- rank_raster(CAZ_area_var)
RAN_area_var_ras <- rank_raster(RAN_area_var)

#species and areas - with weights
CAZ_area_wgt_var_ras <- rank_raster(CAZ_area_wgt_var)
RAN_area_wgt_var_ras <- rank_raster(RAN_area_wgt_var)

#species and areas - with weights and KBA mask
#CAZ_area_wgt_KBA_var_ras <- rank_raster(CAZ_area_wgt_KBA_var)

#species and areas - with weights and inverted KBA mask
CAZ_area_wgt_KBA_inv_var_ras <- rank_raster(CAZ_area_wgt_KBA_inv_var)

CAZ_stack <- stack(CAZ_var_ras,
                   CAZ_wgt_var_ras,
                   CAZ_wgt_KBA_inv_var_ras,
                   RAN_wgt_var_ras,
                   CAZ_area_var_ras,
                   CAZ_area_wgt_var_ras,
                   CAZ_area_wgt_KBA_inv_var_ras,
                   RAN_area_wgt_var_ras)
pairs(CAZ_stack, method = "kendall")

# ##Values
# C_vals <- as.vector(CAZ_var_ras)
# Cw_vals <- as.vector(CAZ_wgt_var_ras)
# CwKBA_vals <- as.vector(CAZ_wgt_KBA_var_ras)
# CA_vals <- as.vector(CAZ_area_var_ras)
# CAw_vals <- as.vector(CAZ_area_wgt_var_ras)
# CAwKBA_vals <- as.vector(CAZ_area_wgt_KBA_var_ras)
# 
# #The following are equivalent
# kencor.1 <- cor.test(C_vals, Cw_vals, alternative = "two.sided", method = "kendall", exact = T) #this performs statistical test
# kencor.2 <- cor.test(C_vals, CwKBA_vals, alternative = "two.sided", method = "kendall", exact = T) #this performs statistical test
# kencor.3 <- cor.test(Cw_vals, CwKBA_vals, alternative = "two.sided", method = "kendall", exact = T) #this performs statistical test
# kencor.4 <- cor.test(CA_vals, CAw_vals, alternative = "two.sided", method = "kendall", exact = T) #this performs statistical test
# kencor.5 <- cor.test(CA_vals, CAwKBA_vals, alternative = "two.sided", method = "kendall", exact = T) #this performs statistical test
# kencor.6 <- cor.test(CAw_vals, CAwKBA_vals, alternative = "two.sided", method = "kendall", exact = T) #this performs statistical test
# cor(CAZ_var_ras[1:791,1:972], CAZ_wgt_var_ras[1:791,1:972], method = "kendall", use = "complete.obs") #this just gives coefficient

######################################################################
##Plotting rank rasters
leg <- zlegend("spectral")

p1 <- rank_plot(CAZ_var_ras)
p2 <- rank_plot(CAZ_wgt_var_ras)
p3 <- rank_plot(CAZ_wgt_KBA_var_ras)
p4 <- rank_plot(CAZ_area_var_ras)
p5 <- rank_plot(CAZ_area_wgt_var_ras)
p6 <- rank_plot(CAZ_area_wgt_KBA_var_ras)

p11 <- rank_plot(RAN_var_ras)
p12 <- rank_plot(RAN_wgt_var_ras)
p13 <- rank_plot(RAN_area_var_ras)
p14 <- rank_plot(RAN_area_wgt_var_ras)

#Couple of examples
species_only_diff <- CAZ_wgt_var_ras - CAZ_var_ras
diff_vals <- values(species_only_diff)
hist(diff_vals)

#plot(species_only_diff, col = coolwarm_hcl) #Red means weighted gave higher priority over unweighted, blue means it gave lower priority
#Comparing rank differences
coolwarm_hcl <- colorspace::diverging_hcl(11,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
species_only_diff_st <- st_as_stars(species_only_diff)
species_only_diff_sf <- st_as_sf(species_only_diff_st)
p8 <- ggplot()+
  geom_sf(data = species_only_diff_sf, aes(fill=layer), 
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


species_area_diff <- CAZ_area_wgt_var_ras - CAZ_area_var_ras
diff_vals2 <- values(species_area_diff)
hist(diff_vals2)
#plot(species_area_diff, col = coolwarm_hcl)
species_area_diff_st <- st_as_stars(species_area_diff)
species_area_diff_sf <- st_as_sf(species_area_diff_st)
p9 <- ggplot()+
  geom_sf(data = species_area_diff_sf, aes(fill=layer), 
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

p7 <- ggarrange(p1,p4,p2,p5, 
                common.legend = T, 
                ncol = 2, nrow = 2,
                labels = c("A","C","B","D"))

p10 <- ggarrange(p8,p9, 
                 common.legend = T, 
                 ncol = 2, nrow = 2,
                 labels = c("E","F"))

ggarrange(p7,p10, 
          common.legend = F, 
          ncol = 1, nrow = 2)