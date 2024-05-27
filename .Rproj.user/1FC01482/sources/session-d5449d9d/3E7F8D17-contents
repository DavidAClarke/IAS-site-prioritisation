#Biodiversity feature performance

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