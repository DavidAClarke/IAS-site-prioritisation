#Results
library(zonator)
library(tidyverse)
library(spatstat)
library(sf)
library(stars)
library(raster)
library(ggpubr)
library(cowplot)
#################################################################################
##Summary information for all species used (n = 5113)
species_info <- read.csv(file.path("Zonation", "maxent_spp_list.csv"))
sum_info <- species_info %>%
  mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern",
                                             "Near Threatened", "Vulnerable",
                                             "Endangered", "Critically Endangered"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", "CR"))) %>%
  group_by(classGroup, redlistCategory) %>%
  summarise(species_per_cat = n()) %>%
  mutate(species_per_group = sum(species_per_cat))

IAS_threats <- species_info %>%
  dplyr::filter(code == "8.1.1" | code == "8.1.2") %>%
  mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern",
                                             "Near Threatened", "Vulnerable",
                                             "Endangered", "Critically Endangered"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", "CR"))) %>%
  group_by(classGroup, redlistCategory) %>%
  summarise(species_per_cat = n()) %>%
  mutate(species_per_group = sum(species_per_cat)) %>%
  mutate(prop.total = species_per_cat/species_per_group)

IAS_species <- species_info %>%
  dplyr::select(ias) %>%
  count(ias) %>%
  drop_na(ias)

#Red List plot
mycols <- c("#808080", "#008000", "#ADFF2F", "#FFFF00", "#FFA500", "#FF0000") #IUCN colours

#All data & all threats
All_species <- ggplot(sum_info, aes(x = redlistCategory, y = species_per_cat, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  facet_wrap(~classGroup, nrow = 2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(trans = "log1p") +
  ylab("Number of species (log + 1)") +
  xlab("IUCN Red List category")

#Species threatened by IAS
All_species <- ggplot(IAS_threats, aes(x = redlistCategory, y = species_per_cat, fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  facet_wrap(~classGroup, nrow = 2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(trans = "log1p") + #could also change to proportion of total
  ylab("Number of species (log + 1)") +
  xlab("IUCN Red List category")

#Alternative stacked plot
All_species <- ggplot(IAS_threats, aes(x = classGroup, y = species_per_cat, fill = redlistCategory)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = mycols, name = "Red List\nCategory") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Proportion of species") +
  xlab("Group")

#Red List index - IAS impact
IAS_impact <- read.csv(file.path("ISO_BL_AUS_aggregated_general.csv"))
ggplot(IAS_impact, aes(x = year, y = rli)) + 
  geom_line(colour = "blue") + 
  geom_ribbon(aes(ymin = qn05, ymax = qn95), linetype = 2, alpha = 0.1) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.75, 1)) +
  #scale_x_discrete(breaks = IAS_impact$year) +
  ylab("Red List Index") +
  xlab("Year")

###################################### Zonation results ############################
##Load in different projects
CAZ <- load_zproject(file.path("Zonation", "species_CAZ_proj"))
CAZ_wgt <- load_zproject(file.path("Zonation", "species_wgt_CAZ_proj"))
#CAZ_wgt_KBA <- load_zproject(file.path("Zonation", "species_wgt_CAZ_KBA_proj"))
CAZ_wgt_KBA_inv <- load_zproject(file.path("Zonation", "species_wgt_CAZ_KBA_proj_inv"))
CAZ_area <- load_zproject(file.path("Zonation", "species_area_CAZ_proj"))
CAZ_area_wgt <- load_zproject(file.path("Zonation", "species_area_wgt_CAZ_proj"))
#CAZ_area_wgt_KBA <- load_zproject(file.path("Zonation", "species_area_wgt_CAZ_KBA_proj"))
CAZ_area_wgt_KBA_inv <- load_zproject(file.path("Zonation", "species_area_wgt_CAZ_KBA_proj_inv"))

##Get all the variants
CAZ_var <- get_variant(CAZ, 1)
CAZ_wgt_var <- get_variant(CAZ_wgt, 1)
#CAZ_wgt_KBA_var <- get_variant(CAZ_wgt_KBA, 1)
CAZ_wgt_KBA_inv_var <- get_variant(CAZ_wgt_KBA_inv, 1)
CAZ_area_var <- get_variant(CAZ_area, 1)
CAZ_area_wgt_var <- get_variant(CAZ_area_wgt, 1)
#CAZ_area_wgt_KBA_var <- get_variant(CAZ_area_wgt_KBA, 1)
CAZ_area_wgt_KBA_inv_var <- get_variant(CAZ_area_wgt_KBA_inv, 1)

#############################################################################
##Rename groups (may need to do two for the variants with and without areas)
species_only <- c(CAZ_var, CAZ_wgt_var, CAZ_wgt_KBA_inv_var)
species_area <- c(CAZ_area_var, CAZ_area_wgt_var, CAZ_area_wgt_KBA_inv_var)

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
#species only - with weights
CAZ_wgt_var_ras <- rank_raster(CAZ_wgt_var)
#species only - with weights and KBA mask
#CAZ_wgt_KBA_var_ras <- rank_raster(CAZ_wgt_KBA_var)
#species only - with weights and inverted KBA mask
CAZ_wgt_KBA_inv_var_ras <- rank_raster(CAZ_wgt_KBA_inv_var)
#species and areas - no weights
CAZ_area_var_ras <- rank_raster(CAZ_area_var)
#species and areas - with weights
CAZ_area_wgt_var_ras <- rank_raster(CAZ_area_wgt_var)
#species and areas - with weights and KBA mask
#CAZ_area_wgt_KBA_var_ras <- rank_raster(CAZ_area_wgt_KBA_var)
#species and areas - with weights and inverted KBA mask
CAZ_area_wgt_KBA_inv_var_ras <- rank_raster(CAZ_area_wgt_KBA_inv_var)
CAZ_stack <- stack(CAZ_var_ras,
                   CAZ_wgt_var_ras,
                   CAZ_wgt_KBA_inv_var_ras,
                   CAZ_area_var_ras,
                   CAZ_area_wgt_var_ras,
                   CAZ_area_wgt_KBA_inv_var_ras)
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
#p1 <- plot(CAZ_var_ras, breaks = leg$values, col = leg$colors)
CAZ_var_ras_st <- st_as_stars(CAZ_var_ras)
CAZ_var_ras_sf <- st_as_sf(CAZ_var_ras_st)
p1 <- ggplot()+
  geom_sf(data = CAZ_var_ras_sf, aes(fill=species_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site sensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

#p2 <- plot(CAZ_wgt_var_ras, breaks = leg$values, col = leg$colors)
CAZ_wgt_var_ras_st <- st_as_stars(CAZ_wgt_var_ras)
CAZ_wgt_var_ras_sf <- st_as_sf(CAZ_wgt_var_ras_st)
p2 <- ggplot()+
  geom_sf(data = CAZ_wgt_var_ras_sf, aes(fill=species_wgt_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site sensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

#p3 <- plot(CAZ_wgt_KBA_var_ras, breaks = leg$values, col = leg$colors)
CAZ_wgt_KBA_var_ras_st <- st_as_stars(CAZ_wgt_KBA_var_ras)
CAZ_wgt_KBA_var_ras_sf <- st_as_sf(CAZ_wgt_KBA_var_ras_st)
p3 <- ggplot()+
  geom_sf(data = CAZ_wgt_KBA_var_ras_sf, aes(fill=species_wgt_CAZ_KBA.CAZ_ME.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site sensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

#p4 <- plot(CAZ_area_var_ras, breaks = leg$values, col = leg$colors)
CAZ_area_var_ras_st <- st_as_stars(CAZ_area_var_ras)
CAZ_area_var_ras_sf <- st_as_sf(CAZ_area_var_ras_st)
p4 <- ggplot()+
  geom_sf(data = CAZ_area_var_ras_sf, aes(fill=species_area_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site sensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())
#p5 <- plot(CAZ_area_wgt_var_ras, breaks = leg$values, col = leg$colors)
CAZ_area_wgt_var_ras_st <- st_as_stars(CAZ_area_wgt_var_ras)
CAZ_area_wgt_var_ras_sf <- st_as_sf(CAZ_area_wgt_var_ras_st)
p5 <- ggplot()+
  geom_sf(data = CAZ_area_wgt_var_ras_sf, aes(fill=species_area_wgt_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site sensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

#p6 <- plot(CAZ_area_wgt_KBA_var_ras, breaks = leg$values, col = leg$colors)
CAZ_area_wgt_KBA_var_ras_st <- st_as_stars(CAZ_area_wgt_KBA_var_ras)
CAZ_area_wgt_KBA_var_ras_sf <- st_as_sf(CAZ_area_wgt_KBA_var_ras_st)
p6 <- ggplot()+
  geom_sf(data = CAZ_area_wgt_KBA_var_ras_sf, aes(fill=species_area_wgt_CAZ_KBA.CAZ_ME.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site sensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

#Comparing rank differences
coolwarm_hcl <- colorspace::diverging_hcl(11,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))

#Couple of examples
species_only_diff <- CAZ_wgt_var_ras - CAZ_var_ras
diff_vals <- values(species_only_diff)
hist(diff_vals)
#plot(species_only_diff, col = coolwarm_hcl) #Red means weighted gave higher priority over unweighted, blue means it gave lower priority
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
######################################################################
#Map for inset (species + weights)
for_inset <- ggplot()+
  geom_sf(data = CAZ_wgt_var_ras_sf, aes(fill=species_wgt_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site sensitivity",
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

                                    ##IAS SDMs##
#Inspect the distribution of priorities within the IAS predicted distribution
#Use binary versions of ensemble committee averaging model
regional_model_path <- "C:/Users/David Clarke.DESKTOP-NNNNVLL/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/IAS_distributions/IAS_regional"

                              ##CAZ with weights##
#Pheidole megacephala
PM_bin <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PM_bin[PM_bin == 0] <- NA
PM_bin <- resample(PM_bin, CAZ_wgt_var_ras) #get equal extent
#PM_hist <- plot_hist(x = CAZ_wgt_var_ras, PM_bin, add.median = TRUE) #approximately left skewed
############################################################
PM_temp_r <- mask(CAZ_wgt_var_ras, PM_bin)
PM_temp_r_values <- values(PM_temp_r)
PM_temp_df <- data.frame(data = PM_temp_r_values)
PM_hist <- ggplot(PM_temp_df, aes(x = PM_temp_r_values)) + 
  geom_histogram(colour = "white",
                 binwidth = 0.05) + 
  scale_x_continuous(breaks = seq(0, 1, 0.25)) + 
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,800)) +
  xlab("Site sensitivity") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 12)) +
  geom_vline(xintercept = median(PM_temp_r_values, #could also use mean
                                 na.rm = T), colour = "red")
#######################################################################
PM_bin_r <- PM_bin
PM_bin_r <- rasterToPolygons(PM_bin_r, dissolve = T)
PM_bin_cr <- crop(PM_bin, PM_bin_r)
PM_CAZ_cr <- crop(CAZ_wgt_var_ras, PM_bin_cr)
PM_CAZ_wgt_st <- st_as_stars(PM_CAZ_cr)
PM_CAZ_wgt_sf <- st_as_sf(PM_CAZ_wgt_st)
PM_bin_sf <- st_as_sf(PM_bin_r)
PM_bin_sf_bb <- st_as_sfc(st_bbox(PM_bin_sf))
PM_plot <- ggplot()+
  geom_sf(data = PM_CAZ_wgt_sf, aes(fill=species_wgt_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0001, 0.25, 0.5, 0.75, 0.9999),
                       labels = as.character(c(0.0, 0.25, 0.5, 0.75, 1))) +
  geom_sf(data = PM_bin_sf, show.legend = F,fill=alpha("black",0.3)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.border = element_rect(colour = "blue", size = 1))
  
#Icerya purchasi
IP_bin <- raster(file.path(regional_model_path, "Icerya.purchasi", "proj_regional", "individual_projections", paste0("Icerya.purchasi","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
IP_bin[IP_bin == 0] <- NA
IP_bin <- resample(IP_bin, CAZ_wgt_var_ras) #get equal extent
############################################################
IP_temp_r <- mask(CAZ_wgt_var_ras, IP_bin)
IP_temp_r_values <- values(IP_temp_r)
IP_temp_df <- data.frame(data = IP_temp_r_values)
IP_hist <- ggplot(IP_temp_df, aes(x = IP_temp_r_values)) + 
  geom_histogram(colour = "white",
                 binwidth = 0.05) + 
  scale_x_continuous(breaks = seq(0, 1, 0.25)) + 
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,12000)) +
  xlab("Site sensitivity") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 12)) +
  geom_vline(xintercept = median(IP_temp_r_values, #could also use mean
                                 na.rm = T), colour = "red")
#######################################################################
IP_bin_r <- IP_bin
IP_bin_r <- rasterToPolygons(IP_bin_r, dissolve = T)
IP_bin_cr <- crop(IP_bin, IP_bin_r)
IP_CAZ_cr <- crop(CAZ_wgt_var_ras, IP_bin_cr)
IP_CAZ_wgt_st <- st_as_stars(IP_CAZ_cr)
IP_CAZ_wgt_sf <- st_as_sf(IP_CAZ_wgt_st)
IP_bin_sf <- st_as_sf(IP_bin_r)
IP_bin_sf_bb <- st_as_sfc(st_bbox(IP_bin_sf))
IP_plot <- ggplot()+
  geom_sf(data = IP_CAZ_wgt_sf, aes(fill=species_wgt_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0001, 0.25, 0.5, 0.75, 0.9999),
                       labels = as.character(c(0.0, 0.25, 0.5, 0.75, 1))) +
  geom_sf(data = IP_bin_sf, show.legend = F,fill=alpha("black",0.3)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.border = element_rect(colour = "black", size = 1))


#Vespula germanica
VG_bin <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VG_bin[VG_bin == 0] <- NA
VG_bin <- resample(VG_bin, CAZ_wgt_var_ras) #get equal extent
#plot_hist(x = CAZ_wgt_var_ras, VG_bin, add.median = TRUE)
############################################################
VG_temp_r <- mask(CAZ_wgt_var_ras, VG_bin)
VG_temp_r_values <- values(VG_temp_r)
VG_temp_df <- data.frame(data = VG_temp_r_values)
VG_hist <- ggplot(VG_temp_df, aes(x = VG_temp_r_values)) + 
  geom_histogram(colour = "white",
                 binwidth = 0.05) + 
  scale_x_continuous(breaks = seq(0, 1, 0.25)) + 
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,7000)) +
  xlab("Site sensitivity") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 12)) +
  geom_vline(xintercept = median(VG_temp_r_values, #could also use mean
                                 na.rm = T), colour = "red")
#######################################################################
VG_bin_r <- VG_bin
VG_bin_r <- rasterToPolygons(VG_bin_r, dissolve = T)
VG_bin_cr <- crop(VG_bin, VG_bin_r)
VG_CAZ_cr <- crop(CAZ_wgt_var_ras, VG_bin_cr)
VG_CAZ_wgt_st <- st_as_stars(VG_CAZ_cr)
VG_CAZ_wgt_sf <- st_as_sf(VG_CAZ_wgt_st)
VG_bin_sf <- st_as_sf(VG_bin_r)
VG_bin_sf_bb <- st_as_sfc(st_bbox(VG_bin_sf))
VG_plot <- ggplot()+
  geom_sf(data = VG_CAZ_wgt_sf, aes(fill=species_wgt_CAZ.CAZ_E.rank.compressed), 
          color=NA, 
          show.legend = T) + 
  scale_fill_gradientn(colours = leg$colors,
                       values = leg$values,
                       name = "Site\nsensitivity",
                       breaks = c(0.0001, 0.25, 0.5, 0.75, 0.9999),
                       labels = as.character(c(0.0, 0.25, 0.5, 0.75, 1))) +
  geom_sf(data = VG_bin_sf, show.legend = F,fill=alpha("black",0.3)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.border = element_rect(colour = "purple", size = 1))

for_inset_all <- for_inset +
  geom_sf(data = PM_bin_sf_bb, fill = NA, col = "blue", size = 1.2) +
  geom_sf(data = IP_bin_sf_bb, fill = NA, col = "black", size = 1.2) +
  geom_sf(data = VG_bin_sf_bb, fill = NA, col = "purple", size = 1.2)

p1 <- for_inset_all
p2 <- ggarrange(PM_plot, IP_plot, VG_plot, legend = "none", ncol = 3, nrow = 1)
p3 <- ggarrange(PM_hist, IP_hist, VG_hist, ncol = 3, nrow = 1)
ggarrange(p1,p2,p3, ncol = 1, nrow = 3)



####################################### KBAs ######################################
#"Inverted" KBA overlap with IAS (KBA threat status 5 more important than 0)
#Pheidole megacephala
PM_bin_KBA <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PM_bin_KBA[PM_bin_KBA == 0] <- NA
PM_bin_KBA <- resample(PM_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent
plot_hist(x = CAZ_wgt_KBA_inv_var_ras, PM_bin_KBA, add.median = TRUE) #approximately left skewed

#Icerya purchasi
IP_bin_KBA <- raster(file.path(regional_model_path, "Icerya.purchasi", "proj_regional", "individual_projections", paste0("Icerya.purchasi","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
IP_bin_KBA[IP_bin_KBA == 0] <- NA
IP_bin_KBA <- resample(IP_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent
plot_hist(x = CAZ_wgt_KBA_inv_var_ras, IP_bin_KBA, add.median = TRUE)

#Vespula germanica
VG_bin_KBA <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VG_bin_KBA[VG_bin_KBA == 0] <- NA
VG_bin_KBA <- resample(VG_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent
plot_hist(x = CAZ_wgt_KBA_inv_var_ras, VG_bin_KBA, add.median = TRUE)

#Comparing distributions
PM_KBA_temp_r <- mask(CAZ_wgt_KBA_inv_var_ras, PM_bin_KBA)
PM_KBA_temp_r_values <- values(PM_KBA_temp_r)
PM_KBA_temp_r_values <- na.omit(PM_KBA_temp_r_values)
PM_temp_r <- mask(CAZ_wgt_var_ras, PM_bin)
PM_temp_r_values <- values(PM_temp_r)
PM_temp_r_values <- na.omit(PM_temp_r_values)

IP_KBA_temp_r <- mask(CAZ_wgt_KBA_inv_var_ras, IP_bin_KBA)
IP_KBA_temp_r_values <- values(IP_KBA_temp_r)
IP_KBA_temp_r_values <- na.omit(IP_KBA_temp_r_values)
IP_temp_r <- mask(CAZ_wgt_var_ras, IP_bin)
IP_temp_r_values <- values(IP_temp_r)
IP_temp_r_values <- na.omit(IP_temp_r_values)

VG_KBA_temp_r <- mask(CAZ_wgt_KBA_inv_var_ras, VG_bin_KBA)
VG_KBA_temp_r_values <- values(VG_KBA_temp_r)
VG_KBA_temp_r_values <- na.omit(VG_KBA_temp_r_values)
VG_temp_r <- mask(CAZ_wgt_var_ras, VG_bin)
VG_temp_r_values <- values(VG_temp_r)
VG_temp_r_values <- na.omit(VG_temp_r_values)

#Kolmogorov-smirnoff tests - KBA vs no KBA
#a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution
#presence of ties a result of rounding
ks.test(PM_KBA_temp_r_values, PM_temp_r_values, alternative = "two.sided")
ks.test(IP_KBA_temp_r_values, IP_temp_r_values, alternative = "two.sided")
ks.test(VG_KBA_temp_r_values, VG_temp_r_values, alternative = "two.sided")

#Difference between species + weight IAS
ks.test(IP_temp_r_values, PM_temp_r_values, alternative = "two.sided")
ks.test(VG_temp_r_values, IP_temp_r_values, alternative = "two.sided")
ks.test(PM_temp_r_values, VG_temp_r_values, alternative = "two.sided")

#Proportion difference (KBA vs no KBA) in number of top sensitive sites
PM_sen_bin <- PM_temp_r_values
PM_sen_bin <- ifelse(PM_sen_bin >= 0.98, 1,0)
sum(PM_sen_bin)/length(PM_sen_bin)*100
PM_KBA_sen_bin <- PM_KBA_temp_r_values
PM_KBA_sen_bin <- ifelse(PM_KBA_sen_bin >= 0.98, 1,0)
sum(PM_KBA_sen_bin)/length(PM_KBA_sen_bin)*100

IP_sen_bin <- IP_temp_r_values
IP_sen_bin <- ifelse(IP_sen_bin >= 0.98, 1,0)
sum(IP_sen_bin)/length(IP_sen_bin)*100
IP_KBA_sen_bin <- IP_KBA_temp_r_values
IP_KBA_sen_bin <- ifelse(IP_KBA_sen_bin >= 0.98, 1,0)
sum(IP_KBA_sen_bin)/length(IP_KBA_sen_bin)*100

VG_sen_bin <- VG_temp_r_values
VG_sen_bin <- ifelse(VG_sen_bin >= 0.98, 1,0)
sum(VG_sen_bin)/length(VG_sen_bin)*100
VG_KBA_sen_bin <- VG_KBA_temp_r_values
VG_KBA_sen_bin <- ifelse(VG_KBA_sen_bin >= 0.98, 1,0)
sum(VG_KBA_sen_bin)/length(VG_KBA_sen_bin)*100

#Top quarter (i.e. >0.75 sensitivity)
PM_sen_bin_75 <- PM_temp_r_values
PM_sen_bin_75 <- ifelse(PM_sen_bin_75 >= 0.75, 1,0)
sum(PM_sen_bin_75)/length(PM_sen_bin_75)*100
PM_KBA_sen_bin_75 <- PM_KBA_temp_r_values
PM_KBA_sen_bin_75 <- ifelse(PM_KBA_sen_bin_75 >= 0.75, 1,0)
sum(PM_KBA_sen_bin_75)/length(PM_KBA_sen_bin_75)*100

IP_sen_bin_75 <- IP_temp_r_values
IP_sen_bin_75 <- ifelse(IP_sen_bin_75 >= 0.75, 1,0)
sum(IP_sen_bin_75)/length(IP_sen_bin_75)*100
IP_KBA_sen_bin_75 <- IP_KBA_temp_r_values
IP_KBA_sen_bin_75 <- ifelse(IP_KBA_sen_bin_75 >= 0.75, 1,0)
sum(IP_KBA_sen_bin_75)/length(IP_KBA_sen_bin_75)*100

VG_sen_bin_75 <- VG_temp_r_values
VG_sen_bin_75 <- ifelse(VG_sen_bin_75 >= 0.75, 1,0)
sum(VG_sen_bin_75)/length(VG_sen_bin_75)*100
VG_KBA_sen_bin_75 <- VG_KBA_temp_r_values
VG_KBA_sen_bin_75 <- ifelse(VG_KBA_sen_bin_75 >= 0.75, 1,0)
sum(VG_KBA_sen_bin_75)/length(VG_KBA_sen_bin_75)*100

#################################### Jaccards similarity #############################
#Variant names
all_names <- c("species_CAZ",
               "species_wgt_CAZ",
               "species_wgt_CAZ_KBA",
               "species_area_CAZ",
               "species_area_wgt_CAZ",
               "species_area_wgt_CAZ_KBA")

calculate_jaccards <- function(rank_stack, x.min, x.max, y.min, y.max, variant_names) {
  
  jaccards <- matrix(nrow = nlayers(rank_stack), ncol = nlayers(rank_stack))
  for (i in 1:nrow(jaccards)) {
    for (j in 1:ncol(jaccards)) {
      if (i == j) {
        jaccards[i, j] <- 1
      }
      else {
        if (is.na(jaccards[j, i])) {
          message(paste0("Calculating Jaccard index between ", 
                         names(rank_stack[[i]]), 
                         " and ", names(rank_stack[[j]])))
          jaccards[i, j] <- jaccard(rank_stack[[i]], rank_stack[[j]], 
                                    x.min = x.min, x.max = x.max, 
                                    y.min = y.min, y.max = y.max)
        }
        else {
          jaccards[i, j] <- NA
        }
      }
    }
  }
  jaccards <- as.data.frame(jaccards)
  colnames(jaccards) <- variant_names
  rownames(jaccards) <- variant_names
  return(jaccards)
}

#Top 2%
top_two <- calculate_jaccards(CAZ_stack, x.min = 0.98, x.max = 1.0, 
                              y.min = 0.98,y.max = 1.0,all_names)

#Top 5%
top_five <- calculate_jaccards(CAZ_stack, x.min = 0.95, x.max = 1.0, 
                               y.min = 0.95,y.max = 1.0,all_names)

#Top 10%
top_ten <- calculate_jaccards(CAZ_stack, x.min = 0.9, x.max = 1.0, 
                              y.min = 0.9,y.max = 1.0,all_names)

#Top 25%
top_twentyfive <- calculate_jaccards(CAZ_stack, x.min = 0.75, x.max = 1.0, 
                                     y.min = 0.75,y.max = 1.0,all_names)

#Top 50%
top_fifty <- calculate_jaccards(CAZ_stack, x.min = 0.5, x.max = 1.0, 
                                y.min = 0.5,y.max = 1.0,all_names)

#Top 80%
top_eighty <- calculate_jaccards(CAZ_stack, x.min = 0.2, x.max = 1.0, 
                                 y.min = 0.2,y.max = 1.0,all_names)

#Total
total <- calculate_jaccards(CAZ_stack, x.min = 0.0, x.max = 1.0, 
                            y.min = 0.0,y.max = 1.0,all_names)

write.csv(top_two, file = file.path("Zonation", "jaccard_two.csv"), na = "-", row.names = T)
write.csv(top_five, file = file.path("Zonation", "jaccard_five.csv"), na = "-", row.names = T)
write.csv(top_ten, file = file.path("Zonation", "jaccard_ten.csv"), na = "-", row.names = T)
write.csv(top_twentyfive, file = file.path("Zonation", "jaccard_twentyfive.csv"), na = "-", row.names = T)
write.csv(top_fifty, file = file.path("Zonation", "jaccard_fifty.csv"), na = "-", row.names = T)
write.csv(top_eighty, file = file.path("Zonation", "jaccard_eighty.csv"), na = "-", row.names = T)
write.csv(total, file = file.path("Zonation", "jaccard_total.csv"), na = "-", row.names = T)

######################################## Distance to coast ##############################
#Prior to doing any kind of distance work, need to project everything to GDA94
CAZ_var_ras_proj <- projectRaster(CAZ_var_ras, res = 5000, crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
CAZ_wgt_var_ras_proj <- projectRaster(CAZ_wgt_var_ras, res = 5000, crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
CAZ_area_var_ras_proj <- projectRaster(CAZ_area_var_ras, res = 5000, crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
CAZ_area_wgt_var_ras_proj <- projectRaster(CAZ_area_wgt_var_ras, res = 5000, crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#For each variant, produce plot and calculate correlation, between top fraction sites and distance to coast
#load distance raster
dist_coast <- raster(file.path("SpatialData", "Raster", "dist_aus_coast.tif"))
dist_coast <- projectRaster(dist_coast, res = 5000, crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
dists <- values(dist_coast)

#extract site priority at each distance
dist_coord <- as.data.frame(coordinates(dist_coast))
pri.list <- list(CAZ_var_ras_proj, CAZ_wgt_var_ras_proj, CAZ_area_var_ras_proj, CAZ_area_wgt_var_ras_proj)
squeeze <- function(pvector){
  for(i in 1:length(pvector)) {
    if(!is.na(pvector[i]) == T & pvector[i] == 1){
      pvector[i] <- pvector[i] - 0.000001
    }
  }
  return(pvector)
}

priority_dist <- data.frame(lapply(pri.list, function(i){
priority <- raster::extract(i, dist_coord)
priority <- squeeze(priority)
names(priority) <- names(i)
return(priority)
}), dist_coord)
colnames(priority_dist) <- c("species_CAZ","species_wgt_CAZ", "species_area_CAZ", "species_area_wgt_CAZ", "longitude", "latitude")

################## Looking for spatial pattern of highest sensitive sites #########
spat_priority_dist <- function(df, n_col){
  for(i in 1:n_col){
    df[,i] <- ifelse(df[,i] >= 0.98, 1,0)}
  return(df)
}
df_new <- spat_priority_dist(priority_dist, 4)
#E.g.variant 1. Do for all 4 in df
df_1 <- data.frame(df_new[,1], df_new[,5], df_new[,6])
df_1 <- df_1 %>% dplyr::filter(df_1[,1] == 1)
p <- ppp(df_1[,2], df_1[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")

df_2 <- data.frame(df_new[,2], df_new[,5], df_new[,6])
df_2 <- df_2 %>% dplyr::filter(df_2[,1] == 1)
p <- ppp(df_2[,2], df_2[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")

df_3 <- data.frame(df_new[,3], df_new[,5], df_new[,6])
df_3 <- df_3 %>% dplyr::filter(df_3[,1] == 1)
p <- ppp(df_3[,2], df_3[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")

df_4 <- data.frame(df_new[,4], df_new[,5], df_new[,6])
df_4 <- df_4 %>% dplyr::filter(df_4[,1] == 1)
p <- ppp(df_4[,2], df_4[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")


#create data frame of priority and distance
#dist_df <- data.frame(priority = priority, distance = dists)

#correlation
cor.test(priority_dist[,1], dists, method = "spearman")
cor.test(priority_dist[,2], dists, method = "spearman")
cor.test(priority_dist[,3], dists, method = "spearman")
cor.test(priority_dist[,4], dists, method = "spearman")

#scatterplot
plot(dists, priority, cex = 0.1)
