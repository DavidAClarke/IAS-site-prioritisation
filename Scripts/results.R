#Results
#Required libraries
library(zonator)
library(tidyverse)
library(spatstat)
library(sf)
library(stars)
library(raster)
library(ggpubr)
library(cowplot)
library(RColorBrewer)

#Call functions
source("Scripts/results_functions.R")

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
IAS_threats <- IAS_threats %>%
  mutate(classGroup = factor(classGroup, levels=c("Mammal", "Bird", "Reptile", "Amphibian", "Fish", "Invertebrate", "Plant")))
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

p1 <- rank_plot(CAZ_var_ras)
p2 <- rank_plot(CAZ_wgt_var_ras)
p3 <- rank_plot(CAZ_wgt_KBA_var_ras)
p4 <- rank_plot(CAZ_area_var_ras)
p5 <- rank_plot(CAZ_area_wgt_var_ras)
p6 <- rank_plot(CAZ_area_wgt_KBA_var_ras)

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
######################################################################

                                    ##IAS SDMs##
regional_model_path <- "F:/PhD/Chapter_3/Data/IAS_distributions/IAS_regional"
spp_list <- c("Digitonthophagus gazella", "Pheidole megacephala",
              "Vespula germanica", "Tetramorium bicarinatum", "Paratrechina longicornis")

#Plots of predicted distributions (ensemble committee averaging)
lapply(spp_list[1:length(spp_list)], FUN = function(i){
  
  IAS_plot(i)
  
})

#Inspect the distribution of priorities within the IAS predicted distribution
#Use binary versions of ensemble committee averaging model


                              ##CAZ with weights##
#Pheidole megacephala
PM_bin <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
PM_bin2 <- PM_bin
PM_bin2[PM_bin2 != 0] <- 1

# Re-class zeros as NA
PM_bin[PM_bin == 0] <- NA
PM_bin <- resample(PM_bin, CAZ_wgt_var_ras) #get equal extent
PM_bin_p <- rasterToPolygons(PM_bin, dissolve = T)
PM_bin_sf <- st_as_sf(PM_bin_p)
PM_bin_sf_bb <- st_as_sfc(st_bbox(PM_bin_sf))


#Vespula germanica
VG_bin <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VG_bin2 <- VG_bin
VG_bin2[VG_bin2 != 0] <- 1

# Re-class zeros as NA
VG_bin[VG_bin == 0] <- NA
VG_bin <- resample(VG_bin, CAZ_wgt_var_ras) #get equal extent
VG_bin_p <- rasterToPolygons(VG_bin, dissolve = T)
VG_bin_sf <- st_as_sf(VG_bin_p)
VG_bin_sf_bb <- st_as_sfc(st_bbox(VG_bin_sf))

#Digitonthophagus gazella
DG_bin <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
DG_bin2 <- DG_bin
DG_bin2[DG_bin2 != 0] <- 1

# Re-class zeros as NA
DG_bin[DG_bin == 0] <- NA
DG_bin <- resample(DG_bin, CAZ_wgt_var_ras) #get equal extent
DG_bin_p <- rasterToPolygons(DG_bin, dissolve = T)
DG_bin_sf <- st_as_sf(DG_bin_p)
DG_bin_sf_bb <- st_as_sfc(st_bbox(DG_bin_sf))

#Tetramorium bicarinatum
TB_bin <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
TB_bin2 <- TB_bin
TB_bin2[TB_bin2 != 0] <- 1

# Re-class zeros as NA
TB_bin[TB_bin == 0] <- NA
TB_bin <- resample(TB_bin, CAZ_wgt_var_ras) #get equal extent
TB_bin_p <- rasterToPolygons(TB_bin, dissolve = T)
TB_bin_sf <- st_as_sf(TB_bin_p)
TB_bin_sf_bb <- st_as_sfc(st_bbox(TB_bin_sf))

#Paratrechina longicornis
PL_bin <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
PL_bin2 <- PL_bin
PL_bin2[PL_bin2 != 0] <- 1

# Re-class zeros as NA
PL_bin[PL_bin == 0] <- NA
PL_bin <- resample(PL_bin, CAZ_wgt_var_ras) #get equal extent
PL_bin_p <- rasterToPolygons(PL_bin, dissolve = T)
PL_bin_sf <- st_as_sf(PL_bin_p)
PL_bin_sf_bb <- st_as_sfc(st_bbox(PL_bin_sf))

#Map for inset (species + weights)
CAZ_wgt_var_ras_st <- st_as_stars(CAZ_wgt_var_ras)
CAZ_wgt_var_ras_sf <- st_as_sf(CAZ_wgt_var_ras_st)
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
        axis.ticks = element_blank(),
        legend.position = "top")
#to add rectangles for bbox of IAS sdm
map_w_borders <- for_inset + 
                  geom_sf(data = PM_bin_sf_bb, fill = NA, color = "#5bb4f2") +
                  geom_sf(data = VG_bin_sf_bb, fill = NA, color = "#9c5200") +
                  geom_sf(data = DG_bin_sf_bb, fill = NA, color = "#b409a7") +
                  geom_sf(data = TB_bin_sf_bb, fill = NA, color = "#040200") +
                  geom_sf(data = PL_bin_sf_bb, fill = NA, color = "#0db02f")

#Histograms of IAS overlap with sensitive sites
PM_hist <- mask_hist(CAZ_wgt_var_ras, PM_bin)
VG_hist <- mask_hist(CAZ_wgt_var_ras, VG_bin)
DG_hist <- mask_hist(CAZ_wgt_var_ras, DG_bin)
TB_hist <- mask_hist(CAZ_wgt_var_ras, TB_bin)
PL_hist <- mask_hist(CAZ_wgt_var_ras, PL_bin)


#Priority plots - overlap of sensitive and susceptible sites (predicted IAS range)
#sdm_col = colour of the predicted range, sdm_brd_col = colour of rectangle delineating range extent
PM_priority <- Priority_plot(PM_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#5bb4f2") #use sdm_brd_col = NA if you don't want a rectangle
VG_priority <- Priority_plot(VG_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#9c5200")
DG_priority <- Priority_plot(DG_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
TB_priority <- Priority_plot(TB_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
PL_priority <- Priority_plot(PL_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")


#Overlap among all three species (turn into function?)
ias_sum <- sum(PM_bin2, VG_bin2, DG_bin2, TB_bin2, PL_bin2)
ias_sum_one <- ias_sum
ias_sum_one[ias_sum_one != 1] <- NA


####################################### KBAs ######################################
#"Inverted" KBA overlap with IAS (KBA threat status 5 more important than 0)

#Pheidole megacephala
PM_bin_KBA <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PM_bin_KBA[PM_bin_KBA == 0] <- NA
PM_bin_KBA <- resample(PM_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Vespula germanica
VG_bin_KBA <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VG_bin_KBA[VG_bin_KBA == 0] <- NA
VG_bin_KBA <- resample(VG_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Digitonthophagus gazella
DG_bin_KBA <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
DG_bin_KBA[DG_bin_KBA == 0] <- NA
DG_bin_KBA <- resample(DG_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Tetramorium bicarinatum
TB_bin_KBA <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
TB_bin_KBA[TB_bin_KBA == 0] <- NA
TB_bin_KBA <- resample(TB_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Paratrechina longicornis
PL_bin_KBA <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PL_bin_KBA[PL_bin_KBA == 0] <- NA
PL_bin_KBA <- resample(PL_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Histograms of IAS overlap with sensitive sites
PM_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, PM_bin_KBA)
VG_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, VG_bin_KBA)
DG_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, DG_bin_KBA)
TB_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, TB_bin_KBA)
PL_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, PL_bin_KBA)

#Get cell values for Kolmogorov-smirnoff tests
PM_vals <- get_msk_vals(CAZ_wgt_var_ras, PM_bin)
PM_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, PM_bin_KBA)

VG_vals <- get_msk_vals(CAZ_wgt_var_ras, VG_bin)
VG_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, VG_bin_KBA)

DG_vals <- get_msk_vals(CAZ_wgt_var_ras, DG_bin)
DG_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, DG_bin_KBA)

TB_vals <- get_msk_vals(CAZ_wgt_var_ras, TB_bin)
TB_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, TB_bin_KBA)

PL_vals <- get_msk_vals(CAZ_wgt_var_ras, PL_bin)
PL_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, PL_bin_KBA)


#Kolmogorov-smirnoff tests - KBA vs no KBA
#a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution
#presence of ties a result of rounding
ks.test(PM_KBA_vals, PM_vals, alternative = "two.sided")
ks.test(VG_KBA_vals, VG_vals, alternative = "two.sided")
ks.test(DG_KBA_vals, DG_vals, alternative = "two.sided")
ks.test(TB_KBA_vals, TB_vals, alternative = "two.sided")
ks.test(PL_KBA_vals, PL_vals, alternative = "two.sided")


#Difference between species + weight IAS
PM_VG <- ks.test(PM_vals, VG_vals, alternative = "two.sided")
PM_DG <- ks.test(PM_vals, DG_vals, alternative = "two.sided")
PM_TB <- ks.test(PM_vals, TB_vals, alternative = "two.sided")
PM_PL <- ks.test(PM_vals, PL_vals, alternative = "two.sided")

VG_DG <- ks.test(VG_vals, DG_vals, alternative = "two.sided")
VG_PL <- ks.test(VG_vals, PL_vals, alternative = "two.sided")
VG_TB <- ks.test(VG_vals, TB_vals, alternative = "two.sided")

DG_TB <- ks.test(DG_vals, TB_vals, alternative = "two.sided")
DG_PL <- ks.test(DG_vals, PL_vals, alternative = "two.sided")

TB_PL <- ks.test(TB_vals, PL_vals, alternative = "two.sided")

v <- c(0,PM_VG$statistic,PM_DG$statistic,PM_TB$statistic, PM_PL$statistic,
       0,0,VG_DG$statistic,VG_PL$statistic,VG_TB$statistic,
       0,0,0,DG_TB$statistic,DG_PL$statistic,
       0,0,0,0,TB_PL$statistic,)
tm <- matrix(v, nrow = 5, ncol = 4)
rownames(tm) <- c("Pm", "Vg", "Dg", "Tb", "Pl")
colnames(tm) <- c("Pm", "Vg", "Dg", "Tb")
corrplot::corrplot(tm, type = "lower", method = "color", 
                   cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
                   tl.srt = 0, tl.col = "black", addCoef.col = "black")

#Proportion difference (KBA vs no KBA) in number of top sensitive sites
#Top two (i.e. >= 0.98 sensitivity)
props <- c(0.98, 0.95, 0.90, 0.75, 0.50, 0.25, 0.00)
diffs <- c()
PM_prop <- multi_props(PM_vals, props)
PM_KBA_prop <- multi_props(PM_KBA_vals, props)

VG_prop <- multi_props(VG_vals, props)
VG_KBA_prop <- multi_props(VG_KBA_vals, props)

DG_prop <- multi_props(DG_vals, props)
DG_KBA_prop <- multi_props(DG_KBA_vals, props)

TB_prop <- multi_props(TB_vals, props)
TB_KBA_prop <- multi_props(TB_KBA_vals, props)

PL_prop <- multi_props(VG_vals, props)
PL_KBA_prop <- multi_props(VG_KBA_vals, props)

Total <- as.data.frame(rbind(PM_prop, VG_prop, DG_prop, TB_prop, PL_prop,
           PM_KBA_prop, VG_KBA_prop, DG_KBA_prop, TB_KBA_prop, PL_KBA_prop))
Type <- c(rep("Unmasked", 5), rep("Masked", 5))
Species <- rep(c("P. megacephala", "V. germanica", "D. gazella",
                 "T. bicarinatum", "P. longicornis"), 2)
Total <- cbind(Species, Type,Total)
Total <- Total %>%
  as_tibble() %>%
  mutate(Species = factor(Species)) %>%
  mutate(Type = factor(Type)) %>%
  pivot_longer(cols = !Species & !Type, 
               names_to = "SiteSensitivity", 
               values_to = "DistributionCoverage") %>%
  mutate(SiteSensitivity = as.double(SiteSensitivity))
cols <- c("#b409a7","#0db02f", "#5bb4f2","#040200","#9c5200")
ggplot(Total, 
       aes(x = SiteSensitivity, 
           y = DistributionCoverage, 
           colour = Species,
           group = interaction(Species,Type))) +
  geom_line(aes(linetype = Type)) +
  scale_colour_manual(values = cols) +
  ylab("Distribution Coverage (%)") +
  xlab("Site Sensitivity") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1,0),
                     trans = "reverse") +
  scale_y_continuous(expand = c(0,0))
  


#################################### Jaccards similarity #############################
#Variant names
all_names <- c("species_CAZ",
               "species_wgt_CAZ",
               "species_wgt_CAZ_KBA",
               "species_area_CAZ",
               "species_area_wgt_CAZ",
               "species_area_wgt_CAZ_KBA")

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
priority_dist <- data.frame(lapply(pri.list, function(i){
  priority <- raster::extract(i, dist_coord)
  priority <- squeeze(priority)
  names(priority) <- names(i)
  return(priority)
}), dist_coord)
colnames(priority_dist) <- c("species_CAZ","species_wgt_CAZ", "species_area_CAZ", "species_area_wgt_CAZ", "longitude", "latitude")

################## Looking for spatial pattern of highest sensitive sites #########
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

# #correlation
# cor.test(priority_dist[,1], dists, method = "spearman")
# cor.test(priority_dist[,2], dists, method = "spearman")
# cor.test(priority_dist[,3], dists, method = "spearman")
# cor.test(priority_dist[,4], dists, method = "spearman")
# 
# #scatterplot
# plot(dists, priority, cex = 0.1)
