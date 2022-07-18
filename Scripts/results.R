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
library(PNWColors)

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

source("Scripts/jaccard_similarities.R")

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
######################################################################

                                    ##IAS SDMs##
regional_model_path <- "G:/Chapter_3/SpatialData/IAS_distributions/IAS_regional" #external hard drive
spp_list <- c("Apis mellifera",  "Monomorium floricola", "Monomorium destructor", 
              "Linepithema humile", "Vespula vulgaris", "Bombus terrestris", "Heteronychus arator", 
              "Digitonthophagus gazella", "Pheidole megacephala", "Vespula germanica", 
              "Tetramorium bicarinatum", "Paratrechina longicornis")
not_run <- c("Solenopsis geminata","Polistes chinensis antennalis","Wasmannia auropunctata", "Apis cerana",
             "Solenopsis invicta","Megachile rotundata")

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
# PM_bin_p <- rasterToPolygons(PM_bin, dissolve = T)
# PM_bin_sf <- st_as_sf(PM_bin_p)
# PM_bin_sf_bb <- st_as_sfc(st_bbox(PM_bin_sf))


#Vespula germanica
VG_bin <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VG_bin2 <- VG_bin
VG_bin2[VG_bin2 != 0] <- 1

# Re-class zeros as NA
VG_bin[VG_bin == 0] <- NA
VG_bin <- resample(VG_bin, CAZ_wgt_var_ras) #get equal extent
# VG_bin_p <- rasterToPolygons(VG_bin, dissolve = T)
# VG_bin_sf <- st_as_sf(VG_bin_p)
# VG_bin_sf_bb <- st_as_sfc(st_bbox(VG_bin_sf))

#Digitonthophagus gazella
DG_bin <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
DG_bin2 <- DG_bin
DG_bin2[DG_bin2 != 0] <- 1

# Re-class zeros as NA
DG_bin[DG_bin == 0] <- NA
DG_bin <- resample(DG_bin, CAZ_wgt_var_ras) #get equal extent
# DG_bin_p <- rasterToPolygons(DG_bin, dissolve = T)
# DG_bin_sf <- st_as_sf(DG_bin_p)
# DG_bin_sf_bb <- st_as_sfc(st_bbox(DG_bin_sf))

#Tetramorium bicarinatum
TB_bin <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
TB_bin2 <- TB_bin
TB_bin2[TB_bin2 != 0] <- 1

# Re-class zeros as NA
TB_bin[TB_bin == 0] <- NA
TB_bin <- resample(TB_bin, CAZ_wgt_var_ras) #get equal extent
# TB_bin_p <- rasterToPolygons(TB_bin, dissolve = T)
# TB_bin_sf <- st_as_sf(TB_bin_p)
# TB_bin_sf_bb <- st_as_sfc(st_bbox(TB_bin_sf))

#Paratrechina longicornis
PL_bin <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
PL_bin2 <- PL_bin
PL_bin2[PL_bin2 != 0] <- 1

# Re-class zeros as NA
PL_bin[PL_bin == 0] <- NA
PL_bin <- resample(PL_bin, CAZ_wgt_var_ras) #get equal extent
# PL_bin_p <- rasterToPolygons(PL_bin, dissolve = T)
# PL_bin_sf <- st_as_sf(PL_bin_p)
# PL_bin_sf_bb <- st_as_sfc(st_bbox(PL_bin_sf))

#Apis mellifera
AM_bin <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
AM_bin2 <- AM_bin
AM_bin2[AM_bin2 != 0] <- 1

# Re-class zeros as NA
AM_bin[AM_bin == 0] <- NA
AM_bin <- resample(AM_bin, CAZ_wgt_var_ras) #get equal extent
# AM_bin_p <- rasterToPolygons(AM_bin, dissolve = T)
# AM_bin_sf <- st_as_sf(AM_bin_p)
# AM_bin_sf_bb <- st_as_sfc(st_bbox(AM_bin_sf))

#Monomorium floricola 
MF_bin <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MF_bin2 <- MF_bin
MF_bin2[MF_bin2 != 0] <- 1

# Re-class zeros as NA
MF_bin[PL_bin == 0] <- NA
MF_bin <- resample(MF_bin, CAZ_wgt_var_ras) #get equal extent
# MF_bin_p <- rasterToPolygons(MF_bin, dissolve = T)
# MF_bin_sf <- st_as_sf(MF_bin_p)
# MF_bin_sf_bb <- st_as_sfc(st_bbox(MF_bin_sf))

#Monomorium destructor
MD_bin <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MD_bin2 <- MD_bin
MD_bin2[MD_bin2 != 0] <- 1

# Re-class zeros as NA
MD_bin[PL_bin == 0] <- NA
MD_bin <- resample(MD_bin, CAZ_wgt_var_ras) #get equal extent
# MD_bin_p <- rasterToPolygons(MD_bin, dissolve = T)
# MD_bin_sf <- st_as_sf(MD_bin_p)
# MD_bin_sf_bb <- st_as_sfc(st_bbox(MD_bin_sf))

#Linepithema humile
LH_bin <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
LH_bin2 <- LH_bin
LH_bin2[LH_bin2 != 0] <- 1

# Re-class zeros as NA
LH_bin[LH_bin == 0] <- NA
LH_bin <- resample(LH_bin, CAZ_wgt_var_ras) #get equal extent
# LH_bin_p <- rasterToPolygons(LH_bin, dissolve = T)
# LH_bin_sf <- st_as_sf(LH_bin_p)
# LH_bin_sf_bb <- st_as_sfc(st_bbox(LH_bin_sf))

#Vespula vulgaris
VV_bin <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VV_bin2 <- VV_bin
VV_bin2[VV_bin2 != 0] <- 1

# Re-class zeros as NA
VV_bin[PL_bin == 0] <- NA
VV_bin <- resample(VV_bin, CAZ_wgt_var_ras) #get equal extent
# VV_bin_p <- rasterToPolygons(VV_bin, dissolve = T)
# VV_bin_sf <- st_as_sf(VV_bin_p)
# VV_bin_sf_bb <- st_as_sfc(st_bbox(VV_bin_sf))

#Megachile rotundata 
# MR_bin <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# MR_bin2 <- MR_bin
# MR_bin2[MR_bin2 != 0] <- 1

# Re-class zeros as NA
# MR_bin[MR_bin == 0] <- NA
# MR_bin <- resample(MR_bin, CAZ_wgt_var_ras) #get equal extent
# MR_bin_p <- rasterToPolygons(MR_bin, dissolve = T)
# MR_bin_sf <- st_as_sf(MR_bin_p)
# MR_bin_sf_bb <- st_as_sfc(st_bbox(MR_bin_sf))

#Bombus terrestris
BT_bin <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
BT_bin2 <- BT_bin
BT_bin2[BT_bin2 != 0] <- 1

# Re-class zeros as NA
BT_bin[BT_bin == 0] <- NA
BT_bin <- resample(BT_bin, CAZ_wgt_var_ras) #get equal extent
# BT_bin_p <- rasterToPolygons(BT_bin, dissolve = T)
# BT_bin_sf <- st_as_sf(BT_bin_p)
# BT_bin_sf_bb <- st_as_sfc(st_bbox(BT_bin_sf))

#Heteronychus arator
HA_bin <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
HA_bin2 <- HA_bin
HA_bin2[HA_bin2 != 0] <- 1

# Re-class zeros as NA
HA_bin[HA_bin == 0] <- NA
HA_bin <- resample(HA_bin, CAZ_wgt_var_ras) #get equal extent
# HA_bin_p <- rasterToPolygons(HA_bin, dissolve = T)
# HA_bin_sf <- st_as_sf(HA_bin_p)
# HA_bin_sf_bb <- st_as_sfc(st_bbox(HA_bin_sf))




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
        legend.position = "right")
#to add rectangles for bbox of IAS sdm
map_w_borders <- for_inset + 
                  geom_sf(data = PM_bin_sf_bb, fill = NA, color = "#5bb4f2") +
                  geom_sf(data = VG_bin_sf_bb, fill = NA, color = "#9c5200") +
                  geom_sf(data = DG_bin_sf_bb, fill = NA, color = "#b409a7") +
                  geom_sf(data = TB_bin_sf_bb, fill = NA, color = "#040200") +
                  geom_sf(data = PL_bin_sf_bb, fill = NA, color = "#0db02f") +
                  geom_sf(data = AM_bin_sf_bb, fill = NA, color = "#5bb4f2") +
                  geom_sf(data = MF_bin_sf_bb, fill = NA, color = "#9c5200") +
                  geom_sf(data = MD_bin_sf_bb, fill = NA, color = "#b409a7") +
                  geom_sf(data = LH_bin_sf_bb, fill = NA, color = "#040200") +
                  geom_sf(data = VV_bin_sf_bb, fill = NA, color = "#0db02f") +
                  geom_sf(data = MR_bin_sf_bb, fill = NA, color = "#0db02f") +
                  geom_sf(data = BT_bin_sf_bb, fill = NA, color = "#0db02f") +
                  geom_sf(data = HA_bin_sf_bb, fill = NA, color = "#0db02f")

#Histograms of IAS overlap with sensitive sites
PM_hist <- mask_hist(CAZ_wgt_var_ras, PM_bin)
VG_hist <- mask_hist(CAZ_wgt_var_ras, VG_bin)
DG_hist <- mask_hist(CAZ_wgt_var_ras, DG_bin)
TB_hist <- mask_hist(CAZ_wgt_var_ras, TB_bin)
PL_hist <- mask_hist(CAZ_wgt_var_ras, PL_bin)
AM_hist <- mask_hist(CAZ_wgt_var_ras, AM_bin)
MF_hist <- mask_hist(CAZ_wgt_var_ras, MF_bin)
MD_hist <- mask_hist(CAZ_wgt_var_ras, MD_bin)
LH_hist <- mask_hist(CAZ_wgt_var_ras, LH_bin)
VV_hist <- mask_hist(CAZ_wgt_var_ras, VV_bin)
#MR_hist <- mask_hist(CAZ_wgt_var_ras, MR_bin)
BT_hist <- mask_hist(CAZ_wgt_var_ras, BT_bin)
HA_hist <- mask_hist(CAZ_wgt_var_ras, HA_bin)


#Priority plots - overlap of sensitive and susceptible sites (predicted IAS range)
#sdm_col = colour of the predicted range, sdm_brd_col = colour of rectangle delineating range extent
PM_priority <- Priority_plot(PM_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#5bb4f2") #use sdm_brd_col = NA if you don't want a rectangle
VG_priority <- Priority_plot(VG_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#9c5200")
DG_priority <- Priority_plot(DG_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
TB_priority <- Priority_plot(TB_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
PL_priority <- Priority_plot(PL_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")
AM_priority <- Priority_plot(AM_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#5bb4f2") #use sdm_brd_col = NA if you don't want a rectangle
MF_priority <- Priority_plot(MF_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#9c5200")
MD_priority <- Priority_plot(MD_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
LH_priority <- Priority_plot(LH_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
VV_priority <- Priority_plot(VV_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")
#MR_priority <- Priority_plot(MR_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
BT_priority <- Priority_plot(BT_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
HA_priority <- Priority_plot(HA_bin, CAZ_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")

                               ##RAN with weights##
#Pheidole megacephala
PM_bin_RAN <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PM_bin_RAN[PM_bin_RAN == 0] <- NA
PM_bin_RAN <- resample(PM_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Vespula germanica
VG_bin_RAN <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VG_bin_RAN[VG_bin_RAN == 0] <- NA
VG_bin_RAN <- resample(VG_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Digitonthophagus gazella
DG_bin_RAN <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
DG_bin_RAN[DG_bin_RAN == 0] <- NA
DG_bin_RAN <- resample(DG_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Tetramorium bicarinatum
TB_bin_RAN <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
TB_bin_RAN[TB_bin_RAN == 0] <- NA
TB_bin_RAN <- resample(TB_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Paratrechina longicornis
PL_bin_RAN <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PL_bin_RAN[PL_bin_RAN == 0] <- NA
PL_bin_RAN <- resample(PL_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Apis mellifera
AM_bin_RAN <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
AM_bin_RAN[AM_bin_RAN == 0] <- NA
AM_bin_RAN <- resample(AM_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Monomorium floricola 
MF_bin_RAN <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MF_bin_RAN[MF_bin_RAN == 0] <- NA
MF_bin_RAN <- resample(MF_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Monomorium destructor
MD_bin_RAN <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MD_bin_RAN[MD_bin_RAN == 0] <- NA
MD_bin_RAN <- resample(MD_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Linepithema humile
LH_bin_RAN <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
LH_bin_RAN[LH_bin_RAN == 0] <- NA
LH_bin_RAN <- resample(LH_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Vespula vulgaris
VV_bin_RAN <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VV_bin_RAN[VV_bin_RAN == 0] <- NA
VV_bin_RAN <- resample(VV_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Megachile rotundata 
# MR_bin_RAN <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# # Re-class zeros as NA
# MR_bin_RAN[MR_bin_RAN == 0] <- NA
# MR_bin_RAN <- resample(MR_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Bombus terrestris
BT_bin_RAN <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
BT_bin_RAN[BT_bin_RAN == 0] <- NA
BT_bin_RAN <- resample(BT_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Heteronychus arator
HA_bin_RAN <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
HA_bin_RAN[HA_bin_RAN == 0] <- NA
HA_bin_RAN <- resample(HA_bin_RAN, RAN_wgt_var_ras) #get equal extent

#Map for inset (random species + weights)
RAN_wgt_var_ras_st <- st_as_stars(RAN_wgt_var_ras)
RAN_wgt_var_ras_sf <- st_as_sf(RAN_wgt_var_ras_st)
for_inset <- ggplot()+
  geom_sf(data = RAN_wgt_var_ras_sf, aes(fill=species_wgt_RAN.RAN_E.rank.compressed), 
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
        legend.position = "right")
#to add rectangles for bbox of IAS sdm
map_w_borders <- for_inset + 
  geom_sf(data = PM_bin_RAN_sf_bb, fill = NA, color = "#5bb4f2") +
  geom_sf(data = VG_bin_RAN_sf_bb, fill = NA, color = "#9c5200") +
  geom_sf(data = DG_bin_RAN_sf_bb, fill = NA, color = "#b409a7") +
  geom_sf(data = TB_bin_RAN_sf_bb, fill = NA, color = "#040200") +
  geom_sf(data = PL_bin_RAN_sf_bb, fill = NA, color = "#0db02f")

#Histograms of IAS overlap with sensitive sites
PM_RAN_hist <- mask_hist(RAN_wgt_var_ras, PM_bin_RAN)
VG_RAN_hist <- mask_hist(RAN_wgt_var_ras, VG_bin_RAN)
DG_RAN_hist <- mask_hist(RAN_wgt_var_ras, DG_bin_RAN)
TB_RAN_hist <- mask_hist(RAN_wgt_var_ras, TB_bin_RAN)
PL_RAN_hist <- mask_hist(RAN_wgt_var_ras, PL_bin_RAN)
AM_RAN_hist <- mask_hist(RAN_wgt_var_ras, AM_bin_RAN)
MF_RAN_hist <- mask_hist(RAN_wgt_var_ras, MF_bin_RAN)
MD_RAN_hist <- mask_hist(RAN_wgt_var_ras, MD_bin_RAN)
LH_RAN_hist <- mask_hist(RAN_wgt_var_ras, LH_bin_RAN)
VV_RAN_hist <- mask_hist(RAN_wgt_var_ras, VV_bin_RAN)
#MR_RAN_hist <- mask_hist(RAN_wgt_var_ras, MR_bin_RAN)
BT_RAN_hist <- mask_hist(RAN_wgt_var_ras, BT_bin_RAN)
HA_RAN_hist <- mask_hist(RAN_wgt_var_ras, HA_bin_RAN)


#Priority plots - overlap of sensitive and susceptible sites (predicted IAS range)
#sdm_col = colour of the predicted range, sdm_brd_col = colour of rectangle delineating range extent
PM_RAN_priority <- Priority_plot(PM_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#5bb4f2") #use sdm_brd_col = NA if you don't want a rectangle
VG_RAN_priority <- Priority_plot(VG_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#9c5200")
DG_RAN_priority <- Priority_plot(DG_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
TB_RAN_priority <- Priority_plot(TB_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
PL_RAN_priority <- Priority_plot(PL_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")
AM_RAN_priority <- Priority_plot(AM_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#5bb4f2") #use sdm_brd_col = NA if you don't want a rectangle
MF_RAN_priority <- Priority_plot(MF_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#9c5200")
MD_RAN_priority <- Priority_plot(MD_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
LH_RAN_priority <- Priority_plot(LH_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
VV_RAN_priority <- Priority_plot(VV_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")
#MR_RAN_priority <- Priority_plot(MR_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
BT_RAN_priority <- Priority_plot(BT_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
HA_RAN_priority <- Priority_plot(HA_bin_RAN, RAN_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")

                          ##CAZ with areas and weights##
#Pheidole megacephala
PM_bin <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# PM_bin2 <- PM_bin
# PM_bin2[PM_bin2 != 0] <- 1

# Re-class zeros as NA
PM_bin[PM_bin == 0] <- NA
PM_bin <- resample(PM_bin, CAZ_area_wgt_var_ras) #get equal extent
# PM_bin_p <- rasterToPolygons(PM_bin, dissolve = T)
# PM_bin_sf <- st_as_sf(PM_bin_p)
# PM_bin_sf_bb <- st_as_sfc(st_bbox(PM_bin_sf))


#Vespula germanica
VG_bin <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# VG_bin2 <- VG_bin
# VG_bin2[VG_bin2 != 0] <- 1

# Re-class zeros as NA
VG_bin[VG_bin == 0] <- NA
VG_bin <- resample(VG_bin, CAZ_area_wgt_var_ras) #get equal extent
# VG_bin_p <- rasterToPolygons(VG_bin, dissolve = T)
# VG_bin_sf <- st_as_sf(VG_bin_p)
# VG_bin_sf_bb <- st_as_sfc(st_bbox(VG_bin_sf))

#Digitonthophagus gazella
DG_bin <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# DG_bin2 <- DG_bin
# DG_bin2[DG_bin2 != 0] <- 1

# Re-class zeros as NA
DG_bin[DG_bin == 0] <- NA
DG_bin <- resample(DG_bin, CAZ_area_wgt_var_ras) #get equal extent
# DG_bin_p <- rasterToPolygons(DG_bin, dissolve = T)
# DG_bin_sf <- st_as_sf(DG_bin_p)
# DG_bin_sf_bb <- st_as_sfc(st_bbox(DG_bin_sf))

#Tetramorium bicarinatum
TB_bin <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# TB_bin2 <- TB_bin
# TB_bin2[TB_bin2 != 0] <- 1

# Re-class zeros as NA
TB_bin[TB_bin == 0] <- NA
TB_bin <- resample(TB_bin, CAZ_area_wgt_var_ras) #get equal extent
# TB_bin_p <- rasterToPolygons(TB_bin, dissolve = T)
# TB_bin_sf <- st_as_sf(TB_bin_p)
# TB_bin_sf_bb <- st_as_sfc(st_bbox(TB_bin_sf))

#Paratrechina longicornis
PL_bin <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# PL_bin2 <- PL_bin
# PL_bin2[PL_bin2 != 0] <- 1

# Re-class zeros as NA
PL_bin[PL_bin == 0] <- NA
PL_bin <- resample(PL_bin, CAZ_area_wgt_var_ras) #get equal extent
# PL_bin_p <- rasterToPolygons(PL_bin, dissolve = T)
# PL_bin_sf <- st_as_sf(PL_bin_p)
# PL_bin_sf_bb <- st_as_sfc(st_bbox(PL_bin_sf))

#Apis mellifera
AM_bin <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# AM_bin2 <- AM_bin
# AM_bin2[AM_bin2 != 0] <- 1

# Re-class zeros as NA
AM_bin[AM_bin == 0] <- NA
AM_bin <- resample(AM_bin, CAZ_area_wgt_var_ras) #get equal extent
# AM_bin_p <- rasterToPolygons(AM_bin, dissolve = T)
# AM_bin_sf <- st_as_sf(AM_bin_p)
# AM_bin_sf_bb <- st_as_sfc(st_bbox(AM_bin_sf))

#Monomorium floricola 
MF_bin <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# MF_bin2 <- MF_bin
# MF_bin2[MF_bin2 != 0] <- 1

# Re-class zeros as NA
MF_bin[PL_bin == 0] <- NA
MF_bin <- resample(MF_bin, CAZ_area_wgt_var_ras) #get equal extent
# MF_bin_p <- rasterToPolygons(MF_bin, dissolve = T)
# MF_bin_sf <- st_as_sf(MF_bin_p)
# MF_bin_sf_bb <- st_as_sfc(st_bbox(MF_bin_sf))

#Monomorium destructor
MD_bin <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# MD_bin2 <- MD_bin
# MD_bin2[MD_bin2 != 0] <- 1

# Re-class zeros as NA
MD_bin[PL_bin == 0] <- NA
MD_bin <- resample(MD_bin, CAZ_area_wgt_var_ras) #get equal extent
# MD_bin_p <- rasterToPolygons(MD_bin, dissolve = T)
# MD_bin_sf <- st_as_sf(MD_bin_p)
# MD_bin_sf_bb <- st_as_sfc(st_bbox(MD_bin_sf))

#Linepithema humile
LH_bin <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# LH_bin2 <- LH_bin
# LH_bin2[LH_bin2 != 0] <- 1

# Re-class zeros as NA
LH_bin[LH_bin == 0] <- NA
LH_bin <- resample(LH_bin, CAZ_area_wgt_var_ras) #get equal extent
# LH_bin_p <- rasterToPolygons(LH_bin, dissolve = T)
# LH_bin_sf <- st_as_sf(LH_bin_p)
# LH_bin_sf_bb <- st_as_sfc(st_bbox(LH_bin_sf))

#Vespula vulgaris
VV_bin <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# VV_bin2 <- VV_bin
# VV_bin2[VV_bin2 != 0] <- 1

# Re-class zeros as NA
VV_bin[PL_bin == 0] <- NA
VV_bin <- resample(VV_bin, CAZ_area_wgt_var_ras) #get equal extent
# VV_bin_p <- rasterToPolygons(VV_bin, dissolve = T)
# VV_bin_sf <- st_as_sf(VV_bin_p)
# VV_bin_sf_bb <- st_as_sfc(st_bbox(VV_bin_sf))

#Megachile rotundata 
# MR_bin <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# MR_bin2 <- MR_bin
# MR_bin2[MR_bin2 != 0] <- 1
# 
# # Re-class zeros as NA
# MR_bin[MR_bin == 0] <- NA
# MR_bin <- resample(MR_bin, CAZ_area_wgt_var_ras) #get equal extent
# MR_bin_p <- rasterToPolygons(MR_bin, dissolve = T)
# MR_bin_sf <- st_as_sf(MR_bin_p)
# MR_bin_sf_bb <- st_as_sfc(st_bbox(MR_bin_sf))

#Bombus terrestris
BT_bin <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# BT_bin2 <- BT_bin
# BT_bin2[BT_bin2 != 0] <- 1

# Re-class zeros as NA
BT_bin[BT_bin == 0] <- NA
BT_bin <- resample(BT_bin, CAZ_area_wgt_var_ras) #get equal extent
# BT_bin_p <- rasterToPolygons(BT_bin, dissolve = T)
# BT_bin_sf <- st_as_sf(BT_bin_p)
# BT_bin_sf_bb <- st_as_sfc(st_bbox(BT_bin_sf))

#Heteronychus arator
HA_bin <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# HA_bin2 <- HA_bin
# HA_bin2[HA_bin2 != 0] <- 1

# Re-class zeros as NA
HA_bin[HA_bin == 0] <- NA
HA_bin <- resample(HA_bin, CAZ_area_wgt_var_ras) #get equal extent
# HA_bin_p <- rasterToPolygons(HA_bin, dissolve = T)
# HA_bin_sf <- st_as_sf(HA_bin_p)
# HA_bin_sf_bb <- st_as_sfc(st_bbox(HA_bin_sf))


                                        ##RAN area with weights##
#Pheidole megacephala
PM_bin_RAN <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PM_bin_RAN[PM_bin_RAN == 0] <- NA
PM_bin_RAN <- resample(PM_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Vespula germanica
VG_bin_RAN <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VG_bin_RAN[VG_bin_RAN == 0] <- NA
VG_bin_RAN <- resample(VG_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Digitonthophagus gazella
DG_bin_RAN <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
DG_bin_RAN[DG_bin_RAN == 0] <- NA
DG_bin_RAN <- resample(DG_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Tetramorium bicarinatum
TB_bin_RAN <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
TB_bin_RAN[TB_bin_RAN == 0] <- NA
TB_bin_RAN <- resample(TB_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Paratrechina longicornis
PL_bin_RAN <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PL_bin_RAN[PL_bin_RAN == 0] <- NA
PL_bin_RAN <- resample(PL_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Apis mellifera
AM_bin_RAN <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
AM_bin_RAN[AM_bin_RAN == 0] <- NA
AM_bin_RAN <- resample(AM_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Monomorium floricola 
MF_bin_RAN <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MF_bin_RAN[MF_bin_RAN == 0] <- NA
MF_bin_RAN <- resample(MF_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Monomorium destructor
MD_bin_RAN <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MD_bin_RAN[MD_bin_RAN == 0] <- NA
MD_bin_RAN <- resample(MD_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Linepithema humile
LH_bin_RAN <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
LH_bin_RAN[LH_bin_RAN == 0] <- NA
LH_bin_RAN <- resample(LH_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Vespula vulgaris
VV_bin_RAN <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VV_bin_RAN[VV_bin_RAN == 0] <- NA
VV_bin_RAN <- resample(VV_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Megachile rotundata 
# MR_bin_RAN <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# # Re-class zeros as NA
# MR_bin_RAN[MR_bin_RAN == 0] <- NA
# MR_bin_RAN <- resample(MR_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Bombus terrestris
BT_bin_RAN <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
BT_bin_RAN[BT_bin_RAN == 0] <- NA
BT_bin_RAN <- resample(BT_bin_RAN, RAN_area_wgt_var_ras) #get equal extent

#Heteronychus arator
HA_bin_RAN <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
HA_bin_RAN[HA_bin_RAN == 0] <- NA
HA_bin_RAN <- resample(HA_bin_RAN, RAN_area_wgt_var_ras) #get equal extent


#Map for inset (species + area + weights)
CAZ_area_wgt_var_ras_st <- st_as_stars(CAZ_area_wgt_var_ras)
CAZ_area_wgt_var_ras_sf <- st_as_sf(CAZ_area_wgt_var_ras_st)
for_inset <- ggplot()+
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
        axis.ticks = element_blank(),
        legend.position = "right")
#to add rectangles for bbox of IAS sdm
map_w_borders <- for_inset + 
  geom_sf(data = PM_bin_sf_bb, fill = NA, color = "#5bb4f2") +
  geom_sf(data = VG_bin_sf_bb, fill = NA, color = "#9c5200") +
  geom_sf(data = DG_bin_sf_bb, fill = NA, color = "#b409a7") +
  geom_sf(data = TB_bin_sf_bb, fill = NA, color = "#040200") +
  geom_sf(data = PL_bin_sf_bb, fill = NA, color = "#0db02f")

#Histograms of IAS overlap with sensitive sites
PM_hist <- mask_hist(CAZ_area_wgt_var_ras, PM_bin)
VG_hist <- mask_hist(CAZ_area_wgt_var_ras, VG_bin)
DG_hist <- mask_hist(CAZ_area_wgt_var_ras, DG_bin)
TB_hist <- mask_hist(CAZ_area_wgt_var_ras, TB_bin)
PL_hist <- mask_hist(CAZ_area_wgt_var_ras, PL_bin)
AM_hist <- mask_hist(CAZ_area_wgt_var_ras, AM_bin)
MF_hist <- mask_hist(CAZ_area_wgt_var_ras, MF_bin)
MD_hist <- mask_hist(CAZ_area_wgt_var_ras, MD_bin)
LH_hist <- mask_hist(CAZ_area_wgt_var_ras, LH_bin)
VV_hist <- mask_hist(CAZ_area_wgt_var_ras, VV_bin)
#MR_hist <- mask_hist(CAZ_area_wgt_var_ras, MR_bin)
BT_hist <- mask_hist(CAZ_area_wgt_var_ras, BT_bin)
HA_hist <- mask_hist(CAZ_area_wgt_var_ras, HA_bin)


#Priority plots - overlap of sensitive and susceptible sites (predicted IAS range)
#sdm_col = colour of the predicted range, sdm_brd_col = colour of rectangle delineating range extent
PM_priority <- Priority_plot(PM_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#5bb4f2") #use sdm_brd_col = NA if you don't want a rectangle
VG_priority <- Priority_plot(VG_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#9c5200")
DG_priority <- Priority_plot(DG_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
TB_priority <- Priority_plot(TB_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
PL_priority <- Priority_plot(PL_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")
AM_priority <- Priority_plot(AM_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#5bb4f2") #use sdm_brd_col = NA if you don't want a rectangle
MF_priority <- Priority_plot(MF_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#9c5200")
MD_priority <- Priority_plot(MD_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
LH_priority <- Priority_plot(LH_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
VV_priority <- Priority_plot(VV_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")
#MR_priority <- Priority_plot(MR_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#b409a7")
BT_priority <- Priority_plot(BT_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#040200")
HA_priority <- Priority_plot(HA_bin, CAZ_area_wgt_var_ras, sdm_col = "black", sdm_brd_col = "#0db02f")


#Overlap among all IAS (turn into function?)
ias_sum <- sum(PM_bin2, VG_bin2, DG_bin2, TB_bin2, PL_bin2, AM_bin2, MF_bin2, MD_bin2, LH_bin2,
               VV_bin2, BT_bin2, HA_bin2)
ias_sum_one <- ias_sum
ias_sum_one[ias_sum_one != 1] <- NA


####################################### KBAs ######################################
#"Inverted" KBA overlap with IAS (KBA threat status 5 more important than 0)
#species + weights

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

#Apis mellifera
AM_bin_KBA <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
AM_bin_KBA[AM_bin_KBA == 0] <- NA
AM_bin_KBA <- resample(AM_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Monomorium floricola 
MF_bin_KBA <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MF_bin_KBA[MF_bin_KBA == 0] <- NA
MF_bin_KBA <- resample(MF_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Monomorium destructor
MD_bin_KBA <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MD_bin_KBA[MD_bin_KBA == 0] <- NA
MD_bin_KBA <- resample(MD_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Linepithema humile
LH_bin_KBA <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
LH_bin_KBA[LH_bin_KBA == 0] <- NA
LH_bin_KBA <- resample(LH_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Vespula vulgaris
VV_bin_KBA <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VV_bin_KBA[VV_bin_KBA == 0] <- NA
VV_bin_KBA <- resample(VV_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Megachile rotundata 
# MR_bin_KBA <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# # Re-class zeros as NA
# MR_bin_KBA[MR_bin_KBA == 0] <- NA
# MR_bin_KBA <- resample(MR_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Bombus terrestris
BT_bin_KBA <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
BT_bin_KBA[BT_bin_KBA == 0] <- NA
BT_bin_KBA <- resample(BT_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Heteronychus arator
HA_bin_KBA <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
HA_bin_KBA[HA_bin_KBA == 0] <- NA
HA_bin_KBA <- resample(HA_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Histograms of IAS overlap with sensitive sites
PM_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, PM_bin_KBA)
VG_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, VG_bin_KBA)
DG_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, DG_bin_KBA)
TB_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, TB_bin_KBA)
PL_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, PL_bin_KBA)
AM_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, AM_bin_KBA)
MF_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, MF_bin_KBA)
MD_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, MD_bin_KBA)
LH_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, LH_bin_KBA)
VV_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, VV_bin_KBA)
#MR_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, MR_bin_KBA)
BT_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, BT_bin_KBA)
HA_KBA_hist <- mask_hist(CAZ_wgt_KBA_inv_var_ras, HA_bin_KBA)

#Get cell values for Kolmogorov-smirnoff tests
PM_vals <- get_msk_vals(CAZ_wgt_var_ras, PM_bin)
PM_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, PM_bin_KBA)
PM_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, PM_bin_RAN)

VG_vals <- get_msk_vals(CAZ_wgt_var_ras, VG_bin)
VG_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, VG_bin_KBA)
VG_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, VG_bin_RAN)

DG_vals <- get_msk_vals(CAZ_wgt_var_ras, DG_bin)
DG_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, DG_bin_KBA)
DG_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, DG_bin_RAN)

TB_vals <- get_msk_vals(CAZ_wgt_var_ras, TB_bin)
TB_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, TB_bin_KBA)
TB_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, TB_bin_RAN)

PL_vals <- get_msk_vals(CAZ_wgt_var_ras, PL_bin)
PL_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, PL_bin_KBA)
PL_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, PL_bin_RAN)

AM_vals <- get_msk_vals(CAZ_wgt_var_ras, AM_bin)
AM_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, AM_bin_KBA)
AM_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, AM_bin_RAN)

MF_vals <- get_msk_vals(CAZ_wgt_var_ras, MF_bin)
MF_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, MF_bin_KBA)
MF_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, MF_bin_RAN)

MD_vals <- get_msk_vals(CAZ_wgt_var_ras, MD_bin)
MD_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, MD_bin_KBA)
MD_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, MD_bin_RAN)

LH_vals <- get_msk_vals(CAZ_wgt_var_ras, LH_bin)
LH_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, LH_bin_KBA)
LH_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, LH_bin_RAN)

VV_vals <- get_msk_vals(CAZ_wgt_var_ras, VV_bin)
VV_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, VV_bin_KBA)
VV_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, VV_bin_RAN)

# MR_vals <- get_msk_vals(CAZ_wgt_var_ras, MR_bin)
# MR_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, MR_bin_KBA)
# MR_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, MR_bin_RAN)

BT_vals <- get_msk_vals(CAZ_wgt_var_ras, BT_bin)
BT_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, BT_bin_KBA)
BT_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, BT_bin_RAN)

HA_vals <- get_msk_vals(CAZ_wgt_var_ras, HA_bin)
HA_KBA_vals <- get_msk_vals(CAZ_wgt_KBA_inv_var_ras, HA_bin_KBA)
HA_RAN_vals <- get_msk_vals(RAN_wgt_var_ras, HA_bin_RAN)

#Multi-panel figures
PM_comb <- Figure("Pheidole megacephala", CAZ_wgt_var_ras, PM_bin)
VG_comb <- Figure("Vespula germanica", CAZ_wgt_var_ras, PM_bin)
DG_comb <- Figure("Digitonthophagus gazella", CAZ_wgt_var_ras, DG_bin)
PL_comb <- Figure("Paratrechina longicornis", CAZ_wgt_var_ras, PL_bin)
TB_comb <- Figure("Tetramorium bicarinatum", CAZ_wgt_var_ras, TB_bin)
AM_comb <- Figure("Apis mellifera", CAZ_wgt_var_ras, AM_bin)
MF_comb <- Figure("Monomorium floricola", CAZ_wgt_var_ras, MF_bin)
MD_comb <- Figure("Monomorium destructor", CAZ_wgt_var_ras, MD_bin)
LH_comb <- Figure("Linepithema humile", CAZ_wgt_var_ras, LH_bin)
VV_comb <- Figure("Vespula vulgaris", CAZ_wgt_var_ras, VV_bin)
#MR_comb <- Figure("Megachile rotundata", CAZ_wgt_var_ras, MR_bin)
BT_comb <- Figure("Bombus terrestris", CAZ_wgt_var_ras, BT_bin)
HA_comb <- Figure("Heteronychus arator", CAZ_wgt_var_ras, HA_bin)

#Kolmogorov-smirnoff tests - KBA vs no KBA
#a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution
#presence of ties a result of rounding
ks.test(PM_KBA_vals, PM_vals, alternative = "two.sided")
ks.test(VG_KBA_vals, VG_vals, alternative = "two.sided")
ks.test(DG_KBA_vals, DG_vals, alternative = "two.sided")
ks.test(TB_KBA_vals, TB_vals, alternative = "two.sided")
ks.test(PL_KBA_vals, PL_vals, alternative = "two.sided")
ks.test(AM_KBA_vals, AM_vals, alternative = "two.sided")
ks.test(MF_KBA_vals, MF_vals, alternative = "two.sided")
ks.test(MD_KBA_vals, MD_vals, alternative = "two.sided")
ks.test(LH_KBA_vals, LH_vals, alternative = "two.sided")
ks.test(VV_KBA_vals, VV_vals, alternative = "two.sided")
#ks.test(MR_KBA_vals, MR_vals, alternative = "two.sided")
ks.test(BT_KBA_vals, BT_vals, alternative = "two.sided")
ks.test(HA_KBA_vals, HA_vals, alternative = "two.sided")


#Difference between species + weight IAS
PM_VG <- ks.test(PM_vals, VG_vals, alternative = "two.sided")
PM_DG <- ks.test(PM_vals, DG_vals, alternative = "two.sided")
PM_TB <- ks.test(PM_vals, TB_vals, alternative = "two.sided")
PM_PL <- ks.test(PM_vals, PL_vals, alternative = "two.sided")
PM_AM <- ks.test(PM_vals, AM_vals, alternative = "two.sided")
PM_MF <- ks.test(PM_vals, MF_vals, alternative = "two.sided")
PM_MD <- ks.test(PM_vals, MD_vals, alternative = "two.sided")
PM_LH <- ks.test(PM_vals, LH_vals, alternative = "two.sided")
PM_VV <- ks.test(PM_vals, VV_vals, alternative = "two.sided")
#PM_MR <- ks.test(PM_vals, MR_vals, alternative = "two.sided")
PM_BT <- ks.test(PM_vals, BT_vals, alternative = "two.sided")
PM_HA <- ks.test(PM_vals, HA_vals, alternative = "two.sided")

VG_DG <- ks.test(VG_vals, DG_vals, alternative = "two.sided")
VG_TB <- ks.test(VG_vals, TB_vals, alternative = "two.sided")
VG_PL <- ks.test(VG_vals, PL_vals, alternative = "two.sided")
VG_AM <- ks.test(VG_vals, AM_vals, alternative = "two.sided")
VG_MF <- ks.test(VG_vals, MF_vals, alternative = "two.sided")
VG_MD <- ks.test(VG_vals, MD_vals, alternative = "two.sided")
VG_LH <- ks.test(VG_vals, LH_vals, alternative = "two.sided")
VG_VV <- ks.test(VG_vals, VV_vals, alternative = "two.sided")
#VG_MR <- ks.test(VG_vals, MR_vals, alternative = "two.sided")
VG_BT <- ks.test(VG_vals, BT_vals, alternative = "two.sided")
VG_HA <- ks.test(VG_vals, HA_vals, alternative = "two.sided")

DG_TB <- ks.test(DG_vals, TB_vals, alternative = "two.sided")
DG_PL <- ks.test(DG_vals, PL_vals, alternative = "two.sided")
DG_AM <- ks.test(DG_vals, AM_vals, alternative = "two.sided")
DG_MF <- ks.test(DG_vals, MF_vals, alternative = "two.sided")
DG_MD <- ks.test(DG_vals, MD_vals, alternative = "two.sided")
DG_LH <- ks.test(DG_vals, LH_vals, alternative = "two.sided")
DG_VV <- ks.test(DG_vals, VV_vals, alternative = "two.sided")
#DG_MR <- ks.test(DG_vals, MR_vals, alternative = "two.sided")
DG_BT <- ks.test(DG_vals, BT_vals, alternative = "two.sided")
DG_HA <- ks.test(DG_vals, HA_vals, alternative = "two.sided")

TB_PL <- ks.test(TB_vals, PL_vals, alternative = "two.sided")
TB_AM <- ks.test(TB_vals, AM_vals, alternative = "two.sided")
TB_MF <- ks.test(TB_vals, MF_vals, alternative = "two.sided")
TB_MD <- ks.test(TB_vals, MD_vals, alternative = "two.sided")
TB_LH <- ks.test(TB_vals, LH_vals, alternative = "two.sided")
TB_VV <- ks.test(TB_vals, VV_vals, alternative = "two.sided")
#TB_MR <- ks.test(TB_vals, MR_vals, alternative = "two.sided")
TB_BT <- ks.test(TB_vals, BT_vals, alternative = "two.sided")
TB_HA <- ks.test(TB_vals, HA_vals, alternative = "two.sided")

PL_AM <- ks.test(PL_vals, AM_vals, alternative = "two.sided")
PL_MF <- ks.test(PL_vals, MF_vals, alternative = "two.sided")
PL_MD <- ks.test(PL_vals, MD_vals, alternative = "two.sided")
PL_LH <- ks.test(PL_vals, LH_vals, alternative = "two.sided")
PL_VV <- ks.test(PL_vals, VV_vals, alternative = "two.sided")
#PL_MR <- ks.test(PL_vals, MR_vals, alternative = "two.sided")
PL_BT <- ks.test(PL_vals, BT_vals, alternative = "two.sided")
PL_HA <- ks.test(PL_vals, HA_vals, alternative = "two.sided")

AM_MF <- ks.test(AM_vals, MF_vals, alternative = "two.sided")
AM_MD <- ks.test(AM_vals, MD_vals, alternative = "two.sided")
AM_LH <- ks.test(AM_vals, LH_vals, alternative = "two.sided")
AM_VV <- ks.test(AM_vals, VV_vals, alternative = "two.sided")
#AM_MR <- ks.test(AM_vals, MR_vals, alternative = "two.sided")
AM_BT <- ks.test(AM_vals, BT_vals, alternative = "two.sided")
AM_HA <- ks.test(AM_vals, HA_vals, alternative = "two.sided")

MF_MD <- ks.test(MF_vals, MD_vals, alternative = "two.sided")
MF_LH <- ks.test(MF_vals, LH_vals, alternative = "two.sided")
MF_VV <- ks.test(MF_vals, VV_vals, alternative = "two.sided")
#MF_MR <- ks.test(MF_vals, MR_vals, alternative = "two.sided")
MF_BT <- ks.test(MF_vals, BT_vals, alternative = "two.sided")
MF_HA <- ks.test(MF_vals, HA_vals, alternative = "two.sided")

MD_LH <- ks.test(MD_vals, LH_vals, alternative = "two.sided")
MD_VV <- ks.test(MD_vals, VV_vals, alternative = "two.sided")
#MD_MR <- ks.test(MD_vals, MR_vals, alternative = "two.sided")
MD_BT <- ks.test(MD_vals, BT_vals, alternative = "two.sided")
MD_HA <- ks.test(MD_vals, HA_vals, alternative = "two.sided")

LH_VV <- ks.test(LH_vals, VV_vals, alternative = "two.sided")
#LH_MR <- ks.test(LH_vals, MR_vals, alternative = "two.sided")
LH_BT <- ks.test(LH_vals, BT_vals, alternative = "two.sided")
LH_HA <- ks.test(LH_vals, HA_vals, alternative = "two.sided")

#VV_MR <- ks.test(VV_vals, MR_vals, alternative = "two.sided")
VV_BT <- ks.test(VV_vals, BT_vals, alternative = "two.sided")
VV_HA <- ks.test(VV_vals, HA_vals, alternative = "two.sided")

#MR_BT <- ks.test(MR_vals, BT_vals, alternative = "two.sided")
#MR_HA <- ks.test(MR_vals, HA_vals, alternative = "two.sided")

BT_HA <- ks.test(BT_vals, HA_vals, alternative = "two.sided")

v <- c(0,PM_VG$statistic,PM_DG$statistic,PM_TB$statistic,PM_PL$statistic,PM_AM$statistic,PM_MF$statistic,PM_MD$statistic,PM_LH$statistic,PM_VV$statistic,PM_BT$statistic,PM_HA$statistic,
       0,0,VG_DG$statistic,VG_TB$statistic,VG_PL$statistic,VG_AM$statistic,VG_MF$statistic,VG_MD$statistic,VG_LH$statistic,VG_VV$statistic,VG_BT$statistic,VG_HA$statistic,
       0,0,0,DG_TB$statistic,DG_PL$statistic,DG_AM$statistic,DG_MF$statistic,DG_MD$statistic,DG_LH$statistic,DG_VV$statistic,DG_BT$statistic,DG_HA$statistic,
       0,0,0,0,TB_PL$statistic,TB_AM$statistic,TB_MF$statistic,TB_MD$statistic,TB_LH$statistic,TB_VV$statistic,TB_BT$statistic,TB_HA$statistic,
       0,0,0,0,0,AM_MF$statistic,AM_MD$statistic,AM_LH$statistic,AM_VV$statistic,AM_BT$statistic,AM_HA$statistic,
       0,0,0,0,0,0,MF_MD$statistic,MF_LH$statistic,MF_VV$statistic,MF_BT$statistic,MF_HA$statistic,
       0,0,0,0,0,0,0,MD_LH$statistic,MD_VV$statistic,MD_BT$statistic,MD_HA$statistic,
       0,0,0,0,0,0,0,0,LH_VV$statistic,LH_BT$statistic,LH_HA$statistic,
       0,0,0,0,0,0,0,0,0,VV_BT$statistic,VV_HA$statistic,
       0,0,0,0,0,0,0,0,0,0,BT_HA$statistic)
       
tm <- matrix(v, nrow = 12, ncol = 11)
rownames(tm) <- c(":italic(Pheidole~~megacephala)", 
                  ":italic(Vespula~~germanica)", 
                  ":italic(Digitonthophagus~~gazella)", 
                  ":italic(Tetramorium~~bicarinatum)", 
                  ":italic(Paratrechina~~longicornis)",
                  ":italic(Apis~~mellifera)",
                  ":italic(Monomorium~~floricola)",
                  ":italic(Monomorium~~destructor)",
                  ":italic(Linepithema~~humile)",
                  ":italic(Vespula~~vulgaris)",
                  ":italic(Bombus~~terrestris)",
                  ":italic(Heteronychus~~arator)")
colnames(tm) <- c("Pm", "Vg", "Dg", "Tb", "Pl", "Am", "Mf", "Md", "Lh", "Vv", "Bt", "Ha")
corrplot::corrplot(tm, type = "lower", method = "color", 
                   cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
                   tl.srt = 0, tl.col = "black", tl.cex = 0.8, addCoef.col = "black",
                   mar = c(0,0,0,0))

#Proportion difference (KBA vs no KBA vs Random) in number of top sensitive sites
#Top two (i.e. >= 0.98 sensitivity)
props <- c(1, 0.98, 0.95, 0.90, 0.75, 0.50, 0.25, 0.00)
diffs <- c()

PM_prop <- multi_props(PM_vals, props)
PM_KBA_prop <- multi_props(PM_KBA_vals, props)
PM_RAN_prop <- multi_props(PM_RAN_vals, props)

VG_prop <- multi_props(VG_vals, props)
VG_KBA_prop <- multi_props(VG_KBA_vals, props)
VG_KBA_prop <- multi_props(VG_RAN_vals, props)

DG_prop <- multi_props(DG_vals, props)
DG_KBA_prop <- multi_props(DG_KBA_vals, props)
DG_KBA_prop <- multi_props(DG_RAN_vals, props)

TB_prop <- multi_props(TB_vals, props)
TB_KBA_prop <- multi_props(TB_KBA_vals, props)
TB_KBA_prop <- multi_props(TB_RAN_vals, props)

PL_prop <- multi_props(PL_vals, props)
PL_KBA_prop <- multi_props(PL_KBA_vals, props)
PL_RAN_prop <- multi_props(PL_RAN_vals, props)

AM_prop <- multi_props(AM_vals, props)
AM_KBA_prop <- multi_props(AM_KBA_vals, props)
AM_RAN_prop <- multi_props(AM_RAN_vals, props)

MF_prop <- multi_props(MF_vals, props)
MF_KBA_prop <- multi_props(MF_KBA_vals, props)
MF_RAN_prop <- multi_props(MF_RAN_vals, props)

MD_prop <- multi_props(MD_vals, props)
MD_KBA_prop <- multi_props(MD_KBA_vals, props)
MD_RAN_prop <- multi_props(MD_RAN_vals, props)

LH_prop <- multi_props(LH_vals, props)
LH_KBA_prop <- multi_props(LH_KBA_vals, props)
LH_RAN_prop <- multi_props(LH_RAN_vals, props)

VV_prop <- multi_props(VV_vals, props)
VV_KBA_prop <- multi_props(VV_KBA_vals, props)
VV_RAN_prop <- multi_props(VV_RAN_vals, props)

# MR_prop <- multi_props(MR_vals, props)
# MR_KBA_prop <- multi_props(MR_KBA_vals, props)
# MR_RAN_prop <- multi_props(MR_RAN_vals, props)

BT_prop <- multi_props(BT_vals, props)
BT_KBA_prop <- multi_props(BT_KBA_vals, props)
BT_RAN_prop <- multi_props(BT_RAN_vals, props)

HA_prop <- multi_props(HA_vals, props)
HA_KBA_prop <- multi_props(HA_KBA_vals, props)
HA_RAN_prop <- multi_props(HA_RAN_vals, props)

Total <- as.data.frame(rbind(PM_prop,VG_prop,DG_prop,TB_prop,PL_prop,AM_prop,MF_prop,MD_prop,LH_prop,VV_prop,BT_prop,HA_prop,
           PM_KBA_prop, VG_KBA_prop, DG_KBA_prop, TB_KBA_prop, PL_KBA_prop,AM_KBA_prop,MF_KBA_prop,MD_KBA_prop,LH_KBA_prop,VV_KBA_prop,BT_KBA_prop,HA_KBA_prop,
           PM_RAN_prop, VG_RAN_prop, DG_RAN_prop, TB_RAN_prop, PL_RAN_prop,AM_RAN_prop,MF_RAN_prop,MD_RAN_prop,LH_RAN_prop,VV_RAN_prop,BT_RAN_prop,HA_RAN_prop))
Type <- c(rep("Unmasked", 12), rep("Masked", 12), rep("Random", 12))
Species <- rep(c("P. megacephala", "V. germanica", "D. gazella",
                 "T. bicarinatum", "P. longicornis", "A. mellifera",
                 "M. floricola", "M. destructor", "L. humile",
                 "V. vulgaris", "B. terrestris", "H. arator"), 3)
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
Total_min <- Total %>%
  dplyr::filter(SiteSensitivity == 1 | SiteSensitivity == 0.98)
df <- data.frame(xmin = c(1,0.98,0.95,0.9,0.75,0.5,0.25),
                 xmax = c(0.98,0.95,0.9,0.75,0.5,0.25,0),
                 ymin = rep(-Inf,7),
                 ymax = rep(Inf,7))
p <- ggplot(Total, 
       aes(x = SiteSensitivity, 
           y = DistributionCoverage, 
           colour = Species,
           group = interaction(Species,Type))) +
  geom_line(stat = "identity", aes(linetype = Type)) +
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
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1,0),
                     trans = "reverse") +
  scale_y_continuous(expand = c(0,0))
p2 <-  p + 
    geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin=ymin, ymax=ymax),
              fill = rev(leg$colors), alpha = 0.2, inherit.aes = F) +
    font("legend.text", face = "italic")

#For inset
ins <- ggplot(Total_min, 
         aes(x = SiteSensitivity, 
             y = DistributionCoverage, 
             colour = Species,
             group = interaction(Species,Type))) +
    geom_line(stat = "identity", aes(linetype = Type)) +
    scale_colour_manual(values = cols) +
    ylab("Distribution Coverage (%)") +
    xlab("Site Sensitivity") +
    theme_bw() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = alpha(leg$colors[7],0.2)),
          axis.text = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank()) +
    scale_x_continuous(expand = c(0,0),
                       limits = c(1,0.98),
                       trans = "reverse") +
    scale_y_continuous(expand = c(0,0))

####################################### KBAs ######################################
#"Inverted" KBA overlap with IAS (KBA threat status 5 more important than 0)
#species + area + weight

#Pheidole megacephala
PM_bin_KBA <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PM_bin_KBA[PM_bin_KBA == 0] <- NA
PM_bin_KBA <- resample(PM_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Vespula germanica
VG_bin_KBA <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VG_bin_KBA[VG_bin_KBA == 0] <- NA
VG_bin_KBA <- resample(VG_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Digitonthophagus gazella
DG_bin_KBA <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
DG_bin_KBA[DG_bin_KBA == 0] <- NA
DG_bin_KBA <- resample(DG_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Tetramorium bicarinatum
TB_bin_KBA <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
TB_bin_KBA[TB_bin_KBA == 0] <- NA
TB_bin_KBA <- resample(TB_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Paratrechina longicornis
PL_bin_KBA <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
PL_bin_KBA[PL_bin_KBA == 0] <- NA
PL_bin_KBA <- resample(PL_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Apis mellifera
AM_bin_KBA <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
AM_bin_KBA[AM_bin_KBA == 0] <- NA
AM_bin_KBA <- resample(AM_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Monomorium floricola 
MF_bin_KBA <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MF_bin_KBA[MF_bin_KBA == 0] <- NA
MF_bin_KBA <- resample(MF_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Monomorium destructor
MD_bin_KBA <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MD_bin_KBA[MD_bin_KBA == 0] <- NA
MD_bin_KBA <- resample(MD_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Linepithema humile
LH_bin_KBA <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
LH_bin_KBA[LH_bin_KBA == 0] <- NA
LH_bin_KBA <- resample(LH_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Vespula vulgaris
VV_bin_KBA <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
VV_bin_KBA[VV_bin_KBA == 0] <- NA
VV_bin_KBA <- resample(VV_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Megachile rotundata 
# MR_bin_KBA <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# # Re-class zeros as NA
# MR_bin_KBA[MR_bin_KBA == 0] <- NA
# MR_bin_KBA <- resample(MR_bin_KBA, CAZ_wgt_KBA_inv_var_ras) #get equal extent

#Bombus terrestris
BT_bin_KBA <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
BT_bin_KBA[BT_bin_KBA == 0] <- NA
BT_bin_KBA <- resample(BT_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Heteronychus arator
HA_bin_KBA <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
HA_bin_KBA[HA_bin_KBA == 0] <- NA
HA_bin_KBA <- resample(HA_bin_KBA, CAZ_area_wgt_KBA_inv_var_ras) #get equal extent

#Histograms of IAS overlap with sensitive sites
PM_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, PM_bin_KBA)
VG_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, VG_bin_KBA)
DG_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, DG_bin_KBA)
TB_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, TB_bin_KBA)
PL_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, PL_bin_KBA)
AM_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, AM_bin_KBA)
MF_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, MF_bin_KBA)
MD_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, MD_bin_KBA)
LH_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, LH_bin_KBA)
VV_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, VV_bin_KBA)
#MR_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, MR_bin_KBA)
BT_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, BT_bin_KBA)
HA_KBA_hist <- mask_hist(CAZ_area_wgt_KBA_inv_var_ras, HA_bin_KBA)

#Get cell values for Kolmogorov-smirnoff tests
PM_vals <- get_msk_vals(CAZ_area_wgt_var_ras, PM_bin)
PM_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, PM_bin_KBA)
PM_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, PM_bin_RAN)

VG_vals <- get_msk_vals(CAZ_area_wgt_var_ras, VG_bin)
VG_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, VG_bin_KBA)
VG_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, VG_bin_RAN)

DG_vals <- get_msk_vals(CAZ_area_wgt_var_ras, DG_bin)
DG_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, DG_bin_KBA)
DG_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, DG_bin_RAN)

TB_vals <- get_msk_vals(CAZ_area_wgt_var_ras, TB_bin)
TB_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, TB_bin_KBA)
TB_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, TB_bin_RAN)

PL_vals <- get_msk_vals(CAZ_area_wgt_var_ras, PL_bin)
PL_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, PL_bin_KBA)
PL_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, PL_bin_RAN)

AM_vals <- get_msk_vals(CAZ_area_wgt_var_ras, AM_bin)
AM_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, AM_bin_KBA)
AM_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, AM_bin_RAN)

MF_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MF_bin)
MF_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MF_bin_KBA)
MF_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, MF_bin_RAN)

MD_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MD_bin)
MD_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MD_bin_KBA)
MD_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, MD_bin_RAN)

LH_vals <- get_msk_vals(CAZ_area_wgt_var_ras, LH_bin)
LH_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, LH_bin_KBA)
LH_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, LH_bin_RAN)

VV_vals <- get_msk_vals(CAZ_area_wgt_var_ras, VV_bin)
VV_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, VV_bin_KBA)
VV_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, VV_bin_RAN)

# MR_vals <- get_msk_vals(CAZ_area_wgt_var_ras, MR_bin)
# MR_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, MR_bin_KBA)
# MR_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, MR_bin_RAN)

BT_vals <- get_msk_vals(CAZ_area_wgt_var_ras, BT_bin)
BT_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, BT_bin_KBA)
BT_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, BT_bin_RAN)

HA_vals <- get_msk_vals(CAZ_area_wgt_var_ras, HA_bin)
HA_KBA_vals <- get_msk_vals(CAZ_area_wgt_KBA_inv_var_ras, HA_bin_KBA)
HA_RAN_vals <- get_msk_vals(RAN_area_wgt_var_ras, HA_bin_RAN)

#Multi-panel figures
PM_comb <- Figure("Pheidole megacephala", CAZ_area_wgt_var_ras, PM_bin)
VG_comb <- Figure("Vespula germanica", CAZ_area_wgt_var_ras, VG_bin)
DG_comb <- Figure("Digitonthophagus gazella", CAZ_area_wgt_var_ras, DG_bin)
PL_comb <- Figure("Paratrechina longicornis", CAZ_area_wgt_var_ras, PL_bin)
TB_comb <- Figure("Tetramorium bicarinatum", CAZ_area_wgt_var_ras, TB_bin)
AM_comb <- Figure("Apis mellifera", CAZ_area_wgt_var_ras, AM_bin)
MF_comb <- Figure("Monomorium floricola", CAZ_area_wgt_var_ras, MF_bin)
MD_comb <- Figure("Monomorium destructor", CAZ_area_wgt_var_ras, MD_bin)
LH_comb <- Figure("Linepithema humile", CAZ_area_wgt_var_ras, LH_bin)
VV_comb <- Figure("Vespula vulgaris", CAZ_area_wgt_var_ras, VV_bin)
#MR_comb <- Figure("Megachile rotundata", CAZ_area_wgt_var_ras, MR_bin)
BT_comb <- Figure("Bombus terrestris", CAZ_area_wgt_var_ras, BT_bin)
HA_comb <- Figure("Heteronychus arator", CAZ_area_wgt_var_ras, HA_bin)

#Kolmogorov-smirnoff tests - KBA vs no KBA
#a two-sample test of the null hypothesis that x and y were drawn from the same continuous distribution
#presence of ties a result of rounding
ks.test(PM_KBA_vals, PM_vals, alternative = "two.sided")
ks.test(VG_KBA_vals, VG_vals, alternative = "two.sided")
ks.test(DG_KBA_vals, DG_vals, alternative = "two.sided")
ks.test(TB_KBA_vals, TB_vals, alternative = "two.sided")
ks.test(PL_KBA_vals, PL_vals, alternative = "two.sided")
ks.test(AM_KBA_vals, AM_vals, alternative = "two.sided")
ks.test(MF_KBA_vals, MF_vals, alternative = "two.sided")
ks.test(MD_KBA_vals, MD_vals, alternative = "two.sided")
ks.test(LH_KBA_vals, LH_vals, alternative = "two.sided")
ks.test(VV_KBA_vals, VV_vals, alternative = "two.sided")
#ks.test(MR_KBA_vals, MR_vals, alternative = "two.sided")
ks.test(BT_KBA_vals, BT_vals, alternative = "two.sided")
ks.test(HA_KBA_vals, HA_vals, alternative = "two.sided")


#Difference between species + area + weight IAS
PM_VG <- ks.test(PM_vals, VG_vals, alternative = "two.sided")
PM_DG <- ks.test(PM_vals, DG_vals, alternative = "two.sided")
PM_TB <- ks.test(PM_vals, TB_vals, alternative = "two.sided")
PM_PL <- ks.test(PM_vals, PL_vals, alternative = "two.sided")
PM_AM <- ks.test(PM_vals, AM_vals, alternative = "two.sided")
PM_MF <- ks.test(PM_vals, MF_vals, alternative = "two.sided")
PM_MD <- ks.test(PM_vals, MD_vals, alternative = "two.sided")
PM_LH <- ks.test(PM_vals, LH_vals, alternative = "two.sided")
PM_VV <- ks.test(PM_vals, VV_vals, alternative = "two.sided")
#PM_MR <- ks.test(PM_vals, MR_vals, alternative = "two.sided")
PM_BT <- ks.test(PM_vals, BT_vals, alternative = "two.sided")
PM_HA <- ks.test(PM_vals, HA_vals, alternative = "two.sided")

VG_DG <- ks.test(VG_vals, DG_vals, alternative = "two.sided")
VG_TB <- ks.test(VG_vals, TB_vals, alternative = "two.sided")
VG_PL <- ks.test(VG_vals, PL_vals, alternative = "two.sided")
VG_AM <- ks.test(VG_vals, AM_vals, alternative = "two.sided")
VG_MF <- ks.test(VG_vals, MF_vals, alternative = "two.sided")
VG_MD <- ks.test(VG_vals, MD_vals, alternative = "two.sided")
VG_LH <- ks.test(VG_vals, LH_vals, alternative = "two.sided")
VG_VV <- ks.test(VG_vals, VV_vals, alternative = "two.sided")
#VG_MR <- ks.test(VG_vals, MR_vals, alternative = "two.sided")
VG_BT <- ks.test(VG_vals, BT_vals, alternative = "two.sided")
VG_HA <- ks.test(VG_vals, HA_vals, alternative = "two.sided")

DG_TB <- ks.test(DG_vals, TB_vals, alternative = "two.sided")
DG_PL <- ks.test(DG_vals, PL_vals, alternative = "two.sided")
DG_AM <- ks.test(DG_vals, AM_vals, alternative = "two.sided")
DG_MF <- ks.test(DG_vals, MF_vals, alternative = "two.sided")
DG_MD <- ks.test(DG_vals, MD_vals, alternative = "two.sided")
DG_LH <- ks.test(DG_vals, LH_vals, alternative = "two.sided")
DG_VV <- ks.test(DG_vals, VV_vals, alternative = "two.sided")
#DG_MR <- ks.test(DG_vals, MR_vals, alternative = "two.sided")
DG_BT <- ks.test(DG_vals, BT_vals, alternative = "two.sided")
DG_HA <- ks.test(DG_vals, HA_vals, alternative = "two.sided")

TB_PL <- ks.test(TB_vals, PL_vals, alternative = "two.sided")
TB_AM <- ks.test(TB_vals, AM_vals, alternative = "two.sided")
TB_MF <- ks.test(TB_vals, MF_vals, alternative = "two.sided")
TB_MD <- ks.test(TB_vals, MD_vals, alternative = "two.sided")
TB_LH <- ks.test(TB_vals, LH_vals, alternative = "two.sided")
TB_VV <- ks.test(TB_vals, VV_vals, alternative = "two.sided")
#TB_MR <- ks.test(TB_vals, MR_vals, alternative = "two.sided")
TB_BT <- ks.test(TB_vals, BT_vals, alternative = "two.sided")
TB_HA <- ks.test(TB_vals, HA_vals, alternative = "two.sided")

PL_AM <- ks.test(PL_vals, AM_vals, alternative = "two.sided")
PL_MF <- ks.test(PL_vals, MF_vals, alternative = "two.sided")
PL_MD <- ks.test(PL_vals, MD_vals, alternative = "two.sided")
PL_LH <- ks.test(PL_vals, LH_vals, alternative = "two.sided")
PL_VV <- ks.test(PL_vals, VV_vals, alternative = "two.sided")
#PL_MR <- ks.test(PL_vals, MR_vals, alternative = "two.sided")
PL_BT <- ks.test(PL_vals, BT_vals, alternative = "two.sided")
PL_HA <- ks.test(PL_vals, HA_vals, alternative = "two.sided")

AM_MF <- ks.test(AM_vals, MF_vals, alternative = "two.sided")
AM_MD <- ks.test(AM_vals, MD_vals, alternative = "two.sided")
AM_LH <- ks.test(AM_vals, LH_vals, alternative = "two.sided")
AM_VV <- ks.test(AM_vals, VV_vals, alternative = "two.sided")
#AM_MR <- ks.test(AM_vals, MR_vals, alternative = "two.sided")
AM_BT <- ks.test(AM_vals, BT_vals, alternative = "two.sided")
AM_HA <- ks.test(AM_vals, HA_vals, alternative = "two.sided")

MF_MD <- ks.test(MF_vals, MD_vals, alternative = "two.sided")
MF_LH <- ks.test(MF_vals, LH_vals, alternative = "two.sided")
MF_VV <- ks.test(MF_vals, VV_vals, alternative = "two.sided")
#MF_MR <- ks.test(MF_vals, MR_vals, alternative = "two.sided")
MF_BT <- ks.test(MF_vals, BT_vals, alternative = "two.sided")
MF_HA <- ks.test(MF_vals, HA_vals, alternative = "two.sided")

MD_LH <- ks.test(MD_vals, LH_vals, alternative = "two.sided")
MD_VV <- ks.test(MD_vals, VV_vals, alternative = "two.sided")
#MD_MR <- ks.test(MD_vals, MR_vals, alternative = "two.sided")
MD_BT <- ks.test(MD_vals, BT_vals, alternative = "two.sided")
MD_HA <- ks.test(MD_vals, HA_vals, alternative = "two.sided")

LH_VV <- ks.test(LH_vals, VV_vals, alternative = "two.sided")
#LH_MR <- ks.test(LH_vals, MR_vals, alternative = "two.sided")
LH_BT <- ks.test(LH_vals, BT_vals, alternative = "two.sided")
LH_HA <- ks.test(LH_vals, HA_vals, alternative = "two.sided")

#VV_MR <- ks.test(VV_vals, MR_vals, alternative = "two.sided")
VV_BT <- ks.test(VV_vals, BT_vals, alternative = "two.sided")
VV_HA <- ks.test(VV_vals, HA_vals, alternative = "two.sided")

#MR_BT <- ks.test(MR_vals, BT_vals, alternative = "two.sided")
#MR_HA <- ks.test(MR_vals, HA_vals, alternative = "two.sided")

BT_HA <- ks.test(BT_vals, HA_vals, alternative = "two.sided")

v <- c(0,PM_VG$statistic,PM_DG$statistic,PM_TB$statistic,PM_PL$statistic,PM_AM$statistic,PM_MF$statistic,PM_MD$statistic,PM_LH$statistic,PM_VV$statistic,PM_BT$statistic,PM_HA$statistic,
       0,0,VG_DG$statistic,VG_TB$statistic,VG_PL$statistic,VG_AM$statistic,VG_MF$statistic,VG_MD$statistic,VG_LH$statistic,VG_VV$statistic,VG_BT$statistic,VG_HA$statistic,
       0,0,0,DG_TB$statistic,DG_PL$statistic,DG_AM$statistic,DG_MF$statistic,DG_MD$statistic,DG_LH$statistic,DG_VV$statistic,DG_BT$statistic,DG_HA$statistic,
       0,0,0,0,TB_PL$statistic,TB_AM$statistic,TB_MF$statistic,TB_MD$statistic,TB_LH$statistic,TB_VV$statistic,TB_BT$statistic,TB_HA$statistic,
       0,0,0,0,0,PL_AM$statistic,PL_MF$statistic,PL_MD$statistic,PL_LH$statistic,PL_VV$statistic,PL_BT$statistic,PL_HA$statistic,
       0,0,0,0,0,0,AM_MF$statistic,AM_MD$statistic,AM_LH$statistic,AM_VV$statistic,AM_BT$statistic,AM_HA$statistic,
       0,0,0,0,0,0,0,MF_MD$statistic,MF_LH$statistic,MF_VV$statistic,MF_BT$statistic,MF_HA$statistic,
       0,0,0,0,0,0,0,0,MD_LH$statistic,MD_VV$statistic,MD_BT$statistic,MD_HA$statistic,
       0,0,0,0,0,0,0,0,0,LH_VV$statistic,LH_BT$statistic,LH_HA$statistic,
       0,0,0,0,0,0,0,0,0,0,VV_BT$statistic,VV_HA$statistic,
       0,0,0,0,0,0,0,0,0,0,0,BT_HA$statistic)

tm <- matrix(v, nrow = 12, ncol = 11)
rownames(tm) <- c(":italic(Pheidole~~megacephala)", 
                  ":italic(Vespula~~germanica)", 
                  ":italic(Digitonthophagus~~gazella)", 
                  ":italic(Tetramorium~~bicarinatum)", 
                  ":italic(Paratrechina~~longicornis)",
                  ":italic(Apis~~mellifera)",
                  ":italic(Monomorium~~floricola)",
                  ":italic(Monomorium~~destructor)",
                  ":italic(Linepithema~~humile)",
                  ":italic(Vespula~~vulgaris)",
                  ":italic(Bombus~~terrestris)",
                  ":italic(Heteronychus~~arator)")
colnames(tm) <- c("Pm", "Vg", "Dg", "Tb", "Pl", "Am", "Mf", "Md", "Lh", "Vv", "Bt")
corrplot::corrplot(tm, type = "lower", method = "color", 
                   cl.pos = "n", col=brewer.pal(n=10, name="Spectral"), 
                   tl.srt = 0, tl.col = "black", tl.cex = 0.8, addCoef.col = "black",
                   mar = c(0,0,0,0))

#Proportion difference (KBA vs no KBA) in number of top sensitive sites
#Top two (i.e. >= 0.98 sensitivity)
props <- c(1, 0.98, 0.95, 0.90, 0.75, 0.50, 0.25, 0.00)
diffs <- c()

PM_prop <- multi_props(PM_vals, props)
PM_KBA_prop <- multi_props(PM_KBA_vals, props)
PM_RAN_prop <- multi_props(PM_RAN_vals, props)

VG_prop <- multi_props(VG_vals, props)
VG_KBA_prop <- multi_props(VG_KBA_vals, props)
VG_RAN_prop <- multi_props(VG_RAN_vals, props)

DG_prop <- multi_props(DG_vals, props)
DG_KBA_prop <- multi_props(DG_KBA_vals, props)
DG_RAN_prop <- multi_props(DG_RAN_vals, props)

TB_prop <- multi_props(TB_vals, props)
TB_KBA_prop <- multi_props(TB_KBA_vals, props)
TB_RAN_prop <- multi_props(TB_RAN_vals, props)

PL_prop <- multi_props(PL_vals, props)
PL_KBA_prop <- multi_props(PL_KBA_vals, props)
PL_RAN_prop <- multi_props(PL_RAN_vals, props)

AM_prop <- multi_props(AM_vals, props)
AM_KBA_prop <- multi_props(AM_KBA_vals, props)
AM_RAN_prop <- multi_props(AM_RAN_vals, props)

MF_prop <- multi_props(MF_vals, props)
MF_KBA_prop <- multi_props(MF_KBA_vals, props)
MF_RAN_prop <- multi_props(MF_RAN_vals, props)

MD_prop <- multi_props(MD_vals, props)
MD_KBA_prop <- multi_props(MD_KBA_vals, props)
MD_RAN_prop <- multi_props(MD_RAN_vals, props)

LH_prop <- multi_props(LH_vals, props)
LH_KBA_prop <- multi_props(LH_KBA_vals, props)
LH_RAN_prop <- multi_props(LH_RAN_vals, props)

VV_prop <- multi_props(VV_vals, props)
VV_KBA_prop <- multi_props(VV_KBA_vals, props)
VV_RAN_prop <- multi_props(VV_RAN_vals, props)

# MR_prop <- multi_props(MR_vals, props)
# MR_KBA_prop <- multi_props(MR_KBA_vals, props)
# MR_RAN_prop <- multi_props(MR_RAN_vals, props)

BT_prop <- multi_props(BT_vals, props)
BT_KBA_prop <- multi_props(BT_KBA_vals, props)
BT_RAN_prop <- multi_props(BT_RAN_vals, props)

HA_prop <- multi_props(HA_vals, props)
HA_KBA_prop <- multi_props(HA_KBA_vals, props)
HA_RAN_prop <- multi_props(HA_RAN_vals, props)

Total <- as.data.frame(rbind(PM_prop,VG_prop,DG_prop,TB_prop,PL_prop,AM_prop,MF_prop,MD_prop,LH_prop,VV_prop,BT_prop,HA_prop,
                             PM_KBA_prop, VG_KBA_prop, DG_KBA_prop, TB_KBA_prop, PL_KBA_prop,AM_KBA_prop,MF_KBA_prop,MD_KBA_prop,LH_KBA_prop,VV_KBA_prop,BT_KBA_prop,HA_KBA_prop,
                             PM_RAN_prop, VG_RAN_prop, DG_RAN_prop, TB_RAN_prop, PL_RAN_prop,AM_RAN_prop,MF_RAN_prop,MD_RAN_prop,LH_RAN_prop,VV_RAN_prop,BT_RAN_prop,HA_RAN_prop))
Type <- c(rep("Unmasked", 12), rep("Masked", 12), rep("Random", 12))
Species <- rep(c("P. megacephala", "V. germanica", "D. gazella",
                 "T. bicarinatum", "P. longicornis", "A. mellifera",
                 "M. floricola", "M. destructor", "L. humile",
                 "V. vulgaris", "B. terrestris", "H. arator"), 3)
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
Total_min <- Total %>%
  dplyr::filter(SiteSensitivity == 1 | SiteSensitivity == 0.98)
df <- data.frame(xmin = c(1,0.98,0.95,0.9,0.75,0.5,0.25),
                 xmax = c(0.98,0.95,0.9,0.75,0.5,0.25,0),
                 ymin = rep(-Inf,7),
                 ymax = rep(Inf,7))
p <- ggplot(Total, 
            aes(x = SiteSensitivity, 
                y = DistributionCoverage, 
                colour = Species,
                group = interaction(Species,Type))) +
  geom_line(stat = "identity", aes(linetype = Type)) +
  scale_linetype_manual(values = c("solid", "dotted", "longdash")) +
  scale_colour_manual(values = pnw_palette("Cascades", 12)) +
  ylab("Distribution Coverage (%)") +
  xlab("Site Sensitivity") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1,0),
                     trans = "reverse") +
  scale_y_continuous(expand = c(0,0))
p2 <-  p + 
  geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin=ymin, ymax=ymax),
            fill = rev(leg$colors), alpha = 0.2, inherit.aes = F) +
  font("legend.text", face = "italic")

#For inset
ins <- ggplot(Total_min, 
              aes(x = SiteSensitivity, 
                  y = DistributionCoverage, 
                  colour = Species,
                  group = interaction(Species,Type))) +
  geom_line(stat = "identity", aes(linetype = Type, size = Type)) +
  scale_linetype_manual(values = c("solid", "dotted", "longdash")) +
  scale_size_manual(values = c(1,1,1)) +
  scale_colour_manual(values = pnw_palette("Cascades", 12)) +
  ylab("Distribution Coverage (%)") +
  xlab("Site Sensitivity") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = alpha(leg$colors[7],0.2)),
        axis.text = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank()) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1,0.98),
                     trans = "reverse") +
  scale_y_continuous(expand = c(0,0))


