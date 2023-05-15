##IAS SDMs##
# spp_list <- c("Apis mellifera",  "Monomorium floricola", "Monomorium destructor", 
#               "Linepithema humile", "Vespula vulgaris", "Bombus terrestris", "Heteronychus arator", 
#               "Digitonthophagus gazella", "Pheidole megacephala", "Vespula germanica", 
#               "Tetramorium bicarinatum", "Paratrechina longicornis")
# not_run <- c("Solenopsis geminata","Polistes chinensis antennalis","Wasmannia auropunctata", "Apis cerana",
#              "Solenopsis invicta","Megachile rotundata")

#Plots of predicted distributions (ensemble committee averaging)
# lapply(spp_list[1:length(spp_list)], FUN = function(i){
#   
#   IAS_plot(i)
#   
# })

#Inspect the distribution of priorities within the IAS predicted distribution
#Use binary versions of ensemble committee averaging model


##CAZ with weights##
#Pheidole megacephala
PM_bin <- raster(file.path(regional_model_path, 
                           "Pheidole.megacephala", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Pheidole.megacephala",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
PM_bin2 <- PM_bin
PM_bin2[PM_bin2 != 0] <- 1

# Re-class zeros as NA
PM_bin[PM_bin == 0] <- NA
#get equal extent
PM_spec_bin <- resample(PM_bin, CAZ_wgt_var_ras, method = "ngb") 
PM_KBA_bin <- resample(PM_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
PM_RAN_bin <- resample(PM_bin, RAN_var_ras, method = "ngb") 
PM_area_bin <- resample(PM_bin, CAZ_area_wgt_var_ras, method = "ngb") 
PM_area_KBA_bin <- resample(PM_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
PM_area_RAN_bin <- resample(PM_bin, RAN_area_var_ras, method = "ngb")

#Vespula germanica
VG_bin <- raster(file.path(regional_model_path, 
                           "Vespula.germanica", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Vespula.germanica",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VG_bin2 <- VG_bin
VG_bin2[VG_bin2 != 0] <- 1

# Re-class zeros as NA
VG_bin[VG_bin == 0] <- NA
#get equal extent
VG_spec_bin <- resample(VG_bin, CAZ_wgt_var_ras, method = "ngb") 
VG_KBA_bin <- resample(VG_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
VG_RAN_bin <- resample(VG_bin, RAN_var_ras, method = "ngb") 
VG_area_bin <- resample(VG_bin, CAZ_area_wgt_var_ras, method = "ngb") 
VG_area_KBA_bin <- resample(VG_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
VG_area_RAN_bin <- resample(VG_bin, RAN_area_var_ras, method = "ngb") 

#Digitonthophagus gazella
DG_bin <- raster(file.path(regional_model_path, 
                           "Digitonthophagus.gazella", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Digitonthophagus.gazella",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
DG_bin2 <- DG_bin
DG_bin2[DG_bin2 != 0] <- 1

# Re-class zeros as NA
DG_bin[DG_bin == 0] <- NA
#get equal extent
DG_spec_bin <- resample(DG_bin, CAZ_wgt_var_ras, method = "ngb") 
DG_KBA_bin <- resample(DG_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
DG_RAN_bin <- resample(DG_bin, RAN_var_ras, method = "ngb") 
DG_area_bin <- resample(DG_bin, CAZ_area_wgt_var_ras, method = "ngb") 
DG_area_KBA_bin <- resample(DG_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
DG_area_RAN_bin <- resample(DG_bin, RAN_area_var_ras, method = "ngb") 

#Tetramorium bicarinatum
TB_bin <- raster(file.path(regional_model_path, 
                           "Tetramorium.bicarinatum", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Tetramorium.bicarinatum",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
TB_bin2 <- TB_bin
TB_bin2[TB_bin2 != 0] <- 1

# Re-class zeros as NA
TB_bin[TB_bin == 0] <- NA
#get equal extent
TB_spec_bin <- resample(TB_bin, CAZ_wgt_var_ras, method = "ngb") 
TB_KBA_bin <- resample(TB_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
TB_RAN_bin <- resample(TB_bin, RAN_var_ras, method = "ngb") 
TB_area_bin <- resample(TB_bin, CAZ_area_wgt_var_ras, method = "ngb") 
TB_area_KBA_bin <- resample(TB_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
TB_area_RAN_bin <- resample(TB_bin, RAN_area_var_ras, method = "ngb") 

#Paratrechina longicornis
PL_bin <- raster(file.path(regional_model_path, 
                           "Paratrechina.longicornis", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Paratrechina.longicornis",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
PL_bin2 <- PL_bin
PL_bin2[PL_bin2 != 0] <- 1

# Re-class zeros as NA
PL_bin[PL_bin == 0] <- NA
#get equal extent
PL_spec_bin <- resample(PL_bin, CAZ_wgt_var_ras, method = "ngb") 
PL_KBA_bin <- resample(PL_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
PL_RAN_bin <- resample(PL_bin, RAN_var_ras, method = "ngb") 
PL_area_bin <- resample(PL_bin, CAZ_area_wgt_var_ras, method = "ngb") 
PL_area_KBA_bin <- resample(PL_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
PL_area_RAN_bin <- resample(PL_bin, RAN_area_var_ras, method = "ngb") 

#Apis mellifera
AM_bin <- raster(file.path(regional_model_path, 
                           "Apis.mellifera", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Apis.mellifera",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
AM_bin2 <- AM_bin
AM_bin2[AM_bin2 != 0] <- 1

# Re-class zeros as NA
AM_bin[AM_bin == 0] <- NA
#get equal extent
AM_spec_bin <- resample(AM_bin, CAZ_wgt_var_ras, method = "ngb") 
AM_KBA_bin <- resample(AM_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
AM_RAN_bin <- resample(AM_bin, RAN_var_ras, method = "ngb") 
AM_area_bin <- resample(AM_bin, CAZ_area_wgt_var_ras, method = "ngb") 
AM_area_KBA_bin <- resample(AM_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
AM_area_RAN_bin <- resample(AM_bin, RAN_area_var_ras, method = "ngb")

#Monomorium floricola 
MF_bin <- raster(file.path(regional_model_path, 
                           "Monomorium.floricola", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Monomorium.floricola",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MF_bin2 <- MF_bin
MF_bin2[MF_bin2 != 0] <- 1

# Re-class zeros as NA
MF_bin[MF_bin == 0] <- NA
#get equal extent
MF_spec_bin <- resample(MF_bin, CAZ_wgt_var_ras, method = "ngb") 
MF_KBA_bin <- resample(MF_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
MF_RAN_bin <- resample(MF_bin, RAN_var_ras, method = "ngb") 
MF_area_bin <- resample(MF_bin, CAZ_area_wgt_var_ras, method = "ngb") 
MF_area_KBA_bin <- resample(MF_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
MF_area_RAN_bin <- resample(MF_bin, RAN_area_var_ras, method = "ngb") 

#Monomorium destructor
MD_bin <- raster(file.path(regional_model_path, 
                           "Monomorium.destructor", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Monomorium.destructor",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MD_bin2 <- MD_bin
MD_bin2[MD_bin2 != 0] <- 1

# Re-class zeros as NA
MD_bin[MD_bin == 0] <- NA
#get equal extent
MD_spec_bin <- resample(MD_bin, CAZ_wgt_var_ras, method = "ngb") 
MD_KBA_bin <- resample(MD_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
MD_RAN_bin <- resample(MD_bin, RAN_var_ras, method = "ngb") 
MD_area_bin <- resample(MD_bin, CAZ_area_wgt_var_ras, method = "ngb") 
MD_area_KBA_bin <- resample(MD_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
MD_area_RAN_bin <- resample(MD_bin, RAN_area_var_ras, method = "ngb") 

#Linepithema humile
LH_bin <- raster(file.path(regional_model_path, 
                           "Linepithema.humile", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Linepithema.humile",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
LH_bin2 <- LH_bin
LH_bin2[LH_bin2 != 0] <- 1

# Re-class zeros as NA
LH_bin[LH_bin == 0] <- NA
#get equal extent
LH_spec_bin <- resample(LH_bin, CAZ_wgt_var_ras, method = "ngb") 
LH_KBA_bin <- resample(LH_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
LH_RAN_bin <- resample(LH_bin, RAN_var_ras, method = "ngb") 
LH_area_bin <- resample(LH_bin, CAZ_area_wgt_var_ras, method = "ngb") 
LH_area_KBA_bin <- resample(LH_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
LH_area_RAN_bin <- resample(LH_bin, RAN_area_var_ras, method = "ngb") 

#Vespula vulgaris
VV_bin <- raster(file.path(regional_model_path, 
                           "Vespula.vulgaris", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Vespula.vulgaris",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VV_bin2 <- VV_bin
VV_bin2[VV_bin2 != 0] <- 1

# Re-class zeros as NA
VV_bin[VV_bin == 0] <- NA
#get equal extent
VV_spec_bin <- resample(VV_bin, CAZ_wgt_var_ras, method = "ngb") 
VV_KBA_bin <- resample(VV_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
VV_RAN_bin <- resample(VV_bin, RAN_var_ras, method = "ngb") 
VV_area_bin <- resample(VV_bin, CAZ_area_wgt_var_ras, method = "ngb") 
VV_area_KBA_bin <- resample(VV_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
VV_area_RAN_bin <- resample(VV_bin, RAN_area_var_ras, method = "ngb") 

#Bombus terrestris
BT_bin <- raster(file.path(regional_model_path, 
                           "Bombus.terrestris", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Bombus.terrestris",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
BT_bin2 <- BT_bin
BT_bin2[BT_bin2 != 0] <- 1

# Re-class zeros as NA
BT_bin[BT_bin == 0] <- NA
#get equal extent
BT_spec_bin <- resample(BT_bin, CAZ_wgt_var_ras, method = "ngb") 
BT_KBA_bin <- resample(BT_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
BT_RAN_bin <- resample(BT_bin, RAN_var_ras, method = "ngb") 
BT_area_bin <- resample(BT_bin, CAZ_area_wgt_var_ras, method = "ngb") 
BT_area_KBA_bin <- resample(BT_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
BT_area_RAN_bin <- resample(BT_bin, RAN_area_var_ras, method = "ngb")

#Heteronychus arator
HA_bin <- raster(file.path(regional_model_path, 
                           "Heteronychus.arator", 
                           "proj_regional", 
                           "individual_projections", 
                           paste0("Heteronychus.arator",
                                  "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
HA_bin2 <- HA_bin
HA_bin2[HA_bin2 != 0] <- 1

# Re-class zeros as NA
HA_bin[HA_bin == 0] <- NA
#get equal extent
HA_spec_bin <- resample(HA_bin, CAZ_wgt_var_ras, method = "ngb") 
HA_KBA_bin <- resample(HA_bin, CAZ_wgt_KBA_inv_var_ras, method = "ngb") 
HA_RAN_bin <- resample(HA_bin, RAN_var_ras, method = "ngb") 
HA_area_bin <- resample(HA_bin, CAZ_area_wgt_var_ras, method = "ngb") 
HA_area_KBA_bin <- resample(HA_bin, CAZ_area_wgt_KBA_inv_var_ras, method = "ngb") 
HA_area_RAN_bin <- resample(HA_bin, RAN_area_var_ras, method = "ngb") 


#Overlap among all IAS
ias_sum <- sum(PM_bin2, VG_bin2, DG_bin2, TB_bin2, PL_bin2, AM_bin2, MF_bin2, 
               MD_bin2, LH_bin2, VV_bin2, BT_bin2, HA_bin2)

ias_sum_one <- ias_sum
ias_sum_one[ias_sum_one != 1] <- NA

ias_sum_sp <- resample(ias_sum, CAZ_wgt_var_ras, method = "ngb")
ias_sum_sp_ar <- resample(ias_sum, CAZ_area_wgt_var_ras, method = "ngb")

ras_st <- st_as_stars(ias_sum_sp_ar)
ras_sf <- st_as_sf(ras_st)
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

non.na <- length(!is.na(ias_sum_sp)) #number of non NA cells
no.cells <- c()
perc.total <- c()

for(i in 0:11){
  
  c <- length(which(values(ias_sum_sp == i)))
  no.cells <- c(no.cells, c)
  p <- c/non.na*100
  perc.total <- c(perc.total, p)
  
}

richness <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
totaltable <- data.frame("Richness" = richness, 
                         "Number of cells" = no.cells, 
                         "Percent total" = round(prop.total,2))

pretty_table <- formattable(totaltable,
            align = c("c", "c", "c"),
            list(Richness = color_tile("#e69b99","#2c6184")))

ggplot(totaltable, aes(log(no.cells), log(prop.total))) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_x_reverse(breaks = log(no.cells), labels = richness) +
  xlab("Alien richness") +
  ylab("Percent of total area (log)")

#SDM evaluations
eval_scores <- read.table(file.path(regional_model_path, "Table_IAS_regional_accuracy.txt"), header = T)
pretty_eval_table <- formattable(eval_scores,
                            align = c("l", "c", "c", "c", "c"))

#Clear console
cat("\014")