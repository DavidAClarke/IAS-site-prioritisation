#Results pertaining to susceptible site identification
##IAS SDMs##
regional_model_path <- "G:/Chapter_3/SpatialData/IAS_distributions/IAS_regional" #external hard drive
spp_list <- c("Apis mellifera",  "Monomorium floricola", "Monomorium destructor", 
              "Linepithema humile", "Vespula vulgaris", "Megachile rotundata", 
              "Bombus terrestris", "Heteronychus arator", "Digitonthophagus gazella", 
              "Pheidole megacephala", "Vespula germanica", "Tetramorium bicarinatum", 
              "Paratrechina longicornis")
not_run <- c("Solenopsis geminata","Polistes chinensis antennalis","Wasmannia auropunctata", "Apis cerana",
             "Solenopsis invicta")

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

#Apis mellifera
AM_bin <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
AM_bin2 <- AM_bin
AM_bin2[AM_bin2 != 0] <- 1

# Re-class zeros as NA
AM_bin[AM_bin == 0] <- NA
AM_bin <- resample(AM_bin, CAZ_wgt_var_ras) #get equal extent
AM_bin_p <- rasterToPolygons(AM_bin, dissolve = T)
AM_bin_sf <- st_as_sf(AM_bin_p)
AM_bin_sf_bb <- st_as_sfc(st_bbox(AM_bin_sf))

#Monomorium floricola 
MF_bin <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MF_bin2 <- MF_bin
MF_bin2[MF_bin2 != 0] <- 1

# Re-class zeros as NA
MF_bin[PL_bin == 0] <- NA
MF_bin <- resample(MF_bin, CAZ_wgt_var_ras) #get equal extent
MF_bin_p <- rasterToPolygons(MF_bin, dissolve = T)
MF_bin_sf <- st_as_sf(MF_bin_p)
MF_bin_sf_bb <- st_as_sfc(st_bbox(MF_bin_sf))

#Monomorium destructor
MD_bin <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MD_bin2 <- MD_bin
MD_bin2[MD_bin2 != 0] <- 1

# Re-class zeros as NA
MD_bin[PL_bin == 0] <- NA
MD_bin <- resample(MD_bin, CAZ_wgt_var_ras) #get equal extent
MD_bin_p <- rasterToPolygons(MD_bin, dissolve = T)
MD_bin_sf <- st_as_sf(MD_bin_p)
MD_bin_sf_bb <- st_as_sfc(st_bbox(MD_bin_sf))

#Linepithema humile
LH_bin <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
LH_bin2 <- LH_bin
LH_bin2[LH_bin2 != 0] <- 1

# Re-class zeros as NA
LH_bin[LH_bin == 0] <- NA
LH_bin <- resample(LH_bin, CAZ_wgt_var_ras) #get equal extent
LH_bin_p <- rasterToPolygons(LH_bin, dissolve = T)
LH_bin_sf <- st_as_sf(LH_bin_p)
LH_bin_sf_bb <- st_as_sfc(st_bbox(LH_bin_sf))

#Vespula vulgaris
VV_bin <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VV_bin2 <- VV_bin
VV_bin2[VV_bin2 != 0] <- 1

# Re-class zeros as NA
VV_bin[PL_bin == 0] <- NA
VV_bin <- resample(VV_bin, CAZ_wgt_var_ras) #get equal extent
VV_bin_p <- rasterToPolygons(VV_bin, dissolve = T)
VV_bin_sf <- st_as_sf(VV_bin_p)
VV_bin_sf_bb <- st_as_sfc(st_bbox(VV_bin_sf))

#Megachile rotundata 
MR_bin <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MR_bin2 <- MR_bin
MR_bin2[MR_bin2 != 0] <- 1

# Re-class zeros as NA
MR_bin[MR_bin == 0] <- NA
MR_bin <- resample(MR_bin, CAZ_wgt_var_ras) #get equal extent
MR_bin_p <- rasterToPolygons(MR_bin, dissolve = T)
MR_bin_sf <- st_as_sf(MR_bin_p)
MR_bin_sf_bb <- st_as_sfc(st_bbox(MR_bin_sf))

#Bombus terrestris
BT_bin <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
BT_bin2 <- BT_bin
BT_bin2[BT_bin2 != 0] <- 1

# Re-class zeros as NA
BT_bin[BT_bin == 0] <- NA
BT_bin <- resample(BT_bin, CAZ_wgt_var_ras) #get equal extent
BT_bin_p <- rasterToPolygons(BT_bin, dissolve = T)
BT_bin_sf <- st_as_sf(BT_bin_p)
BT_bin_sf_bb <- st_as_sfc(st_bbox(BT_bin_sf))

#Heteronychus arator
HA_bin <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
HA_bin2 <- HA_bin
HA_bin2[HA_bin2 != 0] <- 1

# Re-class zeros as NA
HA_bin[HA_bin == 0] <- NA
HA_bin <- resample(HA_bin, CAZ_wgt_var_ras) #get equal extent
HA_bin_p <- rasterToPolygons(HA_bin, dissolve = T)
HA_bin_sf <- st_as_sf(HA_bin_p)
HA_bin_sf_bb <- st_as_sfc(st_bbox(HA_bin_sf))


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
MR_bin_RAN <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
# Re-class zeros as NA
MR_bin_RAN[MR_bin_RAN == 0] <- NA
MR_bin_RAN <- resample(MR_bin_RAN, RAN_wgt_var_ras) #get equal extent

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


                                  ##CAZ with areas and weights##
#Pheidole megacephala
PM_bin <- raster(file.path(regional_model_path, "Pheidole.megacephala", "proj_regional", "individual_projections", paste0("Pheidole.megacephala","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
PM_bin2 <- PM_bin
PM_bin2[PM_bin2 != 0] <- 1

# Re-class zeros as NA
PM_bin[PM_bin == 0] <- NA
PM_bin <- resample(PM_bin, CAZ_area_wgt_var_ras) #get equal extent
PM_bin_p <- rasterToPolygons(PM_bin, dissolve = T)
PM_bin_sf <- st_as_sf(PM_bin_p)
PM_bin_sf_bb <- st_as_sfc(st_bbox(PM_bin_sf))


#Vespula germanica
VG_bin <- raster(file.path(regional_model_path, "Vespula.germanica", "proj_regional", "individual_projections", paste0("Vespula.germanica","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VG_bin2 <- VG_bin
VG_bin2[VG_bin2 != 0] <- 1

# Re-class zeros as NA
VG_bin[VG_bin == 0] <- NA
VG_bin <- resample(VG_bin, CAZ_area_wgt_var_ras) #get equal extent
VG_bin_p <- rasterToPolygons(VG_bin, dissolve = T)
VG_bin_sf <- st_as_sf(VG_bin_p)
VG_bin_sf_bb <- st_as_sfc(st_bbox(VG_bin_sf))

#Digitonthophagus gazella
DG_bin <- raster(file.path(regional_model_path, "Digitonthophagus.gazella", "proj_regional", "individual_projections", paste0("Digitonthophagus.gazella","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
DG_bin2 <- DG_bin
DG_bin2[DG_bin2 != 0] <- 1

# Re-class zeros as NA
DG_bin[DG_bin == 0] <- NA
DG_bin <- resample(DG_bin, CAZ_area_wgt_var_ras) #get equal extent
DG_bin_p <- rasterToPolygons(DG_bin, dissolve = T)
DG_bin_sf <- st_as_sf(DG_bin_p)
DG_bin_sf_bb <- st_as_sfc(st_bbox(DG_bin_sf))

#Tetramorium bicarinatum
TB_bin <- raster(file.path(regional_model_path, "Tetramorium.bicarinatum", "proj_regional", "individual_projections", paste0("Tetramorium.bicarinatum","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
TB_bin2 <- TB_bin
TB_bin2[TB_bin2 != 0] <- 1

# Re-class zeros as NA
TB_bin[TB_bin == 0] <- NA
TB_bin <- resample(TB_bin, CAZ_area_wgt_var_ras) #get equal extent
TB_bin_p <- rasterToPolygons(TB_bin, dissolve = T)
TB_bin_sf <- st_as_sf(TB_bin_p)
TB_bin_sf_bb <- st_as_sfc(st_bbox(TB_bin_sf))

#Paratrechina longicornis
PL_bin <- raster(file.path(regional_model_path, "Paratrechina.longicornis", "proj_regional", "individual_projections", paste0("Paratrechina.longicornis","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
PL_bin2 <- PL_bin
PL_bin2[PL_bin2 != 0] <- 1

# Re-class zeros as NA
PL_bin[PL_bin == 0] <- NA
PL_bin <- resample(PL_bin, CAZ_area_wgt_var_ras) #get equal extent
PL_bin_p <- rasterToPolygons(PL_bin, dissolve = T)
PL_bin_sf <- st_as_sf(PL_bin_p)
PL_bin_sf_bb <- st_as_sfc(st_bbox(PL_bin_sf))

#Apis mellifera
AM_bin <- raster(file.path(regional_model_path, "Apis.mellifera", "proj_regional", "individual_projections", paste0("Apis.mellifera","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
AM_bin2 <- AM_bin
AM_bin2[AM_bin2 != 0] <- 1

# Re-class zeros as NA
AM_bin[AM_bin == 0] <- NA
AM_bin <- resample(AM_bin, CAZ_area_wgt_var_ras) #get equal extent
AM_bin_p <- rasterToPolygons(AM_bin, dissolve = T)
AM_bin_sf <- st_as_sf(AM_bin_p)
AM_bin_sf_bb <- st_as_sfc(st_bbox(AM_bin_sf))

#Monomorium floricola 
MF_bin <- raster(file.path(regional_model_path, "Monomorium.floricola", "proj_regional", "individual_projections", paste0("Monomorium.floricola","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MF_bin2 <- MF_bin
MF_bin2[MF_bin2 != 0] <- 1

# Re-class zeros as NA
MF_bin[PL_bin == 0] <- NA
MF_bin <- resample(MF_bin, CAZ_area_wgt_var_ras) #get equal extent
MF_bin_p <- rasterToPolygons(MF_bin, dissolve = T)
MF_bin_sf <- st_as_sf(MF_bin_p)
MF_bin_sf_bb <- st_as_sfc(st_bbox(MF_bin_sf))

#Monomorium destructor
MD_bin <- raster(file.path(regional_model_path, "Monomorium.destructor", "proj_regional", "individual_projections", paste0("Monomorium.destructor","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MD_bin2 <- MD_bin
MD_bin2[MD_bin2 != 0] <- 1

# Re-class zeros as NA
MD_bin[PL_bin == 0] <- NA
MD_bin <- resample(MD_bin, CAZ_area_wgt_var_ras) #get equal extent
MD_bin_p <- rasterToPolygons(MD_bin, dissolve = T)
MD_bin_sf <- st_as_sf(MD_bin_p)
MD_bin_sf_bb <- st_as_sfc(st_bbox(MD_bin_sf))

#Linepithema humile
LH_bin <- raster(file.path(regional_model_path, "Linepithema.humile", "proj_regional", "individual_projections", paste0("Linepithema.humile","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
LH_bin2 <- LH_bin
LH_bin2[LH_bin2 != 0] <- 1

# Re-class zeros as NA
LH_bin[LH_bin == 0] <- NA
LH_bin <- resample(LH_bin, CAZ_area_wgt_var_ras) #get equal extent
LH_bin_p <- rasterToPolygons(LH_bin, dissolve = T)
LH_bin_sf <- st_as_sf(LH_bin_p)
LH_bin_sf_bb <- st_as_sfc(st_bbox(LH_bin_sf))

#Vespula vulgaris
VV_bin <- raster(file.path(regional_model_path, "Vespula.vulgaris", "proj_regional", "individual_projections", paste0("Vespula.vulgaris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
VV_bin2 <- VV_bin
VV_bin2[VV_bin2 != 0] <- 1

# Re-class zeros as NA
VV_bin[PL_bin == 0] <- NA
VV_bin <- resample(VV_bin, CAZ_area_wgt_var_ras) #get equal extent
VV_bin_p <- rasterToPolygons(VV_bin, dissolve = T)
VV_bin_sf <- st_as_sf(VV_bin_p)
VV_bin_sf_bb <- st_as_sfc(st_bbox(VV_bin_sf))

#Megachile rotundata 
MR_bin <- raster(file.path(regional_model_path, "Megachile.rotundata", "proj_regional", "individual_projections", paste0("Megachile.rotundata","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
MR_bin2 <- MR_bin
MR_bin2[MR_bin2 != 0] <- 1

# Re-class zeros as NA
MR_bin[MR_bin == 0] <- NA
MR_bin <- resample(MR_bin, CAZ_area_wgt_var_ras) #get equal extent
MR_bin_p <- rasterToPolygons(MR_bin, dissolve = T)
MR_bin_sf <- st_as_sf(MR_bin_p)
MR_bin_sf_bb <- st_as_sfc(st_bbox(MR_bin_sf))

#Bombus terrestris
BT_bin <- raster(file.path(regional_model_path, "Bombus.terrestris", "proj_regional", "individual_projections", paste0("Bombus.terrestris","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
BT_bin2 <- BT_bin
BT_bin2[BT_bin2 != 0] <- 1

# Re-class zeros as NA
BT_bin[BT_bin == 0] <- NA
BT_bin <- resample(BT_bin, CAZ_area_wgt_var_ras) #get equal extent
BT_bin_p <- rasterToPolygons(BT_bin, dissolve = T)
BT_bin_sf <- st_as_sf(BT_bin_p)
BT_bin_sf_bb <- st_as_sfc(st_bbox(BT_bin_sf))

#Heteronychus arator
HA_bin <- raster(file.path(regional_model_path, "Heteronychus.arator", "proj_regional", "individual_projections", paste0("Heteronychus.arator","_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri")))
HA_bin2 <- HA_bin
HA_bin2[HA_bin2 != 0] <- 1

# Re-class zeros as NA
HA_bin[HA_bin == 0] <- NA
HA_bin <- resample(HA_bin, CAZ_area_wgt_var_ras) #get equal extent
HA_bin_p <- rasterToPolygons(HA_bin, dissolve = T)
HA_bin_sf <- st_as_sf(HA_bin_p)
HA_bin_sf_bb <- st_as_sfc(st_bbox(HA_bin_sf))