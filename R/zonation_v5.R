## Zonation 5 files

pkgs <- c("tidyverse", "here")
lapply(pkgs, require, character.only = T)

zonation_path <- here("zonation")
species_path <- here("species_lowres")
species_area_path <- here("species_area_lowres")

# dir.create(here(zonation_path, "species_scenarios"))
# dir.create(here(zonation_path, "species_area_scenarios"))

###############################################################################
## Species only scenarios
spp_list <- read.csv(here(zonation_path, "maxent_spp_list_upd.csv")) %>%
    filter(acceptedName != "Endospermum myrmecophilum" & 
           acceptedName != "Ozimops petersi") #files are corrupted

# CAZMAX is equivalent to CAZ in Zonation V4
weight_df <- data.frame(equal = rep(1, nrow(spp_list)), spp_list[,11:15])

for(i in 1:ncol(weight_df)){
  
  nm <- names(weight_df)[i]
  
  dir.create(here(zonation_path, "species_scenarios", paste0("species_",nm)))
  
  variant_path <- here(zonation_path, "species_scenarios", paste0("species_",nm))
  
  feature_list <- data.frame(weight = weight_df[,i],
                             group = spp_list$classGroupNum,
                             filename = list.files(species_path, full.names = T))
  
  write.table(feature_list, 
              file = here(variant_path,"features.txt"), row.names = F)
  
  Settings <- paste("feature list file =", here(variant_path,"features.txt"))
  fileConn <- file(here(variant_path, paste0("species_",nm,".z5")))
  writeLines(Settings, fileConn)
  close(fileConn)
  
}

## With KBA mask
for(i in 1:ncol(weight_df)){
  
  nm <- names(weight_df)[i]
  
  dir.create(here(zonation_path, "species_scenarios", paste0("species_",nm,"_KBA")))
  
  variant_path <- here(zonation_path, "species_scenarios", paste0("species_",nm, "_KBA"))
  
  mask_path <- here("hierarchic_mask.tif")
  
  feature_list <- data.frame(weight = weight_df[,i],
                             group = spp_list$classGroupNum,
                             filename = list.files(species_path, full.names = T))
  
  write.table(feature_list, 
              file = here(variant_path,"features.txt"), row.names = F)
  
  Settings <- paste(
    paste("feature list file =", here(variant_path,"features.txt")),
    paste("hierarchic mask layer =", mask_path), 
    sep = "\n")
  
  fileConn <- file(here(variant_path, paste0("species_",nm,"_KBA.z5")))
  writeLines(Settings, fileConn)
  close(fileConn)
  
}
###############################################################################
## Species area scenarios
spp_list <- read.csv(here(zonation_path, "spp_area_list_upd.csv")) %>%
  filter(acceptedName != "Morelia spilota") #files are corrupted

weight_df <- data.frame(equal = rep(1, nrow(spp_list)), spp_list[,11:15])

for(i in 1:ncol(weight_df)){
  
  nm <- names(weight_df)[i]
  
  dir.create(here(zonation_path, "species_area_scenarios", paste0("species_area_",nm)))
  
  variant_path <- here(zonation_path, "species_area_scenarios", paste0("species_area_",nm))
  
  feature_list <- data.frame(weight = weight_df[,i],
                             group = spp_list$classGroupNum,
                             filename = list.files(species_area_path, full.names = T))
  
  write.table(feature_list, 
              file = here(variant_path,"features.txt"), row.names = F)
  
  Settings <- paste("feature list file =", here(variant_path,"features.txt"))
  fileConn <- file(here(variant_path, paste0("species_area_",nm,".z5")))
  writeLines(Settings, fileConn)
  close(fileConn)
  
}

## With KBA mask
for(i in 1:ncol(weight_df)){
  
  nm <- names(weight_df)[i]
  
  dir.create(here(zonation_path, "species_area_scenarios", paste0("species_area_",nm,"_KBA")))
  
  variant_path <- here(zonation_path, "species_area_scenarios", paste0("species_area_",nm, "_KBA"))
  
  mask_path <- here("hierarchic_mask.tif")
  
  feature_list <- data.frame(weight = weight_df[,i],
                             group = spp_list$classGroupNum,
                             filename = list.files(species_area_path, full.names = T))
  
  write.table(feature_list, 
              file = here(variant_path,"features.txt"), row.names = F)
  
  Settings <- paste(
    paste("feature list file =", here(variant_path,"features.txt")),
    paste("hierarchic mask layer =", mask_path), 
    sep = "\n")
  
  fileConn <- file(here(variant_path, paste0("species_area_",nm,"_KBA.z5")))
  writeLines(Settings, fileConn)
  close(fileConn)
  
}
