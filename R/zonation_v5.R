## Zonation 5 files

rm(list = ls())

pkgs <- c("tidyverse", "here")
lapply(pkgs, require, character.only = T)

zonation_path <- here(dirname(here()), "data", "zonation")
species_path <- here(dirname(here()), "data","species_lowres")
species_area_path <- here(dirname(here()), "data","species_area_lowres")

# dir.create(here(zonation_path, "species_scenarios"))
# dir.create(here(zonation_path, "species_area_scenarios"))

###############################################################################
## Species only scenarios
spp_list <- read.csv(here(dirname(here()), "data", "maxent_spp_list_upd.csv")) %>%
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
spp_list <- read.csv(here(dirname(here()), "data", "spp_area_list_upd.csv")) %>%
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

#################################### NEEDS TESTING #############################
## Create bat files (testing; create loop)
z5_path <- "C:/Program Files/Zonation5/z5.exe"
settings_file <- here(zonation_path, "species_area_scenarios", 
                      "species_area_equal_KBA", "species_area_equal_KBA.z5")
output_dir <- here(zonation_path, "species_area_scenarios", 
                   "species_area_equal_KBA", "output")
z5_options <- "-wg -h -b --mode==CAZMAX"
run1 <- paste(z5_path, z5_options, settings_file, ouput_dir)
run2 <- paste(z5_path, z5_options, settings_file, ouput_dir)
fileConn <- file(here(zonation_path, "species_area_scenarios", "species_area_scenarios.bat"))
writeLines(c(run1,run2), fileConn)
close(fileConn)

## Look at this approach. Perhaps suitable for loop?
# https://sparkbyexamples.com/r-programming/r-write-lines-to-text-file/#:~:text=R%20base%20function%20writeLines(),cat()%20methods%20explained%20below.
# Example 5 - Using cat
cat("I Love R Programming",file="/Users/admin/textFile.txt",sep="\n")
cat("I live in USA",file="/Users/admin/textFile.txt",append=TRUE) #note append=T

# specify file
my_file <- "/Users/admin/textFile.txt"
if(file.exists(my_file)){
  a = T
  cat(run2, file = my_file, append = a)
} else {
  a = F
  cat(run2, file = my_file, append = a)
  
}

# The shell() function can be used to run a .bat file. Maybe use the translate argument (see help)
# Translate = T changes "/" to "\".
# shell(here(zonation_path, "species_area_scenarios", "species_area_scenarios.bat"), translate = F)