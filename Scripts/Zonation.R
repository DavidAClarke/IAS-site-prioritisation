############################################################################################
#                                                                                          #
##          Zonation will not work through R if ANY file paths used have ANY spaces       ##
#                                                                                          #
############################################################################################

#Required Libraries
library(tidyverse)
library(raster)
library(zonator)

#Project path
zonation_path <- "E:/Zonation" #set path

#Input data file path
Input_data_path <- paste0(zonation_path,"/", "Input_data")

#Input data (raster) paths (species or species+area)
species_path <- paste0(Input_data_path,"/", "species")
species_area_path <- paste0(Input_data_path,"/", "species_area")

#Hierarchical mask
KBA_msk <- paste0(Input_data_path, "/","KBA_hiermask.tif") #or _inv.tif
  
#Load species list, class names, threat status and weight
#need to edit to add plant species and add more weight to IAS threatened species
spp_list_sp <- read.csv(paste0(zonation_path,"/","maxent_spp_list.csv"))
spp_list_spa <- read.csv(paste0(zonation_path,"/","spp_area_list.csv"))

#Define variant names 
sp_variant_names <- c("species_CAZ",
                      "species_wgt_CAZ",
                      "species_wgt_CAZ_KBA")
spa_variant_names <- c("species_area_CAZ",
                       "species_area_wgt_CAZ",
                       "species_area_wgt_CAZ_KBA")


                             ##Create new Zonation project ##
#create general settings file
Settings <- "[Settings]
removal rule = 1
warp factor = 100 
edge removal = 1
add edge points = 0
memory save mode = 1"
fileConn <- file(paste0(zonation_path,"/","CAZ_settings.dat"))
writeLines(Settings, fileConn)
close(fileConn)


##create template spp_file's (feature file). 
sp_wgts <- as.numeric(spp_list_sp$weight)
spa_wgts <- as.numeric(spp_list_spa$weight)
create_spp(
  filename = paste0(zonation_path,"/","species.spp"),
  weight = 1, 
  alpha = 1,
  bqp = 1,
  bqp_p = 1,
  cellrem = 0.25,
  spp_file_dir = species_path, 
  recursive = FALSE,
  spp_file_pattern = ".+\\.(tif|img)$",
  override_path = NULL
)
create_spp(
  filename = paste0(zonation_path,"/","species_area.spp"),
  weight = 1, 
  alpha = 1,
  bqp = 1,
  bqp_p = 1,
  cellrem = 0.25,
  spp_file_dir = species_area_path, 
  recursive = FALSE,
  spp_file_pattern = ".+\\.(tif|img)$",
  override_path = NULL
)

#Place groups text file in appropriate place
#species only group
groups_output_sp <- spp_list_sp$classGroupNum
groups_condition_sp <- rep(-1, length(spp_list_sp$X))
groups_retention_sp <- rep(-1, length(spp_list_sp$X))
retention_mode_sp <- rep(-1, length(spp_list_sp$X))
arb_kernal_sp <- rep(-1, length(spp_list_sp$X))
groups_file_sp <- cbind(groups_output_sp, groups_condition_sp, groups_retention_sp, 
                        retention_mode_sp, arb_kernal_sp)
write.table(groups_file_sp, 
            file = paste0(zonation_path,"/","species_proj","/", "species_groups.txt"), 
            sep = " ", 
            col.names = F,
            row.names = F)
groups_path_sp <- "species_groups.txt"

#species and area group
groups_output_spa <- spp_list_spa$classGroupNum
groups_condition_spa <- rep(-1, length(spp_list_spa$X))
groups_retention_spa <- rep(-1, length(spp_list_spa$X))
retention_mode_spa <- rep(-1, length(spp_list_spa$X))
arb_kernal_spa <- rep(-1, length(spp_list_spa$X))
groups_file_spa <- cbind(groups_output_spa, groups_condition_spa, groups_retention_spa, 
                         retention_mode_spa, arb_kernal_spa)
write.table(groups_file_spa, 
            file = paste0(zonation_path,"/","species_area_proj","/", "species_area_groups.txt"), 
            sep = " ", 
            col.names = F,
            row.names = F)
groups_path_spa <- "species_area_groups.txt"

# Create species only Zonation project
sp_sensitive_site_ID <- create_zproject(name = "species_proj", 
                dir = zonation_path,
                variants = sp_variant_names,
                spp_template_dir = species_path,
                spp_template_file = file.path(zonation_path,"species.spp"),
                dat_template_file = file.path(zonation_path, "settings.dat"))

# Create species and area Zonation project
spa_sensitive_site_ID <- create_zproject(name = "species_area_proj", 
                                     dir = zonation_path,
                                     variants = spa_variant_names,
                                     spp_template_dir = species_area_path,
                                     spp_template_file = file.path(zonation_path,"species_area.spp"),
                                     dat_template_file = file.path(zonation_path, "settings.dat"))
                
#Load project in to edit variants
sp_sensitive_site_ID <- load_zproject(file.path(zonation_path, "species_proj"))
spa_sensitive_site_ID <- load_zproject(file.path(zonation_path, "species_area_proj"))

#Edit species only variants
variant <- get_variant(sp_sensitive_site_ID, 1)
# variant_spp_data <- sppdata(variant)
# wgts <- as.numeric(spp_list_sp$weight)
# variant_spp_data$weight <- wgts
# sppdata(variant) <- variant_spp_data
variant <- set_dat_param(variant, "use groups", 1)
variant <- set_dat_param(variant, "groups file", groups_path_sp)
groups(variant) <- as.vector(groups_file_sp[,1])
save_zvariant(variant, dir = paste0(zonation_path,"/", "species_proj"), overwrite = T, debug_msg = F)

variant <- get_variant(sp_sensitive_site_ID, 2)
variant_spp_data <- sppdata(variant)
wgts <- as.numeric(spp_list_sp$weight)
variant_spp_data$weight <- wgts
sppdata(variant) <- variant_spp_data
variant <- set_dat_param(variant, "use groups", 1)
variant <- set_dat_param(variant, "groups file", groups_path_sp)
groups(variant) <- as.vector(groups_file_sp[,1])
save_zvariant(variant, dir = paste0(zonation_path,"/", "species_proj"), overwrite = T, debug_msg = F)

variant <- get_variant(sp_sensitive_site_ID, 3)
variant_spp_data <- sppdata(variant)
wgts <- as.numeric(spp_list_sp$weight)
variant_spp_data$weight <- wgts
sppdata(variant) <- variant_spp_data
variant <- set_dat_param(variant, "use mask", 1)
variant <- set_dat_param(variant, "mask file", KBA_msk)
variant <- set_dat_param(variant, "use groups", 1)
variant <- set_dat_param(variant, "groups file", groups_path_sp)
groups(variant) <- as.vector(groups_file_sp[,1])
save_zvariant(variant, dir = paste0(zonation_path,"/", "species_proj"), overwrite = T, debug_msg = F)

#Edit bat file to add "use threads"
for(i in 1:3){
d <- read.delim(file = paste0(zonation_path, "/", "species_proj","/",variant@name[i],".bat"), sep = " ", header = F)
d$V13 <- apply(d, 1, function(x) paste(x, collapse = " "))
d <- d$V13
d <- paste(d, "--use-threads=4")
write.table(d, file = paste0(zonation_path, "/", "species_proj","/",variant@name[i],".bat"), sep = " ", col.names = F, row.names = F, append = F, quote = F)
}

#run
for(i in 1:3){
run_bat(paste0(zonation_path, "/", "species_proj","/",variant@name[i],".bat"), exe = "zig4")
}

#Edit species and area variants
variant <- get_variant(spa_sensitive_site_ID, 1)
# variant_spp_data <- sppdata(variant)
# wgts <- as.numeric(spp_list_spa$weight)
# variant_spp_data$weight <- wgts
# sppdata(variant) <- variant_spp_data
variant <- set_dat_param(variant, "use groups", 1)
variant <- set_dat_param(variant, "groups file", groups_path_spa)
groups(variant) <- as.vector(groups_file_spa[,1])
save_zvariant(variant, dir = paste0(zonation_path,"/", "species_area_proj"), overwrite = T, debug_msg = F)

variant <- get_variant(sp_sensitive_site_ID, 2)
variant_spp_data <- sppdata(variant)
wgts <- as.numeric(spp_list_spa$weight)
variant_spp_data$weight <- wgts
sppdata(variant) <- variant_spp_data
variant <- set_dat_param(variant, "use groups", 1)
variant <- set_dat_param(variant, "groups file", groups_path_spa)
groups(variant) <- as.vector(groups_file_spa[,1])
save_zvariant(variant, dir = paste0(zonation_path,"/", "species_proj"), overwrite = T, debug_msg = F)

variant <- get_variant(sp_sensitive_site_ID, 3)
variant_spp_data <- sppdata(variant)
wgts <- as.numeric(spp_list_spa$weight)
variant_spp_data$weight <- wgts
sppdata(variant) <- variant_spp_data
variant <- set_dat_param(variant, "use mask", 1)
variant <- set_dat_param(variant, "mask file", KBA_msk)
variant <- set_dat_param(variant, "use groups", 1)
variant <- set_dat_param(variant, "groups file", groups_path_spa)
groups(variant) <- as.vector(groups_file_spa[,1])
save_zvariant(variant, dir = paste0(zonation_path,"/", "species_area_proj"), overwrite = T, debug_msg = F)

#Edit bat file to add "use threads"
for(i in 1:3){
d <- read.delim(file = paste0(zonation_path, "/", "species_area_proj","/",variant@name[i],".bat"), sep = " ", header = F)
d$V13 <- apply(d, 1, function(x) paste(x, collapse = " "))
d <- d$V13
d <- paste(d, "--use-threads=4")
write.table(d, file = paste0(zonation_path, "/", "species_area_proj","/",variant@name[i],".bat"), sep = " ", col.names = F, row.names = F, append = F, quote = F)
}

#run
for(i in 1:3){
  run_bat(paste0(zonation_path, "/", "species_area_proj","/",variant@name[i],".bat"), exe = "zig4")
}

