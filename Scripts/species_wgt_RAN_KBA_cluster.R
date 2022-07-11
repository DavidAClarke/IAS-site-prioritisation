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
#zonation_path <- "C:/Users/david/Documents/PhD/Chapter_3/Zonation"
#zonation_path <- "D:/Chapter_3/Chapter_3/Zonation" #virtual computer
zonation_path <- "/home/dcla0008/nc57_scratch/Zonation"
#zonation_path <- "F:/Zonation" #external

#Input data file path
Input_data_path <- paste0(zonation_path,"/", "Input_data")

#Input data (raster) paths (species or species+area)
species_path <- paste0(Input_data_path,"/", "species_lowres")
#species_area_path <- paste0(Input_data_path,"/", "species_area")

#Hierarchical mask
#KBA_msk <- paste0(Input_data_path, "/","KBA_hiermask.tif")
KBA_msk <- paste0(Input_data_path, "/","KBA_hiermask_inv.tif")

#Load species list, class names, threat status and weight
#need to edit to add plant species and add more weight to IAS threatened species
spp_list <- read.csv(paste0(zonation_path,"/","maxent_spp_list.csv"))
#spp_list <- read.csv(paste0(zonation_path,"/","spp_area_list.csv"))

#For testing
variant_names <- c("species_wgt_RAN_KBA")

#create general dat_files (settings file)
Settings <- "[Settings]
removal rule = 5
warp factor = 100 
edge removal = 0
add edge points = 0
memory save mode = 1"
fileConn <- file(paste0(zonation_path,"/","species_wgt_RAN_KBA.dat"))
writeLines(Settings, fileConn)
close(fileConn)

##create template spp_file (feature file). As is, outputs based on R project.
#Default template. Maybe I only have one feature group and can use get_variant and sppdata to filter features?
wgts <- as.numeric(spp_list$weight)
create_spp(
  filename = paste0(zonation_path,"/","species_wgt_RAN_KBA.spp"),
  weight = wgts, 
  alpha = 1,
  bqp = 1,
  bqp_p = 1,
  cellrem = 0.25,
  spp_file_dir = species_path, #this will change depending on which features; should match filename
  recursive = FALSE,
  spp_file_pattern = ".+\\.(tif|img)$",
  override_path = NULL
)

# Create new Zonation project
sensitive_site_ID <- create_zproject(name = "species_wgt_RAN_KBA_proj", 
                                     dir = zonation_path,
                                     variants = variant_names,
                                     spp_template_dir = species_path,
                                     spp_template_file = paste0(zonation_path,"/", "species_wgt_RAN_KBA.spp"),
                                     dat_template_file = paste0(zonation_path,"/", "species_wgt_RAN_KBA.dat"),
                                     overwrite = T)

#Place groups text file in appropriate place
#load groups file
groups_output <- spp_list$classGroupNum
groups_condition <- rep(-1, length(spp_list$X))
groups_retention <- rep(-1, length(spp_list$X))
retention_mode <- rep(-1, length(spp_list$X))
arb_kernal <- rep(-1, length(spp_list$X))
groups_file <- cbind(groups_output, groups_condition, groups_retention, retention_mode, arb_kernal)
write.table(groups_file, 
            file = paste0(zonation_path,"/","species_wgt_RAN_KBA_proj","/", "species_groups.txt"), 
            sep = " ", 
            col.names = F,
            row.names = F)
groups_path <- "species_groups.txt"

#Load project in to edit variants
sensitive_site_ID <- load_zproject(paste0(zonation_path,"/", "species_wgt_RAN_KBA_proj"))

#Edit variants 
variant <- get_variant(sensitive_site_ID, 1)
variant_spp_data <- sppdata(variant)
wgts <- as.numeric(spp_list$weight)
variant_spp_data$weight <- wgts
sppdata(variant) <- variant_spp_data
variant <- set_dat_param(variant, "use mask", 1)
variant <- set_dat_param(variant, "mask file", KBA_msk)
variant <- set_dat_param(variant, "use groups", 1)
variant <- set_dat_param(variant, "groups file", groups_path)
groups(variant) <- as.vector(groups_file[,1])
save_zvariant(variant, dir = paste0(zonation_path,"/", "species_wgt_RAN_KBA_proj"), overwrite = T, debug_msg = F)

#Edit bat file to add use threads
d <- read.delim(file = paste0(zonation_path, "/", "species_wgt_RAN_KBA_proj","/",variant@name,".bat"), sep = " ", header = F)
d$V13 <- apply(d, 1, function(x) paste(x, collapse = " "))
d <- d$V13
d <- paste(d, "--use-threads=4")
write.table(d, file = paste0(zonation_path, "/", "species_wgt_RAN_KBA_proj","/",variant@name,".bat"), sep = " ", col.names = F, row.names = F, append = F, quote = F)

create_sh_file <- function(x) {
  if (class(x) == "Zvariant") {
    bat_file <- x@bat.file
  } else {
    bat_file <- x
  }
  
  sh_file <- gsub("\\.bat", "\\.sh", bat_file)
  
  cmd_lines <- readLines(bat_file)
  new_cmd_lines <- c("#!/bin/sh")
  
  for (line in cmd_lines) {
    line <- gsub("call ", "", line)
    line <- gsub("\\.exe", "", line)
    new_cmd_lines <- c(new_cmd_lines, line)
  }
  
  file_con <- file(sh_file)
  writeLines(new_cmd_lines, file_con)
  close(file_con)
  Sys.chmod(sh_file)
  return(invisible(TRUE))
}
create_sh_file(variant)

#run
#run_bat(paste0(zonation_path, "/", "species_wgt_ABF_proj","/",variant@name,".bat"), exe = "zig4")
