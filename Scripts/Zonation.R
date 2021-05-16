############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "Zonation.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"),
              local_path = file.path("Scripts", "Zonation.R"))

############################################################################################
#                                                                                          #
##          Zonation will not work through R if ANY file paths used have ANY spaces       ##
#                                                                                          #
############################################################################################

#Create path to a new project
#dir.create("Zonation")
project_path <- file.path("Zonation")

#Define variant names
#E.g. RL_all_ABF_eqw - All RL species, using the ABF rule and equal weights
#E.g. RL_all_ABF_w - All RL species, using the ABF rule and different weights
# variant_names <- c(
#   #All Red List species
#   "RL_all_ABF_eqw",
#   "RL_all_ABF_w",
#   "RL_all_CAZ_eqw",
#   "RL_all_CAZ_w",
#   "RL_all_RAN_eqw",
#   "RL_all_RAN_w",
#   
#   #IAS threatened Red list species
#   "RL_IAS_ABF_eqw",
#   "RL_IAS_ABF_w",
#   "RL_IAS_CAZ_eqw",
#   "RL_IAS_CAZ_w",
#   "RL_IAS_RAN_eqw",
#   "RL_IAS_RAN_w",
#   
#   #All SNES species
#   "SNES_all_ABF_eqw",
#   "SNES_all_ABF_w",
#   "SNES_all_CAZ_eqw",
#   "SNES_all_CAZ_w",
#   "SNES_all_RAN_eqw",
#   "SNES_all_RAN_w",
#   
#   #IAS threatened SNES species
#   "SNES_IAS_ABF_eqw",
#   "SNES_IAS_ABF_w",
#   "SNES_IAS_CAZ_eqw",
#   "SNES_IAS_CAZ_w",
#   "SNES_IAS_RAN_eqw",
#   "SNES_IAS_RAN_w")

#For final review just use one
variant_names <- c("RL_IAS_ABF_eqw")

                             ##Create new Zonation project ##
# See email from Jooma Lehtomaki (and script from his Japan SCP paper) on using different variants
# I may need different projects, one for each species group (therefore 4 projects)

#create dat_file (settings file)
#Test text creation (this works! Perhaps turn into function)
#try different (lower) warp factors for finer (but slower) solutions.
#Also try adding  'add edge points = n'. See bottom of page 103 of user manual
#Will also add a group setting at some point
Settings <- "[Settings] 
removal rule = 2 
warp factor = 1 
edge removal = 1"
fileConn <- file("Zonation/output.dat")
writeLines(Settings, fileConn)
close(fileConn)

#group file could done using a function e.g. lapply over feature list and if name comes from 
#a it is KBA els if from b it is PA, etc.


#There can be no spaces in tif files. This is a way to change file names in bulk
#I should just alter raster function to do this
# old_files <- list.files(file.path("SpatialData", "Input_zonation", "RL_IAS"), pattern = "*.tif", full.names = T)
# new_files <- gsub(" ", "_", old_files)
# file.copy(from = old_files, to = new_files)
# file.remove(old_files)
#example of removing specific files
#to_remove <- old_files[grepl(" ", old_files)]
#file.remove(to_remove)

#Doesn't like me using file.path() for this
#input_raster_dir <- "C:/Users/David Clarke.DESKTOP-NNNNVLL/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/test_data"
#input_raster_dir <- "C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/SpatialData/Input_zonation/RL_IAS" #replace test_data
input_raster_dir <- "C:/Users/david.clarke1/Documents/Chapter_3/Chapter_3/SpatialData/Input_zonation/RL_IAS" #replace test_data

##create template spp_file (feature file). As is, outputs based on R project.
#Default template. (note I will need two different files: 300 species and 70 species)
create_spp(
  filename = "Zonation/featurelist.spp",
  weight = 1, #perhaps use equal weights as default
  alpha = 1,
  bqp = 1,
  bqp_p = 1,
  cellrem = 0.25,
  spp_file_dir = input_raster_dir, #this will change depending on whether using all or only IAS threatened species
  recursive = FALSE,
  spp_file_pattern = ".+\\.(tif|img)$",
  override_path = NULL
)

# Create new Zonation project
FinalReview_project <- create_zproject(name = "FinalReview_proj", 
                dir = project_path,
                variants = variant_names,
                spp_template_dir = input_raster_dir,
                spp_template_file = file.path("Zonation","featurelist.spp"),
                dat_template_file = file.path("Zonation", "output.dat"))#,
                #weight = c(1, 1, 1, 2, 3, 2, 1)) #could already have this info as vector




                                     ## Running Zonation project ##
# I think I just run batch files using run_bat()
#example. It doesn't like me using file.path() for this. Note: this took 7.14 minutes
#run_bat("C:/Users/dcla0008/Dropbox/PhD/Thesis/Data/Chapter_3/Zonation/FinalReview_proj/RL_IAS_ABF_eqw.bat", exe = "zig4")
run_bat("C:/Users/david.clarke1/Documents/Chapter_3/Chapter_3/Zonation/FinalReview_proj/RL_IAS_ABF_eqw.bat", exe = "zig4")

#Another little tid bit I just figured out which may be helpful. 
#Especially since file.path() doesn't work
wd <- getwd()
path <- paste0(wd, "/RL_IAS_ABF_eqw.bat") #could even loop e.g. paste0(wd, i)
#e.g. run_bat(path, exe = "zig4")


                                     ## Working with Zonation results ##
#Look into automatic post-processing:
#landscape comparison option page 138
path <- file.path("Zonation", "FinalReview_proj")
tutorial <- load_zproject(path)

variant <- get_variant(tutorial, "RL_IAS_ABF_eqw")
var_results <- results(variant)
feature_curves <- curves(var_results)
plot(feature_curves) #I should make my own plots instead of using the packages limited plotting functionality



###########################  NOT ACTUALLY REQUIRED (but useful to know) ###############
#rename .spp to .txt (maybe requires function: loop through variant names)
#example for variant 1 (requires RSAGA package)
#spp1 <- file.path("Zonation", "test", "01_variant", "01_variant.spp")
#txt1 <- set.file.extension(spp, extension = ".txt")
#file.rename(spp1,txt1)

#spp2 <- file.path("Zonation", "test", "02_variant", "02_variant.spp")
#txt2 <- set.file.extension(spp, extension = ".txt")
#file.rename(spp2,txt2)

#spp3 <- file.path("Zonation", "test", "03_variant", "03_variant.spp")
#txt3 <- set.file.extension(spp, extension = ".txt")
#file.rename(spp3,txt3)
#########################################################################################
