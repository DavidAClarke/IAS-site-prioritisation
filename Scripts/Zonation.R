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

#Define variant names (maybe don't do random; add mask to species only variants)
variant_names <- c("species_ABF",
                   "species_CAZ",
                   "species_wgt_ABF",
                   "species_wgt_CAZ",
                   "species_wgt_ABF_KBA",
                   "species_wgt_CAZ_KBA",
                   "species_area_wgt_ABF",
                   "species_area_wgt_CAZ",
                   "species_area_wgt_ABF_KBA",
                   "species_area_wgt_CAZ_KBA")

                             ##Create new Zonation project ##
# See email from Jooma Lehtomaki (and script from his Japan SCP paper) on using different variants
# I may need different projects, one for each species group (therefore 4 projects)

#cell removal rules
cell_remove_rule <- c(1,2,5) #CAZ, ABF, Random (double check random number)

#Input data (raster) paths (species or species+area)
species_path <- file.path(Input_data_path, "species")
species_area_path <- file.path(Input_data_path, "species_area")

#### PERHAPS NESTED LOOPS CAN BE USED: 
# FIRST LOOP OVER CELL REMOVAL RULES:
  # ABF:
    # LOOP OVER DATASETS PER CELL REMOVAL RULE:
      # SPECIES, SPECIES+WEIGHTS, SPECIES+WEIGHTS+AREAS (WEIGHTED), SPECIES+WEIGHTS+AREAS (WEIGHTED) + MASK

#create dat_files (settings file)
ABF_settings <- "[Settings] 
removal rule = 2 
warp factor = 1 
edge removal = 1
add edge points = 1000"
fileConn <- file(file.path(project_path,"ABF_settings.dat"))
writeLines(Settings, fileConn)
close(fileConn)

CAZ_settings <- "[Settings] 
removal rule = 1 
warp factor = 1 
edge removal = 1
add edge points = 1000"
fileConn <- file(file.path(project_path,"CAZ_settings.dat"))
writeLines(Settings, fileConn)
close(fileConn)

#get list of dat files
settings <- list.files(project_path, pattern = ".dat$", full.names = T)

#group file could done using a function e.g. lapply over feature list and if name comes from 
#a it is KBA els if from b it is PA, etc. For species could be "class" names


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
#input_raster_dir <- "C:/Users/david.clarke1/Documents/Chapter_3/Chapter_3/SpatialData/Input_zonation/RL_IAS" #replace test_data

##create template spp_file (feature file). As is, outputs based on R project.
#Default template. Maybe I only have one feature group and can use get_variant and sppdata to filter features?
create_spp(
  filename = file.path(project_path,"species.spp"),
  weight = 1, #perhaps use equal weights as default
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
FinalReview_project <- create_zproject(name = "FinalReview_proj", 
                dir = project_path,
                variants = variant_names,
                spp_template_dir = input_raster_dir,
                spp_template_file = file.path(project_path,"species.spp"),
                dat_template_file = file.path("Zonation", "output.dat"))#,
                #weight = c(1, 1, 1, 2, 3, 2, 1)) #could already have this info as vector

#Edit variants (see Lehtomaki script for assistance)
#do some testing to make sure edits work

#try filtering features using sppdata

#Adding red list status, taxonomic class and weights to species
spp_list <- read.csv(file.path(project_path, "spp_list.csv")) #only contains shapefile species at moment
spp_list <- spp_list %>% rename(acceptedName = "x")
redlist_status <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "assessments.csv"))
redlist_taxonomy <- read.csv(file.path("SpeciesData", "redlist_species_data_18052021", "taxonomy.csv"))
TaxInfo <- read.csv(file.path("SpeciesData", "TaxInfo_assess_names.csv"))
spp_list <- spp_list %>%
  left_join(TaxInfo, by = "acceptedName") %>%
  left_join(redlist_status, by= "scientificName") %>%
  left_join(redlist_taxonomy, by = "scientificName") %>%
  dplyr::select(acceptedName, className, redlistCategory) %>%
  dplyr::distinct(acceptedName = acceptedName, .keep_all = T) %>%
  dplyr::mutate(weight = NA)

for(i in 1:length(spp_list$acceptedName)) {
  
  if(spp_list[i,]$redlistCategory == "Least Concern") {
    
    spp_list[i,]$weight <- 1
    
  }
  if(spp_list[i,]$redlistCategory == "Near Threatened") {
    
    spp_list[i,]$weight <- 2
    
  }
  if(spp_list[i,]$redlistCategory == "Vulnerable") {
  
    spp_list[i,]$weight <- 4
    
  }
  if(spp_list[i,]$redlistCategory == "Endangered") {
    
    spp_list[i,]$weight <- 6
    
  }
  if(spp_list[i,]$redlistCategory == "Critically Endangered") {
    
    spp_list[i,]$weight <- 8
    
  }
  if(spp_list[i,]$redlistCategory == "Data Deficient") {
    
    spp_list[i,]$weight <- 2
    
  }} 
  
spp_list <- spp_list %>%
  drop_na(weight) #removing extinct species
write.csv(spp_list,file.path(project_path, "spp_list.csv"))






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
