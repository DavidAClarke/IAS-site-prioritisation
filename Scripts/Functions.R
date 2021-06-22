############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "Functions.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts", "Functions.R"),
              local_path = file.path("Scripts", "Functions.R"), overwrite = T)

###################################################################################################
#~#~# Function loads / installing packages required libraries----
loadLibrary <- function(packageList){
  
  packLoaded <- list.files(.libPaths())
  
  lapply(packageList,function(x){
    
    if(!(x %in% packLoaded)){install.packages(x)}
    
    require(x,character.only = TRUE)
    
  })
  
}

#~#~# Function imports GBIF occurrence data----

Get_gbif_Occ <- function(sppList) {
  
  res.out <- lapply(sppList[1:length(sppList)], function(i) {
    
    print(i)
    
    # Getting GBIF entries with this name: 
    all_occ <- name_backbone(i)
    
    # Getting the species key, unique for the species and its synonimous:
    key_no <- all_occ$speciesKey
    
    # Get the number of occurrences in gbif:
    num_occ <- occ_count(taxonKey=key_no, georeferenced=TRUE)
    
    if (num_occ<=200000) {
      # Search for GBIF occurrences according to 'key_no'
      OccS <- occ_data(taxonKey = key_no,
                       country = "AU",
                       hasCoordinate = T, 
                       hasGeospatialIssue = F, 
                       limit = num_occ)
      if(is_empty(OccS[["data"]]) == FALSE) {
      OccS <- as_tibble(OccS$data)
      OccS <- OccS %>% 
        dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION" | 
                      basisOfRecord == "MACHINE_OBSERVATION" |
                      basisOfRecord == "OBSERVATION" ) %>%
        dplyr::select(scientificName,
                      species,
                      decimalLatitude, 
                      decimalLongitude,
                      coordinateUncertaintyInMeters,
                      geodeticDatum,
                      basisOfRecord, 
                      country)
      } else {NULL}
    } #else {next}
    
  })
  res.out <- do.call(rbind,res.out)
  
  return(res.out)
}

#~#~# Function converting NA to factor level----
NA2fctlvl <- function(df) {
  
  for(j in 1:ncol(df)) {
    
    levels(df[[j]]) <- c(levels(df[[j]]),"NA")
    df[[j]][is.na(df[[j]])] <- "NA"
    
  }
  return(df)
}

#Getting higher level taxonomic information
GetTax <- function(sppList) {
  
  res.out <- lapply(sppList[1:length(sppList)], function(i) {
    
    print(i)
    
    Ords <- gbif_name_usage(name = i, limit = 1)
    
    unlist(lapply(Ords$results, function(e) {
      unlist(e[c("canonicalName", "order", "class")])
    }))
  })
  res.out <- as_tibble(do.call(rbind,res.out))
  
  return(res.out)
}

#Creating rasters for features (in progress)
feature_rst <- function(shp, target_path, feature_name) { 
  
  featureList <- shp[[feature_name]]
  featureList <- unique(featureList) 
  
  #dir.create(file.path("SpatialData", "Input_zonation", target_path))
  
  res.out <- lapply(featureList[1:length(featureList)], function(i) {
    
    print(i)
    
    #sci_name <- shp %>% filter(scientificName) == i
    temp <- shp[,col_name] == i
    temp_num <- as.numeric(paste0(which(temp[, col_name] == T)))
    sci_name <- shp[temp_num,]
    
    new_bb <- c(-1863742, -4840841, 2088051, -1168811) #double check these numbers. Really should be the same as ext when creating raster
    names(new_bb) <- c("xmin", "ymin", "xmax", "ymax")
    attr(new_bb, "class") = "bbox"
    attr(st_geometry(sci_name), "bbox") = new_bb
    
    rst_template <- raster(resolution = 1000,
                           crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", 
                           xmn = -1863742,
                           xmx = 2088051,
                           ymn = -4840841,
                           ymx = -1168811) 
    
    rst_sci_name <- rasterize(st_zm(sci_name), rst_template, field = 1)
    
    writeRaster(rst_sci_name, filename = file.path("SpatialData", "Input_zonation", target_path, paste0(i,".tif")), format = "GTiff")
    
    #drop_upload(file = file.path("SpatialData", "Input_zonation", target_path, paste0(i, ".tif")), 
    #            path = file.path("PhD", "Thesis", "Data", "Chapter 3", "SpatialData", "Input_zonation", target_path))
  })
}


# Obtaining information on whether a given SPRAT web page contains information on "Threats"
# The target_path argument is the folder name for storing the downloaded html documents
Threat_present <- function(urls,species,target_path) {
  
  scientificNames <- as.vector(species)
  ThreatPresent <- vector(mode = "numeric", length = length(urls))
  Profiles <- as.vector(urls)
  
  dir.create(file.path("Data", target_path))
  
  res.out <- lapply(urls[1:length(urls)], function(i = urls) {  
    
    print(i)
    
    #Avoiding a HTTP error 429 (too many requests)
    Sys.sleep(10)
    
    
    download.file(i, destfile = file.path("Data", target_path, paste0(which(urls == i, arr.ind = T), "_scrapedpage.html")), quiet = T)
    spec_sprat <- read_html(file.path("Data", target_path, paste0(which(urls == i, arr.ind = T), "_scrapedpage.html")))
    
    
    ThreatInfo <- spec_sprat %>%
    html_elements(xpath = "//p[preceding::h2[1][contains(.,'Threats')]]") %>%
    html_text()
    
    ThreatPresent <- ifelse(is_empty(ThreatInfo) == F, 1,0)
    
    
    })
  
  res.out <- do.call(rbind,res.out)
  res.out <- data.frame(scientificNames = scientificNames, 
                        ThreatPresent = res.out,
                        Profile = Profiles)
  return(as_tibble(res.out))
}


# Getting threat information. Similar function to above used on those species that return 1
# This needs to be changed to get whole webpage e.g.
#ThreatInfo <- spec_sprat %>%
#     html_element("#bodyContent") %>%
#     html_text()
#Then keyword search on this result i.e. on the whole web page
#Do I save each result to a text file, that I then read in and keyword search later?

Threat_info <- function(urls,species) {
  
  scientificNames <- as.vector(species)
  ThreatInfo <- vector(mode = "character", length = length(urls))
  
  res.out <- lapply(urls[1:length(urls)], function(i = urls) {  
    
    Sys.sleep(10)
    
    download.file(i, destfile = file.path("Data", "rnd2_websites", paste0(which(urls == i, arr.ind = T), "_scrapedpage.html")), quiet = T)
    spec_sprat <- read_html(file.path("Data", "rnd2_websites", paste0(which(urls == i, arr.ind = T), "_scrapedpage.html")))
    
    ThreatInfo <- spec_sprat %>%
           html_element("#bodyContent") %>%
           html_text()
    })
  
  res.out <- do.call(rbind,res.out)
  res.out <- data.frame(scientificNames = scientificNames, 
                        ThreatInfo = res.out)
  return(as_tibble(res.out))
}

#Calculating proportion of total range that is in Australia
#Would use this after filtering out species that have no Australian range 

calc_area <- function(feature_shp, feature_name, country_shp) {
  
  feature_list <- feature_shp[[feature_name]]
  feature_list <- unique(feature_list)
  
  res.out <- lapply(feature_list[1:length(feature_list)], function(i) {
  
    print(i)
  
    spec <- feature_shp %>% filter(BINOMIAL == i)
    spec_dist <- st_intersection(st_buffer(country_shp, dist = 0), st_buffer(spec, dist = 0))
    spec_area <- spec_dist %>% 
      mutate(Area = st_area(spec_dist)) %>%
      st_drop_geometry()
    spec_area <- spec_area %>%
      group_by(BINOMIAL, NAME_0) %>%
      summarise(Country_area = sum(Area)) %>%
      mutate(Perc_Total = Country_area/sum(Country_area)*100) %>%
      rename(scientificName = "BINOMIAL") 
  })

  res.out <- do.call(rbind,res.out)
  return(as_tibble(res.out))
}

#~#~# Function matches co-ordinates to polygons----
#Notes: initial = original data frame, pre = altered version of initial with new name, Layer = shapefile to be matched


xy_match <- function(initial, pre, Layer) {
  pre <- initial
  
  coordinates(pre) <- ~ Long + Lat
  proj4string(pre) <- proj4string(Layer)
  
  final <- cbind(initial, over(pre, Layer))
  return(final)
}
