############################# Uploading and downloading scripts to dropbox #####################
drop_upload(file = file.path("Scripts", "Functions.R"), 
            path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"))

drop_download(path = file.path("PhD", "Thesis", "Data", "Chapter_3", "Scripts"),
              local_path = file.path("Scripts", "Functions.R"))

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
    OccS <- occ_data(scientificName = i, hasCoordinate = T, hasGeospatialIssue = F, limit = 200000)
    OccS <- as_tibble(OccS$data)
    OccS <- OccS %>% dplyr::select(scientificName, decimalLatitude, decimalLongitude,coordinateUncertaintyInMeters,geodeticDatum,basisOfRecord, country)
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
  featureList <- gsub(" | /", "_", featureList)
  
  #dir.create(file.path("SpatialData", "Input_zonation", target_path))
  
  res.out <- lapply(featureList[1:length(featureList)], function(i) {
    
    print(i)
    
    #sci_name <- shp %>% filter(scientificName) == i
    temp <- shp[,feature_name] == i
    temp_num <- as.numeric(paste0(which(temp[, feature_name] == T)))
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
    
    writeRaster(rst_sci_name, filename = file.path("SpatialData", "Input_zonation", target_path, paste0(i,".tif")), format = "GTiff", overwrite = T)
    
    #The following is only required when working from virtual machine
    Sys.sleep(10)
    if(is_empty(drop_acc()) == F){
    drop_upload(file = file.path("SpatialData", "Input_zonation", target_path, paste0(i, ".tif")), 
                path = file.path("PhD", "Thesis", "Data", "Chapter_3", "SpatialData", "Input_zonation", target_path))
    } else if(is_empty(drop_acc()) == T) {
      print("No dropbox account information identified")
    }
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
    Sys.sleep(15)
    
    options(timeout = 200)
    download.file(i, 
                  destfile = file.path("Data", target_path, paste0(which(urls == i, arr.ind = T), "_scrapedpage.html")),
                  quiet = T 
                  )
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

Threat_info <- function(urls, species, target_folder) {
  
  scientificName <- as.vector(species)
  ThreatInfo <- vector(mode = "character", length = length(urls))
  
  dir.create(file.path("Data", target_folder))
  
  res.out <- lapply(urls[1:length(urls)], function(i = urls) {  
    
    print(i)
    
    Sys.sleep(10)
    
    #It helps to have chrome as the default browser
    options(timeout = 300)
    download.file(i, 
                  destfile = file.path("Data", target_folder, paste0(which(urls == i, arr.ind = T), "_scrapedpage.html")), 
                  quiet = T,
                  extra = options(download.file.method = "libcurl")
                  )
    spec_sprat <- read_html(file.path("Data", target_folder, paste0(which(urls == i, arr.ind = T), "_scrapedpage.html")))
    
    ThreatInfo <- spec_sprat %>%
           html_element("#bodyContent") %>%
           html_text()
    })
  
  res.out <- do.call(rbind,res.out)
  res.out <- data.frame(scientificName = scientificName, 
                        ThreatInfo = res.out)
  return(as_tibble(res.out))
}


  
