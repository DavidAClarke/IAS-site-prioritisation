################################################################################
## Script name: 01_functions.R
################################################################################

xy_match <- function(initial, pre, Layer) {
  pre <- initial
  
  coordinates(pre) <- ~ decimalLongitude + decimalLatitude
  crs(pre) <- crs(Layer)
  
  final <- cbind(initial, over(pre, Layer))
  return(final)
}

#Removing points with no environmental data
rm_occs <- function(env, pts) {
  
  oc <- extract(env, pts[,2:3])
  oc <- as.data.frame(oc)
  pts <- cbind(pts, oc)
  pts <- na.omit(pts)
  pts <- pts %>% dplyr::select(scientificName, Longitude, Latitude)
  return(pts)
  
}

# Create a shell (.sh) file
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

#Get scenario variants
scenario_variants <- function(project) {
  
  proj <- load_zproject(file.path("zonation", project))
  proj_var <- get_variant(proj, 1)
  
}


#Create rank raster plots
rank_plot <- function(rank_raster) {
  
  ras_st <- st_as_stars(rank_raster)
  ras_sf <- st_as_sf(ras_st)
  fill_ras_sf <- st_drop_geometry(ras_sf)
  ggplot()+
    geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = leg$colors,
                         values = leg$values,
                         name = "Site\nsensitivity",
                         breaks = c(0.0, 0.2, 0.4, 0.6, 0.8,1)) +
    theme_bw() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank())
}

rank_diff <- function(rank_one, rank_two){
  
  species_only_diff <- rank_one - rank_two
  coolwarm_hcl <- colorspace::diverging_hcl(11,h = c(250, 10), c = 100, 
                                            l = c(37, 88), power = c(0.7, 1.7))
  species_only_diff_sf <- species_only_diff %>%
    st_as_stars() %>%
    st_as_sf()
  
  ggplot()+
    geom_sf(data = species_only_diff_sf, aes(fill=layer), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = rev(coolwarm_hcl),
                         #values = leg$values,
                         name = "Sensitivity\ndifference",
                         breaks = c(-0.5, 0.0, 0.5)) +
    theme_bw() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank())
}

#Priority plots - overlap of sensitive and susceptible sites 
Priority_plot <- function(sdm, rank_raster, sdm_col, sdm_brd_col) {
  
  sdm_r <- sdm
  sdm_r<- rasterToPolygons(sdm_r, dissolve = T)
  sdm_cr <- crop(sdm, sdm_r)
  rank_raster_cr <- crop(rank_raster, sdm_cr)
  rank_raster_cr_st <- st_as_stars(rank_raster_cr)
  rank_raster_cr_sf <- st_as_sf(rank_raster_cr_st)
  sdm_sf <- st_as_sf(sdm_r)
  sdm_sf_bb <- st_as_sfc(st_bbox(sdm_sf))
  fill_rank_raster_cr_sf <- st_drop_geometry(rank_raster_cr_sf)
  g <- ggplot()+
    geom_sf(data = rank_raster_cr_sf, aes(fill=fill_rank_raster_cr_sf[,1]), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = leg$colors,
                         values = leg$values,
                         name = "Site\nsensitivity",
                         breaks = c(0.0001, 0.25, 0.5, 0.75, 0.9999),
                         labels = as.character(c(0.0, 0.25, 0.5, 0.75, 1))) +
    geom_sf(data = sdm_sf, show.legend = F,fill=alpha(sdm_col,0.3)) +
    theme_bw() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          #legend.position = "none",
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.border = element_rect(colour = sdm_brd_col, size = 1.2))
  
  return(g)
}

#Get cell values for Kolmogorov-smirnoff tests
get_msk_vals <- function(rank_raster, mask_file) {
  
  temp_r <- mask(rank_raster, mask_file)
  temp_r_values <- na.omit(values(temp_r))
  
}

#Create histograms
mask_hist <- function(rank_raster, mask_layer) {
  
  vals <- get_msk_vals(rank_raster, mask_layer)
  b <- seq(0.0, 1.0, by = 0.05)
  h <- hist(vals, breaks = b, plot = F)
  up_lim <- max(h$counts)
  up_lim <- round(up_lim + 50, -2)
  vals_df <- data.frame(data = vals)
  
  val_hist <- ggplot(vals_df, aes(x = vals)) + 
    geom_histogram(colour = "white",
                   binwidth = 0.05,
                   boundary = 0) + 
    scale_x_continuous(breaks = seq(0, 1, 0.25)) + 
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,up_lim)) +
    xlab("Site sensitivity") + 
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 12)) +
    geom_vline(xintercept = median(vals, #could also use mean
                                   na.rm = T), colour = "red")
  
  return(val_hist)
  
}

#Proportion difference (KBA vs no KBA) in number of top sensitive sites
#vals = output from get_msk_vals, sens = site sensitivity value
prop_diff <- function(vals, sens) {
  
  vals_bin <- vals
  vals_bin <- ifelse(vals_bin >= sens, 1,0)
  sum(vals_bin)/length(vals_bin)*100
  
}

#Calculate Jaccard similarities
calculate_jaccards <- function(rank_stack, x.min, x.max, y.min, y.max, 
                               variant_names) {
  
  jaccards <- matrix(nrow = nlayers(rank_stack), ncol = nlayers(rank_stack))
  for (i in 1:nrow(jaccards)) {
    for (j in 1:ncol(jaccards)) {
      if (i == j) {
        jaccards[i, j] <- 1
      }
      else {
        if (is.na(jaccards[j, i])) {
          message(paste0("Calculating Jaccard index between ", 
                         names(rank_stack[[i]]), 
                         " and ", names(rank_stack[[j]])))
          jaccards[i, j] <- jaccard(rank_stack[[i]], rank_stack[[j]], 
                                    x.min = x.min, x.max = x.max, 
                                    y.min = y.min, y.max = y.max)
        }
        else {
          jaccards[i, j] <- NA
        }
      }
    }
  }
  jaccards <- as.data.frame(jaccards)
  colnames(jaccards) <- variant_names
  rownames(jaccards) <- variant_names
  return(jaccards)
}

#Taking a little off values of 1
squeeze <- function(pvector){
  for(i in 1:length(pvector)) {
    if(!is.na(pvector[i]) == T & pvector[i] == 1){
      pvector[i] <- pvector[i] - 0.000001
    }
  }
  return(pvector)
}

#Create binary layers based on given site sensitivity
spat_priority_dist <- function(df, n_col){
  for(i in 1:n_col){
    df[,i] <- ifelse(df[,i] >= 0.98, 1,0)}
  return(df)
}

#Create predicted IAS distribution plots
IAS_plot <- function(species){
  species <- gsub(" ", ".", species)
  r <- raster(file.path(regional_model_path, species, "proj_regional", 
                        "individual_projections", 
                        paste0(species,"_EMcaByTSS_mergedAlgo_mergedRun_mergedData.gri")))
  r <- r/1000 #back to a 0-1 scale
  r_st <- st_as_stars(r)
  r_sf <- st_as_sf(r_st)
  r_plot <- ggplot()+
    geom_sf(data = r_sf, aes(fill=layer), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = brewer.pal('YlGnBu', n=9),
                         name = "Probability") +
    theme_bw() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank())
  return(r_plot)
}


#Get proportion differences
multi_props <- function(vals, props){
  
  for(i in props[1:length(props)]){
    
    d <- prop_diff(vals, i)
    diffs <- c(diffs,d)
    
    
  }
  names(diffs) <- c("1.00", "0.98", "0.95", "0.90", "0.75", "0.50", "0.25", 
                    "0.00")
  return(diffs)
}

#Multi-panel figure
Figure <- function(species, rank_raster, binary) {
  
  species <- gsub(" ", ".", species)
  p <- IAS_plot(species)
  hist <- mask_hist(rank_raster, binary)
  priority <- Priority_plot(binary, rank_raster, sdm_col = "black", 
                            sdm_brd_col = NA)
  p2 <- ggarrange(p, priority, nrow = 2, ncol = 1)
  p3 <- ggarrange(p2, hist, nrow = 1, ncol = 2)
  return(p3)
  
}

#Biodiversity feature performance plots
performance_plot <- function(project_variant, feature_groups, brewer_pal) {
  
  if(length(feature_groups) <= 8){
    
    groupnames(project_variant) <- feature_groups
    lost.levels <- seq(0,1, by = 0.05)
    results.caz <- results(project_variant)
    perf <- performance(results.caz, lost.levels,melted = TRUE, groups = T)
    perf <- na.omit(perf)
    perf_wmean <- perf %>% dplyr::filter(str_detect(feature, "w.mean"))
    perf_wmean$feature <- as.factor(perf_wmean$feature)
    perf_min <- perf %>% dplyr::filter(str_detect(feature, "min."))
    perf_min$feature <- as.factor(perf_min$feature)
    feature.groups <- c("Amphibian", "Bird", "Fish", "Fungi", "Invertebrate", 
                        "Mammal", "Plant", "Reptile")
    feat.cols <- RColorBrewer::brewer.pal(length(feature_groups), brewer_pal)
    names(feat.cols) <- levels(perf_wmean$feature)
    colScale <- scale_colour_manual(name = "Feature",
                                    values = feat.cols,
                                    labels = feature.groups)
    ggplot(data = perf, mapping = aes(x = pr.lost, y = perf.levels, 
                                      colour = feature)) +
      geom_line(data = perf_wmean, size = 1) +
      colScale +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 12), 
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            strip.text.x = element_text(size = 12)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      ylab("Distribution remaining") +
      xlab("Proportion of landscape lost")
  } else if(length(feature_groups) > 8) {
    
    groupnames(project_variant) <- feature_groups
    lost.levels <- seq(0,1, by = 0.05)
    results_area.caz <- results(project_variant)
    perf_area <- performance(results_area.caz, lost.levels,melted = TRUE, 
                             groups = T)
    perf_area <- na.omit(perf_area)
    perf_area_wmean <- perf_area %>% dplyr::filter(str_detect(feature, "w.mean"))
    perf_area_wmean$feature <- as.factor(perf_area_wmean$feature)
    #perf_min <- perf %>% dplyr::filter(str_detect(feature, "min."))
    #perf_min$feature <- as.factor(perf_min$feature)
    feature_area.groups <- c("Amphibian", "Bird", "Community", "Ecosystem", 
                             "Fish", "Fungi", "Invertebrate", 
                             "Mammal", "Plant", "Ramsar", "Reptile", "Upstream")
    feat_area.cols <- RColorBrewer::brewer.pal(length(feature_groups), brewer_pal)
    names(feat_area.cols) <- levels(perf_area_wmean$feature)
    colScale <- scale_colour_manual(name = "Feature",
                                    values = feat_area.cols,
                                    labels = feature_area.groups)
    ggplot(data = perf_area, mapping = aes(x = pr.lost, y = perf.levels, 
                                           colour = feature)) +
      geom_line(data = perf_area_wmean, size = 1) +
      colScale +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 12), 
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            strip.text.x = element_text(size = 12)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      ylab("Distribution remaining") +
      xlab("Proportion of landscape lost")
  }
  
}
