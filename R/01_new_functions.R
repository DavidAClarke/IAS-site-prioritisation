################################################################################
## Script name: 01_functions.R
## Author: David Clarke
## Copyright (c) David Clarke, 2024
## Email: david_anthony_clarke@hotmail.com
################################################################################

## Priority rank map----
rank_plot <- function(rankmap){
  
  ras_st <- stars::st_as_stars(rankmap)
  ras_sf <- st_as_sf(ras_st)
  fill_ras_sf <- st_drop_geometry(ras_sf)
  ggplot()+
    geom_sf(data = ras_sf, aes(fill=fill_ras_sf[,1]), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = z_colors_spectral$colors,
                         values = z_colors_spectral$values,
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

## Priority rank difference map----
rank_diff <- function(rankmap1, rankmap2){
  
  species_only_diff <- rankmap1 - rankmap2
  names(species_only_diff) <- "rankmap"
  
  coolwarm_hcl <- colorspace::diverging_hcl(11,h = c(250, 10), c = 100, 
                                            l = c(37, 88), power = c(0.7, 1.7))
  species_only_diff_sf <- species_only_diff %>%
    stars::st_as_stars() %>%
    st_as_sf()
  
  ggplot()+
    geom_sf(data = species_only_diff_sf, aes(fill=rankmap), 
            color=NA, 
            show.legend = T) + 
    scale_fill_gradientn(colours = rev(coolwarm_hcl),
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

## Calculate Jaccard similarities----
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

## Raster correlation matrix----
ras_cor <- function(ras_stack){
  
  rcors <- sapply(full_rank_stack, function(x) 
            sapply(full_rank_stack, function(y) 
              layerCor(c(x,y), "cor", asSample = F)))

  rcors <- rcors[seq(1, length(rcors), 3)] #get every 3rd element (coefs)
  mycors <- c()

  for(i in 1:length(rcors)){
  
    v <- rcors[[i]][2,1]
    mycors <- c(mycors,v)
  
  }

  mymat <- round(matrix(data = mycors, 
                        nrow = sqrt(length(mycors)), 
                        ncol = sqrt(length(mycors)), 
                        byrow = T),2)
  
  return(mymat)

}

## Structural similarity index----
ssim <- function(ras_stack){
  
  ssim_mat <- matrix(nrow = nlyr(ras_stack), ncol = nlyr(ras_stack))
  sim_mat <- matrix(nrow = nlyr(ras_stack), ncol = nlyr(ras_stack))
  siv_mat <- matrix(nrow = nlyr(ras_stack), ncol = nlyr(ras_stack))
  sip_mat <- matrix(nrow = nlyr(ras_stack), ncol = nlyr(ras_stack))
  
  my_combs <- combn(nlyr(ras_stack), 2)
  
  for(k in 1:ncol(my_combs)){
    
    my_ssim <- ssim_raster(ras_stack[[my_combs[1,k]]], ras_stack[[my_combs[2,k]]], global = F)
    
    varnames(my_ssim) <- paste0(names(ras_stack[[my_combs[1,k]]]),"-",
                                      names(ras_stack[[my_combs[2,k]]]))
    
    my_path <- here(dirname(here()), "data", "ssim")
    
    if(file.exists(my_path)){
      
      writeRaster(my_ssim, filename = here(my_path, paste0(varnames(my_ssim), ".tif")))
      
    } else {
      
      dir.create(my_path)
      writeRaster(my_ssim, filename = here(my_path, paste0(varnames(my_ssim), ".tif")))
      
    }
    
    ssim_mat[my_combs[2,k],my_combs[1,k]] <- as.numeric(global(my_ssim[[1]], "mean", na.rm = T))
    sim_mat[my_combs[2,k],my_combs[1,k]] <- as.numeric(global(my_ssim[[2]], "mean", na.rm = T))
    siv_mat[my_combs[2,k],my_combs[1,k]] <- as.numeric(global(my_ssim[[3]], "mean", na.rm = T))
    sip_mat[my_combs[2,k],my_combs[1,k]] <- as.numeric(global(my_ssim[[4]], "mean", na.rm = T))
        
  }
  
 diag(ssim_mat) <- diag(sim_mat) <- diag(siv_mat) <- diag(sip_mat) <- 1
 rownames(ssim_mat) <- rownames(sim_mat) <- rownames(siv_mat) <- rownames(sip_mat) <- names(ras_stack)
 colnames(ssim_mat) <- colnames(sim_mat) <- colnames(siv_mat) <- colnames(sip_mat) <- names(ras_stack)
  
 my_list <- list(ssim_mat = ssim_mat, 
                 sim_mat = sim_mat, 
                 siv_mat = siv_mat, 
                 sip_mat = sip_mat)
 
 return(my_list)
  
}

## Calculate Jaccards----
jaccard <- function(x, y, x.min=0.0, x.max=1.0, y.min=0.0, y.max=1.0,
                    warn.uneven=FALSE, limit.tolerance=4,
                    disable.checks=FALSE) {
  
  if (!disable.checks) {
    # Check the input values
    x.min.value <- round(raster::cellStats(x, stat="min"), limit.tolerance)
    x.max.value <- round(raster::cellStats(x, stat="max"), limit.tolerance)
    y.min.value <- round(raster::cellStats(y, stat="min"), limit.tolerance)
    y.max.value <- round(raster::cellStats(y, stat="max"), limit.tolerance)
    
    if (x.min < x.min.value) {
      stop(paste0("Minimum threshold value for x ("), x.min, ") smaller than
            the computed minimum value in x (", x.min.value, ")")
    }
    if (x.max > x.max.value) {
      stop(paste0("Maximum threshold value for x ("), x.max, ") smaller than
            the computed maximum value in x (", x.max.value, ")")
    }
    if (x.min >= x.max) {
      stop(paste0("Minimum threshold value for x ("), x.min, ") smaller than
             maximum threshold value for x (", x.max, ")")
    }
    if (y.min < y.min.value) {
      stop(paste0("Minimum threshold value for y ("), y.min, ") smaller than
            the computed minimum value in y (", y.min.value, ")")
    }
    if (y.max > y.max.value) {
      stop(paste0("Maximum threshold value for y ("), y.max, ") smaller than
            the computed maximum value in y (", y.max.value, ")")
    }
    if (y.min >= y.max) {
      stop(paste0("Minimum threshold value for y ("), y.min, ") smaller than
             maximum threshold value for y (", y.max, ")")
    }
    
    # Comparisons using just the defaults is probably not feasible
    if (x.min == 0.0 & x.max == 1.0 & y.min == 0.0 & y.max == 1.0) {
      warning("Using all the defaults for x and y ranges")
    }
  } else {
    message("Input limit checks disabled")
  }
  
  # [fixme] - using cellStats(X, "sum") should be safe as we're dealing with
  # binary 0/1 rasters. count() would be preferable, but apparently raster
  # (>= 2.2 at least) doesn't support it anymore.
  
  # Get the values according to the limits provided
  x.bin <- (x >= x.min & x <=x.max)
  y.bin <- (y >= y.min & y <=y.max)
  
  if (warn.uneven) {
    x.size <- raster::cellStats(x.bin, "sum")
    y.size <- raster::cellStats(y.bin, "sum")
    # Sort from smaller to larger
    sizes <- sort(c(x.size, y.size))
    if (sizes[2] / sizes[1] > 20) {
      warning("The extents of raster values above the threshhold differ more",
              "than 20-fold: Jaccard coefficient may not be informative.")
    }
  }
  
  # Calculate the intersection of the two rasters, this is given by adding
  # the binary rasters together -> 2 indicates intersection
  combination <- x.bin + y.bin
  intersection <- combination == 2
  
  # Union is all the area covered by the both rasters
  union <- combination >= 1
  
  return(raster::cellStats(intersection, "sum") / raster::cellStats(union, "sum"))
}
## Susceptible site prep
susc_site_prep <- function(species_name, ras_stack){
  
  sp <- gsub(" ", ".", species_name)

  bin <- rast(raster::raster(here(regional_model_path, 
                             sp, 
                             "proj_regional", 
                             "individual_projections", 
                             paste0(sp,
                                    "_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.gri"))))
  
  bin2 <- bin
  bin2[bin2 != 0] <- 1
  
  # Re-class zeros as NA
  bin[bin == 0] <- NA
  
  bin_list <- list()
  
  for(i in 1:nlyr(ras_stack)){
  
    bin_list[[i]] <- resample(bin, ras_stack[[i]], method = "near") 
    names(bin_list[[i]]) <- names(ras_stack[[i]])
    
  }
  
  return(list(bin_list, bin2))

}

## Get site sensitivity values----
get_msk_vals <- function(rank_raster, mask_file) {
  
  temp_r <- mask(rank_raster, mask_file)
  temp_r_values <- na.omit(values(temp_r))
  
}

## Proportion difference (KBA vs no KBA) in number of top sensitive sites----
#vals = output from get_msk_vals, sens = site sensitivity value
prop_diff <- function(vals, sens) {
  
  vals_bin <- vals
  vals_bin <- ifelse(vals_bin >= sens, 1,0)
  sum(vals_bin)/length(vals_bin)*100
  
}

## Get proportion differences----
multi_props <- function(vals, props){
  
  diffs <- c()
  
  for(i in props[1:length(props)]){
    
    d <- prop_diff(vals, i)
    diffs <- c(diffs,d)
    
    
  }
  names(diffs) <- c("1.00", "0.98", "0.95", "0.90", "0.75", "0.50", "0.25", 
                    "0.00")
  return(diffs)
}
