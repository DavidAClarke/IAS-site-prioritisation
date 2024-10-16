##IAS SDMs##
regional_model_path <- here(dirname(here()), "data", "IAS_regional")

spp_list <- c("Apis mellifera",  "Monomorium floricola",
              "Monomorium destructor","Linepithema humile", "Vespula vulgaris",
              "Bombus terrestris", "Heteronychus arator",
              "Digitonthophagus gazella", "Pheidole megacephala",
              "Vespula germanica","Tetramorium bicarinatum",
              "Paratrechina longicornis")

## Prepare for priority sites results
susceptible_site_prep <- lapply(spp_list[1:length(spp_list)], function(i){
  
  susc_site_prep(i, full_rank_stack)
  
})

## Overlap among all IAS
ias_stack <- rast()

for(i in 1:length(susceptible_site_prep)){
  
  r <- susceptible_site_prep[[i]][[2]]
  ias_stack <- c(ias_stack, r)
  
}

ias_rich <- sum(ias_stack)
ias_rich_one <- ias_rich
ias_rich_one[ias_rich_one != 1] <- NA

ias_sum_sp <- resample(ias_rich, full_rank_stack[[2]], method = "near")
ias_sum_sp_ar <- resample(ias_rich, full_rank_stack[[15]], method = "near")

ras_sf <- st_as_stars(ias_sum_sp_ar) %>%
        st_as_sf()

fill_ras_sf <- st_drop_geometry(ras_sf)
        
cl <- colorRampPalette(c("#e69b99","#2c6184"))

ias_map <- ggplot()+
  geom_sf(data = ras_sf, 
          aes(fill=fill_ras_sf[,1]), 
          color=NA, 
          show.legend = T) +
  scale_fill_gradientn(colours = cl(12),
                       name = "Alien\nrichness",
                       breaks = seq(0,11)) +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

non.na <- as.numeric(global(ias_sum_sp, fun = "notNA")) #number of non NA cells
no.cells <- c()
perc.total <- c()

for(i in 0:11){
  
  c <- length(which(values(ias_sum_sp == i)))
  no.cells <- c(no.cells, c)
  p <- (c/non.na)*100
  perc.total <- c(perc.total, p)
  
}

richness <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
totaltable <- data.frame("Richness" = richness, 
                         "Number of cells" = no.cells, 
                         "Percent total" = round(perc.total,2))

pretty_table <- formattable::formattable(totaltable,
            align = c("c", "c", "c"),
            list(Richness = formattable::color_tile("#e69b99","#2c6184")))

ggplot(totaltable, aes(log(no.cells), log(perc.total))) +
  #geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_x_reverse(breaks = log(no.cells), labels = richness) +
  xlab("Alien richness") +
  ylab("Percent of total area (log)")

#SDM evaluations
eval_scores <- read.table(here(regional_model_path, 
                               "Table_IAS_regional_accuracy.txt"), 
                          header = T)
pretty_eval_table <- formattable::formattable(eval_scores,
                            align = c("l", "c", "c", "c", "c"))

grid.arrange(pretty_table, pretty_eval_table, nrow = 1, ncol = 2)

#Clear console
cat("\014")