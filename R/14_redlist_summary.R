#Summary information for red list

#################################################################################
##Summary information for all species used (n = 5113)
species_info <- read.csv(file.path("zonation", "maxent_spp_list.csv"))
sum_info <- species_info %>%
  mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern",
                                             "Near Threatened", "Vulnerable",
                                             "Endangered", 
                                             "Critically Endangered"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", 
                                             "CR"))) %>%
  group_by(classGroup, redlistCategory) %>%
  summarise(species_per_cat = n()) %>%
  mutate(species_per_group = sum(species_per_cat))

IAS_threats <- species_info %>%
  dplyr::filter(code == "8.1.1" | code == "8.1.2") %>%
  mutate(redlistCategory = factor(redlistCategory, 
                                  levels = c("Data Deficient", "Least Concern",
                                             "Near Threatened", "Vulnerable",
                                             "Endangered", 
                                             "Critically Endangered"),
                                  labels = c("DD", "LC", "NT", "VU", "EN", 
                                             "CR"))) %>%
  group_by(classGroup, redlistCategory) %>%
  summarise(species_per_cat = n()) %>%
  mutate(species_per_group = sum(species_per_cat)) %>%
  mutate(prop.total = species_per_cat/species_per_group)

IAS_species <- species_info %>%
  dplyr::select(ias) %>%
  count(ias) %>%
  drop_na(ias)

#Red List plot
#IUCN colours
mycols <- c("#808080", "#008000", "#ADFF2F", "#FFFF00", "#FFA500", "#FF0000") 

#All data & all threats
All_species <- ggplot(sum_info, aes(x = redlistCategory, y = species_per_cat, 
                                    fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  facet_wrap(~classGroup, nrow = 2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(trans = "log1p") +
  ylab("Number of species (log + 1)") +
  xlab("IUCN Red List category")

#Species threatened by IAS
IAS_threatend_species <- ggplot(IAS_threats, aes(x = redlistCategory, 
                                                 y = species_per_cat, 
                                                 fill = redlistCategory)) +
  geom_bar(stat = "identity") +
  facet_wrap(~classGroup, nrow = 2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(trans = "log1p") + #could also change to proportion of total
  ylab("Number of species (log + 1)") +
  xlab("IUCN Red List category")

#Alternative stacked plot
IAS_threats_stack <- IAS_threats %>%
  mutate(classGroup = factor(classGroup, levels=c("Mammal", "Bird", "Reptile", 
                                                  "Amphibian", "Fish", 
                                                  "Invertebrate", "Plant")))
IAS_threatend_species_stack <- ggplot(IAS_threats_stack, 
                                      aes(x = classGroup, y = species_per_cat, 
                                          fill = redlistCategory)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = mycols, name = "Red List\nCategory") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Proportion of species") +
  xlab("Group")

#Red List index - IAS impact
IAS_impact <- read.csv(file.path("ISO_BL_AUS_aggregated_general.csv"))
RLI_IAS <- ggplot(IAS_impact, aes(x = year, y = rli)) + 
  geom_line(colour = "blue") + 
  geom_ribbon(aes(ymin = qn05, ymax = qn95), linetype = 2, alpha = 0.1) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.75, 1)) +
  #scale_x_discrete(breaks = IAS_impact$year) +
  ylab("Red List Index") +
  xlab("Year")