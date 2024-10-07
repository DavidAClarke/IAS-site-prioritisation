## zonation v5 figures

## Performance
sum_curves <- read.csv(here("zonation", "species_scenarios", "species_weight", 
                            "output", "summary_curves.csv"), sep = " ")
ggplot(sum_curves, aes(x = rank, y = mean)) +
  geom_line() +
  geom_line(aes(x = rank, y = min), linetype = 2) +
  geom_line(aes(x = rank, y = max), linetype = 2)

group_curves <- read.csv(here("zonation", "species_scenarios", "species_weight", 
                            "output", "group_curves.csv"), sep = " ") %>%
  dplyr::select(!X)

group_curves_long <- group_curves %>%
  pivot_longer(cols = !rank, values_to = "Coverage") %>%
  mutate(group = case_when(str_detect(name, "1") ~ "Invertebrate",
                           str_detect(name, "2") ~ "Fish",
                           str_detect(name, "3") ~ "Plant",
                           str_detect(name, "4") ~ "Reptile",
                           str_detect(name, "6") ~ "Mammal",
                           str_detect(name, "7") ~ "Amphibian",
                           str_detect(name, "8") ~ "Fungi",
                           str_detect(name, "9") ~ "Bird")) %>%
  mutate(name = gsub('.{2}$', '', name))

group_curves_long %>%
  dplyr::filter(name == "mean") %>%
  ggplot(aes(x = rank, y = Coverage, color = group)) +
    geom_line(linewidth = 1.3) +
  scale_color_manual(name = "Group",
                     values = RColorBrewer::brewer.pal(8, "Set2")) +
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
    xlab("Priority rank")
