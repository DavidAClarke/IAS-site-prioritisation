#Jaccard similarities among site sensitivities
#Variant names
# all_names <- c("species_CAZ",
#                "species_wgt_CAZ",
#                "species_wgt_CAZ_KBA",
#                "species_wgt_RAN",
#                "species_area_CAZ",
#                "species_area_wgt_CAZ",
#                "species_area_wgt_CAZ_KBA",
#                "species_area_wgt_RAN")

all_names <- names(full_rank_stack)

#Top 2%
top_two <- calculate_jaccards(stack(full_rank_stack), x.min = 0.98, x.max = 1.0, 
                              y.min = 0.98,y.max = 1.0,all_names)

#Top 5%
top_five <- calculate_jaccards(stack(full_rank_stack), x.min = 0.95, x.max = 1.0, 
                               y.min = 0.95,y.max = 1.0,all_names)

#Top 10%
top_ten <- calculate_jaccards(stack(full_rank_stack), x.min = 0.9, x.max = 1.0, 
                              y.min = 0.9,y.max = 1.0,all_names)

#Top 25%
top_twentyfive <- calculate_jaccards(stack(full_rank_stack), x.min = 0.75, x.max = 1.0, 
                                     y.min = 0.75,y.max = 1.0,all_names)

#Top 50%
top_fifty <- calculate_jaccards(stack(full_rank_stack), x.min = 0.5, x.max = 1.0, 
                                y.min = 0.5,y.max = 1.0,all_names)

#Top 80%
top_eighty <- calculate_jaccards(stack(full_rank_stack), x.min = 0.2, x.max = 1.0, 
                                 y.min = 0.2,y.max = 1.0,all_names)

#Total
total <- calculate_jaccards(stack(full_rank_stack), x.min = 0.0, x.max = 1.0, 
                            y.min = 0.0,y.max = 1.0,all_names)

write.csv(top_two, file = here(dirname(here()), "data", "jaccard", "jaccard_two.csv"), na = "-", row.names = T)
write.csv(top_five, file = here(dirname(here()),"data", "jaccard", "jaccard_five.csv"), na = "-", row.names = T)
write.csv(top_ten, file = here(dirname(here()),"data", "jaccard", "jaccard_ten.csv"), na = "-", row.names = T)
write.csv(top_twentyfive, here(dirname(here()),"data", "jaccard", "jaccard_twentyfive.csv"), na = "-", row.names = T)
write.csv(top_fifty, file = here(dirname(here()),"data", "jaccard", "jaccard_fifty.csv"), na = "-", row.names = T)
write.csv(top_eighty, file = here(dirname(here()),"data", "jaccard", "jaccard_eighty.csv"), na = "-", row.names = T)
write.csv(total, file = here(dirname(here()), "data","jaccard", "jaccard_total.csv"), na = "-", row.names = T)