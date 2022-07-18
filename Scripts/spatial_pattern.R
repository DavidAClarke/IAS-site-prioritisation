# Looking for spatial pattern of highest sensitive sites
df_new <- spat_priority_dist(priority_dist, 4)

#E.g.variant 1. Do for all 4 in df
df_1 <- data.frame(df_new[,1], df_new[,5], df_new[,6])
df_1 <- df_1 %>% dplyr::filter(df_1[,1] == 1)
p <- ppp(df_1[,2], df_1[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")

df_2 <- data.frame(df_new[,2], df_new[,5], df_new[,6])
df_2 <- df_2 %>% dplyr::filter(df_2[,1] == 1)
p <- ppp(df_2[,2], df_2[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")

df_3 <- data.frame(df_new[,3], df_new[,5], df_new[,6])
df_3 <- df_3 %>% dplyr::filter(df_3[,1] == 1)
p <- ppp(df_3[,2], df_3[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")

df_4 <- data.frame(df_new[,4], df_new[,5], df_new[,6])
df_4 <- df_4 %>% dplyr::filter(df_4[,1] == 1)
p <- ppp(df_4[,2], df_4[,3], xrange = c(-2131603,2443397),yrange = c(-4947747,-1097747))
clarkevans.test(p, alternative = "clustered", correction = "Donnelly")

#create data frame of priority and distance
#dist_df <- data.frame(priority = priority, distance = dists)

# #correlation
# cor.test(priority_dist[,1], dists, method = "spearman")
# cor.test(priority_dist[,2], dists, method = "spearman")
# cor.test(priority_dist[,3], dists, method = "spearman")
# cor.test(priority_dist[,4], dists, method = "spearman")
# 
# #scatterplot
# plot(dists, priority, cex = 0.1)