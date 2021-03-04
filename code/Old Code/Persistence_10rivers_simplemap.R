## Simple map of sites used in persistence plots

lapply(c("plyr","dplyr","ggplot2","cowplot",
         "rstan", "shinystan","gridExtra","usmap"), require, character.only=T)


site_info <- read.csv("../data/NWIS_Psim_sitedata.csv", header=T)

site_info_map <- site_info[,c("lon","lat","nwis_id")]

sp_transformed <- usmap_transform(site_info_map)

col <- wes_palette("Zissou1", 10, type = "continuous")

plot_usmap() +
  geom_point(data = sp_transformed, aes(x = lon.1, y = lat.1), size=3, color= col)
