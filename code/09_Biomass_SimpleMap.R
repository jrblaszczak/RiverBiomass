## Simple map of sites used in persistence plots

lapply(c("plyr","dplyr","ggplot2","cowplot",
         "rstan", "gridExtra","usmap",
         "wesanderson","ggrepel","ggsn","ggspatial"), require, character.only=T)


site_info <- readRDS("./rds files/NWIS_6siteinfo_subset_SL.rds")
## Change river names to short names
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_02336526"="Proctor Creek, GA",
                                                                               "nwis_01649190"="Paint Branch, MD",
                                                                               "nwis_07191222"="Beaty Creek, OK",
                                                                               "nwis_01608500"="S. Br. Potomac River, WV",
                                                                               "nwis_11044000"="Santa Margarita River, CA",
                                                                               "nwis_08447300"="Pecos River, TX"))

site_info_map <- site_info[,c("lon","lat","nwis_id","short_name","NHD_STREAMORDE")]

sp_transformed <- usmap_transform(site_info_map)

col <- wes_palette("Zissou1", 6, type = "continuous")


## Plot
plot_usmap(exclude = c("AK","HI")) +
  geom_point(data = sp_transformed, aes(x = x, y = y,
                                        fill=short_name),
             shape=21, size=3)+
  geom_label_repel(data = sp_transformed, aes(x = x, y = y, label = short_name,
                                        fill = short_name))+
  theme(legend.position = "none")+
  ggspatial::annotation_scale(location = 'tl')+
  annotation_north_arrow(location = "br")
  