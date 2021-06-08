## Table 1 - Site and time series information

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork"), require, character.only=T)

#################################################
## Import data, site info, and site diagnostics
#################################################

## Source data
source("DataSource_6rivers_StreamLight.R")

## Subset relevant columns
colnames(site_info)
sub <- site_info[,c("short_name","nwis_id","lat","lon","NHD_STREAMORDE","struct.canal_flag","struct.dam_flag","struct.npdes_flag")]


################################################
## Import data on the number of days per year
################################################
saveRDS(site_subset_numdays,"./rds files/NWIS_6site_Ndays_SL.rds")





