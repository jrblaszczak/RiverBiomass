## Table S1 - Site and time series information

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

## save general info
write.csv(sub,"../figures and tables/2022 Tables/TableS_SiteInfo_6sites.csv")


################################################
## Table S2 - Site and time series information
################################################
## Import data on the number of days per year & model diagnostics
site_subset <- readRDS("./rds files/NWIS_6site_subset_SL.rds")
TS_site_subset <- readRDS("./rds files/NWIS_6siteinfo_subset_SL.rds")
site_subset_numdays <- readRDS("./rds files/NWIS_6site_Ndays_SL.rds")



