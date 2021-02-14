## Link Appling sites with NHDv2

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","sf","nhdplusTools"), require, character.only=T)

## Import site data
setwd("../data")
site <- fread("site_data.csv")
names(site)

## Download NHD Plus data
outdir <- "~/GitHub/RiverBiomass/data/NHDPlus"
download_nhdplusv2(outdir,url = paste0("https://s3.amazonaws.com/edap-nhdplus/NHDPlusV21/",
                          "Data/NationalData/NHDPlusV21_NationalData_Seamless",
                          "_Geodatabase_Lower48_07.7zip"))

## Pair sites with NHDplus IDs
test <- list(featureSource = "nwissite", featureID = "USGS-01021050")
discover_nhdplus_id(nldi_feature = test)

char <- discover_nldi_characteristics()




