## Table 1 - Site and time series information

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork"), require, character.only=T)

#################################################
## Import data, site info, and site diagnostics
#################################################

## Source data
source("DataSource_9rivers.R")

# Subset site info
site_info <- site_info[which(site_info$site_name %in% names(df)),]






