## Subset data from hypoxia database that is already linked to NHD

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table"), require, character.only=T)

## Import site data from Appling
setwd("../data")
site <- fread("site_data.csv")
names(site)

## Import hypoxia dataset and subset to Appling
hyp <- fread("GRDO_GEE_HA_NHD_2021_02_07.csv")
PC <- hyp[which(hyp$DB_Source == "PC"),]
names(PC)

## Merge both based on siteID
head(site$site_name); head(PC$SiteID)
colnames(PC)[which(colnames(PC) == "SiteID")] <- "site_name"

df <- left_join(site, PC, by="site_name")

## Export
write.csv(df, "PC_site_attribs.csv")

##########################################
## Re-import for site selection
##########################################
setwd("../data")
df <- fread("PC_site_attribs.csv")
colnames(df)

sub <- df[,c("site_name","ORD_STRA","NHD_STREAMORDE",
              "Start_time","End_time","n_time",
              "US_state")]

## Need multiple years of data
sub$Start_year <- year(sub$Start_time)
sub$End_year <- year(sub$End_time)
sub$n_years <- sub$End_year - sub$Start_year

s <- sub[which(sub$n_years >=3),]

## Compare to Savoy sites
















