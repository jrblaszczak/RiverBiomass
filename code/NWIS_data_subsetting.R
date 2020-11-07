## Wisconsin Data

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2"), require, character.only=T)

###################
## Import data
##################
NWIS <- read.table("../data/daily_predictions.tsv", sep='\t', header = TRUE)
head(NWIS)
NWIS$date <- as.POSIXct(as.character(NWIS$date), format="%Y-%m-%d")

head(NWIS)

## Subset columns and sites
NWIS_sub <- NWIS[,c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                    "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                    "temp.water","discharge","shortwave","velocity")]
colnames(NWIS_sub) <- c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                        "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                        "temp","Q","light","velocity")

## Import site information
NWIS_site_data <- read.csv("../data/site_data.csv", header = TRUE)


#############
## Subset
############
## Wisconsin-specific
WI_locations <- subset(NWIS_site_data, nwis_id > 5400000 & nwis_id < 5410000)
NWIS_WI <- NWIS_sub[which(NWIS_sub$site_name %in% WI_locations$site_name),]
NWIS_WI_list <- split(NWIS_WI, NWIS_WI$site_name)

ggplot(NWIS_WI, aes(date, GPP, color=site_name))+geom_line()
ggplot(NWIS_WI_list$nwis_05406457, aes(date, GPP, color=site_name))+geom_line()

## Export
write.csv(WI_locations,"../data/NWIS_WI_sitedata.csv")
write.csv(NWIS_WI,"../data/NWIS_WI_data.csv")








