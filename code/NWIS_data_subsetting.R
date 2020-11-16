## NWIS Subsetting Time Series from Appling et al. Powell Data Set

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
## Wisconsin-specific ##
WI_locations <- subset(NWIS_site_data, nwis_id > 5400000 & nwis_id < 5410000)
NWIS_WI <- NWIS_sub[which(NWIS_sub$site_name %in% WI_locations$site_name),]
NWIS_WI_list <- split(NWIS_WI, NWIS_WI$site_name)

ggplot(NWIS_WI, aes(date, GPP, color=site_name))+geom_line()
ggplot(NWIS_WI_list$nwis_05406457, aes(date, GPP, color=site_name))+geom_line()

# Export
write.csv(WI_locations,"../data/NWIS_WI_sitedata.csv")
write.csv(NWIS_WI,"../data/NWIS_WI_data.csv")


## Oregon - Clackamas River ##
OR_location <-NWIS_site_data[which(NWIS_site_data$site_name == "nwis_14211010"),]
NWIS_OR <- NWIS_sub[which(NWIS_sub$site_name %in% OR_location$site_name),]

ggplot(NWIS_OR, aes(date, GPP))+geom_point()

# Export
write.csv(OR_location,"../data/NWIS_OR_sitedata.csv")
write.csv(NWIS_OR,"../data/NWIS_OR_data.csv")


## Other sites across the country for P simulation ##
# Which locations have >300 days in a year
NWIS_sub$year <- year(NWIS_sub$date)
NWIS_sub_siteyears <- NWIS_sub %>%
  group_by(site_name, year) %>%
  count()
NWIS_sub_siteyears <- NWIS_sub_siteyears[which(NWIS_sub_siteyears$n > 350),]

#locations from Phil's visualizations
Psim_locations <- NWIS_site_data[which(NWIS_site_data$site_name %in% c("nwis_14211010", #2010
                                                                     "nwis_01608500", #2012
                                                                     "nwis_02156500",
                                                                     "nwis_02168504",
                                                                     "nwis_03298150",
                                                                     "nwis_04121944",
                                                                     "nwis_04136000",
                                                                     "nwis_05435950",
                                                                     "nwis_06711565",
                                                                     "nwis_08181500")),]
NWIS_Psim <- NWIS_sub[which(NWIS_sub$site_name %in% Psim_locations$site_name),]

ggplot(NWIS_Psim, aes(date, GPP))+geom_point()

# Export
write.csv(Psim_locations,"../data/NWIS_Psim_sitedata.csv")
write.csv(NWIS_Psim,"../data/NWIS_Psim_data.csv")



















