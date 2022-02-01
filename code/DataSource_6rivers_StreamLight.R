##==============================================================================
## Script for compiling and formatting data with stream light for stan
## First Year (within sample)
## Code author: J.R. Blaszczak
##==============================================================================

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- readRDS("./rds files/NWIS_6site_subset_SL.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("./rds files/NWIS_6siteinfo_subset_SL.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("./rds files/StreamLight_daily_6riv_all.rds")
colnames(SL)[colnames(SL) == "Date"] <- "date"

## Join data and StreamLight
data <- left_join(data, SL, by=c("site_name", "date"))

## Change river names to short names
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_02336526"="Proctor Creek, GA",
                                                                               "nwis_01649190"="Paint Branch, MD",
                                                                               "nwis_07191222"="Beaty Creek, OK",
                                                                               "nwis_01608500"="S. Br. Potomac River, WV",
                                                                               "nwis_11044000"="Santa Margarita River, CA",
                                                                               "nwis_08447300"="Pecos River, TX"))

## Order for figures by stream order
site_order_list <- c("Proctor Creek, GA",
                     "Paint Branch, MD",
                     "Beaty Creek, OK",
                     "S. Br. Potomac River, WV",
                     "Santa Margarita River, CA",
                     "Pecos River, TX")

## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()
## Select the first of two years
data <- rbind(data[which(data$site_name == "nwis_02336526" & data$year %in% c(2015)),],
              data[which(data$site_name == "nwis_01649190" & data$year %in% c(2010)),],
              data[which(data$site_name == "nwis_07191222" & data$year %in% c(2009)),],
              data[which(data$site_name == "nwis_01608500" & data$year %in% c(2012)),],
              data[which(data$site_name == "nwis_11044000" & data$year %in% c(2015)),],
              data[which(data$site_name == "nwis_08447300" & data$year %in% c(2012)),])


## small: nwis_02336526 2015,2016 (Order 2; PROCTOR CREEK AT JACKSON PARKWAY, AT ATLANTA, GA)
## small: nwis_01649190 2010,2011 (Order 2; PAINT BRANCH NEAR COLLEGE PARK, MD)
## mid: nwis_07191222 2009,2010 (Order 3; Beaty Creek near Jay, OK) 
## mid: nwis_01608500 2012,2013 (Order 5; SOUTH BRANCH POTOMAC RIVER NEAR SPRINGFIELD, WV) 
## large: nwis_11044000 2015,2016 (Order 6; SANTA MARGARITA R NR TEMECULA CA) 
## large: nwis_08447300 2012,2013 (Order 7: Pecos Rv at Brotherton Rh nr Pandale, TX) 

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## split list by ID
l <- split(data, data$site_name)

rel_LQT <- function(x){
  x$light_rel_PPFD <- x$light/max(x$light)
  x$light_rel_PAR <- x$PAR_surface/max(x$PAR_surface)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)

  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))
df <- dat

rm(data, l, SL, data_siteyears, dat)

