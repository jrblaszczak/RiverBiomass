## WI Data source

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- read.csv("../data/NWIS_WI_data.csv", header=T)
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- read.csv("../data/NWIS_WI_sitedata.csv", header=T)

## Figure out how many days of data per site per year
data$year <- year(data$date)
data %>%
  group_by(site_name,year) %>%
  tally()

## subset to only 2010 - two locations along Black Earth Creek
data <- data[which(data$year == "2010"),]

## Change GPP sd column name
colnames(data)[which(colnames(data) == "GPP_daily_sd_name")] <- "GPP_sd"

## Assign IDs
data$ID <- NA
data[which(data$site_name == "nwis_05406457"),]$ID <- "BlackEarth_CrossPlains"
data[which(data$site_name == "nwis_05406500"),]$ID <- "BlackEarth_BlackEarth"

## split list by ID
l <- split(data, data$ID)

rel_LQT <- function(x){
  x$light_rel <- x$light/max(x$light)
  #x$daylength_rel <- x$daylength_hrs/max(x$daylength_hrs)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- (x$Q-mean(x$Q))/sd(x$Q)
  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))

