## OR Data source

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- read.csv("../data/NWIS_OR_data.csv", header=T)
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- read.csv("../data/NWIS_OR_sitedata.csv", header=T)

## Figure out how many days of data per site per year
data$year <- year(data$date)
data %>%
  group_by(site_name,year) %>%
  tally()

## subset to year with more than 300 days
data <- subset(data, year >= 2009 & year <= 2013)

## Assign IDs
data$ID <- NA
data[which(data$site_name == "nwis_14211010"),]$ID <- "Clackamas_OR"

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## Set any GPP < 0 to a small value close to 0
data[which(data$GPP < 0 & data$GPP > -0.5),]$GPP <- sample(exp(-6):exp(-4), 1)
data[which(data$GPP < -0.5),]$GPP <- sample(exp(-6):exp(-4), 1) ## eventually change to NA when figure out how to do so

## Import turbidity data
OR_turb <- read.csv("../data/Clackamas_daily_turb.csv") ## only for 2010 for now
OR_turb$date <- as.POSIXct(as.character(OR_turb$date),format="%Y-%m-%d")
OR_turb$mean_daily_turb <- as.numeric(as.character(OR_turb$mean_daily_turb))
#merge with data
data <- left_join(data, OR_turb, by="date")

## visualize
ggplot(data, aes(date, GPP))+geom_point()+geom_line()

## split list by ID
l <- split(data, data$ID)

rel_LQT <- function(x){

  # Relativize everything by the max
  x$light_rel <- x$light/max(x$light)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)
  
  # Standardize light to mean
  x$std_light <- (x$light-mean(x$light))/sd(x$light)

  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))

rm(data,l, OR_turb)
