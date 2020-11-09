## WI Data source

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

## Change GPP sd column name
colnames(data)[which(colnames(data) == "GPP_daily_sd_name")] <- "GPP_sd"

## Assign IDs
data$ID <- NA
data[which(data$site_name == "nwis_14211010"),]$ID <- "Clackamas_OR"

## Create a GPP SD
data$GPP_sd <- ((data$GPP - data$GPP.lower) + (data$GPP.upper - data$GPP))/2

## Set any GPP < 0 to 0
data[which(data$GPP < 0 & data$GPP > -0.5),]$GPP <- sample(exp(-6):exp(-4), 1)
data[which(data$GPP < -0.5),]$GPP <- sample(exp(-6):exp(-4), 1) ## eventually change to NA when figure out how to do so

## visualize
ggplot(data, aes(date, GPP))+geom_point()+geom_line()

## split list by ID
l <- split(data, data$ID)

rel_LQT <- function(x){
  x$light_rel <- x$light/max(x$light)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)
  
  #x$std_light <- (x$light-mean(x$light))/sd(x$light)
  #x$std_temp <- (x$temp-mean(x$temp))/sd(x$temp)
  #x$tQ <- (x$Q-mean(x$Q))/sd(x$Q)
  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))

