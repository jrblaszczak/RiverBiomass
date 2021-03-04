## 6 rivers data source

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- readRDS("../rds files/NWIS_6site_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("../data/NWIS_6siteinfo_subset.rds")

## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()
## Select the first of two years
data <- rbind(data[which(data$site_name == "nwis_08180700" & data$year %in% c(2010)),],
              data[which(data$site_name == "nwis_10129900" & data$year %in% c(2015)),],
              data[which(data$site_name == "nwis_03058000" & data$year %in% c(2014)),],
              data[which(data$site_name == "nwis_01649500" & data$year %in% c(2012)),],
              data[which(data$site_name == "nwis_14211010" & data$year %in% c(2012)),],
              data[which(data$site_name == "nwis_02234000" & data$year %in% c(2013)),])

## Set any GPP < 0 to a small value close to 0
data[which(data$GPP < 0),]$GPP <- sample(exp(-6):exp(-4), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## visualize
ggplot(data, aes(date, GPP))+
  geom_point()+geom_line()+
  facet_wrap(~site_name,scales = "free_x")

## split list by ID
l <- split(data, data$site_name)

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


rm(data,l, data_siteyears)
