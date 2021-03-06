## 10 rivers data source

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- read.csv("../data/NWIS_Psim_data.csv", header=T)
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- read.csv("../data/NWIS_Psim_sitedata.csv", header=T)

## Figure out how many days of data per site per year
data$year <- year(data$date)
data$site_year <- paste(data$site_name, data$year, sep="_")
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()
# only the year with the most data per site
data_siteyears <- data_siteyears %>%
  group_by(site_name)%>%
  filter(n == max(n))
data_siteyears$site_year <- paste(data_siteyears$site_name, data_siteyears$year, sep="_")
# Tie in nwis_08181500
data_siteyears <- data_siteyears[-which(data_siteyears$site_name == "nwis_08181500" & data_siteyears$year == "2013"),]

## subset to site years
data <- subset(data, site_year %in% data_siteyears$site_year)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## Set any GPP < 0 to a small value close to 0
data[which(data$GPP < 0 & data$GPP > -0.5),]$GPP <- sample(exp(-6):exp(-4), 1)
data[which(data$GPP < -0.5),]$GPP <- sample(exp(-6):exp(-4), 1) ## eventually change to NA when figure out how to do so

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